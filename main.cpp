/*
        2D Finite volume electron MHD solver in Cartesian coordinates
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>
#include <string>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

#include <span>
#include <limits>
#include <list>
#include <filesystem>
namespace fs = std::filesystem;

#include "boost/multi_array.hpp"

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "common.h"
#include "initial_conditions.h"
#include "boundary_conditions.h"
#include "conservation_checks.h"
#include "microphysics.h"
#include "field_evolution.h"

const H5std_string FILE_EXT( ".h5" );

void RK_Step(VectorField & B, VectorField & E, VectorField & J, VectorField & vc, VectorField & B_np1, TransCoeffs & tC, const Domain & domain, const Process & process, const BandBCParams & bparams, double t);

int main(int argc, char **argv)
{
    int world_rank, num_procs;
    /* Initialize the infrastructure necessary for communication */
    MPI_Init(&argc, &argv);
    /* Identify this process */
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    /* Find out how many total processes are active */
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    /*
        Load information from "SimSetup.in"
    */
    SimParams simparams; //object of SimParams class containing simulation parameters (see common.h for full list)
    BandBCParams bparams; //object of BandBCParams bparams class containing initial field and boundary conditions (see common.h for full list)
    load_params(simparams, bparams, world_rank); //load simulation setup parameters from SimSetup.in text file
    const H5std_string FILE_NAME( simparams.OutputFile );

    /*
        Sets number of ghost cells needed at end of each simulation domain. Using
        monotonised central difference limiter (MC) method, so need two at each edge.
    */
    const size_t N_GC = 2;

    //Define problem domain in reduced units
    size_t Nx = simparams.Nx; //number of finite volume cells in x-direction, excluding ghost cells
    size_t Ny = simparams.Ny; //number of finite volume cells in y-direction, excluding ghost cells
    double x_min = simparams.x_min, x_max = simparams.x_max;
    double y_min = simparams.y_min, y_max = simparams.y_max;

    if( (Nx%2 !=0 || Ny%2 !=0) && world_rank == 0 ){
        printf("WARNING: Please choose number of finite volume cells in each direction to be even numbers\n");
    }
    if( (Ny % (size_t)num_procs != 0) && world_rank == 0 ){
        printf("WARNING: Number of finite volume cells cannot be evenly divided among processes\n");
    }

    double Lx = x_max-x_min; //extent of simulation domain in x-direction
    double Ly = y_max-y_min; //extent of simulation domain in y-direction
    double Deltay = Ly/double(Ny); //size of each cell in y-direction in reduced units

    double t = 0.; //current simulation time. Initialize to 0.
    const double t_max = simparams.t_max; //maximum simulation time

    /*
        Define finite volume cell centres and spacing between them.
    */
    std::vector<double> xFull; //Central points of cells over full domain
    std::vector<double> Deltax; //size of each cell in x-direction in reduced units
    Deltax.insert(Deltax.end(), { 0, 0 }); //ghost cell Deltax values. Will fill these with correct values later.
    if(simparams.varying_mesh == true){
        double x_scale = double(Nx)/2; //determines how far apart the cell centres are spaced in x. Larger values are more evenly spaced, while smaller values give more cells near the upper boundary.
        for(size_t i = 0; i<Nx-1; i++){
            if(i == 0){
                xFull.push_back( 0.5*x_min + 0.5*( x_min + Lx*( 1.-exp(-(double(i+1))/x_scale) )/( 1.-exp(-(double(Nx-1))/x_scale) ) ) ); //divides domain into Nx equally-spaced cells with x as their central points
                Deltax.push_back( 2.*( xFull.back() - x_min ) );
            }
            else{
                xFull.push_back( 0.5*( x_min + Lx*( 1.-exp(-(double(i+1))/x_scale) )/( 1.-exp(-(double(Nx-1))/x_scale) ) ) + 0.5*( xFull.back() + Deltax.back()/2. ) );
                Deltax.push_back( 2.*( xFull.back() - xFull[xFull.size()-2] ) - Deltax.back() );
            }
        }
        Deltax.push_back( Deltax.back() );
        xFull.push_back( xFull.back() + Deltax.back() );
    }
    else{
        for(size_t i = 0; i<Nx; i++){
            xFull.push_back( x_min + Lx*(double(i)+0.5)/double(Nx-1) ); //divides domain into Nx equally-spaced cells with x as their central points. Center of the final cell is outside the simulation domain!
            Deltax.push_back( Lx/double(Nx-1) ); //size of each cell in x-direction in reduced units. Subtract off one from denominator to account for non-periodic cells
        }
    }
    //Fill ghost cells in Deltax appropriately i.e., mirrored across the simulation boundaries so Deltax[0]=Deltax[3], Deltax[1]=Deltax[2] for N_GC=2
    for(size_t i=0; i<N_GC; i++){
        Deltax[N_GC-1-i] = Deltax[N_GC+i];
        Deltax.push_back( Deltax[Deltax.size()-3-2*i] );
    }
    std::vector<double> yFull; //Central points of cells over full domain
    yFull = linspace(y_min+Deltay/2., y_max-Deltay/2., Ny); //divides domain into Ny equally-spaced cells with y as their central points

    /*
        Cartesian decomposition of the domain using MPI
    */
    int n_dims = 1; //number of spatial dimensions in the decomposition of the problem (not necessarily the spatial dimension of the problem)
    int dims[1] = {num_procs};
    int isperiodic[1] = {true};
    int reorder = true; //whether to allow MPI to reorder the partial domains
    MPI_Comm comm1D;
    MPI_Cart_create( MPI_COMM_WORLD, n_dims, dims, isperiodic, reorder, &comm1D ); //creates communicator comm1D

    int MyID;
    MPI_Comm_rank( comm1D, &MyID );
    int nbrleft, nbrright;
    MPI_Cart_shift( comm1D, 0, 1, &nbrleft, &nbrright );

    size_t MyS = 0, MyE = 0; //start and end indices of each partial domain, excluding ghost cells. Determined using MPE_Decomp1D
    MPE_Decomp1D(Ny, (size_t)num_procs, (size_t)MyID, MyS, MyE); //decomposes domain into num_procs partial domains in the y-direction with Ny cells

    std::vector<int> starts(num_procs); //Vector of the starting indices (in decomposed domain direction) of each partial domain
    std::vector<int> Ny_locs(num_procs); //Vectors of the extent (in decomposed domain direction) of each partial domain
    //Send starts and Ny_locs to each process.
    MPI_Allgather(&MyS, 1, MPI_INT, starts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    int Ny_loc = int(MyE - MyS);
    MPI_Allgather(&Ny_loc, 1, MPI_INT, Ny_locs.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<double> x(xFull.begin(), xFull.end() ); //Local domain cell mid-points in x-direction. Unsplit
    std::vector<double> y(yFull.begin() + MyS, yFull.begin() + MyE ); //Local domain cell mid-points in y-direction. Split based on the start and end indices computed by MPE_Decomp1D

    /*
        Define physical parameters and simulation domain
    */

    const double rho_cutoff = simparams.rho_cutoff; //cutoff (energy) density in g/cm^3 (lowest density to include in simulation domain).
    DataInterpolation EOS_Interps; //object of DataInterpolation class to contain EOS interpolating functions as function of "radial" coordinate x
    load_EOS(simparams.CrustEOS,rho_cutoff,EOS_Interps);

    ScalarField T(boost::extents[Nx+2*N_GC][MyE-MyS+2*N_GC]);

    //Declare object containing ScalarFields/RadialScalarFields for transport coefficients
    //This object are referenced throughout the code.
    TransCoeffs tC(Nx+2*N_GC, MyE-MyS+2*N_GC);

    ScalarField n_e(boost::extents[Nx+2*N_GC][MyE-MyS+2*N_GC]);
    ScalarField dn_eInvdx(boost::extents[Nx+2*N_GC][MyE-MyS+2*N_GC]);
    ScalarField sigma(boost::extents[Nx+2*N_GC][MyE-MyS+2*N_GC]);

    RadialScalarField A(boost::extents[Nx+2*N_GC]); //mass number
    RadialScalarField Z(boost::extents[Nx+2*N_GC]); //atomic number
    RadialScalarField n_i(boost::extents[Nx+2*N_GC]); //ion number density in fm^-3
    RadialScalarField n_b(boost::extents[Nx+2*N_GC]); //baryon number density in fm^-3
    RadialScalarField rho(boost::extents[Nx+2*N_GC]); //mass-energy density in g/cm^3

    for(size_t i=0; i<Nx+2*N_GC; i++){
        if( i < N_GC ){
            A[i] = gsl_spline_eval(EOS_Interps.A_spline, EOS_Interps.R_cc, EOS_Interps.A_acc);
            Z[i] = gsl_spline_eval(EOS_Interps.Z_spline, EOS_Interps.R_cc, EOS_Interps.Z_acc);
            n_b[i] = gsl_spline_eval(EOS_Interps.n_b_spline, EOS_Interps.R_cc, EOS_Interps.n_b_acc);
            n_i[i] = gsl_spline_eval(EOS_Interps.n_i_spline, EOS_Interps.R_cc, EOS_Interps.n_i_acc);
            rho[i] = gsl_spline_eval(EOS_Interps.rho_spline, EOS_Interps.R_cc, EOS_Interps.rho_acc);
        }
        else if( i > Nx+N_GC-2 ){
            A[i] = gsl_spline_eval(EOS_Interps.A_spline, EOS_Interps.R_rhocutoff, EOS_Interps.A_acc);
            Z[i] = gsl_spline_eval(EOS_Interps.Z_spline, EOS_Interps.R_rhocutoff, EOS_Interps.Z_acc);
            n_b[i] = gsl_spline_eval(EOS_Interps.n_b_spline, EOS_Interps.R_rhocutoff, EOS_Interps.n_b_acc);
            n_i[i] = gsl_spline_eval(EOS_Interps.n_i_spline, EOS_Interps.R_rhocutoff, EOS_Interps.n_i_acc);
            rho[i] = gsl_spline_eval(EOS_Interps.rho_spline, EOS_Interps.R_rhocutoff, EOS_Interps.rho_acc);
        }
        else{
            A[i] = gsl_spline_eval(EOS_Interps.A_spline, EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*(x[i-N_GC]-x_min), EOS_Interps.A_acc);
            Z[i] = gsl_spline_eval(EOS_Interps.Z_spline, EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*(x[i-N_GC]-x_min), EOS_Interps.Z_acc);
            n_b[i] = gsl_spline_eval(EOS_Interps.n_b_spline, EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*(x[i-N_GC]-x_min), EOS_Interps.n_b_acc);
            n_i[i] = gsl_spline_eval(EOS_Interps.n_i_spline, EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*(x[i-N_GC]-x_min), EOS_Interps.n_i_acc);
            rho[i] = gsl_spline_eval(EOS_Interps.rho_spline, EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*(x[i-N_GC]-x_min), EOS_Interps.rho_acc);
        }
        for(size_t j=0; j<MyE-MyS+2*N_GC; j++){
            T[i][j] = simparams.temperature*k_B; //temperature in MeV
            if( i < N_GC ){
//                n_e[i][j] = gsl_spline_eval(EOS_Interps.n_e_spline, EOS_Interps.R_cc, EOS_Interps.n_e_acc);
//                dn_eInvdx[i][j] = 0.;
                n_e[i][j] = gsl_spline_eval(EOS_Interps.n_e_spline, EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*x[N_GC-1-i], EOS_Interps.n_e_acc);
                dn_eInvdx[i][j] = 1./n_e[i][j];
//                dn_eInvdx[i][j] = ( 1./gsl_spline_eval(EOS_Interps.n_e_spline, EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*(x[0]+0.5*Deltax[2]-x_min), EOS_Interps.n_e_acc)
//                                    - 1./gsl_spline_eval(EOS_Interps.n_e_spline, EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*(x[0]-0.5*Deltax[2]-x_min), EOS_Interps.n_e_acc))/Deltax[i];
            }
            else if( i > Nx+2*N_GC-3 ){
                n_e[i][j] = gsl_spline_eval(EOS_Interps.n_e_spline, EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*x[2*Nx+2*N_GC-i-3], EOS_Interps.n_e_acc);
                dn_eInvdx[i][j] = 1./n_e[i][j];
//                n_e[i][j] = gsl_spline_eval(EOS_Interps.n_e_spline, EOS_Interps.R_rhocutoff, EOS_Interps.n_e_acc);
//                dn_eInvdx[i][j] = 0.;
//                dn_eInvdx[i][j] = ( 1./gsl_spline_eval(EOS_Interps.n_e_spline, EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*(x[Nx-1]+0.5*Deltax[Nx+N_GC]-x_min), EOS_Interps.n_e_acc)
//                                    - 1./gsl_spline_eval(EOS_Interps.n_e_spline, EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*(x[Nx-1]-0.5*Deltax[Nx+N_GC]-x_min), EOS_Interps.n_e_acc))/Deltax[i];
            }
            else{
                n_e[i][j] = gsl_spline_eval(EOS_Interps.n_e_spline, EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*x[i-N_GC], EOS_Interps.n_e_acc);
                dn_eInvdx[i][j] = 1./n_e[i][j];
//                n_e[i][j] = gsl_spline_eval(EOS_Interps.n_e_spline, EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*(x[i-N_GC]-x_min), EOS_Interps.n_e_acc);
//                n_e[i][j] = gsl_spline_eval_integ(EOS_Interps.n_e_spline, EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*(x[i-N_GC]-0.5*Deltax[i]-x_min),
//                                                    EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*(x[i-N_GC]+0.5*Deltax[i]-x_min), EOS_Interps.n_e_acc)/Deltax[i]; //electron number density in fm^{-3}
//                dn_eInvdx[i][j] = ( 1./gsl_spline_eval(EOS_Interps.n_e_spline, EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*(x[i-N_GC]+0.5*Deltax[i]-x_min), EOS_Interps.n_e_acc)
//                                    - 1./gsl_spline_eval(EOS_Interps.n_e_spline, EOS_Interps.R_cc + (EOS_Interps.R_rhocutoff-EOS_Interps.R_cc)/Lx*(x[i-N_GC]-0.5*Deltax[i]-x_min), EOS_Interps.n_e_acc))/Deltax[i];
            }
        }
    }

    sigmaCalc(T, n_e, A, Z, n_i, n_b, "Electrical", sigma);

//    double sigma_const = 1e25, ne_const = n_e0;//1e-3;
    double Gamma; //Coulomb coupling parameter
    for(size_t i=0; i<Nx+2*N_GC; i++){
        for(size_t j=N_GC; j<MyE-MyS+N_GC; j++){
//            tC.eta_O[i][j] = c*c*t_0/(4.*pi*sigma_const*L_0*L_0); //Ohmic diffusivity in reduced units
//            tC.eta_H[i][j] = c*B_0*t_0/(4.*pi*unit_e*ne_const*1e39*L_0*L_0); //Hall diffusivity in reduced units. 1e39 factor converts n_e from fm^-3 to cm^-3
//            tC.deta_Hdx[i][j] = 0; //x-derivative of Hall diffusivity in reduced units. 1e39 factor converts n_e from fm^-3 to cm^-3
            tC.eta_O[i][j] = c*c*t_0/(4.*pi*sigma[i][j]*L_0*L_0); //Ohmic diffusivity in reduced units
            tC.eta_H[i][j] = c*B_0*t_0/(4.*pi*unit_e*n_e[i][j]*1e39*L_0*L_0); //Hall diffusivity in reduced units. 1e39 factor converts n_e from fm^-3 to cm^-3
            tC.deta_Hdx[i][j] = c*B_0*t_0/(4.*pi*unit_e*1e39*L_0*L_0)*dn_eInvdx[i][j]; //x-derivative of Hall diffusivity in reduced units. 1e39 factor converts n_e from fm^-3 to cm^-3
//            Gamma = Z[i]*Z[i]*unit_e*unit_e/(k_Bcgs*T[i][j]*T_0*pow(3./(4.*pi*n_i[i]*1e39),1./3.));
            tC.shear_mod[i] = pow( 4.*pi/3.,1./3. )*pow(n_i[i]*1e39,4./3.)*Z[i]*Z[i]*unit_e*unit_e*0.1194/( 1. + 1.781*pow(100./Gamma,2.) ); //shear modulus in dyn/cm^2
            tC.rho[i] = rho[i]; //mass density in g/cm^3
        }
    }

    //Exchange transport coefficients on periodic dimension.
    exchng2Scalar(tC.eta_O,N_GC,comm1D,nbrleft,nbrright);
    exchng2Scalar(tC.eta_H,N_GC,comm1D,nbrleft,nbrright);
    exchng2Scalar(tC.deta_Hdx,N_GC,comm1D,nbrleft,nbrright);

    //Determine time step using cell centre spacing and Hall diffusivity
    double min_Deltax = *std::min_element( Deltax.begin(), Deltax.end() );
    double DeltaL = 1./std::sqrt( 1./(min_Deltax*min_Deltax) + 1./(Deltay*Deltay) ); //spatial step in reduced units
    double max_eta_H = *std::max_element( tC.eta_H.data(), tC.eta_H.data() + tC.eta_H.num_elements() );
    double Deltat = simparams.k_C*DeltaL*DeltaL/max_eta_H; //CFL limit

    //Define structs domain and process, containing data about the overall simulation domain and the current individual process, respectively
    struct Domain domain = {Nx, Ny, N_GC, Lx, Ly, Deltay, Deltax, x, y, Deltat, Ny_locs, starts};
    struct Process process = {simparams.RK_order, world_rank, nbrleft, nbrright, MyS, MyE, comm1D};

    //Compute initial magnetic field and create VectorField objects to hold electric field and updated magnetic field
    VectorField B(boost::extents[3][Nx+2*N_GC][MyE-MyS+2*N_GC]); //Cell-face average values of B across partial domain

    InitializeB(x, y, bparams, domain, N_GC, Deltax, Deltay, B); //Generate initial values of B components by cell-face averaging over initial functional form

    VectorField B_np1(boost::extents[3][Nx+2*N_GC][MyE-MyS+2*N_GC]); //Cell-face average values of B in reduced units at next time step
    VectorField D_np1(boost::extents[3][Nx+2*N_GC][MyE-MyS+2*N_GC]); //Cell-face average values of D in reduced units at next time step
    
    VectorField E(boost::extents[3][Nx+2*N_GC][MyE-MyS+2*N_GC]); //Cell-edge average values of E times c in reduced units (E*c/(B_0*L_0/t_0))
    VectorField J(boost::extents[3][Nx+2*N_GC][MyE-MyS+2*N_GC]); //Cell-edge average values of J in reduced units (J*L_0/B_0)
    VectorField vc(boost::extents[3][Nx+2*N_GC][MyE-MyS+2*N_GC]); //Cell-edge average values of velocity field in reduced units (vc*t_0/L_0)

    //Initializes the D and H fields
    VectorField D(boost::extents[3][Nx+2*N_GC][MyE-MyS+2*N_GC]); //Cell-face average values of D across partial domain
    VectorField H(boost::extents[3][Nx+2*N_GC][MyE-MyS+2*N_GC]); //Cell-edge average values of H

    //Initializes the charge density Rho
    ScalarField Rho(boost::extents[Nx+2*N_GC][MyE-MyS+2*N_GC]);

    //Initializes the velocity field V
    VectorField V(boost::extents[3][Nx+2*N_GC][MyE-MyS+2*N_GC]); 

    B_BoundaryConditions(B, bparams, Ny, N_GC, t, comm1D, world_rank, Ny_locs, starts, nbrleft, nbrright, domain);
    Compute_J(B, J, N_GC, domain);
    exchng2Vector(J, N_GC, comm1D, nbrleft, nbrright); //Exchange J cells in a periodic manner in the y-direction into ghost cells
    Compute_vc(B, vc, J, tC, N_GC, t, domain, bparams, world_rank);
    exchng2Vector(vc, N_GC, comm1D, nbrleft, nbrright); //Exchange vc cells in a periodic manner in the y-direction into ghost cells
    Compute_E(B, B, E, J, vc, tC, N_GC, t, domain, bparams);
    E_BoundaryConditions(E, bparams, N_GC, comm1D, nbrleft, nbrright);

    
    //Create directory for H5 files holding simulation data.
    //Also copy SimSetup.in to this directory, so can see what simulation parameters (e.g., initial conditions) were used for any particular run
    if( world_rank == 0 ){
        std::ostringstream ndir;
        ndir << FILE_NAME + "/" + FILE_NAME + "/";
        std::filesystem::create_directories(ndir.str().c_str());

        fs::path sourceFile = "SimSetup.in";
        std::ostringstream ndir_up;
        ndir_up << FILE_NAME + "/";
        fs::path targetParent = ndir_up.str();
        auto target = targetParent / sourceFile.filename();

        try{
            fs::copy_file(sourceFile, target, fs::copy_options::overwrite_existing);
        }
        catch (std::exception& e){
            std::cout << e.what();
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );

    /*
            Create h5 file
    */
    std::string RANK_NAME = std::to_string(world_rank);
    std::string file_path = FILE_NAME+"/"+FILE_NAME+"/"+FILE_NAME+"_"+RANK_NAME+FILE_EXT;
    // If a previous run left the file present or locked, try to remove it. If remove fails (locked), fall back to a unique filename.
    try{
        if( fs::exists(file_path) ){
            try{
                fs::remove(file_path);
            } catch(const std::exception &e){
                // Couldn't remove (maybe locked). Fall back to using a unique filename to avoid create failure
                std::cerr << "[HDF5] Warning: unable to remove existing file '" << file_path << "': " << e.what() << ". Using unique filename instead." << std::endl;
                std::ostringstream alt;
                alt << FILE_NAME << "/" << FILE_NAME << "/" << FILE_NAME << "_" << RANK_NAME << "_" << getpid() << FILE_EXT;
                file_path = alt.str();
            }
        }
    } catch(const std::exception &e){
        std::cerr << "[HDF5] Warning: filesystem check failed: " << e.what() << std::endl;
    }
    H5::H5File file( file_path, H5F_ACC_TRUNC );

    //Write B to output H5 file and prepare to write data from additional time steps
    hsize_t B_dims[4] = {1,3,Nx,MyE-MyS};
    hsize_t B_maxdims[4] = {H5S_UNLIMITED,H5S_UNLIMITED,H5S_UNLIMITED,H5S_UNLIMITED};
    hsize_t B_size[4] = {1,3,Nx,MyE-MyS}; //Total size of extended array. Increases with every time step. Initially should equal dimsf.
    hsize_t B_offset[4] = {0,0,0,0}; //offset to print new values of u at each timestep. Progressively gets larger in direction that output is extended.
    hsize_t B_dimsext[4] = {1,3,Nx,MyE-MyS}; //size of extended output each time it is extended. Should be the same as the original output i.e. equal to dimsf
    const size_t RANK = 4; //rank (dimension) of arrays to write to H5 file
    H5::DataSpace *dataspace_B = new DataSpace( RANK, B_dims, B_maxdims );

    H5::DSetCreatPropList B_prop;
    hsize_t B_chunk_dims[4] = {1,3,Nx,MyE-MyS};
    B_prop.setChunk(RANK, B_chunk_dims);

    H5::DataSet *B_set = new DataSet(file.createDataSet("B", H5::PredType::NATIVE_DOUBLE, *dataspace_B, B_prop));

    typedef VectorField::index_range range;
    VectorField::array_view<3>::type B_view = B[ boost::indices[range()][range(N_GC,Nx+N_GC)][range(N_GC,MyE-MyS+N_GC)] ];

    VectorField Phys_B(boost::extents[3][Nx][MyE-MyS]); //Cell-face average values of B across partial domain. Excludes ghost cells
    Phys_B = B_view; //copy view of B which removes ghost cells into smaller array which can be written into h5 file
    B_set->write( Phys_B.data(), PredType::NATIVE_DOUBLE ); //does not write ghost cells.

    //Write v_c to output H5 file and prepare to write data from additional time steps
    hsize_t vc_dims[4] = {1,3,Nx,MyE-MyS};
    hsize_t vc_maxdims[4] = {H5S_UNLIMITED,H5S_UNLIMITED,H5S_UNLIMITED,H5S_UNLIMITED};
    hsize_t vc_size[4] = {1,3,Nx,MyE-MyS}; //Total size of extended array. Increases with every time step. Initially should equal dimsf.
    hsize_t vc_offset[4] = {0,0,0,0}; //offset to print new values of u at each timestep. Progressively gets larger in direction that output is extended.
    hsize_t vc_dimsext[4] = {1,3,Nx,MyE-MyS}; //size of extended output each time it is extended. Should be the same as the original output i.e. equal to dimsf
    H5::DataSpace *dataspace_vc = new DataSpace( RANK, vc_dims, vc_maxdims );

    H5::DSetCreatPropList vc_prop;
    hsize_t vc_chunk_dims[4] = {1,3,Nx,MyE-MyS};
    vc_prop.setChunk(RANK, vc_chunk_dims);

    H5::DataSet *vc_set = new DataSet(file.createDataSet("vc", H5::PredType::NATIVE_DOUBLE, *dataspace_vc, vc_prop));

    VectorField::array_view<3>::type vc_view = vc[ boost::indices[range()][range(N_GC,Nx+N_GC)][range(N_GC,MyE-MyS+N_GC)] ];

    VectorField Phys_vc(boost::extents[3][Nx][MyE-MyS]); //Cell-face average values of vc across partial domain. Excludes ghost cells
    Phys_vc = vc_view; //copy view of vc which removes ghost cells into smaller array which can be written into h5 file
    vc_set->write( Phys_vc.data(), PredType::NATIVE_DOUBLE ); //does not write ghost cells.

    // Spatial coordinates (cell centers) for output H5 file
    hsize_t x_dims[2] = {1,Nx};
    hsize_t x_maxdims[2] = {1,Nx};
    const size_t RANK_coords = 2; //rank (dimension) of coordinate arrays to write to H5 file
    H5::DataSpace *dataspace_x = new DataSpace( RANK_coords, x_dims, x_maxdims );

    H5::DSetCreatPropList x_prop;
    hsize_t x_chunk_dims[2] = {1,Nx};
    x_prop.setChunk(RANK_coords, x_chunk_dims);

    H5::DataSet *x_set = new DataSet(file.createDataSet("x", H5::PredType::NATIVE_DOUBLE, *dataspace_x, x_prop));
    x_set->write( x.data(), PredType::NATIVE_DOUBLE );

    hsize_t y_dims[2] = {1,MyE-MyS};
    hsize_t y_maxdims[2] = {1,MyE-MyS};
    H5::DataSpace *dataspace_y = new DataSpace( RANK_coords, y_dims, y_maxdims );
    H5::DSetCreatPropList y_prop;
    hsize_t y_chunk_dims[2] = {1,MyE-MyS};
    y_prop.setChunk(RANK_coords, y_chunk_dims);

    H5::DataSet *y_set = new DataSet(file.createDataSet("y", H5::PredType::NATIVE_DOUBLE, *dataspace_y, y_prop));
    y_set->write( y.data(), PredType::NATIVE_DOUBLE );

    // Time data for output H5 file
    hsize_t t_dims[2] = {1,1};
    hsize_t t_maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
    hsize_t t_size[2] = {1,1};
    hsize_t t_offset[2] = {0,0};
    hsize_t t_dimsext[2] = {1,1};
    H5::DataSpace *dataspace_t = new DataSpace( RANK_coords, t_dims, t_maxdims );

    H5::DSetCreatPropList t_prop;
    hsize_t t_chunk_dims[2] = {1,1};
    t_prop.setChunk(RANK_coords, t_chunk_dims);

    H5::DataSet *t_set = new DataSet(file.createDataSet("t", H5::PredType::NATIVE_DOUBLE, *dataspace_t, t_prop));
    t_set->write( &t, PredType::NATIVE_DOUBLE );

    //Energy conservation data for output H5 file. Only include this in root process H5 file
    //Uses same dataspace and sizes as the time data, so don't need to define separate hsize_t and DataSpace for this
    H5::DataSet *U_B_set = new DataSet(file.createDataSet("U_B", H5::PredType::NATIVE_DOUBLE, *dataspace_t, t_prop));
    H5::DataSet *JH_set = new DataSet(file.createDataSet("JH", H5::PredType::NATIVE_DOUBLE, *dataspace_t, t_prop));
    H5::DataSet *PF_set = new DataSet(file.createDataSet("PF", H5::PredType::NATIVE_DOUBLE, *dataspace_t, t_prop));
    H5::DataSet *DeltaEInt_set = new DataSet(file.createDataSet("DeltaEInt", H5::PredType::NATIVE_DOUBLE, *dataspace_t, t_prop));

    // Add attributes to h5 file
    H5::DataSpace att_space(H5S_SCALAR);
    H5::Attribute att1 = file.createAttribute( "B_0", H5::PredType::NATIVE_DOUBLE, att_space );
    att1.write( H5::PredType::NATIVE_DOUBLE, &B_0 );
    H5::Attribute att2 = file.createAttribute( "t_0", H5::PredType::NATIVE_DOUBLE, att_space );
    att2.write( H5::PredType::NATIVE_DOUBLE, &t_0 );
    H5::Attribute att3 = file.createAttribute( "L_0", H5::PredType::NATIVE_DOUBLE, att_space );
    att3.write( H5::PredType::NATIVE_DOUBLE, &L_0 );
    H5::Attribute att4 = file.createAttribute( "Lx", H5::PredType::NATIVE_DOUBLE, att_space );
    att4.write( H5::PredType::NATIVE_DOUBLE, &Lx );
    H5::Attribute att5 = file.createAttribute( "Ly", H5::PredType::NATIVE_DOUBLE, att_space );
    att5.write( H5::PredType::NATIVE_DOUBLE, &Ly );

    /*
       Run simulation until t >= t_max
    */

    size_t iter = 0; //iteration counter
    const size_t save_cadence = save_cadenceCalc(Deltat,t_max,simparams.saves_number); //cadence for writing current state to H5 file
    const size_t ECons_cadence = simparams.ECons_cadence; //cadence for printing energy conservation information to terminal
    MPI_Barrier( MPI_COMM_WORLD );

    double start_t, end_t;
    if(world_rank == 0){
        start_t = MPI_Wtime();
    }
    double U_BSum0 = 0., U_BSum, JouleSum, PoyntingSum, DeltaEInt = 0., PoyntingFInt;
    std::vector<double> U_BSumVec(num_procs);
    std::vector<double> JouleSumVec(num_procs);
    std::vector<double> PoyntingSumVec(num_procs);

    std::vector<double> JouleSumHist;
    std::vector<double> PoyntingSumHist;
    std::vector<double> DeltaESumHist;

    //Compute_vc(B, vc, tC, N_GC, t, domain, bparams, world_rank);
    //exchng2Vector(vc, N_GC, comm1D, nbrleft, nbrright); //Exchange vc cells in a periodic manner in the y-direction into ghost cells
//    Compute_E(B, B, E, J, vc, tC, N_GC, t, domain, bparams);
//    E_BoundaryConditions(E, bparams, N_GC, comm1D, nbrleft, nbrright);
    EnergyConservation(B, E, J, tC.eta_O, N_GC, t, domain, bparams, U_BSum, JouleSum, PoyntingSum);

    MPI_Allgather(&U_BSum, 1, MPI_DOUBLE, U_BSumVec.data(), 1, MPI_DOUBLE, comm1D);
    MPI_Allgather(&JouleSum, 1, MPI_DOUBLE, JouleSumVec.data(), 1, MPI_DOUBLE, comm1D);
    MPI_Allgather(&PoyntingSum, 1, MPI_DOUBLE, PoyntingSumVec.data(), 1, MPI_DOUBLE, comm1D);

    if(world_rank == 0){
        U_BSum = std::reduce(U_BSumVec.begin(), U_BSumVec.end());
        JouleSum = std::reduce(JouleSumVec.begin(), JouleSumVec.end());
        PoyntingSum = std::reduce(PoyntingSumVec.begin(), PoyntingSumVec.end());

        U_BSum0 = U_BSum;
        JouleSumHist.push_back(JouleSum);
        PoyntingSumHist.push_back(PoyntingSum);
        DeltaESumHist.push_back(JouleSum+PoyntingSum);

        //Add energy conservation data to H5 file
        U_B_set->write( &U_BSum, PredType::NATIVE_DOUBLE );
        JH_set->write( &JouleSum, PredType::NATIVE_DOUBLE );
        PF_set->write( &PoyntingSum, PredType::NATIVE_DOUBLE );
        DeltaEInt_set->write( &DeltaEInt, PredType::NATIVE_DOUBLE );
    }

     ////////////////////////////////Print block/////////////////////////////////////////////////
//    std::fstream Dataout;
//    std::string timestr = std::to_string(t*t_0/yr);
//    if(iter%100 == 0){
//        Dataout.open("eta_HCheck_"+RANK_NAME+"_"+timestr+".dat",std::fstream::out);
//        for(size_t i = 0; i < B.shape()[1]; i++){
//            for(size_t j = 0; j < B.shape()[2]; j++){
//                Dataout << std::setprecision(10) << tC.eta_H[i][j] << "   ";
//            }
//            Dataout << std::endl;
//            Dataout << " " << std::endl;
//        }
//        Dataout.close();
//
//        Dataout.open("eta_OCheck_"+RANK_NAME+"_"+timestr+".dat",std::fstream::out);
//        for(size_t i = 0; i < B.shape()[1]; i++){
//            for(size_t j = 0; j < B.shape()[2]; j++){
//                Dataout << std::setprecision(10) << tC.eta_O[i][j] << "   ";
//            }
//            Dataout << std::endl;
//            Dataout << " " << std::endl;
//        }
//        Dataout.close();

//        Dataout.open("vcxCheck_"+RANK_NAME+"_"+timestr+".dat",std::fstream::out);
//        for(size_t i = 0; i < B.shape()[1]; i++){
//            for(size_t j = 0; j < B.shape()[2]; j++){
//                Dataout << std::setprecision(10) << vc[0][i][j] << "   ";
//            }
//            Dataout << std::endl;
//            Dataout << " " << std::endl;
//        }
//        Dataout.close();
//    }

    while(t <= t_max){

        iter++;

        RK_Step(B, E, J, vc, B_np1, tC, domain, process, bparams, t);

        //Update time and magnetic field
        t = t + Deltat;
        B = B_np1;

        // Check energy conservation. First recompute current density
        Compute_J(B, J, N_GC, domain);
        exchng2Vector(J, N_GC, comm1D, nbrleft, nbrright);
        Compute_vc(B, vc, J, tC, N_GC, t, domain, bparams, world_rank);
        exchng2Vector(vc, N_GC, comm1D, nbrleft, nbrright);
        EnergyConservation(B, E, J, tC.eta_O, N_GC, t, domain, bparams, U_BSum, JouleSum, PoyntingSum);

//        std::fstream Dataout;
//        std::string timestr = std::to_string(t*t_0/yr);
//        if(iter%100 == 0){
//            Dataout.open("BzCheck_"+RANK_NAME+"_"+timestr+".dat",std::fstream::out);
//            for(size_t i = 0; i < B.shape()[1]; i++){
//                for(size_t j = 0; j < B.shape()[2]; j++){
//                    Dataout << std::setprecision(10) << B[2][i][j] << "   ";
//                }
//                Dataout << std::endl;
//                Dataout << " " << std::endl;
//            }
//            Dataout.close();
//
//            Dataout.open("eta_OCheck_"+RANK_NAME+"_"+timestr+".dat",std::fstream::out);
//            for(size_t i = 0; i < B.shape()[1]; i++){
//                for(size_t j = 0; j < B.shape()[2]; j++){
//                    Dataout << std::setprecision(10) << tC.eta_O[i][j] << "   ";
//                }
//                Dataout << std::endl;
//                Dataout << " " << std::endl;
//            }
//            Dataout.close();
//        }

        // energy conservation and print it to terminal
        MPI_Allgather(&U_BSum, 1, MPI_DOUBLE, U_BSumVec.data(), 1, MPI_DOUBLE, comm1D);
        MPI_Allgather(&JouleSum, 1, MPI_DOUBLE, JouleSumVec.data(), 1, MPI_DOUBLE, comm1D);
        MPI_Allgather(&PoyntingSum, 1, MPI_DOUBLE, PoyntingSumVec.data(), 1, MPI_DOUBLE, comm1D);
        if(world_rank == 0){
            U_BSum = std::reduce(U_BSumVec.begin(), U_BSumVec.end());
            JouleSum = std::reduce(JouleSumVec.begin(), JouleSumVec.end());
            PoyntingSum = std::reduce(PoyntingSumVec.begin(), PoyntingSumVec.end());

            JouleSumHist.push_back(JouleSum);
            PoyntingSumHist.push_back(PoyntingSum);
            DeltaESumHist.push_back(JouleSum+PoyntingSum);

            DeltaEInt = TrapezoidIntegrator(DeltaESumHist,Deltat);
            PoyntingFInt = TrapezoidIntegrator(PoyntingSumHist,Deltat);
            if( (iter == 1) || (iter%ECons_cadence == 0) ){
                std::cout << "t = " << std::setprecision(5) << t*t_0/yr << " yr, ΔU_B = " << (U_BSum-U_BSum0)*pow(B_0*L_0,2.) << " erg/cm, ΔE = " << DeltaEInt*pow(B_0*L_0,2.) << " erg/cm";
                std::cout << ", PF Int = " << PoyntingFInt*pow(B_0*L_0,2.) << " erg/cm, % error = " << abs(DeltaEInt/(U_BSum-U_BSum0)-1.)*100. << std::endl;
            }
        }
        if(simparams.divBCheck == true){
            if( (iter == 1) || (iter%ECons_cadence == 0) ) divB_Monitor(B, N_GC, comm1D, world_rank, Ny_locs); //check that divB is still zero
        }
        if( iter % save_cadence == 0 ){
            // Add updated B values to H5 file
            B_size[0] += B_dimsext[0];
            B_size[1] = B_dims[1];
            B_size[2] = B_dims[2];
            B_size[3] = B_dims[3];
            B_set->extend(B_size);

            H5::DataSpace *B_filespace = new H5::DataSpace(B_set->getSpace());
            B_offset[0] += 1;
            B_filespace->selectHyperslab(H5S_SELECT_SET, B_dimsext, B_offset);
            H5::DataSpace *B_memspace = new H5::DataSpace( RANK, B_dimsext, NULL );

            B_view = B[ boost::indices[range()][range(N_GC,Nx+N_GC)][range(N_GC,MyE-MyS+N_GC)] ];
            Phys_B = B_view; //copy view of B which removes ghost cells into smaller array which can be written into h5 file
            B_set->write( Phys_B.data(), PredType::NATIVE_DOUBLE, *B_memspace, *B_filespace );
            std::cerr << "[HDF5] Wrote B block to file (iter=" << iter << ")" << std::endl;

            // Add updated v_c values to H5 file
            vc_size[0] += vc_dimsext[0];
            vc_size[1] = vc_dims[1];
            vc_size[2] = vc_dims[2];
            vc_size[3] = vc_dims[3];
            vc_set->extend(vc_size);

            H5::DataSpace *vc_filespace = new H5::DataSpace(vc_set->getSpace());
            vc_offset[0] += 1;
            vc_filespace->selectHyperslab(H5S_SELECT_SET, vc_dimsext, vc_offset);
            H5::DataSpace *vc_memspace = new H5::DataSpace( RANK, vc_dimsext, NULL );

            vc_view = vc[ boost::indices[range()][range(N_GC,Nx+N_GC)][range(N_GC,MyE-MyS+N_GC)] ];
            Phys_vc = vc_view; //copy view of B which removes ghost cells into smaller array which can be written into h5 file
            vc_set->write( Phys_vc.data(), PredType::NATIVE_DOUBLE, *vc_memspace, *vc_filespace );
                std::cerr << "[HDF5] Wrote vc block to file (iter=" << iter << ")" << std::endl;
                std::cerr << "[HDF5] Debug: completed vc write, about to prepare energy dataset hyperslabs" << std::endl;
                std::cerr.flush();
                // Debug: print current dataset dims before extend/select to diagnose mismatch
                try {
                    H5::DataSpace U_B_cur_space = U_B_set->getSpace();
                    hsize_t U_B_cur_dims[2] = {0,0};
                    U_B_cur_space.getSimpleExtentDims(U_B_cur_dims, NULL);
                    std::cerr << "[HDF5] U_B current dims = (" << U_B_cur_dims[0] << "," << U_B_cur_dims[1] << ")" << std::endl;

                    H5::DataSpace JH_cur_space = JH_set->getSpace();
                    hsize_t JH_cur_dims[2] = {0,0};
                    JH_cur_space.getSimpleExtentDims(JH_cur_dims, NULL);
                    std::cerr << "[HDF5] JH current dims = (" << JH_cur_dims[0] << "," << JH_cur_dims[1] << ")" << std::endl;

                    H5::DataSpace PF_cur_space = PF_set->getSpace();
                    hsize_t PF_cur_dims[2] = {0,0};
                    PF_cur_space.getSimpleExtentDims(PF_cur_dims, NULL);
                    std::cerr << "[HDF5] PF current dims = (" << PF_cur_dims[0] << "," << PF_cur_dims[1] << ")" << std::endl;

                    H5::DataSpace DE_cur_space = DeltaEInt_set->getSpace();
                    hsize_t DE_cur_dims[2] = {0,0};
                    DE_cur_space.getSimpleExtentDims(DE_cur_dims, NULL);
                    std::cerr << "[HDF5] DeltaEInt current dims = (" << DE_cur_dims[0] << "," << DE_cur_dims[1] << ")" << std::endl;
                } catch (const H5::Exception &e) {
                    std::cerr << "[HDF5] Exception querying current dataset spaces: " << e.getCDetailMsg() << std::endl;
                }
            delete B_filespace;
            delete B_memspace;
            delete vc_filespace;
            delete vc_memspace;

            // Add updated t value to H5 file
            t_size[0] += t_dimsext[0];
            t_size[1] = t_dims[0];
            t_set->extend(t_size);

            H5::DataSpace *t_filespace = new H5::DataSpace(t_set->getSpace());
            t_offset[0] += 1;
            t_filespace->selectHyperslab(H5S_SELECT_SET, t_dimsext, t_offset);
            // Use coordinate rank for time dataset memspace
            H5::DataSpace *t_memspace = new H5::DataSpace( RANK_coords, t_dimsext, NULL );
            t_set->write( &t, PredType::NATIVE_DOUBLE, *t_memspace, *t_filespace );
            delete t_filespace;
            delete t_memspace;

            //Add updated energy conservation values to H5 file for root process
            if(world_rank == 0){

                U_B_set->extend(t_size);
                JH_set->extend(t_size);
                PF_set->extend(t_size);
                DeltaEInt_set->extend(t_size);

                H5::DataSpace *U_B_filespace = new H5::DataSpace(U_B_set->getSpace());
                H5::DataSpace *JH_filespace = new H5::DataSpace(JH_set->getSpace());
                H5::DataSpace *PF_filespace = new H5::DataSpace(PF_set->getSpace());
                H5::DataSpace *DeltaEInt_filespace = new H5::DataSpace(DeltaEInt_set->getSpace());
                U_B_filespace->selectHyperslab(H5S_SELECT_SET, t_dimsext, t_offset);
                JH_filespace->selectHyperslab(H5S_SELECT_SET, t_dimsext, t_offset);
                PF_filespace->selectHyperslab(H5S_SELECT_SET, t_dimsext, t_offset);
                DeltaEInt_filespace->selectHyperslab(H5S_SELECT_SET, t_dimsext, t_offset);
                // Use coordinate rank for time/energy 2D datasets (rank = 2)
                H5::DataSpace *U_B_memspace = new H5::DataSpace( RANK_coords, t_dimsext, NULL );
                H5::DataSpace *JH_memspace = new H5::DataSpace( RANK_coords, t_dimsext, NULL );
                H5::DataSpace *PF_memspace = new H5::DataSpace( RANK_coords, t_dimsext, NULL );
                H5::DataSpace *DeltaEInt_memspace = new H5::DataSpace( RANK_coords, t_dimsext, NULL );
                std::cerr << "[HDF5] Writing energy datasets at iter=" << iter << std::endl;
                std::cerr << "[HDF5] t_size = (" << t_size[0] << "," << t_size[1] << "), t_offset = (" << t_offset[0] << "," << t_offset[1] << ")" << std::endl;
                std::cerr << "[HDF5] t_dimsext = (" << t_dimsext[0] << "," << t_dimsext[1] << ")" << std::endl;
                try{
                    std::cerr << "[HDF5] U_B_filespace selected points=" << U_B_filespace->getSelectNpoints() << ", U_B_memspace points=" << U_B_memspace->getSimpleExtentNpoints() << std::endl;
                    U_B_set->write( &U_BSum, PredType::NATIVE_DOUBLE, *U_B_memspace, *U_B_filespace );
                } catch(const H5::Exception &e){
                    std::cerr << "[HDF5] U_B_set write failed: " << e.getCDetailMsg() << std::endl;
                    throw;
                } catch(const std::exception &e){
                    std::cerr << "[HDF5] U_B_set write failed std::exception: " << e.what() << std::endl;
                    throw;
                }
                try{
                    std::cerr << "[HDF5] JH_filespace selected points=" << JH_filespace->getSelectNpoints() << ", JH_memspace points=" << JH_memspace->getSimpleExtentNpoints() << std::endl;
                    JH_set->write( &JouleSum, PredType::NATIVE_DOUBLE, *JH_memspace, *JH_filespace );
                } catch(const H5::Exception &e){
                    std::cerr << "[HDF5] JH_set write failed: " << e.getCDetailMsg() << std::endl;
                    throw;
                } catch(const std::exception &e){
                    std::cerr << "[HDF5] JH_set write failed std::exception: " << e.what() << std::endl;
                    throw;
                }
                try{
                    std::cerr << "[HDF5] PF_filespace selected points=" << PF_filespace->getSelectNpoints() << ", PF_memspace points=" << PF_memspace->getSimpleExtentNpoints() << std::endl;
                    PF_set->write( &PoyntingSum, PredType::NATIVE_DOUBLE, *PF_memspace, *PF_filespace );
                } catch(const H5::Exception &e){
                    std::cerr << "[HDF5] PF_set write failed: " << e.getCDetailMsg() << std::endl;
                    throw;
                } catch(const std::exception &e){
                    std::cerr << "[HDF5] PF_set write failed std::exception: " << e.what() << std::endl;
                    throw;
                }
                try{
                    std::cerr << "[HDF5] DeltaEInt_filespace selected points=" << DeltaEInt_filespace->getSelectNpoints() << ", DeltaEInt_memspace points=" << DeltaEInt_memspace->getSimpleExtentNpoints() << std::endl;
                    DeltaEInt_set->write( &DeltaEInt, PredType::NATIVE_DOUBLE, *DeltaEInt_memspace, *DeltaEInt_filespace );
                } catch(const H5::Exception &e){
                    std::cerr << "[HDF5] DeltaEInt_set write failed: " << e.getCDetailMsg() << std::endl;
                    throw;
                } catch(const std::exception &e){
                    std::cerr << "[HDF5] DeltaEInt_set write failed std::exception: " << e.what() << std::endl;
                    throw;
                }
                std::cerr << "[HDF5] Completed energy dataset writes for iter=" << iter << std::endl;
                delete U_B_filespace;
                delete U_B_memspace;
                delete JH_filespace;
                delete JH_memspace;
                delete PF_filespace;
                delete PF_memspace;
            }
        }
    }

    B_prop.close();
    vc_prop.close();
    x_prop.close();
    y_prop.close();
    t_prop.close();
    delete dataspace_B;
    delete dataspace_vc;
    delete dataspace_t;
    delete dataspace_x;
    delete dataspace_y;
    delete B_set;
    delete vc_set;
    delete x_set;
    delete y_set;
    delete t_set;
    delete U_B_set;
    delete JH_set;
    delete PF_set;
    delete DeltaEInt_set;
    file.close();

    /* Tear down the communication infrastructure */
    MPI_Finalize();

    if(world_rank == 0){
        end_t = MPI_Wtime();
        std::cout << "Solver done. Total run time is " << end_t-start_t << " s" << std::endl;
        std::cout << "Wrote simulation data to " << FILE_NAME << ".h5" << std::endl;
    }

    return 0;
}

/*
        Runge-Kutta timestep. Uses a two-step advance as described in Vigano et al Com. Phys. Commun. 183 (2012), 2042.
        Input: B, E, J, vc: magnetic field, electric field, current density and velocity field in reduced units
               tC: transCoeffs object containing Hall and Ohmic diffusivities and x-derivative of Hall diffusivity in reduced units
               dm, process: Domain and Process objects containing information about the simulation domain and the current process
               bparams: BandBCParams object containing information about the initial magnetic field and boundary conditions
               t: time in reduced units
        Output: B_np1: updated magnetic field
*/
void RK_Step(VectorField & B, VectorField & E, VectorField & J, VectorField & vc, VectorField & B_np1, TransCoeffs & tC, const Domain & dm, const Process & process, const BandBCParams & bparams, double t)
{
    size_t Ny = dm.Ny;
    size_t N_GC = dm.N_GC;
    std::vector<double> Deltax = dm.Deltax;
    double Deltay = dm.Deltay;
    double Deltat = dm.Deltat;
    std::vector<int> Ny_locs = dm.Ny_locs;
    std::vector<int> starts = dm.starts;

    int order = process.order;
    int world_rank = process.world_rank;
    int nbrleft = process.nbrleft;
    int nbrright = process.nbrright;
    MPI_Comm comm1D = process.comm1D;

    size_t i_f = round(bparams.x0_tor*dm.Nx);

    if(order == 2){

        static ScalarField Qx(boost::extents[B.shape()[1]][B.shape()[2]]); //RHS of Faraday's Law (EMF), x-component
        static ScalarField Qy(boost::extents[B.shape()[1]][B.shape()[2]]); //RHS of Faraday's Law (EMF), y-component
        static ScalarField Qz(boost::extents[B.shape()[1]][B.shape()[2]]); //RHS of Faraday's Law, z-component (toroidal field)
        static VectorField B_1(boost::extents[3][B.shape()[1]][B.shape()[2]]); //Updated B after first step

//        Compute_J(B, J, N_GC, dm); //this is already computed outside RKStep
//        exchng2Vector(J, N_GC, comm1D, nbrleft, nbrright); //this is already computed outside RKStep
        B_torEvolve(Qz, B, J, vc, tC, N_GC, t, dm, bparams);
        exchng2Scalar(Qz, N_GC, comm1D, nbrleft, nbrright);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_1[2][i][j] = B[2][i][j] + Deltat*Qz[i][j];
            }
        }
//        Compute_vc(B_1, vc, J, tC, N_GC, t, dm, bparams, world_rank);
//        exchng2Vector(vc, N_GC, comm1D, nbrleft, nbrright);
        Compute_E(B, B_1, E, J, vc, tC, N_GC, t, dm, bparams);
        E_BoundaryConditions(E, bparams, N_GC, comm1D, nbrleft, nbrright);
        Compute_EMF(Qx, Qy, E, N_GC, Deltax, Deltay);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_1[0][i][j] = B[0][i][j] + Deltat*Qx[i][j];
                B_1[1][i][j] = B[1][i][j] + Deltat*Qy[i][j];
            }
        }
        B_BoundaryConditions(B_1, bparams, Ny, N_GC, t, comm1D, world_rank, Ny_locs, starts, nbrleft, nbrright, dm);

//        std::fstream Dataout;
//        std::string timestr = std::to_string(t*t_0/yr);
////        if(iter%100 == 0){
//        Dataout.open("BzCheck_"+timestr+".dat",std::fstream::out);
//        for(size_t i = 0; i < B.shape()[1]; i++){
//            Dataout << i << std::endl;
//            for(size_t j = 0; j < B.shape()[2]; j++){
//                Dataout << std::setprecision(10) << B[2][i][j] << "   ";
//            }
//            Dataout << std::endl;
//            Dataout << " " << std::endl;
//        }
//        Dataout.close();

        Compute_J(B_1, J, N_GC, dm);
        exchng2Vector(J, N_GC, comm1D, nbrleft, nbrright);
        Compute_vc(B_1, vc, J, tC, N_GC, t, dm, bparams, world_rank);
        exchng2Vector(vc, N_GC, comm1D, nbrleft, nbrright);
        B_torEvolve(Qz, B_1, J, vc, tC, N_GC, t, dm, bparams);
        exchng2Scalar(Qz, N_GC, comm1D, nbrleft, nbrright);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_np1[2][i][j] = 0.5*( B[2][i][j] + B_1[2][i][j] + Deltat*Qz[i][j] );
            }
        }
//        Compute_vc(B_1, vc, J, tC, N_GC, t, dm, bparams, world_rank);
//        exchng2Vector(vc, N_GC, comm1D, nbrleft, nbrright);
        Compute_E(B_1, B_np1, E, J, vc, tC, N_GC, t, dm, bparams);
        E_BoundaryConditions(E, bparams, N_GC, comm1D, nbrleft, nbrright);
        Compute_EMF(Qx, Qy, E, N_GC, Deltax, Deltay);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_np1[0][i][j] = 0.5*( B[0][i][j] + B_1[0][i][j] + Deltat*Qx[i][j] );
                B_np1[1][i][j] = 0.5*( B[1][i][j] + B_1[1][i][j] + Deltat*Qy[i][j] );
            }
        }
        B_BoundaryConditions(B_np1, bparams, Ny, N_GC, t, comm1D, world_rank, Ny_locs, starts, nbrleft, nbrright, dm);
        exchng2Vector(B_np1, N_GC, comm1D, nbrleft, nbrright);

//        Dataout.open("Bz_np1Check_"+timestr+".dat",std::fstream::out);
//        for(size_t i = 0; i < B.shape()[1]; i++){
//            Dataout << i << std::endl;
//            for(size_t j = 0; j < B.shape()[2]; j++){
//                Dataout << std::setprecision(10) << B[2][i][j] << "   ";
//            }
//            Dataout << std::endl;
//            Dataout << " " << std::endl;
//        }
//        Dataout.close();
    }
    else if(order == 3){

        static ScalarField Qx(boost::extents[B.shape()[1]][B.shape()[2]]); //RHS of Faraday's Law (EMF), x-component
        static ScalarField Qy(boost::extents[B.shape()[1]][B.shape()[2]]); //RHS of Faraday's Law (EMF), y-component
        static ScalarField Qz(boost::extents[B.shape()[1]][B.shape()[2]]); //RHS of Faraday's Law, z-component (toroidal field)
        static VectorField B_1(boost::extents[3][B.shape()[1]][B.shape()[2]]); //Updated B after first step
        static VectorField B_2(boost::extents[3][B.shape()[1]][B.shape()[2]]); //Updated B after second step

//        Compute_J(B, J, N_GC, dm); //this is already computed outside RKStep
//        exchng2Vector(J, N_GC, comm1D, nbrleft, nbrright); //this is already computed outside RKStep
        Compute_vc(B, vc, J, tC, N_GC, t, dm, bparams, world_rank);
        exchng2Vector(vc, N_GC, comm1D, nbrleft, nbrright);
        B_torEvolve(Qz, B, J, vc, tC, N_GC, t, dm, bparams);
        exchng2Scalar(Qz, N_GC, comm1D, nbrleft, nbrright);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_1[2][i][j] = B[2][i][j] + Deltat*Qz[i][j];
            }
        }
        Compute_E(B, B_1, E, J, vc, tC, N_GC, t, dm, bparams);
        E_BoundaryConditions(E, bparams, N_GC, comm1D, nbrleft, nbrright);
        Compute_EMF(Qx, Qy, E, N_GC, Deltax, Deltay);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_1[0][i][j] = B[0][i][j] + Deltat*Qx[i][j];
                B_1[1][i][j] = B[1][i][j] + Deltat*Qy[i][j];
            }
        }
        B_BoundaryConditions(B_1, bparams, Ny, N_GC, t, comm1D, world_rank, Ny_locs, starts, nbrleft, nbrright, dm);

        Compute_J(B_1, J, N_GC, dm);
        exchng2Vector(J, N_GC, comm1D, nbrleft, nbrright);
        Compute_vc(B, vc, J, tC, N_GC, t, dm, bparams, world_rank);
        exchng2Vector(vc, N_GC, comm1D, nbrleft, nbrright);
        B_torEvolve(Qz, B_1, J, vc, tC, N_GC, t, dm, bparams);
        exchng2Scalar(Qz, N_GC, comm1D, nbrleft, nbrright);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_2[2][i][j] = 0.25*( 3.*B[2][i][j] + B_1[2][i][j] + Deltat*Qz[i][j] );
            }
        }
        Compute_E(B_1, B_2, E, J, vc, tC, N_GC, t, dm, bparams);
        E_BoundaryConditions(E, bparams, N_GC, comm1D, nbrleft, nbrright);
        Compute_EMF(Qx, Qy, E, N_GC, Deltax, Deltay);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_2[0][i][j] = 0.25*( 3.*B[0][i][j] + B_1[0][i][j] + Deltat*Qx[i][j] );
                B_2[1][i][j] = 0.25*( 3.*B[1][i][j] + B_1[1][i][j] + Deltat*Qy[i][j] );
            }
        }

        B_BoundaryConditions(B_2, bparams, Ny, N_GC, t, comm1D, world_rank, Ny_locs, starts, nbrleft, nbrright, dm);
        Compute_J(B_2, J, N_GC, dm);
        exchng2Vector(J, N_GC, comm1D, nbrleft, nbrright);
        Compute_vc(B, vc, J, tC, N_GC, t, dm, bparams, world_rank);
        exchng2Vector(vc, N_GC, comm1D, nbrleft, nbrright);
        B_torEvolve(Qz, B_2, J, vc, tC, N_GC, t, dm, bparams);
        exchng2Scalar(Qz, N_GC, comm1D, nbrleft, nbrright);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_np1[2][i][j] = 1./3.*( B[2][i][j] + 2.*B_2[2][i][j] + 2.*Deltat*Qz[i][j] );
            }
        }
        Compute_E(B_2, B_np1, E, J, vc, tC, N_GC, t, dm, bparams);
        E_BoundaryConditions(E, bparams, N_GC, comm1D, nbrleft, nbrright);
        Compute_EMF(Qx, Qy, E, N_GC, Deltax, Deltay);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_np1[0][i][j] = 1./3.*( B[0][i][j] + 2.*B_2[0][i][j] + 2.*Deltat*Qx[i][j] );
                B_np1[1][i][j] = 1./3.*( B[1][i][j] + 2.*B_2[1][i][j] + 2.*Deltat*Qy[i][j] );
            }
        }
        B_BoundaryConditions(B_np1, bparams, Ny, N_GC, t, comm1D, world_rank, Ny_locs, starts, nbrleft, nbrright, dm);
        exchng2Vector(B_np1, N_GC, comm1D, nbrleft, nbrright);

    }
    else{
        std::cout << "Invalid Runge-Kutta timestepping scheme order" << std::endl;
    }

    return ;

}

        ////////////////////////////////Print block/////////////////////////////////////////////////
//        std::fstream Dataout;
//        std::string timestr = std::to_string(t*t_0/yr);
//        if(iter%100 == 0){
//            Dataout.open("BzCheck_"+RANK_NAME+"_"+timestr+".dat",std::fstream::out);
//            for(size_t i = 0; i < B.shape()[1]; i++){
//                for(size_t j = 0; j < B.shape()[2]; j++){
//                    Dataout << std::setprecision(10) << B[2][i][j] << "   ";
//                }
//                Dataout << std::endl;
//                Dataout << " " << std::endl;
//            }
//            Dataout.close();
//
//            Dataout.open("ByCheck_"+RANK_NAME+"_"+timestr+".dat",std::fstream::out);
//            for(size_t i = 0; i < B.shape()[1]; i++){
//                for(size_t j = 0; j < B.shape()[2]; j++){
//                    Dataout << std::setprecision(10) << B[1][i][j] << "   ";
//                }
//                Dataout << std::endl;
//                Dataout << " " << std::endl;
//            }
//            Dataout.close();
//
//            Dataout.open("BxCheck_"+RANK_NAME+"_nprocs="+PROCS+"_"+timestr+".dat",std::fstream::out);
//            for(size_t i = 0; i < B.shape()[1]; i++){
//                for(size_t j = 0; j < B.shape()[2]; j++){
//                    Dataout << std::setprecision(10) << B[0][i][j] << "   ";
//                }
//                Dataout << std::endl;
//                Dataout << " " << std::endl;
//            }
//            Dataout.close();
//
//            Dataout.open("JyCheck"+RANK_NAME+"_"+timestr+".dat",std::fstream::out);
//            for(size_t i = 0; i < B.shape()[1]; i++){
//                for(size_t j = 0; j < B.shape()[2]; j++){
//                    Dataout << J[1][i][j] << "   ";
//                }
//                Dataout << std::endl;
//                Dataout << " " << std::endl;
//            }
//            Dataout.close();
//        }

        ////////////////////////////////Print block/////////////////////////////////////////////////
//    std::fstream Dataout;
//    Dataout.open("vczCheck_"+RANK_NAME+"_Initial.dat",std::fstream::out);
//    for(size_t i = 0; i < B.shape()[1]; i++){
//        for(size_t j = 0; j < B.shape()[2]; j++){
//            Dataout << std::setprecision(10) << vc[2][i][j] << "   ";
//        }
//        Dataout << std::endl;
//        Dataout << " " << std::endl;
//    }
//    Dataout.close();
////
//    Dataout.open("vcyCheck_"+RANK_NAME+"_Initial.dat",std::fstream::out);
//    for(size_t i = 0; i < B.shape()[1]; i++){
//        for(size_t j = 0; j < B.shape()[2]; j++){
//            Dataout << std::setprecision(10) << vc[1][i][j] << "   ";
//        }
//        Dataout << std::endl;
//        Dataout << " " << std::endl;
//    }
//    Dataout.close();
//
//    Dataout.open("JyCheck_"+RANK_NAME+"_Initial.dat",std::fstream::out);
//    for(size_t i = 0; i < J.shape()[1]; i++){
//        for(size_t j = 0; j < J.shape()[2]; j++){
//            Dataout << std::setprecision(10) << J[1][i][j] << "   ";
//        }
//        Dataout << std::endl;
//        Dataout << " " << std::endl;
//    }
//    Dataout.close();
//
//    Dataout.open("JxCheck_"+RANK_NAME+"_Initial.dat",std::fstream::out);
//    for(size_t i = 0; i < J.shape()[1]; i++){
//        for(size_t j = 0; j < J.shape()[2]; j++){
//            Dataout << std::setprecision(10) << J[0][i][j] << "   ";
//        }
//        Dataout << std::endl;
//        Dataout << " " << std::endl;
//    }
//    Dataout.close();
