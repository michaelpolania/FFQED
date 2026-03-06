#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>
#include <string>
#include <mpi.h>
#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

#include "common.h"

constexpr double pi = 3.141592653589793238463;
constexpr double c = 29979245800; //speed of light in cm/s
constexpr double hbarc = 197.3269804; //hbar times c in MeV*fm
constexpr double unit_e = 4.80320425e-10; //elementary charge in units of statcoulomb
constexpr double yr = 3600*24*365; //1 year in seconds
constexpr double G = 6.67430e-8; //gravitational constant in dyn*cm^2/g^2
constexpr double k_B = 8.61733034e-11; //Boltzmann constant in MeV/K
constexpr double k_Bcgs = 1.380649e-16; //Boltzmann constant in cgs units (erg/K)
constexpr double M_e = 0.51099895000; //electron mass in MeV
constexpr double M_m = 105.66; //muon mass in MeV (105.658375 MeV)
constexpr double M_N = 938.92; //mean nucleon mass in MeV (938.91875434 MeV)
constexpr double M_n = 939.56542052; //neutron mass in MeV
constexpr double MeVtoErg = 1/6.2415e5; //conversion factor from MeV to erg
constexpr double M_solar = 1.98847e33; //solar mass in g
constexpr double alpha_e = 0.0072973525693; //electromagnetic fine structure constant (dimensionless)
constexpr double eB_crit = M_e*M_e; //critical magnetic field times elementary charge in MeV^2
constexpr double ComptonWL = hbarc/M_e; //reduced electron Compton wavelength hbar*c/(m_e*c^2) in fm
constexpr double amu = 931.49410242; //1 amu in MeV
constexpr double microU = 931.49410242*1.0e-6; //1 micro-u (10^{-6} atomic mass units) in MeV
constexpr double gammaEM = 0.577215664901532; //Euler-Mascheroni constant
constexpr double sigma_SB = 5.670374e-5; //Stefan-Boltzmann constant in erg/cm^2/s/K^4

constexpr double GaussConverter = 6.241509074e-8/hbarc; //conversion factor between 1 statCoulomb*Gauss to __ MeV/fm/hbarc = ___ fm^-2
constexpr double GaussConverter2 = 6.241509074e-8*hbarc; //conversion factor between 1 statCoulomb*Gauss to __ MeV/fm*hbarc = ___ MeV^2
constexpr double GaussConverter3 = 4.002719868e16; //conversion factor: 1 sqrt(MeV/fm^3) = 4.002719868e16 G

constexpr double B_0 = 1e13; //characteristic magnetic field (G)
constexpr double n_e0 = 1e-4; //characteristic electron density (fm^{-3})
constexpr double L_0 = 1e5; //characteristic length scale (cm)
constexpr double t_0 = 4.*pi*unit_e*n_e0*1e39*L_0*L_0/(B_0*c); //characteristic timescale (s). Taken as the Hall time-scale with length scale 1 km and field 10^{13} G
constexpr double T_0 = 1e8; //characteristic temperature (K)
constexpr double s_0 = 1e18; //characteristic entropy density (erg/K/cm^3)
constexpr double E_0 = B_0*L_0/(c*t_0); //characteristic electric field (statV/cm)

/*
    Load simulation parameters from SimSetup.in
    Inputs: world_rank: rank of current process
    Output: params: object of class SimParams
            Bparams: object of class BandBCParams
*/
void load_params(SimParams & params, BandBCParams & Bparams, int world_rank)
{
    std::ifstream file("SimSetup.in");
    std::string line, entry;

    std::vector<std::string> labels;
    std::vector<std::string> values;

    while(std::getline(file, line)){
        if( line.rfind("#",0) != 0 && line.size() != 0 ){ //skips comments starting with "#" and blank lines
            std::vector<std::string> current_entries;
            std::string::iterator end_pos = std::remove(line.begin(), line.end(), ' ');
            line.erase(end_pos, line.end()); //removes white space from "line"
            std::stringstream current_line( line );
            while(std::getline(current_line, entry, ':')){
                current_entries.push_back( entry );
            }
            labels.push_back( current_entries[0] ); //label to left of colon
            values.push_back( current_entries[1] ); //value to right of colon
        }
    }
    file.close();

    //Set values of parameters in params object
    for(size_t j=0; j<labels.size(); j++){

        if(labels[j] == "CrustEOS")
            params.CrustEOS  = values[j];
        if(labels[j] == "RK_order")
            params.RK_order = stoi(values[j]);
        if(labels[j] == "varying_mesh"){
            if(values[j] == "true")
                params.varying_mesh = true;
            else params.varying_mesh = false;
        }
        if(labels[j] == "Nx")
            params.Nx = stoul(values[j]);
        if(labels[j] == "Ny")
            params.Ny = stoul(values[j]);
        if(labels[j] == "x_min")
            params.x_min = stod(values[j])/L_0;
        if(labels[j] == "x_max")
            params.x_max = stod(values[j])/L_0;
        if(labels[j] == "y_min")
            params.y_min = stod(values[j])/L_0;
        if(labels[j] == "y_max")
            params.y_max = stod(values[j])/L_0;
        if(labels[j] == "t_max")
            params.t_max = stod(values[j])*yr/t_0;
        if(labels[j] == "k_C")
            params.k_C = stod(values[j]);
        if(labels[j] == "rho_cutoff")
            params.rho_cutoff = stod(values[j]);
        if(labels[j] == "temperature")
            params.temperature = stod(values[j]);
        if(labels[j] == "saves_number")
            params.saves_number = stoul(values[j]);
        if(labels[j] == "ECons_cadence")
            params.ECons_cadence = stoul(values[j]);
        if(labels[j] == "divBCheck"){
            if(values[j] == "true")
                params.divBCheck = true;
            else params.divBCheck = false;
        }
        if(labels[j] == "OutputFile")
            params.OutputFile = values[j];

        // ---- B field initial conditions ----
        if(labels[j] == "B_pol_init")
            Bparams.B_pol_init = stod(values[j])/B_0;
        if(labels[j] == "theta_B")
            Bparams.theta_B = stod(values[j])*pi/180.;
        if(labels[j] == "B_tor_init"){
            if(values[j] == "true")
                Bparams.B_tor_init = true;
            else Bparams.B_tor_init = false;
        }
        if(labels[j] == "B_tor_max")
            Bparams.B_tor_max = stod(values[j])/B_0;
        if(labels[j] == "B_tor_x_center")
            Bparams.B_tor_x_center = stod(values[j])/L_0;
        if(labels[j] == "B_tor_x_width")
            Bparams.B_tor_x_width = stod(values[j])/L_0;
        if(labels[j] == "B_tor_y_center")
            Bparams.B_tor_y_center = stod(values[j])/L_0;
        if(labels[j] == "B_tor_y_width")
            Bparams.B_tor_y_width = stod(values[j])/L_0;

        // ---- D field initial conditions ---- (ADDED)
        if(labels[j] == "D_pol_init")
            Bparams.D_pol_init = stod(values[j])/B_0;
        if(labels[j] == "theta_D")
            Bparams.theta_D = stod(values[j])*pi/180.;
        if(labels[j] == "D_tor_init"){
            if(values[j] == "true")
                Bparams.D_tor_init = true;
            else Bparams.D_tor_init = false;
        }
        if(labels[j] == "D_tor_max")
            Bparams.D_tor_max = stod(values[j])/B_0;
        if(labels[j] == "D_tor_x_center")
            Bparams.D_tor_x_center = stod(values[j])/L_0;
        if(labels[j] == "D_tor_x_width")
            Bparams.D_tor_x_width = stod(values[j])/L_0;
        if(labels[j] == "D_tor_y_center")
            Bparams.D_tor_y_center = stod(values[j])/L_0;
        if(labels[j] == "D_tor_y_width")
            Bparams.D_tor_y_width = stod(values[j])/L_0;

        // ---- B field boundary conditions ----
        if(labels[j] == "B_perp_lower")
            Bparams.B_perp_lower = values[j];
        if(labels[j] == "B_parallel_lower")
            Bparams.B_parallel_lower = values[j];
        if(labels[j] == "B_parallel_shear_type")
            Bparams.B_parallel_shear_type = values[j];
        if(labels[j] == "B_parallel_lower_shear_By")
            Bparams.B_parallel_lower_shear_By = stod(values[j])/B_0;
        if(labels[j] == "B_parallel_lower_shear_Bz")
            Bparams.B_parallel_lower_shear_Bz = stod(values[j])/B_0;
        if(labels[j] == "B_pol_upper")
            Bparams.B_pol_upper = values[j];
        if(labels[j] == "B_tor_upper")
            Bparams.B_tor_upper = values[j];

        // ---- D field boundary conditions ---- (ADDED)
        if(labels[j] == "D_perp_lower")
            Bparams.D_perp_lower = values[j];
        if(labels[j] == "D_parallel_lower")
            Bparams.D_parallel_lower = values[j];
        if(labels[j] == "D_parallel_shear_type")
            Bparams.D_parallel_shear_type = values[j];
        if(labels[j] == "D_parallel_lower_shear_Dy")
            Bparams.D_parallel_lower_shear_Dy = stod(values[j])/B_0;
        if(labels[j] == "D_parallel_lower_shear_Dz")
            Bparams.D_parallel_lower_shear_Dz = stod(values[j])/B_0;
        if(labels[j] == "D_pol_upper")
            Bparams.D_pol_upper = values[j];
        if(labels[j] == "D_tor_upper")
            Bparams.D_tor_upper = values[j];

        // ---- E field boundary conditions ----
        if(labels[j] == "E_perp_lower")
            Bparams.E_perp_lower = values[j];
        if(labels[j] == "E_parallel_lower"){
            Bparams.E_parallel_lower = values[j];
            if(Bparams.B_parallel_lower == "shearing" && Bparams.E_parallel_lower != "continuous"){
                std::cout << "WARNING: inconsistent use of B_parallel_lower and E_parallel_lower boundary conditions" << std::endl;
            }
        }
        if(labels[j] == "E_perp_upper")
            Bparams.E_perp_upper = values[j];
        if(labels[j] == "E_parallel_upper")
            Bparams.E_parallel_upper = values[j];

        // ---- H field boundary conditions ---- (ADDED)
        if(labels[j] == "H_perp_lower")
            Bparams.H_perp_lower = values[j];
        if(labels[j] == "H_parallel_lower")
            Bparams.H_parallel_lower = values[j];
        if(labels[j] == "H_perp_upper")
            Bparams.H_perp_upper = values[j];
        if(labels[j] == "H_parallel_upper")
            Bparams.H_parallel_upper = values[j];

        // ---- Toroidal velocity shear ----
        if(labels[j] == "tor_vel_shear")
            Bparams.tor_vel_shear = values[j];
        if(labels[j] == "t_w")
            Bparams.t_w = stod(values[j])*yr/t_0;
        if(labels[j] == "t_m")
            Bparams.t_m = stod(values[j])*yr/t_0;
        if(labels[j] == "tor_n")
            Bparams.tor_n = stod(values[j]);
        if(labels[j] == "vc_mag")
            Bparams.vc_mag = stod(values[j]);
        if(labels[j] == "xwidth_tor")
            Bparams.xwidth_tor = stod(values[j])/L_0;
        if(labels[j] == "ywidth_tor")
            Bparams.ywidth_tor = stod(values[j])/L_0;
        if(labels[j] == "x0_tor")
            Bparams.x0_tor = stod(values[j])/L_0;
        if(labels[j] == "y0_tor")
            Bparams.y0_tor = stod(values[j])/L_0;
    }

    // ---- Validation block ----
    if( world_rank == 0 ){
        if( params.CrustEOS.empty() )
            std::cout << "Missing crust EOS data table name" << std::endl;
        if( params.RK_order == 0 )
            std::cout << "Missing Runge-Kutta timestepping order" << std::endl;
        if( params.Nx == 0 )
            std::cout << "Missing x-direction resolution" << std::endl;
        if( params.Ny == 0)
            std::cout << "Missing y-direction resolution" << std::endl;
        if( params.x_min < 0. || params.x_max < 0. || params.y_min > 0. || params.y_max < 0. )
            std::cout << "Missing x_min, x_max, y_min or y_max" << std::endl;
        // CHANGED: t_max = 0 is valid for initial condition verification runs, so no warning needed
        if( params.k_C < 1e-20 )
            std::cout << "Missing Courant number" << std::endl;
        if( params.rho_cutoff < 1e-20 )
            std::cout << "Missing cutoff density" << std::endl;
        if( params.temperature < 1e-20 )
            std::cout << "Missing temperature" << std::endl;
        if( params.saves_number == 0 )
            std::cout << "Missing maximum number of snapshots to save" << std::endl;
        if( params.ECons_cadence == 0 )
            std::cout << "Missing energy conservation print-out cadence" << std::endl;
        if( params.OutputFile.empty() )
            std::cout << "Missing output file name" << std::endl;

        // ---- B field BC validation ----
        if( Bparams.B_perp_lower.empty() ){
            Bparams.B_perp_lower = "perfect_conducting";
            std::cout << "Missing lower BC on perpendicular magnetic field; assuming perfect conducting" << std::endl;
        }
        if( Bparams.B_parallel_lower.empty() ){
            Bparams.B_parallel_lower = "perfect_conducting";
            std::cout << "Missing lower BC on parallel magnetic field; assuming perfect conducting" << std::endl;
        }
        if( Bparams.B_pol_upper.empty() ){
            Bparams.B_pol_upper = "vacuum";
            std::cout << "Missing upper BC on poloidal magnetic field; assuming vacuum" << std::endl;
        }
        if( Bparams.B_tor_upper.empty() ){
            Bparams.B_tor_upper = "vacuum";
            std::cout << "Missing upper BC on toroidal magnetic field; assuming vacuum" << std::endl;
        }

        // ---- D field BC validation ---- (ADDED)
        if( Bparams.D_perp_lower.empty() ){
            Bparams.D_perp_lower = "perfect_conducting";
            std::cout << "Missing lower BC on perpendicular displacement field; assuming perfect conducting" << std::endl;
        }
        if( Bparams.D_parallel_lower.empty() ){
            Bparams.D_parallel_lower = "perfect_conducting";
            std::cout << "Missing lower BC on parallel displacement field; assuming perfect conducting" << std::endl;
        }
        if( Bparams.D_pol_upper.empty() ){
            Bparams.D_pol_upper = "vacuum";
            std::cout << "Missing upper BC on poloidal displacement field; assuming vacuum" << std::endl;
        }
        if( Bparams.D_tor_upper.empty() ){
            Bparams.D_tor_upper = "vacuum";
            std::cout << "Missing upper BC on toroidal displacement field; assuming vacuum" << std::endl;
        }

        // ---- E field BC validation ----
        if( Bparams.E_perp_lower.empty() ){
            Bparams.E_perp_lower = "perfect_conducting";
            std::cout << "Missing lower BC on perpendicular electric field; assuming perfect conducting" << std::endl;
        }
        if( Bparams.E_parallel_lower.empty() ){
            Bparams.E_parallel_lower = "perfect_conducting";
            std::cout << "Missing lower BC on parallel electric field; assuming perfect conducting" << std::endl;
        }
        if( Bparams.E_perp_upper.empty() ){
            Bparams.E_perp_upper = "vacuum";
            std::cout << "Missing upper BC on perpendicular electric field; assuming vacuum" << std::endl;
        }
        if( Bparams.E_parallel_upper.empty() ){
            Bparams.E_parallel_upper = "vacuum";
            std::cout << "Missing upper BC on parallel electric field; assuming vacuum" << std::endl;
        }

        // ---- H field BC validation ---- (ADDED)
        if( Bparams.H_perp_lower.empty() ){
            Bparams.H_perp_lower = "perfect_conducting";
            std::cout << "Missing lower BC on perpendicular H field; assuming perfect conducting" << std::endl;
        }
        if( Bparams.H_parallel_lower.empty() ){
            Bparams.H_parallel_lower = "perfect_conducting";
            std::cout << "Missing lower BC on parallel H field; assuming perfect conducting" << std::endl;
        }
        if( Bparams.H_perp_upper.empty() ){
            Bparams.H_perp_upper = "vacuum";
            std::cout << "Missing upper BC on perpendicular H field; assuming vacuum" << std::endl;
        }
        if( Bparams.H_parallel_upper.empty() ){
            Bparams.H_parallel_upper = "vacuum";
            std::cout << "Missing upper BC on parallel H field; assuming vacuum" << std::endl;
        }

        // ---- Toroidal velocity shear validation ----
        if( Bparams.tor_vel_shear != "none" ){
            if( Bparams.t_w < 1e-20 )
                std::cout << "Missing time width of toroidal velocity shear" << std::endl;
            if( Bparams.t_m < 1e-20 )
                std::cout << "Missing centered time of toroidal velocity shear" << std::endl;
            if( Bparams.tor_n < 1e-20 )
                std::cout << "Missing steepness of toroidal velocity shear" << std::endl;
            if( Bparams.vc_mag < 1e-30 )
                std::cout << "Missing magnitude of toroidal velocity shear" << std::endl;
            if( Bparams.xwidth_tor < 1e-20 )
                std::cout << "Missing x-width of toroidal velocity shear" << std::endl;
            if( Bparams.ywidth_tor < 1e-20 )
                std::cout << "Missing y-width of toroidal velocity shear" << std::endl;
        }
    }

    return;
}

/*
    Partitions initial domain into partial domains to be evolved by each process, in 1 dimension
    Input: N: size of domain in direction to be partitioned
           num_procs: number of processes
           MyID: rank of this process
    Output: s, e: start and end indices of partition of array to be evolved by this process
*/
void MPE_Decomp1D( size_t N, int num_procs, int MyID, size_t &s, size_t &e )
{
    size_t nlocal, deficit;

    nlocal = N / num_procs;
    s = MyID * nlocal;
    deficit = N % num_procs;
    s = s + std::min(MyID,(int)deficit);

    if (MyID < (int)deficit){
      nlocal = nlocal;
    }
    e = s + nlocal;
    if( (e > N) || (MyID == num_procs-1)){
     e = N;
    }

    return;
}

/*
    Fills ghost cells by exchanging data between cells - for VectorField objects
    Input: A: local array A
           N_GC: number of ghost cells to exchange
           comm1D: MPI communicator for processes organized into Cartesian grid
           nbrleft, nbrright: which processes are to the left and right of this one
*/
void exchng2Vector(VectorField & A, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright)
{
    size_t N_comps = A.shape()[0];
    size_t Nx = A.shape()[1];
    size_t Ny_loc = A.shape()[2];

    MPI_Datatype stridetype;
    MPI_Type_vector( N_comps*Nx, N_GC, Ny_loc, MPI_DOUBLE, &stridetype);
    MPI_Type_commit( &stridetype );

    MPI_Sendrecv( &(A[0][0][Ny_loc-2*N_GC]), 1, stridetype, nbrright, 0, &(A[0][0][0]), 1, stridetype, nbrleft, 0, comm1D, MPI_STATUS_IGNORE);
    MPI_Sendrecv( &(A[0][0][N_GC]), 1, stridetype, nbrleft, 1, &(A[0][0][Ny_loc-N_GC]), 1, stridetype, nbrright, 1, comm1D, MPI_STATUS_IGNORE);

    MPI_Type_free( &stridetype );

    return;
}

/*
    Fills ghost cells by exchanging data between cells - for ScalarField objects
    Input: A: local array A
           N_GC: number of ghost cells to exchange
           comm1D: MPI communicator for processes organized into Cartesian grid
           nbrleft, nbrright: which processes are to the left and right of this one
*/
void exchng2Scalar(ScalarField & A, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright)
{
    size_t Nx = A.shape()[0];
    size_t Ny_loc = A.shape()[1];

    MPI_Datatype stridetype;
    MPI_Type_vector( Nx, N_GC, Ny_loc, MPI_DOUBLE, &stridetype);
    MPI_Type_commit( &stridetype );

    MPI_Sendrecv( &(A[0][Ny_loc-2*N_GC]), 1, stridetype, nbrright, 0, &(A[0][0]), 1, stridetype, nbrleft, 0, comm1D, MPI_STATUS_IGNORE);
    MPI_Sendrecv( &(A[0][N_GC]), 1, stridetype, nbrleft, 1, &(A[0][Ny_loc-N_GC]), 1, stridetype, nbrright, 1, comm1D, MPI_STATUS_IGNORE);

    MPI_Type_free( &stridetype );

    return;
}

/*
    Computes cadence for saving simulation snapshots to H5 file
    Inputs: Deltat, t_max: timestep and maximum simulation time in reduced units
            saves_number: maximum number of snapshots saved to H5 file
    Output: cadence for saving simulation snapshots to H5 file
*/
size_t save_cadenceCalc(double Deltat, double t_max, size_t saves_number){

    size_t save_cadence = 1;
    size_t N_t = std::ceil( t_max/Deltat );

    if( N_t < saves_number ){
        save_cadence = 1;
    }
    else{
        save_cadence = std::ceil( t_max/Deltat/double(saves_number) );
    }

    return save_cadence;
}

/*
    Minmod function
    Arguments: a, b, c: values to compare
    Output: minmod(a,b,c)
*/
double minmod(double a, double b, double c)
{
    double result = 0.;

    if( a>0 && b>0 && c>0 ){
        result = std::min( a, std::min(b,c) );
    }
    else if( a<0 && b<0 && c<0 ){
        result = std::max( a, std::max(b,c) );
    }

    return result;
}

/*
    Trapezoid rule integration in 1D
    Input: A: vector of integrand sampled in 1D every dx
           dx: spacing between every sample of integrand
    Output: integral
*/
double TrapezoidIntegrator(std::vector<double> & A, double dx)
{
    double integral = 0.;
    for(size_t j=0; j<A.size()-1; j++){
        integral = integral + 0.5*(A[j+1]+A[j])*dx;
    }

    return integral;
}
