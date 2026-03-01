#include <vector>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include "common.h"
#include "boundary_conditions.h"

/*
        Sets boundary condition on magnetic field
        Inputs: B: magnetic field as a vector field
               Bparams: BandBCParams object containing information about boundary conditions
               Ny: extent of combined domain in partitioned direction
               N_GC: number of ghost cells
               t: time in reduced units
               comm1D: MPI communicator for decomposed domain
               world_rank: rank of current process
               Ny_locs and starts: vectors containing the extent and starting indices of the domain in the decomposed direction
               nbrleft and nbrright: the ranks of the processes to the left and right of the current process
               dm: Domain object containing information about the simulation domain
        Output: B with updated boundary values
*/
void B_BoundaryConditions(VectorField & B, const BandBCParams & bparams, size_t Ny, size_t N_GC, double t, MPI_Comm comm1D, int world_rank, std::vector<int> & Ny_locs, std::vector<int> & starts, int nbrleft, int nbrright, const Domain & dm)
{
    //Exchange ghost cells in a periodic manner in the y-direction
    exchng2Vector(B, N_GC, comm1D, nbrleft, nbrright);

    std::vector<double> y = dm.y; //cell centers in the y-direction. Note: no ghost cells.

    std::vector<double> Bx_BC(B.shape()[2]-2*N_GC); //array to hold the Bx values to be Fourier decomposed and then used to compute By
    std::vector<double> By_BC(B.shape()[2]-2*N_GC); //array to hold the By values as computed from Bx according to the boundary condition

    double By_const = 0.; //constant term in By; computed by averaging
    double Byshear, Bzshear; //used for shearing field boundary condition

    //Fill Bx_BC with the values existing at the boundary x=Lx
    for(size_t j=0; j<B.shape()[2]-2*N_GC; j++){
        Bx_BC[j] = B[0][B.shape()[1]-N_GC-1][N_GC+j];
        By_const = By_const + 0.5*( B[1][B.shape()[1]-N_GC-2][N_GC+j] + B[1][B.shape()[1]-N_GC-1][N_GC+j] ); //take average of i=B.shape()[1]-N_GC-2,B.shape()[1]-N_GC-1 because want constant value of By at surface
    }
    By_const = By_const/double(B.shape()[2]-2*N_GC);

    By_BC_Calc(Bx_BC, By_BC, By_const, Ny, Ny_locs, starts, world_rank, comm1D);

    //Impose boundary conditions at non-periodic boundaries
    for(size_t i=0; i<N_GC; i++){
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

            /*
                Lower boundary x=0
            */
            if(bparams.B_perp_lower == "perfect_conducting"){
                B[0][N_GC-1-i][j] = -B[0][N_GC+i+1][j]; //Bx = 0
            }
            else if(bparams.B_perp_lower == "continuous"){
//              B[0][N_GC-1-i][j] = 2.*( Bparams.B_pol_init*sin(Bparams.theta_B) ) - B[0][N_GC+i+1][j]; //Bx = uniform background field
                B[0][N_GC-1-i][j] = B[0][N_GC+i+1][j]; //Bx is continuous
            }

            if(bparams.B_parallel_lower == "perfect_conducting"){
//              B[1][N_GC-1-i][j] = 2.*( Bparams.B_pol_init*cos(Bparams.theta_B) )  - B[1][N_GC+i][j]; //By = uniform background field (i.e., zero)
                B[1][N_GC-1-i][j] = B[1][N_GC+i][j]; //By is continuous
                B[2][N_GC-1-i][j] = B[2][N_GC+i][j]; //Bz is continuous
            }
            else if(bparams.B_parallel_lower == "zero"){
                B[1][N_GC-1-i][j] = -B[1][N_GC+i][j]; //By = 0
                B[2][N_GC-1-i][j] = -B[2][N_GC+i][j]; //Bz = 0
            }
            else if(bparams.B_parallel_lower == "shearing"){
                Bshear(B[1][N_GC+i][j],B[2][N_GC+i][j],y[j-N_GC],t,dm,bparams,Byshear,Bzshear);
                B[1][N_GC-1-i][j] = 2.*Byshear - B[1][N_GC+i][j]; //By is shearing value
                B[2][N_GC-1-i][j] = 2.*Bzshear - B[2][N_GC+i][j]; //Bz is shearing value
//                B[1][N_GC-1][j] = 1.5*Byshear - 0.25*(B[1][N_GC][j]+B[1][N_GC+1][j]); //By is shearing value
//                B[2][N_GC-1][j] = 1.5*Bzshear - 0.25*(B[2][N_GC][j]+B[2][N_GC+1][j]); //Bz is shearing value
//                B[1][N_GC-2][j] = 2.5*Byshear - 0.75*(B[1][N_GC][j]+B[1][N_GC+1][j]); //By is shearing value
//                B[2][N_GC-2][j] = 2.5*Bzshear - 0.75*(B[2][N_GC][j]+B[2][N_GC+1][j]); //Bz is shearing value
            }
            /*
                Upper boundary x=Lx
            */
            if(bparams.B_pol_upper == "vacuum"){
                B[0][B.shape()[1]-N_GC+i][j] = B[0][B.shape()[1]-N_GC-2-i][j]; //Bx is continuous
                //By equals the value related to the Fourier coefficients of Bx such that B is a potential field in the vacuum exterior
                B[1][B.shape()[1]-N_GC-1+i][j] = 2.*By_BC[j-N_GC] - B[1][B.shape()[1]-N_GC-2-i][j];
            }
            else if(bparams.B_pol_upper == "continuous"){
                B[0][B.shape()[1]-N_GC+i][j] = B[0][B.shape()[1]-N_GC-2-i][j]; //Bx is continuous
                B[1][B.shape()[1]-N_GC-1+i][j] = B[1][B.shape()[1]-N_GC-2-i][j]; //By is continuous
            }

            if(bparams.B_tor_upper == "vacuum"){
                B[2][B.shape()[1]-N_GC-1+i][j] = -B[2][B.shape()[1]-N_GC-2-i][j]; //Bz = 0
            }
            else if(bparams.B_tor_upper == "continuous"){
                B[2][B.shape()[1]-N_GC-1+i][j] = B[2][B.shape()[1]-N_GC-2-i][j]; //Bz is continuous
            }

        }
    }

    //Exchange ghost cells in a periodic manner in the y-direction after imposing boundary conditions
    exchng2Vector(B, N_GC, comm1D, nbrleft, nbrright);

    return;

}

/*
        Computes shearing Bz field at lower boundary for "shearing" B_parallel_lower boundary condition
        Inputs: By_init, Bz_init: initial values of By and Bz at boundary
                y: y coordinate in reduced units
                t: time in reduced units
                Bparams: BandBCParams object containing information about boundary conditions
                dm: Domain object containing information about the simulation domain
        Output: Byshear, Bzshear: shearing toroidal magnetic field at crust-core boundary in reduced units
*/
void Bshear(double By_init, double Bz_init, double y, double t, const Domain & dm, const BandBCParams & bparams, double & Byshear, double & Bzshear)
{

    double f_t; //time-dependent part of shearing boundary condition

    static double t_b = 0;//10.*yr/t_0;
    static double t_w = 100.*yr/t_0;
    static double y0 = 0.;
    static double y_w = dm.Ly;
    static double Byshear_max = bparams.B_parallel_lower_shear_By;
    static double Bzshear_max = bparams.B_parallel_lower_shear_Bz;

    if(bparams.B_parallel_shear_type == "uniform"){
        if( t >= t_b ){
            f_t = 1.;//tanh( (t-t_b)/t_w );
            Byshear = Byshear_max*f_t;
            Bzshear = Bzshear_max*f_t;
        }
        else{
            Byshear = -By_init;
            Bzshear = -Bz_init;
        }
    }
    else if(bparams.B_parallel_shear_type == "current_sheet"){
        if( t >= t_b ){
            f_t = 1.;//tanh( (t-t_b)/t_w );
            Byshear = 0.;
            Bzshear = f_t*Bzshear_max*sin(2.*pi*y/dm.Ly)*exp( -(y-y0)*(y-y0)/(2.*y_w*y_w) );
        }
        else{
            Byshear = -By_init;
            Bzshear = -Bz_init;
        }
    }

    return;
}

/*
        Computes By given Bx and the requirement that the external field is a potential field using fast Fourier transforms
        Scatters this information to each process
        Inputs: Bx_BC: vector of values of Bx on the x=Lx boundary
                By_const: mean of By at upper boundary for each process. Combined to give mean across entire domain.
                Ny: length of full simulation domain in y direction excluding ghost cells
                Ny_locs, starts: vectors containing the extent and starting indices of the domain in the decomposed direction
                world_rank: rank of current process
                comm1D: MPI communicator for decomposed domain
        Output: By_BC: vector of values of By on the x=Lx boundary
*/
void By_BC_Calc(std::vector<double> & Bx_BC, std::vector<double> & By_BC, double By_const, size_t Ny, std::vector<int> & Ny_locs, std::vector<int> & starts, int world_rank, MPI_Comm comm1D)
{
    std::vector<double> Bx;
    std::vector<double> By;
    std::vector<double> ByTemp;
    if(world_rank == 0){
        Bx.resize(Ny);
        By.resize(Ny);
        ByTemp.resize(Ny);
    }

    // Gather the values of Bx(x=Lx) from across the different processes
    // Don't take ghost cells
    MPI_Gatherv(&Bx_BC.front(), Ny_locs[world_rank], MPI_DOUBLE, &Bx.front(), Ny_locs.data(), starts.data(), MPI_DOUBLE, 0, comm1D);

    //Gather the values of average of By(x=Lx) from across the different processes
    // Don't take ghost cells
    static int num_procs = std::size(Ny_locs);
    std::vector<double> By_constVec(num_procs);
    MPI_Allgather(&By_const, 1, MPI_DOUBLE, By_constVec.data(), 1, MPI_DOUBLE, comm1D);

    /*
        On root process, perform FFT of Bx(x=Lx), then manipulate the coefficients
        and invert the FFT to obtain By(x=Lx).
    */
    if(world_rank == 0){
        int n = Bx.size();

        gsl_fft_real_wavetable * real;
        gsl_fft_halfcomplex_wavetable * hc;
        gsl_fft_real_workspace * work;

        work = gsl_fft_real_workspace_alloc (n);
        real = gsl_fft_real_wavetable_alloc (n);

        gsl_fft_real_transform (Bx.data(), 1, n, real, work);

        gsl_fft_real_wavetable_free (real);

        /*
            Set the Fourier coefficients of By(x=Lx) using those of Bx(x=Lx)
            Rule: Change cosine series coefficient of Bx to sine series coefficient of By
            Change sine series coefficient of Bx to -cosine series coefficient of By
            GSL real data FFT coefficient storage: constant term in first entry, then
            alternative cosine and sine series coefficients.
            Exception is constant term By[0]: compute this by taking average of
        */
        ByTemp[0] = std::reduce(By_constVec.begin(), By_constVec.end())*double(n)/double(num_procs);
        for(int j = 1; j<n/2; j++)
        {
            ByTemp[2*j-1] = Bx[2*j];
            ByTemp[2*j] = -Bx[2*j-1];
        }

        hc = gsl_fft_halfcomplex_wavetable_alloc (n);

        gsl_fft_halfcomplex_inverse (ByTemp.data(), 1, n, hc, work);

        gsl_fft_halfcomplex_wavetable_free (hc);
        gsl_fft_real_workspace_free (work);

        //Shift By values to left cell edge instead of cell center by averaging cell-centered values on either side of each edge
        //At leftmost cell with periodic boundary conditions, use value at j=0 and j=n-1.
        for(int j = 0; j<n; j++)
        {
            if( j == 0 ) By[j] = 0.5*(ByTemp[n-1] + ByTemp[0]);
            else By[j] = 0.5*(ByTemp[j-1] + ByTemp[j]);
        }

    }

    //Send the values of By(x=Lx) back to the processes.
    MPI_Scatterv(&By.front(), Ny_locs.data(), starts.data(), MPI_DOUBLE, &By_BC.front(), Ny_locs.data()[world_rank], MPI_DOUBLE, 0, comm1D);

    return;
}

/*
        Sets boundary values of the electric field

        Inputs: E: electric field in reduced units at cell edges
                Bparams: BandBCParams object containing information about boundary conditions
                N_GC: number of ghost cells
                comm1D: MPI communicator for the decomposed domain
                nbrleft, nbrright: rank of processes to the left and right of the current process
        Output: E: electric field with boundary values set
*/
void E_BoundaryConditions(VectorField & E, const BandBCParams & bparams, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright)
{

    //Exchange ghost cells in a periodic manner in the y-direction
    exchng2Vector(E, N_GC, comm1D, nbrleft, nbrright);

    //Impose boundary conditions at non-periodic boundaries
    for(size_t i=0; i<N_GC; i++){
        for(size_t j=N_GC; j<E.shape()[2]-N_GC; j++){
            /*
                Lower boundary x=0
            */
            if(bparams.E_perp_lower == "perfect_conducting"){
                E[0][N_GC-1-i][j] = E[0][N_GC+i][j]; //Ex is continuous
//                E[0][N_GC-1][j] = 0.5*(E[0][N_GC][j]+E[0][N_GC+1][j]); //Ex is continuous
//                E[0][N_GC-2][j] = 0.5*(E[0][N_GC+1][j]+E[0][N_GC+2][j]); //Ex is continuous
            }
            else if(bparams.E_perp_lower == "zero"){
                E[0][N_GC-1-i][j] = -E[0][N_GC+i][j]; //Ex is continuous
            }
            if(bparams.E_parallel_lower == "perfect_conducting"){
                E[1][N_GC-1-i][j] = -E[1][N_GC+1+i][j]; //Ey = 0
                E[2][N_GC-1-i][j] = -E[2][N_GC+1+i][j]; //Ez = 0
//                E[1][N_GC][j] = 0.5*E[1][N_GC+1][j]; //to reduce checkerboard oscillation
//                E[2][N_GC][j] = 0.5*E[1][N_GC+1][j]; //to reduce checkerboard oscillation
                E[1][N_GC][j] = 0.;
                E[2][N_GC][j] = 0.;
            }
            else if(bparams.E_parallel_lower == "continuous"){
                E[1][N_GC-1-i][j] = E[1][N_GC+1+i][j]; //Ey is continuous
                E[2][N_GC-1-i][j] = E[2][N_GC+1+i][j]; //Ez is continuous
//                E[1][N_GC][j] = 0.5*E[1][N_GC+1][j]; //to reduced checkerboard oscillation
//                E[2][N_GC][j] = 0.5*E[2][N_GC+1][j]; //to reduced checkerboard oscillation
            }
            /*
                Upper boundary x=Lx
            */
            if(bparams.E_perp_upper == "vacuum"){
                E[0][E.shape()[1]-N_GC-1+i][j] = -E[0][E.shape()[1]-N_GC-2-i][j]; //Ex vanishes
            }
            else if(bparams.E_perp_upper == "continous"){
                E[0][E.shape()[1]-N_GC-1+i][j] = E[0][E.shape()[1]-N_GC-2-i][j]; //Ex is continuous
            }
            if(bparams.E_parallel_upper == "vacuum"){
                E[1][E.shape()[1]-N_GC+i][j] = -E[1][E.shape()[1]-N_GC-2-i][j]; //Ey = 0
                E[2][E.shape()[1]-N_GC+i][j] = -E[2][E.shape()[1]-N_GC-2-i][j]; //Ez = 0
            }
            else if(bparams.E_parallel_upper == "continuous"){
                E[1][E.shape()[1]-N_GC+i][j] = E[1][E.shape()[1]-N_GC-2-i][j]; //Ey is continuous
                E[2][E.shape()[1]-N_GC+i][j] = E[2][E.shape()[1]-N_GC-2-i][j]; //Ez is continuous
            }
        }
    }

    return;

}
