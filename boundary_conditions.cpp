#include <vector>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include "common.h"
#include "boundary_conditions.h"

/*
MP edits below
*/

/*
Initial velocity profile (only along z-direction).
v_z(y) = - v_0 * exp((-(y-y_0))/(2 * sigma_y ** 2)) * sin((2 * pi * f)/(T)) 
*/

double Initialvz(double y, double t, void *driver)
{
    vConfig_params *p = (vConfig_params *) driver;

    double v_max    = p->v_max;
    double y_center = p->y_center;
    double y_width  = p->y_width;
    double f        = p->f;
    
    return -v_max * exp(-pow(y - y_center, 2.0) / (2.0 * y_width * y_width)) * sin(2.0 * M_PI * f * t);   
}



/*
        Sets lower boundary condition on Displacement field
        Inputs: D: magnetic field as a vector field
               bparams: BandBCParams object containing information about boundary conditions
               Ny: extent of combined domain in partitioned direction
               N_GC: number of ghost cells
               t: time in reduced units
               comm1D: MPI communicator for decomposed domain
               world_rank: rank of current process
               Ny_locs and starts: vectors containing the extent and starting indices of the domain in the decomposed direction
               nbrleft and nbrright: the ranks of the processes to the left and right of the current process
               dm: Domain object containing information about the simulation domain
        Output: D with updated boundary values
*/

void LowerBoundary_D(VectorField& D, const VectorField& V, const VectorField& B, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright, double t, vConfig_params& driver, double y_min, double dy)
{
    //Exchange ghost cells
    //exchng2Vector(D, N_GC, comm1D, nbrleft, nbrright);

    //Loops all the way up to last physical cell in the y-direction
    for(size_t j = N_GC; j < D.shape()[2] - N_GC; j++){

        //Compute y coordinate of this cell
        double y_j = y_min + (j - N_GC) * dy;

        //Computes boundary velocity for coordinate y_j
        double Vz = Initialvz(y_j, t, &driver);

        // B components at lower boundary... call B_Boundary condition function?

        // D_BC = -V x B (where v_z is the nonzero component)
        double DBC_x =  Vz * By;   
        double DBC_y = -Vz * B[0][N_GC][j];   
        double DBC_z =  0.0;        

        //D_x already lives on the boundary
        D[0][N_GC][j] = DBC_x;
        for(size_t i = 0; i < N_GC; i++){
            D[0][N_GC-1-i][j] = 2.*DBC_x - D[0][N_GC+1+i][j];
        }

        // Dy is cell-centered in x, do a linear interpolation
        for(size_t i = 0; i < N_GC; i++){
            D[1][N_GC-1-i][j] = 2.*DBC_y - D[1][N_GC+i][j];
        }

        // Dz is cell-centered in x, enforces that D_z = 0 on the boundary
        for(size_t i = 0; i < N_GC; i++){
            D[2][N_GC-1-i][j] = 2.*DBC_z -D[2][N_GC+i][j];
        }
    }

    exchng2Vector(D, N_GC, comm1D, nbrleft, nbrright);
}

/*
        Sets upper boundary condition on Displacement field
        Inputs: D: magnetic field as a vector field
               bparams: BandBCParams object containing information about boundary conditions
               Ny: extent of combined domain in partitioned direction
               N_GC: number of ghost cells
               t: time in reduced units
               comm1D: MPI communicator for decomposed domain
               world_rank: rank of current process
               Ny_locs and starts: vectors containing the extent and starting indices of the domain in the decomposed direction
               nbrleft and nbrright: the ranks of the processes to the left and right of the current process
               dm: Domain object containing information about the simulation domain
        Output: D with updated boundary values
*/

void UpperBoundary_D(VectorField& D, const VectorField& V, const VectorField& B, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright, double t, vConfig_params& driver, double y_min, double dy)
{
    //Exchange ghost cells
    //exchng2Vector(D, N_GC, comm1D, nbrleft, nbrright);

    //Loops all the way up to last physical cell in the y-direction
    for(size_t j = N_GC; j < D.shape()[2] - N_GC; j++){


        // B components at lower boundary... call B_Boundary condition function?

        // D_BC = -V x B = 0 since V = 0 at the upper boundary
        double DBC_x =  0.0;   
        double DBC_y = 0.0;   
        double DBC_z =  0.0;        

        for(size_t i = 0; i < N_GC; i++){
            D[0][D.shape()[1] - N_GC + 1 + i][j] = 2.*DBC_x - D[0][D.shape()[1]-N_GC-1-i][j];
            D[1][D.shape()[1]-N_GC+i][j] = 2.*DBC_y - D[1][D.shape()[1]-N_GC-1-i][j];
            D[2][D.shape()[1]-N_GC+i][j] = 2.*DBC_z -D[2][D.shape()[1]-N_GC-1-i][j];
        }

    }

    exchng2Vector(D, N_GC, comm1D, nbrleft, nbrright);
}





/*
        Sets boundary condition on magnetic field
        Inputs: B: magnetic field as a vector field
               bparams: BandBCParams object containing information about boundary conditions
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

    for(size_t i=0; i<N_GC; i++){
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

            /*
                Lower boundary x=0
            */
            //What are we setting B[0][N_GC][j] to?
            B[0][N_GC][j];
            B[1][N_GC-1-i][j] = B[1][N_GC+i][j]; 
            B[2][N_GC-1-i][j] = B[2][N_GC+i][j];

            /*
                Upper boundary x=Lx
            */
            // What are we setting B[0][B.shape()[1]-N_GC][j]?
            B[0][B.shape()[1] - N_GC][j];
            B[1][B.shape()[1]-N_GC+i][j] = B[1][B.shape()[1]-N_GC-1-i][j];
            B[2][B.shape()[1]-N_GC+i][j] = B[2][B.shape()[1]-N_GC-1-i][j];

            
            }
        }
    

    exchng2Vector(B, N_GC, comm1D, nbrleft, nbrright);

    return;
}

/*
        Sets boundary conditions on electric field
*/
void E_BoundaryConditions(VectorField & E, const BandBCParams & bparams, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright)
{
    exchng2Vector(E, N_GC, comm1D, nbrleft, nbrright);

    for(size_t i=0; i<N_GC; i++){
        for(size_t j=N_GC; j<E.shape()[2]-N_GC; j++){
            /*
                Lower boundary x=0
            */
            if(bparams.E_perp_lower == "perfect_conducting"){
                E[0][N_GC-1-i][j] = E[0][N_GC+i][j];
            }
            else if(bparams.E_perp_lower == "zero"){
                E[0][N_GC-1-i][j] = -E[0][N_GC+i][j];
            }
            if(bparams.E_parallel_lower == "perfect_conducting"){
                E[1][N_GC-1-i][j] = -E[1][N_GC+1+i][j]; //Ey = 0
                E[2][N_GC-1-i][j] = -E[2][N_GC+1+i][j]; //Ez = 0
                E[1][N_GC][j] = 0.;
                E[2][N_GC][j] = 0.;
            }
            else if(bparams.E_parallel_lower == "continuous"){
                E[1][N_GC-1-i][j] = E[1][N_GC+1+i][j];
                E[2][N_GC-1-i][j] = E[2][N_GC+1+i][j];
            }
            /*
                Upper boundary x=Lx
            */
            if(bparams.E_perp_upper == "vacuum"){
                E[0][E.shape()[1]-N_GC-1+i][j] = -E[0][E.shape()[1]-N_GC-2-i][j];
            }
            else if(bparams.E_perp_upper == "continuous"){  // FIX: was "continous" (typo)
                E[0][E.shape()[1]-N_GC-1+i][j] = E[0][E.shape()[1]-N_GC-2-i][j];
            }
            if(bparams.E_parallel_upper == "vacuum"){
                E[1][E.shape()[1]-N_GC+i][j] = -E[1][E.shape()[1]-N_GC-2-i][j];
                E[2][E.shape()[1]-N_GC+i][j] = -E[2][E.shape()[1]-N_GC-2-i][j];
            }
            else if(bparams.E_parallel_upper == "continuous"){
                E[1][E.shape()[1]-N_GC+i][j] = E[1][E.shape()[1]-N_GC-2-i][j];
                E[2][E.shape()[1]-N_GC+i][j] = E[2][E.shape()[1]-N_GC-2-i][j];
            }
        }
    }

    return;
}

/*
        Sets boundary conditions on H field
*/
void H_BoundaryConditions(VectorField & H, const BandBCParams & bparams, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright)
{
    exchng2Vector(H, N_GC, comm1D, nbrleft, nbrright);

    for(size_t i=0; i<N_GC; i++){
        for(size_t j=N_GC; j<H.shape()[2]-N_GC; j++){
            /*
                Lower boundary x=0
            */
            if(bparams.H_perp_lower == "perfect_conducting"){
                H[0][N_GC-1-i][j] = H[0][N_GC+i][j];
            }
            else if(bparams.H_perp_lower == "zero"){
                H[0][N_GC-1-i][j] = -H[0][N_GC+i][j];
            }
            if(bparams.H_parallel_lower == "perfect_conducting"){
                H[1][N_GC-1-i][j] = -H[1][N_GC+1+i][j]; //Hy = 0
                H[2][N_GC-1-i][j] = -H[2][N_GC+1+i][j]; //Hz = 0
                H[1][N_GC][j] = 0.;
                H[2][N_GC][j] = 0.;
            }
            else if(bparams.H_parallel_lower == "continuous"){
                H[1][N_GC-1-i][j] = H[1][N_GC+1+i][j];
                H[2][N_GC-1-i][j] = H[2][N_GC+1+i][j];
            }
            /*
                Upper boundary x=Lx
            */
            if(bparams.H_perp_upper == "vacuum"){
                H[0][H.shape()[1]-N_GC-1+i][j] = -H[0][H.shape()[1]-N_GC-2-i][j];
            }
            else if(bparams.H_perp_upper == "continuous"){  // FIX: was "continous" (typo)
                H[0][H.shape()[1]-N_GC-1+i][j] = H[0][H.shape()[1]-N_GC-2-i][j];
            }
            if(bparams.H_parallel_upper == "vacuum"){
                H[1][H.shape()[1]-N_GC+i][j] = -H[1][H.shape()[1]-N_GC-2-i][j];
                H[2][H.shape()[1]-N_GC+i][j] = -H[2][H.shape()[1]-N_GC-2-i][j];
            }
            else if(bparams.H_parallel_upper == "continuous"){
                H[1][H.shape()[1]-N_GC+i][j] = H[1][H.shape()[1]-N_GC-2-i][j];
                H[2][H.shape()[1]-N_GC+i][j] = H[2][H.shape()[1]-N_GC-2-i][j];
            }
        }
    }

    return;
}

/*
        Sets boundary condition on displacement field D
        FIX: renamed from B_BoundaryConditions to D_BoundaryConditions
*/
void D_BoundaryConditions(VectorField & D, const BandBCParams & bparams, size_t Ny, size_t N_GC, double t, MPI_Comm comm1D, int world_rank, std::vector<int> & Ny_locs, std::vector<int> & starts, int nbrleft, int nbrright, const Domain & dm)
{
    exchng2Vector(D, N_GC, comm1D, nbrleft, nbrright);

    std::vector<double> y = dm.y;

    std::vector<double> Dx_BC(D.shape()[2]-2*N_GC);
    std::vector<double> Dy_BC(D.shape()[2]-2*N_GC);

    double Dy_const = 0.;
    double Dyshear, Dzshear;

    for(size_t j=0; j<D.shape()[2]-2*N_GC; j++){
        Dx_BC[j] = D[0][D.shape()[1]-N_GC-1][N_GC+j];
        Dy_const = Dy_const + 0.5*( D[1][D.shape()[1]-N_GC-2][N_GC+j] + D[1][D.shape()[1]-N_GC-1][N_GC+j] );
    }
    Dy_const = Dy_const/double(D.shape()[2]-2*N_GC);

    // FIX: renamed from Dy_BC_Calc to Dy_BC_Calc (now correctly defined below)
    Dy_BC_Calc(Dx_BC, Dy_BC, Dy_const, Ny, Ny_locs, starts, world_rank, comm1D);

    for(size_t i=0; i<N_GC; i++){
        for(size_t j=N_GC; j<D.shape()[2]-N_GC; j++){

            /*
                Lower boundary x=0
            */
            if(bparams.D_perp_lower == "perfect_conducting"){
                D[0][N_GC-1-i][j] = -D[0][N_GC+i+1][j]; //Dx = 0
            }
            else if(bparams.D_perp_lower == "continuous"){
                // FIX: Bparams -> bparams, theta_B -> theta_D
                D[0][N_GC-1-i][j] = 2.*( bparams.D_pol_init*sin(bparams.theta_D) ) - D[0][N_GC+i+1][j];
            }

            if(bparams.D_parallel_lower == "perfect_conducting"){
                D[1][N_GC-1-i][j] = D[1][N_GC+i][j];
                D[2][N_GC-1-i][j] = D[2][N_GC+i][j];
            }
            else if(bparams.D_parallel_lower == "zero"){
                D[1][N_GC-1-i][j] = -D[1][N_GC+i][j]; //Dy = 0
                D[2][N_GC-1-i][j] = -D[2][N_GC+i][j]; //Dz = 0
            }
            else if(bparams.D_parallel_lower == "shearing"){
                // FIX: Dshear now correctly defined below
                Dshear(D[1][N_GC+i][j],D[2][N_GC+i][j],y[j-N_GC],t,dm,bparams,Dyshear,Dzshear);
                D[1][N_GC-1-i][j] = 2.*Dyshear - D[1][N_GC+i][j];
                D[2][N_GC-1-i][j] = 2.*Dzshear - D[2][N_GC+i][j];
            }

            /*
                Upper boundary x=Lx
            */
            if(bparams.D_pol_upper == "vacuum"){
                D[0][D.shape()[1]-N_GC+i][j] = D[0][D.shape()[1]-N_GC-2-i][j];
                D[1][D.shape()[1]-N_GC-1+i][j] = 2.*Dy_BC[j-N_GC] - D[1][D.shape()[1]-N_GC-2-i][j];
            }
            else if(bparams.D_pol_upper == "continuous"){
                D[0][D.shape()[1]-N_GC+i][j] = D[0][D.shape()[1]-N_GC-2-i][j];
                D[1][D.shape()[1]-N_GC-1+i][j] = D[1][D.shape()[1]-N_GC-2-i][j];
            }

            if(bparams.D_tor_upper == "vacuum"){
                D[2][D.shape()[1]-N_GC-1+i][j] = -D[2][D.shape()[1]-N_GC-2-i][j]; //Dz = 0
            }
            else if(bparams.D_tor_upper == "continuous"){  // FIX: was bparams.B_tor_upper
                D[2][D.shape()[1]-N_GC-1+i][j] = D[2][D.shape()[1]-N_GC-2-i][j];
            }
        }
    }

    exchng2Vector(D, N_GC, comm1D, nbrleft, nbrright);

    return;
}

/*
        Computes shearing Dy/Dz field at lower boundary
*/
void Dshear(double Dy_init, double Dz_init, double y, double t, const Domain & dm, const BandBCParams & bparams, double & Dyshear, double & Dzshear)
{
    double f_t;

    static double t_b = 0.;
    static double t_w = 100.*yr/t_0;
    static double y0 = 0.;
    static double y_w = dm.Ly;
    // FIX: was D_parallel_lower_shear_By/_Bz, now correctly D_parallel_lower_shear_Dy/_Dz
    static double Dyshear_max = bparams.D_parallel_lower_shear_Dy;
    static double Dzshear_max = bparams.D_parallel_lower_shear_Dz;

    if(bparams.D_parallel_shear_type == "uniform"){
        if( t >= t_b ){
            f_t = 1.;
            Dyshear = Dyshear_max*f_t;
            Dzshear = Dzshear_max*f_t;
        }
        else{
            Dyshear = -Dy_init;
            Dzshear = -Dz_init;
        }
    }
    else if(bparams.D_parallel_shear_type == "current_sheet"){
        if( t >= t_b ){
            f_t = 1.;
            Dyshear = 0.;
            Dzshear = f_t*Dzshear_max*sin(2.*pi*y/dm.Ly)*exp( -(y-y0)*(y-y0)/(2.*y_w*y_w) );
        }
        else{
            Dyshear = -Dy_init;
            Dzshear = -Dz_init;
        }
    }

    return;
}

/*
        Computes Dy given Dx at upper boundary using FFT for potential field condition.
        FIX: renamed from By_BC_Calc to Dy_BC_Calc, fixed &Bx.front() -> &Dx.front(),
             &By_BC.front() -> &Dy_BC.front()
*/
void Dy_BC_Calc(std::vector<double> & Dx_BC, std::vector<double> & Dy_BC, double Dy_const, size_t Ny, std::vector<int> & Ny_locs, std::vector<int> & starts, int world_rank, MPI_Comm comm1D)
{
    std::vector<double> Dx;
    std::vector<double> Dy;
    std::vector<double> DyTemp;
    if(world_rank == 0){
        Dx.resize(Ny);
        Dy.resize(Ny);
        DyTemp.resize(Ny);
    }

    // FIX: was &Bx.front(), now correctly &Dx.front()
    MPI_Gatherv(&Dx_BC.front(), Ny_locs[world_rank], MPI_DOUBLE, &Dx.front(), Ny_locs.data(), starts.data(), MPI_DOUBLE, 0, comm1D);

    static int num_procs = std::size(Ny_locs);
    std::vector<double> Dy_constVec(num_procs);
    MPI_Allgather(&Dy_const, 1, MPI_DOUBLE, Dy_constVec.data(), 1, MPI_DOUBLE, comm1D);

    if(world_rank == 0){
        int n = Dx.size();

        gsl_fft_real_wavetable * real;
        gsl_fft_halfcomplex_wavetable * hc;
        gsl_fft_real_workspace * work;

        work = gsl_fft_real_workspace_alloc(n);
        real = gsl_fft_real_wavetable_alloc(n);

        gsl_fft_real_transform(Dx.data(), 1, n, real, work);
        gsl_fft_real_wavetable_free(real);

        DyTemp[0] = std::reduce(Dy_constVec.begin(), Dy_constVec.end())*double(n)/double(num_procs);
        for(int j = 1; j<n/2; j++){
            DyTemp[2*j-1] = Dx[2*j];
            DyTemp[2*j] = -Dx[2*j-1];
        }

        hc = gsl_fft_halfcomplex_wavetable_alloc(n);
        gsl_fft_halfcomplex_inverse(DyTemp.data(), 1, n, hc, work);
        gsl_fft_halfcomplex_wavetable_free(hc);
        gsl_fft_real_workspace_free(work);

        for(int j = 0; j<n; j++){
            if( j == 0 ) Dy[j] = 0.5*(DyTemp[n-1] + DyTemp[0]);
            else Dy[j] = 0.5*(DyTemp[j-1] + DyTemp[j]);
        }
    }

    // FIX: was &By_BC.front(), now correctly &Dy_BC.front()
    MPI_Scatterv(&Dy.front(), Ny_locs.data(), starts.data(), MPI_DOUBLE, &Dy_BC.front(), Ny_locs.data()[world_rank], MPI_DOUBLE, 0, comm1D);

    return;
}

// only need one block for all three components
// make print statements
// vanishing on the sides, lower boundary have perturbation, outer boundary look at paper
// all functions here must match in header file

// MP edits begin below

/*
Defines initial velocity profile at the bottom boundary
(star sturface). v = -v_0 * exp(-(y-y_0)/2 sigma_y^2) * sin(2pi*f/T) z_hat.

*/

struct VelocityParams{

    double v_0; // initial velocity amplitude
    double y0; // center of Gaussian
    double sigma_y; // standard deviation of the Gaussian (width)
    double freq; // f/T parameter
};

double BoundaryVz(double y_boundary, double t, const VelocityParams& vparams){

    return vparams.v_0 * exp((-(y_boundary - vparams.y0))/(2 * vparams.sigma_y * vparams.sigma_y)) * sin(2 * pi * vparams.freq * t)

}

void D_BoundaryConditions(std::vector<double>& x, std::vector<double>& y, 




)

