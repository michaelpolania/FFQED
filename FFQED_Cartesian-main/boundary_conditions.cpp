#include <vector>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include "common.h"
#include "boundary_conditions.h"

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

    std::vector<double> y = dm.y; //cell centers in the y-direction. Note: no ghost cells.

    std::vector<double> Bx_BC(B.shape()[2]-2*N_GC);
    std::vector<double> By_BC(B.shape()[2]-2*N_GC);

    double By_const = 0.;
    double Byshear, Bzshear;

    for(size_t j=0; j<B.shape()[2]-2*N_GC; j++){
        Bx_BC[j] = B[0][B.shape()[1]-N_GC-1][N_GC+j];
        By_const = By_const + 0.5*( B[1][B.shape()[1]-N_GC-2][N_GC+j] + B[1][B.shape()[1]-N_GC-1][N_GC+j] );
    }
    By_const = By_const/double(B.shape()[2]-2*N_GC);

    By_BC_Calc(Bx_BC, By_BC, By_const, Ny, Ny_locs, starts, world_rank, comm1D);

    for(size_t i=0; i<N_GC; i++){
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

            /*
                Lower boundary x=0
            */
            if(bparams.B_perp_lower == "perfect_conducting"){
                B[0][N_GC-1-i][j] = -B[0][N_GC+i+1][j]; //Bx = 0
            }
            else if(bparams.B_perp_lower == "continuous"){
                // FIX: Bparams -> bparams
                B[0][N_GC-1-i][j] = 2.*( bparams.B_pol_init*sin(bparams.theta_B) ) - B[0][N_GC+i+1][j];
            }

            if(bparams.B_parallel_lower == "perfect_conducting"){
                B[1][N_GC-1-i][j] = B[1][N_GC+i][j];
                B[2][N_GC-1-i][j] = B[2][N_GC+i][j];
            }
            else if(bparams.B_parallel_lower == "zero"){
                B[1][N_GC-1-i][j] = -B[1][N_GC+i][j]; //By = 0
                B[2][N_GC-1-i][j] = -B[2][N_GC+i][j]; //Bz = 0
            }
            else if(bparams.B_parallel_lower == "shearing"){
                Bshear(B[1][N_GC+i][j],B[2][N_GC+i][j],y[j-N_GC],t,dm,bparams,Byshear,Bzshear);
                B[1][N_GC-1-i][j] = 2.*Byshear - B[1][N_GC+i][j];
                B[2][N_GC-1-i][j] = 2.*Bzshear - B[2][N_GC+i][j];
            }

            /*
                Upper boundary x=Lx
            */
            if(bparams.B_pol_upper == "vacuum"){
                B[0][B.shape()[1]-N_GC+i][j] = B[0][B.shape()[1]-N_GC-2-i][j];
                B[1][B.shape()[1]-N_GC-1+i][j] = 2.*By_BC[j-N_GC] - B[1][B.shape()[1]-N_GC-2-i][j];
            }
            else if(bparams.B_pol_upper == "continuous"){
                B[0][B.shape()[1]-N_GC+i][j] = B[0][B.shape()[1]-N_GC-2-i][j];
                B[1][B.shape()[1]-N_GC-1+i][j] = B[1][B.shape()[1]-N_GC-2-i][j];
            }

            if(bparams.B_tor_upper == "vacuum"){
                B[2][B.shape()[1]-N_GC-1+i][j] = -B[2][B.shape()[1]-N_GC-2-i][j]; //Bz = 0
            }
            else if(bparams.B_tor_upper == "continuous"){
                B[2][B.shape()[1]-N_GC-1+i][j] = B[2][B.shape()[1]-N_GC-2-i][j];
            }
        }
    }

    exchng2Vector(B, N_GC, comm1D, nbrleft, nbrright);

    return;
}

/*
        Computes shearing By/Bz field at lower boundary
*/
void Bshear(double By_init, double Bz_init, double y, double t, const Domain & dm, const BandBCParams & bparams, double & Byshear, double & Bzshear)
{
    double f_t;

    static double t_b = 0.;
    static double t_w = 100.*yr/t_0;
    static double y0 = 0.;
    static double y_w = dm.Ly;
    static double Byshear_max = bparams.B_parallel_lower_shear_By;
    static double Bzshear_max = bparams.B_parallel_lower_shear_Bz;

    if(bparams.B_parallel_shear_type == "uniform"){
        if( t >= t_b ){
            f_t = 1.;
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
            f_t = 1.;
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
        Computes By given Bx at upper boundary using FFT for potential field condition.
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

    MPI_Gatherv(&Bx_BC.front(), Ny_locs[world_rank], MPI_DOUBLE, &Bx.front(), Ny_locs.data(), starts.data(), MPI_DOUBLE, 0, comm1D);

    static int num_procs = std::size(Ny_locs);
    std::vector<double> By_constVec(num_procs);
    MPI_Allgather(&By_const, 1, MPI_DOUBLE, By_constVec.data(), 1, MPI_DOUBLE, comm1D);

    if(world_rank == 0){
        int n = Bx.size();

        gsl_fft_real_wavetable * real;
        gsl_fft_halfcomplex_wavetable * hc;
        gsl_fft_real_workspace * work;

        work = gsl_fft_real_workspace_alloc(n);
        real = gsl_fft_real_wavetable_alloc(n);

        gsl_fft_real_transform(Bx.data(), 1, n, real, work);
        gsl_fft_real_wavetable_free(real);

        ByTemp[0] = std::reduce(By_constVec.begin(), By_constVec.end())*double(n)/double(num_procs);
        for(int j = 1; j<n/2; j++){
            ByTemp[2*j-1] = Bx[2*j];
            ByTemp[2*j] = -Bx[2*j-1];
        }

        hc = gsl_fft_halfcomplex_wavetable_alloc(n);
        gsl_fft_halfcomplex_inverse(ByTemp.data(), 1, n, hc, work);
        gsl_fft_halfcomplex_wavetable_free(hc);
        gsl_fft_real_workspace_free(work);

        for(int j = 0; j<n; j++){
            if( j == 0 ) By[j] = 0.5*(ByTemp[n-1] + ByTemp[0]);
            else By[j] = 0.5*(ByTemp[j-1] + ByTemp[j]);
        }
    }

    MPI_Scatterv(&By.front(), Ny_locs.data(), starts.data(), MPI_DOUBLE, &By_BC.front(), Ny_locs.data()[world_rank], MPI_DOUBLE, 0, comm1D);

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