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
        double DBC_x =  Vz * B[1][N_GC][j];   
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
        Sets upper boundary condition on Displacement field (enforces value is zero at boundary)
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


        for (size_t i = 0; i < N_GC; i++) {
    
            D[0][D.shape()[1] - 2*N_GC + N_GC + i][j] = -D[0][D.shape()[1] - 2*N_GC + N_GC - 2 - i][j];
            D[1][D.shape()[1] - 2*N_GC + N_GC - 1 + i][j] = -D[1][D.shape()[1] - 2*N_GC + N_GC - 2 - i][j];
            D[2][D.shape()[1] - 2*N_GC + N_GC - 1 + i][j] = -D[2][D.shape()[1] - 2*N_GC + N_GC - 2 - i][j];
}

    }

    exchng2Vector(D, N_GC, comm1D, nbrleft, nbrright);
}


/*
        Sets boundary condition on magnetic field (continous boundary condition)
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
    //exchng2Vector(B, N_GC, comm1D, nbrleft, nbrright);

    for(size_t i=0; i<N_GC; i++){
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

            /*
                Lower boundary x=0
            */
            
            //B[0][N_GC][j];
            B[1][N_GC-1-i][j] = B[1][N_GC+i][j]; 
            B[2][N_GC-1-i][j] = B[2][N_GC+i][j];

            /*
                Upper boundary x=Lx
            */
            
            B[0][B.shape()[1] - 2*N_GC + N_GC + i][j] = B[0][B.shape()[1] - 2*N_GC + N_GC - 2 - i][j];
            B[1][B.shape()[1] - 2*N_GC + N_GC - 1 + i][j] = B[1][B.shape()[1] - 2*N_GC + N_GC - 2 - i][j];
            B[2][B.shape()[1] - 2*N_GC + N_GC - 1 + i][j] = B[2][B.shape()[1] - 2*N_GC + N_GC - 2 - i][j];

            
            }
        }
    

    exchng2Vector(B, N_GC, comm1D, nbrleft, nbrright);

    return;
}
