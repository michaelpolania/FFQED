#include <vector>
#include <iostream>
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

void LowerBoundary_D(std::vector<double> &y, VectorField& D, const VectorField& B, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright, double t, vConfig_params& driver)
{

    //Loops all the way up to last physical cell in the y-direction
    for(size_t j = N_GC; j < D.shape()[2] - N_GC; j++){
        for(size_t i = 0; i < N_GC; i++){

            double y_j = y[j];
            double y_j_avg = 0.5 * (y[j] + y[j-1]);
            double y_j_minus_1 = y[j-1];
            double DBC_x = 0.5 * Initialvz(y_j, t, &driver) * (0.5*(B[1][N_GC][j] + B[1][N_GC][j+1]) + 0.5*(B[1][N_GC-1][j] + B[1][N_GC-1][j+1]));
            double DBC_y = -0.5 * (0.5 * Initialvz(y_j, t, &driver) * (B[0][N_GC][j] + B[0][N_GC + 1][j]) + 0.5 * Initialvz(y_j_minus_1, t, &driver) * (B[0][N_GC][j-1] + B[0][N_GC + 1][j-1]));
            double DBC_z = 0.0;

            //D_x already lives on the boundary
            D[0][N_GC][j] = DBC_x;

            //std::cout << "DBC_y = " << DBC_y << std::endl; 
            D[0][N_GC-1-i][j] = 2.*DBC_x - D[0][N_GC+1+i][j];
            D[1][N_GC-1-i][j] = 2.*DBC_y - D[1][N_GC+i][j];
            D[2][N_GC-1-i][j] = 2.*DBC_z - D[2][N_GC+i][j];
            
        }
    }
    exchng2Vector(D, N_GC, comm1D, nbrleft, nbrright);    
    return;
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

void UpperBoundary_D(VectorField& D, const VectorField& B, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright, double t, vConfig_params& driver)
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

            B[0][N_GC - 1 - i][j] = B[0][N_GC + i + 1][j]; 
            B[1][N_GC - 1 - i][j] = B[1][N_GC + i][j]; 
            B[2][N_GC - 1 - i][j] = B[2][N_GC + i][j];
        

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
