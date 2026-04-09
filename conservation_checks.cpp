#include <vector>
#include <iostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

#include "common.h"
#include "initial_conditions.h"

/*
        Computes the divergence of the magnetic field in each cell and returns its total across current process
        Input: B: vector field containing magnetic field
               N_GC: number of ghost cells at each edge
        Output: mean value of magnetic flux out of each cell. Should be zero to within machine precision.
*/
double divB_Calculator(VectorField & B, size_t N_GC)
{

    double sum = 0.;
    for(size_t i=N_GC; i<B.shape()[1]-N_GC; i++){
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){
            sum += -B[0][i][j] + B[0][i+1][j] - B[1][i][j] + B[1][i][j+1];
        }
    }

    return sum;

}

/*
        Collects the mean div(B) from each process and averages them, then prints out the full domain-averaged div(B) value
        Input: B: vector field containing magnetic field
               N_GC: number of ghost cells at each edge
               comm1D: communicator between processes
               wordl_rank: rank of current process
               Ny_locs: vector containing the extent of the domain in the decomposed direction
*/
void divB_Monitor(VectorField & B, size_t N_GC,  MPI_Comm comm1D, int world_rank, std::vector<int> & Ny_locs)
{

    //Compute divergence of B for each cell, then take average and print out to track this. Should be 0 to within machine precision.
    static int num_procs = std::size(Ny_locs); //gives the number of processes
    std::vector<double> divBVec(num_procs);
    double divB = divB_Calculator(B, N_GC);
    MPI_Allgather(&divB, 1, MPI_DOUBLE, divBVec.data(), 1, MPI_DOUBLE, comm1D);
    double divBSum;
    double num_cells = double( num_procs*(B.shape()[1]-2*N_GC)*(B.shape()[2]-2*N_GC) ); //total number of cells across entire simulation
    if(world_rank == 0){
        divBSum = std::reduce(divBVec.begin(), divBVec.end())/double(num_procs);
        std::cout << "Average div(B) = " << divBSum/num_cells << std::endl;
    }

    return;
}

/*
    Calculates total magnetic field energy, Joule heating rate and outward Poynting flux integrated over decomposed domain
    Inputs: B, E, J: magnetic field, electric field times c, current density, velocity field in reduced units
            N_GC: number of ghost cells
            eta_O: Ohmic diffusivity in reduced units
            t: time in reduced units
            dm: Domain object containing information about the simulation domain
            bparams: BandBCParams object containing information about the initial magnetic field and boundary conditions
    Output: U_B, JouleH, PoyntingF: magnetic field energy, Joule heating rate, Poynting flux integrated over decomposed domain
*/
void EnergyConservation(VectorField & B, VectorField & E, VectorField & J, ScalarField & eta_O, size_t N_GC, double t, const Domain & dm, const BandBCParams & bparams, double & U_B, double & JouleH, double & PoyntingF)
{

    double U_BSum = 0.; //B^2/(8*pi)*cell volume
    double JouleSum = 0.; //J^2/sigma*cell volume
    double PoyntingSum = 0.; //Poynting flux out of cell

    for(size_t i=N_GC; i<B.shape()[1]-N_GC-1; i++){
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

            U_BSum += 1./(8.*pi)*dm.Deltax[i]*dm.Deltay*( 0.5*( pow(B[0][i][j],2.) + pow(B[0][i+1][j],2.) ) + 0.5*( pow(B[1][i][j],2.) + pow(B[1][i][j+1],2.) ) + pow(B[2][i][j],2.) );

            JouleSum += -1./(4.*pi)*dm.Deltax[i]*dm.Deltay*( 0.5*( E[0][i][j]*J[0][i][j] + E[0][i][j+1]*J[0][i][j+1] ) + 0.5*( E[1][i][j]*J[1][i][j] + E[1][i+1][j]*J[1][i+1][j] )
                                                            + 0.25*( E[2][i][j]*J[2][i][j] + E[2][i+1][j]*J[2][i+1][j] + E[2][i][j+1]*J[2][i][j+1] + E[2][i+1][j+1]*J[2][i+1][j+1] ) );

//            PF = -1./(4.*pi)*( -dm.Deltay*( E[1][i][j]*(B[2][i-1][j]+B[2][i][j])/2. - 1./4.*( E[2][i][j]*(B[1][i-1][j]+B[1][i][j]) + E[2][i][j+1]*(B[1][i-1][j+1]+B[1][i][j+1]) ) )
//                            - dm.Deltax[i]*( 1./4.*( E[2][i][j]*(B[0][i][j-1]+B[1][i][j]) + E[2][i+1][j]*(B[0][i+1][j-1]+B[0][i+1][j]) ) - E[0][i][j]*(B[2][i][j-1]+B[2][i][j])/2. ) );

            if( i == N_GC ){ //only compute Poynting flux at non-periodic outer surfaces
                PoyntingSum += -1./(4.*pi)*( -dm.Deltay*( E[1][i][j]*(B[2][i-1][j]+B[2][i][j])/2. - 1./4.*( E[2][i][j]*(B[1][i-1][j]+B[1][i][j]) + E[2][i][j+1]*(B[1][i-1][j+1]+B[1][i][j+1]) ) ) );
            }
            else if( i == B.shape()[1]-N_GC-2 ){
                PoyntingSum += -1./(4.*pi)*( dm.Deltay*( E[1][i+1][j]*(B[2][i][j]+B[2][i+1][j])/2. - 1./4.*( E[2][i+1][j]*(B[1][i][j]+B[1][i+1][j]) + E[2][i+1][j+1]*(B[1][i][j+1]+B[1][i+1][j+1]) ) ) );
            }
        }
    }

    //Set output variables equal to values computed by integrating over simulation domain
    U_B = U_BSum; //magnetic energy density integrated over x-y area; multiply by B_0^2*L_0^2 to get quantity in erg/cm
    JouleH = JouleSum; //Joule heating rate integrated over x-y area: multiply by B_0^2*L_0^2 to get quantity in erg/cm/reduced time
    PoyntingF = PoyntingSum; //Poynting flux through x and y : multiply by B_0^2*L_0^2 to get quantity in erg/cm/reduced time

    return;
}
