#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <utility>
#include <gsl/gsl_poly.h>
#include "boost/multi_array.hpp"

#include "common.h"
#include "initial_conditions.h"
#include "field_evolution.h"

void Compute_Rho(VectorField & D, ScalarField & Rho, const Domain & dm)
{
   

    //Rho is needed in the ghost cells as well (the face-centred current at the first/last interior
    //faces averages the cell-centred current from the neighbouring ghost cell), so compute it
    //everywhere the stencil D[.][i+1][j], D[.][i][j+1] is available.
    for(size_t i=0; i<D.shape()[1]-1; i++){
      for(size_t j=0; j<D.shape()[2]-1; j++){

            Rho[i][j] = (D[0][i+1][j] - D[0][i][j])/dm.Deltax[i] + (D[1][i][j+1] - D[1][i][j])/dm.Deltay;
           // std::cout << "D[0][i][j] = " << D[0][i][j] << std::endl;
            // std::cout << "D[0][i][j] = " << D[0][i][j] << std::endl;
            //std::cout << "D[0][i+1][j] = " << D[0][i+1][j] << std::endl;
            //std::cout << "D[1][i][j] = " << D[1][i][j] << std::endl;
            //std::cout << "D[1][i][j+1] = " << D[1][i][j+1] << std::endl;
            

        }
    }
    return;
}

/*
    Computes cell-centered averaged values for the x and y components of the B or D field.

    Inputs:
        A: Arbitrary VectorField object
        component_index: 0 (x), 1 (y)
        i_offset: offset for i index (e.g. i+1, i-1)
        j_offset: offset for j index (e.g. j+1, j-1)

    Example: If you want to calculate By_i_jminus1, you make the following function call Compute_B_to_cell_center(B, 1, i, j-1)

    Output: Cell centered average value (Bx_bar or By_bar)
        
*/
double Compute_A_to_cell_center(VectorField & A, int component_index, int i_offset, int j_offset) {
    
    // Computes the cell-centered averaged value for Bx
    if (component_index == 0) {
        return 0.5 * (A[0][i_offset][j_offset] + A[0][i_offset + 1][j_offset]);
    }
    
    // Computes the cell-centered averaged value for By
    else if (component_index == 1) {
        return 0.5 * (A[1][i_offset][j_offset] + A[1][i_offset][j_offset + 1]);
    }

    return 0.0;
}

double slope_calc (VectorField & A, int component_index, int i_offset, int j_offset, int slope_direction, double Deltax, double Deltay){
    
    //Compute_A_to_cell_center(A, component_index, i_offset, j_offset);

    if (slope_direction == 0) {

        double val_center = A[component_index][i_offset][j_offset];
        double val_right  = A[component_index][i_offset + 1][j_offset];
        double val_left   = A[component_index][i_offset - 1][j_offset];

        double slope_central = (val_right - val_left) / (2.0 * Deltax);
        double slope_left = 2.0 * (val_center - val_left)/(Deltax);
        double slope_right = 2.0 * (val_right - val_center)/(Deltax);

        return minmod(slope_central, slope_left, slope_right);
    }

    else if (slope_direction == 1){
        
        double val_center = A[component_index][i_offset][j_offset];
        double val_up  = A[component_index][i_offset][j_offset + 1];
        double val_down   = A[component_index][i_offset][j_offset - 1];

        double slope_central = (val_up - val_down) / (2.0 * Deltay);
        double slope_down = 2.0 * (val_center - val_down)/(Deltay);
        double slope_up = 2.0 * (val_up - val_center)/(Deltay);

        return minmod(slope_central, slope_down, slope_up);

    }

    return 0.0;

}

/*
Computes E or H from the cell-faced reconstructions of B and D.

Inputs: 
    i_offset: offset for i index (e.g. i+1, i-1)
    j_offset: offset for j index (e.g. j+1, j-1)
    component_index: 0 (x-dir), 1 (y-dir), 2 (z dir)
    field: "E" for electric field, "H" for auxillarly field
    plus_or_minus: "minus", "plus"
    slope_direction: "x" for x direction or "y" for y direction
    B: Cell-centered average for magnetic field
    D: Cell-centered average for displacement field

*/

double Compute_EH_from_BD(int i_offset, int j_offset, int component_index, const std::string field, const std::string plus_or_minus, const std::string slope_direction, VectorField & B, VectorField & D, const SimParams & params, double Deltax, double Deltay){

    //Ternary operator (if slope_direction == 0 then use Deltax)
    const double delta = (slope_direction == "x") ? Deltax : Deltay;

    const double reconstruction_sign = (plus_or_minus == "minus") ? -1 : 1;

    double B_face[3];
    double D_face[3];

    int slope_dir_idx = (slope_direction == "x") ? 0 : 1;

        // Reconstruct B and D at the same face
    for (int c = 0; c < 3; ++c) {
        double slope_B = slope_calc(B, c, i_offset, j_offset, slope_dir_idx, Deltax, Deltay);
        double slope_D = slope_calc(D, c, i_offset, j_offset, slope_dir_idx, Deltax, Deltay);

        B_face[c] = B[c][i_offset][j_offset] + reconstruction_sign * 0.5 * delta * slope_B;
        D_face[c] = D[c][i_offset][j_offset] + reconstruction_sign * 0.5 * delta * slope_D;
    }

    // 3. Compute with QED Corrections
    if (params.QED_corrections) {
        double D_squared = D_face[0] * D_face[0] + D_face[1] * D_face[1] + D_face[2] * D_face[2];
        double B_squared = B_face[0] * B_face[0] + B_face[1] * B_face[1] + B_face[2] * B_face[2];

        double lambda = Solve_Lambda_Cubic(D_squared, B_squared);

        if (field == "E") {
            return D_face[component_index] / lambda;
        } 
        else if (field == "H") {
            return lambda * B_face[component_index];
        }
    }   
    // 4. Compute standard linear fields
    else {
        if (field == "E") {
            return D_face[component_index];
        }
        else if (field == "H") {
            return B_face[component_index];
        }
    }

    return 0.0;
}

/*
    Computes Ez fluxes using slope limiters to update B_x.

    Inputs: i, j, D, B, Deltax, and Deltay

    Outputs: Ez{i-1/2,j-1/2} and Ez{i-1/2, j+1/2} fluxes computed using Toth 2000 Eq(19).
*/

std::pair<double, double> Ez_Flux_Calculation_Qx(int i, int j, VectorField & D, VectorField & B, const SimParams & params, double Deltax, double Deltay){
    
    //Fluxes for Ez{i,j} aka Ez_{i-1/2, j-1/2}
    //G_1{i,j-1/2} calculation (G_1)

    double Bx_ij_up = B[0][i][j] - 0.5 * Deltay * slope_calc(B, 0, i, j, 1, Deltax, Deltay);
    double Bx_i_jminus1_down = B[0][i][j-1] + 0.5 * Deltay * slope_calc(B, 0, i, j-1, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper 1

    double G_1_up = Compute_EH_from_BD(i, j, 2, "E", "minus", "y", B, D, params, Deltax, Deltay);
    double G_1_down = Compute_EH_from_BD(i, j-1, 2, "E", "plus", "y", B, D, params, Deltax, Deltay);

    double G_1 = 0.5 * (G_1_up + G_1_down) - 0.5 * (Bx_ij_up - Bx_i_jminus1_down);

    //G_1{i-1,j-1/2} calculation (G_2)
    
    //NOTE: B passed in here is already cell-centred (B_avg); do not average it a second time.
    double Bx_iminus1_j_up = B[0][i-1][j] - 0.5 * Deltay * slope_calc(B, 0, i-1, j, 1, Deltax, Deltay);
    double Bx_ijminus1_down = B[0][i-1][j-1] + 0.5 * Deltay * slope_calc(B, 0, i-1, j-1, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper

    double G_2_up = Compute_EH_from_BD(i-1, j, 2, "E", "minus", "y", B, D, params, Deltax, Deltay);
    double G_2_down = Compute_EH_from_BD(i-1, j-1, 2, "E", "plus", "y", B, D, params, Deltax, Deltay);

    double G_2 = 0.5 * (G_2_up + G_2_down) - 0.5 * (Bx_iminus1_j_up - Bx_ijminus1_down);

    //F_2{i-1/2,j} calculation (F_1)

    double By_ij_up = B[1][i][j] - 0.5 * Deltax * slope_calc(B, 1, i, j, 0, Deltax, Deltay);
    double By_iminus1_j_down = B[1][i-1][j] + 0.5 * Deltax * slope_calc(B, 1, i-1, j, 0, Deltax, Deltay);

    double F_1_up = Compute_EH_from_BD(i, j, 2, "E", "minus", "x", B, D, params, Deltax, Deltay);
    double F_1_down = Compute_EH_from_BD(i-1, j, 2, "E", "plus", "x", B, D, params, Deltax, Deltay);

    //Note: Ez is negative in the F matrix in Yu paper

    double F_1 = 0.5 * (-F_1_up - F_1_down) - 0.5 * (By_ij_up - By_iminus1_j_down);
    //F_2{i-1/2,j-1} calculation (F_2)

    double By_i_jminus1_up = B[1][i][j-1] - 0.5 * Deltax * slope_calc(B, 1, i, j-1, 0, Deltax, Deltay);
    double By_ij_minus1_down = B[1][i-1][j-1] + 0.5 * Deltax * slope_calc(B, 1, i-1, j-1, 0, Deltax, Deltay);

    
    double F_2_up = Compute_EH_from_BD(i, j-1, 2, "E", "minus", "x", B, D, params, Deltax, Deltay);
    double F_2_down = Compute_EH_from_BD(i-1, j-1, 2, "E", "plus", "x", B, D, params, Deltax, Deltay);

    //Note: Ez is negative in the F matrix in Yu paper

    double F_2 = 0.5 * (-F_2_up - F_2_down) - 0.5 * (By_i_jminus1_up - By_ij_minus1_down);

    //Fluxes for Ez{i,j+1} aka Ez_{i-1, j+1/2}
    //G_1{i,j+1/2} calculation (G_3)

    double Bx_i_jplus1_up =  B[0][i][j+1] - 0.5 * Deltay * slope_calc(B, 0, i, j+1, 1, Deltax, Deltay);
    double Bx_ij_down = B[0][i][j] + 0.5 * Deltay * slope_calc(B, 0, i, j, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper
    double G_3_up = Compute_EH_from_BD(i, j+1, 2, "E", "minus", "y", B, D, params, Deltax, Deltay);
    double G_3_down = Compute_EH_from_BD(i, j, 2, "E", "plus", "y", B, D, params, Deltax, Deltay);

    double G_3 = 0.5 * (G_3_up + G_3_down) - 0.5 * (Bx_i_jplus1_up - Bx_ij_down);
    
    //G_1{i-1,j+1/2} calculation (G_4)

    double Bx_iminus1_jplus1_up =  B[0][i-1][j+1] - 0.5 * Deltay * slope_calc(B, 0, i-1, j+1, 1, Deltax, Deltay);
    double Bx_iminus1_j_down = B[0][i-1][j] + 0.5 * Deltay * slope_calc(B, 0, i-1, j, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper
    
    double G_4_up = Compute_EH_from_BD(i-1, j+1, 2, "E", "minus", "y", B, D, params, Deltax, Deltay);
    double G_4_down = Compute_EH_from_BD(i-1, j, 2, "E", "plus", "y", B, D, params, Deltax, Deltay);

    double G_4 = 0.5 * (G_4_up + G_4_down) - 0.5 * (Bx_iminus1_jplus1_up - Bx_iminus1_j_down);

    //F_2{i-1/2,j} calculation (F_3) same thing as F_1

    //Note: Ez is negative in the F matrix in Yu paper

    double F_3 = 0.5 * (-F_1_up - F_1_down) - 0.5 * (By_ij_up - By_iminus1_j_down);

    //F_2{i-1/2,j+1} calculation (F_4)

    double By_i_jplus1_up = B[1][i][j+1] - 0.5 * Deltax * slope_calc(B, 1, i, j+1, 0, Deltax, Deltay);
    double By_iminus1_jplus1_down = B[1][i-1][j+1] + 0.5 * Deltax * slope_calc(B, 1, i-1, j+1, 0, Deltax, Deltay);

    double F_4_up = Compute_EH_from_BD(i, j+1, 2, "E", "minus", "x", B, D, params, Deltax, Deltay);
    double F_4_down = Compute_EH_from_BD(i-1, j+1, 2, "E", "plus", "x", B, D, params, Deltax, Deltay); 

    double F_4 = 0.5 * (-F_4_up - F_4_down) - 0.5 * (By_i_jplus1_up - By_iminus1_jplus1_down);

    //Toth 2000 Eq(19)
    double Ez_ij = 0.25 * (-F_1 - F_2 + G_1 + G_2);
    double Ez_i_jplus1 = 0.25 * (-F_3 - F_4 + G_3 + G_4);
    
    return {Ez_ij, Ez_i_jplus1};

}

/*
    Computes Ez fluxes using slope limiters to update B_y.

    Inputs: i, j, D, B, Deltax, and Deltay

    Outputs: Ez{i-1/2,j-1/2} and Ez{i+1/2, j-1/2} fluxes computed using Toth 2000 Eq(19).
*/

std::pair<double, double> Ez_Flux_Calculation_Qy(int i, int j, VectorField & D, VectorField & B, const SimParams & params, double Deltax, double Deltay){

    //Fluxes for Ez{i,j} aka Ez_{i-1/2, j-1/2}
    //G_1{i,j-1/2} calculation (G_1)

    double Bx_ij_up = B[0][i][j] - 0.5 * Deltay * slope_calc(B, 0, i, j, 1, Deltax, Deltay);
    double Bx_i_jminus1_down = B[0][i][j-1] + 0.5 * Deltay * slope_calc(B, 0, i, j-1, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper
    double G_1_up = Compute_EH_from_BD(i, j, 2, "E", "minus", "y", B, D, params, Deltax, Deltay);
    double G_1_down = Compute_EH_from_BD(i, j-1, 2, "E", "plus", "y", B, D, params, Deltax, Deltay);
    
    
    double G_1 = 0.5 * (G_1_up + G_1_down) - 0.5 * (Bx_ij_up - Bx_i_jminus1_down);

    //G_1{i-1,j-1/2} calculation (G_2)

    double Bx_iminus1_j_up = B[0][i-1][j] - 0.5 * Deltay * slope_calc(B, 0, i-1, j, 1, Deltax, Deltay);
    double Bx_ijminus1_down = B[0][i-1][j-1] + 0.5 * Deltay * slope_calc(B, 0, i-1, j-1, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper 

    double G_2_up = Compute_EH_from_BD(i-1, j, 2, "E", "minus", "y", B, D, params, Deltax, Deltay);
    double G_2_down = Compute_EH_from_BD(i-1, j-1, 2, "E", "plus", "y", B, D, params, Deltax, Deltay);

    double G_2 = 0.5 * (G_2_up + G_2_down) - 0.5 * (Bx_iminus1_j_up - Bx_ijminus1_down);


    //F_2{i-1/2,j} calculation (F_1)

    double By_ij_up = B[1][i][j] - 0.5 * Deltax * slope_calc(B, 1, i, j, 0, Deltax, Deltay);
    double By_iminus1_j_down = B[1][i-1][j] + 0.5 * Deltax * slope_calc(B, 1, i-1, j, 0, Deltax, Deltay);

    double F_1_up = Compute_EH_from_BD(i, j, 2, "E", "minus", "x", B, D, params, Deltax, Deltay);
    double F_1_down = Compute_EH_from_BD(i-1, j, 2, "E", "plus", "x", B, D, params, Deltax, Deltay);

    //Note: Ez is negative in the F matrix in Yu paper

    double F_1 = 0.5 * (-F_1_up - F_1_down) - 0.5 * (By_ij_up - By_iminus1_j_down);

    //F_2{i-1/2,j-1} calculation (F_2)

    double By_i_jminus1_up = B[1][i][j-1] - 0.5 * Deltax * slope_calc(B, 1, i, j-1, 0, Deltax, Deltay);
    double By_ij_minus1_down = B[1][i-1][j-1] + 0.5 * Deltax * slope_calc(B, 1, i-1, j-1, 0, Deltax, Deltay);

    double F_2_up = Compute_EH_from_BD(i, j-1, 2, "E", "minus", "x", B, D, params, Deltax, Deltay);
    double F_2_down = Compute_EH_from_BD(i-1, j-1, 2, "E", "plus", "x", B, D, params, Deltax, Deltay);

    double F_2 = 0.5 * (-F_2_up - F_2_down) - 0.5 * (By_i_jminus1_up - By_ij_minus1_down);

    //E_z Fluxes for Ez{i+1,j} aka Ez_{i+1/2, j-1/2}

    //G_1{i,j-1/2} calculation (G_3), same as G_1

    //Note: Ez is positive in the G matrix in Yu paper

    double G_3 = 0.5 * (G_1_up + G_1_down) - 0.5 * (Bx_ij_up - Bx_i_jminus1_down);

    //G_1{i+1, j-1/2}  calculation (G_4)

    double Bx_iplus1_j_up =  B[0][i+1][j] - 0.5 * Deltay * slope_calc(B, 0, i+1, j, 1, Deltax, Deltay);
    double Bx_iplus1_jminus1_down = B[0][i+1][j-1] + 0.5 * Deltay * slope_calc(B, 0, i+1, j-1, 1, Deltax, Deltay);

    double G_4_up = Compute_EH_from_BD(i+1, j, 2, "E", "minus", "y", B, D, params, Deltax, Deltay);
    double G_4_down = Compute_EH_from_BD(i+1, j-1, 2, "E", "plus", "y", B, D, params, Deltax, Deltay);

    double G_4 = 0.5 * (G_4_up + G_4_down) - 0.5 * (Bx_iplus1_j_up - Bx_iplus1_jminus1_down);

    //F_1{i+1/2, j-1} calculation (F_3)

    double By_iplus1_jminus1_up =  B[1][i+1][j-1] - 0.5 * Deltax * slope_calc(B, 1, i+1, j-1, 0, Deltax, Deltay);
    double By_i_jminus1_down =  B[1][i][j-1] + 0.5 * Deltax *  slope_calc(B, 1, i, j-1, 0, Deltax, Deltay);

    double F_3_up = Compute_EH_from_BD(i+1, j-1, 2, "E", "minus", "x", B, D, params, Deltax, Deltay);
    double F_3_down = Compute_EH_from_BD(i, j-1, 2, "E", "plus", "x", B, D, params, Deltax, Deltay);

    double F_3 = 0.5 * (-F_3_up - F_3_down) - 0.5 * (By_iplus1_jminus1_up - By_i_jminus1_down);

    //F_1{i+1/2,j} calculation (F_4)

    double By_iplus1_j_up = B[1][i+1][j] - 0.5 * Deltax * slope_calc(B, 1, i+1, j, 0, Deltax, Deltay);
    double By_ij_down =  B[1][i][j] + 0.5 * Deltax * slope_calc(B, 1, i, j, 0, Deltax, Deltay);  

    double F_4_up = Compute_EH_from_BD(i+1, j, 2, "E", "minus", "x", B, D, params, Deltax, Deltay);
    double F_4_down = Compute_EH_from_BD(i, j, 2, "E", "plus", "x", B, D, params, Deltax, Deltay);

    double F_4 = 0.5 * (-F_4_up - F_4_down) - 0.5 * (By_iplus1_j_up - By_ij_down);

    double Ez_ij = 0.25 * (-F_1 - F_2 + G_1 + G_2);
    double Ez_iplus1_j = 0.25 * (-F_3 - F_4 + G_3 + G_4);

    return {Ez_ij, Ez_iplus1_j};

}

/*
    Computes reconstructed fluxes using slope limiters to update B_z.
    Inputs: i, j, D, B, Deltax, and Deltay

    Outputs: Update for Bz by computing F_3/F_4 (Ey), G_3/G_4 (Ex) fluxes computed using Yu 2011 Eq(9).
*/
double Qz_reconstruction(int i, int j, VectorField & D, VectorField & B, const SimParams & params, double Deltax, double Deltay){

    //G_3{i,j-1/2} calculation G_3

    double G_3_up = Compute_EH_from_BD(i, j, 0, "E", "minus", "y", B, D, params, Deltax, Deltay);
    double G_3_down = Compute_EH_from_BD(i, j-1, 0, "E", "plus", "y", B, D, params, Deltax, Deltay);
    
    double Bz_ij_up_y = B[2][i][j] - 0.5 * Deltay * slope_calc(B, 2, i, j, 1, Deltax, Deltay);
    double Bz_i_jminus1_down_y = B[2][i][j-1] + 0.5 * Deltay * slope_calc(B, 2, i, j-1, 1, Deltax, Deltay);

    double G_3 = 0.5 * (-G_3_up - G_3_down) - 0.5 * (Bz_ij_up_y - Bz_i_jminus1_down_y);

    //G_3{i,j+1/2} calculation G_4 

    double G_4_up = Compute_EH_from_BD(i, j+1, 0, "E", "minus", "y", B, D, params, Deltax, Deltay);
    double G_4_down = Compute_EH_from_BD(i, j, 0, "E", "plus", "y", B, D, params, Deltax, Deltay);

    double Bz_i_jplus1_up_y = B[2][i][j+1] - 0.5 * Deltay * slope_calc(B, 2, i, j+1, 1, Deltax, Deltay);
    double Bz_ij_down_y = B[2][i][j] + 0.5 * Deltay * slope_calc(B, 2, i, j, 1, Deltax, Deltay);

    double G_4 = 0.5 * (-G_4_up - G_4_down) - 0.5 * (Bz_i_jplus1_up_y - Bz_ij_down_y);

    //F_3{i-1/2,j} calculation (F_3) 

    double F_3_up = Compute_EH_from_BD(i, j, 1, "E", "minus", "x", B, D, params, Deltax, Deltay);
    double F_3_down = Compute_EH_from_BD(i-1, j, 1, "E", "plus", "x", B, D, params, Deltax, Deltay);
    
    double Bz_ij_up_x = B[2][i][j] - 0.5 * Deltax * slope_calc(B, 2, i, j, 0, Deltax, Deltay);
    double Bz_iminus1_j_down_x = B[2][i-1][j] + 0.5 * Deltax * slope_calc(B, 2, i-1, j, 0, Deltax, Deltay);

    double F_3 = 0.5 * (F_3_up + F_3_down) - 0.5 * (Bz_ij_up_x - Bz_iminus1_j_down_x);
    
    //F_3{i+1/2,j} calculation (F_4)

    double F_4_up = Compute_EH_from_BD(i+1, j, 1, "E", "minus", "x", B, D, params, Deltax, Deltay);
    double F_4_down = Compute_EH_from_BD(i, j, 1, "E", "plus", "x", B, D, params, Deltax, Deltay);
    
    double Bz_iplus1_j_up_x = B[2][i+1][j] - 0.5 * Deltax * slope_calc(B, 2, i+1, j, 0, Deltax, Deltay);
    double Bz_ij_down_x = B[2][i][j] + 0.5 * Deltax * slope_calc(B, 2, i, j, 0, Deltax, Deltay);

    double F_4 = 0.5 * (F_4_up + F_4_down) - 0.5 * (Bz_iplus1_j_up_x - Bz_ij_down_x);

    //Eq(9) in Yu paper
    double Qz = -(F_4 - F_3)/(Deltax) + (-G_4 + G_3)/(Deltay);

    return Qz;

}

/*
    Computes Hz fluxes using slope limiters to update D_x. Calculation is analagous to the Ez flux calculation for the Bx
    time evolution, hence, copied and pasted the Ez_Flux_Calculation_Qx function.

    Inputs: i, j, B, D, Deltax, and Deltay

    Outputs: Hz{i-1/2,j-1/2} and Hz{i-1/2, j+1/2} fluxes computed using Toth 2000 Eq(19).
*/
std::pair<double, double> Hz_Flux_Calculation_Fx(int i, int j, VectorField & B, VectorField & D, const SimParams & params, double Deltax, double Deltay){

    //Fluxes for Ez{i,j} aka Ez_{i-1/2, j-1/2}
    //G_1{i,j-1/2} calculation (G_1)

    double Dx_ij_up = D[0][i][j] - 0.5 * Deltay * slope_calc(D, 0, i, j, 1, Deltax, Deltay);
    double Dx_i_jminus1_down = D[0][i][j-1] + 0.5 * Deltay * slope_calc(D, 0, i, j-1, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper
    double G_1_up = Compute_EH_from_BD(i, j, 2, "H", "minus", "y", B, D, params, Deltax, Deltay);
    double G_1_down = Compute_EH_from_BD(i, j-1, 2, "H", "plus", "y", B, D, params, Deltax, Deltay);
    
    double G_1 = 0.5 * (-G_1_up - G_1_down) - 0.5 * (Dx_ij_up - Dx_i_jminus1_down);

    //G_1{i-1,j-1/2} calculation (G_2)
    
    double Dx_iminus1_j_up = D[0][i-1][j] - 0.5 * Deltay * slope_calc(D, 0, i-1, j, 1, Deltax, Deltay);
    double Dx_ijminus1_down = D[0][i-1][j-1] + 0.5 * Deltay * slope_calc(D, 0, i-1, j-1, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper 
    double G_2_up = Compute_EH_from_BD(i-1, j, 2, "H", "minus", "y", B, D, params, Deltax, Deltay);
    double G_2_down = Compute_EH_from_BD(i-1, j-1, 2, "H", "plus", "y", B, D, params, Deltax, Deltay);
    
    double G_2 = 0.5 * (-G_2_up - G_2_down) - 0.5 * (Dx_iminus1_j_up - Dx_ijminus1_down);


    //F_2{i-1/2,j} calculation (F_1)

    double Dy_ij_up = D[1][i][j] - 0.5 * Deltax * slope_calc(D, 1, i, j, 0, Deltax, Deltay);
    double Dy_iminus1_j_down = D[1][i-1][j] + 0.5 * Deltax * slope_calc(D, 1, i-1, j, 0, Deltax, Deltay);

    double F_1_up = Compute_EH_from_BD(i, j, 2, "H", "minus", "x", B, D, params, Deltax, Deltay);
    double F_1_down = Compute_EH_from_BD(i-1, j, 2, "H", "plus", "x", B, D, params, Deltax, Deltay);
    

    //Note: Ez is negative in the F matrix in Yu paper

    double F_1 = 0.5 * (F_1_up + F_1_down) - 0.5 * (Dy_ij_up - Dy_iminus1_j_down);

    //F_2{i-1/2,j-1} calculation (F_2)

    double Dy_i_jminus1_up = D[1][i][j-1] - 0.5 * Deltax * slope_calc(D, 1, i, j-1, 0, Deltax, Deltay);
    double Dy_ij_minus1_down = D[1][i-1][j-1] + 0.5 * Deltax * slope_calc(D, 1, i-1, j-1, 0, Deltax, Deltay);

    double F_2_up = Compute_EH_from_BD(i, j-1, 2, "H", "minus", "x", B, D, params, Deltax, Deltay);
    double F_2_down = Compute_EH_from_BD(i-1, j-1, 2, "H", "plus", "x", B, D, params, Deltax, Deltay);

    double F_2 = 0.5 * (F_2_up + F_2_down) - 0.5 * (Dy_i_jminus1_up - Dy_ij_minus1_down);

    //Fluxes for Ez{i,j+1} aka Ez_{i-1, j+1/2}
    //G_1{i,j+1/2} calculation (G_3)

    double Dx_i_jplus1_up = D[0][i][j+1] - 0.5 * Deltay * slope_calc(D, 0, i, j+1, 1, Deltax, Deltay);
    double Dx_ij_down = D[0][i][j] + 0.5 * Deltay * slope_calc(D, 0, i, j, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper
    double G_3_up = Compute_EH_from_BD(i, j+1, 2, "H", "minus", "y", B, D, params, Deltax, Deltay);
    double G_3_down = Compute_EH_from_BD(i, j, 2, "H", "plus", "y", B, D, params, Deltax, Deltay);

    double G_3 = 0.5 * (-G_3_up - G_3_down) - 0.5 * (Dx_i_jplus1_up - Dx_ij_down);
    
    
    //G_1{i-1,j+1/2} calculation (G_4)

    double Dx_iminus1_jplus1_up = D[0][i-1][j+1] - 0.5 * Deltay * slope_calc(D, 0, i-1, j+1, 1, Deltax, Deltay);
    double Dx_iminus1_j_down = D[0][i-1][j] + 0.5 * Deltay * slope_calc(D, 0, i-1, j, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper

    double G_4_up = Compute_EH_from_BD(i-1, j+1, 2, "H", "minus", "y", B, D, params, Deltax, Deltay);
    double G_4_down = Compute_EH_from_BD(i-1, j, 2, "H", "plus", "y", B, D, params, Deltax, Deltay);

    double G_4 = 0.5 * (-G_4_up - G_4_down) - 0.5 * (Dx_iminus1_jplus1_up - Dx_iminus1_j_down);

    //F_2{i-1/2,j} calculation (F_3) same thing as F_1

    //Note: Ez is negative in the F matrix in Yu paper and no dissipation term since Bx lives w/ F_1

    double F_3 = 0.5 * (F_1_up + F_1_down) - 0.5 * (Dy_ij_up - Dy_iminus1_j_down);

    //F_2{i-1/2,j+1} calculation (F_4)

    double Dy_i_jplus1_up = D[1][i][j+1] - 0.5 * Deltax * slope_calc(D, 1, i, j+1, 0, Deltax, Deltay);
    double Dy_iminus1_jplus1_down = D[1][i-1][j+1] + 0.5 * Deltax * slope_calc(D, 1, i-1, j+1, 0, Deltax, Deltay);

    double F_4_up = Compute_EH_from_BD(i, j+1, 2, "H", "minus", "x", B, D, params, Deltax, Deltay);
    double F_4_down = Compute_EH_from_BD(i-1, j+1, 2, "H", "plus", "x", B, D, params, Deltax, Deltay);

    double F_4 = 0.5 * (F_4_up + F_4_down) - 0.5 * (Dy_i_jplus1_up - Dy_iminus1_jplus1_down);

    //Toth 2000 Eq(19)
    double Hz_ij = 0.25 * (F_1 + F_2 - G_1 - G_2);
    double Hz_i_jplus1 = 0.25 * (F_3 + F_4 - G_3 - G_4);

    return {Hz_ij, Hz_i_jplus1};
}

/*
    Computes Hz fluxes using slope limiters to update D_y. Calculation is analagous to the Ez flux calculation for the By
    time evolution, hence, copied and pasted the Ez_Flux_Calculation_Qy function.

    Inputs: i, j, B, D, Deltax, and Deltay

    Outputs: Hz{i-1/2,j-1/2} and Hz{i+1/2, j-1/2} fluxes computed using Toth 2000 Eq(19).
*/

std::pair<double, double> Hz_Flux_Calculation_Fy(int i, int j, VectorField & B, VectorField & D, const SimParams & params, double Deltax, double Deltay){

    //Fluxes for Ez{i,j} aka Ez_{i-1/2, j-1/2}
    //G_1{i,j-1/2} calculation (G_1)

    double Dx_ij_up = D[0][i][j] - 0.5 * Deltay * slope_calc(D, 0, i, j, 1, Deltax, Deltay);
    double Dx_i_jminus1_down = D[0][i][j-1] + 0.5 * Deltay * slope_calc(D, 0, i, j-1, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper

    double G_1_up = Compute_EH_from_BD(i, j, 2, "H", "minus", "y", B, D, params, Deltax, Deltay);
    double G_1_down = Compute_EH_from_BD(i, j-1, 2, "H", "plus", "y", B, D, params, Deltax, Deltay);
    
    double G_1 = 0.5 * (-G_1_up - G_1_down) - 0.5 * (Dx_ij_up - Dx_i_jminus1_down);

    //G_1{i-1,j-1/2} calculation (G_2)

    double Dx_iminus1_j_up = D[0][i-1][j] - 0.5 * Deltay * slope_calc(D, 0, i-1, j, 1, Deltax, Deltay);
    double Dx_ijminus1_down = D[0][i-1][j-1] + 0.5 * Deltay * slope_calc(D, 0, i-1, j-1, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper 
    double G_2_up = Compute_EH_from_BD(i-1, j, 2, "H", "minus", "y", B, D, params, Deltax, Deltay);
    double G_2_down = Compute_EH_from_BD(i-1, j-1, 2, "H", "plus", "y", B, D, params, Deltax, Deltay);

    double G_2 = 0.5 * (-G_2_up - G_2_down) - 0.5 * (Dx_iminus1_j_up - Dx_ijminus1_down);

    //F_2{i-1/2,j} calculation (F_1)

    double Dy_ij_up = D[1][i][j] - 0.5 * Deltax * slope_calc(D, 1, i, j, 0, Deltax, Deltay);
    double Dy_iminus1_j_down = D[1][i-1][j] + 0.5 * Deltax * slope_calc(D, 1, i-1, j, 0, Deltax, Deltay);

    double F_1_up = Compute_EH_from_BD(i, j, 2, "H", "minus", "x", B, D, params, Deltax, Deltay);
    double F_1_down = Compute_EH_from_BD(i-1, j, 2, "H", "plus", "x", B, D, params, Deltax, Deltay);

    //Note: Ez is negative in the F matrix in Yu paper

    double F_1 = 0.5 * (F_1_up + F_1_down) - 0.5 * (Dy_ij_up - Dy_iminus1_j_down);

    //F_2{i-1/2,j-1} calculation (F_2)

    double Dy_i_jminus1_up = D[1][i][j-1] - 0.5 * Deltax * slope_calc(D, 1, i, j-1, 0, Deltax, Deltay);
    double Dy_ij_minus1_down = D[1][i-1][j-1] + 0.5 * Deltax * slope_calc(D, 1, i-1, j-1, 0, Deltax, Deltay);

    double F_2_up = Compute_EH_from_BD(i, j-1, 2, "H", "minus", "x", B, D, params, Deltax, Deltay);
    double F_2_down = Compute_EH_from_BD(i-1, j-1, 2, "H", "plus", "x", B, D, params, Deltax, Deltay);

    double F_2 = 0.5 * (F_2_up + F_2_down) - 0.5 * (Dy_i_jminus1_up - Dy_ij_minus1_down);

    //E_z Fluxes for Ez{i+1,j} aka Ez_{i+1/2, j-1/2}

    //G_1{i,j-1/2} calculation (G_3), same as G_1

    //Note: Ez is positive in the G matrix in Yu paper

    double G_3 = 0.5 * (-G_1_up - G_1_down) - 0.5 * (Dx_ij_up - Dx_i_jminus1_down);

    //G_1{i+1, j-1/2}  calculation (G_4)

    double Dx_iplus1_j_up = D[0][i+1][j] - 0.5 * Deltay * slope_calc(D, 0, i+1, j, 1, Deltax, Deltay);
    double Dx_iplus1_jminus1_down = D[0][i+1][j-1] + 0.5 * Deltay * slope_calc(D, 0, i+1, j-1, 1, Deltax, Deltay);

    double G_4_up = Compute_EH_from_BD(i+1, j, 2, "H", "minus", "y", B, D, params, Deltax, Deltay);
    double G_4_down = Compute_EH_from_BD(i+1, j-1, 2, "H", "plus", "y", B, D, params, Deltax, Deltay);

    double G_4 = 0.5 * (-G_4_up - G_4_down) - 0.5 * (Dx_iplus1_j_up - Dx_iplus1_jminus1_down);

    //F_1{i+1/2, j-1} calculation (F_3)

    double Dy_iplus1_jminus1_up = D[1][i+1][j-1] - 0.5 * Deltax * slope_calc(D, 1, i+1, j-1, 0, Deltax, Deltay);
    double Dy_i_jminus1_down = D[1][i][j-1] + 0.5 * Deltax *  slope_calc(D, 1, i, j-1, 0, Deltax, Deltay);

    double F_3_up = Compute_EH_from_BD(i+1, j-1, 2, "H", "minus", "x", B, D, params, Deltax, Deltay);
    double F_3_down = Compute_EH_from_BD(i, j-1, 2, "H", "plus", "x", B, D, params, Deltax, Deltay);

    double F_3 = 0.5 * (F_3_up + F_3_down) - 0.5 * (Dy_iplus1_jminus1_up - Dy_i_jminus1_down);

    //F_1{i+1/2,j} calculation (F_4)

    double Dy_iplus1_j_up = D[1][i+1][j] - 0.5 * Deltax * slope_calc(D, 1, i+1, j, 0, Deltax, Deltay);
    double Dy_ij_down = D[1][i][j] + 0.5 * Deltax * slope_calc(D, 1, i, j, 0, Deltax, Deltay); 

    double F_4_up = Compute_EH_from_BD(i+1, j, 2, "H", "minus", "x", B, D, params, Deltax, Deltay);
    double F_4_down = Compute_EH_from_BD(i, j, 2, "H", "plus", "x", B, D, params, Deltax, Deltay);

    double F_4 = 0.5 * (F_4_up + F_4_down) - 0.5 * (Dy_iplus1_j_up - Dy_ij_down);

    double Hz_ij = 0.25 * (F_1 + F_2 - G_1 - G_2);
    double Hz_iplus1_j = 0.25 * (F_3 + F_4 - G_3 - G_4);

    return {Hz_ij, Hz_iplus1_j};

}

/*
    Computes reconstructed fluxes using slope limiters to update D_z.
    Inputs: i, j, B, D, Deltax, and Deltay

    Outputs: Update for Dz by computing F_3/F_4 (Hy), G_3/G_4 (Hx) fluxes computed using Yu 2011 Eq(9).
*/

double Fz_reconstruction(int i, int j, VectorField & D, VectorField & B, const SimParams & params, double Deltax, double Deltay){

    //G_3{i,j-1/2} calculation G_3

    //Fix compute_eh to be more explicit

    double G_3_up = Compute_EH_from_BD(i, j, 0, "H", "minus", "y", B, D, params, Deltax, Deltay);
    double G_3_down = Compute_EH_from_BD(i, j-1, 0, "H", "plus", "y", B, D, params, Deltax, Deltay);

    double Dz_ij_up_y = D[2][i][j] - 0.5 * Deltay * slope_calc(D, 2, i, j, 1, Deltax, Deltay); 
    double Dz_i_jminus1_down_y = D[2][i][j-1] + 0.5 * Deltay * slope_calc(D, 2, i, j-1, 1, Deltax, Deltay);

    double G_3 = 0.5 * (G_3_up + G_3_down) - 0.5 * (Dz_ij_up_y - Dz_i_jminus1_down_y);

    //G_3{i,j+1/2} calculation G_4

    double G_4_up = Compute_EH_from_BD(i, j+1, 0, "H", "minus", "y", B, D, params, Deltax, Deltay);
    double G_4_down = Compute_EH_from_BD(i, j, 0, "H", "plus", "y", B, D, params, Deltax, Deltay);

    double Dz_i_jplus1_up_y = D[2][i][j+1] - 0.5 * Deltay * slope_calc(D, 2, i, j+1, 1, Deltax, Deltay);
    double Dz_ij_down_y = D[2][i][j] + 0.5 * Deltay * slope_calc(D, 2, i, j, 1, Deltax, Deltay);

    double G_4 = 0.5 * (G_4_up + G_4_down) - 0.5 * (Dz_i_jplus1_up_y - Dz_ij_down_y);

    //F_3{i-1/2,j} calculation (F_3)

    double F_3_up = Compute_EH_from_BD(i, j, 1, "H", "minus", "x", B, D, params, Deltax, Deltay);
    double F_3_down = Compute_EH_from_BD(i-1, j, 1, "H", "plus", "x", B, D, params, Deltax, Deltay);
    
    double Dz_ij_up_x = D[2][i][j] - 0.5 * Deltax * slope_calc(D, 2, i, j, 0, Deltax, Deltay);
    double Dz_iminus1_j_down_x = D[2][i-1][j] + 0.5 * Deltax * slope_calc(D, 2, i-1, j, 0, Deltax, Deltay);

    double F_3 = 0.5 * (-F_3_up - F_3_down) - 0.5 * (Dz_ij_up_x - Dz_iminus1_j_down_x);
    
    //F_3{i+1/2,j} calculation (F_4)

    double F_4_up = Compute_EH_from_BD(i+1, j, 1, "H", "minus", "x", B, D, params, Deltax, Deltay);
    double F_4_down = Compute_EH_from_BD(i, j, 1, "H", "plus", "x", B, D, params, Deltax, Deltay);
    
    double Dz_iplus1_j_up_x = D[2][i+1][j] - 0.5 * Deltax * slope_calc(D, 2, i+1, j, 0, Deltax, Deltay);
    double Dz_ij_down_x = D[2][i][j] + 0.5 * Deltax * slope_calc(D, 2, i, j, 0, Deltax, Deltay); 

    //Changed from F_4_up+F_4_down
    double F_4 = 0.5 * (-F_4_up - F_4_down) - 0.5 * (Dz_iplus1_j_up_x - Dz_ij_down_x);

    //Eq(9) in Yu paper
    //double Fz =  -(F_4 - F_3)/(Deltax) - (G_4 - G_3)/(Deltay); 

    double Fz =  (-F_4 + F_3)/(Deltax) - (G_4 - G_3)/(Deltay); 

    return Fz;
}




/*
    Computes lambda scalar, where lambda = 1 + 8pi*kappa*(E^2-B^2), and thus, D_vec = lambda * E_vec, H_vec = lambda * B_vec
    This solver uses the Newton-Raphson method to solve for lambda. 
    
    Inputs: 
        D_squared: The cell centered average value for the D field
        B_squared: The cell centered average value for the B field
        kappa: From Petri paper
        max_iter: max iterations
        tol: tolerance used to know when solution is found

    Output:
        Lambda: solution used to solve for each component of E

*/

/*
double Solve_lambda_Newton(double D_squared, double B_squared, double kappa, int max_iter, double tol) {
    
    double lambda = 1.0;

    for (int n = 0; n < max_iter; ++n) {
        double lambda_eq  = lambda*lambda*lambda - (1.0 - 8*M_PI*kappa*B_squared) * lambda*lambda - 8*M_PI*kappa*D_squared;
        double lambda_eq_derivative = 3*lambda*lambda - 2*(1.0 - 8*M_PI*kappa*B_squared) * lambda;
        
        // Computes (f(x_n)/f'(x_n))
        double dlambda = lambda_eq / lambda_eq_derivative; 
        
        // Computes (lambda = lambda - dlambda), the next step for the Newton-Raphson method
        lambda -= dlambda;
        if (std::fabs(dlambda) < tol) break;
    }
    return lambda;
}
*/

/*
    Computes lambda scalar, where lambda = 1 + 8pi*kappa*(E^2-B^2).
    This solver uses the gsl poly method to solve the cubic equation. 
    
    Inputs: 
        D_squared: The cell centered average value for the D field
        B_squared: The cell centered average value for the B field

    Output:
        Lambda0: solution used to solve for each component of E and H

*/
double Solve_Lambda_Cubic(double D_squared, double B_squared){

    double kappa = (4*alpha_e)/(45*4.4*B_0);
    double a = 1;
    double b = 8 * M_PI * kappa * B_squared - 1;
    double c = -8 * M_PI * kappa * D_squared;

    //Defines roots locally
    double lambda0{0.0}, lambda1{0.0}, lambda2{0.0};
    
    int num_roots = gsl_poly_solve_cubic (a, b, c, &lambda0, &lambda1, &lambda2);
    
    double physical_lambda = 0.0;

     if (num_roots == 1) {
        // Clamp the single real root to valid physical bounds [0, 1]
        if (lambda0 >= 0.0 && lambda0 <= 1.0) {
            physical_lambda = lambda0;
        } else {
            physical_lambda = 1.0; 
        }
    }
    
       else if (num_roots == 3) {
        double roots[3] = {lambda0, lambda1, lambda2};
        bool found = false;

       
        physical_lambda = -1.0; 

        for (int k = 0; k < 3; k++) {
            if (roots[k] >= 0.0 && roots[k] <= 1.0) {
               
                if (roots[k] > physical_lambda) {
                    physical_lambda = roots[k];
                    found = true;
                }
            }
        }
    }
      else {

        std::cout << "No real root found." << std::endl;
        return 0.0;

     }
return physical_lambda;

}
    
/*

    Computes E using the cell-centered value of B 
    to use in ExB calculation.

*/

double Compute_E_from_D(int i_offset, int j_offset, int component_index, VectorField & B, VectorField & D, const SimParams & params){

    if (params.QED_corrections){

        double D_squared = D[0][i_offset][j_offset] * D[0][i_offset][j_offset] + D[1][i_offset][j_offset] * D[1][i_offset][j_offset] + D[2][i_offset][j_offset] * D[2][i_offset][j_offset];
        double B_squared = B[0][i_offset][j_offset] * B[0][i_offset][j_offset] + B[1][i_offset][j_offset] * B[1][i_offset][j_offset] + B[2][i_offset][j_offset] * B[2][i_offset][j_offset]; 

        double lambda = Solve_Lambda_Cubic(D_squared, B_squared);

        return D[component_index][i_offset][j_offset]/lambda;
    }

    else{

        return D[component_index][i_offset][j_offset];
    }

} 

    


/*
    Background (reference) fields towards which the absorbing layer relaxes B and D.

    Damping the *full* field to zero (as was done previously) destroys the background poloidal
    field inside the layer: Bx -> 0 there, which creates a large div(B), a tiny B^2 in the
    denominator of the force-free current, and a sharp field gradient at the inner edge of the
    layer that continuously radiates grid-scale noise back into the domain. Only the
    perturbation (B - B_background), (D - D_background) should be absorbed.
*/
static VectorField B_background;
static VectorField D_background;
static bool background_set = false;

void Set_Damping_Background(const VectorField & B, const VectorField & D){
    B_background.resize(boost::extents[B.shape()[0]][B.shape()[1]][B.shape()[2]]);
    D_background.resize(boost::extents[D.shape()[0]][D.shape()[1]][D.shape()[2]]);
    B_background = B;
    D_background = D;
    background_set = true;
}

/*
    Computes the lambda damping term in Cerutti et al. 2015 eq(7) where our x-coordinate
    is the r coordinate.

    Inputs:
        x: physical x-coordinate (reduced units) at which the damped quantity lives
           (cell centre for y/z components, x-face for x components)
        Domain object

    Outputs:
        lambda damping term, eq(7) in Cerutti et al. 2015. Zero outside the absorbing layer.

    Notes:
        - dm.x does NOT contain ghost cells, so callers must pass a coordinate, not a ghosted index.
        - With Heun's RK2 the amplification factor of the pure decay term is 1 - s + s^2/2
          (s = lambda*Deltat), which is minimal at s = 1 and returns to 1 (no damping at all) at
          s = 2. K_abs is therefore limited to 1 so the damping is monotonic across the layer.
*/
double Compute_Damping_Term(double x, const Domain & dm){

    //Computes max x-coordinate
    const double max_dm_x = *std::max_element(dm.x.begin(), dm.x.end());
    const double x_abs = 0.9 * max_dm_x; //inner edge of absorbing layer

    //Numerical parameter that controls damping strength (lambda*Deltat at the outer edge)
    const double K_abs = 1.0;

    if (x <= x_abs) return 0.0;

    double xi = (x - x_abs)/(max_dm_x - x_abs);
    if (xi > 1.0) xi = 1.0;

    return (K_abs)/(dm.Deltat) * xi * xi * xi;
}

void Force_Free_Constraint(VectorField & B, VectorField & D, const Domain & dm, const Process & ps){

    const size_t Sx = D.shape()[1];
    const size_t Sy = D.shape()[2];
    const size_t G = dm.N_GC;

    static VectorField D_new(boost::extents[3][Sx][Sy]);
    static ScalarField s(boost::extents[Sx][Sy]);

    const double EPSILON = 1e-14;

    /*
        1. Remove the component of D parallel to B (D.B = 0). This is done on one extra layer of
           faces on each side of the interior so that the |D| <= |B| check below can be evaluated
           with the *projected* D in every cell adjacent to an interior face.
    */
    for(size_t i = G - 1; i <= Sx - G; i++){
      for(size_t j = G - 1; j <= Sy - G; j++){

            //For Dx component (lives on x-face i-1/2, j)
            double Bx = B[0][i][j];
            double By_avg_x = 0.25 * (B[1][i][j] + B[1][i-1][j] + B[1][i][j+1] + B[1][i-1][j+1]);
            double Bz_avg_x = 0.5 * (B[2][i][j] + B[2][i-1][j]);

            double Dy_avg_x = 0.25 * (D[1][i][j] + D[1][i-1][j] + D[1][i][j+1] + D[1][i-1][j+1]);
            double Dz_avg_x = 0.5 * (D[2][i][j] + D[2][i-1][j]);

            double Bx_mag = std::sqrt((Bx * Bx) + (By_avg_x * By_avg_x) + (Bz_avg_x * Bz_avg_x));

            if (Bx_mag > EPSILON) {
                double Bhat_x = Bx / Bx_mag;
                double Bhat_y = By_avg_x / Bx_mag;
                double Bhat_z = Bz_avg_x / Bx_mag;

                double Ddot_Bhat_x = D[0][i][j] * Bhat_x + Dy_avg_x * Bhat_y + Dz_avg_x * Bhat_z;
                D_new[0][i][j] = D[0][i][j] - (Ddot_Bhat_x * Bhat_x);
            } else {
                D_new[0][i][j] = D[0][i][j]; // Leave unchanged if no B-field exists here
            }

            //For Dy component (lives on y-face i, j-1/2)
            double Bx_avg_y = 0.25 * (B[0][i][j] + B[0][i+1][j] + B[0][i][j-1] + B[0][i+1][j-1]);
            double By = B[1][i][j];
            double Bz_avg_y = 0.5 * (B[2][i][j] + B[2][i][j-1]);

            double Dx_avg_y = 0.25 * (D[0][i][j] + D[0][i+1][j] + D[0][i][j-1] + D[0][i+1][j-1]);
            double Dz_avg_y = 0.5 * (D[2][i][j] + D[2][i][j-1]);

            double By_mag = std::sqrt((Bx_avg_y * Bx_avg_y) + (By * By) + (Bz_avg_y * Bz_avg_y));

            if (By_mag > EPSILON) {
                double Bhat_x_y = Bx_avg_y / By_mag;
                double Bhat_y_y = By / By_mag;
                double Bhat_z_y = Bz_avg_y / By_mag;

                double Ddot_Bhat_y = Dx_avg_y * Bhat_x_y + D[1][i][j] * Bhat_y_y + Dz_avg_y * Bhat_z_y;
                D_new[1][i][j] = D[1][i][j] - (Ddot_Bhat_y * Bhat_y_y);
            } else {
                D_new[1][i][j] = D[1][i][j];
            }

            //For Dz component (lives at cell centre i, j)
            double Bx_avg_z = 0.5 * (B[0][i][j] + B[0][i+1][j]);
            double By_avg_z = 0.5 * (B[1][i][j] + B[1][i][j+1]);
            double Bz = B[2][i][j];

            double Dx_avg_z = 0.5 * (D[0][i][j] + D[0][i+1][j]);
            double Dy_avg_z = 0.5 * (D[1][i][j] + D[1][i][j+1]);

            double Bz_mag = std::sqrt((Bx_avg_z * Bx_avg_z) + (By_avg_z * By_avg_z) + (Bz * Bz));

            if (Bz_mag > EPSILON) {
                double Bhat_x_z = Bx_avg_z / Bz_mag;
                double Bhat_y_z = By_avg_z / Bz_mag;
                double Bhat_z_z = Bz / Bz_mag;

                double Ddot_Bhat_z = Dx_avg_z * Bhat_x_z + Dy_avg_z * Bhat_y_z + D[2][i][j] * Bhat_z_z;
                D_new[2][i][j] = D[2][i][j] - (Ddot_Bhat_z * Bhat_z_z);
            } else {
                D_new[2][i][j] = D[2][i][j];
            }
        }
    }

    /*
        2. Cell-centred limiter factor for |D| <= |B|, evaluated with the projected D
           (previously the un-projected D was used, which over-limits whenever D had a
           parallel component).
    */
    for(size_t i = G - 1; i < Sx - G; i++){
      for(size_t j = G - 1; j < Sy - G; j++){

            double Dx = 0.5 * (D_new[0][i][j] + D_new[0][i+1][j]);
            double Dy = 0.5 * (D_new[1][i][j] + D_new[1][i][j+1]);
            double Dz = D_new[2][i][j];
            double Bxc = 0.5 * (B[0][i][j] + B[0][i+1][j]);
            double Byc = 0.5 * (B[1][i][j] + B[1][i][j+1]);
            double Bzc = B[2][i][j];

            double D_squared = Dx*Dx + Dy*Dy + Dz*Dz;
            double B_squared = Bxc*Bxc + Byc*Byc + Bzc*Bzc;

            s[i][j] = (D_squared > B_squared) ? std::sqrt(B_squared / D_squared) : 1.0;
        }
    }

    /*
        3. Apply the limiter. A face shared by two cells takes the stronger of the two factors,
           so the result does not depend on which side of the face the cell index points to.
    */
    for(size_t i = G; i < Sx - G; i++){
        for(size_t j = G; j < Sy - G; j++){
            D[0][i][j] = D_new[0][i][j] * std::min(s[i-1][j], s[i][j]);
            D[1][i][j] = D_new[1][i][j] * std::min(s[i][j-1], s[i][j]);
            D[2][i][j] = D_new[2][i][j] * s[i][j];
        }
    }

    return;
}

/*
    Reconstructs component c of A (cell-centred) to the face of cell (i,j) in direction dir
    (0 = x, 1 = y). side = -1 for the lower face, +1 for the upper face.
*/
static inline double Face_Value(VectorField & A, int c, int i, int j, int dir, double side, double Deltax, double Deltay){
    const double delta = (dir == 0) ? Deltax : Deltay;
    return A[c][i][j] + side * 0.5 * delta * slope_calc(A, c, i, j, dir, Deltax, Deltay);
}

/*
    Corner electric field Ez_{i-1/2,j-1/2} from the Balsara & Spicer (1999) average of the four
    Rusanov face fluxes adjacent to the corner (Toth 2000 Eq. 19). B and D are cell-centred.

    A single value per corner is used for BOTH the Bx and the By update, which is what makes the
    constrained-transport update preserve div(B) to round-off.
*/
double Ez_Corner(int i, int j, VectorField & D, VectorField & B, const SimParams & params, double Deltax, double Deltay){

    //G = y-flux of Bx = +Ez at y-face (ii, j-1/2)
    auto G_face = [&](int ii){
        double Ez_R = Compute_EH_from_BD(ii, j,   2, "E", "minus", "y", B, D, params, Deltax, Deltay);
        double Ez_L = Compute_EH_from_BD(ii, j-1, 2, "E", "plus",  "y", B, D, params, Deltax, Deltay);
        double Bx_R = Face_Value(B, 0, ii, j,   1, -1.0, Deltax, Deltay);
        double Bx_L = Face_Value(B, 0, ii, j-1, 1, +1.0, Deltax, Deltay);
        return 0.5 * (Ez_R + Ez_L) - 0.5 * (Bx_R - Bx_L);
    };

    //F = x-flux of By = -Ez at x-face (i-1/2, jj)
    auto F_face = [&](int jj){
        double Ez_R = Compute_EH_from_BD(i,   jj, 2, "E", "minus", "x", B, D, params, Deltax, Deltay);
        double Ez_L = Compute_EH_from_BD(i-1, jj, 2, "E", "plus",  "x", B, D, params, Deltax, Deltay);
        double By_R = Face_Value(B, 1, i,   jj, 0, -1.0, Deltax, Deltay);
        double By_L = Face_Value(B, 1, i-1, jj, 0, +1.0, Deltax, Deltay);
        return -0.5 * (Ez_R + Ez_L) - 0.5 * (By_R - By_L);
    };

    return 0.25 * (G_face(i) + G_face(i-1) - F_face(j) - F_face(j-1));
}

/*
    Corner auxiliary field Hz_{i-1/2,j-1/2}, analogous to Ez_Corner for the D update.
    Using one value per corner for both Dx and Dy makes the discrete div(D) (and hence Rho)
    change only through the current, i.e. discrete charge conservation.
*/
double Hz_Corner(int i, int j, VectorField & B, VectorField & D, const SimParams & params, double Deltax, double Deltay){

    //G = y-flux of Dx = -Hz at y-face (ii, j-1/2)
    auto G_face = [&](int ii){
        double Hz_R = Compute_EH_from_BD(ii, j,   2, "H", "minus", "y", B, D, params, Deltax, Deltay);
        double Hz_L = Compute_EH_from_BD(ii, j-1, 2, "H", "plus",  "y", B, D, params, Deltax, Deltay);
        double Dx_R = Face_Value(D, 0, ii, j,   1, -1.0, Deltax, Deltay);
        double Dx_L = Face_Value(D, 0, ii, j-1, 1, +1.0, Deltax, Deltay);
        return -0.5 * (Hz_R + Hz_L) - 0.5 * (Dx_R - Dx_L);
    };

    //F = x-flux of Dy = +Hz at x-face (i-1/2, jj)
    auto F_face = [&](int jj){
        double Hz_R = Compute_EH_from_BD(i,   jj, 2, "H", "minus", "x", B, D, params, Deltax, Deltay);
        double Hz_L = Compute_EH_from_BD(i-1, jj, 2, "H", "plus",  "x", B, D, params, Deltax, Deltay);
        double Dy_R = Face_Value(D, 1, i,   jj, 0, -1.0, Deltax, Deltay);
        double Dy_L = Face_Value(D, 1, i-1, jj, 0, +1.0, Deltax, Deltay);
        return 0.5 * (Hz_R + Hz_L) - 0.5 * (Dy_R - Dy_L);
    };

    return 0.25 * (F_face(j) + F_face(j-1) - G_face(i) - G_face(i-1));
}


/*
    This computes the RHS of the time-evolution equations for B and D.

        dB/dt = -curl(E)
        dD/dt =  curl(H) - 4 pi J
        4 pi J = div(D) (E x B)/B^2 + (B.curl(H) - D.curl(E)) B/B^2

    Staggering (index (i,j) of each array):
        Bx, Dx : x-face  (i-1/2, j)
        By, Dy : y-face  (i, j-1/2)
        Bz, Dz : centre  (i, j)
        Ez, Hz : corner  (i-1/2, j-1/2)

    Inputs:
        Qx, Qy, Qz, Fx, Fy, Fz, Rho, D, B, SimParams object, Domain object, Process Object

    Outputs:
        Updated Qx, Qy, and Qz for B time-evolution equations
        Updated Fx, Fy, and Fz for D time-evolution equations
        J (physical current density, J = (4 pi J)/(4 pi)) at the same locations as D
*/
void Compute_RHS(ScalarField & Qx, ScalarField & Qy, ScalarField & Qz, ScalarField & Fx, ScalarField & Fy, ScalarField & Fz, ScalarField & Rho, VectorField & E, VectorField & H, VectorField & D, VectorField & B, VectorField & J, const SimParams & params, const Domain & dm, const Process & ps, double t)
{
    const size_t Sx = B.shape()[1];
    const size_t Sy = B.shape()[2];
    const size_t G = dm.N_GC;
    const double Deltay = dm.Deltay;

    static VectorField B_avg(boost::extents[3][Sx][Sy]);
    static VectorField D_avg(boost::extents[3][Sx][Sy]);
    static ScalarField Ez_c(boost::extents[Sx][Sy]);   //corner Ez
    static ScalarField Hz_c(boost::extents[Sx][Sy]);   //corner Hz
    static ScalarField Qz_c(boost::extents[Sx][Sy]);   //-curl(E)_z at centres
    static ScalarField Fz_c(boost::extents[Sx][Sy]);   // curl(H)_z at centres
    static VectorField J_c(boost::extents[3][Sx][Sy]); //cell-centred 4 pi J

    /*
        Cell-centred B and D over every cell for which the face average exists.
        (Previously the last row/column was left at zero but still reached by the
        reconstruction stencils at the top/right edge of every MPI sub-domain.)
    */
    for (size_t i = 0; i < Sx; i++){
        for (size_t j = 0; j < Sy; j++){
            D_avg[0][i][j] = (i + 1 < Sx) ? 0.5 * (D[0][i][j] + D[0][i+1][j]) : D[0][i][j];
            D_avg[1][i][j] = (j + 1 < Sy) ? 0.5 * (D[1][i][j] + D[1][i][j+1]) : D[1][i][j];
            D_avg[2][i][j] = D[2][i][j];

            B_avg[0][i][j] = (i + 1 < Sx) ? 0.5 * (B[0][i][j] + B[0][i+1][j]) : B[0][i][j];
            B_avg[1][i][j] = (j + 1 < Sy) ? 0.5 * (B[1][i][j] + B[1][i][j+1]) : B[1][i][j];
            B_avg[2][i][j] = B[2][i][j];
        }
    }

    /*
        Corner EMFs, computed once per corner. The stencil of corner (i,j) reaches cells
        i-2..i+1, j-2..j+1, which stays inside [0, S-2] for corners [G-1, S-G].
    */
    for (size_t i = G - 1; i <= Sx - G; i++){
        for (size_t j = G - 1; j <= Sy - G; j++){
            Ez_c[i][j] = Ez_Corner(i, j, D_avg, B_avg, params, dm.Deltax[i], Deltay);
            Hz_c[i][j] = Hz_Corner(i, j, B_avg, D_avg, params, dm.Deltax[i], Deltay);
        }
    }

    /*
        Cell-centred 4 pi J on cells [G-1, S-G-1]: one ghost layer below the interior is needed
        so that the current can be averaged onto the first interior x- and y-faces.
    */
    for (size_t i = G - 1; i < Sx - G; i++){
        for (size_t j = G - 1; j < Sy - G; j++){

            const double Deltax = dm.Deltax[i];

            //Face-normal curls on the four faces of cell (i,j)
            const double curlH_x_L = (Hz_c[i][j+1]   - Hz_c[i][j])   / Deltay;   //(curl H)_x at (i-1/2, j)
            const double curlH_x_R = (Hz_c[i+1][j+1] - Hz_c[i+1][j]) / Deltay;   //(curl H)_x at (i+1/2, j)
            const double curlH_y_B = -(Hz_c[i+1][j]   - Hz_c[i][j])   / Deltax;  //(curl H)_y at (i, j-1/2)
            const double curlH_y_T = -(Hz_c[i+1][j+1] - Hz_c[i][j+1]) / Deltax;  //(curl H)_y at (i, j+1/2)

            const double curlE_x_L = (Ez_c[i][j+1]   - Ez_c[i][j])   / Deltay;
            const double curlE_x_R = (Ez_c[i+1][j+1] - Ez_c[i+1][j]) / Deltay;
            const double curlE_y_B = -(Ez_c[i+1][j]   - Ez_c[i][j])   / Deltax;
            const double curlE_y_T = -(Ez_c[i+1][j+1] - Ez_c[i][j+1]) / Deltax;

            Fz_c[i][j] = Fz_reconstruction(i, j, D_avg, B_avg, params, Deltax, Deltay); // (curl H)_z
            Qz_c[i][j] = Qz_reconstruction(i, j, D_avg, B_avg, params, Deltax, Deltay); //-(curl E)_z

            const double curlH_z = Fz_c[i][j];
            const double curlE_z = -Qz_c[i][j];

            //B.curl(H) and D.curl(E): each face-normal curl is multiplied by the face-normal
            //field that lives on the same face, then averaged to the cell centre.
            const double B_dot_curlH = 0.5 * (B[0][i][j] * curlH_x_L + B[0][i+1][j] * curlH_x_R)
                                     + 0.5 * (B[1][i][j] * curlH_y_B + B[1][i][j+1] * curlH_y_T)
                                     + B[2][i][j] * curlH_z;
            const double D_dot_curlE = 0.5 * (D[0][i][j] * curlE_x_L + D[0][i+1][j] * curlE_x_R)
                                     + 0.5 * (D[1][i][j] * curlE_y_B + D[1][i][j+1] * curlE_y_T)
                                     + D[2][i][j] * curlE_z;

            const double Bxc = B_avg[0][i][j], Byc = B_avg[1][i][j], Bzc = B_avg[2][i][j];
            const double Exc = Compute_E_from_D(i, j, 0, B_avg, D_avg, params);
            const double Eyc = Compute_E_from_D(i, j, 1, B_avg, D_avg, params);
            const double Ezc = Compute_E_from_D(i, j, 2, B_avg, D_avg, params);

            const double B_squared = Bxc*Bxc + Byc*Byc + Bzc*Bzc;

            if (B_squared > 0.0){
                const double divD = Rho[i][j];
                const double par = B_dot_curlH - D_dot_curlE;

                J_c[0][i][j] = (divD * (Eyc*Bzc - Ezc*Byc) + par * Bxc) / B_squared;
                J_c[1][i][j] = (divD * (Ezc*Bxc - Exc*Bzc) + par * Byc) / B_squared;
                J_c[2][i][j] = (divD * (Exc*Byc - Eyc*Bxc) + par * Bzc) / B_squared;
            }
            else{
                J_c[0][i][j] = J_c[1][i][j] = J_c[2][i][j] = 0.0;
            }
        }
    }

    for (size_t i = G; i < Sx - G; i++){
        for (size_t j = G; j < Sy - G; j++){

            const double Deltax = dm.Deltax[i];

            const double curlE_x = (Ez_c[i][j+1] - Ez_c[i][j]) / Deltay;   //at x-face (i-1/2, j)
            const double curlE_y = -(Ez_c[i+1][j] - Ez_c[i][j]) / Deltax;  //at y-face (i, j-1/2)
            const double curlH_x = (Hz_c[i][j+1] - Hz_c[i][j]) / Deltay;
            const double curlH_y = -(Hz_c[i+1][j] - Hz_c[i][j]) / Deltax;

            //Current averaged onto the location of each D component
            const double Jx = 0.5 * (J_c[0][i-1][j] + J_c[0][i][j]); //x-face (i-1/2, j)
            const double Jy = 0.5 * (J_c[1][i][j-1] + J_c[1][i][j]); //y-face (i, j-1/2)
            const double Jz = J_c[2][i][j];                             //centre

            Qx[i][j] = -curlE_x;
            Qy[i][j] = -curlE_y;
            Qz[i][j] = Qz_c[i][j];

            Fx[i][j] = curlH_x - Jx;
            Fy[i][j] = curlH_y - Jy;
            Fz[i][j] = Fz_c[i][j] - Jz;

            J[0][i][j] = Jx;
            J[1][i][j] = Jy;
            J[2][i][j] = Jz;

            //Absorbing layer: relax the perturbation towards the background state
            const size_t ip = i - G; //index into dm.x, which has no ghost cells
            if (ip < dm.x.size()){
                const double x_c = dm.x[ip];
                const double x_f = x_c - 0.5 * Deltax;
                const double sig_c = Compute_Damping_Term(x_c, dm);
                const double sig_f = Compute_Damping_Term(x_f, dm);

                if (sig_c > 0.0 || sig_f > 0.0){
                    const double Bb0 = background_set ? B_background[0][i][j] : 0.0;
                    const double Bb1 = background_set ? B_background[1][i][j] : 0.0;
                    const double Bb2 = background_set ? B_background[2][i][j] : 0.0;
                    const double Db0 = background_set ? D_background[0][i][j] : 0.0;
                    const double Db1 = background_set ? D_background[1][i][j] : 0.0;
                    const double Db2 = background_set ? D_background[2][i][j] : 0.0;

                    Qx[i][j] -= sig_f * (B[0][i][j] - Bb0);
                    Qy[i][j] -= sig_c * (B[1][i][j] - Bb1);
                    Qz[i][j] -= sig_c * (B[2][i][j] - Bb2);

                    Fx[i][j] -= sig_f * (D[0][i][j] - Db0);
                    Fy[i][j] -= sig_c * (D[1][i][j] - Db1);
                    Fz[i][j] -= sig_c * (D[2][i][j] - Db2);
                }
            }
        }
    }

    return;
}
