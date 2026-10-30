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
   

    for(size_t i=dm.N_GC; i<D.shape()[1]-dm.N_GC; i++){
      for(size_t j=dm.N_GC; j<D.shape()[2]-dm.N_GC; j++){

       //for (size_t i = dm.N_GC - 2; i < D.shape()[1] - dm.N_GC + 2 ; i++){
        //for (size_t j = dm.N_GC - 2; j < D.shape()[2] - dm.N_GC + 2; j++){

            //Rho[i][j] = (D[0][i+1][j] - D[0][i][j])/(4.*pi*dm.Deltax[i]) + (D[1][i][j+1] - D[1][i][j])/(4.*pi*dm.Deltay);
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
    
    double Bx_iminus1_j_up = Compute_A_to_cell_center(B, 0, i-1, j) - 0.5 * Deltay * slope_calc(B, 0, i-1, j, 1, Deltax, Deltay);
    double Bx_ijminus1_down = Compute_A_to_cell_center(B, 0, i-1, j-1) + 0.5 * Deltay * slope_calc(B, 0, i-1, j-1, 1, Deltax, Deltay);

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
    Computes the lambda damping term in Cerutti et al. 2015 eq(7) where our x-coordinate
    is the r coordinate.

    Inputs:
        index i, Domain object
    
    Outputs:
        lambda damping term, eq(7) in Cerutti et al. 2015

*/

double Compute_Damping_Term(int i, const Domain & dm){

    //Computes max x-coordinate
    double max_dm_x = *std::max_element(dm.x.begin(), dm.x.end());

    //Numerical parameter that controls damping strength
    double K_abs = 20.0;

    double damping_term = (K_abs)/(dm.Deltat) * ((dm.x[i] - 0.9 * max_dm_x)/(max_dm_x - 0.9 * max_dm_x)) * ((dm.x[i] - 0.9 * max_dm_x)/(max_dm_x - 0.9 * max_dm_x)) * ((dm.x[i] - 0.9 * max_dm_x)/(max_dm_x - 0.9 * max_dm_x));
    
    return damping_term;
}

void Force_Free_Constraint(VectorField & B, VectorField & D, const Domain & dm, const Process & ps){

   
    //VectorField D_new(boost::extents[3][dm.Nx+2*dm.N_GC][ps.MyE-ps.MyS+2*dm.N_GC]); 

    //VectorField D_new = D;

    static VectorField D_new(boost::extents[3][D.shape()[1]][D.shape()[2]]);


    const double EPSILON = 1e-14;

    for(size_t i = dm.N_GC; i < D.shape()[1] - dm.N_GC; i++){
      for(size_t j = dm.N_GC; j < D.shape()[2] - dm.N_GC; j++){
       //for (size_t i = dm.N_GC - 2; i < D.shape()[1] - dm.N_GC + 2 ; i++){
        //for (size_t j = dm.N_GC - 2; j < D.shape()[2] - dm.N_GC + 2; j++){

            //For Dx component
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
                //D[0][i][j] = D[0][i][j] - (Ddot_Bhat_x * Bhat_x);

            } else {
                

                D_new[0][i][j] = D[0][i][j]; // Leave unchanged if no B-field exists here
            }

            //For Dy component
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

                //D[1][i][j] = D[1][i][j] - (Ddot_Bhat_y * Bhat_y_y);
                
            } else {
                D_new[1][i][j] = D[1][i][j];
                
            }

            //For Dz component
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

                //D[2][i][j] = D[2][i][j] - (Ddot_Bhat_z * Bhat_z_z);
                
            } else {
                D_new[2][i][j] = D[2][i][j];
            }

        }
    }

    for(size_t i = dm.N_GC; i < D.shape()[1] - dm.N_GC; i++){
      for(size_t j = dm.N_GC; j < D.shape()[2] - dm.N_GC; j++){

            double Dx = 0.5 * (D[0][i][j] + D[0][i+1][j]);
            double Dy = 0.5 * (D[1][i][j] + D[1][i][j+1]);
            double Dz = D[2][i][j];
            double Bxc = 0.5 * (B[0][i][j] + B[0][i+1][j]);
            double Byc = 0.5 * (B[1][i][j] + B[1][i][j+1]);
            double Bzc = B[2][i][j];

            double D_squared = Dx*Dx + Dy*Dy + Dz*Dz;
            double B_squared = Bxc*Bxc + Byc*Byc + Bzc*Bzc;

            //If force-free condition is not satisfied.
            if (D_squared > B_squared) {
                double s = std::sqrt(B_squared / D_squared);
                D_new[0][i][j] *= s;
                D_new[1][i][j] *= s;
                D_new[2][i][j] *= s;
            }
        }
    }

    for(size_t c = 0; c < 3; c++){
    for(size_t i = dm.N_GC; i < D.shape()[1] - dm.N_GC; i++){
        for(size_t j = dm.N_GC; j < D.shape()[2] - dm.N_GC; j++){
            D[c][i][j] = D_new[c][i][j];
        }
    }
}
        

    //D = D_new;
    return;
}



/*
    This computes the RHS of the time-evolution equations for B and D. 

    Inputs:
        Qx, Qy, Qz, Fx, Fy, Fz, Rho, D, B, SimParams object, Domain object, Process Object

    Outputs:
        Updated Qx, Qy, and Qz for B time-evolution equations
        Updated Fx, Fy, and Fz for D time-evolution equations
*/
void Compute_RHS(ScalarField & Qx, ScalarField & Qy, ScalarField & Qz, ScalarField & Fx, ScalarField & Fy, ScalarField & Fz, ScalarField & Rho, VectorField & E, VectorField & H, VectorField & D, VectorField & B, VectorField & J, const SimParams & params, const Domain & dm, const Process & ps, double t)
{

    //Compute_EH_from_DB(E, H, D, B, params, dm);
    //exchng2Vector(E, dm.N_GC, ps.comm1D, ps.nbrleft, ps.nbrright);
    //exchng2Vector(H, dm.N_GC, ps.comm1D, ps.nbrleft, ps.nbrright);

    VectorField B_avg(boost::extents[3][dm.Nx+2*dm.N_GC][ps.MyE-ps.MyS+2*dm.N_GC]);
    VectorField D_avg(boost::extents[3][dm.Nx+2*dm.N_GC][ps.MyE-ps.MyS+2*dm.N_GC]);
    
    for (size_t i = dm.N_GC - 2; i < B.shape()[1] - dm.N_GC + 2 ; i++){
        for (size_t j = dm.N_GC - 2; j < B.shape()[2] - dm.N_GC + 2; j++){

            D_avg[0][i][j] = Compute_A_to_cell_center(D, 0, i, j);
            D_avg[1][i][j] = Compute_A_to_cell_center(D, 1, i, j);
            D_avg[2][i][j] = D[2][i][j];

            B_avg[0][i][j] = Compute_A_to_cell_center(B, 0, i, j);
            B_avg[1][i][j] = Compute_A_to_cell_center(B, 1, i, j);
            B_avg[2][i][j] = B[2][i][j];



    }
}

    

    double max_dm_x = *std::max_element(dm.x.begin(), dm.x.end());

    for(size_t i=dm.N_GC; i<D.shape()[1]-dm.N_GC; i++){
        for(size_t j=dm.N_GC; j<D.shape()[2]-dm.N_GC; j++){

            
            
            double Bx_avg_ij = B_avg[0][i][j];
            double Bx_avg_iplus1_j = B_avg[0][i+1][j];
            double Bx_avg_i_jplus1 = B_avg[0][i][j+1];
            double Bx_avg_i_jplus2 = B_avg[0][i][j+2];
            double Bx_avg_iplus2_j = B_avg[0][i+2][j];
            double Bx_avg_iplus1_jplus1 = B_avg[0][i+1][j+1];
            

            double By_avg_ij = B_avg[1][i][j];
            double By_avg_iplus1_j = B_avg[1][i+1][j];
            double By_avg_i_jplus1 = B_avg[1][i][j+1];
            double By_avg_i_jplus2 = B_avg[1][i][j+2];
            double By_avg_iplus1_jplus1 = B_avg[1][i+1][j+1];

            double Dx_avg_ij = D_avg[0][i][j]; 
            double Dx_avg_iplus1_j = D_avg[0][i+1][j];
            double Dx_avg_iplus2_j = D_avg[0][i+2][j];
            double Dx_avg_i_jplus1 = D_avg[0][i][j+1];
            double Dx_avg_i_jplus2 = D_avg[0][i][j+2];
            double Dx_avg_iplus1_jplus1 = D_avg[0][i+1][j+1];

            double Dy_avg_ij = D_avg[1][i][j];
            double Dy_avg_i_jplus1 = D_avg[1][i][j+1];
            double Dy_avg_i_jplus2 = D_avg[1][i][j+2];
            double Dy_avg_iplus1_j = D_avg[1][i+1][j];
            double Dy_avg_iplus1_jplus1 = D_avg[1][i+1][j+1];

   
            //Used for calculating (curl(H) dot B - curl(E) dot D)_ij 
            
            auto [Ez_Qx_left_bottom, Ez_Qx_left_top] = Ez_Flux_Calculation_Qx(i, j, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay); //For i-1/2, j face
            auto [Ez_Qx_right_bottom, Ez_Qx_right_top] = Ez_Flux_Calculation_Qx(i+1, j, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay); //For i+1/2, j face
            
            auto [Ez_Qy_left_bottom, Ez_Qy_right_bottom] = Ez_Flux_Calculation_Qy(i, j, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay); //For i, j-1/2 face
            auto[Ez_Qy_left_top, Ez_Qy_right_top] = Ez_Flux_Calculation_Qy(i, j+1, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay); //For i, j+1/2 face
            
            auto [Hz_Fx_left_bottom, Hz_Fx_left_top] = Hz_Flux_Calculation_Fx(i, j, B_avg, D_avg, params, dm.Deltax[i], dm.Deltay); //For i-1/2, j face
            auto [Hz_Fx_right_bottom, Hz_Fx_right_top] = Hz_Flux_Calculation_Fx(i+1, j, B_avg, D_avg, params, dm.Deltax[i], dm.Deltay); //For i+1/2, j face

            auto [Hz_Fy_left_bottom, Hz_Fy_right_bottom] = Hz_Flux_Calculation_Fy(i, j, B_avg, D_avg, params, dm.Deltax[i], dm.Deltay); //For i, j-1/2 face
            auto[Hz_Fy_left_top, Hz_Fy_right_top] = Hz_Flux_Calculation_Fy(i, j+1, B_avg, D_avg, params, dm.Deltax[i], dm.Deltay); //For i, j+1/2 face

            //Used for calculating (curl(H) dot B - curl(E) dot D)_iplus1_j
            
            auto [Ez_Qx_right2_bottom, Ez_Qx_right2_top] = Ez_Flux_Calculation_Qx(i+2, j, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay); // For i+3/2, j face
            auto [Hz_Qx_right2_bottom, Hz_Qx_right2_top] = Hz_Flux_Calculation_Fx(i+2, j, B_avg, D_avg, params, dm.Deltax[i], dm.Deltay); // For i+3/2, j face

            auto[Ez_Qy_left2_bottom, Ez_Qy_right2_bottom] = Ez_Flux_Calculation_Qy(i+1, j, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay); //For i+1/2, j-1/2 face
            auto[Hz_Fy_left2_bottom, Hz_Fy_right2_bottom] = Hz_Flux_Calculation_Fy(i+1, j, B_avg, D_avg, params, dm.Deltax[i], dm.Deltay); //For i+1/2, j-1/2 face

            auto[Ez_Qy_left2_top, Ez_Qy_right2_top] = Ez_Flux_Calculation_Qy(i+1, j+1, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay); //For i+1/2, j+1/2 face
            auto[Hz_Fy_left2_top, Hz_Fy_right2_top] = Hz_Flux_Calculation_Fy(i+1, j+1, B_avg, D_avg, params, dm.Deltax[i], dm.Deltay); //For i+1/2, j+1/2 face

            //Used for calculating (curl(H) dot B - curl(E) dot D)_i_jplus1

            auto [Ez_Qx_left3_bottom, Ez_Qx_left3_top] = Ez_Flux_Calculation_Qx(i, j+1, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay); //For i-1/2, j+1
            auto [Hz_Fx_left3_bottom, Hz_Fx_left3_top] = Hz_Flux_Calculation_Fx(i, j+1, B_avg, D_avg, params, dm.Deltax[i], dm.Deltay); //For i-1/2, j+1

            auto [Ez_Qx_right3_bottom, Ez_Qx_right3_top] = Ez_Flux_Calculation_Qx(i+1, j+1, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay); //For i+1/2, j+1
            auto [Hz_Fx_right3_bottom, Hz_Fx_right3_top] = Hz_Flux_Calculation_Fx(i+1, j+1, B_avg, D_avg, params, dm.Deltax[i], dm.Deltay); //For i+1/2, j+1

            auto[Ez_Qy_left3_bottom, Ez_Qy_right3_bottom] = Ez_Flux_Calculation_Qy(i, j+1, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay); //For i, j+1/2 face
            auto[Hz_Fy_left3_bottom, Hz_Fy_right3_bottom] = Hz_Flux_Calculation_Fy(i, j+1, B_avg, D_avg, params, dm.Deltax[i], dm.Deltay); //For i, j+1/2 face

            auto[Ez_Qy_left3_top, Ez_Qy_right3_top] = Ez_Flux_Calculation_Qy(i, j+2, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay); //For i, j+3/2 face
            auto[Hz_Fy_left3_top, Hz_Fy_right3_top] = Hz_Flux_Calculation_Fy(i, j+2, B_avg, D_avg, params, dm.Deltax[i], dm.Deltay); //For i, j+3/2 face
            
            //ExB calculations

            //double E_cross_B_x_ij = E[1][i][j] * B[2][i][j] - E[2][i][j] * By_avg_ij;
            //double E_cross_B_x_iplus1_j =  E[1][i+1][j] * B[2][i+1][j] - E[2][i+1][j] * By_avg_iplus1_j;

            //Which Ex, Ey, and Ez do I use? Use compute_E_from_D or what is returned from Flux functions?
            double E_cross_B_x_ij = Compute_E_from_D(i, j, 1, B_avg, D_avg, params) * B[2][i][j] - Compute_E_from_D(i, j, 2, B_avg, D_avg, params) * By_avg_ij;
            double E_cross_B_x_iplus1_j =  Compute_E_from_D(i+1, j, 1, B_avg, D_avg, params) * B[2][i+1][j] - Compute_E_from_D(i+1, j, 2, B_avg, D_avg, params) * By_avg_iplus1_j;



            //double E_cross_B_y_ij = E[2][i][j] * Bx_avg_ij - E[0][i][j] * B_avg[2][i][j];
            //double E_cross_B_y_i_jplus1 =  E[2][i][j+1] * Bx_avg_i_jplus1 - E[0][i][j+1] * B_avg[2][i][j+1];

            //double E_cross_B_z_ij = E[0][i][j] * By_avg_ij - E[1][i][j] * Bx_avg_ij;

            double E_cross_B_y_ij = Compute_E_from_D(i, j, 2, B_avg, D_avg, params) * Bx_avg_ij - Compute_E_from_D(i, j, 0, B_avg, D_avg, params) * B_avg[2][i][j];
            double E_cross_B_y_i_jplus1 =  Compute_E_from_D(i, j+1, 2, B_avg, D_avg, params) * Bx_avg_i_jplus1 - Compute_E_from_D(i, j+1, 0, B_avg, D_avg, params) * B_avg[2][i][j+1];

            double E_cross_B_z_ij = Compute_E_from_D(i, j, 0, B_avg, D_avg, params) * By_avg_ij -  Compute_E_from_D(i, j, 1, B_avg, D_avg, params) * Bx_avg_ij;

            //B^2 calculations

            double B_squared_ij = B[2][i][j] * B[2][i][j] + Bx_avg_ij*Bx_avg_ij + By_avg_ij*By_avg_ij;
            double B_squared_iplus1_j = B[2][i+1][j]*B[2][i+1][j] + Bx_avg_iplus1_j*Bx_avg_iplus1_j + By_avg_iplus1_j*By_avg_iplus1_j;
            double B_squared_i_jplus1 = B[2][i][j+1]*B[2][i][j+1] + Bx_avg_i_jplus1*Bx_avg_i_jplus1 + By_avg_i_jplus1*By_avg_i_jplus1;


            //(curl(H) dot B - curl(E) dot D)_ij calculation
            double curl_H_x_left_ij = (Hz_Fx_left_top - Hz_Fx_left_bottom)/(dm.Deltay); //For i-1/2,j face
            double curl_H_x_right_ij = (Hz_Fx_right_top - Hz_Fx_right_bottom)/(dm.Deltay); //For i+1/2,j face

            double curl_E_x_left_ij = (Ez_Qx_left_top - Ez_Qx_left_bottom)/(dm.Deltay); //For i-1/2,j face
            double curl_E_x_right_ij = (Ez_Qx_right_top - Ez_Qx_right_bottom)/(dm.Deltay); //For i+1/2,j face

            double curl_H_y_bottom_ij = (Hz_Fy_right_bottom - Hz_Fy_left_bottom)/(dm.Deltax[i]);  //For i, j-1/2 face
            double curl_H_y_top_ij = (Hz_Fy_right_top - Hz_Fy_left_top)/(dm.Deltax[i]); //For i, j+1/2 face

            auto[Hz_Fy_left_top_iplus1_jplus1, Hz_Fy_right_top_iplus1_jplus1] = Hz_Flux_Calculation_Fy(i+1, j+1, B_avg, D_avg, params, dm.Deltax[i], dm.Deltay); //For i+1/2, j+1/2 face
            auto[Ez_Qy_left_top_iplus1_jplus1, Ez_Qy_right_top_iplus1_jplus1] = Ez_Flux_Calculation_Qy(i+1, j+1, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay); //For i+1/2, j+1/2 face

            double curl_H_y_top_iplus1_jplus1 = (Hz_Fy_right_top_iplus1_jplus1 - Hz_Fy_left_top_iplus1_jplus1)/(dm.Deltax[i]);

            double curl_E_y_top_iplus1_jplus1 = (Ez_Qy_right_top_iplus1_jplus1 - Ez_Qy_left_top_iplus1_jplus1)/(dm.Deltax[i]);

            
            

            double curl_E_y_bottom_ij = (Ez_Qy_right_bottom - Ez_Qy_left_bottom)/(dm.Deltax[i]);  //For i, j-1/2 face
            double curl_E_y_top_ij = (Ez_Qy_right_top - Ez_Qy_left_top)/(dm.Deltax[i]); //For i, j+1/2 face

            double curl_H_z_ij = Fz_reconstruction(i, j, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay);
            double curl_E_z_ij = Qz_reconstruction(i, j, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay);


            double curl_H_curl_E_ij = 0.5 * (curl_H_x_left_ij * Bx_avg_ij + curl_H_x_right_ij * Bx_avg_iplus1_j) + 0.5 * (curl_H_y_bottom_ij * By_avg_ij + curl_H_y_top_ij * By_avg_i_jplus1) + curl_H_z_ij * B[2][i][j] - 0.5 * (curl_E_x_left_ij * Dx_avg_ij + curl_E_x_right_ij * Dx_avg_iplus1_j) - 0.5 * (curl_E_y_bottom_ij * Dy_avg_ij + curl_E_y_top_ij * Dy_avg_i_jplus1) - (curl_E_z_ij * D[2][i][j]);

            //(curl(H) dot B - curl(E) dot D)_iplus1_j calculation

            double curl_H_x_left_iplus1_j = (Hz_Fx_right_top - Hz_Fx_right_bottom)/(dm.Deltay); //For i+1/2,j face
            double curl_H_x_right_iplus1_j = (Hz_Qx_right2_top - Hz_Qx_right2_bottom)/(dm.Deltay); //For i+3/2,j face

            double curl_E_x_left_iplus1_j = (Ez_Qx_right_top - Ez_Qx_right_bottom)/(dm.Deltay); //For i+1/2,j face
            double curl_E_x_right_iplus1_j =  (Ez_Qx_right2_top - Ez_Qx_right2_bottom)/(dm.Deltay); //For i+3/2,j face

            double curl_H_y_bottom_iplus1_j = (Hz_Fy_right2_bottom - Hz_Fy_left2_bottom)/(dm.Deltax[i]);  //For i+1, j-1/2 face
            double curl_H_y_top_iplus1_j = (Hz_Fy_right2_top - Hz_Fy_left2_top)/(dm.Deltax[i]); //For i+1, j+1/2 face

            double curl_E_y_bottom_iplus1_j = (Ez_Qy_right2_bottom - Ez_Qy_left2_bottom)/(dm.Deltax[i]);  //For i+1, j-1/2 face
            double curl_E_y_top_iplus1_j = (Ez_Qy_right2_top - Ez_Qy_left2_top)/(dm.Deltax[i]); //For i+1, j+1/2 face

            double curl_H_z_iplus1_j = Fz_reconstruction(i+1, j, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay);
            double curl_E_z_iplus1_j = Qz_reconstruction(i+1, j, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay);

            double curl_H_curl_E_iplus1_j = 0.5 * (curl_H_x_left_iplus1_j * Bx_avg_iplus1_j + curl_H_x_right_iplus1_j * Bx_avg_iplus2_j) + 0.5 * (curl_H_y_bottom_iplus1_j * By_avg_iplus1_j + curl_H_y_top_iplus1_j * By_avg_iplus1_jplus1) + curl_H_z_iplus1_j * B[2][i+1][j] - 0.5 * (curl_E_x_left_iplus1_j * Dx_avg_iplus1_j + curl_E_x_right_iplus1_j * Dx_avg_iplus2_j) - 0.5 * (curl_E_y_bottom_iplus1_j * Dy_avg_iplus1_j + curl_E_y_top_iplus1_j * Dy_avg_iplus1_jplus1) - (curl_E_z_iplus1_j * D[2][i+1][j]);

            //double curl_H_curl_E_iplus1_j = 0.5 * (curl_H_x_left_iplus1_j * Bx_avg_iplus1_j + curl_H_x_right_iplus1_j * Bx_avg_iplus2_j) + 0.5 * (curl_H_y_bottom_iplus1_j * By_avg_iplus1_j + curl_H_y_top_iplus1_jplus1 * By_avg_iplus1_jplus1) + curl_H_z_iplus1_j * B[2][i+1][j] - 0.5 * (curl_E_x_left_iplus1_j * Dx_avg_iplus1_j + curl_E_x_right_iplus1_j * Dx_avg_iplus2_j) - 0.5 * (curl_E_y_bottom_iplus1_j * Dy_avg_iplus1_j + curl_E_y_top_iplus1_jplus1 * Dy_avg_iplus1_jplus1) - (curl_E_z_iplus1_j * D[2][i+1][j]);


            //(curl(H) dot B - curl(E) dot D)_i_jplus1 calculation

            double curl_H_x_left_i_jplus1 = (Hz_Fx_left3_top - Hz_Fx_left3_bottom)/(dm.Deltay); //For i-1/2,j+1 face
            double curl_H_x_right_i_jplus1 = (Hz_Fx_right3_top - Hz_Fx_right3_bottom)/(dm.Deltay); //For i+1/2,j+1 face

            double curl_E_x_left_i_jplus1 = (Ez_Qx_left3_top - Ez_Qx_left3_bottom)/(dm.Deltay); //For i-1/2,j+1 face
            double curl_E_x_right_i_jplus1 =  (Ez_Qx_right3_top - Ez_Qx_right3_bottom)/(dm.Deltay); //For i+1/2,j+1 face

            double curl_H_y_bottom_i_jplus1 = (Hz_Fy_right3_bottom - Hz_Fy_left3_bottom)/(dm.Deltax[i]);  //For i, j+1/2 face
            double curl_H_y_top_i_jplus1 = (Hz_Fy_right3_top - Hz_Fy_left3_top)/(dm.Deltax[i]); //For i, j+3/2 face

            double curl_E_y_bottom_i_jplus1 = (Ez_Qy_right3_bottom - Ez_Qy_left3_bottom)/(dm.Deltax[i]);  //For i, j+1/2 face
            double curl_E_y_top_i_jplus1 = (Ez_Qy_right3_top - Ez_Qy_left3_top)/(dm.Deltax[i]); //For i, j+3/2 face

            double curl_H_z_i_jplus1 = Fz_reconstruction(i, j+1, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay);
            double curl_E_z_i_jplus1 = Qz_reconstruction(i, j+1, D_avg, B_avg, params, dm.Deltax[i], dm.Deltay);

            double curl_H_curl_E_i_jplus1 = 0.5 * (curl_H_x_left_i_jplus1 * Bx_avg_i_jplus1 + curl_H_x_right_i_jplus1 * Bx_avg_iplus1_jplus1) + 0.5 * (curl_H_y_bottom_i_jplus1 * By_avg_i_jplus1 + curl_H_y_top_i_jplus1 * By_avg_i_jplus2) + curl_H_z_i_jplus1 * B[2][i][j+1] - 0.5 * (curl_E_x_left_i_jplus1 * Dx_avg_i_jplus1 + curl_E_x_right_i_jplus1 * Dx_avg_iplus1_jplus1) - 0.5 * (curl_E_y_bottom_i_jplus1 * Dy_avg_i_jplus1 + curl_E_y_top_i_jplus1 * Dy_avg_i_jplus2) - curl_E_z_i_jplus1 * D[2][i][j+1];

            //Jx calculation

            double Jx_ij = (Rho[i][j] * E_cross_B_x_ij)/(B_squared_ij) + (curl_H_curl_E_ij * Bx_avg_ij)/(B_squared_ij);
            double Jx_iplus1_j = (Rho[i+1][j] * E_cross_B_x_iplus1_j)/(B_squared_iplus1_j) + (curl_H_curl_E_iplus1_j * Bx_avg_iplus1_j)/(B_squared_iplus1_j);

            double Jx_avg = 0.5 * (Jx_ij + Jx_iplus1_j);

            //Jy calculation

            double Jy_ij = (Rho[i][j] * E_cross_B_y_ij)/(B_squared_ij) + (curl_H_curl_E_ij * By_avg_ij)/(B_squared_ij);
            double Jy_i_jplus1 = (Rho[i][j+1] * E_cross_B_y_i_jplus1)/(B_squared_i_jplus1) + (curl_H_curl_E_i_jplus1 * By_avg_i_jplus1)/(B_squared_i_jplus1);

            double Jy_avg = 0.5 * (Jy_ij + Jy_i_jplus1);

            //Jz calculation

            double Jz_ij = (Rho[i][j] * E_cross_B_z_ij)/(B_squared_ij) + (curl_H_curl_E_ij * B[2][i][j])/(B_squared_ij);

            //std::cout << "Rho[i][j] = " << Rho[i][j] << std::endl;
            if (Jx_avg != 0 || Jy_avg != 0 || Jz_ij !=0){
                //std::cout << "Jx_avg= " << Jx_avg << std::endl;
                //std::cout << "Jy_avg= " << Jy_avg << std::endl;
                //std::cout << "Jz_avg= " << Jz_ij << std::endl;
                //std::cout << "t = " << t << std::endl;


            }
            
            Qx[i][j] = -curl_E_x_left_ij;
            Qy[i][j] = curl_E_y_bottom_ij;
            Qz[i][j] = curl_E_z_ij;

            Fx[i][j] = curl_H_x_left_ij - Jx_avg;
            Fy[i][j] = -curl_H_y_bottom_ij - Jy_avg;
            Fz[i][j] = curl_H_z_ij - Jz_ij;

            //J[0][i][j] = Jx_avg;
            //J[1][i][j] = Jy_avg;
            //J[2][i][j] = Jz_ij;

        }
    }

return;

} 

            /*
            if (dm.x[i] >= 0.9 * max_dm_x){
            
            Qx[i][j] = -(Compute_Damping_Term(i, dm) * B[0][i][j]) - curl_E_x_left_ij;
            Qy[i][j] = -(Compute_Damping_Term(i, dm) * B[1][i][j]) + curl_E_y_bottom_ij;
            Qz[i][j] = -(Compute_Damping_Term(i, dm) * B[2][i][j]) + curl_E_z_ij;

            Fx[i][j] = -(Compute_Damping_Term(i, dm) * D[0][i][j]) + curl_H_x_left_ij - Jx_avg;
            Fy[i][j] = -(Compute_Damping_Term(i, dm) * D[1][i][j]) - curl_H_y_bottom_ij - Jy_avg;
            Fz[i][j] = -(Compute_Damping_Term(i, dm) * D[2][i][j]) + curl_H_z_ij - Jz_ij;

            }

            else{

            Qx[i][j] = -curl_E_x_left_ij;
            Qy[i][j] = curl_E_y_bottom_ij;
            Qz[i][j] = curl_E_z_ij;

            Fx[i][j] = curl_H_x_left_ij - Jx_avg;
            Fy[i][j] = -curl_H_y_bottom_ij - Jy_avg;
            Fz[i][j] = curl_H_z_ij - Jz_ij;
            }
            */
            

            
            

           


            
          
    
  

