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

#include "common.h"
#include "initial_conditions.h"
#include "field_evolution.h"

void Compute_Rho(VectorField & D, ScalarField & Rho, const Domain & dm)
{
    size_t N_GC = dm.N_GC; // in domain 

    for(size_t i=N_GC; i<D.shape()[1]-N_GC; i++){
        for(size_t j=N_GC; j<D.shape()[2]-N_GC; j++){

            Rho[i][j] = (D[0][i+1][j] - D[0][i][j])/(4.*pi*dm.Deltax[i]) + (D[1][i][j+1] - D[1][i][j])/(4.*pi*dm.Deltay);
        }
    }
    return;
}

struct Fields
{
    VectorField & E;
    VectorField & B;
    VectorField & H;
    VectorField & D;
    ScalarField & Rho;
    const Domain & dm;  
};

double compute_A1_x(int i, int j, const Fields & f){

    double rho_avg = 0.5 * (f.Rho[i][j] + f.Rho[i][j-1]);
    
    double Ey = 0.25 * (f.E[1][i][j] + f.E[1][i+1][j] + f.E[1][i][j-1] + f.E[1][i+1][j-1]);   

    double Ez = 0.5 * (f.E[2][i][j] + f.E[2][i+1][j]);

    double Bx_i = 0.5 * (f.B[0][i+1][j] * f.B[0][i+1][j]) + 0.5 *(f.B[0][i][j] * f.B[0][i][j]);
    double Bx_i_minus_1 = 0.5 * (f.B[0][i+1][j-1] * f.B[0][i+1][j-1]) + 0.5 * (f.B[0][i][j-1] * f.B[0][i][j-1]);
    double By_j = 0.5 * (f.B[1][i][j+1] * f.B[1][i][j+1]) + 0.5 * (f.B[1][i][j] * f.B[1][i][j]);
    double By_j_minus_1 = 0.5 * (f.B[1][i][j-1] * f.B[1][i][j-1]) + 0.5*(f.B[1][i][j] * f.B[1][i][j]); 
    double Bz = f.B[2][i][j];
    double Bz_j_minus_1 = f.B[2][i][j-1];

    
    double B_ij = (Bx_i)  + (By_j)  + (Bz * Bz); 
    double B_ij_minus_1 = (Bx_i_minus_1) +(By_j_minus_1) + (Bz_j_minus_1 * Bz_j_minus_1);
    double B_squared = 0.5 * (B_ij + B_ij_minus_1);      
    double Bz_A_1_x = 0.5 * (Bz + Bz_j_minus_1);

    double ExB_x = Ey * Bz_A_1_x  - Ez * f.B[1][i][j];

    return (rho_avg * ExB_x)/(B_squared);
}

double compute_A1_y(int i, int j, const Fields & f){
    // Calculation of A_1_y

    double rho_avg = 0.5 * (f.Rho[i][j] + f.Rho[i-1][j]);
    
    double Ex = 0.25 * (f.E[0][i][j] + f.E[0][i][j+1] + f.E[0][i-1][j] + f.E[0][i-1][j+1]);   

    double Ez = 0.5 * (f.E[2][i][j] + f.E[2][i][j+1]);

    double Bx_i = 0.5 * (f.B[0][i+1][j] + f.B[0][i][j]);
    double Bx_i_minus_1_j = 0.5 * (f.B[0][i-1][j] + f.B[0][i][j]);
    double By_j = 0.5 * (f.B[1][i][j+1] + f.B[1][i][j]);
    double By_i_minus_1_j = 0.5 * (f.B[1][i-1][j+1] + f.B[1][i-1][j]); 
    double Bz = f.B[2][i][j];                          
    double Bz_i_minus_1_j = f.B[2][i-1][j];

    double B_ij = (Bx_i * Bx_i)  + (By_j * By_j)  + (Bz * Bz); 
    double B_i_minus_1_j = (Bx_i_minus_1_j * Bx_i_minus_1_j) +(By_i_minus_1_j * By_i_minus_1_j) + (Bz_i_minus_1_j * Bz_i_minus_1_j);
    double B_squared = 0.5 * (B_ij + B_i_minus_1_j);
            
    double Bz_A_1_y = 0.5 * (Bz + Bz_i_minus_1_j);    
                                       

    double ExB_y = Ez * f.B[0][i][j]  - Ex * Bz_A_1_y;

    return (rho_avg * ExB_y)/(B_squared);
}

double compute_A1_z(int i, int j, const Fields & f){
    // Calculation of (A_1)_z

    double rho_avg = 0.25 * (f.Rho[i][j] + f.Rho[i][j-1] + f.Rho[i-1][j] + f.Rho[i-1][j-1]);

    double Ex_By_A = 0.5 * (f.E[0][i][j] * f.B[1][i][j] + f.E[0][i-1][j] * f.B[1][i-1][j]); // at i-1/2, j-1/2
    double Ey_Bx_A = 0.5 * (f.E[1][i][j] * f.B[0][i][j] + f.E[1][i][j-1] * f.B[0][i][j-1]); // at i-1/2, j-1/2

    double Bx_i = 0.5 * (f.B[0][i+1][j] + f.B[0][i][j]);
    double Bx_i_minus_1_j = 0.5 * (f.B[0][i][j] + f.B[0][i-1][j]);
    double Bx_i_j_minus_1 = 0.5 * (f.B[0][i+1][j-1] + f.B[0][i][j-1]);
    double Bx_ij_minus_1 = 0.5 * (f.B[0][i][j-1] + f.B[0][i-1][j-1]);

    double By_j = 0.5 * (f.B[1][i][j+1] + f.B[1][i][j]);
    double By_j_minus_1 = 0.5 * (f.B[1][i][j-1] + f.B[1][i][j]); 
    double By_i_minus_1_j = 0.5 * (f.B[1][i-1][j+1] + f.B[1][i-1][j]);
    double By_ij_minus_1 = 0.5 * (f.B[1][i-1][j] + f.B[1][i-1][j-1]);  

    double Bz = f.B[2][i][j];
    double Bz_i_minus_1_j = f.B[2][i-1][j];
    double Bz_ij_minus_1 = f.B[2][i-1][j-1];
    double Bz_i_j_minus_1 = f.B[2][i][j-1];

    double B_ij = 0.5*(f.B[0][i][j]*f.B[0][i][j] + f.B[0][i+1][j]*f.B[0][i+1][j]) + 0.5*(f.B[1][i][j]*f.B[1][i][j] + f.B[1][i][j+1]*f.B[1][i][j+1]) + f.B[2][i][j]*f.B[2][i][j];
    double B_i_minus_1_j = 0.5*(f.B[0][i-1][j]*f.B[0][i-1][j] + f.B[0][i][j]*f.B[0][i][j]) + 0.5*(f.B[1][i-1][j]*f.B[1][i-1][j] + f.B[1][i-1][j+1]*f.B[1][i-1][j+1]) + f.B[2][i-1][j]*f.B[2][i-1][j];
    double B_i_j_minus_1 = 0.5*(f.B[0][i][j-1]*f.B[0][i][j-1] + f.B[0][i+1][j-1]*f.B[0][i+1][j-1]) + 0.5*(f.B[1][i][j-1]*f.B[1][i][j-1] + f.B[1][i][j]*f.B[1][i][j]) + f.B[2][i][j-1]*f.B[2][i][j-1];
    double B_ij_minus_1 = 0.5*(f.B[0][i-1][j-1]*f.B[0][i-1][j-1] + f.B[0][i][j-1]*f.B[0][i][j-1]) + 0.5*(f.B[1][i-1][j-1]*f.B[1][i-1][j-1] + f.B[1][i-1][j]*f.B[1][i-1][j]) + f.B[2][i-1][j-1]*f.B[2][i-1][j-1];

    double B_squared = (B_ij + B_i_minus_1_j + B_i_j_minus_1 + B_ij_minus_1) / 4.0;

    return (rho_avg * (Ex_By_A - Ey_Bx_A))/(B_squared);
}

double compute_A2_x(int i, int j, const Fields & f){
        
    // Calculation of (A_2)_x

    // Calculates B dot (curl(H)) at i,j

    double Bx_Hz_A = f.B[0][i][j] * (f.H[2][i][j+1] - f.H[2][i][j])/(f.dm.Deltay); // at i - 1/2, j
    double Bx_Hz_B = f.B[0][i+1][j] * (f.H[2][i+1][j+1] - f.H[2][i+1][j])/(f.dm.Deltay); // at i + 1/2, j

    double Bx_Hz_AB = 0.5 * (Bx_Hz_A + Bx_Hz_B); // at i , j

    double By_Hz_A = f.B[1][i][j+1] * (f.H[2][i+1][j+1] - f.H[2][i][j+1])/(f.dm.Deltax[i]); // at i, j + 1/2
    double By_Hz_B = f.B[1][i][j] * (f.H[2][i+1][j] - f.H[2][i][j])/(f.dm.Deltax[i]); // at i, j - 1/2

    double By_Hz_AB = 0.5 * (By_Hz_A + By_Hz_B); // at i,j

    double Bz_Hy_Hx_AB = f.B[2][i][j] * ((f.H[1][i+1][j] - f.H[1][i][j])/(f.dm.Deltax[i]) + (f.H[0][i][j] - f.H[0][i][j+1])/(f.dm.Deltay)); // at i,j

    double B_dot_curl_H_AB = Bx_Hz_AB - By_Hz_AB + Bz_Hy_Hx_AB; // at i,j

    // Calculates B dot (curl(H)) at i , j-1

    double Bx_Hz_C = f.B[0][i+1][j-1] * (f.H[2][i+1][j] - f.H[2][i+1][j-1])/(f.dm.Deltay); // at i + 1/2, j - 1
    double Bx_Hz_D = f.B[0][i][j-1] * (f.H[2][i][j] - f.H[2][i][j-1])/(f.dm.Deltay); // at i - 1/2, j-1

    double Bx_Hz_CD = 0.5 * (Bx_Hz_C + Bx_Hz_D); // at i , j - 1

    double By_Hz_C = f.B[1][i][j] * (f.H[2][i+1][j] - f.H[2][i][j])/(f.dm.Deltax[i]); // at i, j - 1/2
    double By_Hz_D = f.B[1][i][j-1] * (f.H[2][i+1][j-1] - f.H[2][i][j-1])/(f.dm.Deltax[i]); // at i , j - 3/2

    double By_Hz_CD = 0.5 * (By_Hz_C + By_Hz_D); // at i, j - 1

    double Bz_Hy_Hx_CD = f.B[2][i][j-1] * ((f.H[1][i+1][j-1] - f.H[1][i][j-1])/(f.dm.Deltax[i]) + (f.H[0][i][j-1] - f.H[0][i][j])/(f.dm.Deltay)); // at i,j-1; fixed this on 06/23/26

    double B_dot_curl_H_CD = Bx_Hz_CD - By_Hz_CD + Bz_Hy_Hx_CD; // at i,j-1

    double Bx_i = 0.5 * (f.B[0][i+1][j] + f.B[0][i][j]);
    double Bx_i_minus_1 = 0.5 * (f.B[0][i+1][j-1] + f.B[0][i][j-1]);
    double Bx_i_squared = 0.5 * (f.B[0][i+1][j] * f.B[0][i+1][j]) + 0.5 * (f.B[0][i][j] * f.B[0][i][j]);
    double Bx_i_minus_1_squared = 0.5 * (f.B[0][i+1][j-1]*f.B[0][i+1][j-1]) + 0.5 * (f.B[0][i][j-1] * f.B[0][i][j-1]);
    double By_j_squared = 0.5 * (f.B[1][i][j+1] * f.B[1][i][j+1]) + 0.5 * (f.B[1][i][j] * f.B[1][i][j]);
    double By_j_minus_1_squared = 0.5 * (f.B[1][i][j-1] * f.B[1][i][j-1] ) + 0.5 * (f.B[1][i][j] * f.B[1][i][j]); 
    double Bz = f.B[2][i][j];
    double Bz_j_minus_1 = f.B[2][i][j-1];

    double B_ij = (Bx_i_squared)  + (By_j_squared)  + (Bz * Bz);
    double B_ij_minus_1 = (Bx_i_minus_1_squared) +(By_j_minus_1_squared) + (Bz_j_minus_1 * Bz_j_minus_1);

    return 0.5 * ((B_dot_curl_H_AB * Bx_i)/(B_ij) + (B_dot_curl_H_CD * Bx_i_minus_1)/(B_ij_minus_1));
}

double compute_A2_y(int i, int j, const Fields & f){

    // Calculation of A_2_y

    // Calculates B dot (curl(H)) at i,j 

    double Bx_Hz_A = f.B[0][i][j] * (f.H[2][i][j+1] - f.H[2][i][j])/(f.dm.Deltay); // at i - 1/2, j
    double Bx_Hz_B = f.B[0][i+1][j] * (f.H[2][i+1][j+1] - f.H[2][i+1][j])/(f.dm.Deltay); // at i + 1/2, j

    double Bx_Hz_AB = 0.5 * (Bx_Hz_A + Bx_Hz_B); // at i , j

    double By_Hz_A = f.B[1][i][j+1] * (f.H[2][i+1][j+1] - f.H[2][i][j+1])/(f.dm.Deltax[i]); // at i, j + 1/2
    double By_Hz_B = f.B[1][i][j] * (f.H[2][i+1][j] - f.H[2][i][j])/(f.dm.Deltax[i]); // at i, j - 1/2

    double By_Hz_AB = 0.5 * (By_Hz_A + By_Hz_B); // at i,j

    double Bz_Hy_Hx_AB = f.B[2][i][j] * ((f.H[1][i+1][j] - f.H[1][i][j])/(f.dm.Deltax[i]) + (f.H[0][i][j] - f.H[0][i][j+1])/(f.dm.Deltay)); // at i,j

    double B_dot_curl_H_AB = Bx_Hz_AB - By_Hz_AB + Bz_Hy_Hx_AB; // at i,j

    // Calculates B dot (curl(H)) at i-1 , j

    double Bx_Hz_C = f.B[0][i][j] * (f.H[2][i][j+1] - f.H[2][i][j])/(f.dm.Deltay); // at i - 1/2, j 
    double Bx_Hz_D = f.B[0][i-1][j] * (f.H[2][i-1][j+1] - f.H[2][i-1][j])/(f.dm.Deltay); // at i - 3/2, j chaged on 06/23/26

    double Bx_Hz_CD = 0.5 * (Bx_Hz_C + Bx_Hz_D); // at i-1 , j

    double By_Hz_C = f.B[1][i-1][j+1] * (f.H[2][i][j+1] - f.H[2][i-1][j+1])/(f.dm.Deltax[i]); // at i-1, j + 1/2
    double By_Hz_D = f.B[1][i-1][j] * (f.H[2][i][j] - f.H[2][i-1][j])/(f.dm.Deltax[i]); // at i-1, j - 1/2

    double By_Hz_CD = 0.5 * (By_Hz_C + By_Hz_D); // at i-1, j

    double Bz_Hy_Hx_CD = f.B[2][i-1][j] * ((f.H[1][i][j] - f.H[1][i-1][j])/(f.dm.Deltax[i]) + (f.H[0][i-1][j] - f.H[0][i-1][j+1])/(f.dm.Deltay)); 

    double B_dot_curl_H_CD = Bx_Hz_CD - By_Hz_CD + Bz_Hy_Hx_CD; // at i-1,j

    double By_j = 0.5 * (f.B[1][i][j+1] + f.B[1][i][j]);
    double By_i_minus_1_j = 0.5 * (f.B[1][i-1][j+1] + f.B[1][i-1][j]); 
    
    double Bx_i_squared = 0.5 * (f.B[0][i+1][j] * f.B[0][i+1][j]) + 0.5 * (f.B[0][i][j] * f.B[0][i][j]);
    double Bx_i_minus_1_j_squared = 0.5 * (f.B[0][i][j] * f.B[0][i][j]) + 0.5 *(f.B[0][i-1][j] * f.B[0][i-1][j]);
    double By_j_squared = 0.5 * (f.B[1][i][j+1] * f.B[1][i][j+1]) + 0.5 * (f.B[1][i][j] * f.B[1][i][j]);
    double By_i_minus_1_j_squared = 0.5 * (f.B[1][i-1][j+1] * f.B[1][i-1][j+1]) + 0.5 * (f.B[1][i-1][j] * f.B[1][i-1][j]); 
    double Bz = f.B[2][i][j];
    double Bz_i_minus_1_j = f.B[2][i-1][j];

    double B_ij = (Bx_i_squared)  + (By_j_squared)  + (Bz * Bz); 
    double B_i_minus_1_j = (Bx_i_minus_1_j_squared) +(By_i_minus_1_j_squared) + (Bz_i_minus_1_j * Bz_i_minus_1_j);

    return 0.5 * ((B_dot_curl_H_AB * By_j)/(B_ij) + (B_dot_curl_H_CD * By_i_minus_1_j)/(B_i_minus_1_j));
}

double compute_A2_z(int i, int j, const Fields & f){

    // Calculation of A_2_z

    // Calculates B dot (curl(H)) at i,j 

    double Bx_Hz_A = f.B[0][i][j] * (f.H[2][i][j+1] - f.H[2][i][j])/(f.dm.Deltay); // at i - 1/2, j
    double Bx_Hz_B = f.B[0][i+1][j] * (f.H[2][i+1][j+1] - f.H[2][i+1][j])/(f.dm.Deltay); // at i + 1/2, j

    double Bx_Hz_AB = 0.5 * (Bx_Hz_A + Bx_Hz_B); // at i , j

    double By_Hz_A = f.B[1][i][j+1] * (f.H[2][i+1][j+1] - f.H[2][i][j+1])/(f.dm.Deltax[i]); // at i, j + 1/2
    double By_Hz_B = f.B[1][i][j] * (f.H[2][i+1][j] - f.H[2][i][j])/(f.dm.Deltax[i]); // at i, j - 1/2

    double By_Hz_AB = 0.5 * (By_Hz_A + By_Hz_B); // at i,j

    double Bz_Hy_Hx_AB = f.B[2][i][j] * ((f.H[1][i+1][j] - f.H[1][i][j])/(f.dm.Deltax[i]) + (f.H[0][i][j] - f.H[0][i][j+1])/(f.dm.Deltay)); // at i,j

    double B_dot_curl_H_AB = Bx_Hz_AB - By_Hz_AB + Bz_Hy_Hx_AB; // at i,j

    // Calculates B dot (curl(H)) at i-1 , j

    double Bx_Hz_C = f.B[0][i][j] * (f.H[2][i][j+1] - f.H[2][i][j])/(f.dm.Deltay); // at i - 1/2, j 
    double Bx_Hz_D = f.B[0][i-1][j] * (f.H[2][i-1][j+1] - f.H[2][i-1][j])/(f.dm.Deltay); // at i - 3/2, j, changed on 06/24/26

    double Bx_Hz_CD = 0.5 * (Bx_Hz_C + Bx_Hz_D); // at i-1 , j

    double By_Hz_C = f.B[1][i-1][j+1] * (f.H[2][i][j+1] - f.H[2][i-1][j+1])/(f.dm.Deltax[i]); // at i-1, j + 1/2
    double By_Hz_D = f.B[1][i-1][j] * (f.H[2][i][j] - f.H[2][i-1][j])/(f.dm.Deltax[i]); // at i-1, j - 1/2

    double By_Hz_CD = 0.5 * (By_Hz_C + By_Hz_D); // at i-1, j

    double Bz_Hy_Hx_CD = f.B[2][i-1][j] * ((f.H[1][i][j] - f.H[1][i-1][j])/(f.dm.Deltax[i]) + (f.H[0][i-1][j] - f.H[0][i-1][j+1])/(f.dm.Deltay)); 

    double B_dot_curl_H_CD = Bx_Hz_CD - By_Hz_CD + Bz_Hy_Hx_CD; // at i-1,j

    // Calculates B dot (curl(H)) at i-1, j-1

    double Bx_Hz_E = f.B[0][i][j-1] * (f.H[2][i][j] - f.H[2][i][j-1])/(f.dm.Deltay); // at i - 1/2, j-1 
    double Bx_Hz_F = f.B[0][i-1][j-1] * (f.H[2][i-1][j] - f.H[2][i-1][j-1])/(f.dm.Deltay); // at i - 3/2, j-1

    double Bx_Hz_EF = 0.5 * (Bx_Hz_E + Bx_Hz_F); // at i-1 , j-1

    double By_Hz_E = f.B[1][i-1][j] * (f.H[2][i][j] - f.H[2][i-1][j])/(f.dm.Deltax[i]); // at i-1, j - 1/2
    double By_Hz_F = f.B[1][i-1][j-1] * (f.H[2][i][j-1] - f.H[2][i-1][j-1])/(f.dm.Deltax[i]); // at i-1, j - 3/2

    double By_Hz_EF = 0.5 * (By_Hz_E + By_Hz_F); // at i-1, j-1

    double Bz_Hy_Hx_EF = f.B[2][i-1][j-1] * ((f.H[1][i][j-1] - f.H[1][i-1][j-1])/(f.dm.Deltax[i]) + (f.H[0][i-1][j-1] - f.H[0][i-1][j])/(f.dm.Deltay)); 

    double B_dot_curl_H_EF = Bx_Hz_EF - By_Hz_EF + Bz_Hy_Hx_EF; // at i-1,j-1

    // Calculates B dot (curl(H)) at i, j-1

    double Bx_Hz_G = f.B[0][i+1][j-1] * (f.H[2][i+1][j] - f.H[2][i+1][j-1])/(f.dm.Deltay); // at i + 1/2, j - 1
    double Bx_Hz_H = f.B[0][i][j-1] * (f.H[2][i][j] - f.H[2][i][j-1])/(f.dm.Deltay); // at i - 1/2, j-1

    double Bx_Hz_GH = 0.5 * (Bx_Hz_G + Bx_Hz_H); 

    double By_Hz_G = f.B[1][i][j] * (f.H[2][i+1][j] - f.H[2][i][j])/(f.dm.Deltax[i]); // at i, j - 1/2
    double By_Hz_H = f.B[1][i][j-1] * (f.H[2][i+1][j-1] - f.H[2][i][j-1])/(f.dm.Deltax[i]); // at i, j - 3/2

    double By_Hz_GH = 0.5 * (By_Hz_G + By_Hz_H); 

    double Bz_Hy_Hx_GH = f.B[2][i][j-1] * ((f.H[1][i+1][j-1] - f.H[1][i][j-1])/(f.dm.Deltax[i]) + (f.H[0][i][j-1] - f.H[0][i][j])/(f.dm.Deltay)); // at i,j-1 changed on 06/24/26

    double B_dot_curl_H_GH = Bx_Hz_GH - By_Hz_GH + Bz_Hy_Hx_GH; 

    double Bx_i = 0.5 * (f.B[0][i+1][j] * f.B[0][i+1][j]) + 0.5 * (f.B[0][i][j] * f.B[0][i][j]);
    double Bx_i_minus_1_j = 0.5 * (f.B[0][i][j] * f.B[0][i][j]) + 0.5 * (f.B[0][i-1][j] * f.B[0][i-1][j]);
    double Bx_i_j_minus_1 = 0.5 * (f.B[0][i+1][j-1] * f.B[0][i+1][j-1]) + 0.5 * (f.B[0][i][j-1] * f.B[0][i][j-1]);
    double Bx_ij_minus_1 = 0.5 * (f.B[0][i][j-1] * f.B[0][i][j-1]) + 0.5 * (f.B[0][i-1][j-1] * f.B[0][i-1][j-1]);

    double By_j = 0.5 * (f.B[1][i][j+1] * f.B[1][i][j+1]) + 0.5 * (f.B[1][i][j] * f.B[1][i][j]);
    double By_j_minus_1 = 0.5 * (f.B[1][i][j-1] * f.B[1][i][j-1]) + 0.5 * (f.B[1][i][j] * f.B[1][i][j]); 
    double By_i_minus_1_j = 0.5 * (f.B[1][i-1][j+1] * f.B[1][i-1][j+1]) + 0.5 * (f.B[1][i-1][j] * f.B[1][i-1][j]);
    double By_ij_minus_1 = 0.5 * (f.B[1][i-1][j] * f.B[1][i-1][j]) + 0.5 * (f.B[1][i-1][j-1] * f.B[1][i-1][j-1]);  

    double Bz = f.B[2][i][j];
    double Bz_i_minus_1_j = f.B[2][i-1][j];
    double Bz_ij_minus_1 = f.B[2][i-1][j-1];
    double Bz_i_j_minus_1 = f.B[2][i][j-1];

    double B_ij = (Bx_i)  + (By_j)  + (Bz * Bz); 
    double B_i_minus_1_j = (Bx_i_minus_1_j) +(By_i_minus_1_j) + (Bz_i_minus_1_j * Bz_i_minus_1_j);
    double B_ij_minus_1 = (Bx_ij_minus_1) + (By_ij_minus_1) + (Bz_ij_minus_1 * Bz_ij_minus_1);
    double B_i_j_minus_1 = (Bx_i_j_minus_1) + (By_j_minus_1) + (Bz_i_j_minus_1 * Bz_i_j_minus_1);     
    
    return 0.25 * ((B_dot_curl_H_AB * Bz)/(B_ij) + (B_dot_curl_H_CD * Bz_i_minus_1_j)/(B_i_minus_1_j) + (B_dot_curl_H_EF * Bz_i_j_minus_1)/(B_i_j_minus_1) + (B_dot_curl_H_GH * Bz_ij_minus_1)/(B_ij_minus_1));
}


double compute_A3_x(int i, int j, const Fields & f){

    // Calculation of A_3_x
            
    // Calculates -D dot (curl(E)) at i,j

    double Dx_Ez_A = f.D[0][i][j] * (f.E[2][i][j] - f.E[2][i][j+1])/(f.dm.Deltay); // at i - 1/2, j
    double Dx_Ez_B = f.D[0][i+1][j] * (f.E[2][i+1][j] - f.E[2][i+1][j+1])/(f.dm.Deltay); 

    double Dx_Ez_AB = 0.5 * (Dx_Ez_A + Dx_Ez_B); // at i , j

    double Dy_Ez_A = f.D[1][i][j+1] * (f.E[2][i+1][j+1] - f.E[2][i][j+1])/(f.dm.Deltax[i]); // at i, j + 1/2
    double Dy_Ez_B = f.D[1][i][j] * (f.E[2][i+1][j] - f.E[2][i][j])/(f.dm.Deltax[i]); // at i, j - 1/2

    double Dy_Ez_AB = 0.5 * (Dy_Ez_A + Dy_Ez_B); // FIX: was Dy_Hz_A + Dy_Hz_B

    double Dz_Ey_Ex_AB = f.D[2][i][j] * ((f.E[0][i][j+1] - f.E[0][i][j])/(f.dm.Deltay) + (f.E[1][i][j] - f.E[1][i+1][j])/(f.dm.Deltax[i])); 

    double D_dot_curl_E_AB = Dx_Ez_AB + Dy_Ez_AB + Dz_Ey_Ex_AB; // at i,j

    // Calculates -D dot (curl(E)) at i , j-1

    double Dx_Ez_C = f.D[0][i+1][j-1] * (f.E[2][i+1][j-1] - f.E[2][i+1][j])/(f.dm.Deltay); // at i + 1/2, j - 1
    double Dx_Ez_D = f.D[0][i][j-1] * (f.E[2][i][j-1] - f.E[2][i][j])/(f.dm.Deltay); // at i - 1/2, j-1

    double Dx_Ez_CD = 0.5 * (Dx_Ez_C + Dx_Ez_D); // at i , j - 1

    double Dy_Ez_C = f.D[1][i][j] * (f.E[2][i+1][j] - f.E[2][i][j])/(f.dm.Deltax[i]); // at i, j - 1/2
    double Dy_Ez_D = f.D[1][i][j-1] * (f.E[2][i+1][j-1] - f.E[2][i][j-1])/(f.dm.Deltax[i]); // at i , j - 3/2

    double Dy_Ez_CD = 0.5 * (Dy_Ez_C + Dy_Ez_D); // at i, j - 1

    double Dz_Ey_Ex_CD = f.D[2][i][j-1] * ((f.E[0][i][j] - f.E[0][i][j-1])/(f.dm.Deltay) + (f.E[1][i][j-1] - f.E[1][i+1][j-1])/(f.dm.Deltax[i])); 

    double D_dot_curl_E_CD = Dx_Ez_CD + Dy_Ez_CD + Dz_Ey_Ex_CD; // at i,j-1

    double Bx_i = 0.5 * (f.B[0][i+1][j] + f.B[0][i][j]);
    double Bx_i_minus_1 = 0.5 * (f.B[0][i+1][j-1] + f.B[0][i][j-1]);
    double Bx_i_squared = 0.5 * (f.B[0][i+1][j] * f.B[0][i+1][j]) + 0.5 * (f.B[0][i][j] * f.B[0][i][j]);
    double Bx_i_minus_1_squared = 0.5 * (f.B[0][i+1][j-1]*f.B[0][i+1][j-1]) + 0.5 * (f.B[0][i][j-1] * f.B[0][i][j-1]);
    double By_j_squared = 0.5 * (f.B[1][i][j+1] * f.B[1][i][j+1]) + 0.5 * (f.B[1][i][j] * f.B[1][i][j]);
    double By_j_minus_1_squared = 0.5 * (f.B[1][i][j-1] * f.B[1][i][j-1] ) + 0.5 * (f.B[1][i][j] * f.B[1][i][j]); 
    double Bz = f.B[2][i][j];
    double Bz_j_minus_1 = f.B[2][i][j-1];

    double B_ij = (Bx_i_squared)  + (By_j_squared)  + (Bz * Bz);
    double B_ij_minus_1 = (Bx_i_minus_1_squared) +(By_j_minus_1_squared) + (Bz_j_minus_1 * Bz_j_minus_1);
    
    return 0.5 * ((D_dot_curl_E_AB * Bx_i)/(B_ij) + (D_dot_curl_E_CD * Bx_i_minus_1)/(B_ij_minus_1));
}

double compute_A3_y(int i, int j, const Fields & f){
    
    // Calculation of A_3_y

    // Calculates -D dot (curl(E)) at i,j

    double Dx_Ez_A = f.D[0][i][j] * (f.E[2][i][j] - f.E[2][i][j+1])/(f.dm.Deltay); // at i - 1/2, j
    double Dx_Ez_B = f.D[0][i+1][j] * (f.E[2][i+1][j] - f.E[2][i+1][j+1])/(f.dm.Deltay); // at i + 1/2, j

    double Dx_Ez_AB = 0.5 * (Dx_Ez_A + Dx_Ez_B); // at i , j

    double Dy_Ez_A = f.D[1][i][j+1] * (f.E[2][i+1][j+1] - f.E[2][i][j+1])/(f.dm.Deltax[i]); // at i, j + 1/2
    double Dy_Ez_B = f.D[1][i][j] * (f.E[2][i+1][j] - f.E[2][i][j])/(f.dm.Deltax[i]); // at i, j - 1/2
    
    double Dy_Ez_AB = 0.5 * (Dy_Ez_A + Dy_Ez_B); // FIX: was Dy_Hz_A + Dy_Hz_B
    
    double Dz_Ey_Ex_AB = f.D[2][i][j] * ((f.E[0][i][j+1] - f.E[0][i][j])/(f.dm.Deltay) + (f.E[1][i][j] - f.E[1][i+1][j])/(f.dm.Deltax[i])); 

    double D_dot_curl_E_AB = Dx_Ez_AB + Dy_Ez_AB + Dz_Ey_Ex_AB; // at i,j

    // Calculates -D dot (curl(E)) at i-1 , j

    double Dx_Ez_C = f.D[0][i][j] * (f.E[2][i][j] - f.E[2][i][j+1])/(f.dm.Deltay); // at i - 1/2, j 
    double Dx_Ez_D = f.D[0][i-1][j] * (f.E[2][i-1][j] - f.E[2][i-1][j+1])/(f.dm.Deltay); // at i - 3/2, j fixed on 06/24/26 was E[2][i-1][j-1]

    double Dx_Ez_CD = 0.5 * (Dx_Ez_C + Dx_Ez_D); // FIX: was Dx_Ez_G + Dx_Ez_H

    double Dy_Ez_C = f.D[1][i-1][j+1] * (f.E[2][i][j+1] - f.E[2][i-1][j+1])/(f.dm.Deltax[i]); // at i-1, j + 1/2
    double Dy_Ez_D = f.D[1][i-1][j] * (f.E[2][i][j] - f.E[2][i-1][j])/(f.dm.Deltax[i]); // at i-1, j - 1/2

    double Dy_Ez_CD = 0.5 * (Dy_Ez_C + Dy_Ez_D); // FIX: was Dy_Ez_G + Dy_Ez_H

    double Dz_Ey_Ex_CD = f.D[2][i-1][j] * ((f.E[1][i][j] - f.E[1][i-1][j])/(f.dm.Deltax[i]) + (f.E[0][i-1][j+1] - f.E[0][i-1][j])/(f.dm.Deltay)); 

    double D_dot_curl_E_CD = Dx_Ez_CD + Dy_Ez_CD + Dz_Ey_Ex_CD; // at i-1,j

    double By_j = 0.5 * (f.B[1][i][j+1] + f.B[1][i][j]);
    double By_i_minus_1_j = 0.5 * (f.B[1][i-1][j+1] + f.B[1][i-1][j]); 
    
    double Bx_i_squared = 0.5 * (f.B[0][i+1][j] * f.B[0][i+1][j]) + 0.5 * (f.B[0][i][j] * f.B[0][i][j]);
    double Bx_i_minus_1_j_squared = 0.5 * (f.B[0][i][j] * f.B[0][i][j]) + 0.5 *(f.B[0][i-1][j] * f.B[0][i-1][j]);
    double By_j_squared = 0.5 * (f.B[1][i][j+1] * f.B[1][i][j+1]) + 0.5 * (f.B[1][i][j] * f.B[1][i][j]);
    double By_i_minus_1_j_squared = 0.5 * (f.B[1][i-1][j+1] * f.B[1][i-1][j+1]) + 0.5 * (f.B[1][i-1][j] * f.B[1][i-1][j]); 
    double Bz = f.B[2][i][j];
    double Bz_i_minus_1_j = f.B[2][i-1][j];

    double B_ij = (Bx_i_squared)  + (By_j_squared)  + (Bz * Bz); 
    double B_i_minus_1_j = (Bx_i_minus_1_j_squared) +(By_i_minus_1_j_squared) + (Bz_i_minus_1_j * Bz_i_minus_1_j);
    
    return 0.5 * ((D_dot_curl_E_AB * By_j)/(B_ij) + (D_dot_curl_E_CD * By_i_minus_1_j)/(B_i_minus_1_j));
}

double compute_A3_z(int i, int j, const Fields & f){
    // Calculation of A_3_z

    // Calculates -D dot (curl(E)) at i,j

    double Dx_Ez_A = f.D[0][i][j] * (f.E[2][i][j] - f.E[2][i][j+1])/(f.dm.Deltay); // at i - 1/2, j
    double Dx_Ez_B = f.D[0][i+1][j] * (f.E[2][i+1][j] - f.E[2][i+1][j+1])/(f.dm.Deltay); // at i + 1/2, j

    double Dx_Ez_AB = 0.5 * (Dx_Ez_A + Dx_Ez_B); // at i , j

    double Dy_Ez_A = f.D[1][i][j+1] * (f.E[2][i+1][j+1] - f.E[2][i][j+1])/(f.dm.Deltax[i]); // at i, j + 1/2
    double Dy_Ez_B = f.D[1][i][j] * (f.E[2][i+1][j] - f.E[2][i][j])/(f.dm.Deltax[i]); // at i, j - 1/2

    double Dy_Ez_AB = 0.5 * (Dy_Ez_A + Dy_Ez_B); // FIX: was Dy_Hz_A + Dy_Hz_B

    double Dz_Ey_Ex_AB = f.D[2][i][j] * ((f.E[0][i][j+1] - f.E[0][i][j])/(f.dm.Deltay) + (f.E[1][i][j] - f.E[1][i+1][j])/(f.dm.Deltax[i])); 

    double D_dot_curl_E_AB = Dx_Ez_AB + Dy_Ez_AB + Dz_Ey_Ex_AB; // at i,j

    // Calculates -D dot (curl(E)) at i , j-1

    double Dx_Ez_C = f.D[0][i+1][j-1] * (f.E[2][i+1][j-1] - f.E[2][i+1][j])/(f.dm.Deltay); // at i + 1/2, j - 1
    double Dx_Ez_D = f.D[0][i][j-1] * (f.E[2][i][j-1] - f.E[2][i][j])/(f.dm.Deltay); // at i - 1/2, j-1

    double Dx_Ez_CD = 0.5 * (Dx_Ez_C + Dx_Ez_D); // at i , j - 1

    double Dy_Ez_C = f.D[1][i][j] * (f.E[2][i+1][j] - f.E[2][i][j])/(f.dm.Deltax[i]); // at i, j - 1/2
    double Dy_Ez_D = f.D[1][i][j-1] * (f.E[2][i+1][j-1] - f.E[2][i][j-1])/(f.dm.Deltax[i]); // at i , j - 3/2

    double Dy_Ez_CD = 0.5 * (Dy_Ez_C + Dy_Ez_D); // at i, j - 1

    double Dz_Ey_Ex_CD = f.D[2][i][j-1] * ((f.E[0][i][j] - f.E[0][i][j-1])/(f.dm.Deltay) + (f.E[1][i][j-1] - f.E[1][i+1][j-1])/(f.dm.Deltax[i])); 

    double D_dot_curl_E_CD = Dx_Ez_CD + Dy_Ez_CD + Dz_Ey_Ex_CD; // at i,j-1

    // Calculates -D dot (curl(E)) at i-1, j-1

    double Dx_Ez_E = f.D[0][i][j-1] * (f.E[2][i][j-1] - f.E[2][i][j])/(f.dm.Deltay); // at i - 1/2, j-1 
    double Dx_Ez_F = f.D[0][i-1][j-1] * (f.E[2][i-1][j-1] - f.E[2][i-1][j])/(f.dm.Deltay); // at i - 3/2, j-1

    double Dx_Ez_EF = 0.5 * (Dx_Ez_E + Dx_Ez_F); // FIX: was Dx_Hz_E + Dx_Hz_F

    double Dy_Ez_E = f.D[1][i-1][j] * (f.E[2][i][j] - f.E[2][i-1][j])/(f.dm.Deltax[i]); // at i-1, j - 1/2
    double Dy_Ez_F = f.D[1][i-1][j-1] * (f.E[2][i][j-1] - f.E[2][i-1][j-1])/(f.dm.Deltax[i]); // at i-1, j - 3/2

    double Dy_Ez_EF = 0.5 * (Dy_Ez_E + Dy_Ez_F); // FIX: was Dy_Hz_E + Dy_Hz_F

    double Dz_Ey_Ex_EF = f.D[2][i-1][j-1] * ((f.E[1][i-1][j-1] - f.E[1][i][j-1])/(f.dm.Deltax[i]) + (f.E[0][i-1][j] - f.E[0][i-1][j-1])/(f.dm.Deltay)); 

    double D_dot_curl_E_EF = Dx_Ez_EF + Dy_Ez_EF + Dz_Ey_Ex_EF; // at i-1,j-1

    // Calculates -D dot (curl(E)) at i-1, j

    double Dx_Ez_G = f.D[0][i][j] * (f.E[2][i][j] - f.E[2][i][j+1])/(f.dm.Deltay); // at i - 1/2, j 
    double Dx_Ez_H = f.D[0][i-1][j] * (f.E[2][i-1][j] - f.E[2][i-1][j+1])/(f.dm.Deltay); // at i - 3/2, j

    double Dx_Ez_GH = 0.5 * (Dx_Ez_G + Dx_Ez_H); // at i-1 , j

    double Dy_Ez_G = f.D[1][i-1][j+1] * (f.E[2][i][j+1] - f.E[2][i-1][j+1])/(f.dm.Deltax[i]); // at i-1, j + 1/2
    double Dy_Ez_H = f.D[1][i-1][j] * (f.E[2][i][j] - f.E[2][i-1][j])/(f.dm.Deltax[i]); // at i-1, j - 1/2

    double Dy_Ez_GH = 0.5 * (Dy_Ez_G + Dy_Ez_H); // at i-1, j

    double Dz_Ey_Ex_GH = f.D[2][i-1][j] * ((f.E[1][i-1][j] - f.E[1][i][j])/(f.dm.Deltax[i]) + (f.E[0][i-1][j+1] - f.E[0][i-1][j])/(f.dm.Deltay)); 

    double D_dot_curl_E_GH = Dx_Ez_GH + Dy_Ez_GH + Dz_Ey_Ex_GH; // at i-1,j

    double Bx_i = 0.5 * (f.B[0][i+1][j] + f.B[0][i][j]);
    double Bx_i_minus_1_j = 0.5 * (f.B[0][i][j] + f.B[0][i-1][j]);
    double Bx_i_j_minus_1 = 0.5 * (f.B[0][i+1][j-1] + f.B[0][i][j-1]);
    double Bx_ij_minus_1 = 0.5 * (f.B[0][i][j-1] + f.B[0][i-1][j-1]);

    double By_j = 0.5 * (f.B[1][i][j+1] + f.B[1][i][j]);
    double By_j_minus_1 = 0.5 * (f.B[1][i][j-1] + f.B[1][i][j]); 
    double By_i_minus_1_j = 0.5 * (f.B[1][i-1][j+1] + f.B[1][i-1][j]);
    double By_ij_minus_1 = 0.5 * (f.B[1][i-1][j] + f.B[1][i-1][j-1]);  

    double Bz = f.B[2][i][j];
    double Bz_i_minus_1_j = f.B[2][i-1][j];
    double Bz_ij_minus_1 = f.B[2][i-1][j-1];
    double Bz_i_j_minus_1 = f.B[2][i][j-1];

    double B_ij = 0.5*(f.B[0][i][j]*f.B[0][i][j] + f.B[0][i+1][j]*f.B[0][i+1][j]) + 0.5*(f.B[1][i][j]*f.B[1][i][j] + f.B[1][i][j+1]*f.B[1][i][j+1]) + f.B[2][i][j]*f.B[2][i][j];
    double B_i_minus_1_j = 0.5*(f.B[0][i-1][j]*f.B[0][i-1][j] + f.B[0][i][j]*f.B[0][i][j]) + 0.5*(f.B[1][i-1][j]*f.B[1][i-1][j] + f.B[1][i-1][j+1]*f.B[1][i-1][j+1]) + f.B[2][i-1][j]*f.B[2][i-1][j];
    double B_i_j_minus_1 = 0.5*(f.B[0][i][j-1]*f.B[0][i][j-1] + f.B[0][i+1][j-1]*f.B[0][i+1][j-1]) + 0.5*(f.B[1][i][j-1]*f.B[1][i][j-1] + f.B[1][i][j]*f.B[1][i][j]) + f.B[2][i][j-1]*f.B[2][i][j-1];
    double B_ij_minus_1 = 0.5*(f.B[0][i-1][j-1]*f.B[0][i-1][j-1] + f.B[0][i][j-1]*f.B[0][i][j-1]) + 0.5*(f.B[1][i-1][j-1]*f.B[1][i-1][j-1] + f.B[1][i-1][j]*f.B[1][i-1][j]) + f.B[2][i-1][j-1]*f.B[2][i-1][j-1];

    return 0.25 * ((D_dot_curl_E_AB * Bz)/(B_ij) + (D_dot_curl_E_CD * Bz_i_j_minus_1)/(B_i_j_minus_1) + (D_dot_curl_E_EF * Bz_ij_minus_1)/(B_ij_minus_1) + (D_dot_curl_E_GH * Bz_i_minus_1_j)/(B_i_minus_1_j));
}


/*
    Computes current density J components in reduced units.

    Input: B, E, H, D: magnetic, electric, H, and displacement fields
           Rho: charge density (ScalarField)
           N_GC: number of ghost cells
           dm: Domain object containing information about the simulation domain
    Output: J: current density in reduced units
*/
void Compute_J(VectorField & B, VectorField & E, VectorField & H, VectorField & D, ScalarField & Rho, VectorField & J, size_t N_GC, const Domain & dm){
    Fields f{E, B, H, D, Rho, dm};

    for(size_t i=N_GC; i<B.shape()[1]-N_GC; i++){
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){
            
            J[0][i][j] = compute_A1_x(i, j, f) + compute_A2_x(i, j, f) + compute_A3_x(i, j, f);
            J[1][i][j] = compute_A1_y(i, j, f) + compute_A2_y(i, j, f) + compute_A3_y(i, j, f);
            J[2][i][j] = compute_A1_z(i, j, f) + compute_A2_z(i, j, f) + compute_A3_z(i, j, f);
            double By_ij = 0.5 * (B[1][i][j] + B[1][i][j+1]);
    double By_iminus1_j = 0.5 * (B[1][i-1][j] + B[1][i-1][j+1]);
    double By_iminus1_jminus1 = 0.5 * (B[1][i-1][j] + B[1][i-1][j-1]);
    double By_i_jminus1 = 0.5 * (B[1][i][j] + B[1][i][j-1]);
    double By_iplus1_jminus1 = 0.5 * (B[1][i+1][j] + B[1][i+1][j-1]);
    double By_iplus1_j = 0.5 * (B[1][i+1][j] + B[1][i+1][j+1]);
            

              if (i == 8 && j == 4){
                
                std::cout << slope_calc(E, 2, i, j, 1, dm.Deltax[i], dm.Deltay) << std::endl;
                std::cout << minmod((E[2][i][j+1] - E[2][i][j-1]) / (2 * dm.Deltay),  2 * (E[2][i][j] - E[2][i][j-1]) / dm.Deltay, 2 * (E[2][i][j+1] - E[2][i][j]) / dm.Deltay) << std::endl; 
        
    }
            
            }}
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



/*
    Computes the slope of an arbitrary VectorField object either in the x or y direction using the MCLimiter.

    Inputs:
        A: Arbitrary VectorField object
        component_index: 0 (x), 1 (y), or 2 (z)
        i_offset: offset for i index (e.g. i+1, i-1)
        j_offset: offset for j index (e.g. j+1, j-1)
        slope_direction: 0 for x slope, 1 for y slope
        Deltax
        Deltay

    Output: Slope value in the x or y direction.
*/

double slope_calc (VectorField & A, int component_index, int i_offset, int j_offset, int slope_direction, double Deltax, double Deltay){
    
    Compute_A_to_cell_center(A, component_index, i_offset, j_offset);

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
    Computes Ez fluxes using slope limiters to update B_x.

    Inputs: i, j, E, B, Deltax, and Deltay

    Outputs: Ez{i-1/2,j-1/2} and Ez{i-1/2, j+1/2} fluxes computed using Toth 2000 Eq(19).
*/

std::pair<double, double> Ez_Flux_Calculation_Qx(int i, int j, VectorField & E, VectorField & B, double Deltax, double Deltay){
    
    //Fluxes for Ez{i,j} aka Ez_{i-1/2, j-1/2}
    //G_1{i,j-1/2} calculation (G_1)

    double Bx_ij_up = Compute_A_to_cell_center(B, 0, i, j) - 0.5 * Deltay * slope_calc(B, 0, i, j, 1, Deltax, Deltay);
    double Bx_i_jminus1_down = Compute_A_to_cell_center(B, 0, i, j-1) + 0.5 * Deltay * slope_calc(B, 0, i, j-1, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper
    double G_1_up = E[2][i][j] - 0.5 * slope_calc(E, 2, i, j, 1, Deltax, Deltay) * Deltay;
    double G_1_down = E[2][i][j-1] + 0.5 * slope_calc(E, 2, i, j-1, 1, Deltax, Deltay) * Deltay;
    double G_1 = 0.5 * (G_1_up + G_1_down) - 0.5 * (Bx_ij_up - Bx_i_jminus1_down);

    //G_1{i-1,j-1/2} calculation (G_2)
    
    double Bx_iminus1_j_up = Compute_A_to_cell_center(B, 0, i-1, j) - 0.5 * Deltay * slope_calc(B, 0, i-1, j, 1, Deltax, Deltay);
    double Bx_ijminus1_down = Compute_A_to_cell_center(B, 0, i-1, j-1) + 0.5 * Deltay * slope_calc(B, 0, i-1, j-1, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper 
    double G_2_up = E[2][i-1][j] - 0.5 * slope_calc(E, 2, i-1, j, 1, Deltax, Deltay)  * Deltay;
    double G_2_down = E[2][i-1][j-1] + 0.5 * slope_calc(E, 2, i-1, j-1, 1, Deltax, Deltay) * Deltay;
    double G_2 = 0.5 * (G_2_up + G_2_down) - 0.5 * (Bx_iminus1_j_up - Bx_ijminus1_down);

    //F_2{i-1/2,j} calculation (F_1)

    double By_ij_up = Compute_A_to_cell_center(B, 1, i, j) - 0.5 * Deltax * slope_calc(B, 1, i, j, 0, Deltax, Deltay);
    double By_iminus1_j_down = Compute_A_to_cell_center(B, 1, i-1, j) + 0.5 * Deltax * slope_calc(B, 1, i-1, j, 0, Deltax, Deltay);

    double F_1_up = E[2][i][j] - 0.5 * slope_calc(E, 2, i, j, 0, Deltax, Deltay) * Deltax; 
    double F_1_down = E[2][i-1][j] + 0.5 * slope_calc(E, 2, i-1, j, 0, Deltax, Deltay)  * Deltax;

    //Note: Ez is negative in the F matrix in Yu paper

    double F_1 = 0.5 * (-F_1_up - F_1_down) - 0.5 * (By_ij_up - By_iminus1_j_down);

    //F_2{i-1/2,j-1} calculation (F_2)

    double By_i_jminus1_up = Compute_A_to_cell_center(B, 1, i, j-1) - 0.5 * Deltax * slope_calc(B, 1, i, j-1, 0, Deltax, Deltay);
    double By_ij_minus1_down = Compute_A_to_cell_center(B, 1, i-1, j-1) + 0.5 * Deltax * slope_calc(B, 1, i-1, j-1, 0, Deltax, Deltay);

    double F_2_up = E[2][i][j-1] - 0.5 * slope_calc(E, 2, i, j-1, 0, Deltax, Deltay) * Deltax; 
    double F_2_down = E[2][i-1][j-1] + 0.5 * slope_calc(E, 2, i-1, j-1, 0, Deltax, Deltay) * Deltax;

    //Note: Ez is negative in the F matrix in Yu paper

    double F_2 = 0.5 * (-F_2_up - F_2_down) - 0.5 * (By_i_jminus1_up - By_ij_minus1_down);

    //Fluxes for Ez{i,j+1} aka Ez_{i-1, j+1/2}
    //G_1{i,j+1/2} calculation (G_3)

    double Bx_i_jplus1_up = Compute_A_to_cell_center(B, 0, i, j+1) - 0.5 * Deltay * slope_calc(B, 0, i, j+1, 1, Deltax, Deltay);
    double Bx_ij_down = Compute_A_to_cell_center(B, 0, i, j) + 0.5 * Deltay * slope_calc(B, 0, i, j, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper
    double G_3_up = E[2][i][j+1] - 0.5 * slope_calc(E, 2, i, j+1, 1, Deltax, Deltay) * Deltay; 
    double G_3_down = E[2][i][j] + 0.5 * slope_calc(E, 2, i, j, 1, Deltax, Deltay) * Deltay;
    double G_3 = 0.5 * (G_3_up + G_3_down) - 0.5 * (Bx_i_jplus1_up - Bx_ij_down);
    
    //G_1{i-1,j+1/2} calculation (G_4)

    double Bx_iminus1_jplus1_up = Compute_A_to_cell_center(B, 0, i-1, j+1) - 0.5 * Deltay * slope_calc(B, 0, i-1, j+1, 1, Deltax, Deltay);
    double Bx_iminus1_j_down = Compute_A_to_cell_center(B, 0, i-1, j) + 0.5 * Deltay * slope_calc(B, 0, i-1, j, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper

    double G_4_up = E[2][i-1][j+1] - 0.5 * Deltay * slope_calc(E, 2, i-1, j+1, 1, Deltax, Deltay); 
    double G_4_down = E[2][i-1][j] + 0.5 * Deltay * slope_calc(E, 2, i-1, j, 1, Deltax, Deltay);

    double G_4 = 0.5 * (G_4_up + G_4_down) - 0.5 * (Bx_iminus1_jplus1_up - Bx_iminus1_j_down);

    //F_2{i-1/2,j} calculation (F_3) same thing as F_1

    //Note: Ez is negative in the F matrix in Yu paper

    double F_3 = 0.5 * (-F_1_up - F_1_down) - 0.5 * (By_ij_up - By_iminus1_j_down);

    //F_2{i-1/2,j+1} calculation (F_4)

    double By_i_jplus1_up = Compute_A_to_cell_center(B, 1, i, j+1) - 0.5 * Deltax * slope_calc(B, 1, i, j+1, 0, Deltax, Deltay);
    double By_iminus1_jplus1_down = Compute_A_to_cell_center(B, 1, i-1, j+1) + 0.5 * Deltax * slope_calc(B, 1, i-1, j+1, 0, Deltax, Deltay);

    double F_4_up = E[2][i][j+1] - 0.5 * Deltax * slope_calc(E, 2, i, j+1, 0, Deltax, Deltay); 
    double F_4_down = E[2][i-1][j+1] + 0.5 * Deltax * slope_calc(E, 2, i-1, j+1, 0, Deltax, Deltay);

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

std::pair<double, double> Ez_Flux_Calculation_Qy(int i, int j, VectorField & E, VectorField & B, double Deltax, double Deltay){

    //Fluxes for Ez{i,j} aka Ez_{i-1/2, j-1/2}
    //G_1{i,j-1/2} calculation (G_1)

    double Bx_ij_up = Compute_A_to_cell_center(B, 0, i, j) - 0.5 * Deltay * slope_calc(B, 0, i, j, 1, Deltax, Deltay);
    double Bx_i_jminus1_down = Compute_A_to_cell_center(B, 0, i, j-1) + 0.5 * Deltay * slope_calc(B, 0, i, j-1, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper
    double G_1_up = E[2][i][j] - 0.5 * slope_calc(E, 2, i, j, 1, Deltax, Deltay) * Deltay; 
    double G_1_down = E[2][i][j-1] + 0.5 * slope_calc(E, 2, i, j-1, 1, Deltax, Deltay) * Deltay;
    
    double G_1 = 0.5 * (G_1_up + G_1_down) - 0.5 * (Bx_ij_up - Bx_i_jminus1_down);

    //G_1{i-1,j-1/2} calculation (G_2)

    double Bx_iminus1_j_up = Compute_A_to_cell_center(B, 0, i-1, j) - 0.5 * Deltay * slope_calc(B, 0, i-1, j, 1, Deltax, Deltay);
    double Bx_ijminus1_down = Compute_A_to_cell_center(B, 0, i-1, j-1) + 0.5 * Deltay * slope_calc(B, 0, i-1, j-1, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper 
    double G_2_up = E[2][i-1][j] - 0.5 * slope_calc(E, 2, i-1, j, 1, Deltax, Deltax)  * Deltay;
    double G_2_down = E[2][i-1][j-1] + 0.5 * slope_calc(E, 2, i-1, j-1, 1, Deltax, Deltay) * Deltay;
    double G_2 = 0.5 * (G_2_up + G_2_down) - 0.5 * (Bx_iminus1_j_up - Bx_ijminus1_down);


    //F_2{i-1/2,j} calculation (F_1)

    double By_ij_up = Compute_A_to_cell_center(B, 1, i, j) - 0.5 * Deltax * slope_calc(B, 1, i, j, 0, Deltax, Deltay);
    double By_iminus1_j_down = Compute_A_to_cell_center(B, 1, i-1, j) + 0.5 * Deltax * slope_calc(B, 1, i-1, j, 0, Deltax, Deltay);

    double F_1_up = E[2][i][j] - 0.5 * slope_calc(E, 2, i, j, 0, Deltax, Deltay) * Deltax; 
    double F_1_down = E[2][i-1][j] + 0.5 * slope_calc(E, 2, i-1, j, 0, Deltax, Deltay)  * Deltax;

    //Note: Ez is negative in the F matrix in Yu paper

    double F_1 = 0.5 * (-F_1_up - F_1_down) - 0.5 * (By_ij_up - By_iminus1_j_down);

    //F_2{i-1/2,j-1} calculation (F_2)

    double By_i_jminus1_up = Compute_A_to_cell_center(B, 1, i, j-1) - 0.5 * Deltax * slope_calc(B, 1, i, j-1, 0, Deltax, Deltay);
    double By_ij_minus1_down = Compute_A_to_cell_center(B, 1, i-1, j-1) + 0.5 * Deltax * slope_calc(B, 1, i-1, j-1, 0, Deltax, Deltay);

    double F_2_up = E[2][i][j-1] - 0.5 * slope_calc(E, 2, i, j-1, 0, Deltax, Deltay) * Deltax; 
    double F_2_down = E[2][i-1][j-1] + 0.5 * slope_calc(E, 2, i-1, j-1, 0, Deltax, Deltay) * Deltax;

    double F_2 = 0.5 * (-F_2_up - F_2_down) - 0.5 * (By_i_jminus1_up - By_ij_minus1_down);

    //E_z Fluxes for Ez{i+1,j} aka Ez_{i+1/2, j-1/2}

    //G_1{i,j-1/2} calculation (G_3), same as G_1

    //Note: Ez is positive in the G matrix in Yu paper

    double G_3 = 0.5 * (G_1_up + G_1_down) - 0.5 * (Bx_ij_up - Bx_i_jminus1_down);

    //G_1{i+1, j-1/2}  calculation (G_4)

    double Bx_iplus1_j_up = Compute_A_to_cell_center(B, 0, i+1, j) - 0.5 * Deltay * slope_calc(B, 0, i+1, j, 1, Deltax, Deltay);
    double Bx_iplus1_jminus1_down = Compute_A_to_cell_center(B, 0, i+1, j-1) + 0.5 * Deltay * slope_calc(B, 0, i+1, j-1, 1, Deltax, Deltay);

    double G_4_up = E[2][i+1][j] - 0.5 * Deltay * slope_calc(E, 2, i+1, j, 1, Deltax, Deltay); 
    double G_4_down = E[2][i+1][j-1] + 0.5 * Deltay * slope_calc(E, 2, i+1, j-1, 1, Deltax, Deltay);

    double G_4 = 0.5 * (G_4_up + G_4_down);

    //F_1{i+1/2, j-1} calculation (F_3)

    double By_iplus1_jminus1_up = Compute_A_to_cell_center(B, 1, i+1, j-1) - 0.5 * Deltax * slope_calc(B, 1, i+1, j-1, 0, Deltax, Deltay);
    double By_i_jminus1_down = Compute_A_to_cell_center(B, 1, i, j-1) + 0.5 * Deltax *  slope_calc(B, 1, i, j-1, 0, Deltax, Deltay);
    
    double F_3_up = E[2][i+1][j-1] - 0.5 * slope_calc(E, 2, i+1, j-1, 0, Deltax, Deltay) * Deltax; 
    double F_3_down = E[2][i][j-1] + 0.5 * slope_calc(E, 2, i, j-1, 0, Deltax, Deltay)  * Deltax;

    double F_3 = 0.5 * (F_3_up + F_3_down) - 0.5 * (By_iplus1_jminus1_up - By_i_jminus1_down);

    //F_1{i+1/2,j} calculation (F_4)

    double By_iplus1_j_up = Compute_A_to_cell_center(B, 1, i+1, j) - 0.5 * Deltax * slope_calc(B, 1, i+1, j, 0, Deltax, Deltay);
    double By_ij_down = Compute_A_to_cell_center(B, 1, i, j) + 0.5 * Deltax * slope_calc(B, 1, i, j, 0, Deltax, Deltay); 

    double F_4_up = E[2][i+1][j] - 0.5 * slope_calc(E, 2, i+1, j, 0, Deltax, Deltay) * Deltax;
    double F_4_down = E[2][i][j] + 0.5 * slope_calc(E, 2, i, j, 0, Deltax, Deltay)  * Deltax; 

    double F_4 = 0.5 * (F_4_up + F_4_down) - 0.5 * (By_iplus1_j_up - By_ij_down);

    double Ez_ij = 0.25 * (-F_1 - F_2 + G_1 + G_2);
    double Ez_iplus1_j = 0.25 * (-F_3 - F_4 + G_3 + G_4);

    return {Ez_ij, Ez_iplus1_j};

}

/*
    Computes reconstructed fluxes using slope limiters to update B_z.
    Inputs: i, j, D, B, Deltax, and Deltay

    Outputs: Update for Bz by computing F_3/F_4 (Ey), G_3/G_4 (Ex) fluxes computed using Yu 2011 Eq(9).
*/
double Qz_reconstruction(int i, int j, VectorField & E, VectorField & B, double Deltax, double Deltay){

    //G_3{i,j-1/2} calculation G_3

    double G_3_up = Compute_A_to_cell_center(E, 0, i, j) - 0.5 * Deltay * slope_calc(E, 0, i, j, 1, Deltax, Deltay);
    double G_3_down = Compute_A_to_cell_center(E, 0, i, j-1) + 0.5 * Deltay * slope_calc(E, 0, i, j-1, 1, Deltax, Deltay);

    double Bz_ij_up_y = B[2][i][j] - 0.5 * Deltay * slope_calc(B, 2, i, j, 1, Deltax, Deltay);
    double Bz_i_jminus1_down_y = B[2][i][j-1] + 0.5 * Deltay * slope_calc(B, 2, i, j-1, 1, Deltax, Deltay);

    double G_3 = 0.5 * (-G_3_up - G_3_down) - 0.5 * (Bz_ij_up_y - Bz_i_jminus1_down_y);

    //G_3{i,j+1/2} calculation G_4

    double G_4_up = Compute_A_to_cell_center(E, 0, i, j+1) - 0.5 * Deltay * slope_calc(E, 0, i, j+1, 1, Deltax, Deltay); 
    double G_4_down = Compute_A_to_cell_center(E, 0, i, j) + 0.5 * Deltay * slope_calc(E, 0, i, j, 1, Deltax, Deltay);

    double Bz_i_jplus1_up_y = B[2][i][j+1] - 0.5 * Deltay * slope_calc(B, 2, i, j+1, 1, Deltax, Deltay);
    double Bz_ij_down_y = B[2][i][j] + 0.5 * Deltay * slope_calc(B, 2, i, j, 1, Deltax, Deltay);

    double G_4 = 0.5 * (-G_4_up - G_4_down) - 0.5 * (Bz_i_jplus1_up_y - Bz_ij_down_y);

    //F_3{i-1/2,j} calculation (F_3)

    double F_3_up = Compute_A_to_cell_center(E, 1, i, j) - 0.5 * Deltax * slope_calc(E, 1, i, j, 0, Deltax, Deltay);
    double F_3_down = Compute_A_to_cell_center(E, 1, i-1, j) + 0.5 * Deltax * slope_calc(E, 1, i-1, j, 0, Deltax, Deltay);
    
    double Bz_ij_up_x = B[2][i][j] - 0.5 * Deltax * slope_calc(B, 2, i, j, 0, Deltax, Deltay);
    double Bz_iminus1_j_down_x = B[2][i-1][j] + 0.5 * Deltax * slope_calc(B, 2, i-1, j, 0, Deltax, Deltay);

    double F_3 = 0.5 * (F_3_up + F_3_down) - 0.5 * (Bz_ij_up_x - Bz_iminus1_j_down_x);
    
    //F_3{i+1/2,j} calculation (F_4)

    double F_4_up = Compute_A_to_cell_center(E, 1, i+1, j) - 0.5 * Deltax * slope_calc(E, 1, i+1, j, 0, Deltax, Deltay);
    double F_4_down = Compute_A_to_cell_center(E, 1, i, j) + 0.5 * Deltax * slope_calc(E, 1, i, j, 0, Deltax, Deltay);
    
    double Bz_iplus1_j_up_x = B[2][i+1][j] - 0.5 * Deltax * slope_calc(B, 2, i+1, j, 0, Deltax, Deltay);
    double Bz_ij_down_x = B[2][i][j] + 0.5 * Deltax * slope_calc(B, 2, i, j, 0, Deltax, Deltay);

    double F_4 = 0.5 * (F_4_up + F_4_down) - 0.5 * (Bz_iplus1_j_up_x - Bz_ij_down_x);

    //Eq(9) in Yu paper
    double Qz = - (F_4 - F_3)/(Deltax) - (G_4 - G_3)/(Deltay);

    return Qz;

}

/*
    Computes Hz fluxes using slope limiters to update D_x. Calculation is analagous to the Ez flux calculation for the Bx
    time evolution, hence, copied and pasted the Ez_Flux_Calculation_Qx function.

    Inputs: i, j, B, D, Deltax, and Deltay

    Outputs: Hz{i-1/2,j-1/2} and Hz{i-1/2, j+1/2} fluxes computed using Toth 2000 Eq(19).
*/
std::pair<double, double> Hz_Flux_Calculation_Fx(int i, int j, VectorField & H, VectorField & D, double Deltax, double Deltay){

    //Fluxes for Ez{i,j} aka Ez_{i-1/2, j-1/2}
    //G_1{i,j-1/2} calculation (G_1)

    double Dx_ij_up = Compute_A_to_cell_center(D, 0, i, j) - 0.5 * Deltay * slope_calc(D, 0, i, j, 1, Deltax, Deltay);
    double Dx_i_jminus1_down = Compute_A_to_cell_center(D, 0, i, j-1) + 0.5 * Deltay * slope_calc(D, 0, i, j-1, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper
    double G_1_up = H[2][i][j] - 0.5 * slope_calc(H, 2, i, j, 1, Deltax, Deltay) * Deltay;
    double G_1_down = H[2][i][j-1] + 0.5 * slope_calc(H, 2, i, j-1, 1, Deltax, Deltay) * Deltay;
    double G_1 = 0.5 * (G_1_up + G_1_down) - 0.5 * (Dx_ij_up - Dx_i_jminus1_down);

    //G_1{i-1,j-1/2} calculation (G_2)
    
    double Dx_iminus1_j_up = Compute_A_to_cell_center(D, 0, i-1, j) - 0.5 * Deltay * slope_calc(D, 0, i-1, j, 1, Deltax, Deltay);
    double Dx_ijminus1_down = Compute_A_to_cell_center(D, 0, i-1, j-1) + 0.5 * Deltay * slope_calc(D, 0, i-1, j-1, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper 
    double G_2_up = H[2][i-1][j] - 0.5 * slope_calc(H, 2, i-1, j, 1, Deltax, Deltay)  * Deltay;
    double G_2_down = H[2][i-1][j-1] + 0.5 * slope_calc(H, 2, i-1, j-1, 1, Deltax, Deltay) * Deltay;
    double G_2 = 0.5 * (G_2_up + G_2_down) - 0.5 * (Dx_iminus1_j_up - Dx_ijminus1_down);


    //F_2{i-1/2,j} calculation (F_1)

    double Dy_ij_up = Compute_A_to_cell_center(D, 1, i, j) - 0.5 * Deltax * slope_calc(D, 1, i, j, 0, Deltax, Deltay);
    double Dy_iminus1_j_down = Compute_A_to_cell_center(D, 1, i-1, j) + 0.5 * Deltax * slope_calc(D, 1, i-1, j, 0, Deltax, Deltay);

    double F_1_up = H[2][i][j] - 0.5 * slope_calc(H, 2, i, j, 0, Deltax, Deltay) * Deltax; 
    double F_1_down = H[2][i-1][j] + 0.5 * slope_calc(H, 2, i-1, j, 0, Deltax, Deltay)  * Deltax;

    //Note: Ez is negative in the F matrix in Yu paper

    double F_1 = 0.5 * (-F_1_up - F_1_down) - 0.5 * (Dy_ij_up - Dy_iminus1_j_down);

    //F_2{i-1/2,j-1} calculation (F_2)

    double Dy_i_jminus1_up = Compute_A_to_cell_center(D, 1, i, j-1) - 0.5 * Deltax * slope_calc(D, 1, i, j-1, 0, Deltax, Deltay);
    double Dy_ij_minus1_down = Compute_A_to_cell_center(D, 1, i-1, j-1) + 0.5 * Deltax * slope_calc(D, 1, i-1, j-1, 0, Deltax, Deltay);

    double F_2_up = H[2][i][j-1] - 0.5 * slope_calc(H, 2, i, j-1, 0, Deltax, Deltay) * Deltax; 
    double F_2_down = H[2][i-1][j-1] + 0.5 * slope_calc(H, 2, i-1, j-1, 0, Deltax, Deltay) * Deltax;

    double F_2 = 0.5 * (-F_2_up - F_2_down) - 0.5 * (Dy_i_jminus1_up - Dy_ij_minus1_down);

    //Fluxes for Ez{i,j+1} aka Ez_{i-1, j+1/2}
    //G_1{i,j+1/2} calculation (G_3)

    double Dx_i_jplus1_up = Compute_A_to_cell_center(D, 0, i, j+1) - 0.5 * Deltay * slope_calc(D, 0, i, j+1, 1, Deltax, Deltay);
    double Dx_ij_down = Compute_A_to_cell_center(D, 0, i, j) + 0.5 * Deltax * slope_calc(D, 0, i, j, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper
    double G_3_up = H[2][i][j+1] - 0.5 * slope_calc(H, 2, i, j+1, 1, Deltax, Deltay) * Deltay; 
    double G_3_down = H[2][i][j] + 0.5 * slope_calc(H, 2, i, j, 1, Deltax, Deltay) * Deltay;
    double G_3 = 0.5 * (G_3_up + G_3_down) - 0.5 * (Dx_i_jplus1_up - Dx_ij_down);
    
    
    //G_1{i-1,j+1/2} calculation (G_4)

    double Dx_iminus1_jplus1_up = Compute_A_to_cell_center(D, 0, i-1, j+1) - 0.5 * Deltay * slope_calc(D, 0, i-1, j+1, 1, Deltax, Deltay);
    double Dx_iminus1_j_down = Compute_A_to_cell_center(D, 0, i-1, j) + 0.5 * Deltay * slope_calc(D, 0, i-1, j, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper

    double G_4_up = H[2][i-1][j+1] - 0.5 * Deltay * slope_calc(H, 2, i-1, j+1, 1, Deltax, Deltay); 
    double G_4_down = H[2][i-1][j] + 0.5 * Deltay * slope_calc(H, 2, i-1, j, 1, Deltax, Deltay);

    double G_4 = 0.5 * (G_4_up + G_4_down) - 0.5 * (Dx_iminus1_jplus1_up - Dx_iminus1_j_down);

    //F_2{i-1/2,j} calculation (F_3) same thing as F_1

    //Note: Ez is negative in the F matrix in Yu paper and no dissipation term since Bx lives w/ F_1

    double F_3 = 0.5 * (-F_1_up - F_1_down) - 0.5 * (Dy_ij_up - Dy_iminus1_j_down);

    //F_2{i-1/2,j+1} calculation (F_4)

    double Dy_i_jplus1_up = Compute_A_to_cell_center(D, 1, i, j+1) - 0.5 * Deltax * slope_calc(D, 1, i, j+1, 0, Deltax, Deltay);
    double Dy_iminus1_jplus1_down = Compute_A_to_cell_center(D, 1, i-1, j+1) + 0.5 * Deltax * slope_calc(D, 1, i-1, j+1, 0, Deltax, Deltay);

    double F_4_up = H[2][i][j+1] - 0.5 * Deltax * slope_calc(H, 2, i, j+1, 0, Deltax, Deltay); 
    double F_4_down = H[2][i-1][j+1] + 0.5 * Deltax * slope_calc(H, 2, i-1, j+1, 0, Deltax, Deltay);

    double F_4 = 0.5 * (-F_4_up - F_4_down) - 0.5 * (Dy_i_jplus1_up - Dy_iminus1_jplus1_down);

    //Toth 2000 Eq(19)
    double Hz_ij = 0.25 * (-F_1 - F_2 + G_1 + G_2);
    double Hz_i_jplus1 = 0.25 * (-F_3 - F_4 + G_3 + G_4);

    return {Hz_ij, Hz_i_jplus1};
}

/*
    Computes Hz fluxes using slope limiters to update D_y. Calculation is analagous to the Ez flux calculation for the By
    time evolution, hence, copied and pasted the Ez_Flux_Calculation_Qy function.

    Inputs: i, j, B, D, Deltax, and Deltay

    Outputs: Hz{i-1/2,j-1/2} and Hz{i+1/2, j-1/2} fluxes computed using Toth 2000 Eq(19).
*/

std::pair<double, double> Hz_Flux_Calculation_Fy(int i, int j, VectorField & H, VectorField & D, double Deltax, double Deltay){

    //Fluxes for Ez{i,j} aka Ez_{i-1/2, j-1/2}
    //G_1{i,j-1/2} calculation (G_1)

    double Dx_ij_up = Compute_A_to_cell_center(D, 0, i, j) - 0.5 * Deltay * slope_calc(D, 0, i, j, 1, Deltax, Deltay);
    double Dx_i_jminus1_down = Compute_A_to_cell_center(D, 0, i, j-1) + 0.5 * Deltay * slope_calc(D, 0, i, j-1, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper
    double G_1_up = H[2][i][j] - 0.5 * slope_calc(H, 2, i, j, 1, Deltax, Deltay) * Deltay; 
    double G_1_down = H[2][i][j-1] + 0.5 * slope_calc(H, 2, i, j-1, 1, Deltax, Deltay) * Deltay;
    
    double G_1 = 0.5 * (G_1_up + G_1_down) - 0.5 * (Dx_ij_up - Dx_i_jminus1_down);

    //G_1{i-1,j-1/2} calculation (G_2)

    double Dx_iminus1_j_up = Compute_A_to_cell_center(D, 0, i-1, j) - 0.5 * Deltay * slope_calc(D, 0, i-1, j, 1, Deltax, Deltay);
    double Dx_ijminus1_down = Compute_A_to_cell_center(D, 0, i-1, j-1) + 0.5 * Deltay * slope_calc(D, 0, i-1, j-1, 1, Deltax, Deltay);

    //Note: Ez is positive in the G matrix in Yu paper 
    double G_2_up = H[2][i-1][j] - 0.5 * slope_calc(H, 2, i-1, j, 1, Deltax, Deltax)  * Deltay;
    double G_2_down = H[2][i-1][j-1] + 0.5 * slope_calc(H, 2, i-1, j-1, 1, Deltax, Deltay) * Deltay;
    double G_2 = 0.5 * (G_2_up + G_2_down) - 0.5 * (Dx_iminus1_j_up - Dx_ijminus1_down);

    //F_2{i-1/2,j} calculation (F_1)

    double Dy_ij_up = Compute_A_to_cell_center(D, 1, i, j) - 0.5 * Deltax * slope_calc(D, 1, i, j, 0, Deltax, Deltay);
    double Dy_iminus1_j_down = Compute_A_to_cell_center(D, 1, i-1, j) + 0.5 * Deltax * slope_calc(D, 1, i-1, j, 0, Deltax, Deltay);

    double F_1_up = H[2][i][j] - 0.5 * slope_calc(H, 2, i, j, 0, Deltax, Deltay) * Deltax; 
    double F_1_down = H[2][i-1][j] + 0.5 * slope_calc(H, 2, i-1, j, 0, Deltax, Deltay)  * Deltax;

    //Note: Ez is negative in the F matrix in Yu paper

    double F_1 = 0.5 * (-F_1_up - F_1_down) - 0.5 * (Dy_ij_up - Dy_iminus1_j_down);

    //F_2{i-1/2,j-1} calculation (F_2)

    double Dy_i_jminus1_up = Compute_A_to_cell_center(D, 1, i, j-1) - 0.5 * Deltax * slope_calc(D, 1, i, j-1, 0, Deltax, Deltay);
    double Dy_ij_minus1_down = Compute_A_to_cell_center(D, 1, i-1, j-1) + 0.5 * Deltax * slope_calc(D, 1, i-1, j-1, 0, Deltax, Deltay);

    double F_2_up = H[2][i][j-1] - 0.5 * slope_calc(H, 2, i, j-1, 0, Deltax, Deltay) * Deltax; 
    double F_2_down = H[2][i-1][j-1] + 0.5 * slope_calc(H, 2, i-1, j-1, 0, Deltax, Deltay) * Deltax;

    double F_2 = 0.5 * (-F_2_up - F_2_down) - 0.5 * (Dy_i_jminus1_up - Dy_ij_minus1_down);

    //E_z Fluxes for Ez{i+1,j} aka Ez_{i+1/2, j-1/2}

    //G_1{i,j-1/2} calculation (G_3), same as G_1

    //Note: Ez is positive in the G matrix in Yu paper

    double G_3 = 0.5 * (G_1_up + G_1_down) - 0.5 * (Dx_ij_up - Dx_i_jminus1_down);

    //G_1{i+1, j-1/2}  calculation (G_4)

    double Dx_iplus1_j_up = Compute_A_to_cell_center(D, 0, i+1, j) - 0.5 * Deltay * slope_calc(D, 0, i+1, j, 1, Deltax, Deltay);
    double Dx_iplus1_jminus1_down = Compute_A_to_cell_center(D, 0, i+1, j-1) + 0.5 * Deltay * slope_calc(D, 0, i+1, j-1, 1, Deltax, Deltay);

    double G_4_up = H[2][i+1][j] - 0.5 * Deltay * slope_calc(H, 2, i+1, j, 1, Deltax, Deltay); 
    double G_4_down = H[2][i+1][j-1] + 0.5 * Deltay * slope_calc(H, 2, i+1, j-1, 1, Deltax, Deltay);

    double G_4 = 0.5 * (G_4_up + G_4_down) - 0.5 * (Dx_iplus1_j_up - Dx_iplus1_jminus1_down);

    //F_1{i+1/2, j-1} calculation (F_3)

    double Dy_iplus1_jminus1_up = Compute_A_to_cell_center(D, 1, i+1, j-1) - 0.5 * Deltax * slope_calc(D, 1, i+1, j-1, 0, Deltax, Deltay);
    double Dy_i_jminus1_down = Compute_A_to_cell_center(D, 1, i, j-1) + 0.5 * Deltax *  slope_calc(D, 1, i, j-1, 0, Deltax, Deltay);
    
    double F_3_up = H[2][i+1][j-1] - 0.5 * slope_calc(H, 2, i+1, j-1, 0, Deltax, Deltay) * Deltax; 
    double F_3_down = H[2][i][j-1] + 0.5 * slope_calc(H, 2, i, j-1, 0, Deltax, Deltay)  * Deltax;

    double F_3 = 0.5 * (F_3_up + F_3_down) - 0.5 * (Dy_iplus1_jminus1_up - Dy_i_jminus1_down);

    //F_1{i+1/2,j} calculation (F_4)

    double Dy_iplus1_j_up = Compute_A_to_cell_center(D, 1, i+1, j) - 0.5 * Deltax * slope_calc(D, 1, i+1, j, 0, Deltax, Deltay);
    double Dy_ij_down = Compute_A_to_cell_center(D, 1, i, j) + 0.5 * Deltax * slope_calc(D, 1, i, j, 0, Deltax, Deltay); 

    double F_4_up = H[2][i+1][j] - 0.5 * slope_calc(H, 2, i+1, j, 0, Deltax, Deltay) * Deltax;
    double F_4_down = H[2][i][j] + 0.5 * slope_calc(H, 2, i, j, 0, Deltax, Deltay)  * Deltax; 

    double F_4 = 0.5 * (F_4_up + F_4_down) - 0.5 * (Dy_iplus1_j_up - Dy_ij_down);

    double Ez_ij = 0.25 * (-F_1 - F_2 + G_1 + G_2);
    double Ez_iplus1_j = 0.25 * (-F_3 - F_4 + G_3 + G_4);

    return {Ez_ij, Ez_iplus1_j};

}

/*
    Computes reconstructed fluxes using slope limiters to update D_z.
    Inputs: i, j, B, D, Deltax, and Deltay

    Outputs: Update for Dz by computing F_3/F_4 (Hy), G_3/G_4 (Hx) fluxes computed using Yu 2011 Eq(9).
*/

double Fz_reconstruction(int i, int j, VectorField & E, VectorField & B, double Deltax, double Deltay){

    //G_3{i,j-1/2} calculation G_3

    double G_3_up = Compute_A_to_cell_center(E, 0, i, j) - 0.5 * Deltay *  slope_calc(E, 0, i, j, 1, Deltax, Deltay);
    double G_3_down = Compute_A_to_cell_center(E, 0, i, j-1) + 0.5 * Deltay * slope_calc(E, 0, i, j-1, 1, Deltax, Deltay);

    double Bz_ij_up_y = B[2][i][j] - 0.5 * Deltay * slope_calc(B, 2, i, j, 1, Deltax, Deltay); 
    double Bz_i_jminus1_down_y = B[2][i][j-1] + 0.5 * Deltay * slope_calc(B, 2, i, j-1, 1, Deltax, Deltay);

    double G_3 = 0.5 * (-G_3_up - G_3_down) - 0.5 * (Bz_ij_up_y - Bz_i_jminus1_down_y);

    //G_3{i,j+1/2} calculation G_4

    double G_4_up = Compute_A_to_cell_center(E, 0, i, j+1) - 0.5 * Deltay * slope_calc(E, 0, i, j+1, 1, Deltax, Deltay); 
    double G_4_down = Compute_A_to_cell_center(E, 0, i, j) + 0.5 * Deltay * slope_calc(E, 0, i, j, 1, Deltax, Deltay);

    double Bz_i_jplus1_up_y = B[2][i][j+1] - 0.5 * Deltay * slope_calc(B, 2, i, j+1, 1, Deltax, Deltay);
    double Bz_ij_down_y = B[2][i][j] + 0.5 * Deltay * slope_calc(B, 2, i, j, 1, Deltax, Deltay);

    double G_4 = 0.5 * (-G_4_up - G_4_down) - 0.5 * (Bz_i_jplus1_up_y - Bz_ij_down_y);

    //F_3{i-1/2,j} calculation (F_3)

    double F_3_up = Compute_A_to_cell_center(E, 1, i, j) - 0.5 * Deltax * slope_calc(E, 1, i, j, 0, Deltax, Deltay);
    double F_3_down = Compute_A_to_cell_center(E, 1, i-1, j) + 0.5 * Deltax * slope_calc(E, 1, i-1, j, 0, Deltax, Deltay);
    
    double Bz_ij_up_x = B[2][i][j] - 0.5 * Deltax * slope_calc(B, 2, i, j, 0, Deltax, Deltay);
    double Bz_iminus1_j_down_x = B[2][i-1][j] + 0.5 * Deltax * slope_calc(B, 2, i-1, j, 0, Deltax, Deltay);

    double F_3 = 0.5 * (F_3_up + F_3_down) - 0.5 * (Bz_ij_up_x + Bz_iminus1_j_down_x);
    
    //F_3{i+1/2,j} calculation (F_4)

    double F_4_up = Compute_A_to_cell_center(E, 1, i+1, j) - 0.5 * Deltax * slope_calc(E, 1, i+1, j, 0, Deltax, Deltay);
    double F_4_down = Compute_A_to_cell_center(E, 1, i, j) + 0.5 * Deltax * slope_calc(E, 1, i, j, 0, Deltax, Deltay);
    
    double Bz_iplus1_j_up_x = B[2][i+1][j] - 0.5 * Deltax * slope_calc(B, 2, i+1, j, 0, Deltax, Deltay);
    double Bz_ij_down_x = B[2][i][j] + 0.5 * Deltax * slope_calc(B, 2, i, j, 0, Deltax, Deltay); 

    double F_4 = 0.5 * (F_4_up + F_4_down) - 0.5 * (Bz_iplus1_j_up_x - Bz_ij_down_x);

    //Eq(9) in Yu paper
    double Qz = - (F_4 - F_3)/(Deltax) - (G_4 - G_3)/(Deltay);

    return Qz;
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

    return lambda0;
}


/*
    Computes E and H from the cell-centered average values of D and B.
    Conditional checks whether QED corrections are included.

    Inputs:
        E, H, D, B, SimParams object, Domain object

    Outputs:
        Ex, Ey, Ez, Hx, Hy, Hz for time-evolving B and D
*/
void Compute_EH_from_DB(VectorField & E, VectorField & H, VectorField & D, VectorField & B, const SimParams & params,  const Domain & dm) {

    for (size_t i = dm.N_GC; i < B.shape()[1]-dm.N_GC; i++) {
        for (size_t j = dm.N_GC; j < B.shape()[2]-dm.N_GC; j++) {

            double Dx = Compute_A_to_cell_center(D, 0, i, j);
            double Dy = Compute_A_to_cell_center(D, 1, i, j);
            double Dz = D[2][i][j];

            double Bx = Compute_A_to_cell_center(B, 0, i, j);
            double By = Compute_A_to_cell_center(B, 1, i, j);
            double Bz = B[2][i][j];

            if (params.QED_corrections){

                double D_squared = Dx*Dx + Dy*Dy + Dz*Dz;
                double B_squared = Bx*Bx + By*By + Bz*Bz;


                //double lambda = Solve_lambda_Newton(D_squared, B_squared, kappa, 20, 1e-12);
                double kappa = (4*alpha_e)/(45*4.4*B_0);
                double lambda = Solve_Lambda_Cubic(D_squared, B_squared);

                E[0][i][j] = Dx/lambda;
                E[1][i][j] = Dy/lambda;
                E[2][i][j] = Dz/lambda;

                H[0][i][j] = lambda * Bx;
                H[1][i][j] = lambda * By;
                H[2][i][j] = lambda * Bz;

            }

            //Trivial case (no QED corrections)
            else{

                E[0][i][j] = Dx;
                E[1][i][j] = Dy;
                E[2][i][j] = Dz;

                H[0][i][j] = Bx;
                H[1][i][j] = By;
                H[2][i][j] = Bz;

            }

        }
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
    double K_abs = 40.0;

    double damping_term = (K_abs)/(dm.Deltat) * ((dm.x[i] - 0.9 * max_dm_x)/(max_dm_x - 0.9 * max_dm_x)) * ((dm.x[i] - 0.9 * max_dm_x)/(max_dm_x - 0.9 * max_dm_x)) * ((dm.x[i] - 0.9 * max_dm_x)/(max_dm_x - 0.9 * max_dm_x));
    
    return damping_term;
}


/*
    This computes the RHS of the time-evolution equations for B and D. 

    Inputs:
        Qx, Qy, Qz, Fx, Fy, Fz, Rho, D, B, SimParams object, Domain object, Process Object

    Outputs:
        Updated Qx, Qy, and Qz for B time-evolution equations
        Updated Fx, Fy, and Fz for D time-evolution equations
*/
void Compute_RHS(ScalarField & Qx, ScalarField & Qy, ScalarField & Qz, ScalarField & Fx, ScalarField & Fy, ScalarField & Fz, ScalarField & Rho, VectorField & E, VectorField & H, VectorField & D, VectorField & B, const SimParams & params, const Domain & dm, const Process & ps)
{
    Compute_EH_from_DB(E, H, D, B, params, dm);

    double max_dm_x = *std::max_element(dm.x.begin(), dm.x.end());

    for(size_t i=dm.N_GC; i<D.shape()[1]-dm.N_GC-1; i++){
        for(size_t j=dm.N_GC; j<D.shape()[2]-dm.N_GC-1; j++){
            
            //Change everything from E to D and H to B
            double Bx_avg_ij = Compute_A_to_cell_center(B, 0, i, j);
            double Bx_avg_iplus1_j = Compute_A_to_cell_center(B, 0, i+1, j);
            double Bx_avg_i_jplus1 = Compute_A_to_cell_center(B, 0, i, j+1);
            double Bx_avg_i_jplus2 = Compute_A_to_cell_center(B, 0, i, j+2);
            double Bx_avg_iplus2_j = Compute_A_to_cell_center(B, 0, i+2, j);
            double By_avg_ij = Compute_A_to_cell_center(B, 1, i, j);
            double By_avg_iplus1_j = Compute_A_to_cell_center(B, 1, i+1, j);
            double By_avg_i_jplus1 = Compute_A_to_cell_center(B, 1, i, j+1);
            double By_avg_i_jplus2 = Compute_A_to_cell_center(B, 1, i, j+2);


            double Dx_avg_ij = Compute_A_to_cell_center(D, 0, i, j);
            double Dx_avg_iplus1_j = Compute_A_to_cell_center(D, 0, i+1, j);
            double Dx_avg_iplus2_j = Compute_A_to_cell_center(D, 0, i+2, j);
            double Dx_avg_i_jplus1 = Compute_A_to_cell_center(D, 0, i, j+1);
            double Dx_avg_i_jplus2 = Compute_A_to_cell_center(D, 0, i, j+2);
            double Dy_avg_ij = Compute_A_to_cell_center(D, 1, i, j);
            double Dy_avg_i_jplus1 = Compute_A_to_cell_center(D, 1, i, j+1);
            double Dy_avg_i_jplus2 = Compute_A_to_cell_center(D, 1, i, j+2);

            //Used for calculating (curl(H) dot B - curl(E) dot D)_ij
            auto [Ez_Qx_left_bottom, Ez_Qx_left_top] = Ez_Flux_Calculation_Qx(i, j, D, B, dm.Deltax[i], dm.Deltay); //For i-1/2, j face
            auto [Ez_Qx_right_bottom, Ez_Qx_right_top] = Ez_Flux_Calculation_Qx(i+1, j, D, B, dm.Deltax[i], dm.Deltay); //For i+1/2, j face
            
            auto [Ez_Qy_left_bottom, Ez_Qy_right_bottom] = Ez_Flux_Calculation_Qy(i, j, D, B, dm.Deltax[i], dm.Deltay); //For i, j-1/2 face
            auto[Ez_Qy_left_top, Ez_Qy_right_top] = Ez_Flux_Calculation_Qy(i, j+1, D, B, dm.Deltax[i], dm.Deltay); //For i, j+1/2 face
            
            auto [Hz_Fx_left_bottom, Hz_Fx_left_top] = Hz_Flux_Calculation_Fx(i, j, B, D, dm.Deltax[i], dm.Deltay); //For i-1/2, j face
            auto [Hz_Fx_right_bottom, Hz_Fx_right_top] = Hz_Flux_Calculation_Fx(i+1, j, B, D, dm.Deltax[i], dm.Deltay); //For i+1/2, j face

            auto [Hz_Fy_left_bottom, Hz_Fy_right_bottom] = Hz_Flux_Calculation_Fy(i, j, B, D, dm.Deltax[i], dm.Deltay); //For i, j-1/2 face
            auto[Hz_Fy_left_top, Hz_Fy_right_top] = Hz_Flux_Calculation_Fy(i, j+1, B, D, dm.Deltax[i], dm.Deltay); //For i, j+1/2 face

            //Used for calculating (curl(H) dot B - curl(E) dot D)_iplus1_j
            
            auto [Ez_Qx_right2_bottom, Ez_Qx_right2_top] = Ez_Flux_Calculation_Qx(i+2, j, D, B, dm.Deltax[i], dm.Deltay); // For i+3/2, j face
            auto [Hz_Qx_right2_bottom, Hz_Qx_right2_top] = Hz_Flux_Calculation_Fx(i+2, j, B, D, dm.Deltax[i], dm.Deltay); // For i+3/2, j face

            auto[Ez_Qy_left2_bottom, Ez_Qy_right2_bottom] = Ez_Flux_Calculation_Qy(i+1, j, D, B, dm.Deltax[i], dm.Deltay); //For i+1/2, j-1/2 face
            auto[Hz_Fy_left2_bottom, Hz_Fy_right2_bottom] = Hz_Flux_Calculation_Fy(i+1, j, B, D, dm.Deltax[i], dm.Deltay); //For i+1/2, j-1/2 face

            auto[Ez_Qy_left2_top, Ez_Qy_right2_top] = Ez_Flux_Calculation_Qy(i+1, j+1, D, B, dm.Deltax[i], dm.Deltay); //For i+1/2, j+1/2 face
            auto[Hz_Fy_left2_top, Hz_Fy_right2_top] = Hz_Flux_Calculation_Fy(i+1, j+1, B, D, dm.Deltax[i], dm.Deltay); //For i+1/2, j+1/2 face

            //Used for calculating (curl(H) dot B - curl(E) dot D)_i_jplus1

            auto [Ez_Qx_left3_bottom, Ez_Qx_left3_top] = Ez_Flux_Calculation_Qx(i, j+1, D, B, dm.Deltax[i], dm.Deltay); //For i-1/2, j+1
            auto [Hz_Fx_left3_bottom, Hz_Fx_left3_top] = Hz_Flux_Calculation_Fx(i, j+1, B, D, dm.Deltax[i], dm.Deltay); //For i-1/2, j+1

            auto [Ez_Qx_right3_bottom, Ez_Qx_right3_top] = Ez_Flux_Calculation_Qx(i+1, j+1, D, B, dm.Deltax[i], dm.Deltay); //For i+1/2, j+1
            auto [Hz_Fx_right3_bottom, Hz_Fx_right3_top] = Hz_Flux_Calculation_Fx(i+1, j+1, B, D, dm.Deltax[i], dm.Deltay); //For i+1/2, j+1

            auto[Ez_Qy_left3_bottom, Ez_Qy_right3_bottom] = Ez_Flux_Calculation_Qy(i, j+1, D, B, dm.Deltax[i], dm.Deltay); //For i, j+1/2 face
            auto[Hz_Fy_left3_bottom, Hz_Fy_right3_bottom] = Hz_Flux_Calculation_Fy(i, j+1, B, D, dm.Deltax[i], dm.Deltay); //For i, j+1/2 face

            auto[Ez_Qy_left3_top, Ez_Qy_right3_top] = Ez_Flux_Calculation_Qy(i, j+2, D, B, dm.Deltax[i], dm.Deltay); //For i, j+3/2 face
            auto[Hz_Fy_left3_top, Hz_Fy_right3_top] = Hz_Flux_Calculation_Fy(i, j+2, B, D, dm.Deltax[i], dm.Deltay); //For i, j+3/2 face

            //ExB calculations

            double E_cross_B_x_ij = E[1][i][j] * B[2][i][j] - E[2][i][j] * Compute_A_to_cell_center(B, 1, i, j);
            double E_cross_B_x_iplus1_j =  E[1][i+1][j] * B[2][i+1][j] - E[2][i+1][j] * Compute_A_to_cell_center(B, 1, i+1, j);

            double E_cross_B_y_ij = E[2][i][j] * Compute_A_to_cell_center(B, 0, i, j) - E[0][i][j] * Compute_A_to_cell_center(B, 2, i, j);
            double E_cross_B_y_i_jplus1 =  E[2][i][j+1] * Compute_A_to_cell_center(B, 0, i, j+1) - E[0][i][j+1] * Compute_A_to_cell_center(B, 2, i, j+1);

            double E_cross_B_z_ij = E[0][i][j] * Compute_A_to_cell_center(B, 1, i, j) + E[1][i][j] * Compute_A_to_cell_center(B, 0, i, j);

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

            double curl_E_y_bottom_ij = (Ez_Qy_right_bottom - Ez_Qy_left_bottom)/(dm.Deltax[i]);  //For i, j-1/2 face
            double curl_E_y_top_ij = (Ez_Qy_right_top - Ez_Qy_left_top)/(dm.Deltax[i]); //For i, j+1/2 face

            double curl_H_z_ij = Fz_reconstruction(i, j, B, D, dm.Deltax[i], dm.Deltay);
            double curl_E_z_ij = Qz_reconstruction(i, j, D, B, dm.Deltax[i], dm.Deltay);

            double curl_H_curl_E_ij = 0.5 * (curl_H_x_left_ij * Bx_avg_ij + curl_H_x_right_ij * Bx_avg_iplus1_j) + 0.5 * (curl_H_y_bottom_ij * By_avg_ij + curl_H_y_top_ij * By_avg_i_jplus1) + curl_H_z_ij * B[2][i][j] - 0.5 * (curl_E_x_left_ij * Dx_avg_ij + curl_E_x_right_ij * Dx_avg_iplus1_j) - 0.5 * (curl_E_y_bottom_ij * Dy_avg_ij + curl_E_y_top_ij * Dy_avg_i_jplus1) - curl_E_z_ij * D[2][i][j];

            //(curl(H) dot B - curl(E) dot D)_iplus1_j calculation

            double curl_H_x_left_iplus1_j = (Hz_Fx_right_top - Hz_Fx_right_bottom)/(dm.Deltay); //For i+1/2,j face
            double curl_H_x_right_iplus1_j = (Hz_Qx_right2_top - Hz_Qx_right2_bottom)/(dm.Deltay); //For i+3/2,j face

            double curl_E_x_left_iplus1_j = (Ez_Qx_right_top - Ez_Qx_right_bottom)/(dm.Deltay); //For i+1/2,j face
            double curl_E_x_right_iplus1_j =  (Ez_Qx_right2_top - Ez_Qx_right2_bottom)/(dm.Deltay); //For i+3/2,j face

            double curl_H_y_bottom_iplus1_j = (Hz_Fy_right2_bottom - Hz_Fy_left2_bottom)/(dm.Deltax[i]);  //For i+1, j-1/2 face
            double curl_H_y_top_iplus1_j = (Hz_Fy_right2_top - Hz_Fy_left2_top)/(dm.Deltax[i]); //For i+1, j+1/2 face

            double curl_E_y_bottom_iplus1_j = (Ez_Qy_right2_bottom - Ez_Qy_left2_bottom)/(dm.Deltax[i]);  //For i+1, j-1/2 face
            double curl_E_y_top_iplus1_j = (Ez_Qy_right2_top - Ez_Qy_left2_top)/(dm.Deltax[i]); //For i+1, j+1/2 face

            double curl_H_z_iplus1_j = Fz_reconstruction(i+1, j, B, D, dm.Deltax[i], dm.Deltay);
            double curl_E_z_iplus1_j = Qz_reconstruction(i+1, j, D, B, dm.Deltax[i], dm.Deltay);

            double curl_H_curl_E_iplus1_j = 0.5 * (curl_H_x_left_iplus1_j * Bx_avg_iplus1_j + curl_H_x_right_iplus1_j * Bx_avg_iplus2_j) + 0.5 * (curl_H_y_bottom_iplus1_j * By_avg_i_jplus1 + curl_H_y_top_iplus1_j * By_avg_i_jplus2) + curl_H_z_iplus1_j * B[2][i+1][j] - 0.5 * (curl_E_x_left_iplus1_j * Dx_avg_iplus1_j + curl_E_x_right_iplus1_j * Dx_avg_iplus2_j) - 0.5 * (curl_E_y_bottom_iplus1_j * Dy_avg_i_jplus1 + curl_E_y_top_iplus1_j * Dy_avg_i_jplus2) - curl_E_z_iplus1_j * D[2][i+1][j];

            //(curl(H) dot B - curl(E) dot D)_i_jplus1 calculation

            double curl_H_x_left_i_jplus1 = (Hz_Fx_left3_top - Hz_Fx_left3_bottom)/(dm.Deltay); //For i-1/2,j+1 face
            double curl_H_x_right_i_jplus1 = (Hz_Fx_right3_top - Hz_Fx_right3_bottom)/(dm.Deltay); //For i+1/2,j+1 face

            double curl_E_x_left_i_jplus1 = (Ez_Qx_left3_top - Ez_Qx_left3_bottom)/(dm.Deltay); //For i-1/2,j+1 face
            double curl_E_x_right_i_jplus1 =  (Ez_Qx_right3_top - Ez_Qx_right3_bottom)/(dm.Deltay); //For i+1/2,j+1 face

            double curl_H_y_bottom_i_jplus1 = (Hz_Fy_right3_bottom - Hz_Fy_left3_bottom)/(dm.Deltax[i]);  //For i, j+1/2 face
            double curl_H_y_top_i_jplus1 = (Hz_Fy_right3_top - Hz_Fy_left3_top)/(dm.Deltax[i]); //For i, j+3/2 face

            double curl_E_y_bottom_i_jplus1 = (Ez_Qy_right3_bottom - Ez_Qy_left3_bottom)/(dm.Deltax[i]);  //For i, j+1/2 face
            double curl_E_y_top_i_jplus1 = (Ez_Qy_right3_top - Ez_Qy_left3_top)/(dm.Deltax[i]); //For i, j+3/2 face

            double curl_H_z_i_jplus1 = Fz_reconstruction(i, j+1, B, D, dm.Deltax[i], dm.Deltay);
            double curl_E_z_i_jplus1 = Qz_reconstruction(i, j+1, D, B, dm.Deltax[i], dm.Deltay);

            double curl_H_curl_E_i_jplus1 = 0.5 * (curl_H_x_left_i_jplus1 * Bx_avg_i_jplus1 + curl_H_x_right_i_jplus1 * Bx_avg_i_jplus2) + 0.5 * (curl_H_y_bottom_i_jplus1 * By_avg_i_jplus1 + curl_H_y_top_i_jplus1 * By_avg_i_jplus2) + curl_H_z_i_jplus1 * B[2][i][j+1] - 0.5 * (curl_E_x_left_i_jplus1 * Dx_avg_i_jplus1 + curl_E_x_right_i_jplus1 * Dx_avg_i_jplus2) - 0.5 * (curl_E_y_bottom_i_jplus1 * Dy_avg_i_jplus1 + curl_E_y_top_i_jplus1 * Dy_avg_i_jplus2) - curl_E_z_i_jplus1 * D[2][i][j+1];

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

            if (dm.x[i] >= 0.9 * max_dm_x){
            
            Qx[i][j] = -(Compute_Damping_Term(i, dm) * D[0][i][j]) - curl_E_x_left_ij;
            Qy[i][j] = -(Compute_Damping_Term(i, dm) * D[1][i][j]) - curl_E_y_bottom_ij;
            Qz[i][j] = -(Compute_Damping_Term(i, dm) * D[2][i][j]) - curl_E_z_ij;

            Fx[i][j] = -(Compute_Damping_Term(i, dm) * B[0][i][j]) + curl_H_x_left_ij - Jx_avg;
            Fy[i][j] = -(Compute_Damping_Term(i, dm) * B[1][i][j]) + curl_H_y_bottom_ij - Jy_avg;
            Fz[i][j] = -(Compute_Damping_Term(i, dm) * B[2][i][j]) + curl_H_z_ij - Jz_ij;

            }

            else{

            Qx[i][j] = -curl_E_x_left_ij;
            Qy[i][j] = -curl_E_y_bottom_ij;
            Qz[i][j] = -curl_E_z_ij;

            Fx[i][j] = curl_H_x_left_ij - Jx_avg;
            Fy[i][j] = curl_H_y_bottom_ij - Jy_avg;
            Fz[i][j] = curl_H_z_ij - Jz_ij;
            }
            
    }

    

    return;
}
}
