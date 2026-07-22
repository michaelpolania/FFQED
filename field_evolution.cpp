#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <utility>

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
            
            }}
return;
                
}

/*
std::pair<double, double> Ez_Flux_Calculation(int i, int j, VectorField & D, VectorField & B, std::vector<double> & Deltax, double Deltay){

    //Computes minmod of centered, backward, and forward differences for x at i,j/i-1,j/i,j-1/i-1,j-1/i,j+1/i-1,j+1
    double slope_x_ij = minmod((D[2][i+1][j] - D[2][i-1][j]) / (2 * Deltax[i]),  2 * (D[2][i][j] - D[2][i-1][j]) / Deltax[i], 2 * (D[2][i+1][j] - D[2][i][j]) / Deltax[i]);
    double slope_x_iminus1_j = minmod((D[2][i][j] - D[2][i-2][j]) / (2 * Deltax[i]),  2 * (D[2][i-1][j] - D[2][i-2][j]) / Deltax[i], 2 * (D[2][i][j] - D[2][i-1][j]) / Deltax[i]);
    double slope_x_i_jminus1 = minmod((D[2][i+1][j-1] - D[2][i-1][j-1]) / (2 * Deltax[i]),  2 * (D[2][i][j-1] - D[2][i-1][j-1]) / Deltax[i], 2 * (D[2][i+1][j-1] - D[2][i][j-1]) / Deltax[i]);
    double slope_x_ij_minus1 = minmod((D[2][i][j-1] - D[2][i-2][j-1]) / (2 * Deltax[i]),  2 * (D[2][i-1][j-1] - D[2][i-2][j-1]) / Deltax[i], 2 * (D[2][i][j-1] - D[2][i-1][j-1]) / Deltax[i]);
    double slope_x_i_jplus1 = minmod((D[2][i+1][j+1] - D[2][i-1][j+1]) / (2 * Deltax[i]),  2 * (D[2][i][j+1] - D[2][i-1][j+1]) / Deltax[i], 2 * (D[2][i+1][j+1] - D[2][i][j+1]) / Deltax[i]);
    double slope_x_iminus1_jplus1 = minmod((D[2][i][j+1] - D[2][i-2][j+1]) / (2 * Deltax[i]),  2 * (D[2][i-1][j+1] - D[2][i-2][j+1]) / Deltax[i], 2 * (D[2][i][j+1] - D[2][i-1][j+1]) / Deltax[i]);

    //Computes minmod of centered, backward, and forward differences for y at i,j/i-1,j/i,j-1/i-1,j-1/i,j+1
    double slope_y_ij = minmod((D[2][i][j+1] - D[2][i][j-1]) / (2 * Deltay),  2 * (D[2][i][j] - D[2][i][j-1]) / Deltay, 2 * (D[2][i][j+1] - D[2][i][j]) / Deltay);
    double slope_y_iminus1_j = minmod((D[2][i-1][j+1] - D[2][i-1][j-1]) / (2 * Deltay),  2 * (D[2][i-1][j] - D[2][i-1][j-1]) / Deltay, 2 * (D[2][i-1][j+1] - D[2][i-1][j]) / Deltay);
    double slope_y_i_jminus1 = minmod((D[2][i][j] - D[2][i][j-2]) / (2 * Deltay),  2 * (D[2][i][j-1] - D[2][i][j-2]) / Deltay, 2 * (D[2][i][j] - D[2][i][j-1]) / Deltay);
    double slope_y_ij_minus1 = minmod((D[2][i-1][j] - D[2][i-1][j-2]) / (2 * Deltay),  2 * (D[2][i-1][j-1] - D[2][i-1][j-2]) / Deltay, 2 * (D[2][i-1][j] - D[2][i-1][j-1]) / Deltay);
    double slope_y_i_jplus1 = minmod((D[2][i][j+2] - D[2][i][j]) / (2 * Deltay),  2 * (D[2][i][j+1] - D[2][i][j]) / Deltay, 2 * (D[2][i][j+2] - D[2][i][j+1]) / Deltay);

    //Fluxes for Ez{i,j} aka Ez_{i-1/2, j-1/2}
    //G_1{i,j-1/2} calculation (G_1)

    double G_1_up = D[2][i][j] - 0.5 * slope_y_ij * Deltay;
    double G_1_down = D[2][i][j-1] + 0.5 * slope_y_i_jminus1 * Deltay;
    double G_1 = 0.5 * (G_1_up + G_1_down) - 0.5 * (B[0][i][j] - B[0][i][j-1]);

    //G_1{i-1,j-1/2} calculation (G_2)

    double G_2_up = D[2][i-1][j] - 0.5 * slope_y_iminus1_j * Deltay;
    double G_2_down = D[2][i-1][j-1] + 0.5 * slope_y_ij_minus1 * Deltay;
    double G_2 = 0.5 * (G_2_up + G_2_down) - 0.5 * (B[0][i][j] - B[0][i][j-1]);

    //F_1{i-1/2,j} calculation (F_1)

    double F_1_up = D[2][i][j] - 0.5 * slope_x_ij * Deltax[i];
    double F_1_down = D[2][i-1][j] + 0.5 * slope_x_iminus1_j * Deltax[i];
    double F_1 = 0.5 * (F_1_up + F_1_down) - 0.5 * (B[1][i][j] - B[1][i-1][j]);

    //F_1{i-1/2,j-1} calculation (F_2)

    double F_2_up = D[2][i][j-1] - 0.5 * slope_x_i_jminus1 * Deltax[i];
    double F_2_down = D[2][i-1][j-1] + 0.5 * slope_x_ij_minus1 * Deltax[i];
    double F_2 = 0.5 * (F_2_up + F_2_down) - 0.5 * (B[1][i][j] - B[1][i-1][j]);

    //Fluxes for Ez_{i-1/2,j+1/2}
    //G_3{i,j+1/2} calculation (G_3)
    double G_3_up = D[2][i][j+1] - 0.5 * slope_y_i_jplus1 * Deltay;
    double G_3_down = D[2][i][j] + 0.5 * slope_y_ij * Deltay;
    double G_3 = 0.5 * (G_3_up + G_3_down) - 0.5 * (B[0][i][j+1] - B[0][i][j]);

    //G_4{i-1,j+1/2} calculation (G_4)

    double G_4_up = D[2][i-1][j] - 0.5 * slope_y_iminus1_j * Deltay;
    double G_4_down = D[2][i-1][j-1] + 0.5 * slope_y_ij_minus1 * Deltay;
    double G_4 = 0.5 * (G_4_up + G_4_down) - 0.5 * (B[0][i][j+1] - B[0][i][j]);

    //F_3{i-1/2,j} calculation (F_3)

    double F_3_up = D[2][i][j] - 0.5 * slope_x_ij * Deltax[i];
    double F_3_down = D[2][i-1][j] + 0.5 * slope_x_iminus1_j * Deltax[i];
    double F_3 = 0.5 * (F_3_up + F_3_down) - 0.5 * (B[1][i][j] - B[1][i-1][j]);

    //F_4{i-1/2,j+1} calculation (F_4)

    double F_4_up = D[2][i][j+1] - 0.5 * slope_x_i_jplus1 * Deltax[i]; 
    double F_4_down = D[2][i-1][j+1] + 0.5 * slope_x_iminus1_jplus1 * Deltax[i]; 
    double F_4 = 0.5 * (F_4_up + F_4_down) - 0.5 * (B[1][i][j] - B[1][i-1][j]);

    double Ez_ij = 0.25 * (-F_1 - F_2 + G_1 + G_2);
    double Ez_i_jplus1 = 0.25 * (-F_3 - F_4 + G_3 + G_4);

    return {Ez_ij, Ez_i_jplus1};
}

*/

std::pair<double, double> Ez_Flux_Calculation_Qx(int i, int j, VectorField & E, VectorField & B, std::vector<double> & Deltax, double Deltay){

    //Compute Bx at cell centers by averaging the two Bx face values at i,j/i,j+1/i,j-1/i,j-2/i-1,j/i-1,j-1
    double Bx_ij = 0.5 * (B[0][i][j] + B[0][i+1][j]);
    double Bx_i_jplus1 = 0.5 * (B[0][i][j+1] + B[0][i+1][j+1]);
    double Bx_i_jminus1 = 0.5 * (B[0][i][j-1] + B[0][i+1][j-1]);
    double Bx_i_jminus2 = 0.5 * (B[0][i][j-2] + B[0][i+1][j-2]);
    double Bx_iminus1_j = 0.5 * (B[0][i][j] + B[0][i-1][j]);
    double Bx_iminus1_jplus1 = 0.5 * (B[0][i][j+1] + B[0][i-1][j+1]);
    double Bx_ij_minus1 = 0.5 * (B[0][i][j-1] + B[0][i-1][j-1]);

    double slope_Bx_y_j = minmod((B[0][i][j+1] - B[0][i][j-1]) / (2 * Deltay),  2 * (B[0][i][j] - B[0][i][j-1]) / Deltay, 2 * (B[0][i][j+1] - B[0][i][j]) / Deltay);
    double slope_Bx_y_i_jminus1 = minmod((B[0][i][j] - B[0][i][j-2]) / (2 * Deltay),  2 * (B[0][i][j-1] - B[0][i][j-2]) / Deltay, 2 * (B[0][i][j] - B[0][i][j-1]) / Deltay);
    double slope_Bx_y_i_jplus1 = minmod((B[0][i][j+2] - B[0][i][j]) / (2 * Deltay),  2 * (B[0][i][j+1] - B[0][i][j]) / Deltay, 2 * (B[0][i][j+2] - B[0][i][j+1]) / Deltay);
    double slope_Bx_y_iminus1_j = minmod((B[0][i-1][j+1] - B[0][i-1][j-1]) / (2 * Deltay),  2 * (B[0][i-1][j] - B[0][i-1][j-1]) / Deltay, 2 * (B[0][i-1][j+1] - B[0][i-1][j]) / Deltay);
    double slope_Bx_y_ij_minus1 = minmod((B[0][i-1][j] - B[0][i-1][j-2]) / (2 * Deltay),  2 * (B[0][i-1][j-1] - B[0][i-1][j-2]) / Deltay, 2 * (B[0][i-1][j] - B[0][i-1][j-1]) / Deltay);
    double slope_Bx_y_iminus1_jplus1 = minmod((B[0][i-1][j+2] - B[0][i-1][j]) / (2 * Deltay),  2 * (B[0][i-1][j+1] - B[0][i-1][j]) / Deltay, 2 * (B[0][i-1][j+2] - B[0][i-1][j+1]) / Deltay);

    double slope_Ez_y_j = minmod((E[2][i][j+1] - E[2][i][j-1]) / (2 * Deltay),  2 * (E[2][i][j] - E[2][i][j-1]) / Deltay, 2 * (E[2][i][j+1] - E[2][i][j]) / Deltay);
    double slope_Ez_y_i_jminus1 = minmod((E[2][i][j] - E[2][i][j-2]) / (2 * Deltay),  2 * (E[2][i][j-1] - E[2][i][j-2]) / Deltay, 2 * (E[2][i][j] - E[2][i][j-1]) / Deltay);
    double slope_Ez_y_i_jplus1 = minmod((E[2][i][j+2] - E[2][i][j]) / (2 * Deltay),  2 * (E[2][i][j+1] - E[2][i][j]) / Deltay, 2 * (E[2][i][j+2] - E[2][i][j+1]) / Deltay);
    double slope_Ez_y_iminus1_j = minmod((E[2][i-1][j+1] - E[2][i-1][j-1]) / (2 * Deltay),  2 * (E[2][i-1][j] - E[2][i-1][j-1]) / Deltay, 2 * (E[2][i-1][j+1] - E[2][i-1][j]) / Deltay);
    double slope_Ez_y_ij_minus1 = minmod((E[2][i-1][j] - E[2][i-1][j-2]) / (2 * Deltay),  2 * (E[2][i-1][j-1] - E[2][i-1][j-2]) / Deltay, 2 * (E[2][i-1][j] - E[2][i-1][j-1]) / Deltay);
    double slope_Ez_y_iminus1_jplus1 = minmod((E[2][i-1][j+2] - E[2][i-1][j]) / (2 * Deltay),  2 * (E[2][i-1][j+1] - E[2][i-1][j]) / Deltay, 2 * (E[2][i-1][j+2] - E[2][i-1][j+1]) / Deltay);

    
    double slope_Ez_x_ij = minmod((E[2][i+1][j] - E[2][i-1][j]) / (2 * Deltax[i]),  2 * (E[2][i][j] - E[2][i-1][j]) / Deltax[i], 2 * (E[2][i+1][j] - E[2][i][j]) / Deltax[i]);
    double slope_Ez_x_iminus1_j = minmod((E[2][i][j] - E[2][i-2][j]) / (2 * Deltax[i]),  2 * (E[2][i-1][j] - E[2][i-2][j]) / Deltax[i], 2 * (E[2][i][j] - E[2][i-1][j]) / Deltax[i]);
    double slope_Ez_x_i_jminus1 = minmod((E[2][i+1][j-1] - E[2][i-1][j-1]) / (2 * Deltax[i]),  2 * (E[2][i][j-1] - E[2][i-1][j-1]) / Deltax[i], 2 * (E[2][i+1][j-1] - E[2][i][j-1]) / Deltax[i]);
    double slope_Ez_x_ij_minus1 = minmod((E[2][i][j-1] - E[2][i-2][j-1]) / (2 * Deltax[i]),  2 * (E[2][i-1][j-1] - E[2][i-2][j-1]) / Deltax[i], 2 * (E[2][i][j-1] - E[2][i-1][j-1]) / Deltax[i]);
    double slope_Ez_x_i_jplus1 = minmod((E[2][i+1][j+1] - E[2][i-1][j+1]) / (2 * Deltax[i]),  2 * (E[2][i][j+1] - E[2][i-1][j+1]) / Deltax[i], 2 * (E[2][i+1][j+1] - E[2][i][j+1]) / Deltax[i]);
    double slope_Ez_x_iminus1_jplus1 = minmod((E[2][i][j+1] - E[2][i-2][j+1]) / (2 * Deltax[i]),  2 * (E[2][i-1][j+1] - E[2][i-2][j+1]) / Deltax[i], 2 * (E[2][i][j+1] - E[2][i-1][j+1]) / Deltax[i]);
    
    //Fluxes for Ez{i,j} aka Ez_{i-1/2, j-1/2}
    //G_1{i,j-1/2} calculation (G_1)

    double Bx_ij_up = Bx_ij - 0.5 * Deltay * slope_Bx_y_j;
    double Bx_i_jminus1_down = Bx_i_jminus1 + 0.5 * Deltay * slope_Bx_y_i_jminus1;

    //Note: Ez is positive in the G matrix in Yu paper
    double G_1_up = E[2][i][j] - 0.5 * slope_Ez_y_j * Deltay;
    double G_1_down = E[2][i][j-1] + 0.5 * slope_Ez_y_i_jminus1 * Deltay;
    double G_1 = 0.5 * (G_1_up + G_1_down) - 0.5 * (Bx_ij_up - Bx_i_jminus1_down);

    //G_1{i-1,j-1/2} calculation (G_2)
    
    double Bx_iminus1_j_up = Bx_iminus1_j - 0.5 * Deltay * slope_Bx_y_iminus1_j;
    double Bx_ijminus1_down = Bx_ij_minus1 + 0.5 * Deltay * slope_Bx_y_ij_minus1;

    //Note: Ez is positive in the G matrix in Yu paper
    double G_2_up = E[2][i-1][j] - 0.5 * slope_Ez_y_iminus1_j * Deltay;
    double G_2_down = E[2][i-1][j-1] + 0.5 * slope_Ez_y_ij_minus1 * Deltay;
    double G_2 = 0.5 * (G_2_up + G_2_down) - 0.5 * (Bx_iminus1_j_up - Bx_ijminus1_down);

    //F_2{i-1/2,j} calculation (F_1)

    double F_1_up = E[2][i][j] - 0.5 * slope_Ez_x_ij * Deltax[i];
    double F_1_down = E[2][i-1][j] + 0.5 * slope_Ez_x_iminus1_j  * Deltax[i];

    //Note: Ez is negative in the F matrix in Yu paper and no dissipation term since Bx lives w/ F_1

    double F_1 = 0.5 * (-F_1_up - F_1_down);

    //F_2{i-1/2,j-1} calculation (F_2)

    double F_2_up = E[2][i][j-1] - 0.5 * slope_Ez_x_i_jminus1 * Deltax[i];
    double F_2_down = E[2][i-1][j-1] + 0.5 * slope_Ez_x_ij_minus1  * Deltax[i];

    //Note: Ez is negative in the F matrix in Yu paper and no dissipation term since Bx lives w/ F_1

    double F_2 = 0.5 * (-F_2_up - F_2_down);

    //Fluxes for Ez{i,j+1} aka Ez_{i-1, j+1/2}
    //G_1{i,j+1/2} calculation (G_3)

    double Bx_i_jplus1_up = Bx_i_jplus1 - 0.5 * Deltay * slope_Bx_y_i_jplus1;
    double Bx_ij_down = Bx_ij + 0.5 * Deltay * slope_Bx_y_j;

    //Note: Ez is positive in the G matrix in Yu paper
    double G_3_up = E[2][i][j+1] - 0.5 * slope_Ez_y_i_jplus1 * Deltay;
    double G_3_down = E[2][i][j] + 0.5 * slope_Ez_y_j * Deltay;
    double G_3 = 0.5 * (G_3_up + G_3_down) - 0.5 * (Bx_i_jplus1_up - Bx_ij_down);
    
    //G_1{i-1,j+1/2} calculation (G_4)

    double Bx_iminus1_jplus1_up = Bx_iminus1_jplus1 - 0.5 * Deltay * slope_Bx_y_iminus1_jplus1;
    double Bx_iminus1_j_down = Bx_iminus1_j + 0.5 * Deltay * slope_Bx_y_iminus1_j;

    //Note: Ez is positive in the G matrix in Yu paper

    double G_4_up = E[2][i-1][j+1] - 0.5 * Deltay * slope_Ez_y_iminus1_jplus1;
    double G_4_down = E[2][i-1][j] + 0.5 * Deltay * slope_Ez_y_iminus1_j;

    double G_4 = 0.5 * (G_4_up + G_4_down) - 0.5 * (Bx_iminus1_jplus1_up - Bx_iminus1_j_down);

    //F_2{i-1/2,j} calculation (F_3) same thing as F_1

    //Note: Ez is negative in the F matrix in Yu paper and no dissipation term since Bx lives w/ F_1

    double F_3 = 0.5 * (-F_1_up - F_1_down);

    //F_2{i-1/2,j+1} calculation (F_4)

    double F_4_up = E[2][i][j+1] - 0.5 * Deltax[i] * slope_Ez_x_i_jplus1;
    double F_4_down = E[2][i-1][j+1] + 0.5 * Deltax[i] * slope_Ez_x_iminus1_jplus1;

    double F_4 = 0.5 * (-F_4_up - F_4_down);

    //Toth 2000 Eq(19)
    double Ez_ij = 0.25 * (-F_1 - F_2 + G_1 + G_2);
    double Ez_i_jplus1 = 0.25 * (-F_3 - F_4 + G_3 + G_4);

    return {Ez_ij, Ez_i_jplus1};
}

std::pair<double, double> Ez_Flux_Calculation_Qy(int i, int j, VectorField & E, VectorField & B, std::vector<double> & Deltax, double Deltay){

    //Compute By at cell centers by averaging the two By face values
    double By_ij = 0.5 * (B[1][i][j] + B[1][i][j+1]);
    double By_iminus1_j = 0.5 * (B[1][i-1][j] + B[1][i-1][j+1]);
    double By_iminus1_jminus1 = 0.5 * (B[1][i-1][j] + B[1][i-1][j-1]);
    double By_i_jminus1 = 0.5 * (B[1][i][j] + B[1][i][j-1]);
    double By_iplus1_jminus1 = 0.5 * (B[1][i+1][j] + B[1][i+1][j-1]);
    double By_iplus1_j = 0.5 * (B[1][i+1][j] + B[1][i+1][j+1]);

    double slope_By_x_ij = minmod((B[1][i+1][j] - B[1][i-1][j]) / (2 * Deltax[i]),  2 * (B[1][i][j] - B[1][i-1][j]) / Deltax[i], 2 * (B[1][i+1][j] - B[1][i][j]) / Deltax[i]);
    double slope_By_x_iminus1_j =minmod((B[1][i][j] - B[1][i-2][j]) / (2 * Deltax[i]),  2 * (B[1][i-1][j] - B[1][i-2][j]) / Deltax[i], 2 * (B[1][i][j] - B[1][i-1][j]) / Deltax[i]);
    double slope_By_x_i_jminus1 = minmod((B[1][i+1][j-1] - B[1][i-1][j-1]) / (2 * Deltax[i]),  2 * (B[1][i][j-1] - B[1][i-1][j-1]) / Deltax[i], 2 * (B[1][i+1][j-1] - B[1][i][j-1]) / Deltax[i]);
    double slope_By_x_iminus1_jminus1 = minmod((B[1][i][j-1] - B[1][i-2][j-1]) / (2 * Deltax[i]),  2 * (B[1][i-1][j-1] - B[1][i-2][j-1]) / Deltax[i], 2 * (B[1][i][j-1] - B[1][i-1][j-1]) / Deltax[i]);
    double slope_By_x_iplus1_jminus1 = minmod((B[1][i+2][j-1] - B[1][i][j-1]) / (2 * Deltax[i]),  2 * (B[1][i+1][j-1] - B[1][i][j-1]) / Deltax[i], 2 * (B[1][i+2][j-1] - B[1][i+1][j-1]) / Deltax[i]);
    double slope_By_x_iplus1_j =  minmod((B[1][i+2][j] - B[1][i][j]) / (2 * Deltax[i]),  2 * (B[1][i+1][j] - B[1][i][j]) / Deltax[i], 2 * (B[1][i+2][j] - B[1][i+1][j]) / Deltax[i]);

    double slope_Ez_y_j = minmod((E[2][i][j+1] - E[2][i][j-1]) / (2 * Deltay),  2 * (E[2][i][j] - E[2][i][j-1]) / Deltay, 2 * (E[2][i][j+1] - E[2][i][j]) / Deltay);
    double slope_Ez_y_i_jminus1 = minmod((E[2][i][j] - E[2][i][j-2]) / (2 * Deltay),  2 * (E[2][i][j-1] - E[2][i][j-2]) / Deltay, 2 * (E[2][i][j] - E[2][i][j-1]) / Deltay);
    double slope_Ez_y_i_jplus1 = minmod((E[2][i][j+2] - E[2][i][j]) / (2 * Deltay),  2 * (E[2][i][j+1] - E[2][i][j]) / Deltay, 2 * (E[2][i][j+2] - E[2][i][j+1]) / Deltay);
    double slope_Ez_y_iminus1_j = minmod((E[2][i-1][j+1] - E[2][i-1][j-1]) / (2 * Deltay),  2 * (E[2][i-1][j] - E[2][i-1][j-1]) / Deltay, 2 * (E[2][i-1][j+1] - E[2][i-1][j]) / Deltay);
    double slope_Ez_y_ij_minus1 = minmod((E[2][i-1][j] - E[2][i-1][j-2]) / (2 * Deltay),  2 * (E[2][i-1][j-1] - E[2][i-1][j-2]) / Deltay, 2 * (E[2][i-1][j] - E[2][i-1][j-1]) / Deltay);
    double slope_Ez_y_iminus1_jplus1 = minmod((E[2][i-1][j+2] - E[2][i-1][j]) / (2 * Deltay),  2 * (E[2][i-1][j+1] - E[2][i-1][j]) / Deltay, 2 * (E[2][i-1][j+2] - E[2][i-1][j+1]) / Deltay);
    double slope_Ez_y_iplus1_j = minmod((E[2][i+1][j+1] - E[2][i+1][j-1]) / (2 * Deltay),  2 * (E[2][i+1][j] - E[2][i+1][j-1]) / Deltay, 2 * (E[2][i+1][j+1] - E[2][i+1][j]) / Deltay);
    double slope_Ez_y_iplus1_jminus1 = minmod((E[2][i+1][j] - E[2][i+1][j-2]) / (2 * Deltay),  2 * (E[2][i+1][j-1] - E[2][i+1][j-2]) / Deltay, 2 * (E[2][i+1][j] - E[2][i+1][j-1]) / Deltay);

    double slope_Ez_x_ij = minmod((E[2][i+1][j] - E[2][i-1][j]) / (2 * Deltax[i]),  2 * (E[2][i][j] - E[2][i-1][j]) / Deltax[i], 2 * (E[2][i+1][j] - E[2][i][j]) / Deltax[i]);
    double slope_Ez_x_iminus1_j = minmod((E[2][i][j] - E[2][i-2][j]) / (2 * Deltax[i]),  2 * (E[2][i-1][j] - E[2][i-2][j]) / Deltax[i], 2 * (E[2][i][j] - E[2][i-1][j]) / Deltax[i]);
    double slope_Ez_x_i_jminus1 = minmod((E[2][i+1][j-1] - E[2][i-1][j-1]) / (2 * Deltax[i]),  2 * (E[2][i][j-1] - E[2][i-1][j-1]) / Deltax[i], 2 * (E[2][i+1][j-1] - E[2][i][j-1]) / Deltax[i]);
    double slope_Ez_x_ij_minus1 = minmod((E[2][i][j-1] - E[2][i-2][j-1]) / (2 * Deltax[i]),  2 * (E[2][i-1][j-1] - E[2][i-2][j-1]) / Deltax[i], 2 * (E[2][i][j-1] - E[2][i-1][j-1]) / Deltax[i]);
    double slope_Ez_x_i_jplus1 = minmod((E[2][i+1][j+1] - E[2][i-1][j+1]) / (2 * Deltax[i]),  2 * (E[2][i][j+1] - E[2][i-1][j+1]) / Deltax[i], 2 * (E[2][i+1][j+1] - E[2][i][j+1]) / Deltax[i]);
    double slope_Ez_x_iminus1_jplus1 = minmod((E[2][i][j+1] - E[2][i-2][j+1]) / (2 * Deltax[i]),  2 * (E[2][i-1][j+1] - E[2][i-2][j+1]) / Deltax[i], 2 * (E[2][i][j+1] - E[2][i-1][j+1]) / Deltax[i]);
    double slope_Ez_x_iplus1_jminus1 = minmod((E[2][i+2][j-1] - E[2][i][j-1]) / (2 * Deltax[i]),  2 * (E[2][i+1][j-1] - E[2][i][j-1]) / Deltax[i], 2 * (E[2][i+2][j-1] - E[2][i+1][j-1]) / Deltax[i]);
    double slope_Ez_x_iplus1_j = minmod((E[2][i+2][j] - E[2][i][j]) / (2 * Deltax[i]),  2 * (E[2][i+1][j] - E[2][i][j]) / Deltax[i], 2 * (E[2][i+2][j] - E[2][i+1][j]) / Deltax[i]);

    //Fluxes for Ez{i,j} aka Ez_{i-1/2, j-1/2}
    //G_1{i,j-1/2} calculation (G_1)

    //Note: Ez is positive in the G matrix in Yu paper and no dissipation term since By lives w/ G_1
    double G_1_up = E[2][i][j] - 0.5 * slope_Ez_y_j * Deltay;
    double G_1_down = E[2][i][j-1] + 0.5 * slope_Ez_y_i_jminus1 * Deltay;
    double G_1 = 0.5 * (G_1_up + G_1_down);

    //G_1{i-1,j-1/2} calculation (G_2)

    //Note: Ez is positive in the G matrix in Yu paper and no dissipation term since By lives w/ G_1
    double G_2_up = E[2][i-1][j] - 0.5 * slope_Ez_y_iminus1_j * Deltay;
    double G_2_down = E[2][i-1][j-1] + 0.5 * slope_Ez_y_ij_minus1 * Deltay;
    double G_2 = 0.5 * (G_2_up + G_2_down);

    //F_2{i-1/2,j} calculation (F_1)

    double By_ij_up = By_ij - 0.5 * Deltax[i] * slope_By_x_ij;
    double By_iminus1_j_down = By_iminus1_j + 0.5 * Deltax[i] * slope_By_x_iminus1_j;

    double F_1_up = E[2][i][j] - 0.5 * slope_Ez_x_ij * Deltax[i];
    double F_1_down = E[2][i-1][j] + 0.5 * slope_Ez_x_iminus1_j  * Deltax[i];

    //Note: Ez is negative in the F matrix in Yu paper

    double F_1 = 0.5 * (-F_1_up - F_1_down) - 0.5 * (By_ij_up - By_iminus1_j_down);

    //F_2{i-1/2,j-1} calculation (F_2)

    double By_i_jminus1_up = By_i_jminus1 - 0.5 * Deltax[i] * slope_By_x_i_jminus1;
    double By_iminus1_jminus1_down = By_iminus1_jminus1 + 0.5 * Deltax[i] * slope_By_x_iminus1_jminus1;

    double F_2_up = E[2][i][j-1] - 0.5 * slope_Ez_x_i_jminus1 * Deltax[i];
    double F_2_down = E[2][i-1][j-1] + 0.5 * slope_Ez_x_ij_minus1  * Deltax[i];

    double F_2 = 0.5 * (-F_2_up - F_2_down) - 0.5 * (By_i_jminus1_up - By_iminus1_jminus1_down);

    //E_z Fluxes for Ez{i+1,j} aka Ez_{i+1/2, j-1/2}

    //G_1{i,j-1/2} calculation (G_3), same as G_1

    //Note: Ez is positive in the G matrix in Yu paper and no dissipation term since By lives w/ G_1

    double G_3 = 0.5 * (G_1_up + G_1_down);

    //G_1{i+1, j-1/2}  calculation (G_4)

    double G_4_up = E[2][i+1][j] - 0.5 * Deltay * slope_Ez_y_iplus1_j;
    double G_4_down = E[2][i+1][j-1] + 0.5 * Deltay * slope_Ez_y_iplus1_jminus1;

    double G_4 = 0.5 * (G_4_up + G_4_down);

    //F_1{i+1/2, j-1} calculation (F_3)

    double By_iplus1_jminus1_up = By_iplus1_jminus1 - 0.5 * Deltax[i] * slope_By_x_iplus1_jminus1;  
    double By_i_jminus1_down = By_i_jminus1 + 0.5 * Deltax[i] * slope_By_x_i_jminus1;
    
    double F_3_up = E[2][i+1][j-1] - 0.5 * slope_Ez_x_iplus1_jminus1 * Deltax[i];
    double F_3_down = E[2][i][j-1] + 0.5 * slope_Ez_x_i_jminus1  * Deltax[i];

    double F_3 = 0.5 * (F_3_up + F_3_down) - 0.5 * (By_iplus1_jminus1_up - By_i_jminus1_down);

    //F_1{i+1/2,j} calculation (F_4)

    double By_iplus1_j_up = By_iplus1_j - 0.5 * Deltax[i] * slope_By_x_iplus1_j;
    double By_ij_down = By_ij + 0.5 * Deltax[i] * slope_By_x_ij;

    double F_4_up = E[2][i+1][j] - 0.5 * slope_Ez_x_iplus1_j * Deltax[i];
    double F_4_down = E[2][i][j] + 0.5 * slope_Ez_x_ij  * Deltax[i];

    double F_4 = 0.5 * (F_4_up + F_4_down) - 0.5 * (By_iplus1_j_up - By_ij_down);

    double Ez_ij = 0.25 * (-F_1 - F_2 + G_1 + G_2);
    double Ez_iplus1_j = 0.25 * (-F_3 - F_4 + G_3 + G_4);

    return {Ez_ij, Ez_iplus1_j};

}

std::tuple<double, double, double, double> Ez_Flux_Calculation_Qz(int i, int j, VectorField & E, VectorField & B, std::vector<double> & Deltax, double Deltay){

    double slope_Ex_y_ij = minmod((E[0][i][j+1] - E[0][i][j-1]) / (2 * Deltay),  2 * (E[0][i][j] - E[0][i][j-1]) / Deltay, 2 * (E[0][i][j+1] - E[0][i][j]) / Deltay);
    double slope_Ex_y_i_jminus1 = minmod((E[0][i][j] - E[0][i][j-2]) / (2 * Deltay),  2 * (E[0][i][j-1] - E[0][i][j-2]) / Deltay, 2 * (E[0][i][j] - E[0][i][j-1]) / Deltay);
    double slope_Ex_y_iminus1_j = minmod((E[0][i-1][j+1] - E[0][i-1][j-1]) / (2 * Deltay),  2 * (E[0][i-1][j] - E[0][i-1][j-1]) / Deltay, 2 * (E[0][i-1][j+1] - E[0][i-1][j]) / Deltay);
    double slope_Ex_y_iminus1_jminus1 = minmod((E[0][i-1][j] - E[0][i-1][j-2]) / (2 * Deltay),  2 * (E[0][i-1][j-1] - E[0][i-1][j-2]) / Deltay, 2 * (E[0][i-1][j] - E[0][i-1][j-1]) / Deltay);
   
    double slope_Ey_x_ij = minmod((E[1][i+1][j] - E[1][i-1][j]) / (2 * Deltax[i]),  2 * (E[1][i][j] - E[1][i-1][j]) / Deltax[i], 2 * (E[1][i+1][j] - E[1][i][j]) / Deltax[i]);
    double slope_Ey_x_iminus1_j = minmod((E[1][i][j] - E[1][i-2][j]) / (2 * Deltax[i]),  2 * (E[1][i-1][j] - E[1][i-2][j]) / Deltax[i], 2 * (E[1][i][j] - E[1][i-1][j]) / Deltax[i]);
    double slope_Ey_x_i_jminus1 = minmod((E[1][i+1][j-1] - E[1][i-1][j-1]) / (2 * Deltax[i]),  2 * (E[1][i][j-1] - E[1][i-1][j-1]) / Deltax[i], 2 * (E[1][i+1][j-1] - E[1][i][j-1]) / Deltax[i]);
    double slope_Ey_x_iminus1_jminus1 = minmod((E[1][i][j-1] - E[1][i-2][j-1]) / (2 * Deltax[i]),  2 * (E[1][i-1][j-1] - E[1][i-2][j-1]) / Deltax[i], 2 * (E[1][i][j-1] - E[1][i-1][j-1]) / Deltax[i]);

    double slope_Bz_y_ij = minmod((B[2][i][j+1] - B[2][i][j-1]) / (2 * Deltay),  2 * (B[2][i][j] - B[2][i][j-1]) / Deltay, 2 * (B[2][i][j+1] - B[2][i][j]) / Deltay);
    double slope_Bz_y_i_jminus1 = minmod((B[2][i][j] - B[2][i][j-2]) / (2 * Deltay),  2 * (B[2][i][j-1] - B[2][i][j-2]) / Deltay, 2 * (B[2][i][j] - B[2][i][j-1]) / Deltay);
    double slope_Bz_y_iminus1_j = minmod((B[2][i-1][j+1] - B[2][i-1][j-1]) / (2 * Deltay),  2 * (B[2][i-1][j] - B[2][i-1][j-1]) / Deltay, 2 * (B[2][i-1][j+1] - B[2][i-1][j]) / Deltay); 
    double slope_Bz_y_iminus1_jminus1 = minmod((B[2][i-1][j] - B[2][i-1][j-2]) / (2 * Deltay),  2 * (B[2][i-1][j-1] - B[2][i-1][j-2]) / Deltay, 2 * (B[2][i-1][j] - B[2][i-1][j-1]) / Deltay);
    
    double slope_Bz_x_ij = minmod((B[2][i+1][j] - B[2][i-1][j]) / (2 * Deltax[i]),  2 * (B[2][i][j] - B[2][i-1][j]) / Deltax[i], 2 * (B[2][i+1][j] - B[2][i][j]) / Deltax[i]);
    double slope_Bz_x_iminus1_j = minmod((B[2][i][j] - B[2][i-2][j]) / (2 * Deltax[i]),  2 * (B[2][i-1][j] - B[2][i-2][j]) / Deltax[i], 2 * (B[2][i][j] - B[2][i-1][j]) / Deltax[i]);
    double slope_Bz_x_i_jminus1 = minmod((B[2][i+1][j-1] - B[2][i-1][j-1]) / (2 * Deltax[i]),  2 * (B[2][i][j-1] - B[2][i-1][j-1]) / Deltax[i], 2 * (B[2][i+1][j-1] - B[2][i][j-1]) / Deltax[i]);
    double slope_Bz_x_iminus1_jminus1 = minmod((B[2][i][j-1] - B[2][i-2][j-1]) / (2 * Deltax[i]),  2 * (B[2][i-1][j-1] - B[2][i-2][j-1]) / Deltax[i], 2 * (B[2][i][j-1] - B[2][i-1][j-1]) / Deltax[i]);
    
    // dEx/dy calculation
    
    //Ex{i-1/2,j-1/2}
    //G_1{i,j-1/2} calculation G_1

    double G_1_up = E[0][i][j] - 0.5 * Deltay * slope_Ex_y_ij;
    double G_1_down = E[0][i][j-1] + 0.5 * Deltay * slope_Ex_y_i_jminus1;

    double Bz_ij_up_y = B[2][i][j] - 0.5 * Deltay * slope_Bz_y_ij;
    double Bz_i_jminus1_down_y = B[2][i][j-1] + 0.5 * Deltay * slope_Bz_y_i_jminus1;

    double G_1 = 0.5 * (-G_1_up - G_1_down) - 0.5 * (Bz_ij_up_y - Bz_i_jminus1_down_y);

    //G_1{i-1,j-1/2} calculation G_2

    double G_2_up = E[0][i-1][j] - 0.5 * Deltay * slope_Ex_y_iminus1_j;
    double G_2_down = E[0][i-1][j-1] + 0.5 * Deltay * slope_Ex_y_iminus1_jminus1;

    double Bz_iminus1_j_up_y = B[2][i-1][j] - 0.5 * Deltay * slope_Bz_y_iminus1_j;
    double Bz_iminus1_jminus1_down_y = B[2][i-1][j-1] + 0.5 * Deltay * slope_Bz_y_iminus1_jminus1;

    double G_2 = 0.5 * (-G_2_up - G_2_down) - 0.5 * (Bz_iminus1_j_up_y - Bz_iminus1_jminus1_down_y);

    //F_1{i-1/2, j} calculation F_1

    double F_1_up = E[1][i][j] - 0.5 * Deltax[i] * slope_Ey_x_ij;
    double F_1_down = E[1][i-1][j] + 0.5 * Deltax[i] * slope_Ey_x_iminus1_j;
    
    double Bz_ij_up_x = B[2][i][j] - 0.5 * Deltax[i] * slope_Bz_x_ij;
    double Bz_iminus1_j_down_x = B[2][i-1][j] + 0.5 * Deltax[i] * slope_Bz_x_iminus1_j;

    double F_1 = 0.5 * (F_1_up + F_1_down) - 0.5 * (Bz_ij_up_x + Bz_iminus1_j_down_x);

    //F_1{i-1/2, j-1} calculation F_2

    double F_2_up = E[1][i][j-1] - 0.5 * Deltax[i] * slope_Ey_x_i_jminus1;
    double F_2_down = E[1][i-1][j-1] + 0.5 * Deltax[i] * slope_Ey_x_iminus1_jminus1;

    double Bz_i_jminus1_up_x = B[2][i][j-1] - 0.5 * Deltax[i] * slope_Bz_x_i_jminus1;
    double Bz_iminus1_jminus1_down_x = B[2][i-1][j-1] + 0.5 * Deltax[i] * slope_Bz_x_iminus1_jminus1;

    double F_2 = 0.5 * (F_2_up + F_2_down) - 0.5 * (Bz_i_jminus1_up_x - Bz_iminus1_jminus1_down_x);


    double Ex_ij = 0.25 * (G_1 + G_2 - F_1 - F_2);



}

/*
    Computes electromagnetic force along closed path around cell face. Used to compute time derivative of magnetic field.
    Inputs: E: electric field in reduced units along cell faces
            N_GC: number of ghost cells
            Deltax, Deltay: cell sizes in x and y directions in reduced units
    Output: Qx, Qy, Qz: line integral of E divided by cell face area on bottom cell faces (to update Bx) and left cell faces (to update By) respectively
*/

void Compute_RHS(ScalarField & Qx, ScalarField & Qy, ScalarField & Qz, ScalarField & Fx, ScalarField & Fy, ScalarField & Fz, VectorField & E, VectorField & H, VectorField & J, VectorField & D, VectorField & B, size_t N_GC, std::vector<double> & Deltax, double Deltay)
{

    for(size_t i=N_GC; i<E.shape()[1]-N_GC; i++){
        for(size_t j=N_GC; j<E.shape()[2]-N_GC; j++){
            auto [Ez_Qx_ij, Ez_i_jplus1] = Ez_Flux_Calculation_Qx(i, j, D, B, Deltax, Deltay);
            auto [Ez_Qy_ij, Ez_iplus1_j] = Ez_Flux_Calculation_Qy(i, j, E, B, Deltax, Deltay);
            //Qx[i][j] = -( E[2][i][j+1] - E[2][i][j] )/Deltay;
            Qx[i][j] = -(Ez_i_jplus1 - Ez_Qx_ij)/Deltay;
            
            if( i < E.shape()[1]-N_GC-1 ){
                //Qy[i][j] = -( E[2][i][j] - E[2][i+1][j] )/Deltax[i];
                //Qz[i][j] = ( E[0][i][j+1] - E[0][i][j] )/Deltay -( E[1][i+1][j] - E[1][i][j] )/Deltax[i];

                Qy[i][j] = -(Ez_Qy_ij - Ez_iplus1_j)/Deltax[i];
                //Qz[i][j] = ;

            }
        }
    }

      for(size_t i=N_GC; i<H.shape()[1]-N_GC; i++){
        for(size_t j=N_GC; j<H.shape()[2]-N_GC; j++){

            Fx[i][j] =  ( H[2][i][j+1] - H[2][i][j] )/Deltay - J[0][i][j]; //4pi/c?

            if (i < H.shape()[1] - N_GC - 1) {
                Fy[i][j] = ( H[2][i][j] - H[2][i+1][j] )/Deltax[i] - J[1][i][j];
                Fz[i][j] = -( H[0][i][j+1] - H[0][i][j] )/Deltay + ( H[1][i+1][j] - H[1][i][j] )/Deltax[i] - J[2][i][j];

            }
        }
    }

    return;
}
