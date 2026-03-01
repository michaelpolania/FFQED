#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>

#include "common.h"
#include "initial_conditions.h"
#include "field_evolution.h"

struct Fields
{
    VectorField & E;
    VectorField & B;
    VectorField & H;
    VectorField & D;
    ScalarField & Rho;
    const & dm;

};

double compute_A1_x(int i, int j, const Fields & f){

    double rho_avg = 0.5 * (f.Rho[i][j] + f.Rho[i][j-1]);
    
    double Ey = 0.25 * (f.E[1][i][j] + f.E[1][i+1][j] + f.E[1][i][j-1] + f.E[1][i+1][j-1]);   

    double Ez = 0.5 * (f.E[2][i][j] + f.E[2][i+1][j]);

    double Bx_i = 0.5 * (f.B[0][i+1][j] + f.B[0][i][j]);
    double Bx_i_minus_1 = 0.5 * (f.B[0][i+1][j-1] + f.B[0][i][j-1]);
    double By_j = 0.5 * (f.B[1][i][j+1] + f.B[1][i][j]);
    double By_j_minus_1 = 0.5 * (f.B[1][i][j-1] + f.B[1][i][j]); 
    double Bz = f.B[2][i][j];
    double Bz_j_minus_1 = f.B[2][i][j-1];

    double B_ij = (Bx_i * Bx_i)/4.0  + (By_j * By_j)/4.0  + (Bz * Bz); 
    double B_ij_minus_1 = (Bx_i_minus_1 * Bx_i_minus_1)/4.0 +(By_j_minus_1 * By_j_minus_1)/4.0 + (Bz_j_minus_1 * Bz_j_minus_1);
            
    double Bz_A_1_x = 0.5 * (Bz + Bz_j_minus_1);

    double ExB_x = Ey * Bz_A_1_x  - Ez * f.B[1][i][j];

    return (2.0 * rho_avg * ExB_x)/(B_ij_minus_1 + B_ij);
}

double compute_A1_y(int i, int j, const Fields & f){
    // Calculation of A_1_y

    double rho_avg = 0.5 * (f.Rho[i][j] + f.Rho[i-1][j]);
    
    double Ex = 0.25 * (f.E[0][i][j] + f.E[0][i][j+1] + f.E[0][i-1][j] + f.E[0][i-1][j+1]);   

    double Bz = 0.5 * (f.B[2][i][j] + f.B[2][i-1][j]);

    double Ez = 0.5 * (f.E[2][i][j] + f.E[2][i][j+1]);

    double Bx_i = 0.5 * (f.B[0][i+1][j] + f.B[0][i][j]);
    double Bx_i_minus_1_j = 0.5 * (f.B[0][i-1][j] + f.B[0][i][j]);
    double By_j = 0.5 * (f.B[1][i][j+1] + f.B[1][i][j]);
    double By_i_minus_1_j = 0.5 * (f.B[1][i-1][j+1] + f.B[1][i-1][j]); 
    double Bz = f.B[2][i][j];
    double Bz_i_minus_1_j = f.B[2][i-1][j];

    double B_ij = (Bx_i * Bx_i)/4.0  + (By_j * By_j)/4.0  + (Bz * Bz); 
    double B_ij_minus_1 = (Bx_i_minus_1_j * Bx_i_minus_1_j)/4.0 +(By_i_minus_1_j * By_i_minus_1_j)/4.0 + (Bz_i_minus_1_j * Bz_i_minus_1_j);
            
    double Bz_A_1_x = 0.5 * (Bz + Bz_j_minus_1);

    double ExB_y = Ez_A_1_y * f.B[0][i][j]  - Ex_A_1_y * Bz_A_1_y;

    return (2.0 * rho_avg * ExB_y)/(B_ij_minus_1 + B_ij);
}

double compute_A1_z(int i, int j, const Fields & f){
    // Calculation of (A_1)_z

    double rho_avg = 0.25 * (f.Rho[i][j] + f.Rho[i][j-1] + f.Rho[i-1][j] + f.Rho[i-1][j-1]);

    double Ex_By_A = 0.5 * (f.E[0][i][j] * f.B[1][i][j] + f.E[0][i-1][j] * f.B[1][i-1][j]); // at i-1/2, j-1/2
    double Ey_Bx_A = 0.5 * (f.E[1][i][j] * f.B[0][i][j] + f.E[1][i][j-1] * f.B[0][i][j-1]); // at i-1/2, j-1/2

    double Bx_i = 0.5 * (f.B[0][i+1][j] + f.B[0][i][j]);
    double Bx_i_minus_1_j = 0.5 * (f.B[0][i][j] + f.B[0][i-1][j]);
    double Bx_i_j_minus_1 = 0.5 * (f.B[0][i+1][j-1] + f.B[0][i][j-1]);
    double Bx_ij_minus_1 = 0.5 * (f.B[0][i][j-1] + f.B[0][i-1][j-1]); //fix

    double By_j = 0.5 * (f.B[1][i][j+1] + f.B[1][i][j]);
    double By_j_minus_1 = 0.5 * (f.B[1][i][j-1] + f.B[1][i][j]); 
    double By_i_minus_1_j = 0.5 * (f.B[1][i-1][j+1] + f.B[1][i-1][j]);
    double By_ij_minus_1 = 0.5 * (f.B[1][i-1][j] + f.B[1][i-1][j-1]);  
    

    double Bz = f.B[2][i][j];
    double Bz_i_minus_1_j = f.B[2][i-1][j];
    double Bz_ij_minus_1 = f.B[2][i-1][j-1];
    double Bz_i_j_minus_1 = f.B[2][i][j-1];

    double B_ij = (Bx_i * Bx_i)/4.0  + (By_j * By_j)/4.0  + (Bz * Bz); 
    double B_i_minus_1_j = (Bx_i_minus_1_j * Bx_i_minus_1_j)/4.0 +(By_i_minus_1_j * By_i_minus_1_j)/4.0 + (Bz_i_minus_1_j * Bz_i_minus_1_j);
    double B_ij_minus_1 = (Bx_ij_minus_1 * Bx_ij_minus_1)/4.0 + (By_ij_minus_1 * By_ij_minus_1)/4.0 + (Bz_ij_minus_1 * Bz_ij_minus_1);
    double B_i_j_minus_1 = (Bx_i_j_minus_1 * Bx_i_j_minus_1)/4.0 + (By_j_minus_1 * By_j_minus_1)/4.0 + (Bz_i_j_minus_1 * Bz_i_j_minus_1); 

    double B_squared = 2 * (1)/((B_ij + B_i_minus_1_j + B_ij_minus_1 + B_i_j_minus_1));

    return (4 * rho_avg * (Ex_By_A - Ey_Bx_A))/(B_squared);

}

double compute_A2_x(int i, int j, const Fields & f){
        
        // Calculation of (A_2)_x

        // Calculates B dot (curl(H)) at i,j

        double Bx_Hz_A = f.B[0][i][j] * (f.H[2][i][j+1] - f.H[2][i][j])/(f.dm.Deltay); // at i - 1/2, j, f.dm. might be a problem
        double Bx_Hz_B = f.B[0][i+1][j] * (f.H[2][i+1][j+1] - f.H[2][i+1][j])/(f.dm.Deltay); // at i + 1/2, j

        double Bx_Hz_AB = 0.5 * (Bx_Hz_A + Bx_Hz_B); // at i , j

        double By_Hz_A = f.B[1][i][j+1] * (f.H[2][i+1][j+1] - f.H[2][i][j+1])/(f.dm.Deltax[i]); // at i, j + 1/2
        double By_Hz_B = f.B[1][i][j] * (f.H[2][i+1][j] - f.H[2][i][j])/(f.dm.Deltax[i]); // at i, j - 1/2

        double By_Hz_AB = 0.5 * (By_Hz_A + By_Hz_B); // at i,j

        double Bz_Hy_Hx_AB = f.B[2][i][j] * ((f.H[1][i+1][j] - f.H[1][i][j])/(f.dm.Deltax[i]) + (f.H[0][i][j] - f.H[0][i][j+1])/(f.dm.Deltay)); // at i,j

        double B_dot_curl_H_AB = Bx_Hz_AB -  By_Hz_AB + Bz_Hy_Hx_AB; // at i,j

        // Calculates B dot (curl(H)) at i , j-1

        double Bx_Hz_C = f.B[0][i+1][j-1] * (f.H[2][i+1][j] - f.H[2][i+1][j-1])/(f.dm.Deltay); // at i + 1/2, j - 1
        double Bx_Hz_D = f.B[0][i][j-1] * (f.H[2][i][j] - f.H[2][i][j-1])/(f.dm.Deltay); // at i - 1/2, j-1

        double Bx_Hz_CD = 0.5 * (Bx_Hz_C + Bx_Hz_D); // at i , j - 1

        double By_Hz_C = f.B[1][i][j] * (f.H[2][i+1][j] - f.H[2][i][j])/(f.dm.Deltax[i]); // at i, j - 1/2
        double By_Hz_D = f.B[1][i][j-1] * (f.H[2][i+1][j-1] - f.H[2][i][j-1])/(f.dm.Deltax[i]); // at i , j - 3/2

        double By_Hz_CD = 0.5 * (By_Hz_C + By_Hz_D); // at i, j - 1

        double Bz_Hy_Hx_CD = f.B[2][i][j-1] * ((f.H[1][i+1][j] - f.H[1][i][j])/(f.dm.Deltax[i]) + (f.H[0][i][j-1] - f.H[0][i][j])/(f.dm.Deltay)); // at i,j-1 //continue working

        double B_dot_curl_H_CD = Bx_Hz_CD -  By_Hz_CD + Bz_Hy_Hx_CD; // at i,j-1

        double Bx_i = 0.5 * (f.B[0][i+1][j] + f.B[0][i][j]);
        double Bx_i_minus_1 = 0.5 * (f.B[0][i+1][j-1] + f.B[0][i][j-1]);
        double By_j = 0.5 * (f.B[1][i][j+1] + f.B[1][i][j]);
        double By_j_minus_1 = 0.5 * (f.B[1][i][j-1] + f.B[1][i][j]); 
        double Bz = f.B[2][i][j];
        double Bz_j_minus_1 = f.B[2][i][j-1];

        double B_ij = (Bx_i * Bx_i)/4.0  + (By_j * By_j)/4.0  + (Bz * Bz); 
        double B_ij_minus_1 = (Bx_i_minus_1 * Bx_i_minus_1)/4.0 +(By_j_minus_1 * By_j_minus_1)/4.0 + (Bz_j_minus_1 * Bz_j_minus_1);
        
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

    double B_dot_curl_H_AB = Bx_Hz_AB -  By_Hz_AB + Bz_Hy_Hx_AB; // at i,j

    // Calculates B dot (curl(H)) at i-1 , j

    double Bx_Hz_C = f.B[0][i][j] * (f.H[2][i][j+1] - f.H[2][i][j])/(f.dm.Deltay); // at i - 1/2, j 
    double Bx_Hz_D = f.B[0][i-1][j] * (f.H[2][i-1][j+1] - f.H[2][i-1][j-1])/(f.dm.Deltay); // at i - 3/2, j

    double Bx_Hz_CD = 0.5 * (Bx_Hz_C + Bx_Hz_D); // at i-1 , j

    double By_Hz_C = f.B[1][i-1][j+1] * (f.H[2][i][j+1] - f.H[2][i-1][j+1])/(f.dm.Deltax[i]); // at i, j - 1/2
    double By_Hz_D = f.B[1][i-1][j] * (f.H[2][i][j] - f.H[2][i-1][j])/(f.dm.Deltax[i]); // at i , j - 3/2

    double By_Hz_CD = 0.5 * (By_Hz_C + By_Hz_D); // at i-1, j

    double Bz_Hy_Hx_CD = f.B[2][i-1][j] * ((f.H[1][i][j] - f.H[1][i-1][j])/(f.dm.Deltax[i]) + (H[0][i-1][j] - H[0][i-1][j+1])/(dm.Deltay)); // at i,j-1 //continue working

    double B_dot_curl_H_CD = Bx_Hz_CD -  By_Hz_CD + Bz_Hy_Hx_CD; // at i,j-1

    double Bx_i = 0.5 * (f.B[0][i+1][j] + f.B[0][i][j]);
    double Bx_i_minus_1_j = 0.5 * (f.B[0][i][j] + f.B[0][i-1][j]);
    double By_j = 0.5 * (f.B[1][i][j+1] + f.B[1][i][j]);
    double By_i_minus_1_j = 0.5 * (f.B[1][i-1][j+1] + f.B[1][i-1][j]); 
    double Bz = f.B[2][i][j];
    double Bz_i_minus_1_j = f.B[2][i-1][j];

    double B_ij = (Bx_i * Bx_i)/4.0  + (By_j * By_j)/4.0  + (Bz * Bz); 
    double B_i_minus_1_j = (Bx_i_minus_1_j * Bx_i_minus_1_j)/4.0 +(By_i_minus_1_j * By_i_minus_1_j)/4.0 + (Bz_i_minus_1_j * Bz_i_minus_1_j);
        
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

    double B_dot_curl_H_AB = Bx_Hz_AB -  By_Hz_AB + Bz_Hy_Hx_AB; // at i,j

    // Calculates B dot (curl(H)) at i-1 , j

    double Bx_Hz_C = f.B[0][i][j] * (f.H[2][i][j+1] - f.H[2][i][j])/(f.dm.Deltay); // at i - 1/2, j 
    double Bx_Hz_D = f.B[0][i-1][j] * (f.H[2][i-1][j+1] - f.H[2][i-1][j-1])/(f.dm.Deltay); // at i - 3/2, j

    double Bx_Hz_CD = 0.5 * (Bx_Hz_C + Bx_Hz_D); // at i-1 , j

    double By_Hz_C = f.B[1][i-1][j+1] * (f.H[2][i][j+1] - f.H[2][i-1][j+1])/(f.dm.Deltax[i]); // at i, j - 1/2
    double By_Hz_D = f.B[1][i-1][j] * (f.H[2][i][j] - f.H[2][i-1][j])/(f.dm.Deltax[i]); // at i , j - 3/2

    double By_Hz_CD = 0.5 * (By_Hz_C + By_Hz_D); // at i-1, j

    double Bz_Hy_Hx_CD = f.B[2][i-1][j] * ((f.H[1][i][j] - f.H[1][i-1][j])/(f.dm.Deltax[i]) + (H[0][i-1][j] - H[0][i-1][j+1])/(dm.Deltay)); // at i-1,j

    double B_dot_curl_H_CD = Bx_Hz_CD -  By_Hz_CD + Bz_Hy_Hx_CD; // at i-1,j

    // Calculates B dot (curl(H)) at i-1, j-1

    double Bx_Hz_E = f.B[0][i][j-1] * (f.H[2][i][j] - f.H[2][i][j-1])/(f.dm.Deltay); // at i - 1/2, j-1 
    double Bx_Hz_F = f.B[0][i-1][j-1] * (f.H[2][i-1][j] - f.H[2][i-1][j-1])/(f.dm.Deltay); // at i - 3/2, j-1

    double Bx_Hz_EF = 0.5 * (Bx_Hz_E + Bx_Hz_F); // at i-1 , j-1

    double By_Hz_E = f.B[1][i-1][j] * (f.H[2][i][j] - f.H[2][i-1][j])/(f.dm.Deltax[i]); // at i, j - 1/2
    double By_Hz_F = f.B[1][i-1][j-1] * (f.H[2][i][j-1] - f.H[2][i-1][j-1])/(f.dm.Deltax[i]); // at i , j - 3/2

    double By_Hz_EF = 0.5 * (By_Hz_E + By_Hz_F); // at i-1, j-1

    double Bz_Hy_Hx_EF = f.B[2][i-1][j-1] * ((f.H[1][i][j-1] - f.H[1][i-1][j-1])/(f.dm.Deltax[i]) + (H[0][i-1][j-1] - H[0][i-1][j])/(dm.Deltay)); // at i-1,j

    double B_dot_curl_H_EF = Bx_Hz_EF -  By_Hz_EF + Bz_Hy_Hx_EF; // at i-1,j-1

    // Calculates B dot (curl(H)) at i, j-1

    double Bx_Hz_G = f.B[0][i+1][j-1] * (f.H[2][i+1][j] - f.H[2][i+1][j-1])/(f.dm.Deltay); // at i + 1/2, j - 1
    double Bx_Hz_H = f.B[0][i][j-1] * (f.H[2][i][j] - f.H[2][i][j-1])/(f.dm.Deltay); // at i - 1/2, j-1

    double Bx_Hz_GH = 0.5 * (Bx_Hz_C + Bx_Hz_D); // at i , j - 1

    double By_Hz_G = f.B[1][i][j] * (f.H[2][i+1][j] - f.H[2][i][j])/(f.dm.Deltax[i]); // at i, j - 1/2
    double By_Hz_H = f.B[1][i][j-1] * (f.H[2][i+1][j-1] - f.H[2][i][j-1])/(f.dm.Deltax[i]); // at i , j - 3/2

    double By_Hz_GH = 0.5 * (By_Hz_C + By_Hz_D); // at i, j - 1

    double Bz_Hy_Hx_GH = f.B[2][i][j-1] * ((f.H[1][i+1][j] - f.H[1][i][j])/(f.dm.Deltax[i]) + (f.H[0][i][j-1] - f.H[0][i][j])/(f.dm.Deltay)); // at i,j-1 //continue working
    double B_dot_curl_H_GH = Bx_Hz_GH - By_Hz_GH + Bz_Hy_Hx_GH; 

    double Bx_i = 0.5 * (f.B[0][i+1][j] + f.B[0][i][j]);
    double Bx_i_minus_1_j = 0.5 * (f.B[0][i][j] + f.B[0][i-1][j]);
    double Bx_i_j_minus_1 = 0.5 * (f.B[0][i+1][j-1] + f.B[0][i][j-1]);
    double Bx_ij_minus_1 = 0.5 * (f.B[0][i][j-1] + f.B[0][i-1][j-1]); //fix

    double By_j = 0.5 * (f.B[1][i][j+1] + f.B[1][i][j]);
    double By_j_minus_1 = 0.5 * (f.B[1][i][j-1] + f.B[1][i][j]); 
    double By_i_minus_1_j = 0.5 * (f.B[1][i-1][j+1] + f.B[1][i-1][j]);
    double By_ij_minus_1 = 0.5 * (f.B[1][i-1][j] + f.B[1][i-1][j-1]);  
    

    double Bz = f.B[2][i][j];
    double Bz_i_minus_1_j = f.B[2][i-1][j];
    double Bz_ij_minus_1 = f.B[2][i-1][j-1];
    double Bz_i_j_minus_1 = f.B[2][i][j-1];

    double B_ij = (Bx_i * Bx_i)/4.0  + (By_j * By_j)/4.0  + (Bz * Bz); 
    double B_i_minus_1_j = (Bx_i_minus_1_j * Bx_i_minus_1_j)/4.0 +(By_i_minus_1_j * By_i_minus_1_j)/4.0 + (Bz_i_minus_1_j * Bz_i_minus_1_j);
    double B_ij_minus_1 = (Bx_ij_minus_1 * Bx_ij_minus_1)/4.0 + (By_ij_minus_1 * By_ij_minus_1)/4.0 + (Bz_ij_minus_1 * Bz_ij_minus_1);
    double B_i_j_minus_1 = (Bx_i_j_minus_1 * Bx_i_j_minus_1)/4.0 + (By_j_minus_1 * By_j_minus_1)/4.0 + (Bz_i_j_minus_1 * Bz_i_j_minus_1); 
    
    return 0.25 * ((B_dot_curl_H_AB * Bz)/(B_ij) + (B_dot_curl_H_CD * Bz_i_minus_1_j)/(B_i_minus_1_j) + (B_dot_curl_H_EF * Bz_i_j_minus_1)/(B_i_j_minus_1) + (B_dot_curl_H_GH * Bz_ij_minus_1)/(B_ij_minus_1));
}

double compute_A3_x(int i, int j, const Fields & f){

    //Calculation of A_3_x
            
    // Calculates -D dot (curl(E)) at i,j

    double Dx_Ez_A = f.D[0][i][j] * (f.E[2][i][j] - f.E[2][i][j+1])/(f.dm.Deltay); // at i - 1/2, j
    double Dx_Ez_B = f.D.[0][i+1][j] * (f.E[2][i+1][j] - f.E[2][i+1][j+1])/(f.dm.Deltay); // at i + 1/2, j

    double Dx_Ez_AB = 0.5 * (Dx_Ez_A + Dx_Ez_B); // at i , j

    double Dy_Ez_A = f.D[1][i][j+1] * (f.E[2][i+1][j+1] - f.E[2][i][j+1])/(f.dm.Deltax[i]); // at i, j + 1/2
    double Dy_Ez_B = f.D[1][i][j] * (f.E[2][i+1][j] - f.E[2][i][j])/(f.dm.Deltax[i]); // at i, j - 1/2

    double Dy_Ez_AB = 0.5 * (Dy_Hz_A + Dy_Hz_B); // at i,j

    double Dz_Ey_Ex_AB = f.D[2][i][j] * ((f.E[0][i][j+1] - f.E[0][i][j])/(f.dm.Deltay) + (E[1][i][j] - E[1][i+1][j])/(dm.Deltax[i])); // at i,j

    double D_dot_curl_E_AB = Dx_Ez_AB +  Dy_Ez_AB + Dz_Ey_Ex_AB; // at i,j

    // Calculates -D dot (curl(E)) at i , j-1

    double Dx_Ez_C = f.D[0][i+1][j-1] * (f.E[2][i+1][j-1] - f.E[2][i+1][j])/(f.dm.Deltay); // at i + 1/2, j - 1
    double Dx_Ez_D = f.D[0][i][j-1] * (f.E[2][i][j-1] - f.E[2][i][j])/(f.dm.Deltay); // at i - 1/2, j-1

    double Dx_Ez_CD = 0.5 * (Dx_Ez_C + Dx_Ez_D); // at i , j - 1

    double Dy_Ez_C = f.D[1][i][j] * (f.E[2][i+1][j] - f.E[2][i][j])/(f.dm.Deltax[i]); // at i, j - 1/2
    double Dy_Ez_D = f.D[1][i][j-1] * (f.E[2][i+1][j-1] - f.E[2][i][j-1])/(f.dm.Deltax[i]); // at i , j - 3/2

    double Dy_Ez_CD = 0.5 * (Dy_Ez_C + Dy_Ez_D); // at i, j - 1

    double Dz_Ey_Ex_CD = f.D[2][i][j-1] * ((f.E[0][i][j] - f.E[0][i][j-1])/(f.dm.Deltay) + (E[1][i][j-1] - E[1][i+1][j-1])/(dm.Deltax[i])); // at i,j-1

    double D_dot_curl_E_CD = Dx_Ez_CD +  Dy_Ez_CD + Dz_Ey_Ex_CD; // at i,j-1

    double Bx_i = 0.5 * (f.B[0][i+1][j] + f.B[0][i][j]);
    double Bx_i_minus_1 = 0.5 * (f.B[0][i+1][j-1] + f.B[0][i][j-1]);
    double By_j = 0.5 * (f.B[1][i][j+1] + f.B[1][i][j]);
    double By_j_minus_1 = 0.5 * (f.B[1][i][j-1] + f.B[1][i][j]); 
    double Bz = f.B[2][i][j];
    double Bz_j_minus_1 = f.B[2][i][j-1];

    double B_ij = (Bx_i * Bx_i)/4.0  + (By_j * By_j)/4.0  + (Bz * Bz); 
    double B_ij_minus_1 = (Bx_i_minus_1 * Bx_i_minus_1)/4.0 +(By_j_minus_1 * By_j_minus_1)/4.0 + (Bz_j_minus_1 * Bz_j_minus_1);
    
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

    double Dy_Ez_AB = 0.5 * (Dy_Hz_A + Dy_Hz_B); // at i,j

    double Dz_Ey_Ex_AB = f.D[2][i][j] * ((f.E[0][i][j+1] - f.E[0][i][j])/(f.dm.Deltay) + (E[1][i][j] - E[1][i+1][j])/(dm.Deltax[i])); // at i,j

    double D_dot_curl_E_AB = Dx_Ez_AB +  Dy_Ez_AB + Dz_Ey_Ex_AB; // at i,j

    // Calculates -D dot (curl(E)) at i-1 , j

    double Dx_Ez_C = f.D[0][i][j] * (f.E[2][i][j] - f.E[2][i][j+1])/(f.dm.Deltay); // at i - 1/2, j 
    double Dx_Ez_D = f.D[0][i-1][j] * (f.E[2][i-1][j-1] - f.E[2][i-1][j+1])/(f.dm.Deltay); // at i - 3/2, j

    double Dx_Ez_CD = 0.5 * (Dx_Ez_G + Dx_Ez_H); // at i-1 , j

    double Dy_Ez_C = f.D[1][i-1][j+1] * (f.E[2][i][j+1] - f.E[2][i-1][j+1])/(f.dm.Deltax[i]); // at i, j - 1/2
    double Dy_Ez_D = f.D[1][i-1][j] * (f.E[2][i][j] - f.E[2][i-1][j])/(f.dm.Deltax[i]); // at i , j - 3/2

    double Dy_Ez_CD = 0.5 * (Dy_Ez_G + Dy_Ez_H); // at i-1, j

    double Dz_Ey_Ex_CD = f.D[2][i-1][j] * ((f.E[1][i-1][j] - f.E[1][i][j])/(f.dm.Deltax[i]) + (E[0][i-1][j] - E[0][i-1][j+1])/(dm.Deltay)); // at i-1,j

    double D_dot_curl_E_CD = Dx_Ez_CD + Dy_Ez_CD + Dz_Ey_Ex_CD; // at i-1,j

    double Bx_i = 0.5 * (f.B[0][i+1][j] + f.B[0][i][j]);
    double Bx_i_minus_1_j = 0.5 * (f.B[0][i][j] + f.B[0][i-1][j]);

    double By_j = 0.5 * (f.B[1][i][j+1] + f.B[1][i][j]);
    double By_i_minus_1_j = 0.5 * (f.B[1][i-1][j+1] + f.B[1][i-1][j]);
 
    double Bz = f.B[2][i][j];
    double Bz_i_minus_1_j = f.B[2][i-1][j];

    double B_ij = (Bx_i * Bx_i)/4.0  + (By_j * By_j)/4.0  + (Bz * Bz); 
    double B_i_minus_1_j = (Bx_i_minus_1_j * Bx_i_minus_1_j)/4.0 +(By_i_minus_1_j * By_i_minus_1_j)/4.0 + (Bz_i_minus_1_j * Bz_i_minus_1_j);
    
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

    double Dy_Ez_AB = 0.5 * (Dy_Hz_A + Dy_Hz_B); // at i,j

    double Dz_Ey_Ex_AB = f.D[2][i][j] * ((f.E[0][i][j+1] - f.E[0][i][j])/(f.dm.Deltay) + (E[1][i][j] - E[1][i+1][j])/(dm.Deltax[i])); // at i,j

    double D_dot_curl_E_AB = Dx_Ez_AB +  Dy_Ez_AB + Dz_Ey_Ex_AB; // at i,j

    // Calculates -D dot (curl(E)) at i , j-1

    double Dx_Ez_C = f.D[0][i+1][j-1] * (f.E[2][i+1][j-1] - f.E[2][i+1][j])/(f.dm.Deltay); // at i + 1/2, j - 1
    double Dx_Ez_D = f.D[0][i][j-1] * (f.E[2][i][j-1] - f.E[2][i][j])/(f.dm.Deltay); // at i - 1/2, j-1

    double Dx_Ez_CD = 0.5 * (Dx_Ez_C + Dx_Ez_D); // at i , j - 1

    double Dy_Ez_C = f.D[1][i][j] * (f.E[2][i+1][j] - f.E[2][i][j])/(f.dm.Deltax[i]); // at i, j - 1/2
    double Dy_Ez_D = f.D[1][i][j-1] * (f.E[2][i+1][j-1] - f.E[2][i][j-1])/(f.dm.Deltax[i]); // at i , j - 3/2

    double Dy_Ez_CD = 0.5 * (Dy_Ez_C + Dy_Ez_D); // at i, j - 1

    double Dz_Ey_Ex_CD = f.D[2][i][j-1] * ((f.E[0][i][j] - f.E[0][i][j-1])/(f.dm.Deltay) + (E[1][i][j-1] - E[1][i+1][j-1])/(dm.Deltax[i])); // at i,j-1

    double D_dot_curl_E_CD = Dx_Ez_CD +  Dy_Ez_CD + Dz_Ey_Ex_CD; // at i,j-1

    // Calculates -D dot (curl(E)) at i-1, j-1

    double Dx_Ez_E = f.D[0][i][j-1] * (f.E[2][i][j-1] - f.E[2][i][j])/(f.dm.Deltay); // at i - 1/2, j-1 
    double Dx_Ez_F = f.D[0][i-1][j-1] * (f.E[2][i-1][j-1] - f.E[2][i-1][j])/(f.dm.Deltay); // at i - 3/2, j-1

    double Dx_Ez_EF = 0.5 * (Dx_Hz_E + Dx_Hz_F); // at i-1 , j-1

    double Dy_Ez_E = f.D[1][i-1][j] * (f.E[2][i][j] - f.E[2][i-1][j])/(f.dm.Deltax[i]); // at i, j - 1/2
    double Dy_Ez_F = f.D[1][i-1][j-1] * (f.E[2][i][j-1] - f.E[2][i-1][j-1])/(f.dm.Deltax[i]); // at i , j - 3/2

    double Dy_Ez_EF = 0.5 * (Dy_Ez_E + Dy_Ez_F); // at i-1, j-1

    double Dz_Ey_Ex_EF = f.D[2][i-1][j-1] * ((f.E[1][i][j-1] - f.E[1][i-1][j-1])/(f.dm.Deltax[i]) + (E[0][i-1][j-1] - E[0][i-1][j])/(dm.Deltay)); // at i-1,j

    double D_dot_curl_E_EF = Dx_Ez_EF + Dy_Ez_EF + Dz_Ey_Ex_EF; // at i-1,j-1

    // Calculates -D dot (curl(E)) at i-1, j

    double Dx_Ez_G = f.D[0][i][j] * (f.E[2][i][j] - f.E[2][i][j+1])/(f.dm.Deltay); // at i - 1/2, j 
    double Dx_Ez_H = f.D[0][i-1][j] * (f.E[2][i-1][j-1] - f.E[2][i-1][j+1])/(f.dm.Deltay); // at i - 3/2, j

    double Dx_Ez_GH = 0.5 * (Dx_Ez_G + Dx_Ez_H); // at i-1 , j

    double Dy_Ez_G = f.D[1][i-1][j+1] * (f.E[2][i][j+1] - f.E[2][i-1][j+1])/(f.dm.Deltax[i]); // at i, j - 1/2
    double Dy_Ez_H = f.D[1][i-1][j] * (f.E[2][i][j] - f.E[2][i-1][j])/(f.dm.Deltax[i]); // at i , j - 3/2

    double Dy_Ez_GH = 0.5 * (Dy_Ez_G + Dy_Ez_H); // at i-1, j

    double Dz_Ey_Ex_GH = f.D[2][i-1][j] * ((f.E[1][i-1][j] - f.E[1][i][j])/(f.dm.Deltax[i]) + (E[0][i-1][j] - E[0][i-1][j+1])/(dm.Deltay)); // at i-1,j

    double D_dot_curl_E_GH = Dx_Ez_GH + Dy_Ez_GH + Dz_Ey_Ex_GH; // at i-1,j


    double Bx_i = 0.5 * (f.B[0][i+1][j] + f.B[0][i][j]);
    double Bx_i_minus_1_j = 0.5 * (f.B[0][i][j] + f.B[0][i-1][j]);
    double Bx_i_j_minus_1 = 0.5 * (f.B[0][i+1][j-1] + f.B[0][i][j-1]);
    double Bx_ij_minus_1 = 0.5 * (f.B[0][i][j-1] + f.B[0][i-1][j-1]); //fix

    double By_j = 0.5 * (f.B[1][i][j+1] + f.B[1][i][j]);
    double By_j_minus_1 = 0.5 * (f.B[1][i][j-1] + f.B[1][i][j]); 
    double By_i_minus_1_j = 0.5 * (f.B[1][i-1][j+1] + f.B[1][i-1][j]);
    double By_ij_minus_1 = 0.5 * (f.B[1][i-1][j] + f.B[1][i-1][j-1]);  
    

    double Bz = f.B[2][i][j];
    double Bz_i_minus_1_j = f.B[2][i-1][j];
    double Bz_ij_minus_1 = f.B[2][i-1][j-1];
    double Bz_i_j_minus_1 = f.B[2][i][j-1];

    double B_ij = (Bx_i * Bx_i)/4.0  + (By_j * By_j)/4.0  + (Bz * Bz); 
    double B_i_minus_1_j = (Bx_i_minus_1_j * Bx_i_minus_1_j)/4.0 +(By_i_minus_1_j * By_i_minus_1_j)/4.0 + (Bz_i_minus_1_j * Bz_i_minus_1_j);
    double B_ij_minus_1 = (Bx_ij_minus_1 * Bx_ij_minus_1)/4.0 + (By_ij_minus_1 * By_ij_minus_1)/4.0 + (Bz_ij_minus_1 * Bz_ij_minus_1);
    double B_i_j_minus_1 = (Bx_i_j_minus_1 * Bx_i_j_minus_1)/4.0 + (By_j_minus_1 * By_j_minus_1)/4.0 + (Bz_i_j_minus_1 * Bz_i_j_minus_1); 
    
    return 0.25 * ((D_dot_curl_E_AB * Bz)/(B_ij) + (D_dot_curl_E_CD * Bz_i_j_minus_1)/(B_i_j_minus_1) + (D_dot_curl_E_EF * Bz_ij_minus_1)/(B_ij_minus_1) + (D_dot_curl_E_GH * Bz_i_minus_1_j)/(B_i_minus_1_j));
}


/*
    Computes current divided by 4*pi/c components along cell edges in reduced units i.e., c*t_0/(L_0*B_0)*E

    Input: B: magnetic field in reduced units at cell centers
           N_GC: number of ghost cells
           dm: Domain object containing information about the simulation domain
    Output: E, J: electric field and current density in reduced units
*/
void Compute_J(VectorField & B, VectorField & J, size_t N_GC, const Domain & dm)
{
    double Deltax_n1;
    for(size_t i=N_GC; i<B.shape()[1]-N_GC; i++){
        Deltax_n1 = 0.5*(dm.Deltax[i-1]+dm.Deltax[i]);
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

            J[0][i][j] = compute_A1_x(i, j, Fields) + compute_A2_x(i, j, Fields) + compute_A3_x(i, j, Fields);
            J[1][i][j] = compute_A1_y(i, j, Fields) + compute_A2_y(i, j, Fields) + compute_A3_y(i, j, Fields);
            J[2][i][j] = compute_A1_z(i, j, Fields) + compute_A2_z(i, j, Fields) + compute_A3_z(i, j, Fields);
        }
    }

    return;
}

void Compute_J(VectorField & B, VectorField & J, size_t N_GC, const Domain & dm)
{
    double Deltax_n1;
    for(size_t i=N_GC; i<B.shape()[1]-N_GC; i++){
        Deltax_n1 = 0.5*(dm.Deltax[i-1]+dm.Deltax[i]);
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

            auto J_vec = compute_current_density(B, J, H, D, Rho, i, j)
            


            J[0][i][j] = J_vec.Jx;
            J[1][i][j] = J_vec.Jy;                                                                                                                                                                                                                                                                                                               
            J[2][i][j] = J_vec.Jz;
        }
    }

    return;
}

double rn_uniform() {
    // Making rng static ensures that it stays the same
    // Between different invocations of the function
    static std::default_random_engine rng;

    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
}

double rn_weibull() {
    // Making rng static ensures that it stays the same
    // Between different invocations of the function
    static std::default_random_engine rng;

    std::weibull_distribution<double> dist(0.5,0.5);

    return dist(rng);
}

/*
    Computes c*electric field components along cell edges in reduced units i.e., c*t_0/(L_0*B_0)*E

    Inputs: B, J, vc: magnetic field in reduced units at cell face centers, current density and velocity field at cell face edges
            Bn: updated intermediate values of B, in which only toroidal field is included
            N_GC: number of ghost cells
            tC: transCoeffs object containing Ohmic and Hall diffusivities in reduced units and shear modulus and mass density in cgs units
            Deltax, Deltay: cell sizes in x and y directions in reduced units
            t: time in reduced units
            dm: Domain object containing information about the simulation domain
            bparams: BandBCParams object containing information about the initial magnetic field and boundary conditions
    Output: E: electric field on cell face edges in reduced units
*/
void Compute_E(VectorField & B, VectorField & Bn, VectorField & E, VectorField & J, VectorField & vc, TransCoeffs & tC, size_t N_GC, double t, const Domain & dm, const BandBCParams & bparams)
{

    std::vector<double> Deltax = dm.Deltax;
    double Deltay = dm.Deltay;

    std::vector<double> x = dm.x; //cell centers in the x-direction. Note: no ghost cells.
    std::vector<double> y = dm.y; //cell centers in the y-direction. Note: no ghost cells.

    double vx_av, vy_av, vz_av; //averaged values of the velocity = -current density
    double Bx_up, By_up, Bz_up; //upwinded field values

    double alphaz_jm1, alphaz_j, alphaz_im1, alphaz_i;
    double alphay_im1, alphay_i, alphax_jm1, alphax_j;
    double Deltax_n2, Deltax_n1, Deltax_1; //spacing between cell centres of i-2 and i-1, i-1 and i, i and i+1 cells, respectively.

    //Compute electric field at edges using magnetic field
    //E[0][i][j] = E_x^{i,j-1/2}
    //E[1][i][j] = E_y^{i-1/2,j}
    //E[2][i][j] = E_z^{i-1/2,j-1/2}
    for(size_t i=N_GC; i<B.shape()[1]-N_GC; i++){
        Deltax_n2 = 0.5*(Deltax[i-2]+Deltax[i-1]);
        Deltax_n1 = 0.5*(Deltax[i-1]+Deltax[i]);
        Deltax_1 = 0.5*(Deltax[i]+Deltax[i+1]);
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

            J[0][i][j] = ( Bn[2][i][j] - Bn[2][i][j-1] )/Deltay;
            J[1][i][j] = ( Bn[2][i-1][j] - Bn[2][i][j] )/Deltax_n1;
            if(i < B.shape()[1]-N_GC-1){
                vy_av = 0.25*( -0.5*(tC.eta_H[i][j]+tC.eta_H[i-1][j])*J[1][i][j] - 0.5*(tC.eta_H[i][j]+tC.eta_H[i+1][j])*( Bn[2][i][j] - Bn[2][i+1][j] )/Deltax_1
                                - 0.5*(tC.eta_H[i-1][j-1]+tC.eta_H[i][j-1])*( Bn[2][i-1][j-1] - Bn[2][i][j-1] )/Deltax_n1 - 0.5*(tC.eta_H[i][j-1]+tC.eta_H[i+1][j-1])*( Bn[2][i][j-1] - Bn[2][i+1][j-1] )/Deltax_1 );
                if(vy_av > 0){
                    alphaz_jm1 = minmod( (Bn[2][i][j]-Bn[2][i][j-2])/(2*Deltay), 2*(Bn[2][i][j]-Bn[2][i][j-1])/Deltay, 2*(Bn[2][i][j-1]-Bn[2][i][j-2])/Deltay );
                    Bz_up = Bn[2][i][j-1] + alphaz_jm1*(Deltay/2.);
                }
                else{
                    alphaz_j = minmod( (Bn[2][i][j+1]-Bn[2][i][j-1])/(2*Deltay), 2*(Bn[2][i][j+1]-Bn[2][i][j])/Deltay, 2*(Bn[2][i][j]-Bn[2][i][j-1])/Deltay );
                    Bz_up = Bn[2][i][j] - alphaz_j*(Deltay/2.);
                }
                vz_av = 0.5*( -0.25*(tC.eta_H[i][j]+tC.eta_H[i-1][j]+tC.eta_H[i][j-1]+tC.eta_H[i-1][j-1])*J[2][i][j] - 0.25*(tC.eta_H[i][j]+tC.eta_H[i][j-1]+tC.eta_H[i+1][j-1]+tC.eta_H[i+1][j])*J[2][i+1][j] );
                //Can use J[2][i+1][j] above and J[2][i][j+1] below because these reference already-up-to-date, computed currents
                By_up = B[1][i][j];
                E[0][i][j] = 0.5*(tC.eta_O[i][j]+tC.eta_O[i][j-1])*J[0][i][j] - ( vy_av*Bz_up - vz_av*By_up )
                            - ( 0.25*(vc[1][i][j]+vc[1][i+1][j]+vc[1][i][j-1]*vc[1][i+1][j-1])*0.5*(B[2][i][j]+B[2][i][j-1]) - 0.5*(vc[2][i][j]+vc[2][i+1][j])*B[1][i][j] );
            }

            vz_av = 0.5*( -0.25*(tC.eta_H[i][j]+tC.eta_H[i-1][j]+tC.eta_H[i][j-1]+tC.eta_H[i-1][j-1])*J[2][i][j] - 0.25*(tC.eta_H[i][j]+tC.eta_H[i-1][j]+tC.eta_H[i][j+1]+tC.eta_H[i-1][j+1])*J[2][i][j+1] );
            Bx_up = B[0][i][j];
            vx_av = 0.25*( -0.5*(tC.eta_H[i][j]+tC.eta_H[i][j-1])*J[0][i][j] - 0.5*(tC.eta_H[i-1][j]+tC.eta_H[i-1][j-1])*( Bn[2][i-1][j] - Bn[2][i-1][j-1] )/Deltay
                            - 0.5*(tC.eta_H[i][j+1]+tC.eta_H[i][j])*( Bn[2][i][j+1] - Bn[2][i][j] )/Deltay - 0.5*(tC.eta_H[i-1][j+1]+tC.eta_H[i-1][j])*( Bn[2][i-1][j+1] - Bn[2][i-1][j] )/Deltay );
            if(vx_av > 0){
                alphaz_im1 = minmod( (Bn[2][i][j]-Bn[2][i-2][j])/(Deltax_n2+Deltax_n1), 2*(Bn[2][i][j]-Bn[2][i-1][j])/Deltax_n1, 2*(Bn[2][i-1][j]-Bn[2][i-2][j])/Deltax_n2 );
                Bz_up = Bn[2][i-1][j] + alphaz_im1*(Deltax[i-1]/2.);
            }
            else{
                alphaz_i = minmod( (Bn[2][i+1][j]-Bn[2][i-1][j])/(Deltax_n1+Deltax_1), 2*(Bn[2][i+1][j]-Bn[2][i][j])/Deltax_1, 2*(Bn[2][i][j]-Bn[2][i-1][j])/Deltax_n1 );
                Bz_up = Bn[2][i][j] - alphaz_i*(Deltax[i]/2.);
            }
            E[1][i][j] = 0.5*(tC.eta_O[i][j]+tC.eta_O[i-1][j])*J[1][i][j] - ( vz_av*Bx_up - vx_av*Bz_up )
                        - ( 0.5*(vc[2][i][j]+vc[2][i][j+1])*B[0][i][j] - 0.25*(vc[0][i][j]+vc[0][i][j+1]+vc[0][i-1][j]+vc[0][i-1][j+1])*0.5*(B[2][i][j]+B[2][i-1][j]) );

            vx_av = 0.5*( -0.5*(tC.eta_H[i][j]+tC.eta_H[i][j-1])*J[0][i][j] - 0.5*(tC.eta_H[i-1][j]+tC.eta_H[i-1][j-1])*( Bn[2][i-1][j] - Bn[2][i-1][j-1] )/Deltay );
            if(vx_av > 0){
                alphay_im1 = minmod( (B[1][i][j]-B[1][i-2][j])/(Deltax_n2+Deltax_n1), 2*(B[1][i][j]-B[1][i-1][j])/Deltax_n1, 2*(B[1][i-1][j]-B[1][i-2][j])/Deltax_n2 );
                By_up = B[1][i-1][j] + alphay_im1*(Deltax[i-1]/2.);
            }
            else{
                alphay_i = minmod( (B[1][i+1][j]-B[1][i-1][j])/(Deltax_n1+Deltax_1), 2*(B[1][i+1][j]-B[1][i][j])/Deltax_1, 2*(B[1][i][j]-B[1][i-1][j])/Deltax_n1 );
                By_up = B[1][i][j] - alphay_i*(Deltax[i]/2.);
            }
            vy_av = 0.5*( -0.5*(tC.eta_H[i][j]+tC.eta_H[i-1][j])*J[1][i][j] - 0.5*(tC.eta_H[i][j-1]+tC.eta_H[i-1][j-1])*( Bn[2][i-1][j-1] - Bn[2][i][j-1] )/Deltax_n1 );
            if(vy_av > 0){
                alphax_jm1 = minmod( (B[0][i][j]-B[0][i][j-2])/(2*Deltay), 2*(B[0][i][j]-B[0][i][j-1])/Deltay, 2*(B[0][i][j-1]-B[0][i][j-2])/Deltay );
                Bx_up = B[0][i][j-1] + alphax_jm1*(Deltay/2.);
            }
            else{
                alphax_j = minmod( (B[0][i][j+1]-B[0][i][j-1])/(2*Deltay), 2*(B[0][i][j+1]-B[0][i][j])/Deltay, 2*(B[0][i][j]-B[0][i][j-1])/Deltay );
                Bx_up = B[0][i][j] - alphax_j*(Deltay/2.);
            }
            E[2][i][j] = 0.25*(tC.eta_O[i][j] + tC.eta_O[i-1][j] + tC.eta_O[i][j-1] + tC.eta_O[i-1][j-1])*J[2][i][j] - ( vx_av*By_up - vy_av*Bx_up )
                        - ( 0.5*(vc[0][i][j]*B[1][i][j]+vc[0][i-1][j]*B[1][i-1][j]) - 0.5*(vc[1][i][j]*B[0][i][j]+vc[1][i][j-1]*B[0][i][j-1]) );
        }
    }

    return;
}

/*
    Computes electromagnetic force along closed path around cell face. Used to compute time derivative of magnetic field.
    Inputs: E: electric field in reduced units along cell faces
            N_GC: number of ghost cells
            Deltax, Deltay: cell sizes in x and y directions in reduced units
    Output: Qx, Qy: line integral of E divided by cell face area on bottom cell faces (to update Bx) and left cell faces (to update By) respectively
*/
void Compute_EMF(ScalarField & Qx, ScalarField & Qy, VectorField & E, size_t N_GC, std::vector<double> & Deltax, double Deltay)
{

    for(size_t i=N_GC; i<E.shape()[1]-N_GC; i++){
        for(size_t j=N_GC; j<E.shape()[2]-N_GC; j++){

            Qx[i][j] = -( E[2][i][j+1] - E[2][i][j] )/Deltay;

            if( i < E.shape()[1]-N_GC-1 ){
                Qy[i][j] = -( E[2][i][j] - E[2][i+1][j] )/Deltax[i];
            }
        }
    }

    return;
}
