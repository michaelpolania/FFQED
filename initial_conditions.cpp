#include <vector>
#include "boost/multi_array.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <iostream>

#include "common.h"
#include "initial_conditions.h"

/*
    Initializes the cell face-averaged magnetic field
    Inputs: x, y: vectors of the cell mid-points in each direction
    domain: Domain class object containing the size of the full domain as Lx and Ly
    N_GC: number of ghost cells added in each direction
    Deltax, Deltay: grid spacing in each direction
    Output: B: 3D multi_array whose first dimension corresponds to the three field components and the final two dimensions are the grid points
*/
void InitializeB(std::vector<double> & x, std::vector<double> & y, const BandBCParams& bparams, const Domain & domain, size_t N_GC, std::vector<double> & Deltax, double Deltay, VectorField & B)
{

    double Lx = domain.Lx; //Defines total length of the simulation in x
    double Ly = domain.Ly; 

    double B_pol_max = bparams.B_pol_init;
    double theta_B = bparams.theta_B;
    double B_tor_max = bparams.B_tor_max;
    double x_center = bparams.B_tor_x_center; // Center of Gaussian flux tube in the x-direction
    double x_width = bparams.B_tor_x_width; // Standard deviation of the Gaussian, aka width of flux tube in x
    double y_center = bparams.B_tor_y_center;
    double y_width = bparams.B_tor_y_width;

    struct BConfig_params params = {B_pol_max, theta_B, B_tor_max, x_center, x_width, y_center, y_width};

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double result1, result2, error;

    //Wraps function and parameters into gsl function F_x
    gsl_function Fx;
    Fx.function = &InitialBx;
    Fx.params = &params;

    gsl_function Fy;
    Fy.function = &InitialBy;
    Fy.params = &params;

    //x-dependent part of the Bz component (Gaussian)
    gsl_function Fz1;
    Fz1.function = &InitialBz1;
    Fz1.params = &params;

    gsl_function Fz2;
    Fz2.function = &InitialBz2;
    Fz2.params = &params;

    //Compute cell face-averaged value of Bx, By, Bz for every cell except for the top, given the functions InitialBx, InitialBy, InitialBz1-2, which computes the analytic initial profile for Bx, By, Bz
   
    //iterates over every physical cell, not ghost cells
    for(size_t i=0; i<x.size(); i++){
        for(size_t j=0; j<y.size(); j++){

            //Computes integration for Bx, integrates over y
            gsl_integration_qag(&Fx, y[j]-Deltay/2., y[j]+Deltay/2., 0, 1e-7, 1000, 3, w, &result1, &error);
            B[0][i+N_GC][j+N_GC] = result1/Deltay;
            
            //Ask Peter?
            //For every cell except the last one, compute integration normally (as done for x)
            if(i < x.size()-1){
                gsl_integration_qag(&Fy, x[i]-Deltax[i]/2., x[i]+Deltax[i]/2., 0, 1e-7, 1000, 3, w, &result1, &error);
                B[1][i+N_GC][j+N_GC] = result1/Deltax[i];
                gsl_integration_qag(&Fz1, x[i]-Deltax[i]/2., x[i]+Deltax[i]/2., 0, 1e-7, 1000, 3, w, &result1, &error);
                gsl_integration_qag(&Fz2, y[j]-Deltay/2., y[j]+Deltay/2., 0, 1e-7, 1000, 3, w, &result2, &error);
                B[2][i+N_GC][j+N_GC] = result1*result2/(Deltax[i]*Deltay);
            }
            
            //If at the last one, define By and Bz to be the values computed at the previous cell until BCs are defined
            else{
                B[1][i+N_GC][j+N_GC] = B[1][i+N_GC-1][j+N_GC]; //initialize By continuous across the outer boundary. This gets changed by B_BoundaryConditions, but is needed to get the By BC correct initially.
                B[2][i+N_GC][j+N_GC] = B[2][i+N_GC-1][j+N_GC]; //initialize Bz continuous across the outer boundary. This gets changed by B_BoundaryConditions.
            }
        }
    }

    gsl_integration_workspace_free (w);

    //If initial toroidal field is turned off, then set this component equal to zero everywhere
    if(bparams.B_tor_init == false){
        for(size_t i=0; i<x.size(); i++){
            for(size_t j=0; j<y.size(); j++){
                B[2][i+N_GC][j+N_GC] = 0.;
            }
        }
    }

    return;
}

/*
        Initializes x-component of B field
        Arguments: y: coordinates
                  params: parameters for function (needed by GSL even if empty)
        Output: function value Bx(x,y)
*/
double InitialBx(double y, void * params)
{
    struct BConfig_params *p = (struct BConfig_params *) params;
    double B_pol_max = p -> B_pol_max;
    double theta_B = p -> theta_B;

    return B_pol_max*sin(theta_B);
}

/*
        Initializes y-component of B field
        Arguments: x: coordinate
                  params: parameters for function (needed by GSL even if empty)
        Output: function value By(x,y)
*/
double InitialBy(double x, void * params)
{
    struct BConfig_params *p = (struct BConfig_params *) params;
    double B_pol_max = p -> B_pol_max;
    double theta_B = p -> theta_B;

    return B_pol_max*cos(theta_B);
}

/*
        Initializes z-component of B field
        Arguments: x: coordinates
                  params: parameters for function (needed by GSL even if empty)
        Output: function value B1(x) where Bz(x,y) = B1(x)*B2(y)
*/
double InitialBz1(double x, void * params)
{
    struct BConfig_params *p = (struct BConfig_params *) params;
    double B_tor_max = p -> B_tor_max;
    double x_center = p -> x_center;
    double x_width = p -> x_width;

    return B_tor_max*exp(-pow(x-x_center,2.)/(2.*x_width*x_width));
}

/*
        Initializes z-component of B field
        Arguments: y: coordinates
                  params: parameters for function (needed by GSL even if empty)
        Output: function value Bz2(y) where Bz(x,y) = Bz1(x)*Bz2(y)
*/
double InitialBz2(double y, void * params)
{
    struct BConfig_params *p = (struct BConfig_params *) params;
    double B_tor_max = p -> B_tor_max;
    double y_center = p -> y_center;
    double y_width = p -> y_width;

    return B_tor_max * exp(-pow(y-y_center,2.)/(2.*y_width*y_width));

}

/*
Initializes velocity field (V)

Arguments: x, y coordinates

Output: three velocity components across domain
*/

void InitializeV(std::vector<double> &x, std::vector<double> &y, size_t N_GC, VectorField &V)
{
    for (size_t i = 0; i < x.size(); i++) {
        for (size_t j = 0; j < y.size(); j++) {

            V[0][i+N_GC][j+N_GC] = 0.0; // vx
            V[1][i+N_GC][j+N_GC] = 0.0; // vy
            V[2][i+N_GC][j+N_GC] = 0.0; // vz

        }
    }
}

/*
Initializes displacement field (D)

Arguments: x, y coordinates

Output: three displacement field components across domain
*/
void InitializeD(std::vector<double> &x, std::vector<double> &y, size_t N_GC, VectorField &D, const Domain & dm, double x_min, double y_min)
{
    double E1 = 1.0;
    double E2 = 1.0;

    for (size_t i = 0; i < x.size(); i++) {
        
        for (size_t j = 0; j < y.size(); j++) {

            //int ii = static_cast<int>(i);
            //int jj = static_cast<int>(j);
            //int ng = static_cast<int>(N_GC);

            //double x_i = x_min + (ii - ng) * dm.Deltax[i];
            //double y_j = y_min + (jj - ng) * dm.Deltay;
            
            double x_i = x[i];
            double y_j = y[j];
            
            
            //double y_j = y_min + (j - N_GC) * dm.Deltay;
            //double x_i = x_min + (i - N_GC) * dm.Deltax[0];
            
            D[0][i+N_GC][j+N_GC] = 0.0; // Dx
            D[1][i+N_GC][j+N_GC] = E1 * y_j; // Dy
            D[2][i+N_GC][j+N_GC] = E2 * x_i; // Dz  

            std::cout << "y[0]= " << x[0]<< std::endl;
            
            
            //std::cout << "x_i = " << x_i
          //<< " y_j = " << y_j << std::endl;
            //std::cout << "D = "
            //<< D[0][i+N_GC][j+N_GC] << " "
            //<< D[1][i+N_GC][j+N_GC] << " "
            //<< D[2][i+N_GC][j+N_GC] << std::endl;   
           // if (i == 2 && j == 0) {
             //   std::cout << "D = "
               // << D[0][i+N_GC][j+N_GC] << " "
               // << D[1][i+N_GC][j+N_GC] << " "
                //<< D[2][i+N_GC][j+N_GC] << std::endl;

            //}
            // Print

        }
    }
}

/*
Initializes electric field (E) including ghost cells.

Arguments: x, y coordinates

Output: three electric field components across domain
*/

void InitializeE(std::vector<double> &x, std::vector<double> &y, size_t N_GC, VectorField &E, VectorField &D)
{
    for (size_t i = 0; i < D.shape()[1]; i++) {
        for (size_t j = 0; j < D.shape()[2]; j++) {

            E[0][i][j] = D[0][i][j]; // Ex
            E[1][i][j] = D[1][i][j]; // Ey
            E[2][i][j] = D[2][i][j]; // Ez

            //std::cout << "D = "
            //<< E[0][i+N_GC][j+N_GC] << " "
            //<< E[1][i+N_GC][j+N_GC] << " "
            //<< E[2][i+N_GC][j+N_GC] << std::endl;  

        }
    }
}

/*
Initializes auxillary field (H)

Arguments: x, y coordinates

Output: three auxillary field components across domain
*/

void InitializeH(std::vector<double> &x, std::vector<double> &y, size_t N_GC, const VectorField &B, VectorField &H)
{
    for (size_t i = 0; i < B.shape()[1]; i++) {
        for (size_t j = 0; j < B.shape()[2]; j++) {

            H[0][i][j] = B[0][i][j]; // Hx
            H[1][i][j] = B[1][i][j]; // Hy
            H[2][i][j] = B[2][i][j]; // Hz

        }
    }
}







