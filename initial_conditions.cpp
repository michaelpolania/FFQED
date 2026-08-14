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
                //gsl_integration_qag(&Fz1, x[i]-Deltax[i]/2., x[i]+Deltax[i]/2., 0, 1e-7, 1000, 3, w, &result1, &error);
                //gsl_integration_qag(&Fz2, y[j]-Deltay/2., y[j]+Deltay/2., 0, 1e-7, 1000, 3, w, &result2, &error);
                //B[2][i+N_GC][j+N_GC] = result1*result2/(Deltax[i]*Deltay);
                B[2][i+N_GC][j+N_GC] = 5.;
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
    //if(bparams.B_tor_init == false){
      //  for(size_t i=0; i<x.size(); i++){
        //    for(size_t j=0; j<y.size(); j++){
          //      B[2][i+N_GC][j+N_GC] = 0.;
            //}
        //}
    //}

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
    //struct BConfig_params *p = (struct BConfig_params *) params;
    //double B_pol_max = p -> B_pol_max;
    //double theta_B = p -> theta_B;

    //return B_pol_max*sin(theta_B);
    return 0.;
}

/*
        Initializes y-component of B field
        Arguments: x: coordinate
                  params: parameters for function (needed by GSL even if empty)
        Output: function value By(x,y)
*/
double InitialBy(double x, void * params)
{
    //struct BConfig_params *p = (struct BConfig_params *) params;
    //double B_pol_max = p -> B_pol_max;
    //double theta_B = p -> theta_B;

    //return B_pol_max*cos(theta_B);
    return 3.*x;
}

/*
        Initializes z-component of B field
        Arguments: x: coordinates
                  params: parameters for function (needed by GSL even if empty)
        Output: function value B1(x) where Bz(x,y) = B1(x)*B2(y)
*/
double InitialBz1(double x, void * params)
{
    //struct BConfig_params *p = (struct BConfig_params *) params;
    //double B_tor_max = p -> B_tor_max;
    //double x_center = p -> x_center;
    //double x_width = p -> x_width;

    //return B_tor_max*exp(-pow(x-x_center,2.)/(2.*x_width*x_width));
    return 5;
}

/*
        Initializes z-component of B field
        Arguments: y: coordinates
                  params: parameters for function (needed by GSL even if empty)
        Output: function value Bz2(y) where Bz(x,y) = Bz1(x)*Bz2(y)
*/
double InitialBz2(double y, void * params)
{
    //struct BConfig_params *p = (struct BConfig_params *) params;
    //double B_tor_max = p -> B_tor_max;
    //double y_center = p -> y_center;
    //double y_width = p -> y_width;

    //return B_tor_max * exp(-pow(y-y_center,2.)/(2.*y_width*y_width));
    return 0;

}

void InitializeB_test(size_t N_GC, VectorField &B)
{
    for(size_t i = N_GC; i < B.shape()[1]-N_GC; i++){
        for(size_t j = N_GC; j < B.shape()[2]-N_GC; j++){
            B[0][i][j] = 10.0;  // Bx = 10
            B[1][i][j] = 0.0;   // By = 0
            B[2][i][j] = 0.0;   // Bz = 0
        }
    }
    return;
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

struct EConfig_params {
    double E1;
    double E2;
};

double InitialDy(double y, void * params)
{
    struct EConfig_params *p = (struct EConfig_params *) params;
    return p->E1 * y;
}

// Struct to pass the fixed y-coordinate into the x-integration routine
struct DyIntegrationParams {
    double fixed_y;
    EConfig_params *config;
};

// Integrand for Dy: GSL will vary 'x', but we evaluate using the fixed 'y'
double Dy_over_x_integrand(double x, void *params) {
    DyIntegrationParams *p = (DyIntegrationParams *)params;
    double y = p->fixed_y;
    double E1 = p->config->E1;
    
    return E1 * y; 
}

// Integrand for Dz: GSL varies 'x'
double InitialDz(double x, void *params) {
    EConfig_params *p = (EConfig_params *)params;
    return p->E2 * x;
}

void InitializeD(std::vector<double> &x, std::vector<double> &y, size_t N_GC, VectorField &D, const Domain & dm, std::vector<double> &Deltax, double Deltay)
{
    double E1 = 1.0;
    double E2 = 1.0;
    EConfig_params config = {E1, E2};

    // Define Dz GSL function (direct 1D integration over x)
    gsl_function Fz;
    Fz.function = &InitialDz;
    Fz.params = &config;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
    double result1, result2, error;

    for(size_t i = 0; i < x.size(); i++){
        for(size_t j = 0; j < y.size(); j++){
            D[0][i+N_GC][j+N_GC] = 0.0;

            if(i < x.size()-1){
                double x_min = x[i] - Deltax[i]/2.0;
                double x_max = x[i] + Deltax[i]/2.0;

                // 1. Compute Dy: Setup the struct with the current cell's y-coordinate
                DyIntegrationParams dy_params = { y[j]- Deltay*0.5, &config };
                gsl_function Fy;
                Fy.function = &Dy_over_x_integrand;
                Fy.params = &dy_params;

                // Integrate over x
                gsl_integration_qag(&Fy, x_min, x_max, 0, 1e-7, 1000, 3, w, &result1, &error);
                D[1][i+N_GC][j+N_GC] = result1 / Deltax[i];

                // 2. Compute Dz: Integrate over x
                gsl_integration_qag(&Fz, x_min, x_max, 0, 1e-7, 1000, 3, w, &result2, &error);
                D[2][i+N_GC][j+N_GC] = result2 / Deltax[i];
            }
            else{
                D[1][i+N_GC][j+N_GC] = D[1][i+N_GC-1][j+N_GC]; 
                D[2][i+N_GC][j+N_GC] = D[2][i+N_GC-1][j+N_GC]; 
            }

             if (i == 98 && j == 98){
       
       // std::cout << "A3x_avg = " << (y[j] + y[j-1])/(20) << std::endl;
        //std::cout << "x[97]= " << x[i-1] << std::endl;
        //std::cout << "x[98]= " << x[i] << std::endl;
        //std::cout << "y[97]= " << y[j-1] << std::endl;
       // std::cout << "y[98]= " << y[j] << std::endl;
        //std::cout << Deltax[i] << std::endl;
        //std::cout << Deltay << std::endl;

        //std::cout << "A1z_avg = " << -(1/(40*M_PI)) * 0.5 *(y[j] + y[j-1]) << std::endl;

    }
        }
    }

    gsl_integration_workspace_free(w);

    // Print statements matching your verification logic
    //std::cout << "D[0][5][5] = " << D[0][100][100] << std::endl;  
    //std::cout << "D[1][5][5] = " << D[1][100][100] << std::endl;  
    //std::cout << "D[2][5][5] = " << D[2][100][100] << std::endl;  
    //std::cout << "x[3] " << x[3] << std::endl;
    //std::cout << "y[2] " << y[2] << std::endl;
    //std::cout << "y[3] " << y[3] << std::endl;
    //std::cout << "E1 " << E1 << std::endl;

   
}


/*
Initializes displacement field (D)

Arguments: x, y coordinates

Output: three displacement field components across domain
*/
/*
void InitializeD(std::vector<double> &x, std::vector<double> &y, size_t N_GC, VectorField &D, const Domain & dm, double x_min, double y_min)
{
    double E1 = 1.0;
    double E2 = 1.0;

    //iterates over every physical cell, not ghost cells
    for(size_t i=0; i<x.size(); i++){
        for(size_t j=0; j<y.size(); j++){
    


            //Computes integration for Bx, integrates over y
            gsl_integration_qag(&Fx, y[j]-Deltay/2., y[j]+Deltay/2., 0, 1e-7, 1000, 3, w, &result1, &error);
            B[0][i+N_GC][j+N_GC] = result1/Deltay;



    
    
    
}   }}
   
    for (size_t i = N_GC; i < D.shape()[1]-N_GC; i++) {
        
        for (size_t j = N_GC; j < D.shape()[2]-N_GC; j++) {

            //int ii = static_cast<int>(i);
            //int jj = static_cast<int>(j);
            //int ng = static_cast<int>(N_GC);

            //double x_i = x_min + (ii - ng) * dm.Deltax[i];
            //double y_j = y_min + (jj - ng) * dm.Deltay;
            
            double x_i = x[i-N_GC];
            double y_j = y[j-N_GC];
            
            
            //double y_j = y_min + (j - N_GC) * dm.Deltay;
            //double x_i = x_min + (i - N_GC) * dm.Deltax[0];
            
            D[0][i][j] = 0.0; // Dx
            D[1][i][j] = E1 * y_j; // Dy
            D[2][i][j] = E2 * x_i; // Dz  
            //D[2][i][j] = 0.0;
            //std::cout << "y[98]= " << x[98]<< std::endl;
            std::cout << "y[97]= " << x.size() << std::endl;
            
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

*/
/*
Initializes electric field (E) excluding ghost cells.

Arguments: x, y coordinates

Output: three electric field components across domain
*/

void InitializeE(std::vector<double> &x, std::vector<double> &y, size_t N_GC, VectorField &E, VectorField &D)
{
    for (size_t i = N_GC; i < D.shape()[1]-N_GC; i++) {
        for (size_t j = N_GC; j < D.shape()[2]-N_GC; j++) {

            //Update from 06/29/26
            //Recall, E=D, however, E and D do not live at the same locations, so must interpolate E to cell edges
            //E[0][i][j] =  0.25*(D[0][i][j] + D[0][i+1][j] + D[0][i][j-1] + D[0][i+1][j-1]);
            //E[1][i][j] = 0.25*(D[1][i][j] + D[1][i-1][j] + D[1][i-1][j+1] + D[1][i][j+1]);
            //E[2][i][j] = 0.25*(D[2][i][j] + D[2][i-1][j] + D[2][i][j-1] + D[2][i-1][j-1]);

            //Predominantly cell-centered formalism (E, H live at cell centers)

            E[0][i][j] = 0.5 * (D[0][i][j] + D[0][i+1][j]);
            E[1][i][j] = 0.5 * (D[1][i][j] + D[1][i][j+1]);
            E[2][i][j] = D[2][i][j];
        
        }
    }
}

/*
Initializes auxillary field (H)

Arguments: x, y coordinates

Output: three auxillary field components across domain
*/

void InitializeH(std::vector<double> &x, std::vector<double> &y, size_t N_GC, VectorField &B, VectorField &H)
{
    for (size_t i = N_GC; i < B.shape()[1]-N_GC; i++) {
        for (size_t j = N_GC; j < B.shape()[2]-N_GC; j++) {

            //Recall, H=B, however, H and B do not live at the same locations, so must interpolate H to cell edges
            //H[0][i][j] =  0.25*(B[0][i][j] + B[0][i+1][j] + B[0][i][j-1] + B[0][i+1][j-1]);
            //H[1][i][j] = 0.25*(B[1][i][j] + B[1][i-1][j] + B[1][i-1][j+1] + B[1][i][j+1]);
            //H[2][i][j] = 0.25*(B[2][i][j] + B[2][i-1][j] + B[2][i][j-1] + B[2][i-1][j-1]);

            //Predominantly cell-centered formalism (E, H live at cell centers)

            H[0][i][j] = 0.5 * (B[0][i][j] + B[0][i+1][j]);
            H[1][i][j] = 0.5 * (B[1][i][j] + B[1][i][j+1]);
            H[2][i][j] = B[2][i][j];

        }
    }
}







