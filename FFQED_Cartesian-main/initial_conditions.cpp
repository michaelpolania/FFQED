#include <vector>
#include "boost/multi_array.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

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

    double Lx = domain.Lx;
    double Ly = domain.Ly;

    double B_pol_max = bparams.B_pol_init;
    double theta_B = bparams.theta_B;
    double B_tor_max = bparams.B_tor_max;
    double x_center = bparams.B_tor_x_center;
    double x_width = bparams.B_tor_x_width;
    double y_center = bparams.B_tor_y_center;
    double y_width = bparams.B_tor_y_width;

    struct BConfig_params params = {B_pol_max, theta_B, B_tor_max, x_center, x_width, y_center, y_width};

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double result1, result2, error;

    gsl_function Fx;
    Fx.function = &InitialBx;
    Fx.params = &params;

    gsl_function Fy;
    Fy.function = &InitialBy;
    Fy.params = &params;

    gsl_function Fz1;
    Fz1.function = &InitialBz1;
    Fz1.params = &params;

    gsl_function Fz2;
    Fz2.function = &InitialBz2;
    Fz2.params = &params;

    //Compute cell face-averaged value of Bx, By, Bz for every cell except for the top, given the functions InitialBx, InitialBy, InitialBz1-2, which computes the analytic initial profile for Bx, By, Bz
    for(size_t i=0; i<x.size(); i++){
        for(size_t j=0; j<y.size(); j++){
            gsl_integration_qag(&Fx, y[j]-Deltay/2., y[j]+Deltay/2., 0, 1e-7, 1000, 3, w, &result1, &error);
            B[0][i+N_GC][j+N_GC] = result1/Deltay;
            if(i < x.size()-1){
                gsl_integration_qag(&Fy, x[i]-Deltax[i]/2., x[i]+Deltax[i]/2., 0, 1e-7, 1000, 3, w, &result1, &error);
                B[1][i+N_GC][j+N_GC] = result1/Deltax[i];
                gsl_integration_qag(&Fz1, x[i]-Deltax[i]/2., x[i]+Deltax[i]/2., 0, 1e-7, 1000, 3, w, &result1, &error);
                gsl_integration_qag(&Fz2, y[j]-Deltay/2., y[j]+Deltay/2., 0, 1e-7, 1000, 3, w, &result2, &error);
                B[2][i+N_GC][j+N_GC] = result1*result2/(Deltax[i]*Deltay);
            }
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
    double y_center = p -> y_center;
    double y_width = p -> y_width;

    return exp(-pow(y-y_center,2.)/(2.*y_width*y_width));
}
