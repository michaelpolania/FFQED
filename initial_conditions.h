#ifndef INITIAL_CONDITIONS_H_INCLUDED
#define INITIAL_CONDITIONS_H_INCLUDED

struct BConfig_params{
    double B_pol_max, theta_B, B_tor_max, x_center, x_width, y_center, y_width;
};

void InitializeB(std::vector<double> & x, std::vector<double> & y, const BandBCParams& bparams, const Domain & domain, size_t N_GC, std::vector<double> & Deltax, double Deltay, VectorField & B);
double InitialBx(double y, void * params);
double InitialBy(double x, void * params);
double InitialBz1(double x, void * params);
double InitialBz2(double y, void * params);
void InitializeV(std::vector<double> &x, std::vector<double> &y, size_t N_GC, VectorField &V);
//void InitializeD(std::vector<double> &x, std::vector<double> &y, size_t N_GC, VectorField &D, const Domain & dm, double x_min, double y_min);
void InitializeD(std::vector<double> &x, std::vector<double> &y, size_t N_GC, VectorField &D, const Domain & dm, std::vector<double> &Deltax, double Deltay);
void InitializeE(std::vector<double> &x, std::vector<double> &y, size_t N_GC, VectorField &E, VectorField &D);
void InitializeH(std::vector<double> &x, std::vector<double> &y, size_t N_GC, const VectorField &B, VectorField &H);
void InitializeB_test(size_t N_GC, VectorField &B);
struct EConfig_params;
struct DyIntegrationParams;
double Dy_over_x_integrand(double x, void *params);
double InitialDy(double y, void * params);
double InitialDz(double x, void * params);



#endif // INITIAL_CONDITIONS_H_INCLUDED
