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

#endif // INITIAL_CONDITIONS_H_INCLUDED
