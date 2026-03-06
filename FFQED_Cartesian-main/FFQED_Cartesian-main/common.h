#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#include <mpi.h>
#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

#include "boost/multi_array.hpp"

extern const double pi;
extern const double c; //speed of light in cm/s
extern const double hbarc; //hbar times c in MeV*fm
extern const double unit_e; //elementary charge in units of statcoulomb
extern const double yr; //1 year in seconds
extern const double G; //Gravitational extern constant in dyn*cm^2/g^2
extern const double k_B; //Boltzmann extern constant in MeV/K
extern const double k_Bcgs; //Boltzmann extern constant in cgs units (erg/K)
extern const double M_e; //electron mass in MeV
extern const double M_m; //Muon mass in MeV (105.658375 MeV)
extern const double M_N; //nucleon mass in MeV
extern const double M_n; //neutron mass in MeV
extern const double MeVtoErg; //conversion factor from MeV to erg
extern const double M_solar; //solar mass in g
extern const double alpha_e; //electromagnetic fine structure extern constant (dimensionless)
extern const double eB_crit; //critical magnetic field times elementary charge in MeV^2
extern const double Cstruc; //lattice structure extern constant for bcc lattice
extern const double ComptonWL; //Reduced electron Compton wavelength hbar*c/(m_e*c^2) in fm
extern const double amu; //1 amu in MeV
extern const double microU; //1 micro-u (10^{-6} atomic mass units) in MeV
extern const double gammaEM; //Euler-Mascheroni extern constant
extern const double sigma_SB; //Stefan-Boltzmann extern constant in erg/cm^2/s/K^4

extern const double GaussConverter; //Conversion factor between 1 statCoulomb*Gauss to __ MeV/fm/hbarc = ___ fm^-2
extern const double GaussConverter2; //Conversion factor between 1 statCoulomb*Gauss to __ MeV/fm*hbarc = ___ MeV^2
extern const double GaussConverter3; //Conversion factor: 1 sqrt(MeV/fm^3) = 4.002719868e16 G
extern const double GaussConverter4; //Conversion factor: 1 MeV^2 = 1.444027592e13 G

extern const double B_0; //characteristic magnetic field (G)
extern const double n_e0; //characteristic electron density (fm^{-3})
extern const double L_0; //characteristic length scale (cm)
extern const double t_0; //characteristic timescale (s)
extern const double T_0; //characteristic temperature (K)
extern const double s_0; //characteristic entropy density (erg/K/cm^3)
extern const double E_0; //characteristic electric field (statV/cm)

typedef boost::multi_array<double, 3> VectorField; //type definition for three-component vector fields defined over two spatial dimensions. First index is components (0=x, 1=y, 2=z), second and third indices are x and y coordinates
typedef boost::multi_array_types::index_range range; //range of indices used in array slicing
typedef boost::multi_array<double, 2> ScalarField; //type definition for scalar fields defined over two spatial dimensions. Indices are x and y coordinates
typedef boost::multi_array<double, 1> RadialScalarField; //type definition for scalar fields defined over one ("radial" or x) spatial dimension. Index is x coordinate

struct Process {
    int order, world_rank, nbrleft, nbrright;
    size_t MyS, MyE;
    MPI_Comm comm1D;
};

struct Domain {
    size_t Nx, Ny, N_GC;
    double Lx, Ly, Deltay;
    std::vector<double> Deltax;
    std::vector<double> x, y; //cell centers in reduced units
    double Deltat;
    std::vector<int> Ny_locs, starts;
};

struct SimParams {

    std::string CrustEOS; //File name for crust EOS data table
    int RK_order = 0; //order of Runge-Kutta timestep method (2 or 3)
    bool varying_mesh = false; //true for varying cell-center spacing in "radial"- (x-)direction or false otherwise
    size_t Nx = 0, Ny = 0; //resolution in x-and y-directions
    double x_min = -1., x_max = -1., y_min = 1., y_max = -1.; //limits of simulation domain in reduced units
    double t_max = 0.; //maximum time in reduced units
    double k_C = 0.; //Courant number
    double rho_cutoff = 0.; //cutoff (energy) density in g/cm^3
    double temperature = 0.; //temperature of crust in K
    size_t saves_number = 0; //maximum number of snapshots to save to H5 file
    size_t ECons_cadence = 0; //cadence to print energy conservation information to terminal
    std::string OutputFile; //File name for output H5 file (excluding file extension)
    bool divBCheck = false; //Explicitly show div(B) check every ECons_cadence. Defaults to false if unspecified.

};

struct BandBCParams {

    // ---- B field initial conditions ----
    double B_pol_init = 1., theta_B = 1.5707963; //initial uniform poloidal field in reduced units and angle with respect to y-axis in radians (pi/2 = x-direction)
    bool B_tor_init = true; //whether to include an initial toroidal magnetic field
    double B_tor_max = 0.1; //maximum value of toroidal magnetic field in reduced units (default 1e12 G / B_0 = 0.1)
    double B_tor_x_center = 0., B_tor_x_width = 0., B_tor_y_center = 0., B_tor_y_width = 0.; //center and width of initial toroidal magnetic field in reduced units

    // ---- B field boundary conditions ----
    std::string B_perp_lower, B_parallel_lower, B_pol_upper, B_tor_upper;
    std::string B_parallel_shear_type = "uniform"; //shearing BC type at lower boundary
    double B_parallel_lower_shear_By = 0., B_parallel_lower_shear_Bz = 0.; //magnitude of shearing field at lower boundary in reduced units

    // ---- D field initial conditions ----
    double D_pol_init = 0., theta_D = 1.5707963; //FIX: was theta_B (duplicate). Initial uniform poloidal displacement field in reduced units and angle with respect to y-axis in radians
    bool D_tor_init = false; //whether to include an initial toroidal displacement field
    double D_tor_max = 0.; //FIX: was 0/B_0 (division by zero). Maximum value of toroidal displacement field in reduced units
    double D_tor_x_center = 0., D_tor_x_width = 0., D_tor_y_center = 0., D_tor_y_width = 0.; //center and width of initial toroidal displacement field in reduced units

    // ---- D field boundary conditions ----
    std::string D_perp_lower, D_parallel_lower, D_pol_upper, D_tor_upper;
    std::string D_parallel_shear_type = "uniform"; //FIX: added missing shearing BC type for D
    double D_parallel_lower_shear_Dy = 0., D_parallel_lower_shear_Dz = 0.; //FIX: added missing shearing magnitudes for D

    // ---- E field boundary conditions ----
    std::string E_perp_lower, E_parallel_lower, E_perp_upper, E_parallel_upper;

    // ---- H field boundary conditions ----
    std::string H_perp_lower, H_parallel_lower, H_perp_upper, H_parallel_upper;

    // ---- Toroidal velocity shear ----
    std::string tor_vel_shear = "none"; //which toroidal lattice velocity shear to include. Defaults to "none"
    bool pol_vel = false; //whether to include a poloidal velocity field. Defaults to false
    double t_w = 0., t_m = 0., tor_n = 0.; //time width, centered time, steepness of toroidal velocity shear
    double vc_mag = 0., xwidth_tor = 0., ywidth_tor = 0., x0_tor = 0.5, y0_tor = 0.; //properties of toroidal velocity shear

};

/*
    TransCoeffs class

    Members are ScalarFields containing transport coefficients

    Class objects are instantiated with x and y Nx and Ny sizes of the ScalarFields
*/
class TransCoeffs {
    public:
        size_t Nx, Ny;

    ScalarField eta_H;
    ScalarField deta_Hdx;
    ScalarField eta_O;
    RadialScalarField shear_mod;
    RadialScalarField rho;

    TransCoeffs(size_t Nx, size_t Ny) :
        eta_H(boost::extents[Nx][Ny]),
        deta_Hdx(boost::extents[Nx][Ny]),
        eta_O(boost::extents[Nx][Ny]),
        shear_mod(boost::extents[Nx]),
        rho(boost::extents[Nx])
    {}
};

template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

void load_params(SimParams & params, BandBCParams & bparams, int world_rank);
void MPE_Decomp1D(size_t N, int num_procs, int MyID, size_t & s, size_t & e);
void exchng2Vector(VectorField & A, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright);
void exchng2Scalar(ScalarField & A, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright);
size_t save_cadenceCalc(double Deltat, double t_max, size_t saves_number);
double minmod(double a, double b, double c);
double TrapezoidIntegrator(std::vector<double> & A, double dx);

#endif // COMMON_H_INCLUDED
