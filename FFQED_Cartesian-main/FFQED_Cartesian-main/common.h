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
extern const double Cstruc; //lattice structure extern constant for bcc lattice (from Baiko, Potekhin and Yakovlev 2001: take their K_m and multiply by geometric factor (4*pi/3)^{1/3})
extern const double ComptonWL; //Reduced electron Compton wavelength hbar*c/(m_e*c^2) in fm (this is the Compton wavelength appearing in Chamel and Stoyanov 2020)
extern const double amu; //1 amu in MeV
extern const double microU; //1 micro-u (10^{-6} atomic mass units) in MeV
extern const double gammaEM; //Euler-Mascheroni extern constant
extern const double sigma_SB; //Stefan-Boltzmann extern constant in erg/cm^2/s/K^4

extern const double GaussConverter; //Conversion factor between 1 statCoulomb*Gauss to __ MeV/fm/hbarc = ___ fm^-2
extern const double GaussConverter2; //Conversion factor between 1 statCoulomb*Gauss to __ MeV/fm*hbarc = ___ MeV^2
extern const double GaussConverter3; //Conversion factor: 1 sqrt(MeV/fm^3) = 4.002719868e16 G
extern const double GaussConverter4; //Conversion factor: 1 MeV^2 = 4.002719868e16/(197.3269804)^(3/2) G = 1.444027592e13 G (hbarc=c=1 all energy units)

extern const double B_0; //characteristic magnetic field (G)
extern const double n_e0; //characteristic electron density (fm^{-3})
extern const double L_0; //characteristic length scale (cm)
extern const double t_0; //characteristic timescale (s). Taken as the Hall time-scale with length scale 10 m and field 10^{15} G
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
    bool varying_mesh = false; //true for varying cell-center spacing in "radial"- (x-)direction or false otherwise. Defaults to false if unspecified.
    size_t Nx = 0, Ny = 0; //resolution in x-and y-directions
    double x_min = -1., x_max = -1., y_min = 1., y_max = -1.; //limits of simulation domain in reduced units
    double t_max = 0.; //maximum time in reduced units
    double k_C = 0.; //Courant number
    double rho_cutoff = 0.; //cutoff (energy) density in g/cm^3 (lowest density to include in simulation domain)
    double temperature = 0.; //temperature of crust in K (used to compute electrical conductivity)
    std::string B_perp_lower, B_parallel_lower, B_pol_upper, B_tor_upper; //boundary conditions for the magnetic field
    size_t saves_number = 0; //maximum number of snapshots to save to H5 file
    size_t ECons_cadence = 0; //cadence to print energy conservation information to terminal
    std::string OutputFile; //File name for output H5 file (excluding file extension)
    bool divBCheck = false; //Explicitly show div(B) check every ECons_cadence. Defaults to false if unspecified.

};

struct BandBCParams {

    double B_pol_init = 1., theta_B = pi/2.; //initial uniform poloidal field in reduced units and angle of this field with respect to the y-axis in radians.
    bool B_tor_init = true; //whether to include an initial toroidal magnetic field in the simulation. Defaults to true if unspecified.
    double B_tor_max = 1e12/B_0; //maximum value of toroidal magnetic field if turned on, in reduced units. Defaults to 1e12 G if unspecified.
    double B_tor_x_center, B_tor_x_width, B_tor_y_center, B_tor_y_width; //center position and width of initial toroidal magnetic field in reduced units
    std::string B_perp_lower, B_parallel_lower, B_pol_upper, B_tor_upper; //boundary conditions for the magnetic field
    std::string E_perp_lower, E_parallel_lower, E_perp_upper, E_parallel_upper; //lower boundary conditions for electric field
    std::string B_parallel_shear_type = "uniform"; //which shearing magnetic field to use at lower boundary if B_parallel_lower = "shearing". Defaults to "uniform" if unspecified.
    double B_parallel_lower_shear_By = 0., B_parallel_lower_shear_Bz = 0.; //magnitude of shearing field at lower boundary in reduced units. Defaults to zero if unspecified.
    std::string tor_vel_shear = "none"; //which toroidal lattice velocity shear to include in the simulation. Defaults to "none"
    bool pol_vel = false; //whether to include a poloidal velocity field in the simulation. Defaults to false if unspecified.
    double t_w = 0., t_m = 0., tor_n = 0.; //time width, centered time, steepness of toroidal velocity shear
    double vc_mag = 0., xwidth_tor = 0., ywidth_tor = 0., x0_tor = 5e4/L_0, y0_tor = 0.; //properties of toroidal velocity shear

};

/*
    TransCoeffs class

    Members are ScalarFields containing transport coefficients

    Class objects are instantiated with x and y Nx and Ny sizes of the ScalarFields
*/
class TransCoeffs {
    public:
        size_t Nx, Ny;
        //bool cond_aniso_In; //input variable for conductivity_anisotropy boolean
        //bool conductivity_anisotropy; //permanent variable for conductivity_anisotropy boolean

    ScalarField eta_H; //
    ScalarField deta_Hdx; //
    ScalarField eta_O; //
    RadialScalarField shear_mod;
    RadialScalarField rho;
//    ScalarField eta_O_delta; //
//    ScalarField eta_O_perp; //
//    ScalarField kappa; //
//    ScalarField kappa_delta; //
//    ScalarField kappa_perp; //
//    ScalarField kappa_H; //

    TransCoeffs(size_t Nx, size_t Ny) :
        eta_H(boost::extents[Nx][Ny]), //
        deta_Hdx(boost::extents[Nx][Ny]), //
        eta_O(boost::extents[Nx][Ny]), //
        shear_mod(boost::extents[Nx]), //
        rho(boost::extents[Nx])
        //eta_O_delta(boost::extents[Nr][Nth]), //
        //eta_O_perp(boost::extents[Nr][Nth]), //
        //kappa(boost::extents[Nr][Nth]), //
        //kappa_delta(boost::extents[Nr][Nth]),//
        //kappa_perp(boost::extents[Nr][Nth]),//
        //kappa_H(boost::extents[Nr][Nth]) //
    {
//        conductivity_anisotropy = cond_aniso_In;
//        if(conductivity_anisotropy == true){
//            //resize isotropic conductivities to zero if unused
//            eta_O.resize(boost::extents[0][0]);
//            kappa.resize(boost::extents[0][0]);
//        }
//        else{
//            //resize anisotropic conductivities to zero if unused
//            eta_O_delta.resize(boost::extents[0][0]);
//            eta_O_perp.resize(boost::extents[0][0]);
//            kappa_delta.resize(boost::extents[0][0]);
//            kappa_perp.resize(boost::extents[0][0]);
//            kappa_H.resize(boost::extents[0][0]);
//        }
    }
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
