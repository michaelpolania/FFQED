#ifndef MICROPHYSICS_H_INCLUDED
#define MICROPHYSICS_H_INCLUDED

#include <gsl/gsl_spline.h>

class DataInterpolation{
	public:
	const gsl_spline* n_e_spline;
	const gsl_spline* A_spline;
	const gsl_spline* Z_spline;
	const gsl_spline* n_i_spline;
	const gsl_spline* n_b_spline;
	const gsl_spline* rho_spline;
	const gsl_spline* n_nf_spline;
	const gsl_spline* mu_nf_spline;
	const gsl_spline* gtt_spline;
	const gsl_spline* grr_spline;
	gsl_interp_accel* n_e_acc;
	gsl_interp_accel* A_acc;
	gsl_interp_accel* Z_acc;
	gsl_interp_accel* n_i_acc;
	gsl_interp_accel* n_b_acc;
	gsl_interp_accel* rho_acc;
	gsl_interp_accel* n_nf_acc;
	gsl_interp_accel* mu_nf_acc;
	gsl_interp_accel* gtt_acc;
	gsl_interp_accel* grr_acc;
	double R_cc; //radius of crust-core transition in reduced units
	double RStar; //radius of star in reduced units
	double R_rhocutoff; //radius of star at cutoff density (minimum density considered in crust) in reduced units
};

void load_EOS(std::string CrustEOS_filename, double rho_cutoff, DataInterpolation & Interpolators);

void sigmaCalc(ScalarField & T, ScalarField & n_e, RadialScalarField & A, RadialScalarField & Z, RadialScalarField & n_i, RadialScalarField & n_b, std::string EorT, ScalarField & sigma);

double expE1_rp(double x);

#endif // MICROPHYSICS_H_INCLUDED
