#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_expint.h>

#include "common.h"
#include "microphysics.h"

/*
    Loads crust equation of state file and generates interpolating function for the electron number density
    Inputs: filename: name of EOS file, including file extension
            rho_cutoff: low density cutoff for EOS in g/cm^3. xmax will thus correspond to the radius at which rho=rho_cutoff
    Output: Interpolators: object of DataInterpolation class to contain EOS interpolators and required auxiliary objects
*/
void load_EOS(std::string CrustEOS_filename, double rho_cutoff, DataInterpolation & Interpolators)
{
    double km_to_cm = 1e5; //converts radius in km to cm

	//Determines number of rows in data table
	std::ostringstream filename;
	filename << CrustEOS_filename;
	std::fstream Datafile;
    Datafile.open(filename.str().c_str(),std::ifstream::in);
	std::string line, dummyLine;
    size_t N_dat = 0; //number of rows in data table. Initialize to zero
	getline(Datafile, dummyLine); //skips header (1 line)
	while (std::getline(Datafile, line))
		++N_dat;

	Datafile.close();

    //vectors to store values from data table
	std::vector<double> r_InterpVec; //radius in cm. Vector; will convert to array
	std::vector<double> gtt_InterpVec; //Metric factor g_{tt}=-e^{nu}. Vector; will convert to array
	std::vector<double> grr_InterpVec; //Metric factor g_{rr}=e^{lambda}. Vector; will convert to array
	std::vector<double> n_b_InterpVec; //Baryon number density in fm^{-3}. Vector; will convert to array
    std::vector<double> rho_InterpVec; //mass density in g/cm^3. Vector; will convert to array. Note: not equal to the mass-energy density called rho in the data table!
	std::vector<double> n_i_InterpVec; //ion number density in fm^{-3}. Vector; will convert to array
	std::vector<double> A_InterpVec; //Mass number of crustal nuclei (excluding dripped neutrons and protons). Vector; will convert to array
	std::vector<double> Z_InterpVec; //Atomic number of crustal nuclei (excluding dripped protons). Vector; will convert to array
	std::vector<double> n_e_InterpVec; //Electron number density in fm^{-3}. Vector; will convert to array
	std::vector<double> n_nf_InterpVec; //Number density of free (dripped) neutrons in fm^{-3}. Vector; will convert to array
    std::vector<double> mu_nf_InterpVec; //Free (dripped) neutron chemical potential in MeV. Vector; will convert to array

    //temporary variables to hold data table values
    double r_dt, gtt_dt, grr_dt, n_b_dt, A_dt, Z_eq_dt, Z_dt, n_e_dt, n_nf_dt, mu_nf_dt;
    double dummy; //dummy variable to avoid saving all quantities from loaded data table

	Datafile.open(filename.str().c_str(),std::fstream::in);
	std::getline(Datafile, dummyLine); //skips header (1 line)
	//Gets background TOV solution values from data table
	for(size_t i=0; i<N_dat; i++)
	{
		Datafile >> std::setprecision(15) >> r_dt >> dummy >> dummy >> dummy >> gtt_dt >> grr_dt >> n_b_dt >> dummy >> A_dt >> Z_eq_dt >> Z_dt >> n_e_dt >> n_nf_dt >> mu_nf_dt >> dummy >> dummy;

		r_InterpVec.push_back(r_dt*km_to_cm/L_0);
		gtt_InterpVec.push_back(gtt_dt);
		grr_InterpVec.push_back(grr_dt);
		n_b_InterpVec.push_back(n_b_dt);
		rho_InterpVec.push_back(M_n*MeVtoErg*pow(1e13,3)/pow(c,2)*n_b_dt);
		n_i_InterpVec.push_back(n_e_dt/Z_eq_dt);
		A_InterpVec.push_back(A_dt);
		Z_InterpVec.push_back(Z_dt);
		n_e_InterpVec.push_back(n_e_dt);
		n_nf_InterpVec.push_back(n_nf_dt);
		mu_nf_InterpVec.push_back(mu_nf_dt);
	}
	Datafile.close();

    //std::cout << r_InterpVec[0] << std::endl;
    //std::cout << r_InterpVec.back() << std::endl; 

	//converts vectors of data to arrays to use in interpolation
	//double* r_Interp = &r_InterpVec[0]; //radius in cm

    double r_Interp[N_dat];
    int i = 0;

    //std::cout << "Calculation of r_Interp[i] - r_Interp[i-1]" << std::endl;
    for(int i=0; i<N_dat; i++){
        r_Interp[i] = r_InterpVec[i];
        //std::cout << std::setprecision(15) << i << "      " << r_Interp[i] - r_Interp[i-1] << std::endl;
    }

    //std::cout << "End of r_Interp[i] - r_Interp[i-1] calculation" <<  std::endl;
    //std::cout << "Begin printing r_InterpVec values" << std::endl;

    for(int i=0; i<N_dat; i++){
        r_Interp[i] = r_InterpVec[i];
        //std::cout << std::setprecision(15) << i << "      " << r_Interp[i] << std::endl;
    }
    

	//double* r_Interp = r_InterpVec.data(); //radius in cm
    double* gtt_Interp = &gtt_InterpVec[0]; //Metric factor g_{tt}=-e^{nu}
	double* grr_Interp = &grr_InterpVec[0]; //Metric factor g_{rr}=e^{lambda}
	double* n_b_Interp = &n_b_InterpVec[0]; //Baryon number density in fm^{-3}
    double* rho_Interp = &rho_InterpVec[0]; //mass density in g/cm^3
    double* n_i_Interp = &n_i_InterpVec[0]; //ion number density in fm^{-3}
	double*	A_Interp = &A_InterpVec[0]; //Mass number of crustal nuclei
	double* Z_Interp = &Z_InterpVec[0]; //Atomic number of crustal nuclei
	double* n_e_Interp = &n_e_InterpVec[0]; //Electron number density in fm^{-3}
	double* n_nf_Interp= &n_nf_InterpVec[0]; //Number density of free (dripped) neutrons in fm^{-3}
    double* mu_nf_Interp= &mu_nf_InterpVec[0]; //Chemical potential of free (dripped) neutrons in MeV

    //std::cout << r_InterpVec.back() << std::endl;
    //std::cout << r_Interp[N_dat - 1] << std::endl;
    //std::cout << N_dat << std::endl;

	//Generates spline and acc objects for each quantity to interpolate them and assigns them to DataInterpolation object
	gsl_interp_accel *acc_gtt = gsl_interp_accel_alloc ();
	gsl_spline *spline_gtt = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_gtt, r_Interp, gtt_Interp, N_dat);
	Interpolators.gtt_spline = spline_gtt;
	Interpolators.gtt_acc = acc_gtt;

	gsl_interp_accel *acc_grr = gsl_interp_accel_alloc ();
	gsl_spline *spline_grr = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_grr, r_Interp, grr_Interp, N_dat);
	Interpolators.grr_spline = spline_grr;
	Interpolators.grr_acc = acc_grr;

	gsl_interp_accel *acc_n_b = gsl_interp_accel_alloc ();
	gsl_spline *spline_n_b = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_n_b, r_Interp, n_b_Interp, N_dat);
	Interpolators.n_b_spline = spline_n_b;
	Interpolators.n_b_acc = acc_n_b;

    gsl_interp_accel *acc_rho = gsl_interp_accel_alloc ();
	gsl_spline *spline_rho = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_rho, r_Interp, rho_Interp, N_dat);
	Interpolators.rho_spline = spline_rho;
	Interpolators.rho_acc = acc_rho;

    gsl_interp_accel *acc_n_i = gsl_interp_accel_alloc ();
	gsl_spline *spline_n_i = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_n_i, r_Interp, n_i_Interp, N_dat);
	Interpolators.n_i_spline = spline_n_i;
	Interpolators.n_i_acc = acc_n_i;

    gsl_interp_accel *acc_A = gsl_interp_accel_alloc ();
	gsl_spline *spline_A = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_A, r_Interp, A_Interp, N_dat);
	Interpolators.A_spline = spline_A;
	Interpolators.A_acc = acc_A;

    gsl_interp_accel *acc_Z = gsl_interp_accel_alloc ();
	gsl_spline *spline_Z = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_Z, r_Interp, Z_Interp, N_dat);
	Interpolators.Z_spline = spline_Z;
	Interpolators.Z_acc = acc_Z;

    gsl_interp_accel *acc_n_e = gsl_interp_accel_alloc ();
	gsl_spline *spline_n_e = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_n_e, r_Interp, n_e_Interp, N_dat);
	Interpolators.n_e_spline = spline_n_e;
	Interpolators.n_e_acc = acc_n_e;

    gsl_interp_accel *acc_n_nf = gsl_interp_accel_alloc ();
	gsl_spline *spline_n_nf = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_n_nf, r_Interp, n_nf_Interp, N_dat);
	Interpolators.n_nf_spline = spline_n_nf;
	Interpolators.n_nf_acc = acc_n_nf;

    gsl_interp_accel *acc_mu_nf = gsl_interp_accel_alloc ();
	gsl_spline *spline_mu_nf = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_mu_nf, r_Interp, mu_nf_Interp, N_dat);
	Interpolators.mu_nf_spline = spline_mu_nf;
	Interpolators.mu_nf_acc = acc_mu_nf;

    //std::cout << "Here 1" << std::endl;
    Interpolators.R_cc = r_InterpVec[0]; //coordinate radius of the crust-core transition in reduced units
	Interpolators.RStar = r_InterpVec.back(); //coordinate radius of the star in reduced units

	std::reverse(rho_InterpVec.begin(), rho_InterpVec.end()); //flip order of enegry density vector (need to do this since rho decreases with r)
	std::reverse(r_InterpVec.begin(), r_InterpVec.end()); //flip order of radial coordinate vector (need to do this since rho decreases with r)
	double* r_Interp_rev = &r_InterpVec[0]; //radius in cm
	double* rho_Interp_rev = &rho_InterpVec[0]; //energy density in g/cm^3

    // (no monotonicity check here - allow GSL to report errors if data are invalid)

    gsl_interp_accel *acc_r = gsl_interp_accel_alloc ();
    gsl_spline *spline_r = gsl_spline_alloc (gsl_interp_steffen, N_dat);
    gsl_spline_init(spline_r, rho_Interp_rev, r_Interp_rev, N_dat);

	Interpolators.R_rhocutoff = gsl_spline_eval(spline_r, rho_cutoff, acc_r); //coordinate radius at cutoff density rho_cutoff in reduced units

    gsl_spline_free (spline_r);
    gsl_interp_accel_free (acc_r);

    //std::cout << "Here 2" << std::endl;

	return;
}


/*
    Computes electrical conductivity in s^{-1}
    Uses results from Potekhin A. Y., Baiko D. A., Haensel P. and Yakovlev D. G., 1999, A&A 346, 345,
    Potekhin A. Y., 1999, A&A 351, 787, Gnedin, O. Y. et al, 2001, MNRAS 324, 725
    and reaction rates and transport in neutron stars by A. Schmitt and P. Shternin
    Use Pearson, J. M. et al, 2018, MNRAS 481, 2994, Eq. (C21) and Table C10 for computing x_nuc=x_n

    Inputs: T: temperature in MeV
            n_e: electron number density in fm^{-3}
            A: mass number of nuclei (excludes dripped neutrons and protons)
            Z: atomic number of nuclei (excludes dripped protons)
            n_i: ion number density in fm^{-3}
            n_b: baryon number density in fm^{-3}
            EorT: "Electrical" (for electrical conductivity) or "Thermal" (for thermal conductivity)
    Output: sigma: ScalarField object
*/
void sigmaCalc(ScalarField & T, ScalarField & n_e, RadialScalarField & A, RadialScalarField & Z, RadialScalarField & n_i, RadialScalarField & n_b, std::string EorT, ScalarField & sigma)
{
    static double un1 = 3.;
    static double un2 = 13.0; //n = -2 frequency moment of bcc Coulombic lattice
    static double hbar = 6.582119569e-22; //reduced Planck constant in MeV*s
    static double amu = 931.49410242; //1 amu in MeV

    static double k_F, mu_e, a_i, Gamma, T_plasma, eta, eta_0;
    static double v_F, beta, q_D, q_i2, d2Pdmu2, k_TF2, w0, x_nuc, s, w1, D;
    static double G_s, w, sw, sw1, sinv, svF2, expi_sw, expi_sww, expi_sw1, expi_sw1w1, Lambda_1, Lambda_2, Lambda0_1, Lambda0_2;
    static double CoulombLogarithmVal, nu_ei;

    for(size_t i=0; i<sigma.shape()[0]; i++){
        for(size_t j=0; j<sigma.shape()[1]; j++){

            k_F = pow(3.*pi*pi*n_e[i][j],1./3.)*hbarc; //Fermi wave number in MeV
            mu_e = sqrt(k_F*k_F+M_e*M_e); //electron chemical potential in MeV

            a_i = pow(3./(4.*pi*n_i[i]),1./3.); //ion sphere radius in fm
            Gamma = alpha_e*hbarc*Z[i]*Z[i]/(T[i][j]*a_i); //Coulomb plasma parameter (dimensionless)
            T_plasma = sqrt(4.*pi*alpha_e*Z[i]*Z[i]*n_i[i]*pow(hbarc,3)/(amu*A[i])); //plasma temperature in MeV
            eta = T[i][j]/T_plasma; //ratio of temperature to plasma temperature
            eta_0 = 0.19/pow(Z[i],1./6.); //

            v_F = sqrt( 1 - pow(M_e/mu_e,2) ); //Fermi velocity in units of c
            beta = pi*alpha_e*Z[i]*v_F; //dimensionless beta parameter
            q_D = sqrt(3.*Gamma)*hbarc/a_i; //inverse Debye screening length in MeV
            q_i2 = pow(q_D,2)*( 1. + 0.06*Gamma )*exp(-sqrt(Gamma)); //q_i^2 in MeV^2
            d2Pdmu2 = sqrt(mu_e*mu_e-M_e*M_e)*mu_e/(pi*pi); //in MeV^2
            k_TF2 = 4*pi*alpha_e*d2Pdmu2; //Thomas-Fermi wave number squared in MeV^2
            w0 = un2*pow(2.*k_F/q_D,2)*( 1. + beta/3. ); //w(2*k_TF2); called "w" in e.g. Potekhin et al 1999
            x_nuc = ( 0.1035 + 1.944*pow(n_b[i],0.5717) )/( 1. + 608.*pow(n_b[i],3.143) ) + 0.0225*Z[i]*pow(n_b[i],1.26); //(0.1045+2.09*n_b**0.586)/(1-513*n_b**3)
            s = ( q_i2 + k_TF2 )/(4*k_F*k_F)*exp(-beta);
            w1 = 14.73*pow(x_nuc,2)*( 1. + sqrt(x_nuc)*Z[i]/13. )*( 1. + beta/3. );
            D = exp(-0.42*sqrt(k_F/(M_e*A[i]*Z[i]))*un1*exp(-9.1*eta));

            if(EorT == "Electrical" ){
                G_s = eta/sqrt(eta*eta + eta_0*eta_0)*( 1. + 0.122*beta*beta ); //G_sigma
            }
            else if(EorT == "Thermal" ){
                G_s = eta/sqrt(eta*eta + eta_0*eta_0)*( 1. + 0.122*beta*beta ) + 0.0105*( 1.-1./Z[i] )*( 1. + pow(v_F,3)*beta )*( 1. + pow(x_nuc,2)*sqrt(2.*Z[i]) )*eta/pow( eta*eta + 0.0081,1.5 ); //G_kappa
            }
            else{
                std::cout << "Invalid value for EorT in sigmaCalc" << std::endl;
            }

            w = w0 + w1;
            sw = s*w; //s*w
            sw1 = s*w1; //s*w1
            sinv = 1./s;
            svF2 = v_F*v_F*s; //s times (v_F/c)^2 = s*(1-(M_e/mu_e^2)

            expi_sw = gsl_sf_expint_Ei(-sw);
            expi_sww = gsl_sf_expint_Ei(-sw-w);
            expi_sw1 = gsl_sf_expint_Ei(-sw1);
            expi_sw1w1 = gsl_sf_expint_Ei(-sw1-w1);

            if(sw < 10.){
                Lambda_1 = ( log(1.+sinv) + 1./(1.+sinv)*(1.-exp(-w)) - (1.+sw)*exp(sw)*(-expi_sw+expi_sww) )/2.;
                Lambda_2 = ( v_F*v_F*(exp(-w)-1.+w)/w - svF2/(1.+sinv)*(1.-exp(-w)) - 2.*svF2*log(1.+sinv) + svF2*(2.+sw)*exp(sw)*(-expi_sw+expi_sww) )/2.;
            }
            else{
                Lambda_1 = ( log(1.+sinv) + 1./(1.+sinv)*(1.-exp(-w)) - (1.+sw)*( expE1_rp(sw)-exp(-w)*expE1_rp(sw+w) ) )/2.;
                Lambda_2 = ( v_F*v_F*(exp(-w)-1.+w)/w - svF2/(1.+sinv)*(1.-exp(-w)) - 2.*svF2*log(1.+sinv) + svF2*(2.+sw)*( expE1_rp(sw)-exp(-w)*expE1_rp(sw+w) ) )/2.;
            }
            if(sw1 < 10.){
                Lambda0_1 = ( log(1.+sinv) + 1./(1.+sinv)*(1.-exp(-w1)) - (1.+sw1)*exp(sw1)*(-expi_sw1+expi_sw1w1) )/2.;
                Lambda0_2 = ( v_F*v_F*(exp(-w1)-1.+w1)/w1 - svF2/(1.+sinv)*(1.-exp(-w1)) - 2.*svF2*log(1.+sinv) + svF2*(2.+sw1)*exp(sw1)*(-expi_sw1+expi_sw1w1) )/2.;
            }
            else{
                Lambda0_1 = ( log(1.+sinv) + 1./(1.+sinv)*(1.-exp(-w1)) - (1.+sw1)*( expE1_rp(sw1)-exp(-w1)*expE1_rp(sw1+w1) ) )/2.;
                Lambda0_2 = ( v_F*v_F*(exp(-w1)-1.+w1)/w1 - svF2/(1.+sinv)*(1.-exp(-w1)) - 2.*svF2*log(1.+sinv) + svF2*(2.+sw1)*( expE1_rp(sw1)-exp(-w1)*expE1_rp(sw1+w1) ) )/2.;
            }

            CoulombLogarithmVal = ((Lambda_1 - Lambda_2)-(Lambda0_1-Lambda0_2))*G_s*D; //Fitted formula for Coulomb logarithm from Gnedin et al 2001. Note we absorb the (v_F/c)^2 factor into the definition of Lambda_2

            nu_ei = 4.*pi*pow(Z[i]*alpha_e,2)*n_i[i]*pow(hbarc,3)/(mu_e*mu_e*pow(v_F,3))*CoulombLogarithmVal/hbar; //effective electron-ion collision frequency in s^{-1}

            sigma[i][j] = ( pow(unit_e,2)*n_e[i][j]*(1e39)/(mu_e*MeVtoErg/pow(c,2)*nu_ei) ); //returns conductivity in s^{-1}. Numerical factors here are converting n_e to cm^{-3} and mu_e to erg/c^2 = g

        }
    }

    return;
}

/*
    Computes exponential integral exp(x)*E1(x)=exp(x)*int_x^{inf}e^{-y}/y dy using rational-polynomial approximation
    Eq. 5.1.55 and 5.1.56 from Abramowitz and Stegun, Handbook of Mathematical Functions (pg. 231)
*/
double expE1_rp(double x)
{
    double expE1;

    static double a_1 = 8.5733287401;
    static double a_2 = 18.0590169730;
    static double a_3 = 8.6347608925;
    static double a_4 = 0.2677737343;
    static double b_1 = 9.5733223454;
    static double b_2 = 25.6329561486;
    static double b_3 = 21.0996530827;
    static double b_4 = 3.9584969228;

    static double c_0 = -0.57721566;
    static double c_1 = 0.99999193;
    static double c_2 = -0.24991055;
    static double c_3 = 0.05519968;
    static double c_4 = -0.00976004;
    static double c_5 = 0.00107857;

    if( x > 1. ){
        expE1 = ( pow(x,4) + a_1*pow(x,3) + a_2*pow(x,2) + a_3*x + a_4 )/( pow(x,4) + b_1*pow(x,3) + b_2*pow(x,2) + b_3*x + b_4 )/x;
    }
    else{
        expE1 = exp(x)*( -log(x) + c_0 + c_1*x + c_2*pow(x,2) + c_3*pow(x,3.) + c_4*pow(x,4) + c_5*pow(x,5) );
    }

    return expE1;
}
