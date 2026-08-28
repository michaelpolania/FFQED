#ifndef FIELD_EVOLUTION_H_INCLUDED
#define FIELD_EVOLUTION_H_INCLUDED

struct Fields;
double compute_A1_x(int i, int j, const Fields & f);
double compute_A1_y(int i, int j, const Fields & f);
double compute_A1_z(int i, int j, const Fields & f);
double compute_A2_x(int i, int j, const Fields & f);
double compute_A2_y(int i, int j, const Fields & f);
double compute_A2_z(int i, int j, const Fields & f);
double compute_A3_x(int i, int j, const Fields & f);
double compute_A3_y(int i, int j, const Fields & f);
double compute_A3_z(int i, int j, const Fields & f);
void Compute_Rho(VectorField & D, ScalarField & Rho, const Domain & dm);
void Compute_J(VectorField & B, VectorField & E, VectorField & H, VectorField & D, ScalarField & Rho, VectorField & J, size_t N_GC, const Domain & dm);
void Compute_E(VectorField & B, VectorField & Bn, VectorField & E, VectorField & vc, VectorField & J, TransCoeffs & tC, size_t N_GC, double t, const Domain & dm, const BandBCParams & bparams);
std::pair<double, double> Ez_Flux_Calculation_Qx(int i, int j, VectorField & E, VectorField & B, std::vector<double> & Deltax, double Deltay);
std::pair<double, double> Ez_Flux_Calculation_Qy(int i, int j, VectorField & E, VectorField & B, double Deltax, double Deltay);
double Qz_reconstruction(int i, int j, VectorField & E, VectorField & B, std::vector<double> & Deltax, double Deltay);
std::pair<double, double> Hz_Flux_Calculation_Fx(int i, int j, VectorField & H, VectorField & D, double Deltax, double Deltay);
std::pair<double, double> Hz_Flux_Calculation_Fy(int i, int j, VectorField & H, VectorField & D, double Deltax, double Deltay);
double Fz_reconstruction(int i, int j, VectorField & E, VectorField & B, std::vector<double> & Deltax, double Deltay);
//std::pair<double, double> Ez_Flux_Calculation(int i, int j, VectorField & D, VectorField & B, std::vector<double> & Deltax, double Deltay);
void Compute_RHS(ScalarField & Qx, ScalarField & Qy, ScalarField & Qz, ScalarField & Fx, ScalarField & Fy, ScalarField & Fz, ScalarField & Rho, VectorField & E, VectorField & H, VectorField & D, VectorField & B, const SimParams & params, const Domain & dm, const Process & ps);
double Compute_A_to_cell_center(VectorField & A, int component_index, int i_offset, int j_offset);
double Solve_lambda_Newton(double D_squared, double B_squared, double kappa, int max_iter, double tol);
double Solve_Lambda_Cubic(double D_squared, double B_squared);
void Compute_EH_from_DB(VectorField & E, VectorField & H, VectorField & D, VectorField & B, const SimParams & params,  const Domain & dm);
double slope_calc (VectorField & A, int component_index, int i_offset, int j_offset, int slope_direction, double Deltax, double Deltay);
double Compute_Damping_Term(int i, const Domain & dm);
#endif