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
void Compute_J(VectorField & B, VectorField & E, VectorField & H, VectorField & D, ScalarField & Rho, VectorField & J, size_t N_GC, const Domain & dm);
void Compute_vc(VectorField & B, VectorField & vc, VectorField & J, TransCoeffs & tC, size_t N_GC, double t, const Domain & dm, const BandBCParams & bparams, int world_rank);
void B_torEvolve(ScalarField & Qz, VectorField & B, VectorField & J, VectorField & vc, TransCoeffs & tC, size_t N_GC, double t, const Domain & dm, const BandBCParams & bparams);
void Compute_E(VectorField & B, VectorField & Bn, VectorField & E, VectorField & vc, VectorField & J, TransCoeffs & tC, size_t N_GC, double t, const Domain & dm, const BandBCParams & bparams);
void Compute_EMF(ScalarField & Qx, ScalarField & Qy, VectorField & E, size_t N_GC, std::vector<double> & Deltax, double Deltay);

#endif // BOUNDARY_CONDITIONS_H_INCLUDED

