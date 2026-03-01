#ifndef FIELD_EVOLUTION_H_INCLUDED
#define FIELD_EVOLUTION_H_INCLUDED

void Compute_J(VectorField & B, VectorField & J, size_t N_GC, const Domain & dm);
double rn_uniform();
double rn_weibull();
void Compute_vc(VectorField & B, VectorField & vc, VectorField & J, TransCoeffs & tC, size_t N_GC, double t, const Domain & dm, const BandBCParams & bparams, int world_rank);
void B_torEvolve(ScalarField & Qz, VectorField & B, VectorField & J, VectorField & vc, TransCoeffs & tC, size_t N_GC, double t, const Domain & dm, const BandBCParams & bparams);
void Compute_E(VectorField & B, VectorField & Bn, VectorField & E, VectorField & vc, VectorField & J, TransCoeffs & tC, size_t N_GC, double t, const Domain & dm, const BandBCParams & bparams);
void Compute_EMF(ScalarField & Qx, ScalarField & Qy, VectorField & E, size_t N_GC, std::vector<double> & Deltax, double Deltay);

#endif // BOUNDARY_CONDITIONS_H_INCLUDED

