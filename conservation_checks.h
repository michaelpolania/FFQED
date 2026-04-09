#ifndef CONSERVATION_CHECKS_H_INCLUDED
#define CONSERVATION_CHECKS_H_INCLUDED

double divB_Calculator(VectorField & B, size_t N_GC);
void divB_Monitor(VectorField & B, size_t N_GC,  MPI_Comm comm1D, int world_rank, std::vector<int> & Ny_locs);
void EnergyConservation(VectorField & B, VectorField & E, VectorField & J, ScalarField & eta_O, size_t N_GC, double t, const Domain & dm, const BandBCParams & bparams, double & U_B, double & JouleH, double & PoyntingF);

#endif // CONSERVATION_CHECKS_H_INCLUDED
