#ifndef BOUNDARY_CONDITIONS_H_INCLUDED
#define BOUNDARY_CONDITIONS_H_INCLUDED

void B_BoundaryConditions(VectorField & B, const BandBCParams & bparams, size_t Ny, size_t N_GC, double t, MPI_Comm comm1D, int world_rank, std::vector<int> & Ny_locs, std::vector<int> & starts, int nbrleft, int nbrright, const Domain & dm);
void Bshear(double By_init, double Bz_init, double y, double t, const Domain & dm, const BandBCParams & bparams, double & Byshear, double & Bzshear);
void By_BC_Calc(std::vector<double> & Bx_BC, std::vector<double> & By_BC, double By_const, size_t Ny, std::vector<int> & Ny_locs, std::vector<int> & starts, int world_rank, MPI_Comm comm1D);
void E_BoundaryConditions(VectorField & E, const BandBCParams & bparams, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright);

#endif // BOUNDARY_CONDITIONS_H_INCLUDED
