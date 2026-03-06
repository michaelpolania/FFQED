#ifndef BOUNDARY_CONDITIONS_H_INCLUDED
#define BOUNDARY_CONDITIONS_H_INCLUDED

// B field
void B_BoundaryConditions(VectorField & B, const BandBCParams & bparams, size_t Ny, size_t N_GC, double t, MPI_Comm comm1D, int world_rank, std::vector<int> & Ny_locs, std::vector<int> & starts, int nbrleft, int nbrright, const Domain & dm);
void Bshear(double By_init, double Bz_init, double y, double t, const Domain & dm, const BandBCParams & bparams, double & Byshear, double & Bzshear);
void By_BC_Calc(std::vector<double> & Bx_BC, std::vector<double> & By_BC, double By_const, size_t Ny, std::vector<int> & Ny_locs, std::vector<int> & starts, int world_rank, MPI_Comm comm1D);

// D field
void D_BoundaryConditions(VectorField & D, const BandBCParams & bparams, size_t Ny, size_t N_GC, double t, MPI_Comm comm1D, int world_rank, std::vector<int> & Ny_locs, std::vector<int> & starts, int nbrleft, int nbrright, const Domain & dm);
void Dshear(double Dy_init, double Dz_init, double y, double t, const Domain & dm, const BandBCParams & bparams, double & Dyshear, double & Dzshear);
void Dy_BC_Calc(std::vector<double> & Dx_BC, std::vector<double> & Dy_BC, double Dy_const, size_t Ny, std::vector<int> & Ny_locs, std::vector<int> & starts, int world_rank, MPI_Comm comm1D);

// E field
void E_BoundaryConditions(VectorField & E, const BandBCParams & bparams, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright);

// H field
void H_BoundaryConditions(VectorField & H, const BandBCParams & bparams, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright);

#endif // BOUNDARY_CONDITIONS_H_INCLUDED