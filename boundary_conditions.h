#ifndef BOUNDARY_CONDITIONS_H_INCLUDED
#define BOUNDARY_CONDITIONS_H_INCLUDED

// Velocity field
struct vConfig_params {
    double v_max;
    double y_center;
    double y_width;
    double f;
};

double Initialvz(double y, double t, void *driver);



// B field
void B_BoundaryConditions(VectorField & B, const BandBCParams & bparams, size_t Ny, size_t N_GC, double t, MPI_Comm comm1D, int world_rank, std::vector<int> & Ny_locs, std::vector<int> & starts, int nbrleft, int nbrright, const Domain & dm);


// D field
void LowerBoundary_D(VectorField& D, const VectorField& V, const VectorField& B, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright, double t, vConfig_params& driver, double y_min, double dy)
void UpperBoundary_D(VectorField& D, const VectorField& V, const VectorField& B, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright, double t, vConfig_params& driver, double y_min, double dy)

// E field
void E_BoundaryConditions(VectorField & E, const BandBCParams & bparams, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright);

// H field
void H_BoundaryConditions(VectorField & H, const BandBCParams & bparams, size_t N_GC, MPI_Comm comm1D, int nbrleft, int nbrright);

#endif // BOUNDARY_CONDITIONS_H_INCLUDED