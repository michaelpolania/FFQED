#CXX=g++
CXX=mpicxx
CFLAGS=-O3 -std=c++17
#CFLAGS=-O2 -O3

# Use CXXFLAGS for C++ code (-std=c++17 belongs here, not in CFLAGS)
CXXFLAGS = -O3 -std=c++17
#CXXFLAGS = -g -O0 -std=c++17

ifdef MPI_DEBUG
CXXFLAGS += -DMPI_DEBUG
endif
# Use Homebrew's stable, un-versioned symlink paths
LDHDIR = -I/opt/homebrew/opt/open-mpi/include \
         -I/opt/homebrew/opt/hdf5-mpi/include \
         -I/opt/homebrew/opt/boost/include \
         -I/opt/homebrew/opt/gsl/include \
         -I/opt/homebrew/opt/openssl@3/include \
         -I/opt/homebrew/opt/zlib/include

LDLDIR = -L/opt/homebrew/opt/open-mpi/lib \
         -L/opt/homebrew/opt/hdf5-mpi/lib \
         -L/opt/homebrew/opt/boost/lib \
         -L/opt/homebrew/opt/gsl/lib \
         -L/opt/homebrew/opt/openssl@3/lib \
         -L/opt/homebrew/opt/zlib/lib

LDLIBS = -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 \
         -lgsl -lgslcblas \
         -lmpi -lcrypto -lcurl -lpthread -lz -lm

# Combine the include folders into your compilation flags
CXXFLAGS += $(LDHDIR)
# Use Homebrew paths dynamically
#MPI_PREFIX=$(shell brew --prefix open-mpi)
#HDF5_PREFIX=$(shell brew --prefix hdf5-mpi)
#BOOST_PREFIX = $(shell brew --prefix boost)
#GSL_PREFIX = $(shell brew --prefix gsl)
# Include directories
#LDHDIR=-I/opt/homebrew/Cellar/open-mpi/5.0.8/include -I/opt/homebrew/Cellar/gsl/2.8/include -I/opt/homebrew/Cellar/boost/1.89.0/include -I/opt/homebrew/Cellar/hdf5-mpi/1.14.6/include -I/opt/homebrew/opt/openssl@3/include -I/opt/homebrew/opt/zlib/include

# Library directories
#LDLDIR=-L/opt/homebrew/Cellar/open-mpi/5.0.8/lib -L/opt/homebrew/Cellar/hdf5-mpi/1.14.6/lib -L/opt/homebrew/Cellar/gsl/2.8/lib -L/opt/homebrew/opt/openssl@3/lib -L/opt/homebrew/opt/zlib/lib

# Libraries to link
#LDLIBS=-lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 \
       -lgsl -lgslcblas \
       -lmpi -lcrypto -lcurl -lpthread -lz -lm



# Optional defines (kept from your original)
LDDEFINES=-DOLD_HEADER_FILENAME -DHDF_NO_NAMESPACE -DNO_STATIC_CAST

# Quick check for HDF5 C++ header. If it's missing, stop with a helpful message.
H5_HEADER=$(HDF5_PREFIX)/include/H5Cpp.h
H5_FOUND=$(wildcard $(H5_HEADER))
ifeq ($(H5_FOUND),)
$(warning HDF5 header not found at $(H5_HEADER))
$(warning If you haven't installed HDF5, try: brew install hdf5 )
$(warning If HDF5 is installed in a different prefix, set HDF5_PREFIX in this makefile.)
endif

DEPS= common.h initial_conditions.h boundary_conditions.h conservation_checks.h microphysics.h field_evolution.h
OBJ= main.o common.o initial_conditions.o boundary_conditions.o conservation_checks.o microphysics.o field_evolution.o

all: main ;

%.o: %.cpp $(DEPS)
	# CHANGED: $(CFLAGS) replaced with $(CXXFLAGS)
	$(CXX) $(CXXFLAGS) $(LDDEFINES) -c -o $@ $<

main: $(OBJ)
	# CHANGED: $(CFLAGS) replaced with $(CXXFLAGS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLDIR) $(LDLIBS)

clean:
	rm -f main *.o