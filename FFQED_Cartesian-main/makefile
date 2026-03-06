CXX=g++
CFLAGS=-O3 -std=c++17
#CFLAGS=-O2 -O3

# Use Homebrew paths dynamically
MPI_PREFIX=$(shell brew --prefix open-mpi)
HDF5_PREFIX=$(shell brew --prefix hdf5-mpi)
BOOST_PREFIX = $(shell brew --prefix boost)
GSL_PREFIX = $(shell brew --prefix gsl)
# Include directories
LDHDIR=-I/opt/homebrew/Cellar/open-mpi/5.0.8/include -I/opt/homebrew/Cellar/gsl/2.8/include -I/opt/homebrew/Cellar/boost/1.89.0/include -I/opt/homebrew/Cellar/hdf5-mpi/1.14.6/include -I/opt/homebrew/opt/openssl@3/include -I/opt/homebrew/opt/zlib/include

# Library directories
LDLDIR=-L/opt/homebrew/Cellar/open-mpi/5.0.8/lib -L/opt/homebrew/Cellar/hdf5-mpi/1.14.6/lib -L/opt/homebrew/Cellar/gsl/2.8/lib -L/opt/homebrew/opt/openssl@3/lib -L/opt/homebrew/opt/zlib/lib

# Libraries to link
LDLIBS=-lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 \
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
	# Compile: only pass include flags and defines here. Library -L flags are for the link stage
	$(CXX) $(CFLAGS) $(LDHDIR) $(LDDEFINES) -c -o $@ $<

main: $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $^ $(LDLDIR) $(LDLIBS)

clean:
	rm -f main *.o