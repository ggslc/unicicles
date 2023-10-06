ICE_DIR=/home/robin/off/data-plane/UKESM_ICE
BIKE_DIR=$(ICE_DIR)/UNICICLES/trunk/BISICLES/code/cdriver
HDF_DIR=$(ICE_DIR)/hdf5/parallel
CDF_DIR=$(ICE_DIR)/netcdf/parallel
Z_DIR=$(ICE_DIR)/zlib
GLIM_DIR=$(ICE_DIR)/UNICICLES/trunk//cism-parallel

CDF_INC=-I$(CDF_DIR)/include
CDF_LIBS=-L$(CDF_DIR)/lib -L$(HDF_DIR)/lib -L$(Z_DIR)/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz 
BIKE_LIBS=-L$(BIKE_DIR) -L/usr/lib64 -lBisicles2d.Linux.64.mpicxx.mpif90.DEBUG.OPT.MPI -lChomboLibs2d.Linux.64.mpicxx.mpif90.DEBUG.OPT.MPI -lstdc++ -lpython2.6

GLIM_INC=-I$(GLIM_DIR)/include
GLIM_LIBS=-L$(GLIM_DIR)/lib -lglint -lglide -lglimmer -lglimmer-solve -lglimmer-IO

all:	wrapper

OBJS=wrapper_mod.o gl_mod.o wrapper_main.o

FC=mpif90
FCFLAGS=-DNO_RESCALE

wrapper:  $(OBJS)
	$(FC)  -o wrapper $(OBJS) $(GLIM_LIBS) $(BIKE_LIBS) $(CDF_LIBS)

gl_mod.o:gl_mod.f90
	$(FC) $(FCFLAGS) -g -free $(GLIM_INC) $(CDF_INC) -c gl_mod.f90

wrapper_mod.o:wrapper_mod.f90
	$(FC) $(FCFLAGS) -g -free $(GLIM_INC) $(CDF_INC) -c wrapper_mod.f90

wrapper_main.o:wrapper_main.f90
	$(FC) $(FCFLAGS) -g -free -c wrapper_main.f90

clean:
	rm -f *.o *.mod
