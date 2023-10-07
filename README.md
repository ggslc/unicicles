# UNICICLES

## Build notes

### Ubuntu 22.04 (serial, debug)

Assumes BISICLES at $BISICLES_HOME/bisicles-uob and Chombo at $BISICLES_HOME/Chombo.
Both built with DEBUG=TRUE, OPT=FALSE, MPI=FALSE, USE_PETSC=FALSE

```

cd $BISICLES_HOME
git clone git@github.com:ggslc/unicicles.git
cd unicicles

# Build glimmer-cism
cd glimmer-cism
./bootstrap
cd ..
mkdir glimbike-serial
cd glimbike-serial

BIKE=$BISICLES_HOME/bisicles-uob/
BIKE_CONFIG=2d.Linux.64.g++.gfortran.DEBUG
HDF5=/usr/lib/x86_64-linux-gnu/hdf5/serial/
NETCDF=/usr

PYTHON=python3 FC=gfortran FCFLAGS="-fno-range-check -ffree-line-length-0 -DBISICLES_CDRIVER -DNO_RESCALE -g -I$BIKE/code/src " LDFLAGS="-L$BIKE/code/lib -lBisicles$BIKE_CONFIG -lChomboLibs$BIKE_CONFIG -lpython3.10 -L$HDF5 -lhdf5 -lz " ../glimmer-cism/configure --with-netcdf=$NETCDF --with-hdf5=$HDF5 --prefix=$PWD --disable-python
make
make install

# build wrappers/ukesm-ice_NETCDF
cd wrappers/ukesm-ice_NETCDF
make clean -f Makefile.ubuntu22.04
make -f Makefile.ubuntu22.04

```

### Ubuntu 22.04 (parallel, opt, petsc, debug)

Assumes BISICLES at $BISICLES_HOME/bisicles-uob and Chombo at $BISICLES_HOME/Chombo.
Both built with DEBUG=TRUE, OPT=TRUE, MPI=TRUE, USE_PETSC=TRUE

```

cd $BISICLES_HOME
git clone git@github.com:ggslc/unicicles.git
cd unicicles

# Build glimmer-cism
cd glimmer-cism
./bootstrap
cd ..
mkdir glimbike-parallel
cd glimbike-parallel

BIKE=$BISICLES_HOME/bisicles-uob/
BIKE_CONFIG=2d.Linux.64.mpiCC.mpif90.DEBUG.OPT.MPI.PETSC
HDF5=/usr/lib/x86_64-linux-gnu/hdf5/openmpi/
NETCDF=/usr

PYTHON=python3 FC=gfortran FCFLAGS="-fno-range-check -ffree-line-length-0 -DBISICLES_CDRIVER -DNO_RESCALE -g -I$BIKE/code/src " LDFLAGS="-L$BIKE/code/lib -lBisicles$BIKE_CONFIG -lChomboLibs$BIKE_CONFIG -lpython3.10 -L$HDF5 -lhdf5 -lz " ../glimmer-cism/configure --with-netcdf=$NETCDF --with-hdf5=$HDF5 --prefix=$PWD --disable-python
make
make install

# build wrappers/ukesm-ice_NETCDF
cd wrappers/ukesm-ice_NETCDF
make clean -f Makefile.ubuntu22.04_opt_mpi_petsc
make -f Makefile.ubuntu22.04_opt_mpi_petsc

```




### Ubuntu 20.04 (serial, debug)

```
cd unicicles

# Build glimmer-cism
cd glimmer-cism
./bootstrap
cd ..
mkdir glimbike-serial
cd glimbike-serial
BIKE=$BISICLES_HOME/bisicles-uob/
HDF5=/usr/lib/x86_64-linux-gnu/hdf5/serial/
FC=gfortran FCFLAGS="-fno-range-check -ffree-line-length-0 -DBISICLES_CDRIVER -DNO_RESCALE -g -I$BIKE/code/src " LDFLAGS="-L$BIKE/code/lib -lBisicles2d.Linux.64.g++.gfortran.DEBUG -lChomboLibs2d.Linux.64.g++.gfortran.DEBUG -lpython3.8 -L$HDF5 -lhdf5 -lz " ../glimmer-cism/configure --with-netcdf=/usr --with-hdf5=$HDF5 --prefix=$PWD
make
make install


```


