# UNICICLES

## Build notes


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


