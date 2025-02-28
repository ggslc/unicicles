import cf
import numpy as np

nclonlatfile="eORCA025_LONLAT-FILL_v2.3.nc"
ncgridfile="cf_gridfile_eORCA025_v2.3.nc"

f=cf.read(nclonlatfile)

#northern_lim_j=500
#lat_g=cf.Data(f[1].array[0:northern_lim_j,:],units='degrees_north')
#lon_g=cf.Data(f[0].array[0:northern_lim_j,:],units='degrees_east')
lat_g=cf.Data(f[1].array,units='degrees_north')
lon_g=cf.Data(f[0].array,units='degrees_east')

nx=np.shape(lat_g.array)[1]
ny=np.shape(lat_g.array)[0]

g = cf.Field()

x_1d=np.arange(nx)
dimx = cf.DimensionCoordinate(data=cf.Data(x_1d, 'm'), properties={'axis': 'X'})
g.insert_dim(dimx)
y_1d=np.arange(ny)
dimy = cf.DimensionCoordinate(data=cf.Data(y_1d, 'm'), properties={'axis': 'Y'})
g.insert_dim(dimy)

lats = cf.AuxiliaryCoordinate(data=lat_g, properties={'standard_name': 'latitude'})
g.insert_aux(lats, axes=('Y', 'X'))
lons = cf.AuxiliaryCoordinate(data=lon_g, properties={'standard_name': 'longitude'})
g.insert_aux(lons, axes=('Y', 'X'))

g.insert_data(cf.Data(np.zeros([ny,nx])),axes=('Y', 'X'))
cf.write(g,ncgridfile,fmt="NETCDF4")
