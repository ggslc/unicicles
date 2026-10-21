import cf
import numpy as np

nclonlatfile="mesh_mask.nc"
ncgridfile="cf_gridfile.nc"

f=cf.read(nclonlatfile)

#northern_lim_j=500
#lat_g=cf.Data(f[1].array[0:northern_lim_j,:],units='degrees_north')
#lon_g=cf.Data(f[0].array[0:northern_lim_j,:],units='degrees_east')
lat_g=cf.Data(f.select_field("ncvar%gphit").squeeze().array,units='degrees_north')
lon_g=cf.Data(f.select_field("ncvar%glamt").squeeze().array,units='degrees_east')
e1t_g=cf.Data(f.select_field("ncvar%e1t").squeeze().array,units='metres')
e2t_g=cf.Data(f.select_field("ncvar%e2t").squeeze().array,units='metres')

nx=np.shape(lat_g.array)[1]
ny=np.shape(lat_g.array)[0]

g = cf.Field()

x_1d=np.arange(nx)
g.set_construct(cf.DomainAxis(size=len(x_1d)),key="X")
dimx = cf.DimensionCoordinate(data=cf.Data(x_1d, 'm'))
g.set_construct(dimx,axes="X")

y_1d=np.arange(ny)
g.set_construct(cf.DomainAxis(size=len(y_1d)),key="Y")
dimy = cf.DimensionCoordinate(data=cf.Data(y_1d, 'm'))
g.set_construct(dimy,axes="Y")

lats = cf.AuxiliaryCoordinate(data=lat_g, properties={'standard_name': 'latitude'})
g.set_construct(lats, axes=('Y', 'X'))
lons = cf.AuxiliaryCoordinate(data=lon_g, properties={'standard_name': 'longitude'})
g.set_construct(lons, axes=('Y', 'X'))

e1t = cf.AuxiliaryCoordinate(data=e1t_g, properties={'standard_name': 'cell_x_length'})
g.set_construct(e1t, axes=('Y', 'X'))
e2t = cf.AuxiliaryCoordinate(data=e2t_g, properties={'standard_name': 'cell_y_length'})
g.set_construct(e2t, axes=('Y', 'X'))

g.set_data(cf.Data(np.zeros([ny,nx])),axes=('Y', 'X'))
cf.write(g,ncgridfile,fmt="NETCDF4")
