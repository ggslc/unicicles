import numpy as np
from amrfile import io as amrio
from pyproj import Transformer
import sys

def parse_commandline(argv):
  import os.path
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument("--input",      help="input, name of BISICLES sample input")
  parser.add_argument("--level",      help="input, degree of BISICLES refinement")
  args = parser.parse_args()

  good_files=0
  for file in [args.input]:
    if file:
      if  os.path.isfile(file):
        good_files+=1
      else:
          print("")
          print("ERROR: specified file does not exist:",file)

  if good_files==1:
    input=args.input
  else:
    print("")
    print("ERROR: problem with one or more input files. I want ALL of them")
    print("")
    parser.print_help()
    sys.exit(2)

  return input, int(args.level)


hdfplotfile, level = parse_commandline(sys.argv)
ncgridfile="cf_gridfile_BISICLES_lev"+str(level)+".nc"

xdumpfile="x_bike.dump"
ydumpfile="y_bike.dump"
 
order=0 # interpolation order, 0 for piecewise constant, 1 for linear

amrID = amrio.load(hdfplotfile)
lo,hi = amrio.queryDomainCorners(amrID, level)

nx=hi[0]-lo[0]+1
ny=hi[1]-lo[1]+1

x_1d,y_1d,bas_t = amrio.readBox2D(amrID, level, lo, hi, "Z_base", order)
x_1d,y_1d,bot_t = amrio.readBox2D(amrID, level, lo, hi, "Z_bottom", order)
x_1d,y_1d,top_t = amrio.readBox2D(amrID, level, lo, hi, "Z_surface", order)
x_1d,y_1d,thk_t = amrio.readBox2D(amrID, level, lo, hi, "thickness", order)

print(top_t.shape)

dx=x_1d[1]-x_1d[0]
dy=y_1d[1]-y_1d[0]

print(dx,x_1d[0],x_1d[-1])
print(dy,y_1d[0],y_1d[-1])

y_2d,x_2d=np.mgrid[y_1d[0]:y_1d[-1]+dy:dy, x_1d[0]:x_1d[-1]+dx:dx]
y_2db1,x_2db1=np.mgrid[y_1d[0]-dy/2:y_1d[-1]+dy-dy/2:dy, x_1d[0]-dx/2:x_1d[-1]+dx-dx/2:dx]
y_2db2,x_2db2=np.mgrid[y_1d[0]-dy/2:y_1d[-1]+dy-dy/2:dy, x_1d[0]+dx/2:x_1d[-1]+dx+dx/2:dx]
y_2db3,x_2db3=np.mgrid[y_1d[0]+dy/2:y_1d[-1]+dy+dy/2:dy, x_1d[0]+dx/2:x_1d[-1]+dx+dx/2:dx]
y_2db4,x_2db4=np.mgrid[y_1d[0]+dy/2:y_1d[-1]+dy+dy/2:dy, x_1d[0]-dx/2:x_1d[-1]+dx-dx/2:dx]

x_2d.dump(xdumpfile)
y_2d.dump(ydumpfile)

#GrIS
inProj = Proj(init = 'epsg:3413') ## ESPG number for Polar Sterographic North (70 degN, 45 degW) projection.
x0 = -654650.0 # sez Steph/amrtocf
y0 = -3385950.0


#BEDMAP2 uses this
##inProj = Proj(init='epsg:3031') #(WGS 84 / Antarctic Polar Stereographic)
#inProj = 'epsg:3031'

#ALL AIS
#x0=-3072000 # sez Steph
#y0 =-3072000

#ASE. Tony sez
#echo -1831000 -904000 | cs2cs +proj=stere +lat_0=-90 +lon_0=0 +datum=WGS84 +x_0=0 +y_0=0 +lat_ts=-71 +units=m +no_defs +to +proj=longlat
#116d16'35.189"W 71d21'55.712"S 0.000
#I matched the ASE subset to bedmap2 coordinate system (ASE origin at -1831000 -904000 in bedmap2 coordinates) and converted to lat/long. This lat/long looks about right.
#is ~4deg off compared to ferret overplot of land, but so is the NEMO tmask vs NEMO GLAMT,GPHIT, so I think OK
#x0=-1831000 
#y0=-904000

#outProj = Proj(init='epsg:4326')
outProj = 'epsg:4326'

##Steph says he applied a x-mirror to make it match up with bedmap2!
##mfact=-1.
mfact=1.
t = Transformer.from_crs(inProj,outProj,always_xy=True)
xl_2d,yl_2d = t.transform((x_2d + x0)*mfact,(y_2d + y0))

xl_2db1,yl_2db1 = t.transform((x_2db1+x0)*mfact,(y_2db1+y0))
xl_2db2,yl_2db2 = t.transform((x_2db2+x0)*mfact,(y_2db2+y0))
xl_2db3,yl_2db3 = t.transform((x_2db3+x0)*mfact,(y_2db3+y0))
xl_2db4,yl_2db4 = t.transform((x_2db4+x0)*mfact,(y_2db4+y0))

##Ferret, at least, says it's all off by ~90 lon as well as reflected in x!
##steph_offset=90.
steph_offset=0.
xl_2d=xl_2d+steph_offset
xl_2db1=xl_2db1+steph_offset
xl_2db2=xl_2db2+steph_offset
xl_2db3=xl_2db3+steph_offset
xl_2db4=xl_2db4+steph_offset

#there's a clash between cf and amrfile HDF5 shared libraries - do cf second
import cf
f = cf.Field()
#
dimx = cf.DimensionCoordinate(data = cf.Data(x_1d, 'm'),properties = {'axis': 'X'})
f.set_construct(cf.DomainAxis(size=len(x_1d)),key="X")
f.set_construct(dimx,axes="X")

#
dimy = cf.DimensionCoordinate(data = cf.Data(y_1d, 'm'),properties = {'axis': 'Y'})
f.set_construct(cf.DomainAxis(size=len(y_1d)),key="Y")
f.set_construct(dimy,axes="Y")

print(y_1d.shape, x_1d.shape)
#
lats = cf.AuxiliaryCoordinate(data=cf.Data(yl_2d, 'degrees_north'), properties={'standard_name': 'latitude'})
lats_bounds = np.empty(lats.shape+(4,), dtype=float)
lats_bounds[...,0] = yl_2db1
lats_bounds[...,1] = yl_2db2
lats_bounds[...,2] = yl_2db3
lats_bounds[...,3] = yl_2db4
lats.set_bounds(cf.Bounds(data=lats_bounds))
f.set_construct(lats, axes=('Y', 'X'))
#
lons = cf.AuxiliaryCoordinate(data=cf.Data(xl_2d, 'degrees_east') , properties={'standard_name': 'longitude'})
lons_bounds = np.empty(lons.shape+(4,), dtype=float)
lons_bounds[...,0] = xl_2db1
lons_bounds[...,1] = xl_2db2
lons_bounds[...,2] = xl_2db3
lons_bounds[...,3] = xl_2db4
lons.set_bounds(cf.Bounds(data=lons_bounds))
f.set_construct(lons, axes=('Y', 'X'))

f.set_data(cf.Data(top_t,units='m'), axes=('Y', 'X'))
#f.set_data(cf.Data(bas_t,units='m'), axes=('Y', 'X'))
#f.set_data(cf.Data(bot_t,units='m'), axes=('Y', 'X'))
#f.set_data(cf.Data(thk_t,units='m'), axes=('Y', 'X'))

cf.write(f,ncgridfile)
