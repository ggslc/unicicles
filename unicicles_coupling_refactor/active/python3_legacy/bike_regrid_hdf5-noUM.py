import numpy as np
import sys
import os
from arg_to_file_exist import arg_to_file_exist
import cf

def validate_warn(self,filename=None, warn=True):
  self.validate_errors(filename=filename, warn=True)
  print("(no more than one warning is noted)")
  return

def disable_mule_validators(umf):
  umf.validate_errors=umf.validate
  print("OVERWRITING MULE'S VALIDATION FUNCTION")
  print("(the original is saved as self.validate_errors)")
  funcType = type(umf.validate)
  print("VALIDATION ERRORS DOWNGRADED TO WARNINGS")
  umf.validate=funcType(validate_warn, umf)
  return umf


def parse_commandline(argv):
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument("--hdf5_input",    help = "input, name of unicicles hdf5")
  parser.add_argument("--hdf5_output",   help = "output, name of unicicles melange sink hdf5")
  parser.add_argument("--nc_output",     help = "output, name of nemo iceberg seed nc")
  parser.add_argument("--coupling_period", help = "climate-ice coupling period in seconds")
  parser.add_argument("--ice_tstep_multiply", help = "is ice accelerated?")
  parser.add_argument("--region",        help = "GrIS or AIS?")
  args = parser.parse_args()

  err = 0
  hdf5_input,err  = arg_to_file_exist(args.hdf5_input,err=err)
  nc_output,err   = arg_to_file_exist(args.nc_output,io="out",err=err)
  hdf5_output,err = arg_to_file_exist(args.hdf5_output,io="out",err=err)

  if args.coupling_period != None:
    coupling_period = float(args.coupling_period)
  else:
    err = 4
    print("ERROR: specifiy a climate-ice coupling period")

  if args.ice_tstep_multiply != None:
    ice_tstep_multiply = float(args.ice_tstep_multiply)
  else:
    ice_tstep_multiply = 1.


  region="GrIS"
  if args.region != None:
    region = args.region

  if err > 0:
    parser.print_help()
    sys.exit(2)

  return hdf5_input, hdf5_output, nc_output, coupling_period, ice_tstep_multiply, region 


def extract_bisicles_topography(hdf5_input, region="GrIS"):
  from amrfile import io as amrio
  #from pyproj import Proj, transform #deprecated proj1 syntax now commented out and replaced
  from pyproj import Transformer

  level = 0 # level of grid refinement
  order = 0 # interpolation order, 0 for piecewise constant, 1 for linear

  #load BISICLES hdf5
  amrID  =  amrio.load(hdf5_input)
  lo,hi = amrio.queryDomainCorners(amrID, level)
  nx = hi[0] - lo[0] + 1
  ny = hi[1] - lo[1] + 1

  #extract cartesian grid and variables. x and y are the box centres distance from origin 
  x_1d,y_1d,z_surf = amrio.readBox2D(amrID, level, lo, hi, "Z_surface", order)
  x_1d,y_1d,calving = amrio.readBox2D(amrID, level, lo, hi, "melangeThickness", order)
  x_1d,y_1d,active_mask = amrio.readBox2D(amrID, level, lo, hi, "stableSourcesMask", order)

  #APPLY ACTIVE MASK
  inactive = np.where(active_mask < 0.5)
  calving[inactive] = 0.

  #MAKE 2D MESHES OF X,Y LIKE WITH THE OLD MISOMIP-TYP INITIAL SCRIPT!!!!
  dx = x_1d[1] - x_1d[0]
  dy = y_1d[1] - y_1d[0]
  y_2d,x_2d = np.mgrid[y_1d[0]:y_1d[-1] + dy:dy, x_1d[0]:x_1d[-1] + dx:dx]

  #corners not needed yet
  #y_2db1,x_2db1 = np.mgrid[y_1d[0]-dy/2:y_1d[-1]+dy-dy/2:dy, x_1d[0]-dx/2:x_1d[-1]+dx-dx/2:dx]
  #y_2db2,x_2db2 = np.mgrid[y_1d[0]-dy/2:y_1d[-1]+dy-dy/2:dy, x_1d[0]+dx/2:x_1d[-1]+dx+dx/2:dx]
  #y_2db3,x_2db3 = np.mgrid[y_1d[0]+dy/2:y_1d[-1]+dy+dy/2:dy, x_1d[0]+dx/2:x_1d[-1]+dx+dx/2:dx]
  #y_2db4,x_2db4 = np.mgrid[y_1d[0]+dy/2:y_1d[-1]+dy+dy/2:dy, x_1d[0]-dx/2:x_1d[-1]+dx-dx/2:dx]

  #outProj = Proj(init = 'epsg:4326')
  outProj = 'epsg:4326'

  print("RSS region is ",region)

  if region == "GrIS":
    #inProj = Proj(init = 'epsg:3413') ## ESPG number for Polar Sterographic North (70 degN, 45 degW) projection.
    inProj = 'epsg:3413'
    x0 = -654650.0 # sez Steph/amrtocf
    y0 = -3385950.0
  elif region == "AIS":
    #inProj = Proj(init='epsg:3031') #(WGS 84 / Antarctic Polar Stereographic)
    inProj = 'epsg:3031'
    x0=-3072000 # sez Steph
    y0 =-3072000
  elif region == "ASE":
    #inProj = Proj(init='epsg:3031') #(WGS 84 / Antarctic Polar Stereographic)
    inProj = 'epsg:3031'
    x0=-1831000 # sez Tony
    y0 =-904000
  else:
    print("extract_bisicles_topography: region unknown ",region)

  #Transform from the BISICLES grid info to lons and lats
  #Steph says he applied a x-mirror to make it match up with bedmap2! Needed? (Antarctica)
  mfact = 1.
  #if region == "AIS": mfact = -1

  #xl_2d,yl_2d = transform(inProj,outProj,(x_2d + x0)*mfact,(y_2d + y0))
  t = Transformer.from_crs(inProj,outProj,always_xy=True)
  xl_2d,yl_2d = t.transform((x_2d + x0)*mfact,(y_2d + y0))

  #corners not needed yet
  #xl_2db1,yl_2db1 = transform(inProj,outProj,(x_2db1+x0)*mfact,(y_2db1+y0))
  #xl_2db2,yl_2db2 = transform(inProj,outProj,(x_2db2+x0)*mfact,(y_2db2+y0))
  #xl_2db3,yl_2db3 = transform(inProj,outProj,(x_2db3+x0)*mfact,(y_2db3+y0))
  #xl_2db4,yl_2db4 = transform(inProj,outProj,(x_2db4+x0)*mfact,(y_2db4+y0))

  #Ferret, at least, says it's all off by ~90 lon as well as reflected in x! Or not? (Antarctica)
  steph_offset = 0.
  #if region == "AIS": steph_offset=90.

  xl_2d = xl_2d + steph_offset
  #corners not needed yet
  #xl_2db1 = xl_2db1+steph_offset
  #xl_2db2 = xl_2db2+steph_offset
  #xl_2db3 = xl_2db3+steph_offset
  #xl_2db4 = xl_2db4+steph_offset

  return z_surf,calving,active_mask,xl_2d,yl_2d,x_1d,y_1d

def construct_CF(field,lon,lat,x_1d,y_1d):


  #construct a minimal CF field to let cf-python do regridding
  field_CF = cf.Field()
  dimx = cf.DimensionCoordinate(data = cf.Data(x_1d, 'm'), properties = {'axis': 'X'})
  #dimx = cf.DimensionCoordinate(data = cf.Data(x_1d, 'm'))
  field_CF.set_construct(cf.DomainAxis(size=len(x_1d)),key="X")
  field_CF.set_construct(dimx,axes="X")
  dimy = cf.DimensionCoordinate(data = cf.Data(y_1d, 'm'), properties = {'axis': 'Y'})
  #dimy = cf.DimensionCoordinate(data = cf.Data(y_1d, 'm'))
  field_CF.set_construct(cf.DomainAxis(size=len(y_1d)),key="Y")
  field_CF.set_construct(dimy,axes="Y")
  lats = cf.AuxiliaryCoordinate(data = cf.Data(lat, 'degrees_north'), properties = {'standard_name': 'latitude'})
  field_CF.set_construct(lats, axes = ('Y', 'X'))
  lons = cf.AuxiliaryCoordinate(data = cf.Data(lon, 'degrees_east') , properties = {'standard_name': 'longitude'})
  field_CF.set_construct(lons, axes = ('Y', 'X'))

  field_ma=np.ma.masked_equal(field,0.)
  field_CF.set_data(cf.Data(field_ma,units = 'm'), axes = ('Y', 'X'))
  field_CF.set_property('_FillValue',-99.)


  return field_CF


def make_BISICLES_melange_sink(bike,x,y, coupling_interval_secs, calv_hdf):

  from netCDF4 import Dataset

  #turn the melange mass (metres of ice) we've just taken into a flux (m/s) to EXTRACT
  #from BISICLES over the next coupling period. 
  melange_sink = np.copy(bike)*-1 / coupling_interval_secs

  inactive = np.where(melange_sink > 0.)
  melange_sink[inactive] = 0.

  #turn the data into a BISICLES hdf5 via netcdf (can amrio do this directly, better?)
  ncfile_out = Dataset("calv_temp.nc",'w',format = 'NETCDF3_CLASSIC')
  ncfile_out.createDimension('x',x.shape[0])
  ncfile_out.createDimension('y',y.shape[0])

  x_nc = ncfile_out.createVariable('x',np.dtype('float64').char,('x'))
  y_nc = ncfile_out.createVariable('y',np.dtype('float64').char,('y'))
  data_nc = ncfile_out.createVariable('calved_melange',np.dtype('float64').char,('y','x'))

  x_nc[:] = x
  y_nc[:] = y
  data_nc[:] = melange_sink
  ncfile_out.close()

  cmd='./nctoamr2d.ex calv_temp.nc ' +  calv_hdf + ' calved_melange'
  os.system(cmd)
  cmd='rm calv_temp.nc'
  os.system(cmd)

  return

def regrid_one_stage_to_NEMO(field_CF,nemo_ncgridfile):

  #NEMO3.6 actually wants a volume of ice at a point - not a cell average (?)
  #NEMO4.2 has changed the input units back to mass
  #For volume multiply by the cell areas here, on the nice regular grid, and do a
  #simple nearest-neighbour point map across

  #python3 cf returns lon/lat with the old syntax of asking for X and Y, we want to work
  #out cell lengths here
  #key_dy=list(field_CF.dimension_coordinates().keys())[0]
  #key_dx=list(field_CF.dimension_coordinates().keys())[1]
  #dx = field_CF.dimension_coordinate(key_dx).array[1] - field_CF.dimension_coordinate(key_dx).array[0]
  #dy = field_CF.dimension_coordinate(key_dy).array[1] - field_CF.dimension_coordinate(key_dy).array[0]
  #field_CF = field_CF * dx * dy

  #python3 cf masking on field_CF was screwing up the regrid
  f1=np.asarray(field_CF.array)
  field_CF.set_data(cf.Data(f1))

  g = cf.read(nemo_ncgridfile)[0]
  #field CF is not identifying axes in a way regrids can identify. Both axis and aux_coord are
  #identifying as 'X'? Add src_axes=['latitude','longitude']?
  field_REGRID = field_CF.regrids(g, dst_cyclic = True, method = 'nearest_dtos', src_axes={'X': 'key%X', 'Y': 'key%Y'})


  return field_REGRID

def reroute_calving_to_ocean(calv_cf,mask_cf,route_file,coupling_period,ice_tstep_multiply,orig_calving,output_file):

  # extract the numpy array, with 0s instead of missing_data
  calv_ma = calv_cf.array
  np.ma.set_fill_value(calv_ma,0)
  calv = np.ma.filled(calv_ma)

  mask = mask_cf.array

  #NEMO3.6 this is m3 per coupling interval. We want km^3/yr
  #what are NEMO4.2 units?
  #Assume 360_day calendar
  calv = calv / 1e9 / coupling_period * 360. * 24. * 60. *60. 

  #if ice has been accelerated, pull back the iceberg magnitude given to NEMO
  calv = calv / ice_tstep_multiply

  #set up the routing map
  mapping = cf.read(route_file)
  try:
    iindex = mapping.select_field("ncvar%xindex").array
  except:
    iindex = mapping[1].array
  try:
    jindex = mapping.select_field("ncvar%yindex").array
  except:
    jindex = mapping[0].array

  nj = np.shape(iindex)[0]
  ni = np.shape(iindex)[1]

  #loop over the routed area and move everything to its destination
  calv_routed = np.zeros([nj,ni])
  for j in np.arange(nj):
    for i in np.arange(ni):

      if mask[j,i] > 0:

        iiindex = int(iindex[j,i])
        ijindex = int(jindex[j,i])

        if (ijindex == 0) & (iiindex == 0):  #if routfile is not appropriate
          calv_routed[j,i] = calv[j,i]
        else: 
          calv_routed[ijindex,iiindex] = calv_routed[ijindex,iiindex] + calv[j,i]

  #this is just the regional field - will construct a unified global field elsewhere
  h = cf.read(orig_calving)

  original_seed = h.select_field("ncvar%calvingmask")

  new_seed=np.copy(original_seed.array)
  for i in range(np.shape(new_seed)[0]):
    new_values = np.where(calv_routed > 0.)
    new_seed[i,:,:][new_values]=calv_routed[new_values]
    #new_seed[i,:,:]=calv[:,:]
  
  original_seed.set_data(cf.Data(new_seed))

  cf.write(h,output_file,fmt = "NETCDF4")

  return

if __name__ == "__main__":

  hdf5in_file, hdf5out_file, ncout_file, coupling_period, ice_tstep_multiply, region  = parse_commandline(sys.argv)
  
  #get what we need from the BISICLES plot file
  print("extract fields from BISICLES hdf5")
  surface_height,calving,active_mask,lon,lat,x_dist,y_dist = extract_bisicles_topography(hdf5in_file,region=region)
  
  #make calving source/sink files for NEMO/BISICLES. These grid/template files are linked to 
  #appropriate sources by the suite
  print("deal with the iceberg/melange files")

  nemo_ncgridfile = "cf_gridfile.nc"
  route_file = "cf_routfile.nc"
  orig_calving = "orig_calving.nc"
  
  make_BISICLES_melange_sink(calving,x_dist,y_dist,coupling_period, hdf5out_file)
  
  calving_CF = construct_CF(calving,lon,lat,x_dist,y_dist)
  calving_REGRID = regrid_one_stage_to_NEMO(calving_CF,nemo_ncgridfile)
  active_mask_CF = construct_CF(active_mask,lon,lat,x_dist,y_dist)
  active_mask_REGRID = regrid_one_stage_to_NEMO(active_mask_CF,nemo_ncgridfile)
  reroute_calving_to_ocean(calving_REGRID,active_mask_REGRID,route_file,coupling_period,ice_tstep_multiply,orig_calving,ncout_file)
