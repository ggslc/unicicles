import numpy as np
import sys
import os
from arg_to_file_exist import arg_to_file_exist

def parse_commandline(argv):
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument("--um_input",      help = "input, name of UM ancil to be modified")
  parser.add_argument("--hdf5_input",    help = "input, name of unicicles hdf5")
  parser.add_argument("--hdf5_output",   help = "output, name of unicicles melange sink hdf5")
  parser.add_argument("--nc_output",     help = "output, name of nemo iceberg seed nc")
  parser.add_argument("--um_output",     help = "output, name of modified UM ancil")
  parser.add_argument("--coupling_period", help = "climate-ice coupling period in seconds")
  parser.add_argument("--ice_tstep_multiply", help = "is ice accelerated?")
  parser.add_argument("--region",        help = "GrIS or AIS?")
  args = parser.parse_args()

  err = 0
  um_input,err    = arg_to_file_exist(args.um_input,err=err)
  hdf5_input,err  = arg_to_file_exist(args.hdf5_input,err=err)
  um_output,err   = arg_to_file_exist(args.um_output, io="out",err=err)
  nc_output,err   = arg_to_file_exist(args.nc_output,io="out",err=err)
  hdf5_output,err = arg_to_file_exist(args.hdf5_output,io="out",err=err)

  if args.coupling_period != None:
    coupling_period = np.float(args.coupling_period)
  else:
    err = 4
    print "ERROR: specifiy a climate-ice coupling period"

  if args.ice_tstep_multiply != None:
    ice_tstep_multiply = np.float(args.ice_tstep_multiply)
  else:
    ice_tstep_multiply = 1.


  region="GrIS"
  if args.region != None:
    region = args.region

  if err > 0:
    parser.print_help()
    sys.exit(2)

  return um_input, hdf5_input, um_output, hdf5_output, nc_output, coupling_period, ice_tstep_multiply, region 


def extract_bisicles_topography(hdf5_input, region="GrIS"):
  from amrfile import io as amrio
  from pyproj import Proj, transform

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
  #x_1d,y_1d,calving = amrio.readBox2D(amrID, level, lo, hi, "accumCalvedIceThickness", order)

  #MAKE 2D MESHES OF X,Y LIKE WITH THE OLD MISOMIP-TYP INITIAL SCRIPT!!!!
  dx = x_1d[1] - x_1d[0]
  dy = y_1d[1] - y_1d[0]
  y_2d,x_2d = np.mgrid[y_1d[0]:y_1d[-1] + dy:dy, x_1d[0]:x_1d[-1] + dx:dx]

  #corners not needed yet
  #y_2db1,x_2db1 = np.mgrid[y_1d[0]-dy/2:y_1d[-1]+dy-dy/2:dy, x_1d[0]-dx/2:x_1d[-1]+dx-dx/2:dx]
  #y_2db2,x_2db2 = np.mgrid[y_1d[0]-dy/2:y_1d[-1]+dy-dy/2:dy, x_1d[0]+dx/2:x_1d[-1]+dx+dx/2:dx]
  #y_2db3,x_2db3 = np.mgrid[y_1d[0]+dy/2:y_1d[-1]+dy+dy/2:dy, x_1d[0]+dx/2:x_1d[-1]+dx+dx/2:dx]
  #y_2db4,x_2db4 = np.mgrid[y_1d[0]+dy/2:y_1d[-1]+dy+dy/2:dy, x_1d[0]-dx/2:x_1d[-1]+dx-dx/2:dx]

  outProj = Proj(init = 'epsg:4326')

  if region == "GrIS":
    inProj = Proj(init = 'epsg:3413') ## ESPG number for Polar Sterographic North (70 degN, 45 degW) projection.
    x0 = -654650.0 # sez Steph/amrtocf
    y0 = -3385950.0
  elif region == "AIS":
    inProj = Proj(init='epsg:3031') #(WGS 84 / Antarctic Polar Stereographic)
    x0=-3072000 # sez Steph
    y0 =-3072000
  else:
    print "extract_bisicles_topography: region unknown ",region

  #Transform from the BISICLES grid info to lons and lats
  #Steph says he applied a x-mirror to make it match up with bedmap2! Needed? (Antarctica)
  mfact = 1.
  #if region == "AIS": mfact = -1

  xl_2d,yl_2d = transform(inProj,outProj,(x_2d + x0)*mfact,(y_2d + y0))

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

  return z_surf,calving,xl_2d,yl_2d,x_1d,y_1d

def construct_CF(field,lon,lat,x_1d,y_1d):

  import cf

  #construct a minimal CF field to let cf-python do regridding
  field_CF = cf.Field()
  dimx = cf.DimensionCoordinate(data = cf.Data(x_1d, 'm'), properties = {'axis': 'X'})
  field_CF.insert_dim(dimx)
  dimy = cf.DimensionCoordinate(data = cf.Data(y_1d, 'm'), properties = {'axis': 'Y'})
  field_CF.insert_dim(dimy)
  lats = cf.AuxiliaryCoordinate(data = cf.Data(lat, 'degrees_north'), properties = {'standard_name': 'latitude'})
  field_CF.insert_aux(lats, axes = ('Y', 'X'))
  lons = cf.AuxiliaryCoordinate(data = cf.Data(lon, 'degrees_east') , properties = {'standard_name': 'longitude'})
  field_CF.insert_aux(lons, axes = ('Y', 'X'))

  field_ma=np.ma.masked_equal(field,0.)
  field_CF.insert_data(cf.Data(field_ma,units = 'm'), axes = ('Y', 'X'))
  field_CF.setprop('_FillValue',-99.)

  return field_CF

def regrid_two_stage_to_GLOBE30(field_CF,region="GrIS"):
  
  import cf

  #cf-python ESMF regridding is very quick for the basic transform onto a regular lat-lon grid.
  #But memory limitations mean it can't handle going to the 30min GLOBE grid we want to be able 
  #to feed into CAP. We do a coarser regrid in cf-python, then call cdo (which is slow for the
  #first regrid) to refine it. Honestly found to be the only solution with acceptable run times
  #in testing...

  #construct the coarse regular grid in CF
  if region == "GrIS":
    ymin = 45.
    ymax = 90.
    xmin = 250.
    xmax = 360.
    xcyclic = False
  elif region == "AIS":
    ymin = -90.
    ymax = -55.
    xmin = 0.
    xmax = 360.
    xcyclic = True

  dy = 0.083333
  dx = 0.083333

  x_1d = np.arange(xmin + dx/2.,xmax - dx/2.,dx)
  y_1d = np.arange(ymin + dx/2.,ymax - dx/2.,dy)

  dimx = cf.DimensionCoordinate(data = cf.Data(x_1d, 'degrees_east'), properties = {'axis': 'X','standard_name': 'longitude'})
  dimy = cf.DimensionCoordinate(data = cf.Data(y_1d, 'degrees_north'), properties = {'axis': 'Y','standard_name': 'latitude'})
  g = cf.Field()
  g.insert_dim(dimx)
  g.insert_dim(dimy)
  g.insert_data(cf.Data(np.zeros([dimy.shape[0],dimx.shape[0]])), axes = ('Y', 'X'))

  #regrid from BISICLES to the coarse, regular grid
  field_REGRID_1 = field_CF.regrids(g,dst_cyclic = xcyclic, method = 'bilinear')
  #cf.write(field_CF,'cf_polarstereofield.nc', fmt = "NETCDF3_CLASSIC")
  cf.write(field_REGRID_1,'cf_0.083deglatlonfield.nc', fmt = "NETCDF3_CLASSIC")

  #write out cdo grid definitions for coarse and fine grids
  gf = open('gridfile-Reg0.083','w')
  gf.write('#0.083 lon lat grid file \n')
  gf.write('gridtype  =  lonlat\n')
  gf.write('xsize ='+str(len(x_1d))+'\n' )
  gf.write('ysize ='+str(len(y_1d))+'\n' )
  gf.write('xfirst ='+str(xmin+dx/2.)+'\n' )
  gf.write('xinc ='+str(dx)+'\n' )
  gf.write('yfirst ='+str(ymin+dy/2.)+'\n' )
  gf.write('yinc ='+str(dy)+'\n' )
  gf.close()

  resfact = 10.
  gf = open('gridfile-Reg0.0083','w')
  gf.write('#0.0083 res lon lat grid file \n')
  gf.write('gridtype = lonlat\n')
  gf.write('xsize ='+str(resfact*len(x_1d))+'\n' )
  gf.write('ysize ='+str(resfact*len(y_1d))+'\n' )
  gf.write('xfirst ='+str(xmin+dx/(2.*resfact))+'\n' )
  gf.write('xinc ='+str(dx/10)+'\n' )
  gf.write('yfirst ='+str(ymin+dy/(2.*resfact))+'\n' )
  gf.write('yinc ='+str(dy/resfact)+'\n' )
  gf.close()

  #get cdo to refine by a factor of 10
  cmd="cdo remapbil,gridfile-Reg0.0083 -setgrid,gridfile-Reg0.083 cf_0.083deglatlonfield.nc cdo_0.0083deglatlonfield.nc"
  os.system(cmd)


  field_REGRID_2 = cf.read("cdo_0.0083deglatlonfield.nc")[0].array
  cmd = 'rm cdo_0.0083deglatlonfield.nc gridfile-Reg0.0083 gridfile-Reg0.083 cf_0.083deglatlonfield.nc'
  os.system(cmd)

  return field_REGRID_2

def splice_into_ancil_template(topog_data,ancil_template,ancil_output):
  import mule
  import mule_rss

  fieldg = mule.load_umfile(ancil_template)
  globe30 = fieldg.fields[0].get_data()

#flip N-S
  topog_data = topog_data[::-1,:]

#update where we have data
  update = np.where(topog_data.mask == False)
  globe30[update] = topog_data[update]

  fieldg.fields[0] = mule_rss.overwrite_data(fieldg.fields[0],globe30)
  fieldg.to_file(ancil_output)

  return

def make_BISICLES_melange_sink(bike,x,y, coupling_interval_secs, calv_hdf):

  from netCDF4 import Dataset

  #turn the melange mass (metres of ice) we've just taken into a flux (m/s) to EXTRACT
  #from BISICLES over the next coupling period. 
  melange_sink = np.copy(bike)*-1 / coupling_interval_secs

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

  import cf

  #NEMO actually wants a volume of ice at a point - not a cell average (?)
  #Multiply by the cell areas here, on the nice regular grid, and do a
  #simple nearest-neighbour point map across
  dx = field_CF.item('X').array[1] - field_CF.item('X').array[0]
  dy = field_CF.item('Y').array[1] - field_CF.item('Y').array[0]
  field_CF = field_CF * dx * dy

  g = cf.read_field(nemo_ncgridfile)
  field_REGRID = field_CF.regrids(g,dst_cyclic = True,method = 'nearest_dtos')

  return field_REGRID

def reroute_calving_to_ocean(calv_cf,route_file,coupling_period,ice_tstep_multiply,orig_calving,output_file):

  import cf

  # extract the numpy array, with 0s instead of missing_data
  calv_ma = calv_cf.array
  np.ma.set_fill_value(calv_ma,0)
  calv = np.ma.filled(calv_ma)

  #this is m3 per coupling interval. We want km^3/yr for NEMO.
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
      if (jindex[j,i] == 0) & (iindex[j,i] == 0):  #if routfile is not appropriate
        calv_routed[j,i] = calv[j,i]
      else: 
        calv_routed[jindex[j,i],iindex[j,i]] = calv_routed[jindex[j,i],iindex[j,i]] + calv[j,i]

  #this is just the regional field - will construct a unified global field elsewhere
  h = cf.read(orig_calving)

  original_seed = h.select_field("ncvar%calvingmask")

  new_seed=np.copy(original_seed.array)
  for i in range(np.shape(new_seed)[0]):
    new_seed[i,:,:]=calv_routed[:,:]
  
  original_seed.insert_data(cf.Data(new_seed))

  cf.write(h,output_file,fmt = "NETCDF4")

  return

if __name__ == "__main__":

  ancilin_file, hdf5in_file, ancilout_file, hdf5out_file, ncout_file, coupling_period, ice_tstep_multiply, region  = parse_commandline(sys.argv)
  
  #get what we need from the BISICLES plot file
  print "extract fields from BISICLES hdf5"
  surface_height,calving,lon,lat,x_dist,y_dist = extract_bisicles_topography(hdf5in_file,region=region)
  
  #make topography ancil for CAP
  print "regrid, splice CAP fields"
  surface_height_CF = construct_CF(surface_height,lon,lat,x_dist,y_dist)
  surface_height_REGRID = regrid_two_stage_to_GLOBE30(surface_height_CF,region=region)
  splice_into_ancil_template(surface_height_REGRID,ancilin_file,ancilout_file)
  
  #make calving source/sink files for NEMO/BISICLES. These grid/template files are linked to 
  #appropriate sources by the suite

  nemo_ncgridfile = "cf_gridfile.nc"
  route_file = "cf_routfile.nc"
  orig_calving = "orig_calving.nc"
  
  make_BISICLES_melange_sink(calving,x_dist,y_dist,coupling_period, hdf5out_file)
  
  calving_CF = construct_CF(calving,lon,lat,x_dist,y_dist)
  calving_REGRID = regrid_one_stage_to_NEMO(calving_CF,nemo_ncgridfile)
  reroute_calving_to_ocean(calving_REGRID,route_file,coupling_period,ice_tstep_multiply,orig_calving,ncout_file)
