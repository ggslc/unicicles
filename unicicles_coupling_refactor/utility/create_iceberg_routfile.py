import numpy as np
from netCDF4 import Dataset
import sys
#from ferretlook import *

def arg_to_file_exist(arg, mandatory=True, io="in",err=0):
  import os

  filename = None
  if arg == None:
    if mandatory:
      print "ERROR: mandatory argument missing"
      err = 1
  else:
      filename=arg
      exists = os.path.isfile(filename)
      if (exists) & (io == "out"):
        print "ERROR: output file already exists ",filename
        err = 2
      if (not exists) & (io == "in"):
        if mandatory:
          print "ERROR: input file does not exist ",filename
          err = 3
        else:
          print "WARNING: input file does not exist ",filename
          filename = None

  return filename, err

def parse_commandline(argv):
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument("--input",  help="input, name of nemo meshmask file")
  parser.add_argument("--output", help="output, name of iceberg routing file")
  args = parser.parse_args()

  err = 0
  input,err  = arg_to_file_exist(args.input,err=err)
  output,err = arg_to_file_exist(args.output,err=err,io="out")

  if err > 0:
    parser.print_help()
    sys.exit(2)

  return input, output

if __name__ == "__main__":

  bathyfile, routfile = parse_commandline(sys.argv)

  try:
      tmask = Dataset(bathyfile).variables['tmask'][:].squeeze()[0,...]
  except:
      print "ERROR: tmask variable not found in input file"
      exit()

  ny = np.shape(tmask)[0]
  nx = np.shape(tmask)[1]

  res = None
  if ny == 332:  res = "eORCA1"
  if ny == 1207: res = "eORCA025"
  
  #bounds of boxes within which we'll say BISICLES might produce
  #land-locked icebergs - only find routes for land inside these areas
  #below s_lim is Antarctica
  #above n_lim and between e_lim and w_lim is ~Greenland
  if res == "eORCA025":
    s_lim = 400
    n_lim = 950
    e_lim = 940
    w_lim = 1100
  elif res == "eORCA1":
    s_lim = 100
    n_lim = 270
    e_lim = 235
    w_lim = 275
  else:
    print "ERROR: resolution of input meshmask not recognised - can only do eORCA1 or eORCA025"
    exit()
  
  #indices of land/shelf points
  land_y,land_x = np.where(tmask < 1)
  
  #rough find of indices of icesheet points
  ais  = np.where(land_y < s_lim)
  gris = np.where( (land_y > n_lim) & (land_x > e_lim) & (land_x < w_lim) )
  
  #mask of useful things: 0 is ocean, 1 is land outside our regions, 99 is ice to route from
  mask = np.zeros_like(tmask)
  mask[land_y,land_x] = 9
  mask[land_y[ais],land_x[ais]] = 99
  mask[land_y[gris],land_x[gris]] = 99

  #ferretlook(mask)
  
  #indices of open ocean targets and potential ice source points
  open_y,open_x = np.where(mask < 1)
  ice_y,ice_x   = np.where(mask > 10)
  
  #set all the non-ice points to route to themselves
  x_index = np.zeros_like(tmask, dtype=np.int16)
  y_index = np.zeros_like(tmask, dtype=np.int16)
  x_index[open_y,open_x] = open_x
  y_index[open_y,open_x] = open_y
  x_index[land_y,land_x] = land_x
  y_index[land_y,land_x] = land_y
  
  #for each ice point find the nearest open ocean point
  print " "
  print "finding nearest wet cells to ice"
  for i in range(len(ice_y)): 
    print '\r', i,":",len(ice_y),'         ',
    dx = np.abs(ice_x[i] - open_x)
    dy = np.abs(ice_y[i] - open_y)
    dist = dx*dx + dy*dy
    min_pos = np.argmin(dist)
    x_min = open_x[min_pos]
    y_min = open_y[min_pos]
  
    x_index[ice_y[i],ice_x[i]] = x_min
    y_index[ice_y[i],ice_x[i]] = y_min
  
  #write an output file     
  outfile = Dataset(routfile,'w')
  outfile.createDimension('x',nx)
  outfile.createDimension('y',ny)
  
  outfield = outfile.createVariable('lon',np.float,('y','x'))
  outfield.setncattr('units','degrees_east')
  outfield.setncattr('long_name','longitude')
  outfield[:] = Dataset(bathyfile).variables['nav_lon'][:].squeeze()
  
  outfield = outfile.createVariable('lat',np.float,('y','x'))
  outfield.setncattr('units','degrees_north')
  outfield.setncattr('long_name','latitude')
  outfield[:] = Dataset(bathyfile).variables['nav_lat'][:].squeeze()
  
  outfield = outfile.createVariable('xindex',np.float,('y','x'))
  outfield.setncattr('long_name','x index of nearest wet surface cell to ice')
  outfield[:] = x_index
  
  outfield = outfile.createVariable('yindex',np.float,('y','x'))
  outfield.setncattr('long_name','y index of nearest wet surface cell to ice')
  outfield[:] = y_index
  
  outfile.close()
