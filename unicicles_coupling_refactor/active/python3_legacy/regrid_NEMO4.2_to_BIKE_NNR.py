import os
import numpy as np
import sys
import cf
#import mapping_class as mp


def parse_commandline(argv):
  import os.path
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument("--input",      help="input, name of nemo input nc file")
  parser.add_argument("--output",     help="input, name of bisicles output hdf file")
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

  if args.output != None:
    output=args.output
  else:
    print("")
    print("ERROR: explicitly specifiy an output file")
    print("")
    parser.print_help()
    sys.exit(2)

  return input,output


nctgridfile, regrid_hdf5file=parse_commandline(sys.argv)

regrid_ncfile  =os.path.splitext(os.path.basename(regrid_hdf5file))[0]+".nc"

bike_ncgridfile="cf_bikegridfile.nc"
nemo_ncgridfile="cf_gridfile.nc"

#NEMO melt rate is positive for melting and negative for freezing, BISICLES has the opposite sign convention, so multiply NEMO melt by -1 
water_unit_factor=-60*60*24*360/1e3
heat_unit_factor=60*60*24*360

h=cf.read(nctgridfile)
g=cf.read(bike_ncgridfile)[0]

melt_water=h.select_field('long_name=Ice shelf fresh water flux ( from isf to oce )').array.squeeze()
melt_water = np.ma.masked_equal(melt_water,0.)
melt_heat =h.select_field('long_name=Ice shelf ocean  heat flux ( from isf to oce )').array.squeeze()
melt_heat = np.ma.masked_equal(melt_heat,0.)

print(" regridding melt field with NNR info to GL!!!")
f = cf.read(nemo_ncgridfile)[0]
f.set_data(cf.Data(melt_water, units='kg/m2/s'), axes=('Y', 'X'))

w_regrid_nnr = f.regrids(g, src_cyclic=True, method='nearest_stod')
w_regrid_lin = f.regrids(g, src_cyclic=True, method='linear')
melt_lin = w_regrid_lin.array
melt_nnr = w_regrid_nnr.array

melt_lin.fill_value = 0.
melt = melt_lin.filled() * water_unit_factor
melt[np.where(melt==0.)]=melt_nnr[np.where(melt==0.)]*water_unit_factor
melt = np.where(melt > 1e20, 0., melt)
#mask = np.ma.getmaskarray(melt_lin)
#melt = np.where(mask, melt_nnr.data, melt_lin.data)
#melt = np.where(melt > 1e20, 0., melt)*water_unit_factor

print(" regridding heat field with NNR info to GL!!!")
f = cf.read(nemo_ncgridfile)[0]
f.set_data(cf.Data(melt_heat, units='J/m2/s'), axes=('Y', 'X'))
h_regrid_nnr = f.regrids(g, src_cyclic=True, method='nearest_stod')
h_regrid_lin = f.regrids(g, src_cyclic=True, method='linear')
heat_lin = h_regrid_lin.array
heat_nnr = h_regrid_nnr.array

heat_lin.fill_value = 0.
heat = heat_lin.filled() * heat_unit_factor
heat[np.where(heat==0.)]=heat_nnr[np.where(heat==0.)]*heat_unit_factor
heat = np.where(heat > 1e20, 0., heat)
#mask = np.ma.getmaskarray(heat_lin)
#heat = np.where(mask, heat_nnr.data, heat_lin.data)
#heat = np.where(heat > 1e20, 0., heat)*heat_unit_factor

from netCDF4 import Dataset
#like BIKEtoNEMO some domains needs a rollaxis (and I don't get why)
#we need to switch x and y here to get it the way BISICLES wants it
#heat=np.swapaxes(h_regrid.array,0,1)
#melt=np.swapaxes(w_regrid.array,0,1)
#heat=h_regrid.array
#melt=w_regrid.array

x_bike=np.load("bike_xcoords.dump",allow_pickle=True,fix_imports=True,encoding='latin1')
y_bike=np.load("bike_ycoords.dump",allow_pickle=True,fix_imports=True,encoding='latin1')

ncfile_out = Dataset(regrid_ncfile,'w',format='NETCDF3_CLASSIC')
ncfile_out.createDimension('x',x_bike.shape[1])
ncfile_out.createDimension('y',x_bike.shape[0])

x_nc=ncfile_out.createVariable('x',np.dtype('float64').char,('x'))
y_nc=ncfile_out.createVariable('y',np.dtype('float64').char,('y'))

x_nc[:]=x_bike[0,:]
y_nc[:]=y_bike[:,0]

heat_nc=ncfile_out.createVariable('melt_heat',np.dtype('float64').char,('y','x'))
melt_nc=ncfile_out.createVariable('melt_water',np.dtype('float64').char,('y','x'))

heat_nc[:]=heat
melt_nc[:]=melt

ncfile_out.close()

print(" calling nctoamr")
os.system('./nctoamr2d.ex %(regrid_ncfile)s %(regrid_hdf5file)s melt_water melt_heat' % locals())
