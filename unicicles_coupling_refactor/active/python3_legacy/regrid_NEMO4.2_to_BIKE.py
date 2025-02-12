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

#melt_water=h.select('long_name=Ice shelf melting')[-1].array.squeeze()
#melt_heat =h.select('long_name=Ice shelf heat content flux')[-1].array.squeeze()
melt_water=h.select_field('long_name=Ice shelf fresh water flux ( from isf to oce )').array.squeeze()
#melt_heat =h.select_field('long_name=Ice shelf heat content flux of injected water ( from isf to oce )').array.squeeze()
melt_heat =h.select_field('long_name=Ice shelf ocean  heat flux ( from isf to oce )').array.squeeze()


#if you mask out grounded areas, cf python interpolates missing_data blanks into the cavity 
#quite a way. Better to use 0s and over estimate the useful area we can cover
##melt_water=np.ma.masked_equal(melt_water,0)
##melt_heat=np.ma.masked_equal(melt_heat,0)

#need to do better than this - diagnostics on ORCA grids may have junk in the halos/overlap points, it seems

#time average
if melt_water.ndim==3:
  print("Doing time average",np.shape(melt_water))
  print("I shouldn't have to do this - have I picked up a 3d volume field by accident?!")
  exit(1)
  melt_water=np.mean(melt_water,axis=0)
  melt_heat =np.mean(melt_heat,axis=0)
  print(np.shape(melt_water))

g=cf.read(bike_ncgridfile)[0]

print(" regrid 1")
f=cf.read(nemo_ncgridfile)[0]
f.set_data(cf.Data(melt_water,units='kg/m2/s'))
w_regrid=f.regrids(g,src_cyclic=True,method='linear')
#w_regrid=g

print(" regrid 2")
f=cf.read(nemo_ncgridfile)[0]
f.set_data(cf.Data(melt_heat,units='kg/m2/s'), axes=('Y', 'X'))
h_regrid=f.regrids(g,src_cyclic=True,method='linear')
#h_regrid=g


#i=cf.FieldList()
#i.append(w_regrid)
#i.append(h_regrid)
#cf.write(i,regrid_ncfile)

from netCDF4 import Dataset
#like BIKEtoNEMO needs a rollaxis (and I don't get why)
#we need to switch x and y here to get it the way BISICLES wants it
#heat=np.swapaxes(h_regrid.array,0,1)
#melt=np.swapaxes(w_regrid.array,0,1)
heat=h_regrid.array
melt=w_regrid.array

heat.fill_value=0.
heat=heat.filled()*heat_unit_factor
melt.fill_value=0.
melt=melt.filled()*water_unit_factor

x_bike=np.load("bike_xcoords.dump",allow_pickle=True,fix_imports=True,encoding='latin1')
y_bike=np.load("bike_ycoords.dump",allow_pickle=True,fix_imports=True,encoding='latin1')
#touch commit test #2

#--------the beginning of the sub-script to fill GL meltrate and the conservation of total melt flux (a.siahaan)
#--------I think this is just conservation, I don't think it's working reliably
#--------is it the calculation, or that this map file has some points mapping to multiple NEMO boxes?
#bisicles_to_nemo_mapping_file="bikegridmap_file.map2d"
#bn_map=mp.load(bisicles_to_nemo_mapping_file)
#ny_n=100
#nx_n=np.shape(bn_map.x)[1]
#
#for j in range(ny_n):
#  for i in range(nx_n):
#     if (np.abs(melt_water[j,i]) > 0):
#        ncontrib = np.int(bn_map.nmap[j,i])
#        meltsub = np.zeros(ncontrib)
#        heatsub = np.zeros(ncontrib)
#        for jcont in range(ncontrib):
#            ji, jj = bn_map.x[j,i,jcont],bn_map.y[j,i,jcont]
#            meltsub[jcont] = melt[ji,jj]
#            heatsub[jcont] = heat[ji,jj]
#        avgmeltsub = np.sum(meltsub)/np.sum(np.abs(meltsub)>1e-10)
#        avgheatsub = np.sum(heatsub)/np.sum(np.abs(heatsub)>1e-10)
#        for jcont in range(ncontrib):
#            ji, jj = bn_map.x[j,i,jcont],bn_map.y[j,i,jcont]
#            melt[ji,jj] = melt[ji,jj] * melt_water[j,i] * water_unit_factor/avgmeltsub
#            heat[ji,jj] = heat[ji,jj] * melt_heat[j,i] * heat_unit_factor/avgheatsub

#--------the end of the sub-script to fill GL meltrate and the conservation of total melt flux------


#add Steph's calving-edge-enforced-by-melt criterion - how now?
#melt[np.where()]=-1.0e+4

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
