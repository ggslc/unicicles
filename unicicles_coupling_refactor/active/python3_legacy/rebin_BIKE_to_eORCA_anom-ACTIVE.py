import numpy as np
import sys
from amrfile import io as amrio
import mapping_class as mp
from arg_to_file_exist import arg_to_file_exist
from netCDF4 import Dataset


thresh=1e-3
mindep_isf=1e-3
mindep_bath=1e-3
mincav=1e-3

def parse_commandline(argv):
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument("--hdf5_input",  help="input, name of bisicles input file")
  parser.add_argument("--nc_output",   help="output, name of nemo output file")
  parser.add_argument("--nc_template", help="input, name of nemo original file")
  parser.add_argument("--verbose" )
  args = parser.parse_args()

  err = 0
  l_verb = False
  if args.verbose: l_verb = True

  hdf5_input, err  = arg_to_file_exist(args.hdf5_input, mandatory=True, err=err)
  nc_template, err = arg_to_file_exist(args.nc_template, mandatory=True, err=err)
  nc_output, err   = arg_to_file_exist(args.nc_output, mandatory=True, io="out", err=err)

  if err > 0:
    parser.print_help()
    sys.exit(2)

  return hdf5_input, nc_template, nc_output, l_verb

def enforce_minima(array,thresh,mindep):

  shallow=np.where( (array < mindep) & (array > thresh) )
  array[shallow]=mindep
  zero=np.where( (array < thresh))
  array[zero]=0.

  return array

def make_surface_mask(bathy,isf,mindep_isf,mindep_bath,thresh):

  mask=np.zeros_like(bathy)
  #mask[ (isf > mindep_isf-thresh) | (bathy < mindep_bath+thresh) ] = 1
  mask[ (isf > thresh) | (bathy < thresh) ] = 1

  return mask


def isf_from_bisicles_median(bathy_in, isf_in, n_mapping, x_mapping, y_mapping, bathy_thresh,l_verb=False):

  n_contrib = 0.
  n_cav = 0
  n_open = 0
  bathy_acc = 0.
  isf_acc = 0.

  for n in range(int(n_mapping)):

    bathy_in_point = bathy_in[int(y_mapping[n]),int(x_mapping[n])]
    isf_in_point   = isf_in[int(y_mapping[n]),int(x_mapping[n])]

    if l_verb: print(n,bathy_in_point,isf_in_point)

    if (bathy_in_point > bathy_thresh):
      n_contrib = n_contrib + 1.
      bathy_acc = bathy_acc + bathy_in_point
      isf_acc   = isf_acc   + isf_in_point

      if   (np.abs(isf_in_point - bathy_in_point)  > thresh) \
         & (np.abs(isf_in_point) > thresh): n_cav = n_cav + 1 #not grounded, not open ocean

      if np.abs(isf_in_point) < thresh: n_open = n_open + 1   #open ocean


  if l_verb: print(n_contrib, n_cav, n_open)
  if l_verb: print(isf_acc, bathy_acc)

  isf_out = -99999. #use for detection of out-of-domain points later
  bathy_out = -99999. 

  if n_contrib > 0: 

    bathy_out = bathy_acc / n_contrib

    if (n_open+n_cav) >= n_contrib/2.: #enough points are not grounded to do something
      #if n_open > n_cav: 
      if n_open == n_contrib:
          isf_out = 0.
      else:
          isf_out = np.abs((bathy_acc - isf_acc)) / n_contrib #enough points are ungrounded 
                                                              #and not open to have a cavity
    else:
      isf_out = -1. #grounded

    if l_verb: print(isf_acc/n_contrib, isf_out)

  return isf_out, bathy_out

hdfplotfile, nemo_templatefile, rebin_ncfile, l_verb = parse_commandline(sys.argv)

bisicles_to_nemo_mapping_file="bikegridmap_file.map2d"

#bisicles bathymetry that we're testing is -ve (below SL) and +ve (above SL)
bathy_threshold=-5000

level=0 # level of grid refinement
order=1 # interpolation order, 0 for piecewise constant, 1 for linear

if (l_verb): print('reading hdf file')
amrID = amrio.load(hdfplotfile)
lo,hi = amrio.queryDomainCorners(amrID, level)
_,_,isf = amrio.readBox2D(amrID, level, lo, hi, "Z_bottom", order)
_,_,bathy = amrio.readBox2D(amrID, level, lo, hi, "Z_base", order)
_,_,active_mask = amrio.readBox2D(amrID, level, lo, hi, "stableSourcesMask", order)

#APPLY ACTIVE MASK
inactive = np.where(active_mask < 0.5)
isf[inactive] = -9999.
bathy[inactive] = -9999.

print(bathy.shape, isf.shape)

##somewhere in the tangle reading the hdf5 this way gets me arrays transposed compared
##to the way they were read/used when generating the lats and lons 
#bathy=np.transpose(bathy)
#isf=np.transpose(isf)
##

if (l_verb): print('restoring mapping files')
bn_map=mp.load(bisicles_to_nemo_mapping_file)

ny_n=np.shape(bn_map.x)[0]
nx_n=np.shape(bn_map.x)[1]

isf_cav=np.zeros([ny_n,nx_n])
bathy_ism=np.zeros([ny_n,nx_n])

print(np.sum(isf))


template = Dataset(nemo_templatefile,'r')
bathy_tem_v = template.variables['Bathymetry_isf']
isf_tem_v = template.variables['isf_draft']
bathy_tem = bathy_tem_v[:]
isf_tem = isf_tem_v[:]

test_i=128
test_j=37

print(test_i,test_j, bathy_tem[test_j,test_i],  isf_tem[test_j,test_i])

if (l_verb): print("remap geometry to NEMO")
for j in range(ny_n):
  if (l_verb): print('\r',j,'  ', end=' ')
  sys.stdout.flush()
  for i in range(nx_n):

    isf_cav[j,i], bathy_ism[j,i] = isf_from_bisicles_median(bathy,isf, \
                    bn_map.nmap[j,i],bn_map.x[j,i,:],bn_map.y[j,i,:], \
                    bathy_threshold,False)

    if (i == test_i) & (j==test_j) :
        print(test_i,test_j, bn_map.nmap[j,i], bathy_tem[j,i],  isf_tem[j,i])
        print(test_i,test_j, bathy[int(bn_map.y[j,i,0]),int(bn_map.x[j,i,0])],isf[int(bn_map.y[j,i,0]),int(bn_map.x[j,i,0])])
        print(test_i,test_j, bathy[int(bn_map.y[j,i,1]),int(bn_map.x[j,i,1])],isf[int(bn_map.y[j,i,1]),int(bn_map.x[j,i,1])])
        print(test_i,test_j, bathy[int(bn_map.y[j,i,2]),int(bn_map.x[j,i,2])],isf[int(bn_map.y[j,i,2]),int(bn_map.x[j,i,2])])

    #if isf_cav[j,i] < -90000: #no/bad data in mapping

print(test_i,test_j, isf_cav[test_j,test_i])

ny_f=np.shape(isf_tem)[0]
nx_f=np.shape(isf_tem)[1]

bathy = bathy_tem

isf_new=np.zeros([ny_f,nx_f]) # default open

update=np.where(isf_cav > thresh) # do non-open, non-grounded points
isf_new[0:ny_n,0:nx_n][update] = bathy[0:ny_n,0:nx_n][update] - isf_cav[update]

grounded  = np.where( (isf_cav < 0) & (isf_cav > -9000) ) #ground where the flag says
isf_new[0:ny_n,0:nx_n][grounded] = bathy[0:ny_n,0:nx_n][grounded]

isf_new = np.maximum(isf_new,0) #sanity, safety

isf = isf_new

print(test_i,test_j, bathy[test_j,test_i],isf[test_j,test_i])

isf = enforce_minima(isf,thresh,mindep_isf)

mask_tem = make_surface_mask(bathy_tem,isf_tem,mindep_isf,mindep_bath,thresh)
mask     = make_surface_mask(bathy,isf,mindep_isf,mindep_bath,thresh)
fillmask = mask-mask_tem

##old data is open, new isn't; remove shelf, possibly deepen bathy
closed=np.where((fillmask > thresh))
if l_verb:  print(np.shape(closed)[1], "were open ocean, now closed - opening")
isf[closed] = 0

closed=np.where((fillmask > thresh) & (bathy<mindep_bath))
if l_verb:  print(np.shape(closed)[1],"of which now have land, not ice")
#if l_bath: bathy[closed] = mindep_bath
#if l_bath: bathy[closed] = 12.

##new data is open, old isn't, bathymetry is shallow ; create shelf and ground it
open_shallow=np.where((fillmask < -thresh) & (bathy<mindep_isf+mincav) & (isf_cav > -9000) )
if l_verb:  print(np.shape(open_shallow)[1],"were closed ocean but shallow, BISICLES wants a deeper cavity than bathymetry allows, general mismatch - creating grounded ice")
print(np.shape(open_shallow)[1],"were closed ocean, now open and shallow - creating grounded ice")
isf[open_shallow] = bathy[open_shallow]

##new data is open, old isn't, bathymetry is shallow ; create shelf and ground it
open_shallow=np.where((fillmask < -thresh) & (bathy<mindep_isf+mincav) )
if l_verb:  print(np.shape(open_shallow)[1],"were closed ocean, now open and shallow - creating grounded ice")
isf[open_shallow] = bathy[open_shallow]

##new data is open, old isn't, bathymetry is deep ; create shelf
open_deep=np.where((fillmask < -thresh) & (bathy>=mindep_isf+mincav) )
if l_verb:  print(np.shape(open_deep)[1],"were closed ocean, now open and deep - creating a cavity")
#isf[open_deep] = mindep_isf
isf[open_deep] = 12.

#new data has thin cavity; ground it
thin_cavity=np.where( (isf > thresh) & (np.abs(bathy-isf)>thresh) & (np.abs(bathy-isf)<mincav) & (isf_cav > -9000) )
if l_verb:  print(np.shape(thin_cavity)[1],"have a cavity, but I'm grounding as it's thinner than ",mincav," m")
isf[thin_cavity]=bathy[thin_cavity]

print(test_i,test_j, bathy[test_j,test_i],isf[test_j,test_i])

#write to netCDF file
output= Dataset(rebin_ncfile,'w',format='NETCDF3_CLASSIC')
for dimname in template.dimensions:
   dim=template.dimensions[dimname]

   output.createDimension(dimname,len(dim))

bathyout=output.createVariable('Bathymetry_isf',bathy_tem_v.dtype,bathy_tem_v.dimensions)
isfout  =output.createVariable('isf_draft',  isf_tem_v.dtype,  isf_tem_v.dimensions)

#APPLY ACTIVE MASK
inactive = np.where(isf_cav < -9000)
active = np.where(isf_cav > -1)

isf[inactive] = isf_tem[inactive]

print(test_i,test_j, bathy[test_j,test_i],isf[test_j,test_i])

bathyout[:]=bathy[:]
isfout[:]=isf[:]

output.close()
