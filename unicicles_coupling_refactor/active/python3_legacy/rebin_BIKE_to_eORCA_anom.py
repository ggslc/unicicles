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
  parser.add_argument("--new_bathy",   help="write a new NEMO bathymetry as well as the icedraft")
  parser.add_argument("--verbose" )
  args = parser.parse_args()

  err = 0
  l_verb = False
  if args.verbose: l_verb = True

  l_bath = False
  if args.new_bathy: l_bath = True


  hdf5_input, err  = arg_to_file_exist(args.hdf5_input, mandatory=True, err=err)
  nc_template, err = arg_to_file_exist(args.nc_template, mandatory=True, err=err)
  nc_output, err   = arg_to_file_exist(args.nc_output, mandatory=True, io="out", err=err)

  if err > 0:
    parser.print_help()
    sys.exit(2)

  return hdf5_input, nc_template, nc_output, l_verb, l_bath

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

  isf_out = 99999. #use for detection of out-of-domain points later
  bathy_out = 99999. 

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
      isf_out = -9999. #grounded

    if l_verb: print(isf_acc/n_contrib, isf_out)
  
  return isf_out, bathy_out

def splice_field_distance(full_bath,ice_bath, ice_isf, merge_mask=None):
    from scipy import ndimage

    if merge_mask == None:
      merge_mask = np.zeros_like(ice_isf)
      is_ice = np.where( (ice_isf > thresh) & (ice_isf < 1e4) )
      merge_mask[is_ice] = 1.

    #we could square this? Or multiply by a factor? Should it drop off /slower/?
    ice_dist_int = ndimage.distance_transform_edt(1-merge_mask)
    ice_dist = ice_dist_int.astype(np.float)
    is_ice = np.where(ice_dist < 1.)
    ice_dist[is_ice] = 1.
    #ice_weights = 1. / ice_dist  # drops off too fast?
    linear_distance_threshold=10. # linear dropoff over specified distance
    ice_weights = np.maximum( (1-(ice_dist-1.)/linear_distance_threshold) , 0)

    bad_ice_bath = np.where(ice_bath > 1e4)
    ice_weights[bad_ice_bath] = 0.

    ny = np.shape(ice_bath)[0]
    nx = np.shape(ice_bath)[1]

    full_bath_new = np.copy(full_bath)

    full_bath_new[:ny,:nx] = full_bath[:ny,:nx]*(1-ice_weights) + ice_bath*ice_weights*-1

    full_bath_new = np.maximum(full_bath_new,0.)

    return full_bath_new


hdfplotfile, nemo_templatefile, rebin_ncfile, l_verb, l_bath = parse_commandline(sys.argv)

#nemo_templatefile="BIKE_geom_rebin_eORCA025_grid.nc_match_12.0_30.nc"
bisicles_to_nemo_mapping_file="bikegridmap_file.map2d"

bathy_threshold=-5000
isf_threshold=-9999

level=2 # level of grid refinement
order=1 # interpolation order, 0 for piecewise constant, 1 for linear

if (l_verb): print('reading hdf file')
amrID = amrio.load(hdfplotfile)
lo,hi = amrio.queryDomainCorners(amrID, level)
_,_,isf = amrio.readBox2D(amrID, level, lo, hi, "Z_bottom", order)
_,_,bathy = amrio.readBox2D(amrID, level, lo, hi, "Z_base", order)

##somewhere in the tangle reading the hdf5 this way gets me arrays transposed compared
##to the way they were read/used when generating the lats and lons 
bathy=np.transpose(bathy)
isf=np.transpose(isf)
##

if (l_verb): print('restoring mapping files')
bn_map=mp.load(bisicles_to_nemo_mapping_file)

ny_n=np.shape(bn_map.x)[0]
nx_n=np.shape(bn_map.x)[1]

isf_cav=np.zeros([ny_n,nx_n])
bathy_ism=np.zeros([ny_n,nx_n])

if (l_verb): print("remap geometry to NEMO")
for j in range(ny_n):
  if (l_verb): print('\r',j,'  ', end=' ')
  sys.stdout.flush()
  for i in range(nx_n):

    isf_cav[j,i], bathy_ism[j,i] = isf_from_bisicles_median(bathy,isf, \
                    bn_map.nmap[j,i],bn_map.x[j,i,:],bn_map.y[j,i,:], \
                    bathy_threshold,False)

    if isf_cav[j,i] > 90000: #no/bad data in mapping
        isf_cav[j,i] =0.     #default open ocean
        #300 is representative of ORCA025, but will work for ORCA1 too
        if j<300: isf_cav[j,i] = -9999. #this far south, any bad data is eORCA grid weirdness
                                        #shouldn't matter, but nicer if we ground it
                                        #to whatever the bathymetry array has



template = Dataset(nemo_templatefile,'r')
bathy_tem_v = template.variables['Bathymetry_isf']
isf_tem_v = template.variables['isf_draft']
bathy_tem = bathy_tem_v[:]
isf_tem = isf_tem_v[:]

ny_f=np.shape(isf_tem)[0]
nx_f=np.shape(isf_tem)[1]


if (l_bath): 
  bathy = splice_field_distance(bathy_tem,bathy_ism, isf_cav)
  bathy[:,-2:]=bathy[:,0:2] #wrap for haloe
else:
  bathy = bathy_tem
  bathy[:,-2:]=bathy[:,0:2] #wrap for haloe

isf_new=np.zeros([ny_f,nx_f]) # default open

update=np.where(isf_cav > thresh) # do non-open, non-grounded points
isf_new[0:ny_n,0:nx_n][update] = bathy[0:ny_n,0:nx_n][update] - isf_cav[update]

grounded  = np.where(isf_cav < 0) #ground where the flag says
isf_new[0:ny_n,0:nx_n][grounded] = bathy[0:ny_n,0:nx_n][grounded]

isf_new[:,-2:]=isf_new[:,0:2] #wrap for haloes
isf_new = np.maximum(isf_new,0) #sanity, safety

isf = isf_new

if l_bath: bathy = enforce_minima(bathy,thresh,mindep_bath)
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
if l_bath: bathy[closed] = 12.

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
thin_cavity=np.where( (isf > thresh) & (np.abs(bathy-isf)>thresh) & (np.abs(bathy-isf)<mincav) )
if l_verb:  print(np.shape(thin_cavity)[1],"have a cavity, but I'm grounding as it's thinner than ",mincav," m")
isf[thin_cavity]=bathy[thin_cavity]

#wrap the haloes
bathy[:,-2:] = bathy[:,0:2]
isf[:,-2:]   = isf[:,0:2]

#write to netCDF file
output= Dataset(rebin_ncfile,'w',format='NETCDF3_CLASSIC')
for dimname in template.dimensions:
   dim=template.dimensions[dimname]

   output.createDimension(dimname,len(dim))

bathyout=output.createVariable('Bathymetry_isf',bathy_tem_v.dtype,bathy_tem_v.dimensions)
isfout  =output.createVariable('isf_draft',  isf_tem_v.dtype,  isf_tem_v.dimensions)

bathyout[:]=bathy[:]
isfout[:]=isf[:]

output.close()




