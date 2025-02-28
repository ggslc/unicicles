import sys
import numpy as np
import mule_rss
from params_and_constants import *
#import cfplot as cfp

def um_to_np_2d(um_dump, stash):
  try:
    index=mule_rss.findindex(um_dump,stash)[0]
  except IndexError:
    print ""
    print "ERROR in um_to_np_2d. Field not found for STASH: ",stash
    print ""
    sys.exit(2)

  return um_dump.fields[index].get_data()


def find_pslevs(um_dump,indices):
  pslev=[]
  for i in indices:
     pslev.append(um_dump.fields[i].lbuser5)

  pslev=np.asarray(pslev)
  return pslev


def select_elev_indices(um_dump,indices):
  pslev=find_pslevs(um_dump,indices)

  indices=indices[np.where( ((pslev >= start_tile_id_elev_ice)   \
                             & (pslev < start_tile_id_snowlayers)) | \
                            (pslev >= start_tile_id_elev_ice_snowlayers)  \
                          )]

  pslev=pslev[np.where( ((pslev >= start_tile_id_elev_ice)   \
                             & (pslev < start_tile_id_snowlayers)) | \
                        (pslev >= start_tile_id_elev_ice_snowlayers)  \
                          )]

  return indices, pslev


def um_to_np_3d(um_dump, stash, elev=False):
  try:
    indices=mule_rss.findindex(um_dump,stash)
  except IndexError:
    print ""
    print "ERROR in um_to_np_3d. Field not found for STASH: ",stash
    print ""
    sys.exit(2)

  pslev =  find_pslevs(um_dump,indices)
  if elev: indices, pslev = select_elev_indices(um_dump,indices)
    
  return mule_rss.get_data3d(um_dump,indices),pslev


def um_to_np_4d(um_dump, stash, elev=False):
  try:
    indices=mule_rss.findindex(um_dump,stash)
  except IndexError:
    print ""
    print "ERROR in um_to_np_4d. Field not found for STASH: ",stash
    print ""
    sys.exit(2)

  pslev =  find_pslevs(um_dump,indices)
  if elev: indices, pslev = select_elev_indices(um_dump,indices)
    
  field3d = mule_rss.get_data3d(um_dump,indices)

  nsl=max(pslev-(pslev//1000)*1000)
  nt =len(pslev)/nsl
  ny =np.shape(field3d)[0]
  nx =np.shape(field3d)[1]

  field4d = np.zeros([ny,nx,nt,nsl])

  for i in range(len(pslev)):
      tile_id=pslev[i]//1000
      snowlayer_index=pslev[i]-(tile_id*1000)-1
      tile_index=np.mod(i,nt)
      field4d[:,:,tile_index,snowlayer_index]=field3d[:,:,i]

  return field4d,pslev


def np_to_um_2d(um_dump, stash, array):
  try:
    index=mule_rss.findindex(um_dump,stash)[0]
  except IndexError:
    print ""
    print "ERROR in np_to_um_2d. Field not found for STASH: ",stash
    print ""
    sys.exit(2)

  um_dump.fields[index]=mule_rss.overwrite_data(um_dump.fields[index],array)

  return um_dump


def np_to_um_3d(um_dump, stash, array, elev=False):
  try:
    indices=mule_rss.findindex(um_dump,stash)
  except IndexError:
    print ""
    print "ERROR in np_to_um_3d. Field not found for STASH: ",stash
    print ""
    sys.exit(2)

  pslev =  find_pslevs(um_dump,indices)
  if elev: indices, pslev = select_elev_indices(um_dump,indices)

  if len(indices) != np.shape(array)[2]:
    print "ERROR in np_to_um_3d. Array doesn't match 3d axis length"
    print ""
    sys.exit(2)
  else:
    #not sure this function shouldn't be returning the um_dump
    #for modifications to register?!
    um_dump = mule_rss.put_data3d(um_dump,indices,array)

  return um_dump


def np_to_um_4d(um_dump, stash, field4d, elev=False):
  try:
    indices=mule_rss.findindex(um_dump,stash)
  except IndexError:
    print ""
    print "ERROR in np_to_um_4d. Field not found for STASH: ",stash
    print ""
    sys.exit(2)

  pslev =  find_pslevs(um_dump,indices)
  if elev: indices, pslev = select_elev_indices(um_dump,indices)

  nsl_dump=max(pslev-(pslev//1000)*1000)
  nt_dump =len(pslev)/nsl_dump
  ny =np.shape(field4d)[0]
  nx =np.shape(field4d)[1]
  nt =np.shape(field4d)[2]
  nsl=np.shape(field4d)[3]

  array = np.zeros([ny,nx,nt*nsl])

  for i in range(len(pslev)):
      tile_id=pslev[i]//1000
      snowlayer_index=pslev[i]-(tile_id*1000)-1
      tile_index=np.mod(i,nt)
      array[:,:,i]=field4d[:,:,tile_index,snowlayer_index]
 
  if len(indices) != np.shape(array)[2]:
    print "ERROR in np_to_um_4d. Array doesn't match 3d axis length"
    print ""
    sys.exit(2)
  else:
    #not sure this function shouldn't be returning the um_dump
    #for modifications to register?!
    um_dump = mule_rss.put_data3d(um_dump,indices,array)

  return um_dump


def ice_to_np(ice_file, ncvar_name,ns_flip=True):
  try: 
    field=ice_file.variables[ncvar_name][:].squeeze()
  except KeyError: 
    print ""
    print "ERROR in ice_to_np: field ",ncvar_name," not found in file"
    print ""
    sys.exit(2)

  if field.ndim>2:
    # use UM convention: ice_file arrays are [(nsnow),ntile,ny,nx], um dumps are [ny,nx,ntile]
    field=shift_tile_dim(field)
  if field.ndim>3:
#   the above has done the snow. Now move the tile fraction back past the 2D
    field=shift_tile_dim(field,end=0)

  #The glint-produced ice_file arrays are N-S flipped wrt UM
  #The atmos.iceceouple file fed *to* glint is the right way up though
  if ns_flip: field=field[::-1,...]

  return field


def shift_tile_dim(field,end=1):
  # shuffle a leading-dimension (ie tiles) to the back, eg
  # [ntile,ny,nx] -> [ny,nx,ntile] (this is my mule convention)
  # this swaps position 1 -> 3
  field=np.swapaxes(field,0,1)
  field=np.swapaxes(field,1,2)

  # a 4D array needs one more swap to go right to the back
  if (field.ndim==4) & (end==1):
    field=np.swapaxes(field,2,3)

  return field


def get_ice_um_grid_overlap(ice_file):
  #pull 3d land fraction from ice. We need this for the possibly 
  #fractional overlap of the ice GCM grids coverage at the edges of the 
  #ice domain (the sum of land_frac_ice over elevs for each grid box tells us this)
  land_frac_ice=ice_to_np(ice_file, 'tile_land_fraction')
  fmask_ice=np.sum(land_frac_ice,axis=2)
  
  return fmask_ice


def merge_ice_um_arrays(um, ice, overlap, points=None, noFrac=False):

  if (noFrac):
    overlap[np.where(np.abs(overlap-1.0)< 1e-4) ] = 1.
    overlap[np.where(np.abs(overlap-1.0)>=1e-4) ] = 0.

  #splice in the new field from the ice, weighted by the UM/ice grid overlap
  if points==None:
    um = ice*overlap + um*(1.-overlap)
  else:
    um[points] = ice[points]*overlap[points] + um[points]*(1.-overlap[points])

    #testing
    #um[points] = -9999.

  return um
