import numpy as np
from params_and_constants import *
from um_to_np_utility import um_to_np_2d,             \
                             get_ice_um_grid_overlap, \
                             ice_to_np,               \
                             merge_ice_um_arrays,     \
                             np_to_um_2d

def splice_mean_orog(um_dump,um_ref_dump,ice_file):

  #pull UM coastal fraction
  coast_um=um_to_np_2d(um_ref_dump,stashcode_land_frac)

  #get overlap between UM and ice grids
  fmask_ice=get_ice_um_grid_overlap(ice_file)

  #update where the ice grid has data and the UM mask says there's some land
  #we need to include the coastal fraction weighting in here too?
  ice_updates=np.where((fmask_ice > 0.) & (coast_um >0))

  #pull orography as a 2d np array
  orog_um=um_to_np_2d(um_ref_dump,stashcode_orog)
  #pull 2d GBM orog from ice
  orog_ice=ice_to_np(ice_file, 'surface_elevation')

  #splice in the new orography, weighted by the UM/ice grid overlap and the UM land fraction
  #orog_new=np.copy(orog_um) # REDUNDANT?
  orog_new=merge_ice_um_arrays(orog_um,orog_ice,fmask_ice*coast_um,points=ice_updates,noFrac=True)

  um_dump=np_to_um_2d(um_dump, stashcode_orog, orog_new)

  return um_dump
