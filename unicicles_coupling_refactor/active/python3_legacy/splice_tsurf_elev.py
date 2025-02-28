import numpy as np
from params_and_constants import *
from um_to_np_utility import ice_to_np,               \
                             get_ice_um_grid_overlap, \
                             um_to_np_3d,             \
                             merge_ice_um_arrays,     \
                             np_to_um_3d

def splice_tsurf_elev(um_dump,um_ref_dump,ice_file):

  #Glimmer deals in degC, UM in K
  tsurf_elev_ice=ice_to_np(ice_file, 'ice_stemp')+273.15
  #tsurf_elev_ice=ice_to_np(ice_file, 'snow_ice_hflux')
  #get BISICLES orog on the UM grid - this is the best thing to use as a mask
  #for where the ISM domain has real data worth adding
  #Could use tile_ice_fraction, levle by level, better?
  orog_ism=ice_to_np(ice_file, 'surface_elevation')

  tsurf_elev_um,_ =um_to_np_3d(um_ref_dump, stashcode_tsurf_elev_surft, elev=True)

  frac_um,pslev   =um_to_np_3d(um_ref_dump, stashcode_frac_surf_type,   elev=True)
  fmask_um2d=np.sum(frac_um, axis=2)

  ##check how many elevation classes we have, and what the max is
  max_elev_id=np.max(pslev)
  nelev_um=np.mod(max_elev_id, max_number_elevs)

  #get overlap between UM and ice grids
  fmask_ice2d=get_ice_um_grid_overlap(ice_file)

  #only update where there is some ice grid coverage and the UM allows elevated tiles 
  frac_updates=np.where((orog_ism > 0.) & (fmask_um2d > 0.5))

  #only want to update the ice elevs, not rock subsurfaces

  #as for orog, mix in the new ice-derived fields with the old UM ones weighted by the
  #partial grid coverage of the ice grid at the edges (fmask_ice)
  tsurf_elev_new=np.copy(tsurf_elev_um)
  for elev in range(nelev_um):
      tsurf_elev_ice_1lev=tsurf_elev_ice[:,:,elev]
      tsurf_elev_um_1lev=tsurf_elev_um[:,:,elev]

      tsurf_elev_um_1lev=merge_ice_um_arrays(tsurf_elev_um_1lev,  \
                                             tsurf_elev_ice_1lev, \
                                             fmask_ice2d,         \
                                             points=frac_updates, \
                                             noFrac=True)

      tsurf_elev_um[:,:,elev]=tsurf_elev_um_1lev

  um_dump=np_to_um_3d(um_dump, stashcode_tsurf_elev_surft, \
                      tsurf_elev_um, elev=True)

  return um_dump
