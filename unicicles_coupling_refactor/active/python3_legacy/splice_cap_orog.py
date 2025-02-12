import numpy as np
from params_and_constants import *
from um_to_np_utility import um_to_np_2d, \
                             ice_to_np, \
                             get_ice_um_grid_overlap, \
                             merge_ice_um_arrays, \
                             np_to_um_2d

def splice_cap_orog(um_dump,um_ref_dump,cap_file,ice_file, x_offset=0, y_offset=0):

  #pull UM coastal fraction
  coast_um=um_to_np_2d(um_ref_dump,stashcode_land_frac)

  #get overlap between UM and ice grids
  fmask_ice=get_ice_um_grid_overlap(ice_file)

  #get BISICLES orog on the UM grid - this is the best thing to use as a mask
  #for where the ISM domain has real data worth adding
  orog_ism=ice_to_np(ice_file, 'surface_elevation')

  #update where the ice grid has data and the UM mask says there's some land
  #we need to include the coastal fraction weighting in here too?
  ice_updates=np.where((orog_ism > 0.) & (coast_um >0))

  for stash in [stashcode_orog,           \
                stashcode_orog_var,       \
                stashcode_orog_x_grad,    \
                stashcode_orog_y_grad,    \
                stashcode_orog_gdxx,      \
                stashcode_orog_gdxy,      \
                stashcode_orog_gdyy,      \
                stashcode_sil_orog_rough, \
                stashcode_hlf_pk_to_trf_ht]:

      field_um_ref  =um_to_np_2d(um_ref_dump, stash)
      field_cap_reg =um_to_np_2d(cap_file,stash)

      #modern CAP doesn't recalculate these.
      #NEED TO GET ANSWER FROM STUART WEBSTER(?)
      #old CAP does HPTH=SD, SOR=SD*res_dep_factor
      #res_dep_factor =5.9e-5 for N24
      #res_dep_factor =8.4e-5 for N48
      #guess:
      #res_dep_factor=1e-4

      #eyeballing G'land and Ant cf N96 orig ancil, compromising
      #doesn't get the highs, over estimates the lows
      #in general, factors are too small for Gland and too high for Ant
      SOR_res_dep_factor=5e-5
      HPT_res_dep_factor=0.15

      if stash==stashcode_sil_orog_rough:
         field_cap_reg =um_to_np_2d(cap_file,stashcode_orog_var)*SOR_res_dep_factor

      if stash==stashcode_hlf_pk_to_trf_ht:
         field_cap_reg =um_to_np_2d(cap_file,stashcode_orog_var)*HPT_res_dep_factor

      #CAP fields may be on a regional subgrid
      ny_reg=np.shape(field_cap_reg)[0]
      nx_reg=np.shape(field_cap_reg)[1]
      field_cap_glob=np.zeros_like(field_um_ref)

      field_cap_glob[y_offset:y_offset+ny_reg,x_offset:x_offset+nx_reg] \
                    = field_cap_reg

  #splice in the CAP-processed fields, weighted by the UM/ice grid overlap and the UM land fraction
      field_spliced=merge_ice_um_arrays(field_um_ref,       \
                                        field_cap_glob,     \
                                        fmask_ice*coast_um, \
                                        points=ice_updates,  \
                                        noFrac=True)

      um_dump=np_to_um_2d(um_dump, stash, field_spliced)

  return um_dump
