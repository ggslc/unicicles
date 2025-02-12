import numpy as np
from params_and_constants import *
from calculate_areas      import calculate_areas
from um_to_np_utility     import ice_to_np,      \
                                 um_to_np_3d,    \
                                 shift_tile_dim

def conservation_of_ice(um_dump,toice_file, fromice_file, frac_elev0, snicemass_elev0):

  smb_ice          =ice_to_np(toice_file, 'ice_smb', ns_flip=False)
  snicemass_elev,_ =um_to_np_3d(um_dump, stashcode_snow_tile,      elev=True)
  frac_elev,pslev  =um_to_np_3d(um_dump, stashcode_frac_surf_type, elev=True)

  ice_elevs=np.where( (pslev >= start_tile_id_elev_ice)  \
                    & (pslev <  start_tile_id_elev_rock) \
                    )

  ice_calving      =ice_to_np(fromice_file, 'cell_calving_flux')
  delta_ice_vol    =ice_to_np(fromice_file, 'change_in_ice_volume')

  fracmask_old2d=np.sum(frac_elev0,axis=2)
  fracmask_new2d=np.sum(frac_elev, axis=2)

  nt=np.shape(frac_elev)[2]

  fracmask_old3d =shift_tile_dim(np.tile(fracmask_old2d,[nt,1,1]))
  fracmask_new3d =shift_tile_dim(np.tile(fracmask_new2d,[nt,1,1]))

  icepts_old2d=np.where(fracmask_old2d > 0.5)
  icepts_old3d=np.where(fracmask_old3d > 0.5)
  icepts_new2d=np.where(fracmask_new2d > 0.5)
  icepts_new3d=np.where(fracmask_new3d > 0.5)

  tileareas_old,_ = calculate_areas(um_dump,frac_elev0)
  tileareas_new,_ = calculate_areas(um_dump,frac_elev)

  #temp=np.empty_like(snicemass_elev)
  #temp[icemask_old3d]=snicemass_elev0[icemask_old3d]*frac_elev0[icemask_old3d]
  #snicemass_gbm0[icemask_old2d]=np.sum(temp,axis=2)[icemask_old2d]
  #print "old",np.sum(snicemass_elev0[icemask_old3d]*frac_elev0[icemask_old3d])

  #temp=np.empty_like(smb_ice)
  #temp[icemask_old3d]=smb_ice[icemask_old3d]*frac_elev0[icemask_old3d]
  #smb_gbm[icemask_old2d]=np.sum(temp,axis=2)[icemask_old2d]
  #print "flux",np.sum(smb_ice[icemask_old3d]*frac_elev0[icemask_old3d])

  #temp=np.empty_like(snicemass_elev)
  #temp[icemask_new3d]=snicemass_elev[icemask_new3d]*frac_elev[icemask_new3d]
  #snicemass_gbm[icemask_new2d]=np.sum(temp,axis=2)[icemask_new2d]
  #print "new",np.sum(snicemass_elev[icemask_new3d]*frac_elev[icemask_new3d])
