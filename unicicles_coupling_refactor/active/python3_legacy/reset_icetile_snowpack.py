import numpy as np
from params_and_constants import *
from um_to_np_utility import shift_tile_dim, \
                             um_to_np_2d,    \
                             um_to_np_3d,    \
                             um_to_np_4d,    \
                             ice_to_np,      \
                             np_to_um_2d,    \
                             np_to_um_3d,    \
                             np_to_um_4d

def reset_icetile_snowpack(um_dump, ice_file, coupling_period_in_secs,real_tile_only=False):

#GBM
  soil_vsat_2d   =um_to_np_2d(um_dump, stashcode_vol_smc_sat)
  snicemass_gbm  =um_to_np_2d(um_dump, stashcode_snow_amount)
#TILES
  nsnice_elev,pslev    =um_to_np_3d(um_dump, stashcode_nsnow_layrs_tiles, elev=True)
  snicemass_elev,_     =um_to_np_3d(um_dump, stashcode_snow_tile,         elev=True)
  snicedep_elev,_      =um_to_np_3d(um_dump, stashcode_snowdep_grd_tile,  elev=True)
  snicerho_elev,_      =um_to_np_3d(um_dump, stashcode_snowpack_bk_dens,  elev=True)
  frac_elev,_          =um_to_np_3d(um_dump, stashcode_frac_surf_type,    elev=True)
#LAYERS
  snice_layer_mass,_   =um_to_np_4d(um_dump, stashcode_snow_ice_tile,     elev=True)
  snice_layer_dep,_    =um_to_np_4d(um_dump, stashcode_snow_laythk_tiles, elev=True)
  snice_layer_den,_    =um_to_np_4d(um_dump, stashcode_snow_laydns_tiles, elev=True)

  nt=len(pslev)
  max_elev_id=np.max(pslev)
  nelev=np.mod(max_elev_id, max_number_elevs)

  soil_vsat_3d   =shift_tile_dim(np.tile(soil_vsat_2d,[nt,1,1]))
  #pick out all ice gridboxes
  icemask2d=np.where(soil_vsat_2d[:,:] == 0.)
  #pick out the *ice* tiles of the elevated ones - only 0:nelev are looked at!
  icemask3d=np.where(soil_vsat_3d[:,:,0:nelev] == 0.)
  icemask3d_alltype=np.where(soil_vsat_3d[:,:,0:2*nelev] == 0.)

  if real_tile_only: icemask3d = np.where(frac_elev[:,:,0:nelev] > 1e-3)

  #need to convert smb kg/m2/s to the total mass that was given to the
  #icesheet for this cycle. ns_flip becase the toice file is the same NS
  #as the UM, unlike the fromice
  smb_ice=ice_to_np(ice_file, 'ice_smb', ns_flip=False)*coupling_period_in_secs
  
#sanity check - do we have nsmax snow on all ice tiles?
  if np.any(nsnice_elev[icemask3d] < nsmax):
    print("3d:nsnow is no longer 10 on all the ice tiles I'm looking at!")
    print("if you *know* this is only on virtual tiles and you don't mind that,")
    print("set real_tile_only=True calling this script to ignore virtual tiles")
    print("if real_tile_only is already True and you're reading this, you've got real problems")
    exit(1)

  snice_layer_mass[:,:,:,-1][icemask3d] = snice_layer_mass[:,:,:,-1][icemask3d] - smb_ice[icemask3d]
  snice_layer_dep[:,:,:,-1][icemask3d]  = snice_layer_dep[:,:,:,-1][icemask3d]  - (smb_ice[icemask3d]/rho_ice_glimmer)
  snice_layer_den[:,:,:,-1][icemask3d]  = snice_layer_mass[:,:,:,-1][icemask3d]/snice_layer_dep[:,:,:,-1][icemask3d]

#make the same adjustment to the tile total mass
  snicemass_elev[icemask3d]=snicemass_elev[icemask3d]-smb_ice[icemask3d]
#make the same adjustment to the tile total depth
  snicedep_elev[icemask3d]=snicedep_elev[icemask3d]-(smb_ice[icemask3d]/rho_ice_glimmer)
#make the same adjustment to the tile mean density
  snicerho_elev[icemask3d]=snicemass_elev[icemask3d]/snicedep_elev[icemask3d]

  temp=np.empty_like(snicemass_elev)
  temp[icemask3d_alltype]=snicemass_elev[icemask3d_alltype]*frac_elev[icemask3d_alltype]
#To do the GBM, the above needs to looks at ice and rock elevations, not just icemask_3d. 
  snicemass_gbm[icemask2d]=np.sum(temp,axis=2)[icemask2d]

#GBM
  um_dump = np_to_um_2d(um_dump, stashcode_snow_amount,       snicemass_gbm)
#TILES
  um_dump = np_to_um_3d(um_dump, stashcode_snow_tile,         snicemass_elev, elev=True)
  um_dump = np_to_um_3d(um_dump, stashcode_snowdep_grd_tile,  snicedep_elev,  elev=True)
  um_dump = np_to_um_3d(um_dump, stashcode_snowpack_bk_dens,  snicerho_elev,  elev=True)
#LAYERS
  um_dump = np_to_um_4d(um_dump, stashcode_snow_ice_tile,     snice_layer_mass, elev=True)
  um_dump = np_to_um_4d(um_dump, stashcode_snow_laythk_tiles, snice_layer_dep,  elev=True)
  um_dump = np_to_um_4d(um_dump, stashcode_snow_laydns_tiles, snice_layer_den,  elev=True)

  return um_dump
