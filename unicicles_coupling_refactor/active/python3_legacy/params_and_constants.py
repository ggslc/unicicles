import numpy as np

#this matches the value in Glimmer
rho_ice_glimmer=910.

#these matche the value in the UM/JULES run
nsmax  =10
dzsnow=np.asarray([0.04,0.12,0.34,0.5,0.5,0.5,1.,1.,2.,2.])

planet_radius = 6371229.0
pi            = 3.14159265358979323846


#these match jules/src/control/shared/max_dimensions.F90
max_number_elevs=25

#these match jules/src/params/um/land_tile_ids.F90
start_tile_id_elev_ice=901
start_tile_id_elev_rock=926
start_tile_id_snowlayers=1*1000
start_tile_id_elev_ice_snowlayers=start_tile_id_elev_ice*1000

#these match um/src/control/stash/um_stashcode_mod.F90
stashcode_orog_x_grad       =  5
stashcode_orog_y_grad       =  6
stashcode_sil_orog_rough    = 17
stashcode_hlf_pk_to_trf_ht  = 18
stashcode_snow_amount       = 23
stashcode_lsm               = 30
stashcode_orog              = 33
stashcode_orog              = 33
stashcode_orog_var          = 34
stashcode_orog_gdxx         = 35
stashcode_orog_gdxy         = 36
stashcode_orog_gdyy         = 37
stashcode_vol_smc_sat       = 43
stashcode_frac_surf_type    = 216
stashcode_rgrain            = 231
stashcode_snow_tile         = 240
stashcode_snowdep_grd_tile  = 376
stashcode_snowpack_bk_dens  = 377
stashcode_nsnow_layrs_tiles = 380
stashcode_snow_laythk_tiles = 381
stashcode_snow_ice_tile     = 382
stashcode_snow_liq_tile     = 383
stashcode_snow_t_tile       = 384
stashcode_snow_laydns_tiles = 385
stashcode_snow_grnsiz_tiles = 386
stashcode_land_frac         = 505
stashcode_tsurf_elev_surft  = 576
