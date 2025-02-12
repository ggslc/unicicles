import numpy as np
from params_and_constants import *
from um_to_np_utility import ice_to_np,   \
                             shift_tile_dim, \
                             um_to_np_2d, \
                             um_to_np_3d, \
                             um_to_np_4d, \
                             np_to_um_2d, \
                             np_to_um_3d, \
                             np_to_um_4d

from ctypes import *
libsnow = cdll.LoadLibrary("libsnowpack_manipulations.so")

adjust_fort=libsnow.adjust_stdalone_
adjust_fort.argtypes = [POINTER(c_int)   \
                       ,POINTER(c_int)   \
                       ,POINTER(c_int)   \
                       ,POINTER(c_int)   \
                       ,POINTER(c_float) \
                       ,POINTER(c_int)   \
                       ,POINTER(c_int)   \
                       ,POINTER(c_float) \
                       ,POINTER(c_float) \
                       ,POINTER(c_int)   \
                       ,POINTER(c_float) \
                       ,POINTER(c_float) \
                       ,POINTER(c_float) \
                       ,POINTER(c_float) \
                       ,POINTER(c_float) \
                       ,POINTER(c_float) \
                       ,POINTER(c_float) \
                       ,POINTER(c_float) \
                       ,POINTER(c_float) \
                       ,POINTER(c_float) \
                       ,POINTER(c_float) \
                       ,POINTER(c_float) \
                       ]

def np_to_fortpoint(array,dtype=c_float):
 array_f =np.asfortranarray(array,dtype=dtype)
 return array_f, array_f.ctypes.data_as(POINTER(dtype))

def fort_to_np(array,dtype=np.float):
 return np.ascontiguousarray(array,dtype=dtype)

def adjust_snow(um_dump,fromice_file,tile_frac_old):

#should be no need to change the gridbox mean mass field - the whole idea
#is to conserve that. NO, THAT's ONLY FOR WHEN THE JUST THE FRACTIONS CHANGE NOW. 
#IF GLIMMER HAS MADE/ELIMINATED ICE, THE GBM WILL CHANGE!

#GBM
  soil_vsat_2d   =um_to_np_2d(um_dump, stashcode_vol_smc_sat)
  snicemass_gbm  =um_to_np_2d(um_dump, stashcode_snow_amount)

# on just rock elevs. Changes Glint has made
  nonice_snowdepth_ice = ice_to_np(fromice_file, 'nonice_snowdepth')

# on ice and rock elevs
  tile_frac,pslev =um_to_np_3d(um_dump, stashcode_frac_surf_type,    elev=True)
  nsnow,_         =um_to_np_3d(um_dump, stashcode_nsnow_layrs_tiles, elev=True)
  snow_tile,_     =um_to_np_3d(um_dump, stashcode_snow_tile,         elev=True)
  snowdepth,_     =um_to_np_3d(um_dump, stashcode_snowdep_grd_tile,  elev=True)
  rho_snow_grnd,_ =um_to_np_3d(um_dump, stashcode_snowpack_bk_dens,  elev=True)
  rgrain,_        =um_to_np_3d(um_dump, stashcode_rgrain,            elev=True)

# on snowlayers, on ice and rock elevs
  rgrainl,_  =um_to_np_4d(um_dump, stashcode_snow_grnsiz_tiles, elev=True)
  ds,_       =um_to_np_4d(um_dump, stashcode_snow_laythk_tiles, elev=True)
  sice,_     =um_to_np_4d(um_dump, stashcode_snow_ice_tile,     elev=True)
  sliq,_     =um_to_np_4d(um_dump, stashcode_snow_liq_tile,     elev=True)
  rho_snow,_ =um_to_np_4d(um_dump, stashcode_snow_laydns_tiles, elev=True)
  tsnow,_    =um_to_np_4d(um_dump, stashcode_snow_t_tile,       elev=True)

  nx=np.shape(tile_frac)[0] #this label is misleading, arrays actually [lat,lon,tile]. Use is
  ny=np.shape(tile_frac)[1] #consistent from here on through the fortran despite this
  ntiles=np.shape(tile_frac)[2]
  #nsmax in params_and_constants
  max_elev_id=np.max(pslev)
  elev_levels=np.mod(max_elev_id, max_number_elevs)
#we've only pulled off the elevated tiles of these fields
  elev_start=1

  soil_vsat_3d   =shift_tile_dim(np.tile(soil_vsat_2d,[ntiles,1,1]))
  icemask2d=np.where(soil_vsat_2d[:,:] == 0.)
  icemask3d_alltype=np.where(soil_vsat_3d[:,:,:] == 0.)

# the library actually expects *everything* in on ntiles - need to expand field
# with the first nelevs (the ice_elev parts) zeroed
  nonice_snowdepth = np.zeros([nx,ny,ntiles])
  nonice_snowdepth[:,:,elev_start+elev_levels-1:]=nonice_snowdepth_ice

  smb_from_nisnow = np.zeros([nx,ny,ntiles])

  #make fortran-ordered versions of the numpy arrays, with
  #pointers for passing
  #IN
  dzsnow_f,dzsnow_fp                    =np_to_fortpoint(dzsnow)
  tile_frac_old_f,tile_frac_old_fp      =np_to_fortpoint(tile_frac_old)

  tile_frac_f,tile_frac_fp              =np_to_fortpoint(tile_frac)
  nonice_snowdepth_f,nonice_snowdepth_fp=np_to_fortpoint(nonice_snowdepth)
  
  #INOUT, on tiles
  nsnow_f,nsnow_fp                      =np_to_fortpoint(nsnow,dtype=c_int)
  snow_tile_f,snow_tile_fp              =np_to_fortpoint(snow_tile)
  snowdepth_f,snowdepth_fp              =np_to_fortpoint(snowdepth)
  rho_snow_grnd_f,rho_snow_grnd_fp      =np_to_fortpoint(rho_snow_grnd)
  rgrain_f,rgrain_fp                    =np_to_fortpoint(rgrain)
  
  #INOUT, on snowlayers
  ds_f,ds_fp                            =np_to_fortpoint(ds)
  sliq_f,sliq_fp                        =np_to_fortpoint(sliq)
  sice_f,sice_fp                        =np_to_fortpoint(sice)
  rho_snow_f,rho_snow_fp                =np_to_fortpoint(rho_snow)
  rgrainl_f,rgrainl_fp                  =np_to_fortpoint(rgrainl)
  tsnow_f,tsnow_fp                      =np_to_fortpoint(tsnow)

  #OUT, on tiles
  smb_from_nisnow_f,smb_from_nisnow_fp = np_to_fortpoint(smb_from_nisnow)
 
  print "calling fortran lib"
  dummy=adjust_fort(c_int(nx)           \
                   ,c_int(ny)           \
                   ,c_int(ntiles)       \
                   ,c_int(nsmax)        \
                   ,dzsnow_fp           \
                   ,c_int(elev_start)   \
                   ,c_int(elev_levels)  \
                   ,tile_frac_fp        \
                   ,tile_frac_old_fp    \
                   ,nsnow_fp            \
                   ,nonice_snowdepth_fp \
                   ,snow_tile_fp        \
                   ,snowdepth_fp        \
                   ,rho_snow_grnd_fp    \
                   ,rgrain_fp           \
                   ,sice_fp             \
                   ,sliq_fp             \
                   ,rgrainl_fp          \
                   ,tsnow_fp            \
                   ,rho_snow_fp         \
                   ,ds_fp               \
                   ,smb_from_nisnow_fp  \
                   )
  print "back from fortran lib"
  
  #convert the modified fortran-style arrays back to 
  #regular numpy ones
  #OUT, on tiles
  nsnow           = fort_to_np(nsnow_f)
  snow_tile       = fort_to_np(snow_tile_f)
  snowdepth       = fort_to_np(snowdepth_f)
  rho_snow_grnd   = fort_to_np(rho_snow_grnd_f)
  rgrain          = fort_to_np(rgrain_f)
  smb_from_nisnow = fort_to_np(smb_from_nisnow_f)

  #smb_from_nisnow needs to be dimensioned (10,144,192) to go back out to the next year's coupling
  smb_from_nisnow=np.rollaxis(smb_from_nisnow,2,0)
  smb_from_nisnow=smb_from_nisnow[0:10,:,:]
  
  #OUT, on snowlayers
  ds       = fort_to_np(ds_f)
  sliq     = fort_to_np(sliq_f)
  sice     = fort_to_np(sice_f)
  rho_snow = fort_to_np(rho_snow_f)
  rgrainl  = fort_to_np(rgrainl_f)
  tsnow    = fort_to_np(tsnow_f)

  #splicing isn't needed here, I don't think. The tile fractions we're altering
  #by have already been done, and the change in ni_snow coming out of glimmer is
  #weighted according the gridbox frac it occupies too
 
  #AS LONG AS AREAS OUTSIDE OF THE UNICICLES DOMAIN AREN'T BEING AFFECTED. THINK OK.?

  #reconstruct GBM
  temp=np.empty_like(snow_tile)
  temp[icemask3d_alltype]=snow_tile[icemask3d_alltype]*tile_frac[icemask3d_alltype]
#To do the GBM, the above needs to looks at ice and rock elevations, not just icemask_3d. 
  snicemass_gbm[icemask2d]=np.sum(temp,axis=2)[icemask2d]

  #GBM
  um_dump = np_to_um_2d(um_dump, stashcode_snow_amount,       snicemass_gbm)

  um_dump = np_to_um_3d(um_dump, stashcode_nsnow_layrs_tiles, nsnow,         elev=True)
  um_dump = np_to_um_3d(um_dump, stashcode_snow_tile,         snow_tile,     elev=True)
  um_dump = np_to_um_3d(um_dump, stashcode_snowdep_grd_tile, snowdepth,     elev=True)
  um_dump = np_to_um_3d(um_dump, stashcode_snowpack_bk_dens,  rho_snow_grnd, elev=True)
  um_dump = np_to_um_3d(um_dump, stashcode_rgrain,            rgrain,        elev=True)

  um_dump = np_to_um_4d(um_dump, stashcode_snow_laythk_tiles, ds,       elev=True)
  um_dump = np_to_um_4d(um_dump, stashcode_snow_ice_tile,     sice,     elev=True)
  um_dump = np_to_um_4d(um_dump, stashcode_snow_liq_tile,     sliq,     elev=True)
  um_dump = np_to_um_4d(um_dump, stashcode_snow_laydns_tiles, rho_snow, elev=True)
  um_dump = np_to_um_4d(um_dump, stashcode_snow_grnsiz_tiles, rgrainl,  elev=True)
  um_dump = np_to_um_4d(um_dump, stashcode_snow_t_tile,       tsnow,    elev=True)
  
  return um_dump, smb_from_nisnow
