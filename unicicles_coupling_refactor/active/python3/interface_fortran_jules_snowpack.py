"""
interface_fortran_jules_snowpack

interface between python and fortran so that unicicles_cap_global_to_um can use
original JULES multilayer snowpack code to calculate modifications to snow
prognostics implied by changes in the ice sheet

Called by adjust_snow in unicicles_cap_global_to_um
"""

from ctypes import c_float, c_int, cdll, POINTER
import numpy as np
from common_params_and_constants import dzsnow, max_number_elevs, \
    nsmax, number_elevs, \
    stashcode_frac_surf_type, stashcode_nsnow_layrs_tiles, \
    stashcode_rgrain, stashcode_snow_amount, \
    stashcode_snow_grnsiz_tiles, stashcode_snow_ice_tile, \
    stashcode_snow_laydns_tiles, stashcode_snow_laythk_tiles, \
    stashcode_snow_liq_tile, stashcode_snow_t_tile, \
    stashcode_snow_tile, stashcode_snowdep_grd_tile, \
    stashcode_snowpack_bk_dens, stashcode_vol_smc_sat
from common_um_to_np import ice_to_np, np_to_um_2d, np_to_um_3d, np_to_um_4d, \
    shift_tile_dim, um_to_np_2d, um_to_np_3d, um_to_np_4d
libsnow = cdll.LoadLibrary("libsnowpack_manipulations.so")

adjust_fort = libsnow.adjust_stdalone_
adjust_fort.argtypes = [POINTER(c_int),
                        POINTER(c_int),
                        POINTER(c_int),
                        POINTER(c_int),
                        POINTER(c_float),
                        POINTER(c_int),
                        POINTER(c_int),
                        POINTER(c_float),
                        POINTER(c_float),
                        POINTER(c_int),
                        POINTER(c_float),
                        POINTER(c_float),
                        POINTER(c_float),
                        POINTER(c_float),
                        POINTER(c_float),
                        POINTER(c_float),
                        POINTER(c_float),
                        POINTER(c_float),
                        POINTER(c_float),
                        POINTER(c_float),
                        POINTER(c_float),
                        POINTER(c_float)]


def np_to_fortpoint(array, dtype=c_float):
    #cast numpy array as fortran
    array_f = np.asfortranarray(array, dtype=dtype)
    return array_f, array_f.ctypes.data_as(POINTER(dtype))


def fort_to_np(array, dtype=float):
    #cast fortran array as numpy
    return np.ascontiguousarray(array, dtype=dtype)


def adjust_snow(um_dump, fromice_file, tile_frac_old):

    # gather inputs

    # get simple lon,lat fields
    soil_vsat_2d = um_to_np_2d(um_dump, stashcode_vol_smc_sat)
    snicemass_gbm = um_to_np_2d(um_dump, stashcode_snow_amount)

    # get field from Glint whose third dimension is just rock elev tiles. The Glint file
    # source contains *changes* Glint has made to the rock elev snow depth
    nonice_snowdepth_ice = ice_to_np(fromice_file, 'nonice_snowdepth')

    # get fields whose third dimension is ice and rock elev tiles
    tile_frac, pslev = um_to_np_3d(um_dump, stashcode_frac_surf_type,
                                   elev=True)
    nsnow, _ = um_to_np_3d(um_dump, stashcode_nsnow_layrs_tiles, elev=True)
    snow_tile, _ = um_to_np_3d(um_dump, stashcode_snow_tile, elev=True)
    snowdepth, _ = um_to_np_3d(um_dump, stashcode_snowdep_grd_tile, elev=True)
    rho_snow_grnd, _ = um_to_np_3d(um_dump, stashcode_snowpack_bk_dens,
                                   elev=True)
    rgrain, _ = um_to_np_3d(um_dump, stashcode_rgrain, elev=True)

    # get fields with snowlayers, on ice and rock elevs
    rgrainl, _ = um_to_np_4d(um_dump, stashcode_snow_grnsiz_tiles, elev=True)
    ds, _ = um_to_np_4d(um_dump, stashcode_snow_laythk_tiles, elev=True)
    sice, _ = um_to_np_4d(um_dump, stashcode_snow_ice_tile, elev=True)
    sliq, _ = um_to_np_4d(um_dump, stashcode_snow_liq_tile, elev=True)
    rho_snow, _ = um_to_np_4d(um_dump, stashcode_snow_laydns_tiles, elev=True)
    tsnow, _ = um_to_np_4d(um_dump, stashcode_snow_t_tile, elev=True)

    # Derive the dimensions of teh arrays we're dealing with
    # These labels is misleading, arrays actually [lat,lon,tile]. Use is
    # consistent from here on through the fortran despite this
    nx = np.shape(tile_frac)[0]
    ny = np.shape(tile_frac)[1]
    ntiles = np.shape(tile_frac)[2]
    # nsmax in params_and_constants
    max_elev_id = np.max(pslev)
    elev_levels = np.mod(max_elev_id, max_number_elevs)
    # We've only pulled off the elevated tiles of these fields
    elev_start = 1

    # Make masks (with and without a tile dimensions)  of where the pre-defined ice 
    # gridboxes are in the UM grid. Checking for vol_smc_sat=0. is a common approach.
    soil_vsat_3d = shift_tile_dim(np.tile(soil_vsat_2d, [ntiles, 1, 1]))
    icemask2d = np.where(soil_vsat_2d[:, :] == 0.)
    icemask3d_alltype = np.where(soil_vsat_3d[:, :, :] == 0.)

    # The JULES library actually expects *everything* in on ntiles - need to
    # expand field with the first nelevs (the ice_elev parts) zeroed
    nonice_snowdepth = np.zeros([nx, ny, ntiles])
    nonice_snowdepth[:, :, elev_start + elev_levels - 1:] \
        = nonice_snowdepth_ice

    smb_from_nisnow = np.zeros([nx, ny, ntiles])

    # Make fortran-ordered versions of the numpy arrays, with pointers for passing
    # IN
    dzsnow_f, dzsnow_fp = np_to_fortpoint(dzsnow)
    tile_frac_old_f, tile_frac_old_fp = np_to_fortpoint(tile_frac_old)

    tile_frac_f, tile_frac_fp = np_to_fortpoint(tile_frac)
    nonice_snowdepth_f, nonice_snowdepth_fp = np_to_fortpoint(nonice_snowdepth)

    # INOUT, on tiles
    nsnow_f, nsnow_fp = np_to_fortpoint(nsnow, dtype=c_int)
    snow_tile_f, snow_tile_fp = np_to_fortpoint(snow_tile)
    snowdepth_f, snowdepth_fp = np_to_fortpoint(snowdepth)
    rho_snow_grnd_f, rho_snow_grnd_fp = np_to_fortpoint(rho_snow_grnd)
    rgrain_f, rgrain_fp = np_to_fortpoint(rgrain)

    # INOUT, on snowlayers
    ds_f, ds_fp = np_to_fortpoint(ds)
    sliq_f, sliq_fp = np_to_fortpoint(sliq)
    sice_f, sice_fp = np_to_fortpoint(sice)
    rho_snow_f, rho_snow_fp = np_to_fortpoint(rho_snow)
    rgrainl_f, rgrainl_fp = np_to_fortpoint(rgrainl)
    tsnow_f, tsnow_fp = np_to_fortpoint(tsnow)

    # OUT, on tiles
    smb_from_nisnow_f, smb_from_nisnow_fp = np_to_fortpoint(smb_from_nisnow)

    print("calling adjustsnow fortran lib")
    dummy = adjust_fort(c_int(nx),
                        c_int(ny),
                        c_int(ntiles),
                        c_int(nsmax),
                        dzsnow_fp,
                        c_int(elev_start),
                        c_int(elev_levels),
                        tile_frac_fp,
                        tile_frac_old_fp,
                        nsnow_fp,
                        nonice_snowdepth_fp,
                        snow_tile_fp,
                        snowdepth_fp,
                        rho_snow_grnd_fp,
                        rgrain_fp,
                        sice_fp,
                        sliq_fp,
                        rgrainl_fp,
                        tsnow_fp,
                        rho_snow_fp,
                        ds_fp,
                        smb_from_nisnow_fp)
    print("back from fortran lib")

    # Convert the modified fortran-style arrays back to regular numpy ones
    # OUT, on tiles
    nsnow = fort_to_np(nsnow_f)
    snow_tile = fort_to_np(snow_tile_f)
    snowdepth = fort_to_np(snowdepth_f)
    rho_snow_grnd = fort_to_np(rho_snow_grnd_f)
    rgrain = fort_to_np(rgrain_f)
    smb_from_nisnow = fort_to_np(smb_from_nisnow_f)

    # smb_from_nisnow needs to be dimensioned (10,144,192) to go back out to
    # the next year's coupling
    smb_from_nisnow = np.rollaxis(smb_from_nisnow, 2, 0)
    smb_from_nisnow = smb_from_nisnow[0:number_elevs, :, :]

    # OUT, on snowlayers
    ds = fort_to_np(ds_f)
    sliq = fort_to_np(sliq_f)
    sice = fort_to_np(sice_f)
    rho_snow = fort_to_np(rho_snow_f)
    rgrainl = fort_to_np(rgrainl_f)
    tsnow = fort_to_np(tsnow_f)

    # Reconstruct gridbox mean snow
    tmp = np.empty_like(snow_tile)
    tmp[icemask3d_alltype] = snow_tile[icemask3d_alltype] * \
                             tile_frac[icemask3d_alltype]
    snicemass_gbm[icemask2d] = np.sum(tmp, axis=2)[icemask2d]

    #Put the altered fields back in the UM dump

    # put simple lon,lat fields
    um_dump = np_to_um_2d(um_dump, stashcode_snow_amount, snicemass_gbm)

    # get fields whose third dimension is ice and rock elev tiles
    um_dump = np_to_um_3d(um_dump, stashcode_nsnow_layrs_tiles,
                          nsnow, elev=True)
    um_dump = np_to_um_3d(um_dump, stashcode_snow_tile, snow_tile, elev=True)
    um_dump = np_to_um_3d(um_dump, stashcode_snowdep_grd_tile,
                          snowdepth, elev=True)
    um_dump = np_to_um_3d(um_dump, stashcode_snowpack_bk_dens,
                          rho_snow_grnd, elev=True)
    um_dump = np_to_um_3d(um_dump, stashcode_rgrain, rgrain, elev=True)

    # put fields with snowlayers, on ice and rock elevs
    um_dump = np_to_um_4d(um_dump, stashcode_snow_laythk_tiles, ds, elev=True)
    um_dump = np_to_um_4d(um_dump, stashcode_snow_ice_tile, sice, elev=True)
    um_dump = np_to_um_4d(um_dump, stashcode_snow_liq_tile, sliq, elev=True)
    um_dump = np_to_um_4d(um_dump, stashcode_snow_laydns_tiles,
                          rho_snow, elev=True)
    um_dump = np_to_um_4d(um_dump, stashcode_snow_grnsiz_tiles,
                          rgrainl, elev=True)
    um_dump = np_to_um_4d(um_dump, stashcode_snow_t_tile, tsnow, elev=True)

    return um_dump, smb_from_nisnow
