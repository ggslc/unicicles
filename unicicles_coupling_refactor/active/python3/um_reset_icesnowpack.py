"""
um_reset_icesnowpack

Reset the snowpack on the ice sheet tiles in UM/JULES by adding back in
the SMB that has been extracted from them over the last couping cycle, and
adjust all dependent UM snow prognostics to match

Called by unicicles_cap_global_to_um
"""

import sys
import numpy as np
from common_params_and_constants import max_number_elevs, nsmax, \
    rho_ice_glimmer, stashcode_frac_surf_type, stashcode_nsnow_layrs_tiles, \
    stashcode_snow_amount, stashcode_snow_ice_tile, \
    stashcode_snow_laydns_tiles, stashcode_snow_laythk_tiles, \
    stashcode_snow_tile, stashcode_snowdep_grd_tile, \
    stashcode_snowpack_bk_dens, stashcode_vol_smc_sat
from common_um_to_np import ice_to_np, np_to_um_2d, np_to_um_3d, np_to_um_4d, \
    shift_tile_dim, um_to_np_2d, um_to_np_3d, um_to_np_4d


def reset_icetile_snowpack(um_dump, ice_file, coupling_period_in_secs,
                           real_tile_only=False):
    """
    Reset the snowpack on the ice sheet tiles in UM/JULES by adding back in
    the SMB that has been extracted from them over the last couping cycle, and
    adjust all dependent UM snow prognostics to match
    """
    # get arrays for gridbox mean variables
    soil_vsat_2d = um_to_np_2d(um_dump, stashcode_vol_smc_sat)
    snicemass_gbm = um_to_np_2d(um_dump, stashcode_snow_amount)

    # get arrays for variables on elevated tiles
    nsnice_elev, pslev = um_to_np_3d(um_dump, stashcode_nsnow_layrs_tiles,
                                     elev=True)
    snicemass_elev, _ = um_to_np_3d(um_dump, stashcode_snow_tile, elev=True)
    snicedep_elev, _ = um_to_np_3d(um_dump, stashcode_snowdep_grd_tile,
                                   elev=True)
    snicerho_elev, _ = um_to_np_3d(um_dump, stashcode_snowpack_bk_dens,
                                   elev=True)
    frac_elev, _ = um_to_np_3d(um_dump, stashcode_frac_surf_type, elev=True)

    # get arrays for variables on individual snow levels of the elevated tiles
    snice_layer_mass, _ = um_to_np_4d(um_dump, stashcode_snow_ice_tile,
                                      elev=True)
    snice_layer_dep, _ = um_to_np_4d(um_dump, stashcode_snow_laythk_tiles,
                                     elev=True)
    snice_layer_den, _ = um_to_np_4d(um_dump, stashcode_snow_laydns_tiles,
                                     elev=True)

    # get the number of elevated tiles, the highest tile_id and the number of different heights
    nt = len(pslev)
    max_elev_id = np.max(pslev)
    nelev = np.mod(max_elev_id, max_number_elevs)

    # make global masks of which gridboxes and tiles are ice_or rock_ elevated, based on
    # them not having the regular soil property vol_smc_sat
    soil_vsat_3d = shift_tile_dim(np.tile(soil_vsat_2d, [nt, 1, 1]))
    # Pick out all ice gridboxes
    icemask2d = np.where(soil_vsat_2d[:, :] == 0.)
    # Pick out the *ice* tiles of the elevated ones - only 0:nelev are looked
    # at!
    icemask3d = np.where(soil_vsat_3d[:, :, 0:nelev] == 0.)
    icemask3d_alltype = np.where(soil_vsat_3d[:, :, 0:2 * nelev] == 0.)

    # usually we look at all possible elevated tiles, real_tile_only setting is for only 
    # considering tiles currently active in the prognostic model - see sanity check below
    if real_tile_only:
        icemask3d = np.where(frac_elev[:, :, 0:nelev] > 1e-3)

    # Need to convert smb kg/m2/s to the total mass that was given to the
    # icesheet for this cycle. No ns_flip becase the toice file is the same north-south
    # array ordering as the UM, unlike the fromice
    smb_ice = ice_to_np(ice_file, 'ice_smb', ns_flip=False) * \
              coupling_period_in_secs

    # Sanity check - do we have nsmax snow on all ice tiles? Ice tiles should *always* have a 
    # very large depth of snowpack or the coupling framework has broken somewhere
    if np.any(nsnice_elev[icemask3d] < nsmax):
        print("3d:nsnow is no longer 10 on all the ice tiles I'm looking at!")
        print("if you *know* this is only on virtual tiles and you don't ",
              "mind that,")
        print("Set real_tile_only=True calling this script to ignore ",
              "virtual tiles")
        print("If real_tile_only is already True and you're reading this, ",
              "you've got real problems")
        sys.exit(1)

    # Put the variables on individual snow layers back to their original value by adding /
    # /subtracting the relevant amount of mass as solid ice from the bottom layer
    snice_layer_mass[:, :, :, -1][icemask3d] \
        = snice_layer_mass[:, :, :, -1][icemask3d] - smb_ice[icemask3d]
    snice_layer_dep[:, :, :, -1][icemask3d] \
        = snice_layer_dep[:, :, :, -1][icemask3d] - \
        (smb_ice[icemask3d] / rho_ice_glimmer)
    snice_layer_den[:, :, :, -1][icemask3d] \
        = snice_layer_mass[:, :, :, -1][icemask3d] / \
        snice_layer_dep[:, :, :, -1][icemask3d]

    # Make the same adjustment to the variables for tile total mass, depth and density
    snicemass_elev[icemask3d] = snicemass_elev[icemask3d] - smb_ice[icemask3d]
    snicedep_elev[icemask3d] = snicedep_elev[icemask3d] - \
                               (smb_ice[icemask3d] / rho_ice_glimmer)
    snicerho_elev[icemask3d] = snicemass_elev[icemask3d] / \
                               snicedep_elev[icemask3d]

    # Do the gridbox mean snowpack mass, looking at snow on ice /and/ rock tiles, not
    # just icemask_3d.
    tmp = np.empty_like(snicemass_elev)
    tmp[icemask3d_alltype] = snicemass_elev[icemask3d_alltype] * \
                             frac_elev[icemask3d_alltype]
    snicemass_gbm[icemask2d] = np.sum(tmp, axis=2)[icemask2d]

    # put gridbox mean fields back into the UM dump
    um_dump = np_to_um_2d(um_dump, stashcode_snow_amount, snicemass_gbm)

    # put the elevated tile fields back into the UM dump
    um_dump = np_to_um_3d(um_dump, stashcode_snow_tile, snicemass_elev,
                          elev=True)
    um_dump = np_to_um_3d(um_dump, stashcode_snowdep_grd_tile,
                          snicedep_elev, elev=True)
    um_dump = np_to_um_3d(um_dump, stashcode_snowpack_bk_dens,
                          snicerho_elev, elev=True)

    # put the fields on individual snow layers  back into the UM dump
    um_dump = np_to_um_4d(um_dump, stashcode_snow_ice_tile,
                          snice_layer_mass, elev=True)
    um_dump = np_to_um_4d(um_dump, stashcode_snow_laythk_tiles,
                          snice_layer_dep, elev=True)
    um_dump = np_to_um_4d(um_dump, stashcode_snow_laydns_tiles,
                          snice_layer_den, elev=True)

    return um_dump
