"""
cap_orog_to_um

subroutine used by unicicles_cap_global_to_um to insert the modified surface
topography fields into the UM restart dump

Called by unicicles_cap_global_to_um
"""

import numpy as np
from common_params_and_constants import stashcode_hlf_pk_to_trf_ht, \
    stashcode_land_frac, stashcode_orog, stashcode_orog_gdxx, \
    stashcode_orog_gdxy, stashcode_orog_gdyy, stashcode_orog_var, \
    stashcode_orog_x_grad, stashcode_orog_y_grad, stashcode_sil_orog_rough
from common_um_to_np import get_ice_um_grid_overlap, ice_to_np, \
    merge_ice_um_arrays, np_to_um_2d, um_to_np_2d

def insert_cap_ice_orog(um_dump, um_ref_dump, cap_file, ice_file):
    """
    Insert the modified surface topography fields into the UM restart dump.
    """
    # def splice_cap_orog(um_dump, um_ref_dump, cap_file, ice_file,
    #                    x_offset=0, y_offset=0):

    # Pull UM coastal fraction
    coast_um = um_to_np_2d(um_ref_dump, stashcode_land_frac)

    # Get overlap between UM and ice grids
    fmask_ice = get_ice_um_grid_overlap(ice_file)

    # Get BISICLES orog on the UM grid - this is the best thing to use as
    # a mask for where the ISM domain has real data worth adding
    orog_ism = ice_to_np(ice_file, 'surface_elevation')

    # Update where the ice grid has data and the UM mask says there's some
    # land we need to include the coastal fraction weighting in here too?
    ice_updates = np.where((orog_ism > 0.) & (coast_um > 0))

    for stash in [stashcode_orog,
                  stashcode_orog_var,
                  stashcode_orog_x_grad,
                  stashcode_orog_y_grad,
                  stashcode_orog_gdxx,
                  stashcode_orog_gdxy,
                  stashcode_orog_gdyy,
                  stashcode_sil_orog_rough,
                  stashcode_hlf_pk_to_trf_ht]:

        #get fields for each stashcode from the UM reference source 
        #(what we will change) and from the global CAP ancil 
        #(what we will change it with)
        field_um_ref = um_to_np_2d(um_ref_dump, stash)
        field_cap_glob = um_to_np_2d(cap_file, stash)

        # Splice in the CAP-processed fields, weighted by the UM/ice grid
        # overlap and the UM land fraction
        field_spliced = merge_ice_um_arrays(field_um_ref,
                                            field_cap_glob,
                                            fmask_ice * coast_um,
                                            points=ice_updates,
                                            noFrac=True)

        #put the altered field back into the UM dump.
        um_dump = np_to_um_2d(um_dump, stash, field_spliced)

    return um_dump
