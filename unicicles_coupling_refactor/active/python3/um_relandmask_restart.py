"""
um_relandmask_restart

EXPERIMENTAL: part of being able to change landsea mask in the ice coupling
workflow. Changes the land sea mask in a UM dump and then resizes *all* fields 
that care to match the new definition

Called by unicicles_cap_global_to_um
"""

import numpy as np
import common_mule_rss as mule_rss
from common_params_and_constants import stashcode_land_frac, \
    stashcode_land_frac_diag, stashcode_lsm

def remask_var(array, missing_flag, newlsm):

    """
    simple routine for filling out any new land areas with their nearest (in array space, not lon-lat)
    pre-existing neighbour
    """

    from scipy import ndimage

    # http://stackoverflow.com/questions/3662361/fill-in-missing-values-with-nearest-neighbour-in-python-numpy-masked-arrays
    array_masked = np.ma.masked_equal(array, missing_flag)

    sea_pts = np.where(newlsm < -0.5)
    # Print len(sea_pts[0])," sea points in new mask"
    ind = ndimage.distance_transform_edt(array_masked.mask,
                                         return_distances=False,
                                         return_indices=True)

    a = array[tuple(ind)]

    if len(sea_pts[0]) > 0:
        a[sea_pts] = missing_flag

    return a


def new_lsm(umf, lsm_file):

    # apply mask from lsm_file netCDF to an already open UM dump

    #lists of mask stashcodes in the dump
    lsm_stash_bin = [stashcode_lsm]
    lsm_stash_frac = [stashcode_land_frac, stashcode_land_frac_diag]

    #read in the new (fractional?) land sea mask, make an index of points with
    #enough area to be considered land
    fland = lsm_file.variables['lsmask'][:].squeeze()
    fland = np.ma.filled(fland, fill_value=0.)
    minval = .01
    new_landindex = np.where(fland > minval)

    #grab (first instance) of the UM binary land sea mask
    index_LSM = mule_rss.findindex(umf, lsm_stash_bin[0])[0]
    LSM_old = mule_rss.get_data2d(umf, index_LSM)

    #grab (first instance) of the UM fractional land sea mask
    index_fLSM = mule_rss.findindex(umf, lsm_stash_frac[0])[0]
    fLSM_old = mule_rss.get_data2d(umf, index_fLSM)

    #make new arrays with the new mask values
    LSM_new = LSM_old.copy()
    LSM_new[:, :] = 0.
    LSM_new[new_landindex] = 1

    fLSM_new = fLSM_old.copy()
    fLSM_new[:, :] = 0.
    fLSM_new[new_landindex] = fland[new_landindex]


    # To resize all the land fields in the dump we simply need to loop 
    # over them, reading them in from the dump and then putting them back
    # once the new binary mask is in place
    last_stash = -99
    count = 0

    for index in range(len(umf.fields)):

        # Check stash for this field in the loop
        this_stash = umf.fields[index].lbuser4
        missing_flag = umf.fields[index].bmdi

        # Handy to know which level we're on
        if this_stash == last_stash:
            count = count + 1
        if this_stash > last_stash:
            count = 1

        last_stash = this_stash

        # Read field from dump
        field = mule_rss.get_data2d(umf, index)

        if this_stash in lsm_stash_bin:
            print("putting binary mask direct into ", this_stash)
            field = LSM_new

        # My UM4.5 version did specific things to sea-ice too so there were
        # valid SSTs/albedo where we made new ocean. Essential?
        # open sea sst
        # if this_stash == 507: field = ?

        missing = np.where(field == missing_flag)
        valid = np.where(field != missing_flag)

        # Identify land fields that need resizing
        # Seems safe to assume all fields that appear with missing data and are the same
        # size as the global mask array are land and can be masked? Incl diagnostics
        if (len(missing[0]) > 0) & (len(field) == len(LSM_new)) > 0:

            if this_stash in lsm_stash_frac:
                #treat the fractional mask specially
                field = fLSM_new
                print("Putting fractional mask direct into ", this_stash)
            else:
                print("Simple near-neighbor fill for stash ", this_stash,
                      count, field.shape)
                field = remask_var(field, missing_flag, LSM_new)

        # write field back out again
        umf = mule_rss.put_data2d(umf, index, field)

    return umf
