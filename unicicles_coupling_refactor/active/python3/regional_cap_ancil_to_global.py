"""
regional_cap_ancil_to_global

Subroutine for merging topography information calculated by CAP from the
regional ice sheet domains into a unified global field

Used by unicicles_cap_global_to_um
"""

import sys
import mule
import numpy as np
from common_params_and_constants import assig_n216, assig_n96, h_n216, h_n96, \
    las_a, las_bigsig, las_cd, las_k, \
    ncols_AIS_n216, ncols_AIS_n96, ncols_GrIS_n216, ncols_GrIS_n96, \
    stashcode_hlf_pk_to_trf_ht, stashcode_orog, stashcode_orog_gdxx, \
    stashcode_orog_gdxy, stashcode_orog_gdyy, stashcode_orog_x_grad, \
    stashcode_orog_y_grad, stashcode_orog_var, stashcode_sil_orog_rough, \
    xoff_AIS_n216, xoff_AIS_n96, xoff_GrIS_n216, xoff_GrIS_n96, \
    yoff_AIS_n216, yoff_AIS_n96, yoff_GrIS_n216, yoff_GrIS_n96
from common_arg_to_file_exist import arg_to_file_exist


def parse_commandline():
    """ Read the line command line arguments """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_GrIS",
                        help="input, name of regional Greenland ancil")
    parser.add_argument("--input_AIS",
                        help="input, name of regional Antarctic ancil")
    parser.add_argument("--input_global",
                        help="input, name of background global ancil")
    parser.add_argument("--output_global",
                        help="output, name of modified global ancil")
    args = parser.parse_args()

    err = 0

    GrIS_filename, err = arg_to_file_exist(args.input_GrIS, mandatory=False,
                                           err=err)
    if err != 0:
        parser.print_help()
        sys.exit(err)
    AIS_filename, err = arg_to_file_exist(args.input_AIS, mandatory=False,
                                          err=err)
    if err != 0:
        parser.print_help()
        sys.exit(err)
    global_filename_in, err = arg_to_file_exist(args.input_global, err=err)
    if err != 0:
        parser.print_help()
        sys.exit(err)
    global_filename_out, err = arg_to_file_exist(args.output_global, io="out",
                                                 err=err)
    if err != 0:
        parser.print_help()
        sys.exit(err)

    return GrIS_filename, AIS_filename, global_filename_in, global_filename_out


def do_splice(um_dump, cap_file, x_offset=0, y_offset=0, assig_coef=assig_n96,
              h_coef=h_n96, las_linear=True):

    """
    put the values from the regional domain CAP-generated ancil into the global
    domain fields in the UM dump
    """

    from common_um_to_np import np_to_um_2d, um_to_np_2d

    # Silhouette orographic roughness and the half-peak-to-trough aren't
    # calculated direct from the topography by the CAP anymore, so our
    # coupling doesn't have updated info for these. Use the older generation 
    # CAP subroutine to derive something usable from the standard deviation
    orog_sd = um_to_np_2d(cap_file, stashcode_orog_var)

    orog_as, orog_h = assig_4p4(orog_sd, assig_coef=assig_coef,
                                h_coef=h_coef, las_linear=las_linear)

    # loop over the orography fields and splice each in turn
    for stash in [stashcode_orog,
                  stashcode_orog_var,
                  stashcode_orog_x_grad,
                  stashcode_orog_y_grad,
                  stashcode_orog_gdxx,
                  stashcode_orog_gdxy,
                  stashcode_orog_gdyy,
                  stashcode_sil_orog_rough,
                  stashcode_hlf_pk_to_trf_ht]:

        # global field from the UM
        field_um = um_to_np_2d(um_dump, stash)

        # regional field from the CAP calculation
        field_cap_reg = um_to_np_2d(cap_file, stash)

        # these two fields don't have info in the CAP-derived file, we 
        # worked out some values ourselves in the subroutien call above.
        if stash == stashcode_sil_orog_rough:
            field_cap_reg = orog_as
        if stash == stashcode_hlf_pk_to_trf_ht:
            field_cap_reg = orog_h

        # CAP fields may be on a regional subgrid
        ny_reg = np.shape(field_cap_reg)[0]
        nx_reg = np.shape(field_cap_reg)[1]

        field_cap_glob = field_um.copy()

        # use the known offsets to put the regional field into the right place
        # in the global domain
        field_cap_glob[y_offset:y_offset + ny_reg,
                       x_offset:x_offset + nx_reg] = field_cap_reg

        um_dump = np_to_um_2d(um_dump, stash, field_cap_glob)

    return um_dump


def splice_cap_ancil(GrIS_filename, AIS_filename, global_filename_in):

    res = None
    found_file = False

    # check for existence and resolution of input for Greenland domain
    if GrIS_filename is not None:
        in_GrIS = mule.load_umfile(GrIS_filename)
        ncols = in_GrIS.integer_constants.num_cols
        found_file = True
        if ncols == ncols_GrIS_n216:
            res = "N216"
        elif ncols == ncols_GrIS_n96:
            res = "N96"

    # check for existence and resolution of input for Antarctic domain
    if AIS_filename is not None:
        in_AIS = mule.load_umfile(AIS_filename)
        ncols = in_AIS.integer_constants.num_cols
        found_file = True
        if ncols == ncols_AIS_n216:
            res = "N216"
        elif ncols == ncols_AIS_n96:
            res = "N96"

    if not found_file:
        print("Need either a valid UM Greenland or an Antarctic file")
        sys.exit(2)

    # load resolution dependent parameter sets
    if res == "N216":
        assig_coeff = assig_n216
        h_coef = h_n216
        x_offset_gris = xoff_GrIS_n216
        y_offset_gris = yoff_GrIS_n216
        x_offset_ais = xoff_AIS_n216
        y_offset_ais = yoff_AIS_n216
    elif res == "N96":
        assig_coeff = assig_n96
        h_coef = h_n96
        x_offset_gris = xoff_GrIS_n96
        y_offset_gris = yoff_GrIS_n96
        x_offset_ais = xoff_AIS_n96
        y_offset_ais = yoff_AIS_n96
    else:
        print("The resolution/size of the input files has not been recognised")
        sys.exit(2)

    # load the UM file with the global fields where going to modify
    in_global = mule.load_umfile(global_filename_in)

    # if there's Greenland information, use it to modify the global field
    if GrIS_filename is not None:
        final = do_splice(in_global, in_GrIS, x_offset=x_offset_gris,
                          y_offset=y_offset_gris, assig_coef=assig_coeff,
                          h_coef=h_coef, las_linear=True)

    # if there's Antartic information, choose which are the relevant global
    # fields to modify and then modify them
    if AIS_filename is not None:
        if GrIS_filename is None:
            global_background = in_global
        else:
            global_background = final

        final = do_splice(global_background, in_AIS, x_offset=x_offset_ais,
                          y_offset=y_offset_ais, assig_coef=assig_coeff,
                          h_coef=h_coef, las_linear=True)

    # some versions of mule won't write accept ancil-format fields that are
    # missing this header 
    # final.level_dependent_constants \
    #    = mule.ff.FF_LevelDependentConstants(
    #      np.zeros([final.integer_constants.num_levels,4]))

    return final


def assig_4p4(orog_sd, assig_coef=assig_n96, h_coef=h_n96, las_linear=True):
    """

    # Subroutine to calculate the A/S parameter from sigma-h

    # Based on the pre-GLOBE CAP routine (here vn4.4)
    # We also return the half-peak-to_trough field - the CAP just set
    # this to the standard deviation, optionally scaled by a namelist
    # parameter.
    # Eyeballing G'land and Ant cf N96 orig ancil and compromising,
    # simple linear scalings for the earliest UKESM1-IS were as
    # the defaults above. Doesn't get the highs, overestimates the lows
    # In general, factors are too small for Gland and too high for Ant
    ### N216 eyeballed: assig_coef=1e-4 h_coef=0.3###


    *from UM CAP@vn4.4*
    *
    *
    * There are two possible methods
    *
    * 1) A simple linear relationship (LAS_LINEAR)
    *
    *  A/S =a * sigma-h
    *
    * where a is a coefficent (argument ASSIG_COEF)
    *
    * a is model grid dependent and can be calculated externally
    * and fed to the program using the namelist or calculated
    * internally
    #
    # we can't do this internal calculation - relies on there being
    # pre-existing regions where we know both orog_sd and orog_as
    # and we don't have those here
    *
    *  Standard values are  (10' data)          (5'data)
    *  oper glob            a=0.0001679         a=0.000156687
    *  oper LAM             a=0.0002986         a=0.00027718
    *  climate              a=0.000084454       a=0.0000743820
    # in vn4.5 UM-Glimmer coupling, we used
    #  HadCM3               a=8.4e-5
    #  FAMOUS               a=5.9e-5
    *
    * 2) A more complicated method
    *
    *   A/S=  2k**2/Cd/(ln(zc/zo))**2
    *
    *  k=von Karman's constant (0.4)
    *  Cd= drag coeffiecent =0.3
    *  zc=SQRT(2)SIGMA
    *  zo=aSIGMA**2   for sigma <= BIGSIG
    *    =aBIGSIG(2SIGMA-BIGSIG) otherwise
    * where
    *  BIGSIG=100    a=5E-4

    """

    valid = np.where(orog_sd > 0.)

    orog_h = np.zeros_like(orog_sd)
    orog_h[valid] = orog_sd[valid] * h_coef

    orog_as = np.zeros_like(orog_sd)

    if las_linear:

        orog_as[valid] = orog_sd[valid] * assig_coef

        orog_as = np.maximum(orog_as, 0.)
        orog_as = np.minimum(orog_as, 0.2)

    else:

        k = las_k
        cd = las_cd
        a = las_a
        bigsig = las_bigsig

        big = np.where(orog_sd > bigsig)

        zc = np.sqrt(2.0) * orog_sd
        z0 = a * (orog_sd)**2.0
        z0[big] = a * bigsig * (2.0 * orog_sd[big] - bigsig)

        as1 = 2.0 * k * k / cd
        as2 = (np.log(zc[valid] / z0[valid]))**2.0

        orog_as[valid] = as1 / as2

    return orog_as, orog_h


if __name__ == "__main__":

    # Read the command line arguments
    GrIS_filename, AIS_filename, global_filename_in, global_filename_out \
        = parse_commandline()

    # put the regional ancil info into the global UM fields
    final = splice_cap_ancil(GrIS_filename, AIS_filename, global_filename_in)

    # write out the final results
    final.to_file(global_filename_out)
