"""
regional_iceberg_seeds_to_global

Subroutine for merging iceberg discharge information from BISICLES
regional ice sheet domains into a unified global field of iceberg seeding

Used by bisicles_global_to_nemo
"""

import sys
import cf
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
    AIS_filename, err = arg_to_file_exist(args.input_AIS, mandatory=False,
                                          err=err)
    global_filename_in, err = arg_to_file_exist(args.input_global, err=err)
    global_filename_out, err = arg_to_file_exist(args.output_global,
                                                 io="out", err=err)

    if err > 0:
        parser.print_help()
        sys.exit(err)

    return GrIS_filename, AIS_filename, global_filename_in, global_filename_out


def do_splice(global_background, regional_field, region=None):

    """
    put the berg array for either the northern or southern hemisphere into
    the global CF Field
    """

    import numpy as np

    global_array = global_background.array

    # find the index of equatorial latitude so we can split the world
    # into N and S hemispheres
    half_y = int(np.shape(global_array)[-2] / 2)

    # avoid overwriting changes that may have been made by the other ice sheet
    # by only updating the NH part of the array if we have GrIS input, and only
    # the SH part if we have AIS
    if region == "GrIS":
        global_array[..., half_y:, :] = regional_field.array[..., half_y:, :]
    elif region == "AIS":
        global_array[..., :half_y, :] = regional_field.array[..., :half_y, :]
    else:
        print("unknown region specified:", region)
        sys.exit(1)

    # put the array into the CF Field
    global_background.set_data(cf.Data(global_array))

    return global_background


def splice_iceberg_seeds(GrIS_filename, AIS_filename, global_filename_in):

    """
    put the berg info for the northern and/or southern hemispheres into
    the global CF Field
    """

    found_file = False
    calving_string = 'ncvar%calvingmask'


    # we might have been given info for either GrIS, AIS or both

    if GrIS_filename is not None:
        in_GrIS = cf.read(GrIS_filename).select_field(calving_string)
        found_file = True

    if AIS_filename is not None:
        in_AIS = cf.read(AIS_filename).select_field(calving_string)
        found_file = True

    if not found_file:
        print("Need either a Greenland or an Antarctic file")
        sys.exit(1)

    # read in the template calving field. It already contains a calving
    # pattern for both NH and SH, so we'll just be overwriting a hemisphere
    # if we have information from the ice sheet for it
    in_global = cf.read(global_filename_in).select_field(calving_string)

    # if there's Greenland information, use it to modify the global field
    if GrIS_filename is not None:
        final = do_splice(in_global, in_GrIS, region="GrIS")
    
    # if there's Antartic information, choose whether to modify the original
    # template or the version modified by GrIS, if it was done above
    if AIS_filename is not None:
        if GrIS_filename is None:
            global_background = in_global
        else:
            global_background = final

        final = do_splice(global_background, in_AIS, region="AIS")

    return final


if __name__ == "__main__":

    # Read the command line arguments
    GrIS_filename, AIS_filename, global_filename_in, global_filename_out \
        = parse_commandline()
    
    # put the regional iceberg info into the global NEMO field
    final = splice_iceberg_seeds(GrIS_filename, AIS_filename,
                                 global_filename_in)

    # write out the final results
    cf.write(final, global_filename_out)
