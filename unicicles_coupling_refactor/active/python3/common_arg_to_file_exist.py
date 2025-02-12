"""
common_arg_to_file_exist

subroutine used by all scripts to check whether field specified on the
commandline for input or output really/already exist
"""

import os

def arg_to_file_exist(arg, mandatory=True, io="in", err=0):
    """
    Check whether field specified on command line for input or output exists
    """
    filename = None
    if arg is None:
        if mandatory:
            print("ERROR: mandatory argument missing")
            err = 1
    else:
        filename = arg
        exists = os.path.isfile(filename)

        if (exists) & (io == "out"):
            print("WARNING: output file already exists ", filename)
        if (not exists) & (io == "in"):
            if mandatory:
                print("ERROR: input file does not exist ", filename)
                err = 3
            else:
                print("WARNING: input file does not exist ", filename)
                filename = None

    return filename, err
