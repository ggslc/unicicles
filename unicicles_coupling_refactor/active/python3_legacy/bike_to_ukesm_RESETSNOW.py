import sys
from mule import load_umfile
import numpy as np
from netCDF4 import Dataset
from save_originals import *
from splice_cap_orog import *
from splice_tile_fracs import *
from reset_icetile_snowpack import *
from splice_tsurf_elev import *
from conservation_of_ice import *
from adjust_snow  import *

def validate_warn(self,filename=None, warn=True):
  self.validate_errors(filename=filename, warn=True)
  print("(no more than one warning is noted)")
  return

def disable_mule_validators(umf):
  umf.validate_errors=umf.validate
  print("OVERWRITING MULE'S VALIDATION FUNCTION")
  print("(the original is saved as self.validate_errors)")
  funcType = type(umf.validate)
  print("VALIDATION ERRORS DOWNGRADED TO WARNINGS")
  umf.validate=funcType(validate_warn, umf)
  return umf

                         

def parse_commandline(argv):
  import os.path
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument("--um_input",      help="input, name of UM dump to be modified")
  parser.add_argument("--toice_input",   help="input, name of file going to Unicicles")
  parser.add_argument("--um_output",     help="output, name of modified UM dump")
  parser.add_argument("--coupling_period",  help="climate-ice coupling period in secs")
  args = parser.parse_args()

  good_files=0
  for file in [args.um_input,args.toice_input]:
    if file:
      if  os.path.isfile(file):
        good_files+=1
      else:
          print("")
          print("ERROR: specified file does not exist:",file)

  if good_files==2:
    um_input=args.um_input
    toice_input=args.toice_input
  else:
    print("")
    print("ERROR: problem with one or more input files. I want ALL of them")
    print("")
    parser.print_help()
    sys.exit(2)

  if args.um_output != None:
    um_output=args.um_output
  else:
    print("")
    print("ERROR: specifiy an output file")
    print("")
    parser.print_help()
    sys.exit(2)

  if args.coupling_period != None:
    coupling_period=np.float(args.coupling_period)
  else:
    print("")
    print("ERROR: specifiy a climate-ice coupling period")
    print("")
    parser.print_help()
    sys.exit(2)

  return um_input,toice_input,um_output,coupling_period

##process filenames etc from the command line
if __name__ == "__main__":
  um_input, toice_input, um_output, coupling_period = parse_commandline(sys.argv)

##read in the UM dump to be modified
  print("Reading UM input ", um_input)
  ##beware, some version of MULE don't like the extra diagnostics in a full UKESM dump - two landmask entries?!
  ##running model needs all 152 fields though, not just 153 prognostics, or won't restart. This is fixed in MULE now?
  umf=load_umfile(um_input)
  um_dump=umf.copy()
  um_dump.fields=[]
  um_dump = disable_mule_validators(um_dump)
  for i in range(umf.fixed_length_header.raw[152]):
        um_dump.fields.append(umf.fields[i])

##read in the wrapper-processed output from unicicles
  print("Reading toice input ", toice_input)
  toice_file = Dataset(toice_input,'r')

##4. elev_ice snowpacks
  print("Resetting ice tile snow mass")
  um_dump=reset_icetile_snowpack(um_dump,toice_file,coupling_period,real_tile_only=True)

#output the modified dump
  print("Writing UM output")
  um_dump.to_file(um_output)
