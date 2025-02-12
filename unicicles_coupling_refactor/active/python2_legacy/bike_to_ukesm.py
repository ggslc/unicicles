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
from new_landseamask import *

def validate_warn(self,filename=None, warn=True):
  self.validate_errors(filename=filename, warn=True)
  print "(no more than one warning is noted)"
  return

def disable_mule_validators(umf):
  umf.validate_errors=umf.validate
  print "OVERWRITING MULE'S VALIDATION FUNCTION"
  print "(the original is saved as self.validate_errors)"
  funcType = type(umf.validate)
  print "VALIDATION ERRORS DOWNGRADED TO WARNINGS"
  umf.validate=funcType(validate_warn, umf)
  return umf

                         

def parse_commandline(argv):
  import os.path
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument("--um_input",      help="input, name of UM dump to be modified")
  parser.add_argument("--um_ref_input",  help="input, name of source of UM reference fields")
  parser.add_argument("--toice_input",   help="input, name of file going to Unicicles")
  parser.add_argument("--fromice_input", help="input, name of file output by Unicicles")
  parser.add_argument("--fromcap_input", help="input, name of file output by CAP")
  parser.add_argument("--um_output",     help="output, name of modified UM dump")
  parser.add_argument("--coupling_period",  help="climate-ice coupling period in secs")
  parser.add_argument("--tonextice_output", help="output, additonal contribution to icesheet SMB next time")
  parser.add_argument("--do_snow_shuffling", help="(re)enable non-ice snow shuffling, in coupling, experimental", action="store_true")
  parser.add_argument("--new_landseamask", help="is there a new ice shelf land sea mask?")
  args = parser.parse_args()

  good_files=0
  for file in [args.um_input,args.um_ref_input,args.toice_input,args.fromice_input, args.fromcap_input]:
    if file:
      if  os.path.isfile(file):
        good_files+=1
      else:
          print ""
          print "ERROR: specified file does not exist:",file

  if good_files==5:
    um_input=args.um_input
    um_ref_input=args.um_ref_input
    toice_input=args.toice_input
    fromice_input=args.fromice_input
    fromcap_input=args.fromcap_input
  else:
    print ""
    print "ERROR: problem with one or more input files. I want ALL of them"
    print ""
    parser.print_help()
    sys.exit(2)

  if args.um_output != None:
    um_output=args.um_output
  else:
    print ""
    print "ERROR: specify an output file"
    print ""
    parser.print_help()
    sys.exit(2)

  if args.tonextice_output != None:
    tonextice_output=args.tonextice_output
  else:
    print ""
    print "ERROR: specify an output file"
    print ""
    parser.print_help()
    sys.exit(2)

  if args.coupling_period != None:
    coupling_period=np.float(args.coupling_period)
  else:
    print ""
    print "ERROR: specify a climate-ice coupling period"
    print ""
    parser.print_help()
    sys.exit(2)

  do_snow_shuffling = False
  if args.do_snow_shuffling: do_snow_shuffling = True

  if args.new_landseamask != None:
    if  os.path.isfile(args.new_landseamask):
      new_landseamask = args.new_landseamask
    else:
      print ""
      print "ERROR: specified new land-sea mask file doesn't exist"
      print ""
      parser.print_help()
      sys.exit(2)
  else:
    new_landseamask = None


  return um_input,um_ref_input,toice_input,fromice_input,fromcap_input,um_output, tonextice_output, coupling_period, do_snow_shuffling, new_landseamask

##process filenames etc from the command line
if __name__ == "__main__":
  um_input, um_ref_input, toice_input, fromice_input, fromcap_input, um_output, tonextice_output, coupling_period, do_snow_shuffling, new_landseamask = parse_commandline(sys.argv)

##read in the UM dump to be modified
  print "Reading UM input ", um_input
  ##beware, some version of MULE don't like the extra diagnostics in a full UKESM dump - two landmask entries?!
  ##running model needs all 152 fields though, not just 153 prognostics, or won't restart. This is fixed in MULE now?
  umf=load_umfile(um_input)
  um_dump=umf.copy()
  um_dump.fields=[]
  um_dump = disable_mule_validators(um_dump)
  for i in range(umf.fixed_length_header.raw[152]):
        um_dump.fields.append(umf.fields[i])

##read in the UM reference fields
  print "Reading UM reference input ", um_ref_input
  umf=load_umfile(um_ref_input)
  um_ref_dump=umf.copy()
  um_ref_dump.fields=[]
  for i in range(umf.fixed_length_header.raw[152]):
        um_ref_dump.fields.append(umf.fields[i])

##read in the wrapper-processed output from unicicles
  print "Reading toice input ", toice_input
  toice_file = Dataset(toice_input,'r')

  print "Reading fromice input ", fromice_input
  fromice_file = Dataset(fromice_input,'r')

##read in the CAP-processed output from BISICLES
  print "Reading CAP input ",fromcap_input
  fromcap_file = load_umfile(fromcap_input)

## modify:
##0. save the originals of some fields for conservation calcs later
  frac_elev0,snicemass_elev0 = save_originals(um_dump)

##do the new landsea mask
  if new_landseamask != None:
    print "Putting in new land sea mask"
    print "Reading mask", new_landseamask
    lsm_file = Dataset(new_landseamask,'r')
    um_dump=new_lsm(um_dump,lsm_file)

##1. mean orography 
###just splice in the field from the wrapper
#TAKE FROM CAP OUTPUT INSTEAD => CAN HAVE FILTERING etc IF REQ'D
#  print "Handling mean orography"
#  um_dump=splice_mean_orog(um_dump,um_ref_dump,fromice_file)

##2. subgrid orog stats
###just splice in the fields from the CAP output
###Right now, CAP fields are on a regional subgrid. Some care needs to
###be taken around resolutions, and whether the ancil and dump count
###latitudes N-S or S-N
  print "Putting in CAP orog fields"
  um_dump=splice_cap_orog(um_dump,um_ref_dump,       \
                          fromcap_file,fromice_file)

##3. tile fractions
  print "Update tile fractions"
  um_dump=splice_tile_fracs(um_dump,um_ref_dump,fromice_file)

##4. elev_ice snowpacks
  print "Resetting ice tile snow mass"
  um_dump=reset_icetile_snowpack(um_dump,toice_file,coupling_period)

##5. elev_rock snowpacks
  if do_snow_shuffling == True:
    print "SNOW SHUFFLING ON"
    print "Adjusting the rock elev snow on the new fracs"
    um_dump, smb_from_nisnow = adjust_snow(um_dump,fromice_file,frac_elev0)

    if np.sum(smb_from_nisnow) > 0.0 :
      print "Write file for addition to next year's ice SMB"
      import cf
      f=cf.read(toice_input)[0]
      g=f.copy()
      g.setprop('long_name', 'smb_from_nisnow')
      g.id = 'smb_from_nisnow'
      for prop in ['stash_code','um_stash_source']:
        if g.hasprop(prop): g.delprop(prop)

      g.insert_data(cf.Data(smb_from_nisnow/coupling_period, units='kg/m2/s'))
      cf.write(g,tonextice_output,fmt='NETCDF4')

  else:
    print "SNOW SHUFFLING OFF"
    print "NOT ADJUSTING ELEV_ROCK SNOWPACKS"

##6. elev_tile subsurface temperatures
# BISICLES not using temps at the moment - just getting 0s back
  print "NOT updating sub-icesnowpack boundary conditions, BISICLES all 0"
#  print "Updating sub-icesnowpack boundary conditions"
#  um_dump=splice_tsurf_elev(um_dump,um_ref_dump,fromice_file)

##7. check conservation. Want change in ice volume, and calving too really
#  TOO FLAKY RIGHT NOW, NOT ATTEMPTING TO MONITOR/ALTER ONLINE
#
#  print "Assessing conservation"
#  conservation_of_ice(um_dump,toice_file,fromice_file, \
#                      frac_elev0,snicemass_elev0)

##8. put the Ice->ocean flux into NEMO dump
#    NOW DONE EARLIER VIA ICEBERG SEEDING FILE IN BIKE_regrid_to_hdf.py

#output the modified dump
  print "Writing UM output"
  um_dump.to_file(um_output)

