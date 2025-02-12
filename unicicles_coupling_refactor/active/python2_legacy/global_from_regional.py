import sys
import cf
import mule
from splice_cap_ancil import splice_cap_ancil
from splice_iceberg_seeds import splice_iceberg_seeds
from splice_icecouple import splice_icecouple
from arg_to_file_exist import arg_to_file_exist

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
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument("--input_berg_GrIS",    help="input, name of regional Greenland ancil")
  parser.add_argument("--input_berg_AIS",     help="input, name of regional Antarctic ancil")
  parser.add_argument("--input_berg_global",  help="input, name of background global ancil")
  parser.add_argument("--output_berg_global", help="output, name of modified global ancil")
  parser.add_argument("--input_orog_GrIS",    help="input, name of regional Greenland ancil")
  parser.add_argument("--input_orog_AIS",     help="input, name of regional Antarctic ancil")
  parser.add_argument("--input_orog_global",  help="input, name of background global ancil")
  parser.add_argument("--output_orog_global", help="output, name of modified global ancil")
  parser.add_argument("--input_icecouple_GrIS",    help="input, name of regional Greenland icecouple file")
  parser.add_argument("--input_icecouple_AIS",     help="input, name of regional Antarctic icecouple file")
  parser.add_argument("--output_icecouple_global_out", help="output, name of modified global icecouple file")
  args = parser.parse_args()

  err = 0

  GrIS_berg_filename,err       = arg_to_file_exist(args.input_berg_GrIS, mandatory=False, err=err)
  AIS_berg_filename,err        = arg_to_file_exist(args.input_berg_AIS, mandatory=False, err=err)
  global_berg_filename_in,err  = arg_to_file_exist(args.input_berg_global, err=err)
  global_berg_filename_out,err = arg_to_file_exist(args.output_berg_global,io="out", err=err)
  GrIS_orog_filename,err       = arg_to_file_exist(args.input_orog_GrIS, mandatory=False, err=err)
  AIS_orog_filename,err        = arg_to_file_exist(args.input_orog_AIS, mandatory=False, err=err)
  global_orog_filename_in,err  = arg_to_file_exist(args.input_orog_global, err=err)
  global_orog_filename_out,err = arg_to_file_exist(args.output_orog_global,io="out", err=err)
  GrIS_icecouple_filename, err = arg_to_file_exist(args.input_icecouple_GrIS, mandatory=False)
  AIS_icecouple_filename, err  = arg_to_file_exist(args.input_icecouple_AIS, mandatory=False)
  global_icecouple_filename_out, err = arg_to_file_exist(args.output_icecouple_global_out,io="out")

  if err > 0:
    parser.print_help()
    sys.exit(2)

  return GrIS_berg_filename, AIS_berg_filename, global_berg_filename_in, global_berg_filename_out, GrIS_orog_filename, AIS_orog_filename, global_orog_filename_in, global_orog_filename_out, GrIS_icecouple_filename, AIS_icecouple_filename, global_icecouple_filename_out

if __name__ == "__main__":

    GrIS_berg_filename, AIS_berg_filename, global_berg_filename_in, global_berg_filename_out, GrIS_orog_filename, AIS_orog_filename, global_orog_filename_in, global_orog_filename_out, GrIS_icecouple_filename, AIS_icecouple_filename, global_icecouple_filename_out = parse_commandline(sys.argv)

    final_berg = splice_iceberg_seeds(GrIS_berg_filename, AIS_berg_filename, global_berg_filename_in)

    cf.write(final_berg,global_berg_filename_out)

    final_orog = splice_cap_ancil(GrIS_orog_filename, AIS_orog_filename, global_orog_filename_in)

    final_orog = disable_mule_validators(final_orog)

    final_orog.to_file(global_orog_filename_out)

    splice_icecouple(GrIS_icecouple_filename, AIS_icecouple_filename, global_icecouple_filename_out)
