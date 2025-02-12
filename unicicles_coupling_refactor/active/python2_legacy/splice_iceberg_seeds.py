import sys
import cf
from arg_to_file_exist import arg_to_file_exist


def parse_commandline(argv):
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument("--input_GrIS",    help="input, name of regional Greenland ancil")
  parser.add_argument("--input_AIS",     help="input, name of regional Antarctic ancil")
  parser.add_argument("--input_global",  help="input, name of background global ancil")
  parser.add_argument("--output_global", help="output, name of modified global ancil")
  args = parser.parse_args()

  err = 0

  GrIS_filename, err    = arg_to_file_exist(args.input_GrIS, mandatory=False, err=err)
  AIS_filename, err     = arg_to_file_exist(args.input_AIS, mandatory=False, err=err)
  global_filename_in, err  = arg_to_file_exist(args.input_global, err=err)
  global_filename_out, err = arg_to_file_exist(args.output_global,io="out", err=err)

  if err > 0:
    parser.print_help()
    sys.exit(2)

  return GrIS_filename, AIS_filename, global_filename_in, global_filename_out

def do_splice(global_background,regional_field,region=None):

    import numpy as np

    global_array=global_background.array
    half_y = np.shape(global_array)[-2]/2

    if region == "GrIS":
        global_array[...,half_y:,:]=regional_field.array[...,half_y:,:]
    elif region == "AIS":
        global_array[...,:half_y,:]=regional_field.array[...,:half_y,:]
    else:
        print "unknown region specified:", region 
        exit()

    global_background.insert_data(cf.Data(global_array))

    return global_background

def splice_iceberg_seeds(GrIS_filename, AIS_filename, global_filename_in):

    found_file=False
    calving_string = 'ncvar%calvingmask'

    if GrIS_filename != None:
         in_GrIS=cf.read(GrIS_filename).select_field(calving_string)
         found_file=True
    
    if AIS_filename != None:
         in_AIS=cf.read(AIS_filename).select_field(calving_string)
         found_file=True
    
    if not(found_file):
         print "Need either a Greenland or an Antarctic file"
         exit()
    
    in_global=cf.read(global_filename_in).select_field(calving_string)
    
    if GrIS_filename != None:
          final = do_splice(in_global,in_GrIS, region="GrIS")

    if AIS_filename != None:
      if GrIS_filename == None:
          global_background = in_global
      else:
          global_background = final

      final = do_splice(global_background, in_AIS,region="AIS")

    return final

if __name__ == "__main__":

    GrIS_filename, AIS_filename, global_filename_in, global_filename_out = parse_commandline(sys.argv)

    final = splice_iceberg_seeds(GrIS_filename, AIS_filename, global_filename_in)

    cf.write(final,global_filename_out)
