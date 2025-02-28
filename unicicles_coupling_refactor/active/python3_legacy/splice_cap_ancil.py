from params_and_constants import *
import sys
import mule
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


def do_splice(um_dump,cap_file, x_offset=0, y_offset=0, assig_coef=5e-5, h_coef=0.15, las_linear=True):

  import numpy as np
  from um_to_np_utility import um_to_np_2d, \
                               np_to_um_2d
  from cap_assig import assig_4p4 as assig

  orog_sd = um_to_np_2d(cap_file,stashcode_orog_var)

  orog_as, orog_h = assig(orog_sd, assig_coef=assig_coef, h_coef=h_coef, las_linear=las_linear)

  for stash in [stashcode_orog,           \
                stashcode_orog_var,       \
                stashcode_orog_x_grad,    \
                stashcode_orog_y_grad,    \
                stashcode_orog_gdxx,      \
                stashcode_orog_gdxy,      \
                stashcode_orog_gdyy,      \
                stashcode_sil_orog_rough, \
                stashcode_hlf_pk_to_trf_ht]:

      field_um      =um_to_np_2d(um_dump, stash)
      field_cap_reg =um_to_np_2d(cap_file,stash)

      if stash==stashcode_sil_orog_rough: field_cap_reg = orog_as
      if stash==stashcode_hlf_pk_to_trf_ht: field_cap_reg = orog_h

      #CAP fields may be on a regional subgrid
      ny_reg=np.shape(field_cap_reg)[0]
      nx_reg=np.shape(field_cap_reg)[1]

      field_cap_glob=field_um.copy()

      field_cap_glob[y_offset:y_offset+ny_reg,x_offset:x_offset+nx_reg] \
          = field_cap_reg

      um_dump=np_to_um_2d(um_dump, stash, field_cap_glob)

  return um_dump


def splice_cap_ancil(GrIS_filename, AIS_filename, global_filename_in):

    res=None
    found_file=False
    if GrIS_filename != None:
        in_GrIS=mule.load_umfile(GrIS_filename)
        ncols=in_GrIS.integer_constants.num_cols
        found_file=True
        if ncols == 122:
            res="N216"
        elif ncols == 54:
            res="N96"
        
    if AIS_filename != None:
        in_AIS=mule.load_umfile(AIS_filename)
        ncols=in_AIS.integer_constants.num_cols
        found_file=True
        if ncols == 432:
            res="N216"
        elif ncols == 192:
            res="N96"

    if not(found_file):
        print("Need either a valid UM Greenland or an Antarctic file")
        exit(2)

    if res == "N216":
            assig_coeff=1e-4
            h_coef=0.3
            x_offset_gris=309
            y_offset_gris=266
            x_offset_ais=0
            y_offset_ais=1
    elif res == "N96":
            assig_coeff=5e-5
            h_coef=0.15
            x_offset_gris=137
            y_offset_gris=118
            x_offset_ais=0
            y_offset_ais=1
    else: 
        print("The resolution/size of the input files has not been recognised")
        exit(2)

    in_global=mule.load_umfile(global_filename_in)

    if GrIS_filename != None:
      final = do_splice(in_global,in_GrIS,x_offset=x_offset_gris,y_offset=y_offset_gris,assig_coef=assig_coeff,h_coef=h_coef,las_linear=True)

    if AIS_filename != None:
      if GrIS_filename == None:
        global_background = in_global
      else:
        global_background = final

      final = do_splice(global_background,in_AIS,x_offset=x_offset_ais,y_offset=y_offset_ais,assig_coef=assig_coeff,h_coef=h_coef,las_linear=True)

    #final.level_dependent_constants=mule.ff.FF_LevelDependentConstants(np.zeros([final.integer_constants.num_levels,4]))

    return final

if __name__ == "__main__":

 GrIS_filename, AIS_filename, global_filename_in, global_filename_out = parse_commandline(sys.argv)
 
 final = splice_cap_ancil(GrIS_filename, AIS_filename, global_filename_in)
 
 final.to_file(global_filename_out)
