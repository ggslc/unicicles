from arg_to_file_exist import arg_to_file_exist
import numpy as np
import os

def parse_commandline(argv):
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument("--input_GrIS",    help="input, name of regional Greenland ancil")
  parser.add_argument("--input_AIS",     help="input, name of regional Antarctic ancil")
  parser.add_argument("--output_global", help="output, name of modified global ancil")
  args = parser.parse_args()

  err = 0

  GrIS_filename, err    = arg_to_file_exist(args.input_GrIS, mandatory=False, err=err)
  AIS_filename, err     = arg_to_file_exist(args.input_AIS, mandatory=False, err=err)
  global_filename_out, err = arg_to_file_exist(args.output_global,io="out", err=err)

  if err > 0:
    parser.print_help()
    sys.exit(2)

  return GrIS_filename, AIS_filename, global_filename_out

def splice_icecouple(GrIS_filename, AIS_filename, global_filename_out):
    from netCDF4 import Dataset

    found_file = 0
    fields_to_splice = ["surface_elevation","cell_calving_flux","tile_land_fraction","tile_ice_fraction","nonice_snowdepth","ice_stemp"]

    if GrIS_filename != None: 
      gris_file = Dataset(GrIS_filename)
      found_file += 1
    
    if AIS_filename != None:
      ais_file = Dataset(AIS_filename)
      found_file += 2

    if found_file == 0: 
      print "must specify at least one of GrIS, AIS i.icecouple files", GrIS_filename, AIS_filename
      exit()

    if found_file == 1:  os.system("cp "+GrIS_filename+" "+global_filename_out)

    if found_file == 2:  os.system("cp "+AIS_filename+"  "+global_filename_out)

    if found_file == 3:  
      out_file=Dataset(global_filename_out,'w')

#echo global dimensions
      for dimname in gris_file.dimensions:
        dim = gris_file.dimensions[dimname]
        out_file.createDimension(dimname,len(dim))

      coverage_gris = np.sum(gris_file.variables['tile_land_fraction'][:],axis=0)
      coverage_ais  = np.sum(ais_file.variables['tile_land_fraction'][:],axis=0)
      valid = np.where( (coverage_gris > 0) | (coverage_ais > 0) )

#echo individual variables, with attributes
      for varname in gris_file.variables:
        var      = gris_file.variables[varname]
        var_ais  = ais_file.variables[varname]

        if varname == "total_ice_volume":
          outfield=out_file.createVariable(varname,var.dtype,var.dimensions)
          outfield[:] = var[:] + var_ais[:]
          outfield=out_file.createVariable('total_ice_volume_GrIS',var.dtype,var.dimensions)
          outfield[:] = var[:]
          outfield=out_file.createVariable('total_ice_volume_AIS',var.dtype,var.dimensions)
          outfield[:] = var_ais[:]

        elif varname == "change_in_ice_volume":
          outfield=out_file.createVariable(varname,var.dtype,var.dimensions)
          outfield[:] = var[:] + var_ais[:]
          outfield=out_file.createVariable('change_in_ice_volume_GrIS',var.dtype,var.dimensions)
          outfield[:] = var[:]
          outfield=out_file.createVariable('change_in_ice_volume_AIS',var.dtype,var.dimensions)
          outfield[:] = var_ais[:]

        elif varname in fields_to_splice:
          outfield=out_file.createVariable(varname,var.dtype,var.dimensions)
          merge = np.zeros_like(var[:])

          l_tile = len(np.shape(var[:])) == 3
          if l_tile: 
            for tile in range(np.shape(var[:])[0]):
              var_2d     = var[:][tile,:,:]
              var_ais_2d = var_ais[:][tile,:,:]
          
              merge_2d = np.zeros_like(var_2d)
              merge_2d[valid] = ( (var_2d[valid]*coverage_gris[valid]) + (var_ais_2d[valid]*coverage_ais[valid]) )/(coverage_gris[valid]+coverage_ais[valid])
              merge[tile,:,:] = merge_2d

          else:
            merge[valid] = ( (var[:][valid]*coverage_gris[valid]) + (var_ais[:][valid]*coverage_ais[valid]) )/(coverage_gris[valid]+coverage_ais[valid])

          outfield[:] = merge
        else:
          outfield=out_file.createVariable(varname,var.dtype,var.dimensions)
          outfield[:] = var

      out_file.close()

if __name__ == "__main__":
    import sys

    GrIS_filename, AIS_filename, global_filename_out = parse_commandline(sys.argv)

    splice_icecouple(GrIS_filename, AIS_filename, global_filename_out)
