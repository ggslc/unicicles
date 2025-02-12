import sys
import os
import numpy as np
import cf

def parse_commandline(argv):
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument("--input",      help="input, name of UM forcing data")
  parser.add_argument("--output",     help="output, name of input to uncicles")
  parser.add_argument("--gbm_smb",    help="aggregate individual tile SMB - eg lose height dependence")
  parser.add_argument("--gbm_lapse",  help="use Edwards et al '14 SMB lapse rate to fake height dependence")
  parser.add_argument("--exclude_code",  action="append",help="don't expect these fields in input - output them as 0")
  parser.add_argument("--do_snow_shuffling",  action="append",help="(re)enable non-ice snow shuffling, in coupling, experimental")
  args = parser.parse_args()

  good_files=0
  for file in [args.input]:
    if file:
      if  os.path.isfile(file):
        good_files += 1
      else:
          print ""
          print "ERROR: specified file does not exist:",file

  if good_files == 1:
    input=args.input
  else:
    print ""
    print "ERROR: problem with one or more input files. I want ALL of them"
    print ""
    parser.print_help()
    sys.exit(2)

  if args.output != None:
    output=args.output
  else:
    print ""
    print "ERROR: explicitly specifiy an output file"
    print ""
    parser.print_help()
    sys.exit(2)

  gbm_smb = False
  if args.gbm_smb != None: 
     if args.gbm_smb != "False": gbm_smb = True

  gbm_lapse = False
  if args.gbm_lapse != None: 
     if args.gbm_lapse != "False": gbm_lapse = True

  do_snow_shuffling = False
  if args.do_snow_shuffling != None: 
     if args.do_snow_shuffling != "False": do_snow_shuffling = True

  exclude_codes = []
  if type(args.exclude_code) == type([]): exclude_codes = list(map(int, args.exclude_code))

  return input,output,gbm_smb,gbm_lapse,exclude_codes, do_snow_shuffling

def extract_fields(all_fields,extract_codes,exclude_codes):

  extract_list=cf.FieldList()

  for code in extract_codes: 
    read_code = code

    print code, "in extract_fields, excluding", exclude_codes

    if code in exclude_codes:
      #for now, assume we're only faking tile fields
      print "faking a tile field for",code 
      read_code = 3317

    select_string = 'stash_code:^' + str(read_code) + '$'
    select_field = all_fields.select_field(select_string).copy()

    if code in exclude_codes:
      select_field.setprop("stash_code",str(code))
      select_field.delprop("long_name")
      empty = np.zeros_like(select_field.array)
      select_field.insert_data(cf.Data(empty))
      #is this enough metadata hacking?

    if len(select_field) > 1:
      print "*******"
      print "WARNING: ",select_string," returned multiple fields?"
      print "*******"
    elif len(select_field) == 0:
      print "*******"
      print "WARNING: ",select_string," not found in input file?"
      print "*******"
    else:
      print "extracting ",code
      if code == 8236:
          select_field = select_field/910
      if code == 8576:
          select_field = select_field-273.15

      if select_field.shape[0] < 50:
        if select_field.shape[0] == 9: #elevated tiles not on. Make a one-way coupling field
          print "dummy elevated tiles from old 9 tile config (ice and bare soil) ",select_string
          nelev=10
          dummy=cf.Field()
          for attr in ['long_name','stash_code','um_stash_source','_FillValue']:
            dummy.setprop(attr,select_field.getprop(attr))

          dummy.insert_dim(select_field.dim('dim0'))
          dummy.insert_dim(select_field.dim('dim2'))
          dummy.insert_dim(select_field.dim('dim3'))

          fakeelevs=cf.DimensionCoordinate()
          fakeelevs.setprop('long_name','pseudolevel')
          if code == 8236:
            fakeelevs.insert_data(cf.Data(np.arange(nelev)+926))
            array = select_field.subspace[select_field.coord('long_name:pseudolevel') == 8, ...].array.squeeze()
          else:
            fakeelevs.insert_data(cf.Data(np.arange(nelev)+901))
            array = select_field.subspace[select_field.coord('long_name:pseudolevel') == 9, ...].array.squeeze()
            if code == 3317:
              array[:,:]=1./np.float(nelev)
            if code == 8576:
              array[:,:]=0.

          select_string = 'stash_code:^3317$'
          select_field = all_fields.select_field(select_string)
          icemask = select_field.subspace[select_field.coord('long_name:pseudolevel') == 9, ...].array.squeeze()
          array[np.where(np.abs(icemask - 0.) < 0.1)] = 0.
          fakearray=np.ma.zeros([nelev,np.shape(array)[0],np.shape(array)[1]])
          for n in range(nelev):
            fakearray[n,...]=array

          dummy.insert_dim(fakeelevs)
          dummy.insert_data(cf.Data(fakearray))
          
          select_field = dummy

        else: 
          print "subspace elevated pseudo levels ", code
          if (code == 8236) & (read_code == 8236):
            select_field=select_field.subspace[select_field.coord('long_name:pseudolevel') > 924, ...]
          else:
            select_field=select_field.subspace[select_field.coord('long_name:pseudolevel') > 899, ...]
            select_field=select_field.subspace[select_field.coord('long_name:pseudolevel') < 924, ...]

      extract_list.append(select_field)

  return extract_list

def calculate_area(field_list):

    planet_radius = 6371229.0
    pi            = 3.14159265358979323846

    r_theta = field_list.select_field('stash_code:^33$') + planet_radius
    coastal_fraction = field_list.select_field('stash_code:^3395$')
    tile_fraction = field_list.select_field('stash_code:^3317$')

    latitude_1d = r_theta.coord('latitude').array   * pi / 180.
    longitude_1d = r_theta.coord('longitude').array * pi / 180.
    delta_phi = latitude_1d[1] - latitude_1d[0]
    delta_lambda = longitude_1d[1] - longitude_1d[0]
    latitude_2d = np.zeros([latitude_1d.size,longitude_1d.size])
    for j in range(latitude_1d.size): 
      latitude_2d[j,:] = latitude_1d[j]

    cos_theta_latitude = r_theta.copy()
    cos_theta_latitude.insert_data(cf.Data(np.cos(latitude_2d)))

    #calculate gridbox areas as in um/src/control/coupling/ice_sheet_mass.F90
    cell_area = r_theta * r_theta              \
          * delta_lambda * delta_phi           \
          * cos_theta_latitude                 \
          * coastal_fraction

    tile_area_np = tile_fraction.array
    for k in range(tile_area_np.shape[0]):
      tile_area_np[k,:,:] = tile_area_np[k,:,:] * cell_area.array

    tile_area = tile_fraction.copy()
    tile_area.insert_data(cf.Data(tile_area_np))

    tile_area.setprop('long_name','tile surface area')
    tile_area.id = 'tile_surface_area'
    for attr in ['standard_name','stash_code','um_stash_source']:
      if tile_area.hasprop(attr): tile_area.delprop(attr)

    return tile_area

def change_names(out_list):

  field = out_list.select_field('stash_code:^8236$')
  field.id = 'nonice_snowdepth'
  field.coord('long_name:pseudolevel').setprop('standard_name','tile_id')

  field = out_list.select_field('stash_code:^8576$')
  field.id = 'ice_stemp'
  field.coord('long_name:pseudolevel').setprop('standard_name','tile_id')

  field = out_list.select_field('stash_code:^8578$')
  field.id = 'ice_smb'
  field.coord('long_name:pseudolevel').setprop('standard_name','tile_id')

  field = out_list.select_field('stash_code:^8577$')
  field.id = 'snow_ice_hflux'
  field.coord('long_name:pseudolevel').setprop('standard_name','tile_id')

  field = out_list.select_field('stash_code:^3317$')
  field.id = 'tile_frac'
  field.coord('long_name:pseudolevel').setprop('standard_name','tile_id')

  field = out_list.select_field('long_name:tile surface area')
  field.coord('long_name:pseudolevel').setprop('standard_name','tile_id')

  return out_list

def aggregate_tiles(tile_field,field_list,scaling_method=0):

  print "aggregating a field on elevation tiles:",tile_field.name()

  tile_fraction = field_list.select_field('stash_code:^3317$')

  tile_array = tile_field.array * tile_fraction.array
  aggregate_array = np.sum(tile_array,axis=0)

  for tile in range(np.shape(tile_array)[0]):
      tile_array[tile,:,:] = aggregate_array[:,:]

  tile_field.insert_data(cf.Data(tile_array))

  return tile_field

def fakelapse_tiles(tile_field,field_list):

  nor_sou_boundary = 77.

  lapse_pos_nor = 0.09
  lapse_pos_sou = 0.07
  lapse_neg_nor = 0.56
  lapse_neg_sou = 1.91

  days_in_year = 360.
  secs_in_year = 86400.* days_in_year

  elevations = np.asarray([31.96,301.57,550.11,850.85,1156.89,1453.16,1806.74,2268.25,2753.70,3341.4])

#to do Edwards lapse rate we need: current orog. orog of the levels. Identification of +/- value. Identification of latitude.
  print "faking SMB on elevation tiles using Edwards'14 lapse rates"

  orog_array = field_list.select_field('stash_code:^33$').array.squeeze()
  tile_fraction_array = field_list.select_field('stash_code:^3317$').array.squeeze()
  tile_array_orig = tile_field.array.squeeze()
  latitude_array = tile_field.dim('latitude').array.squeeze()

  tile_array_new = np.zeros_like(tile_array_orig)

  if np.abs(np.sum(tile_array_orig[0,:,:] - tile_array_orig[9,:,:])) > 1e-3:
    print "faking the lapse rates only make sense if we're working with grid box means"
    exit()

  for j in range(np.shape(tile_array_orig)[1]):
    l_pos = False
    if (latitude_array[j] > nor_sou_boundary) | (latitude_array[j] < 0.) : l_nor = True

    for i in range(np.shape(tile_array_orig)[2]):
      #mask on ice boxes only. This field is now *elevation* tile fractions,  may be real or fake
      if np.sum(tile_fraction_array[:,j,i]) > 0.:
        l_pos = False
        if (tile_array_orig[0,j,i] > 0): l_pos = True
  
        lapse = lapse_neg_sou
        if l_nor & l_pos: 
          lapse = lapse_pos_nor
        elif l_nor:
          lapse = lapse_neg_nor
        else:
          lapse = lapse_pos_sou
  
        for tile in range(np.shape(tile_array_orig)[0]):
          dh = elevations[tile] - orog_array[j,i]
          tile_array_new[tile,j,i] = tile_array_orig[tile,j,i] + (lapse * dh)/secs_in_year

  tile_field.insert_data(cf.Data(tile_array_new))

  return tile_field

if __name__ == "__main__":

  #parse command line
  file_in, file_out, gbm_smb, gbm_lapse, exclude_codes, do_snow_shuffling = parse_commandline(sys.argv)

  #extract the fields we want
  all_fields = cf.read(file_in)
  
  process_codes = [33,3395,3317,8578,8236]
  echo_codes = [3317,8576,8577]
  
  process_list = extract_fields(all_fields,process_codes,exclude_codes)
  echo_list = extract_fields(all_fields,echo_codes,exclude_codes)

  #set up the new list of fields to output
  out_list = cf.FieldList()
  
  #create tile areas
  tile_area = calculate_area(process_list)
  out_list.append(tile_area)
  
  ice_smb=process_list.select_field('stash_code:^8578$')

  if  gbm_smb == True:   ice_smb = aggregate_tiles(ice_smb,process_list)

  if  gbm_lapse == True: ice_smb = fakelapse_tiles(ice_smb,process_list)

  # we may already have a file_out with an SMB contribution from last year's nisnow - add in
  if do_snow_shuffling == True:
    if  os.path.isfile(file_out):
      print "SNOW SHUFFLING ON"
      print "ADDING A CONTRIBUTION FROM LAST YEARS's NISNOW CONTRIBUTION TO SMB"
      print "THIS IS GOOD FOR CONSERVATION, BUT EXPERIMENTAL AND HAS LED TO UNPHYSICAL"
      print "SNOW/SMB THINGS IN PAST RUNS"

      f=cf.read(file_out)

      try: 
        nisnow_smb=f.select_field('long_name:smb_from_nisnow')
        print "existing smb", np.sum(ice_smb.array)
        print "addition to add", np.sum(nisnow_smb.array)
        ice_smb.data=ice_smb.data+nisnow_smb.data
        ##add the nisnow back to the saved file, for reference
        ##getting some odd errors to do with the time axis if I add this back in....
        #out_list.append(nisnow_smb)
        ##write to its own file for reference
        nisnow_file=file_out.rsplit("_",1)[0]+"-nisnow_"+file_out.rsplit("_",1)[1]
        cf.write(nisnow_smb,nisnow_file, fmt = "NETCDF4")

        print "new smb", np.sum(ice_smb.array)

      except:
        print "something went wrong looking for the nisnow field in the file - it's probably not there"
        print "skipping adding this contribution in after all"

  else:

    print "SNOW SHUFFLING OFF"
    print "NOT ADDING PREVIOUS YEAR's NISNOW POTENTIAL CONTRIBUTION TO SMB"

  out_list.append(ice_smb)

  nisnow=process_list.select_field('stash_code:^8236$')

  if  gbm_smb == True: nisnow = aggregate_tiles(nisnow,process_list)

  if do_snow_shuffling == False:
    print "SNOW SHUFFLING OFF"
    print "ZEROING NISNOW PASSED INTO BISICLES"
    nisnow.data=nisnow.data-nisnow.data

  out_list.append(nisnow)

  #add the rest of the fields to the output
  for field in echo_list: out_list.append(field)
  
  #change netCDF_out names to human-readable, common format with JULES
  out_list = change_names(out_list)

  #write new file, overwriting the nisnow stuff if necessary
  cf.write(out_list,file_out,fmt = "NETCDF4", overwrite=True)
