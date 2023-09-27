import numpy as np
import cf
import os.path

def cmdline_args():
  import sys

  if len(sys.argv) != 2:
    print 'Usage: extract_fields.py <inputfile>'
    sys.exit("wrong command line options")

  fname=sys.argv[1]

  if not os.path.isfile(fname):
    sys.exit("extract_fields.py: input file does not exist")

  return fname

def extract_fields(all_fields,extract_codes):

  ExtractList=cf.FieldList()

  for code in extract_codes: 
    select_string='stash_code:^'+str(code)+'$'
    select_field=all_fields.select(select_string)
    #select_field.subspace[select_field.coord('long_name:pseudolevel')>900]
    if len(select_field) > 1:
      print "*******"
      print "WARNING: ",select_string," returned multiple fields?"
      print "*******"
    elif len(select_field) == 0:
      print "*******"
      print "WARNING: ",select_string," not found in input file?"
      print "*******"
    else:
      print "extracting ",select_string
      if code == 8236:
          select_field=select_field/910
      if code == 8576:
          select_field=select_field-273.15
      ExtractList.append(select_field)

  return ExtractList

def calculate_area(fieldlist):

    import numpy as np

    planet_radius = 6371229.0
    pi            = 3.14159265358979323846

    r_theta=fieldlist.select('stash_code:^33$')+planet_radius
    coastal_fraction=fieldlist.select('stash_code:^3395$')
    tile_fraction=fieldlist.select('stash_code:^3317$')

    latitude_1d=r_theta.coord('latitude').array   * pi/180.
    longitude_1d=r_theta.coord('longitude').array * pi/180.
    delta_phi=latitude_1d[1]-latitude_1d[0]
    delta_lambda=longitude_1d[1]-longitude_1d[0]
    latitude_2d=np.zeros([latitude_1d.size,longitude_1d.size])
    for j in range(latitude_1d.size): 
      latitude_2d[j,:]=latitude_1d[j]

    cos_theta_latitude=r_theta.copy()
    cos_theta_latitude.insert_data(cf.Data(np.cos(latitude_2d)))

    #calculate gridbox areas as in um/src/control/coupling/ice_sheet_mass.F90
    cell_area = r_theta * r_theta                  \
          * delta_lambda * delta_phi           \
          * cos_theta_latitude                 \
          * coastal_fraction

    tile_area_np=tile_fraction.array
    for k in range(tile_area_np.shape[0]):
      tile_area_np[k,:,:]=tile_area_np[k,:,:]*cell_area.array

    tile_area=tile_fraction.copy()
    tile_area.insert_data(cf.Data(tile_area_np))

    tile_area.setprop('long_name','tile surface area')
    tile_area.id='tile_surface_area'
    for attr in ['standard_name','stash_code','um_stash_source']:
      if tile_area.hasprop(attr): tile_area.delprop(attr)

    return tile_area

def fixfield(good,bad):
  fixed=good.copy()
  fixed.insert_data(cf.Data(bad.array))
  fixed.id=bad.id
  for attr in ['add_offset','cell_methods','comment','Conventions','_FillValue', \
               'flag_masks','flag_meanings','flag_values','history','long_name', \
               'missing_value','references','scale_factor','stash_code',         \
               'standard_error_multiplier','standard_name','title','units',      \
               'valid_max','valid_min','valid_range','um_stash_source']:
    if fixed.hasprop(attr):
      if bad.hasprop(attr):
        fixed.setprop(attr,bad.getprop(attr))
      else:
        fixed.setprop(attr,'none')

  return fixed

def change_names(OutList):

  field=OutList.select('stash_code:^8236$')
  field.id='nonice_snowdepth'
  field.coord('pseudolevel').setprop('standard_name','tile_id')

  field=OutList.select('stash_code:^8576$')
  field.id='ice_stemp'
  field.coord('pseudolevel').setprop('standard_name','tile_id')

  field=OutList.select('stash_code:^8578$')
  field.id='ice_smb'
  field.coord('pseudolevel').setprop('standard_name','tile_id')

  field=OutList.select('stash_code:^8577$')
  field.id='snow_ice_hflux'
  field.coord('pseudolevel').setprop('standard_name','tile_id')

  field=OutList.select('stash_code:^3317$')
  field.id='tile_frac'
  field.coord('pseudolevel').setprop('standard_name','tile_id')

  field=OutList.select('long_name:tile surface area')
  field.coord('pseudolevel').setprop('standard_name','tile_id')

  return OutList

good_codes=[33,3395,3317,8236]
bad_codes=[8576,8577,8578]
template_code=3317

#parse command line
file_in=cmdline_args()

file_out=os.path.splitext(file_in)[0]+'.IceCoupling.nc'

#extract the fields we want
all_fields=cf.read(file_in)

GoodList=extract_fields(all_fields,good_codes)
BadList=extract_fields(all_fields,bad_codes)

#create tile areas
tile_area=calculate_area(GoodList)
tile_fraction=GoodList.select('stash_code:^3317$')

OutList=cf.FieldList()
for field in [tile_area,tile_fraction \
              ,GoodList.select('stash_code:^8236$') \
              ]:
  OutList.append(field)

#work around the weird dimensioning cf-python has made
#from the incorrect lbvc pp_header on some fields
template_field=GoodList.select('stash_code:'+str(template_code))
for field in BadList:
  OutList.append(fixfield(template_field,field))

#change netCDF_out names to human-readable, common format with JULES
OutList=change_names(OutList)

#write new file
cf.write(OutList,file_out)
