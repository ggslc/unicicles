import numpy as np
import mule
import mule_rss

import sys

def parse_commandline(argv):
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument("--region",    help="input, which ice sheet")
  parser.add_argument("--flip_array",     help="input, round fractional mask to binary",action="store_true")
  parser.add_argument("--binary",          help="input, flip array N-S (not meta)",action="store_true")
  parser.add_argument("--input_global",    help="input, name of global file")
  parser.add_argument("--output_regional", help="output, name of regional file")

  args = parser.parse_args()

  for arg in vars(args):
    value = getattr(args, arg)
    if value != None: 
      ##ooh sketchy. You can make a variable named like the input string
      globals()[arg] = value
    else:
      ##between this and the fact that argparse will get upset with any arguments
      ##not listed above this assignment isn't massively unsafe(?!)
      print("give a value for commandline arg:",arg)
      parser.print_help(sys.stderr)
      exit()

  return input_global, output_regional, region, flip_array, binary

infile, outfile, region, flip, binary = parse_commandline(sys.argv)

fieldg = mule.load_umfile(infile)
fieldg.level_dependent_constants=mule.ff.FF_LevelDependentConstants(np.zeros([fieldg.integer_constants.num_levels,4]))

#report what's there
nx_orig=fieldg.integer_constants.raw[6]
ny_orig=fieldg.integer_constants.raw[7]

if region == "AIS":
#Cutout Ant from N96e
  nx=192
  ny=25
  x_offset=0
  y_offset=0
elif region == "GrIS":
#Cutout GrIS from N96e (-1 at each edge to keep CAP domain from failing/wrapping?)
  nx=54
  ny=25
  x_offset=137
  y_offset=118
else: 
  print("Don't know parameters for region",region," exiting")
  exit(1)

#Cutout GrIS+ from GLOBE30
#nx=13200
#ny=5400
#x_offset=30000
#y_offset=0
#Cutout Ant from GLOBE30
#nx=43200
#ny=4200
#x_offset=0
#y_offset=17400

dx=fieldg.fields[0].raw[62]
dy=fieldg.fields[0].raw[60]

print("X",nx,dx,fieldg.real_constants.raw[4] + (x_offset*dx), fieldg.real_constants.raw[4] + (x_offset*dx) + nx*dx)
print("Y",ny,dy,fieldg.real_constants.raw[3] + (y_offset*dy), fieldg.real_constants.raw[3] + (y_offset*dy) + ny*dy)


if (nx < fieldg.integer_constants.raw[6]) or (ny < fieldg.integer_constants.raw[7]):
#change from global to non-global latlon grid without wrap
  print("setting to non-global grid type")
  fieldg.fixed_length_header.raw[4]=3

# switch to newer version of coding the UM version
fieldg.fixed_length_header.raw[12]=901
#1 should always be MDI
fieldg.fixed_length_header.raw[1]=-32768

fieldg.integer_constants.raw[6]=nx
fieldg.integer_constants.raw[7]=ny

fieldg.real_constants.raw[3]=fieldg.real_constants.raw[3] + (y_offset*dy)
fieldg.real_constants.raw[4]=fieldg.real_constants.raw[4] + (x_offset*dx)


for i in range(len(fieldg.fields)):
  fieldg.fields[i].raw[17]=3
  fieldg.fields[i].raw[18]=ny
  fieldg.fields[i].raw[19]=nx
  fieldg.fields[i].raw[59]=fieldg.fields[i].raw[59]  + (y_offset*dy)
  fieldg.fields[i].raw[61]=fieldg.fields[i].raw[61]  + (x_offset*dx)

  datag=fieldg.fields[i].get_data()
  if flip:   
    print("flipping input array NS, apparently it was generated wrong")
    datag = datag[-1:0:-1,:]
  if binary: 
    print("I think I've been given a fraction landmask, make it binary")
    #by eyeballing what I currently use the binary mask is around this cf fractional!?
    fieldg.fields[0].raw[23]=38
    threshold=0.5
    datag[np.where(datag >= threshold)] = 1.
    datag[np.where(datag < threshold)] = 0.
    datag = datag.astype(int)
    fieldg.fields[i].raw[42] = 30

#write as logical data
  fieldg.fields[i].raw[38] = 9011111
  fieldg.fields[0].raw[39] = 3

  fieldg.fields[i]=mule_rss.overwrite_data(fieldg.fields[i],datag[y_offset:y_offset+ny,x_offset:x_offset+nx])

fieldg.to_file(outfile)
