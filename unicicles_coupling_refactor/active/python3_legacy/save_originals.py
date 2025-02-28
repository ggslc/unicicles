from params_and_constants import *
from um_to_np_utility import um_to_np_3d

def save_originals(um_dump):

  frac_elev0,_      =um_to_np_3d(um_dump, stashcode_frac_surf_type,   elev=True)
  snicemass_elev0,_ =um_to_np_3d(um_dump, stashcode_snow_tile,        elev=True)

  return frac_elev0,snicemass_elev0
