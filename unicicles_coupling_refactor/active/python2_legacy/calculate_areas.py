import mule_rss
import numpy as np
from params_and_constants import *
from um_to_np_utility import um_to_np_2d, \
                             um_to_np_3d

def calculate_areas(um_dump, frac_um=None):

    """
    calculate actual gridbox and tile-fraction areas. 
    cell area calc as in /src/control/coupling/ice_sheet_mass.F90
    """

    r_theta=um_to_np_2d(um_dump,stashcode_orog)+planet_radius
    coastal_fraction=um_to_np_2d(um_dump,stashcode_land_frac)
    if frac_um == None:
      frac_um,_ =um_to_np_3d(um_dump, stashcode_frac_surf_type, elev=True)

    nx=np.shape(coastal_fraction)[0]
    ny=np.shape(coastal_fraction)[1]
    ntiles=np.shape(frac_um)[2]

    index=mule_rss.findindex(um_dump,stashcode_orog)
    start_phi=um_dump.fields[index].bzy    * pi/180.
    delta_phi=um_dump.fields[index].bdy    * pi/180.
    delta_lambda=um_dump.fields[index].bdx * pi/180.

    latitude_1d=start_phi+np.arange(ny,dtype=float)*delta_phi
    latitude_2d=np.zeros_like(coastal_fraction)
    for j in range(ny): latitude_2d[:,j]=latitude_1d[j]

    cos_theta_latitude=np.cos(latitude_2d)

    cell_area = r_theta * r_theta                  \
          * delta_lambda * delta_phi           \
          * cos_theta_latitude                 \
          * coastal_fraction

    tile_area=np.zeros_like(frac_um)
    for k in range(ntiles):
      tile_area[:,:,k]=frac_um[:,:,k]*cell_area

    return tile_area, cell_area
