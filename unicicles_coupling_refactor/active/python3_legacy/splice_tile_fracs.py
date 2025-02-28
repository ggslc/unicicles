import numpy as np
from params_and_constants import *
from um_to_np_utility import ice_to_np,               \
                             get_ice_um_grid_overlap, \
                             shift_tile_dim,          \
                             um_to_np_3d,             \
                             merge_ice_um_arrays,     \
                             np_to_um_3d

def fix_polarrow_tilefrac(ice_file):

    #tile_land_frac from Glint is the amount of land it knows about in an elevation class
    #depending on Glint mask settings it may discriminate between land and sea, or simply
    #whether a point is in the regional domain Glint knows about from the GCM
    #Either way, isolated gridpoints on the South polar rows of UKESM with a tile-sum less
    #than 1 are incorrect, some artefact of rebinning Glints polar stereographic to the the UM
    #lat-lon. Results in small artefacts near the poles when we merge Glint info back into UM

    field=ice_file.variables['tile_land_fraction'][:]

    #only two polar rows affected
    for index in range(2):
      #Glint files are N first, need last two for S 
      row = field.shape[1]-1-index
      for col in range(field.shape[2]):
        tilesum = np.sum(field[:,row,col],axis=0)
        #make sure we end with a tilefrac sum of 1
        if ( tilesum < 1e-6 ):
          field[-1,row,col] = 1.
        elif ( tilesum < 1 ):
          field[:,row,col] = field[:,row,col] / tilesum

    ice_file.variables['tile_land_fraction'][:] = field[:]

   #Do the same with ice fraction?

    field=ice_file.variables['tile_ice_fraction'][:]

    #only two polar rows affected
    for index in range(2):
      #Glint files are N first, need last two for S 
      row = field.shape[1]-1-index
      for col in range(field.shape[2]):
        tilesum = np.sum(field[:,row,col],axis=0)
        #make sure we end with a tilefrac sum of 1
        if ( tilesum < 1e-6 ):
          field[-1,row,col] = 1.
        elif ( tilesum < 1 ):
          field[:,row,col] = field[:,row,col] / tilesum

    ice_file.variables['tile_ice_fraction'][:] = field[:]

    return ice_file

def splice_tile_fracs(um_dump,um_ref_dump,ice_file):

  ##for BIKE, land fraction is frac of BIKE cells for that GCM gridbox 
  ##that is in a given GCM elevation
  ##ice fraction is the proportion of the land fraction that is ice
  land_frac_ice=ice_to_np(ice_file, 'tile_land_fraction')
  ice_frac_ice =ice_to_np(ice_file, 'tile_ice_fraction')
  #get BISICLES orog on the UM grid - this is the best thing to use as a mask
  #for where the ISM domain has real data worth adding
  orog_ism=ice_to_np(ice_file, 'surface_elevation')

  nelev_ice=np.shape(land_frac_ice)[2]

  #get overlap between UM and ice grids
  fmask_ice2d=get_ice_um_grid_overlap(ice_file)

  #to get to the UM fraction convention, we need to do icefrac*landfrac
  #The fmask_ice normalises out the partial coverage weighting from the edge of
  #the ice grid in these quantities - makes things clearer later
  fmask_ice3d   =shift_tile_dim(np.tile(fmask_ice2d,[nelev_ice,1,1]))
  fmask_good=np.where(fmask_ice3d > 0)

  ice_frac_gcm  =np.zeros_like(ice_frac_ice)
  land_frac_gcm =np.zeros_like(ice_frac_ice)

  ice_frac_gcm[fmask_good]=( ice_frac_ice[fmask_good] *  \
                             land_frac_ice[fmask_good] ) \
                           / fmask_ice3d[fmask_good]


  land_frac_gcm[fmask_good]=( (1-ice_frac_ice[fmask_good]) * \
                              land_frac_ice[fmask_good] )    \
                            / fmask_ice3d[fmask_good]


  frac_um, pslev =um_to_np_3d(um_ref_dump, stashcode_frac_surf_type, elev=True)

  ##check how many elevation classes we have, and what the max is
  max_elev_id=np.max(pslev)
  nelev_um=np.mod(max_elev_id, max_number_elevs)

  if nelev_um != nelev_ice:
     print("ERROR in splice_tile_fracs: nelev_um != nelev_ice - OOPS")
     sys.exit(2)

  #update where there is some ice domain coverage and the UM allows elevated tiles 
  fmask_um2d=np.sum(frac_um,axis=2)
  frac_updates=np.where((orog_ism > 0.) & (fmask_um2d > 0.5))

  #index elev is the ice-type for an elevation band
  #index elev+nelev_um is the rock-type for that elevation band
  #as for orog, mix in the new ice-derived fields with the old UM ones weighted by the
  #partial grid coverage of the ice grid at the edges (fmask_ice)
  frac_new=np.copy(frac_um)
  for elev in range(nelev_um):
    ice_frac_gcm_1lev=ice_frac_gcm[:,:,elev]
    land_frac_gcm_1lev=land_frac_gcm[:,:,elev]

    #we have rock and ice types in the UM
    if max_elev_id >= (start_tile_id_elev_ice+max_number_elevs):
      ice_frac_um_1lev=frac_um[:,:,elev]
      rock_frac_um_1lev=frac_um[:,:,elev+nelev_um]

      ice_frac_um_1lev=merge_ice_um_arrays(ice_frac_um_1lev,    \
                                           ice_frac_gcm_1lev,   \
                                           fmask_ice2d,         \
                                           points=frac_updates, \
                                           noFrac=True)

      rock_frac_um_1lev=merge_ice_um_arrays(rock_frac_um_1lev,  \
                                            land_frac_gcm_1lev, \
                                            fmask_ice2d,        \
                                            points=frac_updates, \
                                            noFrac=True)

      frac_new[:,:,elev]          = ice_frac_um_1lev
      frac_new[:,:,elev+nelev_um] = rock_frac_um_1lev
    #we only have ice types in the UM
    else:
      ice_frac_um_1lev=frac_um[:,:,elev]
      ice_frac_um_1lev=merge_ice_um_arrays(ice_frac_um_1lev,                     \
                                           land_frac_gcm_1lev+ice_frac_gcm_1lev, \
                                           fmask_ice2d,                          \
                                           points=frac_updates, \
                                           noFrac=True)

      frac_new[:,:,elev]=ice_frac_um_1lev

  um_dump=np_to_um_3d(um_dump, stashcode_frac_surf_type, frac_new, elev=True)
 
  return um_dump
