import numpy as np
import mule
import mule_rss

def remask_var(array,missing_flag,newlsm):

  from scipy import ndimage

  #http://stackoverflow.com/questions/3662361/fill-in-missing-values-with-nearest-neighbour-in-python-numpy-masked-arrays
  array_masked = np.ma.masked_equal(array,missing_flag)

  sea_pts=np.where(newlsm < -0.5)
  #print len(sea_pts[0])," sea points in new mask"
  ind = ndimage.distance_transform_edt(array_masked.mask, \
                                       return_distances=False, \
                                       return_indices=True)

  a=array[tuple(ind)]

  if len(sea_pts[0]) > 0:
    a[sea_pts]=missing_flag

  return a


def new_lsm(umf,lsm_file):

  # should use the local module with codes in, keep
  # everything in one place
  lsm_stash_bin=[30]
  lsm_stash_frac=[505,3395]

  fland = lsm_file.variables['lsmask'][:].squeeze()
  fland = np.ma.filled(fland,fill_value=0.)
  minval = .01
  new_landindex = np.where(fland > minval)

  index_LSM = mule_rss.findindex(umf,lsm_stash_bin[0])[0]
  LSM_old   = mule_rss.get_data2d(umf,index_LSM)
  LSM_new   = LSM_old.copy()

  index_fLSM = mule_rss.findindex(umf,lsm_stash_frac[0])[0]
  fLSM_old   = mule_rss.get_data2d(umf,index_fLSM)
  fLSM_new   = fLSM_old.copy()

  LSM_new[:,:]=0.
  LSM_new[new_landindex] = 1

  fLSM_new[:,:]=0.
  fLSM_new[new_landindex] = fland[new_landindex]

  #THIS ISN'T ENOUGH! WE NEED TO PROCESS ALL LAND FIELDS INTO GLOBAL ONES (USING THE OLD MASK) 
  #AND BACK AGAIN (WITH THE NEW MASK) OR THEY'LL GET MANGLED 
  #umf = mule_rss.put_data2d(umf,index_LSM,LSM_new)
  #umf = mule_rss.put_data2d(umf,index_fLSM,fLSM_new)

  #initilise handy counters
  last_stash = -99
  count = 0

  for index in range(len(umf.fields)):
  #for index in range(umf.fixed_length_header.raw[153]):
  
    #check stash for this field in the loop
    this_stash   = umf.fields[index].lbuser4
    missing_flag = umf.fields[index].bmdi

    #handy to know which level we're on
    if this_stash == last_stash: count = count + 1
    if this_stash > last_stash: count = 1

    last_stash = this_stash

    #read field from dump
    field = mule_rss.get_data2d(umf,index)

    if this_stash in lsm_stash_bin: 
      print "putting binary mask direct into ",this_stash
      field = LSM_new
    #my UM4.5 did specific things to sea-ice too so there were valid SSTs/albedo
    #where we made new ocean. Essential?
    #open sea sst
    #if this_stash == 507: field = ?

    missing = np.where(field == missing_flag)
    valid   = np.where(field != missing_flag)

    #is it safe to assume all fields that have missing and are the same
    #size as the mask are land and can be masked? Incl diagnostics?
    if (len(missing[0]) > 0) & (len(field) == len(LSM_new)) > 0:

      #print this_stash, "starts with missing, avg", len(missing[0]), np.mean(field[valid])
      if this_stash in lsm_stash_frac: 

        field = fLSM_new
        print "putting fractional mask direct into ",this_stash

      else:

        print "crude near-neighbor fill for stash ",this_stash, count, field.shape
        field = remask_var(field, missing_flag, LSM_new)

        #missing = np.where(field == missing_flag)
        #valid   = np.where(field != missing_flag)
        #print this_stash, "ends with missing, avg", len(missing[0]), np.mean(field[valid])

    umf = mule_rss.put_data2d(umf,index,field)

  #need to check how landmask is treated later on - is it got from here, or implied from the 
  #fractional areas etc from previous tile fileds in UM or ice?

  return umf

  
