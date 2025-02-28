import numpy as np

def assig_4p4(orog_sd, assig_coef=5E-5, h_coef=0.15, las_linear=True):
  """

* Subroutine to calculate the A/S parameter from sigma-h

# Based on the pre-GLOBE CAP routine (here vn4.4)
# We also return the half-peak-to_trough field - the CAP just set 
# this to the standard deviation, optionally scaled by a namelist
# parameter.
# Eyeballing G'land and Ant cf N96 orig ancil and compromising,
# simple linear scalings for the earliest UKESM1-IS were as 
# the defaults above. Doesn't get the highs, overestimates the lows
# In general, factors are too small for Gland and too high for Ant
### N216 eyeballed: assig_coef=1e-4 h_coef=0.3###


*from UM CAP@vn4.4*
*
*
* There are two possible methods
*
* 1) A simple linear relationship (LAS_LINEAR)
*
*  A/S =a * sigma-h
*
* where a is a coefficent (argument ASSIG_COEF)
*
* a is model grid dependent and can be calculated externally
* and fed to the program using the namelist or calculated
* internally
#
# we can't do this internal calculation - relies on there being
# pre-existing regions where we know both orog_sd and orog_as
# and we don't have those here 
*
*  Standard values are  (10' data)          (5'data)
*  oper glob            a=0.0001679         a=0.000156687
*  oper LAM             a=0.0002986         a=0.00027718
*  climate              a=0.000084454       a=0.0000743820
# in vn4.5 UM-Glimmer coupling, we used
#  HadCM3               a=8.4e-5
#  FAMOUS               a=5.9e-5
*
* 2) A more complicated method
*
*   A/S=  2k**2/Cd/(ln(zc/zo))**2
*
*  k=von Karman's constant (0.4)
*  Cd= drag coeffiecent =0.3
*  zc=SQRT(2)SIGMA
*  zo=aSIGMA**2   for sigma <= BIGSIG
*    =aBIGSIG(2SIGMA-BIGSIG) otherwise
* where
*  BIGSIG=100    a=5E-4

  """

  valid = np.where(orog_sd >0.)

  orog_h = np.zeros_like(orog_sd)
  orog_h[valid]  = orog_sd[valid] * h_coef

  orog_as = np.zeros_like(orog_sd)

  if las_linear==True:

    orog_as[valid] = orog_sd[valid] * assig_coef

    orog_as = np.maximum(orog_as,0.)
    orog_as = np.minimum(orog_as,0.2)

  else:

    k=0.4
    cd=0.3
    a=5e-4
    bigsig=100.

    big   = np.where(orog_sd > bigsig)

    zc = np.sqrt(2.0)*orog_sd
    z0 = a * (orog_sd)**2.0
    z0[big] = a * bigsig * ( 2.0 * orog_sd[big] - bigsig )

    as1=2.0*k*k/cd
    as2=(np.log( zc[valid]/z0[valid] ))**2.0

    orog_as[valid]=as1/as2

  return orog_as, orog_h
