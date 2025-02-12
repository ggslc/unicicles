#BIKEtoNEMO
#
import numpy as np
import sys
from amrfile import io as amrio
from ferretlook import *
hdfplotfile="plot.ant-meltx-tsai.2lev.4sub.000000.2d.hdf5"
regrid_ncfile="regridded.nc"
bike_ncgridfile="cf_gridfile_BISICLES_lev0.nc"
nemo_ncgridfile="cf_gridfile_eORCA025_v2.3.nc"
bathy_threshold=-5000
isf_threshold=-9999
level=0 # level of grid refinement
order=0 # interpolation order, 0 for piecewise constant, 1 for linear

amrID = amrio.load(hdfplotfile)
lo,hi = amrio.queryDomainCorners(amrID, level)

variable="thickness"

_,_,isf = amrio.readBox2D(amrID, level, lo, hi, variable, order)
import cf

g=cf.read(nemo_ncgridfile)
f=cf.read(bike_ncgridfile)
f.insert_data(cf.Data(isf,units='m'), axes=('X', 'Y'))
f_regrid=f.regrids(g,dst_cyclic=True,method='bilinear')
remap_isf_np=np.rollaxis(f_regrid.array,1)
ferretlook(remap_isf_np)

#NEMOtoBIKE
#
import os
import numpy as np
import sys
import cf
from ferretlook import *
nctgridfile="eORCA025_bathymetry_isf_v2.3.nc"
regrid_hdf5file="regridded.hdf5"
regrid_ncfile  =os.path.splitext(os.path.basename(regrid_hdf5file))[0]+".nc"
bike_ncgridfile="cf_gridfile_BISICLES_lev0.nc"
nemo_ncgridfile="cf_gridfile_eORCA025_v2.3.nc"

water_unit_factor=60*60*24*360/1e3
heat_unit_factor=60*60*24*360
h=cf.read(nctgridfile)
variable="ncvar%isf_draft"
melt_water=h.select(variable).array.squeeze()
melt_water=np.ma.masked_equal(melt_water,0)
g=cf.read(bike_ncgridfile)
f=cf.read(nemo_ncgridfile)
f.insert_data(cf.Data(melt_water,units='kg/m2/s'), axes=('Y', 'X'))
w_regrid=f.regrids(g,src_cyclic=True,method='bilinear',i=True)
melt=np.swapaxes(w_regrid.array,0,1)
ferretlook(melt)
