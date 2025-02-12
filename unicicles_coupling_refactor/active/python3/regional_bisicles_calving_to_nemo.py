#!/usr/bin/env python

"""
regional_bisicles_calving_to_nemo

The plot.*hdf5 output directly from BISICLES is read and the field of
accumulated calved ice is processed.  The field is nearest-neighbour
regridded (it is not an area-average) onto the NEMO grid and routed to the
nearest (treating the NEMO array as uniform) open-ocean gridbox.
The same BISICLES field is used to produce a file containing a sink term
that will remove the calved ice from BISICLES over the next cycle.

Called by suite as regional_ice_sheet_calving_to_ocean
"""

import sys
import os
import numpy as np
from common_arg_to_file_exist import arg_to_file_exist
from common_params_and_constants import days_in_year, \
    epsg_AIS, epsg_GrIS, epsg_global, secs_in_day, \
    x0_AIS, x0_GrIS, y0_AIS, y0_GrIS


def parse_commandline():
    """ Read the line command line arguments """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--hdf5_input", help="input, name of unicicles hdf5")
    parser.add_argument("--hdf5_output",
                        help="output, name of unicicles melange sink hdf5")
    parser.add_argument("--nc_output",
                        help="output, name of nemo iceberg seed nc")
    parser.add_argument("--coupling_period",
                        help="climate-ice coupling period in seconds")
    parser.add_argument("--ice_tstep_multiply", help="is ice accelerated?")
    parser.add_argument("--region", help="GrIS or AIS?")
    args = parser.parse_args()

    err = 0
    hdf5_input, err = arg_to_file_exist(args.hdf5_input, err=err)
    if err > 0:
        parser.print_help()
        sys.exit(err)
    nc_output, err = arg_to_file_exist(args.nc_output, io="out", err=err)
    if err > 0:
        parser.print_help()
        sys.exit(err)
    hdf5_output, err = arg_to_file_exist(args.hdf5_output, io="out", err=err)
    if err > 0:
        parser.print_help()
        sys.exit(err)

    if args.coupling_period is not None:
        coupling_period = float(args.coupling_period)
    else:
        err = err + 1
        print("ERROR: specifiy a climate-ice coupling period")

    if args.ice_tstep_multiply is not None:
        ice_tstep_multiply = float(args.ice_tstep_multiply)
    else:
        ice_tstep_multiply = 1.

    region = "GrIS"
    if args.region is not None:
        region = args.region

    if err > 0:
        parser.print_help()
        sys.exit(err)

    return hdf5_input, hdf5_output, nc_output, coupling_period, \
        ice_tstep_multiply, region


def extract_bisicles_calving(hdf5_input, region="GrIS"):
    from amrfile import io as amrio
    from pyproj import Transformer

    """
    Pull the amount of calved ice out of the BISICLES plot files and return it
    along with the longitude and latitudes of the grid centres
    """

    level = 0  # level of grid refinement, 0 = lowest
    order = 0  # interpolation order, 0 for piecewise constant, 1 for linear

    # Load BISICLES hdf5
    amrID = amrio.load(hdf5_input)
    lo, hi = amrio.queryDomainCorners(amrID, level)

    # Extract cartesian grid and variables. x and y are the box centres
    # distance from origin. 
    #
    # We actually use the inventory of calved ice from the melange model
    # in BISICLES rather than the diagnostic of calving rate itself
    x_1d, y_1d, calving = amrio.readBox2D(amrID, level, lo, hi,
                                          "melangeThickness", order)

    # Make 2D meshes of x,y as coordinate arrays to match the data
    dx = x_1d[1] - x_1d[0]
    dy = y_1d[1] - y_1d[0]
    y_2d, x_2d = np.mgrid[y_1d[0]:y_1d[-1] + dy:dy, x_1d[0]:x_1d[-1] + dx:dx]

    # select projection defintions for pyproj to use
    outProj = epsg_global

    # If you set up a new BISICLES domain, check that coastlines etc look 
    # right either side of the transformation - these offsets have changed 
    # in the past, odd reflections have been applied etc.
    if region == "GrIS":
        inProj = epsg_GrIS
        x0 = x0_GrIS
        y0 = y0_GrIS
    elif region == "AIS":
        inProj = epsg_AIS
        x0 = x0_AIS
        y0 = y0_AIS
    else:
        print("extract_bisicles_calving: region unknown ", region)

    # Transform from the BISICLES grid info to lons and lats.

    # Original AIS definition needed longitude reflection
    # if region == "AIS": mfact = -1
    mfact = 1.

    t = Transformer.from_crs(inProj, outProj, always_xy=True)
    xl_2d, yl_2d = t.transform((x_2d + x0) * mfact, (y_2d + y0))

    # Original AIS definition needed *additional* longitude offset
    # if region == "AIS": steph_offset=90.
    steph_offset = 0.
    xl_2d = xl_2d + steph_offset

    return calving, xl_2d, yl_2d, x_1d, y_1d


def construct_CF(field, lon, lat, x_1d, y_1d):

    import cf

    """
    Construct a minimal CF Field from a numpy array and some lon-lat data
    so that cf-python can then use it in regridding routines
    """
    field_CF = cf.Field()

    dimx = cf.DimensionCoordinate(data=cf.Data(x_1d, 'm'))
    field_CF.set_construct(cf.DomainAxis(size=len(x_1d)), key="X")
    field_CF.set_construct(dimx, axes="X")

    dimy = cf.DimensionCoordinate(data=cf.Data(y_1d, 'm'))
    field_CF.set_construct(cf.DomainAxis(size=len(y_1d)), key="Y")
    field_CF.set_construct(dimy, axes="Y")
    lats = cf.AuxiliaryCoordinate(data=cf.Data(lat, 'degrees_north'),
                                  properties={'standard_name': 'latitude'})
    field_CF.set_construct(lats, axes=('Y', 'X'))
    lons = cf.AuxiliaryCoordinate(data=cf.Data(lon, 'degrees_east'),
                                  properties={'standard_name': 'longitude'})
    field_CF.set_construct(lons, axes=('Y', 'X'))

    field_ma = np.ma.masked_equal(field, 0.)
    field_CF.set_data(cf.Data(field_ma, units='m'), axes=('Y', 'X'))
    field_CF.set_property('_FillValue', -99.)

    return field_CF


def make_BISICLES_melange_sink(bike, x, y, coupling_interval_secs, calv_hdf):

    from netCDF4 import Dataset

    """
    Make a BISICLES hdf5 file to use as a boundary condition to extract
    the amount of calved ice we've just given to NEMO from the BISICLES
    inventory over the next coupling period
    """

    # Turn the melange mass (metres of ice) we've just taken into a
    # flux (m/s) to EXTRACT from BISICLES
    melange_sink = np.copy(bike) * -1 / coupling_interval_secs

    # Turn the data into a BISICLES hdf5 via netcdf and the filesystem utility
    ncfile_out = Dataset("calv_tmp.nc", 'w', format='NETCDF3_CLASSIC')
    ncfile_out.createDimension('x', x.shape[0])
    ncfile_out.createDimension('y', y.shape[0])

    x_nc = ncfile_out.createVariable('x', np.dtype('float64').char, ('x'))
    y_nc = ncfile_out.createVariable('y', np.dtype('float64').char, ('y'))
    data_nc = ncfile_out.createVariable('calved_melange',
                                        np.dtype('float64').char, ('y', 'x'))

    x_nc[:] = x
    y_nc[:] = y
    data_nc[:] = melange_sink
    ncfile_out.close()

    # call filesystem utility
    cmd = './nctoamr2d.ex calv_tmp.nc ' + calv_hdf + ' calved_melange'
    status = os.system(cmd)
    if status != 0:
        print("ERROR: failed to run nctoamr2d.ex")
        sys.exit(1)
    cmd = 'rm calv_tmp.nc'
    status = os.system(cmd)
    if status != 0:
        print("ERROR: failed to remove calv_tmp.nc")
        sys.exit(1)

    return


def regrid_to_NEMO(field_CF, nemo_ncgridfile):

    import cf

    """
    regrid the BISICLES cartesian iceberg field onto the NEMO orca grid
    """

    # NEMO actually wants a volume of ice at a point - NOT a cell average
    # THIS MAY CHANGE IN FUTURE RELEASES?
    # Multiply field by the cell areas here, on the nice regular grid, and do a
    # simple nearest-neighbour point map to the orca grid
    #python3 cf returns lon/lat with the old syntax of asking for X and Y, we want to work
    #out cell lengths here
    key_dy=list(field_CF.dimension_coordinates().keys())[0]
    key_dx=list(field_CF.dimension_coordinates().keys())[1]

    dx = field_CF.dimension_coordinate(key_dx).array[1] - field_CF.dimension_coordinate(key_dx).array[0]
    dy = field_CF.dimension_coordinate(key_dy).array[1] - field_CF.dimension_coordinate(key_dy).array[0]
    field_CF = field_CF * dx * dy

    #masking in python3 cf is screwing up the regrid
    f1=np.asarray(field_CF.array)
    field_CF.set_data(cf.Data(f1))

    g = cf.read(nemo_ncgridfile)[0]
    # Different versions of cf-python seem to label (and identify) our too-simple
    # constructed field axes in different ways. Adding explicit identifiers
    # to src_axes for the cf version you have seems to be the way.
    # src_axes=['latitude','longitude']?
    field_REGRID = field_CF.regrids(g, dst_cyclic=True, method='nearest_dtos',
                                    src_axes={'X': 'key%X', 'Y': 'key%Y'})

    return field_REGRID


def reroute_calving_to_ocean(calv_cf, route_file, coupling_period,
                             ice_tstep_multiply, orig_calving, output_file):

    """
    Coastlines of the BISICLES and orca grids to do not match perfectly, so
    some of the remapped BISICLES calving will be located on orca land cells.
    Use a precalculated map to warp all our remapped BISICLES data to the nearest 
    open ocean point
    """

    import cf

    # Extract the numpy array, with 0s instead of missing_data
    calv_ma = calv_cf.array
    np.ma.set_fill_value(calv_ma, 0)
    calv = np.ma.filled(calv_ma)

    # This is m3 per coupling interval. We want km^3/yr for NEMO.
    # THIS MAY CHANGE IN FUTURE RELEASES?
    # Assume 360_day calendar
    calv = calv / 1e9 / coupling_period * secs_in_day * days_in_year

    # If ice has been accelerated with asynchronous coupling to climate model, 
    # pull back the iceberg magnitude given to NEMO to a characteristic per year
    # value for the next 1 year of climate
    calv = calv / ice_tstep_multiply

    # Set up the routing map
    mapping = cf.read(route_file)
    try:
        iindex = mapping.select_field("ncvar%xindex").array
    except BaseException:
        iindex = mapping[1].array
    try:
        jindex = mapping.select_field("ncvar%yindex").array
    except BaseException:
        jindex = mapping[0].array

    nj = np.shape(iindex)[0]
    ni = np.shape(iindex)[1]

    # Loop over the routed area and move everything to its destination
    calv_routed = np.zeros([nj, ni])
    for j in np.arange(nj):
        for i in np.arange(ni):

            iiindex = int(iindex[j, i])
            ijindex = int(jindex[j, i])

            if (ijindex == 0) & (iiindex == 0):
                # If routfile is not appropriate
                calv_routed[j, i] = calv[j, i]
            else:
                calv_routed[ijindex, iiindex] \
                    = calv_routed[ijindex, iiindex] + calv[j, i]

    # This is just the regional field - will construct a unified global field
    # elsewhere
    h = cf.read(orig_calving)

    original_seed = h.select_field("ncvar%calvingmask")

    # At NEMO4.2 the halos have been stripped from calving reference file, but
    # not from calving mask. Hence, the use of i_start, i_end and j_end below
    # to strip the halos from calv_routed. Keep an eye on this convention as 
    # UKESM configs with NEMO4.2 evolve.
    new_seed = np.copy(original_seed.array)
    new_seed_shape = np.shape(new_seed)
    j_end = new_seed_shape[1]
    if new_seed_shape[2] == ni:
        i_start = 0
    else:
        i_start = 1
    i_end = new_seed_shape[2] + i_start
    for k in range(new_seed_shape[0]):
        new_seed[k, :, :] = calv_routed[0:j_end, i_start:i_end]

    original_seed.set_data(cf.Data(new_seed))

    # write the final, rerouted iceberg info out to the file NEMO will use.
    cf.write(h, output_file, fmt="NETCDF4")

    return


if __name__ == "__main__":

    # Read the command line arguments
    hdf5in_file, hdf5out_file, ncout_file, coupling_period, \
        ice_tstep_multiply, region = parse_commandline()

    # Get what we need from the BISICLES plot file
    print("Extract fields from BISICLES hdf5")
    calving, lon, lat, x_dist, y_dist \
        = extract_bisicles_calving(hdf5in_file, region=region)

    # Make calving source/sink files for NEMO/BISICLES. These grid/template
    # files are linked to appropriate sources by the suite
    print("Deal with the iceberg/melange files")

    nemo_ncgridfile = "cf_gridfile.nc"
    route_file = "cf_routfile.nc"
    orig_calving = "orig_calving.nc"

    # make sink file for BISICLES
    make_BISICLES_melange_sink(calving, x_dist, y_dist, coupling_period,
                               hdf5out_file)

    # make a CF Field from the calving info
    calving_CF = construct_CF(calving, lon, lat, x_dist, y_dist)

    # remap the CF Field to orca grid
    calving_REGRID = regrid_to_NEMO(calving_CF, nemo_ncgridfile)

    # warp all the icebergs to a wet cell in NEMO and make the source file for NEMO
    reroute_calving_to_ocean(calving_REGRID, route_file, coupling_period,
                             ice_tstep_multiply, orig_calving, ncout_file)
