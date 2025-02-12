#!/usr/bin/env python

"""
um_to_unicicles

The annual mean pp output from the UM is read and a netCDF file for input
to UniCiCles is made. Mostly just juggling of variable names and units.
of HadGEM without ice sheet tiles or the multilayer snowpack model

Called by suite as atmos_to_ice_sheets
"""

import sys
import os
import numpy as np
import cf
from common_arg_to_file_exist import arg_to_file_exist
from common_params_and_constants import days_in_year, deg_to_rad, \
    degC_to_Kelvin, max_number_elevs, number_elevs, planet_radius, \
    rho_ice_glimmer, secs_in_day, start_tile_id_elev_ice, \
    start_tile_id_elev_rock, stashcode_land_frac_diag, stashcode_orog, \
    stashcode_tile_fractions, stashcode_tile_smb, stashcode_tile_snow_htflux, \
    stashcode_tile_snow_mass, stashcode_tile_snow_temp, tile_elevations


def parse_commandline():
    """ 
    Read the line command line arguments 
    "snow_shuffling" tries to make adjustments to the UM snowpacks away from the ice sheet
    so that accumulation can incept new solid ice points and maintain water mass conservation 
    as ice sheet tile area changes. This has never proven technically robust in practice
    "exclude code" lets the routine be fed with archived UM output that doesn't contain
    all the stashcodes the ice sheet coupling expects, but enough to still make an SMB forcing
    """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="input, name of UM forcing data")
    parser.add_argument("--output", help="output, name of input to uncicles")
    parser.add_argument("--do_snow_shuffling",
                        help="experimental. If True, will pass UM non-ice "
                        "snow to UniCiCles to try to incept new ice points")
    parser.add_argument("--exclude_code", action="append",
                        help="don't expect these fields in input - output "
                        "them as 0")
    args = parser.parse_args()

    err = 0
    input, err = arg_to_file_exist(args.input, mandatory=True, err=err)
    output, err = arg_to_file_exist(args.output, io="out", err=err)

    do_snow_shuffling = False
    if args.do_snow_shuffling is not None:
        if args.do_snow_shuffling == "True":
            do_snow_shuffling = True

    exclude_codes = []
    if isinstance(args.exclude_code, type([])):
        exclude_codes = list(map(int, args.exclude_code))

    if err > 0:
        parser.print_help()
        sys.exit(err)

    return input, output, exclude_codes, do_snow_shuffling

def extract_fields(all_fields, extract_codes, exclude_codes):

    """
    get the fields we need to pass to UniCiCles from the UM output, allowing for
    some non-essential ones to be supplied empty if they aren't available
    """

    extract_list = cf.FieldList()

    for code in extract_codes:
        read_code = code

        print(code, "in extract_fields, excluding", exclude_codes)

        if code in exclude_codes:
            # For now, assume we're only faking fields dimensioned
            # on elevated tiles. Read in the tile fraction field
            # for that purpose
            print("faking a tile field for", code)
            read_code = stashcode_tile_fractions

        select_string = 'stash_code=' + str(read_code)
        select_field = all_fields.select_field(select_string).copy()

        # make up the necessary metadata (and clear the actual data) for a faked field
        if code in exclude_codes:
            select_field.set_property("stash_code", str(code))
            select_field.del_property("long_name")
            empty = np.zeros_like(select_field.array)
            select_field.set_data(cf.Data(empty))

        print("extracting ", code)

        # change units for certain fields
        if code == stashcode_tile_snow_mass:
            select_field = select_field / rho_ice_glimmer
        if code == stashcode_tile_snow_temp:
            select_field = select_field + degC_to_Kelvin

        #we want EITHER ice_ OR rock_ elevated tile for some fields
        #subspace to just those levels
        if select_field.shape[0] < max_number_elevs * 2:
            if select_field.shape[0] == 9:
            # Elevated tiles not on in source UM.
                print("9 tile config (ice and bare soil) dected, can't make coupling field")
                exit(1)
            else:
                print("subspace elevated pseudo levels ", code)
                if (code == stashcode_tile_snow_mass) & \
                   (read_code == stashcode_tile_snow_mass):
                    #we want the NON-ICE snowpack mass
                    select_field = select_field.subspace(
                        **{'long_name=pseudolevel':
                           cf.ge(start_tile_id_elev_rock)})
                else:
                    #we want the ice tile quantities of everything else
                    select_field = select_field.subspace(
                        **{'long_name=pseudolevel':
                           cf.ge(start_tile_id_elev_ice)})
                    select_field = select_field.subspace(
                        **{'long_name=pseudolevel':
                           cf.lt(start_tile_id_elev_rock)})

        extract_list.append(select_field)

    return extract_list


def calculate_area(field_list):

    """
    calculate areas of active elevated tiles - the coupler can use this for conservation
    """

    #get heights (from the centre of the Earth) and tile and coastal area fractions
    select_string = 'stash_code=' + str(stashcode_orog)
    r_theta = field_list.select_field(select_string) + planet_radius
    select_string = 'stash_code=' + str(stashcode_land_frac_diag)
    coastal_fraction = field_list.select_field(select_string)
    select_string = 'stash_code=' + str(stashcode_tile_fractions)
    tile_fraction = field_list.select_field(select_string)

    #get latitudes and longitudes
    latitude_1d = r_theta.coord('latitude').array * deg_to_rad
    longitude_1d = r_theta.coord('longitude').array * deg_to_rad
    delta_phi = latitude_1d[1] - latitude_1d[0]
    delta_lambda = longitude_1d[1] - longitude_1d[0]
    latitude_2d = np.zeros([latitude_1d.size, longitude_1d.size])
    for j in range(latitude_1d.size):
        latitude_2d[j, :] = latitude_1d[j]

    cos_theta_latitude = r_theta.copy()
    cos_theta_latitude.set_data(cf.Data(np.cos(latitude_2d)))

    # Calculate gridbox areas as in um/src/control/coupling/ice_sheet_mass.F90
    cell_area = r_theta * r_theta *                  \
                delta_lambda * delta_phi *           \
                cos_theta_latitude *                 \
                coastal_fraction

    # Then areas of each tile in the gridboxes
    tile_area_np = tile_fraction.array
    for k in range(tile_area_np.shape[0]):
        tile_area_np[k, :, :] = tile_area_np[k,:,:] * cell_area.array

    # write back out as a cf Field
    tile_area = tile_fraction.copy()
    tile_area.set_data(cf.Data(tile_area_np))

    tile_area.set_property('long_name', 'tile surface area')
    tile_area.id = 'tile_surface_area'
    tile_area.nc_set_variable('tile_surface_area')
    for attr in ['standard_name', 'stash_code', 'um_stash_source']:
        if tile_area.has_property(attr):
            tile_area.del_property(attr)

    return tile_area


def change_names(out_list):

    # change certain field names and metadata so they will be correctly identified 
    # by the UniCiCles wrapper
    select_string = 'stash_code=' + str(stashcode_tile_snow_mass)
    field = out_list.select_field(select_string)
    field.id = 'nonice_snowdepth'
    field.nc_set_variable('nonice_snowdepth')
    field.coord('long_name=pseudolevel').set_property('standard_name',
                                                      'tile_id')

    select_string = 'stash_code=' + str(stashcode_tile_snow_temp)
    field = out_list.select_field(select_string)
    field.id = 'ice_stemp'
    field.nc_set_variable('ice_stemp')
    field.coord('long_name=pseudolevel').set_property('standard_name',
                                                      'tile_id')

    select_string = 'stash_code=' + str(stashcode_tile_smb)
    field = out_list.select_field(select_string)
    field.id = 'ice_smb'
    field.nc_set_variable('ice_smb')
    field.coord('long_name=pseudolevel').set_property('standard_name',
                                                      'tile_id')

    select_string = 'stash_code=' + str(stashcode_tile_snow_htflux)
    field = out_list.select_field(select_string)
    field.id = 'snow_ice_hflux'
    field.nc_set_variable('snow_ice_hflux')
    field.coord('long_name=pseudolevel').set_property('standard_name',
                                                      'tile_id')

    select_string = 'stash_code=' + str(stashcode_tile_fractions)
    field = out_list.select_field(select_string)
    field.id = 'tile_frac'
    field.nc_set_variable('tile_frac')
    field.coord('long_name=pseudolevel').set_property('standard_name',
                                                      'tile_id')

    field = out_list.select_field('long_name=tile surface area')
    field.coord('long_name=pseudolevel').set_property('standard_name',
                                                      'tile_id')

    return out_list

if __name__ == "__main__":

    # Parse command line
    file_in, file_out, exclude_codes, do_snow_shuffling \
        = parse_commandline()

    # Extract the fields we want
    all_fields = cf.read(file_in)

    # We need these fields to calculate compound quantities
    process_codes = [stashcode_orog, stashcode_land_frac_diag,
                     stashcode_tile_fractions, stashcode_tile_smb,
                     stashcode_tile_snow_mass]

    process_list = extract_fields(all_fields, process_codes, exclude_codes)

    # These fields wil just be copied over as they are
    echo_codes = [stashcode_tile_fractions, stashcode_tile_snow_temp,
                  stashcode_tile_snow_htflux]

    echo_list = extract_fields(all_fields, echo_codes, exclude_codes)

    # Set up the new list of fields to output
    out_list = cf.FieldList()

    # Create tile areas
    tile_area = calculate_area(process_list)
    out_list.append(tile_area)

    # Copy SMB to output
    select_string = 'stash_code=' + str(stashcode_tile_smb)
    ice_smb = process_list.select_field(select_string)

    out_list.append(ice_smb)

    # Non-ice snow - check if the coupler is trying to shuffle ni_snow around?
    select_string = 'stash_code=' + str(stashcode_tile_snow_mass)
    nisnow = process_list.select_field(select_string)

    if not do_snow_shuffling:
        print("SNOW SHUFFLING OFF")
        print("ZEROING NISNOW PASSED INTO BISICLES")
        print("FOR CONSISTENCY WITH NOT ATTEMPTING TO REARANGE IT IN THE ",
              "UM AFTERWARDS")
        print("(see unicicles_cap_to_um)")
        nisnow.data = nisnow.data -nisnow.data

    out_list.append(nisnow)

    # Add the rest of the fields to the output
    for field in echo_list:
        out_list.append(field)

    # Change netCDF_out names to human-readable, common format with JULES
    out_list = change_names(out_list)

    # Write new file, overwriting the nisnow stuff if necessary
    cf.write(out_list, file_out, fmt = "NETCDF4", overwrite=True)
