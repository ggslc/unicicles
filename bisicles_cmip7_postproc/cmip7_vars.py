"""
CMIP7/CF variable definitions and BISICLES-to-CMIP7 field mapping.

Provides metadata tables for:
  - 2D spatial field mapping (BISICLES name -> CMIP7 name, CF attributes, unit conversions)
  - Derived spatial fields (computed from mask + iceFrac)
  - Scalar diagnostic mapping (diagnostics CSV -> CMIP7 timeseries variables)

Variable names and metadata are taken from the CMIP7 CMOR landIce table
(table_id: "landIce", table_date: 2026-03-08, Conventions: CF-1.12 CMIP-7.0).

Fields marked ``cmip7_compliant: False`` have no equivalent ``out_name`` in
the CMIP7 landIce table but are retained because they are useful diagnostics
produceable from BISICLES output.

BISICLES mask values (from IceConstants.H):
  GROUNDEDMASKVAL = 1
  FLOATINGMASKVAL = 2
  OPENSEAMASKVAL  = 4
  OPENLANDMASKVAL = 8
"""

# Physical constants used for unit conversions
ICE_DENSITY = 918.0        # kg m-3
WATER_DENSITY = 1028.0     # kg m-3
GRAVITY = 9.81             # m s-2
SECS_PER_YEAR = 31556926.0  # s a-1  (SECONDS_PER_TROPICAL_YEAR from BISICLES IceConstants.H)

# BISICLES mask integer values
GROUNDED_MASK_VAL = 1
FLOATING_MASK_VAL = 2
OPEN_SEA_MASK_VAL = 4
OPEN_LAND_MASK_VAL = 8

# Conversion factors from BISICLES units to CMIP7/SI units
_M_A_TO_M_S = 1.0 / SECS_PER_YEAR                    # m a-1        -> m s-1
_M_A_ICE_TO_KG_M2_S = ICE_DENSITY / SECS_PER_YEAR    # m a-1 (ice)  -> kg m-2 s-1
_KG_M2_A_TO_KG_M2_S = 1.0 / SECS_PER_YEAR            # kg m-2 a-1   -> kg m-2 s-1


# ---------------------------------------------------------------------------
# 2D spatial field mapping
#
# Keys are BISICLES internal field names as they appear in plot HDF5 files.
# "conversion_factor" multiplies the raw BISICLES value to reach cmip7_units.
# "cell_methods" follows the CF convention and matches the CMIP7 table entry
#   for snapshot (time-independent) output; temporal averaging is not applied
#   by this code, which processes individual plot files.
# "cmip7_compliant: False" marks fields with no CMIP7 landIce out_name.
# ---------------------------------------------------------------------------
FIELD_MAPPING = {
    # ------------------------------------------------------------------
    # CMIP7-compliant fields
    # ------------------------------------------------------------------
    "thickness": {
        "cmip7_name": "lithk",
        "standard_name": "land_ice_thickness",
        "long_name": "Ice Sheet Thickness",
        "cmip7_units": "m",
        "bisicles_units": "m",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },
    "Z_surface": {
        "cmip7_name": "orog",
        "standard_name": "surface_altitude",
        "long_name": "Surface Altitude",
        "cmip7_units": "m",
        "bisicles_units": "m",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },
    "Z_base": {
        "cmip7_name": "topg",
        "standard_name": "bedrock_altitude",
        "long_name": "Bedrock Altitude",
        "cmip7_units": "m",
        "bisicles_units": "m",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean where grounded_ice_sheet (mask=sfgrlf)",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },
    "activeSurfaceThicknessSource": {
        "cmip7_name": "acabf",
        "standard_name": "land_ice_surface_specific_mass_balance_flux",
        "long_name": "Ice Sheet Surface Mass Balance Flux",
        "cmip7_units": "kg m-2 s-1",
        "bisicles_units": "m a-1",
        "conversion_factor": _M_A_ICE_TO_KG_M2_S,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
        "comment": "Positive values indicate accumulation (gain of ice mass).",
    },
    "activeBasalThicknessSource": {
        "cmip7_name": "libmassbf",
        "standard_name": "land_ice_basal_specific_mass_balance_flux",
        "long_name": "Land Ice Basal Specific Mass Balance Flux",
        "cmip7_units": "kg m-2 s-1",
        "bisicles_units": "m a-1",
        "conversion_factor": _M_A_ICE_TO_KG_M2_S,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
        "comment": (
            "Combined grounded and floating basal mass balance. "
            "The CMIP7 table defines separate grounded (mask=sftgrf) and "
            "floating (mask=sftflf) variants; this field covers both. "
            "Positive values indicate accumulation (freezing)."
        ),
    },
    "calvingFlux": {
        "cmip7_name": "licalvf",
        "standard_name": "land_ice_specific_mass_flux_due_to_calving",
        "long_name": "Land Ice Calving Flux",
        "cmip7_units": "kg m-2 s-1",
        "bisicles_units": "m a-1",
        "conversion_factor": _M_A_ICE_TO_KG_M2_S,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
        "comment": "Positive values indicate calving (loss of ice mass).",
    },
    "dragCoef": {
        "cmip7_name": "strbasemag",
        "standard_name": "land_ice_basal_drag",
        "long_name": "Land Ice Basal Drag",
        "cmip7_units": "Pa",
        "bisicles_units": "Pa",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },
    # xVel / yVel are depth-averaged (columnar mean) velocities in BISICLES
    # (SSA / L1L2 approximation).  The closest CMIP7 landIce variables are
    # xvelmean / yvelmean (land_ice_vertical_mean_x/y_velocity).
    # Note: xvelsurf / xvelbase are distinct CMIP7 variables (surface/basal)
    # that require a multi-layer model.
    "xVel": {
        "cmip7_name": "xvelmean",
        "standard_name": "land_ice_vertical_mean_x_velocity",
        "long_name": "X-Component of Land Ice Vertical Mean Velocity",
        "cmip7_units": "m s-1",
        "bisicles_units": "m a-1",
        "conversion_factor": _M_A_TO_M_S,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },
    "yVel": {
        "cmip7_name": "yvelmean",
        "standard_name": "land_ice_vertical_mean_y_velocity",
        "long_name": "Y-Component of Land Ice Vertical Mean Velocity",
        "cmip7_units": "m s-1",
        "bisicles_units": "m a-1",
        "conversion_factor": _M_A_TO_M_S,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },

    # ------------------------------------------------------------------
    # Non-CMIP7 fields (useful BISICLES diagnostics with no CMIP7 out_name
    # in the landIce table; retained for scientific completeness)
    # ------------------------------------------------------------------

    # iceFrac (total land ice area fraction) has no CMIP7 landIce out_name.
    # The CMIP7 table only defines sftgrf (grounded %) and sftflf (floating %)
    # as derived fields.  iceFrac is retained here under the CF standard_name
    # land_ice_area_fraction with units "1" (fraction, not percentage).
    "iceFrac": {
        "cmip7_name": "sftgif",
        "standard_name": "land_ice_area_fraction",
        "long_name": "Land Ice Area Fraction",
        "cmip7_units": "1",
        "bisicles_units": "1",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean",
        "modeling_realm": "landIce",
        "cmip7_compliant": False,
        "comment": (
            "Total land ice area fraction (grounded + floating). "
            "Not a CMIP7 landIce out_name; sftgrf and sftflf are the CMIP7 "
            "equivalents expressed as percentages."
        ),
    },
    # land_ice_base_elevation is not in the CMIP7 landIce table.
    # Z_bottom = z_surface - thickness (bottom surface of the ice column).
    "Z_bottom": {
        "cmip7_name": "base",
        "standard_name": "land_ice_base_elevation",
        "long_name": "Land Ice Base Elevation",
        "cmip7_units": "m",
        "bisicles_units": "m",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": False,
        "comment": "Not a CMIP7 landIce out_name; retained as a useful diagnostic.",
    },
    # tendency_of_land_ice_thickness is not in the CMIP7 landIce table.
    "dThickness/dt": {
        "cmip7_name": "dlithkdt",
        "standard_name": "tendency_of_land_ice_thickness",
        "long_name": "Tendency of Land Ice Thickness",
        "cmip7_units": "m s-1",
        "bisicles_units": "m a-1",
        "conversion_factor": _M_A_TO_M_S,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": False,
        "comment": "Not a CMIP7 landIce out_name; retained as a useful diagnostic.",
    },
    # Flux divergence — no CMIP7 equivalent.
    "divergenceThicknessFlux": {
        "cmip7_name": "divflux",
        "standard_name": "tendency_of_land_ice_thickness_due_to_flow",
        "long_name": "Divergence of Ice Thickness Flux",
        "cmip7_units": "m s-1",
        "bisicles_units": "m a-1",
        "conversion_factor": _M_A_TO_M_S,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": False,
        "comment": "Not a CMIP7 landIce out_name; retained as a useful diagnostic.",
    },
    # Till water depth — no CMIP7 equivalent.
    "tillWaterDepth": {
        "cmip7_name": "tillwatd",
        "standard_name": "till_water_thickness",
        "long_name": "Till Water Depth",
        "cmip7_units": "m",
        "bisicles_units": "m",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean where grounded_ice_sheet (mask=sfgrlf)",
        "modeling_realm": "landIce",
        "cmip7_compliant": False,
        "comment": "Not a CMIP7 landIce out_name; retained as a useful diagnostic.",
    },
    # Melange thickness — no CMIP7 equivalent.
    "melangeThickness": {
        "cmip7_name": "melangeThickness",
        "standard_name": "melange_thickness",
        "long_name": "Melange Thickness",
        "cmip7_units": "m",
        "bisicles_units": "m",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": False,
        "comment": (
            "Thickness of ice mélange at the calving front. "
            "Not a CMIP7 landIce out_name; retained as a useful diagnostic."
        ),
    },
}

# Reverse lookup: CMIP7 name -> BISICLES name
CMIP7_TO_BISICLES = {v["cmip7_name"]: k for k, v in FIELD_MAPPING.items()}


# ---------------------------------------------------------------------------
# CF-plot file field mapping
#
# When BISICLES is configured to write CF-output files (plot.CF-*.hdf5), it
# uses CMIP7-compatible variable names directly (lithk, acabf, sftgrf, …)
# rather than its internal names (thickness, activeSurfaceThicknessSource, …).
# The metadata written by AmrIceIO.cpp is incomplete for CMIP7 compliance —
# units are in per-year rather than per-second, cell_methods, modeling_realm
# and other required attributes are absent.
#
# This table maps each BISICLES CF output name to the full CMIP7/CF metadata
# and the conversion factor needed to reach SI units.
#
# BISICLES CF-output units (as written by AmrIceIO.cpp):
#   velocity fields:     m yr-1
#   mass-balance fluxes: kg m-2 yr-1
#   temperature:         K
#   thickness/elevation: m
#   area fractions:      1  (dimensionless fraction, 0–1)
# ---------------------------------------------------------------------------
CF_FIELD_MAPPING = {
    # ------------------------------------------------------------------
    # CMIP7-compliant fields
    # ------------------------------------------------------------------
    "lithk": {
        "cmip7_name": "lithk",
        "standard_name": "land_ice_thickness",
        "long_name": "Ice Sheet Thickness",
        "cmip7_units": "m",
        "bisicles_units": "m",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },
    "orog": {
        "cmip7_name": "orog",
        "standard_name": "surface_altitude",
        "long_name": "Surface Altitude",
        "cmip7_units": "m",
        "bisicles_units": "m",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },
    "topg": {
        "cmip7_name": "topg",
        "standard_name": "bedrock_altitude",
        "long_name": "Bedrock Altitude",
        "cmip7_units": "m",
        "bisicles_units": "m",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean where grounded_ice_sheet (mask=sfgrlf)",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },
    "xvelbase": {
        "cmip7_name": "xvelbase",
        "standard_name": "land_ice_basal_x_velocity",
        "long_name": "X-Component of Land Ice Basal Velocity",
        "cmip7_units": "m s-1",
        "bisicles_units": "m yr-1",
        "conversion_factor": _M_A_TO_M_S,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },
    "yvelbase": {
        "cmip7_name": "yvelbase",
        "standard_name": "land_ice_basal_y_velocity",
        "long_name": "Y-Component of Land Ice Basal Velocity",
        "cmip7_units": "m s-1",
        "bisicles_units": "m yr-1",
        "conversion_factor": _M_A_TO_M_S,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },
    "xvelsurf": {
        "cmip7_name": "xvelsurf",
        "standard_name": "land_ice_surface_x_velocity",
        "long_name": "X-Component of Land Ice Surface Velocity",
        "cmip7_units": "m s-1",
        "bisicles_units": "m yr-1",
        "conversion_factor": _M_A_TO_M_S,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },
    "yvelsurf": {
        "cmip7_name": "yvelsurf",
        "standard_name": "land_ice_surface_y_velocity",
        "long_name": "Y-Component of Land Ice Surface Velocity",
        "cmip7_units": "m s-1",
        "bisicles_units": "m yr-1",
        "conversion_factor": _M_A_TO_M_S,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },
    "litemptop": {
        "cmip7_name": "litemptop",
        "standard_name": "temperature_at_top_of_ice_sheet_model",
        "long_name": "Temperature at Top of Ice Sheet Model",
        "cmip7_units": "K",
        "bisicles_units": "K",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },
    "litempbot": {
        "cmip7_name": "litempbot",
        "standard_name": "temperature_at_base_of_ice_sheet_model",
        "long_name": "Temperature at Base of Ice Sheet Model",
        "cmip7_units": "K",
        "bisicles_units": "K",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },
    # CMIP7 requires sftgrf and sftflf in % (0–100); BISICLES CF output writes
    # them as dimensionless fractions (0–1), so conversion_factor = 100.
    "sftgrf": {
        "cmip7_name": "sftgrf",
        "standard_name": "grounded_ice_sheet_area_fraction",
        "long_name": "Grounded Ice Sheet Area Percentage",
        "cmip7_units": "%",
        "bisicles_units": "1",
        "conversion_factor": 100.0,
        "cell_methods": "area: time: mean",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
        "comment": (
            "Percentage of cell area covered by grounded ice. "
            "Directly output by BISICLES (not derived from flotation criterion)."
        ),
    },
    "sftflf": {
        "cmip7_name": "sftflf",
        "standard_name": "floating_ice_shelf_area_fraction",
        "long_name": "Floating Ice Shelf Area Percentage",
        "cmip7_units": "%",
        "bisicles_units": "1",
        "conversion_factor": 100.0,
        "cell_methods": "area: time: mean",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
        "comment": (
            "Percentage of cell area covered by floating ice. "
            "Directly output by BISICLES (not derived from flotation criterion)."
        ),
    },
    "acabf": {
        "cmip7_name": "acabf",
        "standard_name": "land_ice_surface_specific_mass_balance_flux",
        "long_name": "Ice Sheet Surface Mass Balance Flux",
        "cmip7_units": "kg m-2 s-1",
        "bisicles_units": "kg m-2 yr-1",
        "conversion_factor": _KG_M2_A_TO_KG_M2_S,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
        "comment": "Positive values indicate accumulation (gain of ice mass).",
    },
    "libmassbf": {
        "cmip7_name": "libmassbf",
        "standard_name": "land_ice_basal_specific_mass_balance_flux",
        "long_name": "Land Ice Basal Specific Mass Balance Flux",
        "cmip7_units": "kg m-2 s-1",
        "bisicles_units": "kg m-2 yr-1",
        "conversion_factor": _KG_M2_A_TO_KG_M2_S,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
        "comment": "Positive values indicate basal accumulation (freezing).",
    },
    "licalvf": {
        "cmip7_name": "licalvf",
        "standard_name": "land_ice_specific_mass_flux_due_to_calving",
        "long_name": "Land Ice Calving Flux",
        "cmip7_units": "kg m-2 s-1",
        "bisicles_units": "kg m-2 yr-1",
        "conversion_factor": _KG_M2_A_TO_KG_M2_S,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
        "comment": "Positive values indicate calving (loss of ice mass).",
    },
    "strbasemag": {
        "cmip7_name": "strbasemag",
        "standard_name": "land_ice_basal_drag",
        "long_name": "Land Ice Basal Drag",
        "cmip7_units": "Pa",
        "bisicles_units": "Pa",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },
    "hfgeoubed": {
        "cmip7_name": "hfgeoubed",
        "standard_name": "upward_geothermal_heat_flux_at_ground_level_in_land_ice",
        "long_name": "Geothermal Heat Flux at Ice Sheet Base",
        "cmip7_units": "W m-2",
        "bisicles_units": "W m-2",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean where grounded_ice_sheet (mask=sfgrlf)",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },

    # ------------------------------------------------------------------
    # Non-CMIP7 fields present in BISICLES CF output
    # Retained because they are useful diagnostics.
    # ------------------------------------------------------------------

    # sftgif (total land ice area fraction) has no CMIP7 landIce out_name.
    "sftgif": {
        "cmip7_name": "sftgif",
        "standard_name": "land_ice_area_fraction",
        "long_name": "Land Ice Area Fraction",
        "cmip7_units": "1",
        "bisicles_units": "1",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean",
        "modeling_realm": "landIce",
        "cmip7_compliant": False,
        "comment": (
            "Total land ice area fraction (grounded + floating). "
            "Not a CMIP7 landIce out_name; sftgrf and sftflf are the CMIP7 "
            "equivalents expressed as percentages."
        ),
    },
    # Grounded/floating basal mass balance splits — retain but mark non-CMIP7.
    "libmassbfgr": {
        "cmip7_name": "libmassbfgr",
        "standard_name": "land_ice_basal_specific_mass_balance_flux",
        "long_name": "Basal Mass Balance Flux of Grounded Ice Sheet",
        "cmip7_units": "kg m-2 s-1",
        "bisicles_units": "kg m-2 yr-1",
        "conversion_factor": _KG_M2_A_TO_KG_M2_S,
        "cell_methods": "area: time: mean where grounded_ice_sheet (mask=sfgrlf)",
        "modeling_realm": "landIce",
        "cmip7_compliant": False,
        "comment": "Not a CMIP7 landIce out_name; retained as a useful diagnostic.",
    },
    "libmassbffl": {
        "cmip7_name": "libmassbffl",
        "standard_name": "land_ice_basal_specific_mass_balance_flux",
        "long_name": "Basal Mass Balance Flux of Floating Ice Shelf",
        "cmip7_units": "kg m-2 s-1",
        "bisicles_units": "kg m-2 yr-1",
        "conversion_factor": _KG_M2_A_TO_KG_M2_S,
        "cell_methods": "area: time: mean where floating_ice_shelf (mask=sftflf)",
        "modeling_realm": "landIce",
        "cmip7_compliant": False,
        "comment": "Not a CMIP7 landIce out_name; retained as a useful diagnostic.",
    },
    # dlithkdt — not in CMIP7 landIce table.
    "dlithkdt": {
        "cmip7_name": "dlithkdt",
        "standard_name": "tendency_of_land_ice_thickness",
        "long_name": "Tendency of Land Ice Thickness",
        "cmip7_units": "m s-1",
        "bisicles_units": "m yr-1",
        "conversion_factor": _M_A_TO_M_S,
        "cell_methods": "area: time: mean where ice_sheet",
        "modeling_realm": "landIce",
        "cmip7_compliant": False,
        "comment": "Not a CMIP7 landIce out_name; retained as a useful diagnostic.",
    },
    # Basal temperature splits — retain but mark non-CMIP7.
    "litempbotgr": {
        "cmip7_name": "litempbotgr",
        "standard_name": "temperature_at_base_of_ice_sheet_model",
        "long_name": "Temperature at Base of Grounded Ice Sheet",
        "cmip7_units": "K",
        "bisicles_units": "K",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean where grounded_ice_sheet (mask=sfgrlf)",
        "modeling_realm": "landIce",
        "cmip7_compliant": False,
        "comment": "Not a CMIP7 landIce out_name; retained as a useful diagnostic.",
    },
    "litempbotfl": {
        "cmip7_name": "litempbotfl",
        "standard_name": "temperature_at_base_of_ice_sheet_model",
        "long_name": "Temperature at Base of Floating Ice Shelf",
        "cmip7_units": "K",
        "bisicles_units": "K",
        "conversion_factor": 1.0,
        "cell_methods": "area: time: mean where floating_ice_shelf (mask=sftflf)",
        "modeling_realm": "landIce",
        "cmip7_compliant": False,
        "comment": "Not a CMIP7 landIce out_name; retained as a useful diagnostic.",
    },
}


# ---------------------------------------------------------------------------
# Derived 2D spatial fields
#
# Computed in Python from combinations of BISICLES fields.
# "conversion_factor" is applied when writing to NetCDF; internal computation
# always uses fraction (0–1) values.
#
# CMIP7 landIce table specifies sftgrf and sftflf in units of "%" (percentage),
# i.e. values 0–100.  The conversion_factor of 100 converts from the
# internally-computed fraction (0–1) to the required percentage.
# ---------------------------------------------------------------------------
DERIVED_FIELDS = {
    "sftgrf": {
        "standard_name": "grounded_ice_sheet_area_fraction",
        "long_name": "Grounded Ice Sheet Area Percentage",
        "cmip7_units": "%",
        "conversion_factor": 100.0,   # fraction (0–1) -> percentage (0–100)
        "cell_methods": "area: time: mean",
        "modeling_realm": "landIce",
        "source_fields": ["iceFrac", "mask"],
        "comment": (
            "Percentage of cell area covered by grounded ice. "
            "Computed as 100 * iceFrac where mask == GROUNDEDMASKVAL (1), else 0. "
            "iceFrac and/or mask may be derived from geometry if absent from the "
            "plot file (see flatten.compute_flotation_mask)."
        ),
    },
    "sftflf": {
        "standard_name": "floating_ice_shelf_area_fraction",
        "long_name": "Floating Ice Shelf Area Percentage",
        "cmip7_units": "%",
        "conversion_factor": 100.0,   # fraction (0–1) -> percentage (0–100)
        "cell_methods": "area: time: mean",
        "modeling_realm": "landIce",
        "source_fields": ["iceFrac", "mask"],
        "comment": (
            "Percentage of cell area covered by floating ice. "
            "Computed as 100 * iceFrac where mask == FLOATINGMASKVAL (2), else 0. "
            "iceFrac and/or mask may be derived from geometry if absent from the "
            "plot file (see flatten.compute_flotation_mask)."
        ),
    },
}


# ---------------------------------------------------------------------------
# Grid-geometry fields
#
# These are not read from BISICLES plot files but computed from the output
# grid coordinates.  Each is a static (time-invariant) spatial field.
# ---------------------------------------------------------------------------
GRID_GEOMETRY_FIELDS = {
    "modelcellareai": {
        "standard_name": "cell_area",
        "long_name": "Model Grid-Cell Area for Ice-Sheet Variables",
        "cmip7_units": "m2",
        "cell_methods": "area: sum",
        "modeling_realm": "landIce",
        "cmip7_compliant": True,
    },
}


# ---------------------------------------------------------------------------
# Scalar diagnostic variable mapping
#
# Maps (region, quantity) pairs from the diagnostics CSV output to CMIP7
# scalar timeseries variable names.
#
# diagnostics CSV units:
#   volume quantities: m3
#   area quantities:   m2
#   flux quantities:   m3 a-1  (volume of ice per year)
#
# "conversion_factor" converts from the diagnostics units to cmip7_units.
# "cmip7_compliant: False" marks entries with no CMIP7 landIce out_name.
# ---------------------------------------------------------------------------
SCALAR_MAPPING = {
    # ------------------------------------------------------------------
    # CMIP7-compliant scalar fields
    # ------------------------------------------------------------------

    # lim — land_ice_mass — CMIP7: lim_tavg-u-hm-is
    ("entire", "volume"): {
        "cmip7_name": "lim",
        "standard_name": "land_ice_mass",
        "long_name": "Ice Sheet Mass",
        "cmip7_units": "kg",
        "conversion_factor": ICE_DENSITY,   # m3 * kg m-3 = kg
        "cell_methods": "area: sum where ice_sheet time: mean",
        "cmip7_compliant": True,
        "ismip7_type":"ST",
        "comment": "Total volume of land ice multiplied by ice density.",
    },

    # limnsw — land_ice_mass_not_displacing_sea_water — CMIP7: limnsw_tavg-u-hm-is
    ("entire", "volumeAbove"): {
        "cmip7_name": "limnsw",
        "standard_name": "land_ice_mass_not_displacing_sea_water",
        "long_name": "Ice Sheet Mass That Does not Displace Sea Water",
        "cmip7_units": "kg",
        "conversion_factor": ICE_DENSITY,
        "cell_methods": "area: sum where ice_sheet time: mean",
        "cmip7_compliant": True,
        "ismip7_type":"ST",
        "comment": "Ice mass above flotation, relevant to sea-level change.",
    },

    # iareagr — grounded_ice_sheet_area — CMIP7: iareagr_tavg-u-hm-gis
    ("grounded", "area"): {
        "cmip7_name": "iareagr",
        "standard_name": "grounded_ice_sheet_area",
        "long_name": "Area Covered by Grounded Ice Sheet",
        "cmip7_units": "m^2",
        "conversion_factor": 1.0,
        "cell_methods": "area: sum where grounded_ice_sheet (mask=sfgrlf) time: mean",
        "cmip7_compliant": True,
        "ismip7_type":"ST",
    },

    # iareafl — floating_ice_shelf_area — CMIP7: iareafl_tavg-u-hm-fis
    ("floating", "area"): {
        "cmip7_name": "iareafl",
        "standard_name": "floating_ice_shelf_area",
        "long_name": "Area Covered by Floating Ice Shelves",
        "cmip7_units": "m^2",
        "conversion_factor": 1.0,
        "cell_methods": "area: sum where floating_ice_shelf (mask=sftflf) time: mean",
        "cmip7_compliant": True,
        "ismip7_type":"ST",
    },

    # tendacabf — CMIP7: tendacabf_tavg-u-hm-is
    ("ice", "SMB"): {
        "cmip7_name": "tendacabf",
        "standard_name": "tendency_of_land_ice_mass_due_to_surface_mass_balance",
        "long_name": "Total Surface Mass Balance Flux",
        "cmip7_units": "kg s-1",
        "conversion_factor": ICE_DENSITY / SECS_PER_YEAR,
        "cell_methods": "area: sum where ice_sheet time: mean",
        "cmip7_compliant": True,
        "ismip7_type":"FL",
    },

    # tendlibmassbf — CMIP7: tendlibmassbf_tavg-u-hm-is
    ("ice", "BMB"): {
        "cmip7_name": "tendlibmassbf",
        "standard_name": "tendency_of_land_ice_mass_due_to_basal_mass_balance",
        "long_name": "Total Basal Mass Balance Flux",
        "cmip7_units": "kg s-1",
        "conversion_factor": ICE_DENSITY / SECS_PER_YEAR,
        "cell_methods": "area: sum where ice_sheet time: mean",
        "cmip7_compliant": True,
        "ismip7_type":"FL",
    },
   
    # tendlibmassbffl — CMIP7: tendlibmassbf_tavg-u-hm-is
    ("floating", "BMB"): {
        "cmip7_name": "tendlibmassbffl",
        "standard_name": "tendency_of_land_ice_mass_due_to_basal_mass_balance",
        "long_name": "Total Basal Mass Balance Flux",
        "cmip7_units": "kg s-1",
        "conversion_factor": ICE_DENSITY / SECS_PER_YEAR,
        "cell_methods": "area: sum where ice_sheet time: mean",
        "cmip7_compliant": True,
        "ismip7_type":"FL",
    },
    
    # tendlibmassbffl — CMIP7: tendlibmassbf_tavg-u-hm-is
    ("grounded", "BMB"): {
    "cmip7_name": "tendlibmassbfgr",
    "standard_name": "tendency_of_land_ice_mass_due_to_basal_mass_balance",
    "long_name": "Total Basal Mass Balance Flux",
    "cmip7_units": "kg s-1",
    "conversion_factor": ICE_DENSITY / SECS_PER_YEAR,
    "cell_methods": "area: sum where ice_sheet time: mean",
    "cmip7_compliant": True,
    "ismip7_type":"FL",
},

    
    
    # tendlicalvf — CMIP7: tendlicalvf_tavg-u-hm-is
    ("ice", "calving"): {
        "cmip7_name": "tendlicalvf",
        "standard_name": "tendency_of_land_ice_mass_due_to_calving",
        "long_name": "Total Calving Flux",
        "cmip7_units": "kg s-1",
        "conversion_factor": ICE_DENSITY / SECS_PER_YEAR,
        "cell_methods": "area: sum where ice_sheet time: mean",
        "cmip7_compliant": True,
        "comment": "Positive indicates mass loss by calving.",
        "ismip7_type":"FL",
    },

    # ------------------------------------------------------------------
    # Non-CMIP7 scalar diagnostics
    # Retained because the diagnostics tool computes them and they are
    # scientifically useful, but they have no out_name in the CMIP7 table.
    # ------------------------------------------------------------------

    # Total ice-fraction-weighted area — no CMIP7 scalar equivalent.
    ("ice", "fracArea"): {
        "cmip7_name": "iarea",
        "standard_name": "land_ice_area",
        "long_name": "Total Land Ice Area",
        "cmip7_units": "m^2",
        "conversion_factor": 1.0,
        "cell_methods": "area: sum where ice_sheet time: mean",
        "cmip7_compliant": False, "ismip7_type":"ST",
        "comment": (
            "Not a CMIP7 landIce out_name. "
            "Sum of ice-fraction-weighted cell areas over all ice-covered cells."
        ),
    },
}
