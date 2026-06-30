"""
Wrapper for the BISICLES flatten file tool and post-processing utilities for
writing CMIP7/CF-compliant 2D spatial NetCDF files from BISICLES plot HDF5 files.

Workflow
--------
1. Call the flatten executable to collapse the multi-level AMR hierarchy onto a
   single uniform grid and write a preliminary NetCDF file.
2. Open that NetCDF file in Python, apply unit conversions, rename variables to
   CMIP7 names, compute derived fields (sftgrf, sftflf from mask + iceFrac),
   and write a new, fully CF-compliant NetCDF file with all required metadata.

The flatten tool usage::

    flatten2d.*.ex <input.2d.hdf5> <output.nc> <level> [x0 [y0]]

where ``level`` is the refinement level to flatten onto (0 = coarsest grid,
1 = first refined level, etc.).

Grounded/floating mask fallback
---------------------------------
The derived fields ``sftgrf`` and ``sftflf`` require knowledge of which cells
contain grounded versus floating ice.  BISICLES plot files may or may not
include an explicit ``mask`` variable or ``iceFrac`` field.  When these are
absent, the module falls back to computing them from ice geometry:

* The **flotation criterion** determines whether ice at a given cell is
  grounded or floating based solely on thickness (*h*) and bed topography
  (*b*)::

      flotation_thickness = max(0, -b * rho_w / rho_i)
      floating  if  h > h_min  AND  b < 0  AND  h <= flotation_thickness
      grounded  if  h > h_min  AND  NOT floating

* When ``iceFrac`` is absent a binary approximation is used:
  ``iceFrac = 1`` where ``h > h_min``, else ``0``.  This overestimates ice
  extent at sub-grid margins; the ``bisicles_name`` attribute on the written
  variable records how the field was derived.

The function :func:`compute_flotation_mask` is also available as a public API
for use in other contexts.

Output file organisation
------------------------
Both :func:`process_plotfile` (single file) and :func:`process_directory`
(directory of files) write **one output NetCDF per CMIP7 variable**, named
``{cmip7_name}.nc``, inside the specified ``output_dir``.  When
processing a directory the time axis of each file spans all input timesteps,
giving a multi-year timeseries per variable.
"""

import os
import subprocess
import tempfile
import re
from pathlib import Path

import numpy as np

from .cmip7_vars import (
    FIELD_MAPPING,
    CF_FIELD_MAPPING,
    DERIVED_FIELDS,
    GRID_GEOMETRY_FIELDS,
    GROUNDED_MASK_VAL,
    FLOATING_MASK_VAL,
    OPEN_SEA_MASK_VAL,
    OPEN_LAND_MASK_VAL,
)
from .cf_utils import (
    CF_CONVENTIONS,
    FILL_VALUE,
    ISMIP7_FILL_VALUE,
    UKESM_GRID_ORIGINS,
    get_global_attributes,
    period_to_cmip_frequency,
    exact_days_since,
    add_time_variable,
    add_time_bounds,
    add_xy_variables,
    add_crs_variable,
    compute_latlon_arrays,
    add_latlon_variables,
    _ismip7_drs_filename,
)
from .filename_parser import parse_bisicles_filename

# Standalone BISICLES plot file patterns (no UKESM suite ID in the name)
_RE_STANDALONE_CFMEAN   = re.compile(r'plot\.CF\.\d+\.2d\.hdf5$',  re.IGNORECASE)
_RE_STANDALONE_SNAPSHOT = re.compile(r'plot\.\d+\.2d\.hdf5$',       re.IGNORECASE)


# ---------------------------------------------------------------------------
# Running the flatten tool
# ---------------------------------------------------------------------------

def run_flatten(
    input_hdf5,
    output_nc,
    exe_path,
    level=0,
    x0=None,
    y0=None,
):
    """
    Run the BISICLES flatten tool to project an AMR hierarchy onto a uniform grid.

    Parameters
    ----------
    input_hdf5 : str or Path
        Path to the BISICLES plot HDF5 file.
    output_nc : str or Path
        Path for the output NetCDF file (must end in .nc).
    exe_path : str or Path
        Full path to the flatten executable.
    level : int
        AMR level to flatten onto.  0 = coarsest grid; higher values give finer
        resolution.
    x0 : float, optional
        X-coordinate origin override (metres).
    y0 : float, optional
        Y-coordinate origin override (metres).

    Raises
    ------
    RuntimeError
        If the flatten tool exits with a non-zero return code.
    """
    if level < 0:
        raise ValueError(
            f"level must be >= 0 (0 = coarsest grid). Got: {level}"
        )
    if exe_path is None:
        raise ValueError(
            "exe_path is required. Provide the full path to the flatten executable, e.g.:\n"
            "  --exe-path /path/to/filetools/flatten2d.Linux.64.g++.gfortran.OPT.ex"
        )
    p = Path(exe_path)
    if p.is_dir():
        raise ValueError(
            f"exe_path points to a directory, not an executable: {exe_path}\n"
            "Provide the full path to the executable file, e.g.:\n"
            "  /path/to/filetools/flatten2d.Linux.64.g++.gfortran.OPT.ex"
        )
    if not p.is_file():
        raise FileNotFoundError(f"Flatten executable not found: {exe_path}")
    exe = str(p)
    args = [exe, str(input_hdf5), str(output_nc), str(level)]
    if x0 is not None and y0 is not None:
        args += [str(x0), str(y0)]

    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"flatten tool failed (return code {result.returncode}):\n"
            f"stdout: {result.stdout}\n"
            f"stderr: {result.stderr}"
        )


# ---------------------------------------------------------------------------
# Flotation-based mask computation
# ---------------------------------------------------------------------------

def compute_flotation_mask(
    thickness,
    topography,
    ice_density=918.0,
    water_density=1028.0,
    h_min=1.0,
):
    """
    Compute an integer ice-type mask from ice geometry using the flotation criterion.

    A cell is classified as floating when the ice is thin enough that the
    ocean can support it hydrostatically:

    .. code-block:: text

        flotation_thickness = max(0, -topography * rho_w / rho_i)
        floating  if  thickness > h_min  AND  topography < 0
                                         AND  thickness <= flotation_thickness
        grounded  if  thickness > h_min  AND  NOT floating

    This replicates the logic in BISICLES ``SigmaCSF.ChF`` / ``LevelSigmaCS``.

    Parameters
    ----------
    thickness : ndarray
        Ice thickness in metres.
    topography : ndarray
        Bed topography (positive up) in metres.  Negative values are below
        present-day sea level.
    ice_density : float
        Ice density in kg m-3 (default 918.0).
    water_density : float
        Ocean water density in kg m-3 (default 1028.0).
    h_min : float
        Minimum ice thickness threshold in metres below which a cell is
        treated as ice-free (default 1.0 m).

    Returns
    -------
    mask : ndarray of int, same shape as *thickness*
        Integer mask with values:

        * ``GROUNDED_MASK_VAL`` (1) – grounded ice
        * ``FLOATING_MASK_VAL`` (2) – floating ice
        * ``OPEN_SEA_MASK_VAL`` (4) – ice-free ocean (bed < 0)
        * ``OPEN_LAND_MASK_VAL`` (8) – ice-free land (bed >= 0)
    """
    # Thickness of ice needed just to ground: if bed is above sea level,
    # ice is always grounded regardless of thickness.
    flotation_thickness = np.maximum(
        0.0, -topography * (water_density / ice_density)
    )

    has_ice = thickness > h_min
    below_sea = topography < 0.0

    is_floating = has_ice & below_sea & (thickness <= flotation_thickness)
    is_grounded = has_ice & ~is_floating
    is_open_sea = ~has_ice & below_sea

    # Start with everything as open land, then overwrite
    mask = np.full(thickness.shape, OPEN_LAND_MASK_VAL, dtype=np.int32)
    mask = np.where(is_open_sea, OPEN_SEA_MASK_VAL, mask)
    mask = np.where(is_floating, FLOATING_MASK_VAL, mask)
    mask = np.where(is_grounded, GROUNDED_MASK_VAL, mask)
    return mask



def _regrid(flat_data, xygrid):
    """
    Interpolate flat_data fields from the grid defined by flat_data['x'],
    flat_data['y'] to the grid defined by x, y = xygrid.

    Points on the target grid that fall outside the source domain are filled
    with NaN; the caller is responsible for replacing these with a fill value.
    """

    from scipy.interpolate import RegularGridInterpolator

    result = {}
    result["x"], result["y"] = xygrid
    result["time"] = flat_data["time"]
    result["time_units"] = flat_data["time_units"]
    result["epsg"] = flat_data["epsg"]
    result["variables"] = {}
    xy = np.meshgrid(result["y"], result["x"])

    for key, var in flat_data["variables"].items():
        rg = RegularGridInterpolator(
            (flat_data["y"], flat_data["x"]),
            var,
            bounds_error=False,
            fill_value=np.nan,
        )
        result["variables"][key] = rg(xy).T

    return result
    

# ---------------------------------------------------------------------------
# Reading the flatten tool output and applying CMIP7 conversions
# ---------------------------------------------------------------------------

def _read_flatten_nc(nc_path):
    """
    Read a NetCDF file produced by the flatten tool and return its contents.

    Returns
    -------
    dict with keys:
      'x', 'y'       : 1-D coordinate arrays (m)
      'time'         : scalar simulation time (years)
      'epsg'         : EPSG code (int) or None
      'variables'    : dict of {bisicles_name: ndarray (ny, nx)}
      'raw_attrs'    : dict of global attributes from the flatten NC
    """
    from netCDF4 import Dataset

    result = {"variables": {}, "x": None, "y": None, "time": None, "epsg": None}

    with Dataset(str(nc_path), "r") as ds:
        result["raw_attrs"] = {k: ds.getncattr(k) for k in ds.ncattrs()}

        # Coordinate variables
        if "x" in ds.variables:
            result["x"] = ds.variables["x"][:].data
        if "y" in ds.variables:
            result["y"] = ds.variables["y"][:].data
        if "time" in ds.variables:
            tv = ds.variables["time"]
            t_data = tv[:]
            # flatten tool may write a 1-D array; take the last (or only) value
            result["time"] = float(np.asarray(t_data).flat[-1])
            # Try to extract simulation time in years from units attribute
            if hasattr(tv, "units"):
                result["time_units"] = tv.units
        if "crs" in ds.variables:
            crs_v = ds.variables["crs"]
            if hasattr(crs_v, "epsg_code"):
                try:
                    result["epsg"] = int(str(crs_v.epsg_code).replace("EPSG:", ""))
                except ValueError:
                    pass

        # Time bounds — read if the flatten tool wrote them (CF output files)
        for _bnds_name in ("time_bnds", "time_bounds"):
            if _bnds_name in ds.variables:
                _bnds = np.asarray(ds.variables[_bnds_name][:]).ravel()
                if len(_bnds) >= 2:
                    result["time_start"] = float(_bnds[0])
                    result["time_end"]   = float(_bnds[-1])
                break

        # Data variables (skip coordinate/metadata variables)
        skip = {"x", "y", "time", "time_bnds", "time_bounds", "crs", "level"}
        for name in ds.variables:
            if name in skip:
                continue
            v = ds.variables[name]
            if v.ndim < 2:
                continue
            arr = v[:]
            if hasattr(arr, "data"):
                arr = arr.data.astype(float)
                if hasattr(v, "_FillValue"):
                    arr[arr == float(v._FillValue)] = np.nan
            result["variables"][name] = arr

    return result


def _compute_derived_fields(
    variables,
    ice_density=918.0,
    water_density=1028.0,
    h_min=1.0,
):
    """
    Compute derived CMIP7 area-fraction fields (sftgrf, sftflf) in-place.

    Prioritises model-output fields where available, falling back to
    geometry-based computation when they are absent:

    * **mask** – used directly if present; otherwise computed from
      ``thickness`` and ``Z_base`` via the flotation criterion
      (:func:`compute_flotation_mask`).
    * **iceFrac** – used directly if present; otherwise approximated as a
      binary field (1 where ``thickness > h_min``, else 0).

    The fallback paths are recorded on the returned ``_derived_sources`` dict
    which ``write_cmip7_per_variable_netcdfs`` uses to annotate the output
    variables.

    Parameters
    ----------
    variables : dict
        Mutable dict of {bisicles_name: ndarray}.  Modified in-place.
    ice_density : float
        Ice density in kg m-3, used only when computing the flotation mask.
    water_density : float
        Ocean water density in kg m-3, used only when computing the flotation mask.
    h_min : float
        Ice thickness threshold in metres, used when deriving iceFrac and
        when computing the flotation mask.

    Returns
    -------
    sources : dict
        Maps derived variable name -> provenance string describing how the
        field was computed.  Empty if no derived fields could be produced.
    """
    sources = {}

    # ------------------------------------------------------------------
    # Resolve ice fraction
    # ------------------------------------------------------------------
    if "iceFrac" in variables:
        ice_frac = variables["iceFrac"]
        icefrac_source = "iceFrac from plot file"
    elif "thickness" in variables:
        ice_frac = np.where(variables["thickness"] > h_min, 1.0, 0.0)
        icefrac_source = (
            f"binary approximation (1 where thickness > {h_min} m, else 0); "
            "sub-grid margin fractions are not represented"
        )
    else:
        return sources  # Cannot proceed without at least thickness

    # ------------------------------------------------------------------
    # Resolve grounded/floating mask
    # ------------------------------------------------------------------
    if "mask" in variables:
        mask = variables["mask"]
        mask_source = "mask from plot file"
    elif "thickness" in variables and "Z_base" in variables:
        mask = compute_flotation_mask(
            variables["thickness"],
            variables["Z_base"],
            ice_density=ice_density,
            water_density=water_density,
            h_min=h_min,
        )
        mask_source = (
            f"computed from flotation criterion using thickness and Z_base "
            f"(rho_ice={ice_density} kg/m3, rho_water={water_density} kg/m3, "
            f"h_min={h_min} m)"
        )
    else:
        return sources  # Cannot determine grounded/floating without geometry

    # ------------------------------------------------------------------
    # Compute sftgrf and sftflf
    # ------------------------------------------------------------------
    ice_frac_clean = np.where(np.isnan(ice_frac), 0.0, ice_frac)
    mask_int = np.round(mask).astype(np.int32)

    variables["sftgrf"] = np.where(mask_int == GROUNDED_MASK_VAL, ice_frac_clean, 0.0)
    variables["sftflf"] = np.where(mask_int == FLOATING_MASK_VAL, ice_frac_clean, 0.0)

    provenance = f"ice_frac: {icefrac_source}; mask: {mask_source}"
    sources["sftgrf"] = provenance
    sources["sftflf"] = provenance
    return sources


# ---------------------------------------------------------------------------
# Module-level variable-writing helper
# ---------------------------------------------------------------------------

def _write_2d_var(ds, out_name, arr_3d, mapping, time_cell_method, grid_mapping,
                  bisicles_name=None, coordinates="", fill_value=FILL_VALUE):
    """
    Write one ``(time, y, x)`` data variable with full CF/CMIP7 metadata.

    Parameters
    ----------
    ds : netCDF4.Dataset
        Open, writable Dataset; must already have ``time``, ``y``, ``x``
        dimensions defined.
    out_name : str
        Variable name in the output file (CMIP7 name).
    arr_3d : ndarray, shape (T, ny, nx)
        Data array in BISICLES units; the mapping's ``conversion_factor`` is
        applied inside this function.
    mapping : dict
        Entry from :data:`FIELD_MAPPING` or :data:`CF_FIELD_MAPPING`.
    time_cell_method : str
        ``"time: mean"`` or ``"time: point"``.  Controls whether the stored
        ``cell_methods`` attribute reflects averaging or instantaneous output.
    grid_mapping : str or None
        Name of the CRS variable, or ``None`` if no projection is defined.
    bisicles_name : str, optional
        BISICLES internal field name; written as ``bisicles_name`` attribute
        when supplied.
    coordinates : str
        Space-separated list of *auxiliary* coordinate variable names to
        record in the CF ``coordinates`` attribute.  Per CF-1.12, only
        auxiliary coordinates (e.g. ``"lat lon"``) should be listed here —
        dimension coordinates (``x``, ``y``, ``time``) are self-identifying
        and must NOT appear.  Pass an empty string to omit the attribute.
    fill_value : float
        Fill value to use for missing data.  Use :data:`FILL_VALUE` (1e20,
        default) for CMIP7/UKESM output and :data:`ISMIP7_FILL_VALUE` for
        ISMIP7 submissions.
    """
    arr = arr_3d * mapping["conversion_factor"]
    arr = np.where(np.isnan(arr), fill_value, arr)
    cell_methods = mapping["cell_methods"]
    if time_cell_method == "time: point":
        cell_methods = cell_methods.replace("time: mean", "time: point")
    var = ds.createVariable(out_name, "f4", ("time", "y", "x"), fill_value=fill_value)
    var[:] = arr
    var.standard_name = mapping["standard_name"]
    var.long_name = mapping["long_name"]
    var.units = mapping["cmip7_units"]
    var.missing_value = np.float32(fill_value)
    var.cell_methods = cell_methods
    var.modeling_realm = mapping["modeling_realm"]
    if coordinates:
        var.coordinates = coordinates
    if grid_mapping:
        var.grid_mapping = grid_mapping
    if "comment" in mapping:
        var.comment = mapping["comment"]
    if bisicles_name is not None:
        var.bisicles_name = bisicles_name
    var.bisicles_units = mapping["bisicles_units"]


# ---------------------------------------------------------------------------
# Writing per-variable CMIP7 NetCDF files
# ---------------------------------------------------------------------------

def write_cmip7_per_variable_netcdfs(
    all_data,
    output_dir,
    epsg_code=None,
    x0=None,
    y0=None,
    reference_year=1850,
    calendar="standard",
    cmip7_only=False,
    ice_density=918.0,
    water_density=1028.0,
    h_min=1.0,
    institution="",
    source="BISICLES adaptive mesh refinement ice sheet model",
    experiment="",
    variant_label="",
    ice_sheet="",
    extra_attrs=None,
    source_files=None,
    frequency="",
    ismip7_mode=False,
    model_id="",
    member_id="",
    # ISMIP7 DRS filename components
    source_id="",
    ism_id="",
    ism_member_id="m001",
    esm_id="standalone",
    forcing_member_id="f001",
    set_counter="C001",
    # ISMIP7 mandatory global attributes
    group="",
    contact_name="",
    contact_email="",
    grid_spec=None,
):
    """
    Write one CMIP7/CF-compliant NetCDF per variable from a collection of
    flatten outputs.

    Each output file is named ``{cmip7_name}.nc`` inside *output_dir*
    and contains a single data variable with a ``time`` dimension spanning all
    supplied timesteps.

    Parameters
    ----------
    all_data : list of (dict, BISICLESFileInfo or None)
        Each element is ``(flatten_data, file_info)`` where *flatten_data* is
        the dict returned by :func:`_read_flatten_nc` and *file_info* is the
        result of :func:`~.filename_parser.parse_bisicles_filename` (or
        ``None`` when filename parsing is not applicable).  One entry per
        input timestep.
    output_dir : str or Path
        Directory for output files.  Created if it does not exist.
    epsg_code : int, optional
        EPSG code for the projection.  Overrides the value read from the
        flatten NetCDF files when supplied.
    x0 : float, optional
        X-coordinate of the lower-left corner of the domain (metres).
    y0 : float, optional
        Y-coordinate of the lower-left corner of the domain (metres).
    reference_year : int
        Reference year for the CF time axis (default 1850).
    calendar : str
        CF calendar name: ``"gregorian"`` (default) or ``"360_day"``.
    cmip7_only : bool
        If True, only write CMIP7-standard variables.  If False (default),
        also write unmapped BISICLES variables with their original names.
    ice_density : float
        Ice density in kg m-3 for the flotation-based mask fallback.
    water_density : float
        Ocean water density in kg m-3 for the flotation-based mask fallback.
    h_min : float
        Minimum ice thickness in metres for the mask/iceFrac fallback.
    institution : str
        Written to the global ``institution`` attribute.
    source : str
        Written to the global ``source`` attribute.
    experiment : str
        Experiment identifier.
    variant_label : str
        CMIP variant label.
    ice_sheet : str
        Ice sheet identifier (e.g. ``'GrIS'``, ``'AIS'``).
    extra_attrs : dict, optional
        Additional global attributes.
    ismip7_mode : bool
        If True, write ISMIP7-compliant output: DRS filenames, ``Conventions``
        set to ``"CF-1.12"``, and ISMIP7 mandatory global attributes.
    model_id : str
        Model identifier for the CMIP7-coupled workflow (global attribute).
    member_id : str
        Member identifier for the CMIP7-coupled workflow (global attribute).
    source_id : str
        ISMIP7 DRS: modelling group name (e.g. ``"BristolGlaciology"``).
    ism_id : str
        ISMIP7 DRS: ISM name and version (e.g. ``"BISICLES"``).
    ism_member_id : str
        ISMIP7 DRS: ISM choice variant identifier (default ``"m001"``).
    esm_id : str
        ISMIP7 DRS: CMIP ESM used for forcing (e.g. ``"CESM2-WACCM"`` or
        ``"standalone"`` for idealised forcing).
    forcing_member_id : str
        ISMIP7 DRS: forcing choice variant identifier (default ``"f001"``).
    set_counter : str
        ISMIP7 DRS: set counter (e.g. ``"C001"``).
    group : str
        ISMIP7 mandatory global attribute: modelling group name.
    contact_name : str
        ISMIP7 mandatory global attribute: contact person name(s).
    contact_email : str
        ISMIP7 mandatory global attribute: contact person email(s).

    Returns
    -------
    list of Path
        Paths of output NetCDF files written (one per variable).
    """
    from netCDF4 import Dataset

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not all_data:
        raise ValueError("all_data is empty; nothing to write.")

    # -----------------------------------------------------------------------
    # Phase 1: accumulate time metadata and per-variable 2-D arrays
    # -----------------------------------------------------------------------
    times_list = []
    is_time_mean_list = []
    time_start_list = []
    time_end_list = []
    time_cell_method_list = []

    # Accumulators: {name: [arr_at_t0, arr_at_t1, ...]}
    bisicles_arrays = {}
    cf_arrays = {}
    derived_arrays = {}
    derived_src = {}     # {derived_name: provenance_str} – from last timestep
    unmapped_arrays = {}

    x = None
    y = None
    epsg = epsg_code  # overridden from first file when not supplied

    for flatten_data, file_info in all_data:
        # Coordinates and projection from first file (assumed constant)
        if x is None:
            x = flatten_data.get("x")
            y = flatten_data.get("y")
            if epsg is None:
                epsg = flatten_data.get("epsg")

        # Time metadata
        if file_info is not None:
            t = file_info.time_years
            is_mean = file_info.is_time_mean
            t_start = file_info.start_time_years
            t_end = file_info.end_time_years
            tcm = file_info.cell_methods_time
        else:
            t_raw = flatten_data.get("time") or 0.0
            is_cf_mean = flatten_data.get("_is_cf_mean", False)
            is_mean = is_cf_mean
            if is_cf_mean:
                # Use time bounds from the flatten NC if the flatten tool wrote them,
                # otherwise infer the averaging period from the time value:
                #   N.5  → mid-year convention: averaging year N, period [N, N+1]
                #   N.0  → end-of-period convention: averaging year N-1, period [N-1, N]
                if "time_start" in flatten_data and "time_end" in flatten_data:
                    t_start = flatten_data["time_start"]
                    t_end   = flatten_data["time_end"]
                else:
                    frac = t_raw - int(t_raw)
                    if abs(frac - 0.5) < 0.1:
                        year = int(t_raw)           # midpoint: 2015.5 → year 2015
                    else:
                        year = int(round(t_raw)) - 1  # end-of-period: 2016.0 → year 2015
                    t_start = float(year)
                    t_end   = float(year + 1)
                t = 0.5 * (t_start + t_end)
                tcm = "time: mean"
            else:
                t = t_raw
                t_start = None
                t_end = None
                tcm = "time: point"

        times_list.append(t)
        is_time_mean_list.append(is_mean)
        time_start_list.append(t_start)
        time_end_list.append(t_end)
        time_cell_method_list.append(tcm)

        # Work on a mutable copy so _compute_derived_fields can add sftgrf/sftflf
        variables = dict(flatten_data.get("variables", {}))
        sources = _compute_derived_fields(
            variables,
            ice_density=ice_density,
            water_density=water_density,
            h_min=h_min,
        )
        derived_src.update(sources)

        # BISICLES internal name -> CMIP7
        for bname in FIELD_MAPPING:
            if bname in variables:
                bisicles_arrays.setdefault(bname, []).append(variables[bname])

        # CF-output name -> CMIP7 (plot.CF-*.hdf5 files)
        written_cf = set()
        for cf_name in CF_FIELD_MAPPING:
            if cf_name in variables:
                cf_arrays.setdefault(cf_name, []).append(variables[cf_name])
                written_cf.add(cf_name)

        # Derived fields (sftgrf, sftflf) – skip if CF file already provides them
        for dname in DERIVED_FIELDS:
            if dname not in written_cf and dname in variables:
                derived_arrays.setdefault(dname, []).append(variables[dname])

        # Unmapped BISICLES variables
        if not cmip7_only:
            known = set(FIELD_MAPPING) | set(CF_FIELD_MAPPING) | set(DERIVED_FIELDS)
            for bname, arr in variables.items():
                if bname not in known and arr.ndim >= 2:
                    unmapped_arrays.setdefault(bname, []).append(arr)

    if x is None or y is None:
        raise ValueError("Flatten data is missing x or y coordinate arrays.")

    # -----------------------------------------------------------------------
    # Phase 2: sort all per-timestep data by ascending time
    # -----------------------------------------------------------------------
    sort_idx = sorted(range(len(times_list)), key=lambda i: times_list[i])
    times_sorted = [times_list[i] for i in sort_idx]
    is_mean_sorted = [is_time_mean_list[i] for i in sort_idx]
    t_start_sorted = [time_start_list[i] for i in sort_idx]
    t_end_sorted = [time_end_list[i] for i in sort_idx]
    tcm_sorted = [time_cell_method_list[i] for i in sort_idx]

    # Use the cell_method from the first timestep for the whole file.
    # Within a single run all files are consistently either time-mean or snapshot.
    time_cell_method = tcm_sorted[0] if tcm_sorted else "time: point"
    has_time_bounds = any(is_mean_sorted)
    time_arr = np.asarray(times_sorted)

    # -----------------------------------------------------------------------
    # Compute ISMIP7-correct time coordinates using exact calendar arithmetic:
    #   Annual means  → timestamp = YYYY-07-01, bounds = [YYYY-01-01, (YYYY+1)-01-01]
    #   Snapshots     → timestamp = YYYY-01-01 (integer year rounded)
    # This avoids the ±1-3 day error in the 365.25-days/year approximation.
    # -----------------------------------------------------------------------
    _cal_key = "standard" if calendar == "gregorian" else calendar

    def _to_exact_days(t, is_mean, t_start, t_end):
        """Return (time_days, bound_start_days, bound_end_days) for one timestep."""
        if is_mean and t_start is not None and t_end is not None:
            s = round(t_start)
            e = round(t_end)
            if abs(t_start - s) < 0.01 and abs(t_end - e) < 0.01 and e - s == 1:
                # Annual mean with integer Jan-1 boundaries → exact July-1 timestamp
                return (
                    exact_days_since(int(s), 7, 1, reference_year, _cal_key),
                    exact_days_since(int(s), 1, 1, reference_year, _cal_key),
                    exact_days_since(int(e), 1, 1, reference_year, _cal_key),
                )
        # Snapshot or non-annual mean → round to nearest integer year, use Jan-1
        yr = int(round(t))
        d = exact_days_since(yr, 1, 1, reference_year, _cal_key)
        return d, d, d  # bounds ignored for snapshots

    _exact_t  = []
    _exact_bs = []
    _exact_be = []
    for _i in range(len(times_sorted)):
        _td, _bs, _be = _to_exact_days(
            times_sorted[_i], is_mean_sorted[_i],
            t_start_sorted[_i], t_end_sorted[_i],
        )
        _exact_t.append(_td)
        _exact_bs.append(_bs)
        _exact_be.append(_be)
    _exact_t  = np.array(_exact_t)
    _exact_bs = np.array(_exact_bs) if has_time_bounds else None
    _exact_be = np.array(_exact_be) if has_time_bounds else None

    def _sort(arr_list):
        """Return arr_list reordered by sort_idx."""
        return [arr_list[i] for i in sort_idx]

    # -----------------------------------------------------------------------
    # Phase 3: shared file-setup and write helpers
    # -----------------------------------------------------------------------
    ny, nx_size = len(y), len(x)
    crs_name = "crs"
    grid_mapping = crs_name if epsg is not None else None

    # Attempt to compute 2-D lat/lon auxiliary coordinates (requires pyproj).
    # lat/lon are required by CF-1.12 / CMOR for projected-coordinate grids.
    lat_2d = lon_2d = None
    if epsg is not None and epsg != 4326:
        try:
            lat_2d, lon_2d = compute_latlon_arrays(x, y, epsg)
        except ImportError as _exc:
            import warnings
            warnings.warn(
                f"pyproj is not installed; lat/lon auxiliary coordinates will "
                f"not be written (CF/CMOR compliance requires them). "
                f"Install pyproj with: pip install pyproj\n"
                f"Original error: {_exc}",
                stacklevel=2,
            )
    has_latlon = lat_2d is not None
    coords_str = "lat lon" if has_latlon else ""

    if not frequency:
        first_info = next((fi for _, fi in all_data if fi is not None), None)
        if first_info is not None and first_info.period:
            frequency = period_to_cmip_frequency(first_info.period)

    # ISMIP7 requires f4 precision and the netCDF4 default fill value;
    # UKESM/CMIP7 output uses f8 for time/geometry and the conventional 1e20 fill value.
    _fill_value = ISMIP7_FILL_VALUE if ismip7_mode else FILL_VALUE
    _time_dtype = "f4" if ismip7_mode else "f8"
    _geom_dtype = "f4" if ismip7_mode else "f8"

    # In ISMIP7 mode, populate group/model from DRS fields when not set explicitly
    _group = group or source_id
    _crs_str = f"epsg:{epsg}" if (ismip7_mode and epsg is not None) else ""

    # Derive human-readable grid description and nominal resolution from x/y spacing
    _dx = abs(float(x[1] - x[0])) if len(x) > 1 else 0.0
    _dy = abs(float(y[1] - y[0])) if len(y) > 1 else 0.0
    _res_m = max(_dx, _dy)
    if _res_m >= 1000.0:
        _res_str = f"{_res_m / 1000:.4g} km"
    else:
        _res_str = f"{_res_m:.4g} m"

    _epsg_names = {
        3031: "Antarctic Polar Stereographic",
        3413: "NSIDC Sea Ice Polar Stereographic North",
    }
    _proj_name = _epsg_names.get(epsg, f"EPSG:{epsg}") if epsg is not None else "unknown projection"
    _grid_desc = f"{_proj_name} (EPSG:{epsg}), {_res_str}" if epsg is not None else ""

    # Build regrid provenance metadata when a target grid was supplied
    _extra_history = ""
    _grid_label = "gn"
    if grid_spec is not None:
        x_min, x_max, nx, y_min, y_max, ny = grid_spec
        _extra_history = (
            f"Regridded to x=[{x_min},{x_max}] nx={nx}, "
            f"y=[{y_min},{y_max}] ny={ny} "
            f"using scipy.interpolate.RegularGridInterpolator (linear)"
        )
        if not ismip7_mode:
            _grid_label = "gr"

    global_attrs = get_global_attributes(
        institution=institution,
        source=source,
        experiment=experiment,
        variant_label=variant_label,
        ice_sheet=ice_sheet,
        source_files=source_files,
        frequency=frequency,
        conventions="CF-1.12" if ismip7_mode else CF_CONVENTIONS,
        model_id=model_id,
        member_id=member_id,
        group=_group if ismip7_mode else group,
        contact_name=contact_name,
        contact_email=contact_email,
        crs=_crs_str,
        extra_history=_extra_history,
        grid_label=_grid_label,
        grid=_grid_desc,
        nominal_resolution=_res_str,
    )
    if ismip7_mode and ism_id:
        global_attrs["model"] = ism_id
    global_attrs["external_variables"] = "modelcellareai"
    if extra_attrs:
        global_attrs.update(extra_attrs)

    # Pre-build time-bounds lists once (used inside _setup_ds)
    if has_time_bounds:
        _tb_start = [
            t_start_sorted[i] if is_mean_sorted[i] else times_sorted[i]
            for i in range(len(times_sorted))
        ]
        _tb_end = [
            t_end_sorted[i] if is_mean_sorted[i] else times_sorted[i]
            for i in range(len(times_sorted))
        ]
    else:
        _tb_start = _tb_end = None

    def _setup_ds(ds):
        """Populate dimensions, coordinates, and global attributes."""
        ds.setncatts(global_attrs)
        ds.createDimension("time", None)  # unlimited record dimension
        ds.createDimension("y", ny)
        ds.createDimension("x", nx_size)
        add_xy_variables(ds, x, y, epsg_code=epsg)
        add_time_variable(ds, time_arr, reference_year=reference_year,
                          calendar=calendar, dtype=_time_dtype,
                          time_days=_exact_t)
        if has_time_bounds:
            add_time_bounds(ds, _tb_start, _tb_end,
                            reference_year=reference_year, calendar=calendar,
                            dtype=_time_dtype,
                            start_days=_exact_bs,
                            end_days=_exact_be)
        if epsg is not None:
            add_crs_variable(ds, epsg, x0=x0, y0=y0)
        if has_latlon:
            add_latlon_variables(ds, lat_2d, lon_2d)

    def _title(long_name):
        """Return the NetCDF title attribute for a variable."""
        if ismip7_mode:
            return long_name
        return f"UniCiCles (BISICLES) output from UKESM: {long_name}"

    def _out_path(varname, mask_no=0):
        """Return the output Path for a variable."""
        if ismip7_mode:
            fname = _ismip7_drs_filename(
                varname, ice_sheet, source_id, ism_id, ism_member_id,
                esm_id, forcing_member_id, experiment, set_counter,
                times_sorted, mask_no=mask_no,
            )
        else:
            fname = f"{varname}.nc"
        return output_dir / fname

    output_files = []

    # -----------------------------------------------------------------------
    # Phase 4: write one file per variable
    # -----------------------------------------------------------------------

    # 1. BISICLES internal name -> CMIP7 (standard plot files)
    for bisicles_name, mapping in FIELD_MAPPING.items():
        if bisicles_name not in bisicles_arrays:
            continue
        arr_stack = np.stack(_sort(bisicles_arrays[bisicles_name]), axis=0)
        out_name = mapping["cmip7_name"]
        out_path = _out_path(out_name)
        with Dataset(str(out_path), "w", format="NETCDF4") as ds:
            _setup_ds(ds)
            ds.variable_id = out_name
            ds.variable_name = out_name
            ds.title = _title(mapping["long_name"])
            _write_2d_var(ds, out_name, arr_stack, mapping, time_cell_method,
                          grid_mapping, bisicles_name=bisicles_name,
                          coordinates=coords_str, fill_value=_fill_value)
        output_files.append(out_path)

    # 2. CF-output name -> CMIP7 (plot.CF-*.hdf5 files)
    for cf_name, mapping in CF_FIELD_MAPPING.items():
        if cf_name not in cf_arrays:
            continue
        arr_stack = np.stack(_sort(cf_arrays[cf_name]), axis=0)
        out_name = mapping["cmip7_name"]
        out_path = _out_path(out_name)
        with Dataset(str(out_path), "w", format="NETCDF4") as ds:
            _setup_ds(ds)
            ds.variable_id = out_name
            ds.variable_name = out_name
            ds.title = _title(mapping["long_name"])
            _write_2d_var(ds, out_name, arr_stack, mapping, time_cell_method,
                          grid_mapping, coordinates=coords_str,
                          fill_value=_fill_value)
        output_files.append(out_path)

    # 3. Derived fields (sftgrf, sftflf) – skip if CF file already provided them
    for derived_name, dmeta in DERIVED_FIELDS.items():
        if derived_name in cf_arrays:
            continue
        if derived_name not in derived_arrays:
            continue
        arr_stack = np.stack(_sort(derived_arrays[derived_name]), axis=0)
        arr_stack = arr_stack * dmeta["conversion_factor"]
        arr_stack = np.where(np.isnan(arr_stack), _fill_value, arr_stack)
        cell_methods = dmeta["cell_methods"]
        if time_cell_method == "time: point":
            cell_methods = cell_methods.replace("time: mean", "time: point")
        out_path = _out_path(derived_name)
        with Dataset(str(out_path), "w", format="NETCDF4") as ds:
            _setup_ds(ds)
            ds.variable_id = derived_name
            ds.variable_name = derived_name
            ds.title = _title(dmeta["long_name"])
            var = ds.createVariable(
                derived_name, "f4", ("time", "y", "x"), fill_value=_fill_value
            )
            var[:] = arr_stack
            var.standard_name = dmeta["standard_name"]
            var.long_name = dmeta["long_name"]
            var.units = dmeta["cmip7_units"]
            var.missing_value = np.float32(_fill_value)
            var.cell_methods = cell_methods
            var.modeling_realm = dmeta["modeling_realm"]
            if coords_str:
                var.coordinates = coords_str
            if grid_mapping:
                var.grid_mapping = grid_mapping
            if "comment" in dmeta:
                var.comment = dmeta["comment"]
            if derived_name in derived_src:
                var.derived_from = derived_src[derived_name]
        output_files.append(out_path)

    # 4. Unmapped BISICLES variables (preserved unless cmip7_only=True)
    if not cmip7_only:
        for bname, arr_list in unmapped_arrays.items():
            arr_stack = np.stack(_sort(arr_list), axis=0)
            arr_stack = np.where(np.isnan(arr_stack), _fill_value, arr_stack)
            safe_name = bname.replace("/", "_")
            out_path = output_dir / f"{safe_name}.nc"
            with Dataset(str(out_path), "w", format="NETCDF4") as ds:
                _setup_ds(ds)
                ds.variable_id = safe_name
                ds.variable_name = safe_name
                ds.title = _title(bname)
                var = ds.createVariable(
                    safe_name, "f4", ("time", "y", "x"), fill_value=_fill_value
                )
                var[:] = arr_stack
                var.long_name = bname
                var.missing_value = np.float32(FILL_VALUE)
                if coords_str:
                    var.coordinates = coords_str
                var.note = "Unmapped BISICLES field; units as in original plot file."
                if grid_mapping:
                    var.grid_mapping = grid_mapping
            output_files.append(out_path)

    # 5. Grid-geometry fields (modelcellareai) — computed from x/y, no time dim
    dx = abs(float(x[1] - x[0])) if nx_size > 1 else 1.0
    dy = abs(float(y[1] - y[0])) if ny > 1 else 1.0
    cell_area = np.full((ny, nx_size), dx * dy, dtype=np.float64)

    for geom_name, gmeta in GRID_GEOMETRY_FIELDS.items():
        out_path = output_dir / f"{geom_name}.nc"
        with Dataset(str(out_path), "w", format="NETCDF4") as ds:
            ds.setncatts(global_attrs)
            ds.variable_id = geom_name
            ds.variable_name = geom_name
            ds.title = _title(gmeta["long_name"])
            ds.createDimension("y", ny)
            ds.createDimension("x", nx_size)
            add_xy_variables(ds, x, y, epsg_code=epsg)
            if epsg is not None:
                add_crs_variable(ds, epsg, x0=x0, y0=y0)
            if has_latlon:
                add_latlon_variables(ds, lat_2d, lon_2d)
            var = ds.createVariable(geom_name, _geom_dtype, ("y", "x"))
            var[:] = cell_area
            var.standard_name = gmeta["standard_name"]
            var.long_name = gmeta["long_name"]
            var.units = gmeta["cmip7_units"]
            var.cell_methods = gmeta["cell_methods"]
            var.modeling_realm = gmeta["modeling_realm"]
            if coords_str:
                var.coordinates = coords_str
            if grid_mapping:
                var.grid_mapping = grid_mapping
        output_files.append(out_path)

    return output_files


# ---------------------------------------------------------------------------
# Internal helper: run flatten tool and read result for one plot file
# ---------------------------------------------------------------------------

def _flatten_plot_file(
    plot_file,
    exe_path,
    level=0,
    x0=None,
    y0=None,
    xygrid=None,
    keep_intermediate=False,
    intermediate_nc=None,
    verbose=False,
):
    """
    Run the flatten executable on *plot_file* and return the parsed data.

    Parameters
    ----------
    plot_file : Path
        Input BISICLES plot HDF5 file.
    exe_path : str or Path
        Full path to the flatten executable.
    level : int
        AMR level to flatten onto.
    x0, y0 : float, optional
        Grid origin overrides.
    keep_intermediate : bool
        If True, keep the intermediate NetCDF file (written next to
        *intermediate_nc* or alongside *plot_file*).
    intermediate_nc : Path, optional
        Explicit path for the intermediate file.  A temporary file is used
        when not supplied.
    verbose : bool
        Print progress messages.

    Returns
    -------
    tuple (flatten_data, file_info)
        *flatten_data* is the dict from :func:`_read_flatten_nc`;
        *file_info* is from :func:`~.filename_parser.parse_bisicles_filename`
        (may be ``None`` for non-UKESM filenames).
    """
    if keep_intermediate and intermediate_nc is not None:
        tmp_nc = intermediate_nc
        _tmp_path = None
    else:
        fd, _tmp = tempfile.mkstemp(suffix=".nc", prefix="bisicles_flatten_")
        os.close(fd)
        tmp_nc = Path(_tmp)
        _tmp_path = tmp_nc

    file_info = parse_bisicles_filename(plot_file)

    # For standalone BISICLES files that don't match UKESM naming, detect
    # whether this is a CF time-mean file (plot.CF.NNNNNN.2d.hdf5) or a
    # snapshot (plot.NNNNNN.2d.hdf5) from the filename alone.
    _standalone_cf = (
        file_info is None and _RE_STANDALONE_CFMEAN.search(plot_file.name) is not None
    )

    try:
        if verbose:
            print(f"  Flattening {plot_file.name} onto level {level}...")
        run_flatten(plot_file, tmp_nc, exe_path=exe_path, level=level, x0=x0, y0=y0)

        if verbose:
            print(f"  Reading flattened data...")
        flatten_data = _read_flatten_nc(tmp_nc)
        if xygrid is not None:
            flatten_data = _regrid(flatten_data, xygrid)
    finally:
        if _tmp_path is not None and _tmp_path.exists():
            os.unlink(_tmp_path)

    # Store the CF-mean flag so write_cmip7_per_variable_netcdfs can set
    # correct time encoding without needing to re-examine the filename.
    flatten_data["_is_cf_mean"] = _standalone_cf

    return flatten_data, file_info


# ---------------------------------------------------------------------------
# High-level convenience functions
# ---------------------------------------------------------------------------

def process_plotfile(
    plot_file,
    output_dir,
    exe_path,
    level=0,
    epsg_code=None,
    x0=None,
    y0=None,
    reference_year=1850,
    calendar="gregorian",
    cmip7_only=False,
    ice_density=918.0,
    water_density=1028.0,
    h_min=1.0,
    keep_intermediate=False,
    grid_spec=None,
    verbose=True,
    **nc_kwargs,
):
    """
    Flatten a BISICLES plot HDF5 file and write one CMIP7/CF-compliant NetCDF
    per variable into *output_dir*.

    Each output file is named ``{cmip7_name}.nc`` and contains a single
    data variable with a time dimension of size 1 (the single input timestep).

    Parameters
    ----------
    plot_file : str or Path
        Input BISICLES plot HDF5 file.
    output_dir : str or Path
        Directory for output NetCDF files.  Created if it does not exist.
    exe_path : str or Path
        Full path to the flatten executable.
    level : int
        AMR level to flatten onto (0 = coarsest grid).
    epsg_code : int, optional
        EPSG code for the projection.  Overrides value from HDF5 metadata.
    x0 : float, optional
        X-coordinate of the lower-left corner of the domain (metres).
        Defaults to the standard BISICLES value for the given ``epsg_code``
        (see :data:`cf_utils.UKESM_GRID_ORIGINS`).
    y0 : float, optional
        Y-coordinate of the lower-left corner of the domain (metres).
        Defaults to the standard BISICLES value for the given ``epsg_code``
        (see :data:`cf_utils.UKESM_GRID_ORIGINS`).
     grid_spec : tuple (x_min, x_max, nx, y_min, y_max, ny)
        specify an output grid on different to the BISICLES grid but on the same
        projection . This will be a regular nx by ny grid,
        with x_min <= x <= x_max & similar for y
    reference_year : int
        Reference year for the CF time axis.
    calendar : str
        CF calendar name: ``"gregorian"`` (default) or ``"360_day"`` for
        UKESM-coupled runs.
    cmip7_only : bool
        If True, only write CMIP7-named variables.
    ice_density : float
        Ice density in kg m-3 for the flotation-based mask fallback (default 918.0).
    water_density : float
        Ocean water density in kg m-3 for the flotation-based mask fallback
        (default 1028.0).
    h_min : float
        Minimum ice thickness threshold in metres for the mask/iceFrac fallback
        (default 1.0 m).
    keep_intermediate : bool
        If True, keep the intermediate flatten NetCDF file (written to
        *output_dir* with suffix ``_flatten_raw.nc``).
    verbose : bool
        Print progress messages.
    **nc_kwargs
        Keyword arguments for :func:`write_cmip7_per_variable_netcdfs`
        (institution, source, experiment, variant_label, ice_sheet, extra_attrs).

    Returns
    -------
    list of Path
        Paths of output NetCDF files written.
    """
    plot_file = Path(plot_file)
    output_dir = Path(output_dir)

    intermediate_nc = None
    if keep_intermediate:
        output_dir.mkdir(parents=True, exist_ok=True)
        intermediate_nc = output_dir / (plot_file.stem.replace(".2d", "") + "_flatten_raw.nc")

    # Parse simulation time from the filename (UKESM BISICLES internal
    # time is always 0, so the filename is the authoritative source).
    file_info = parse_bisicles_filename(plot_file)
    if file_info is not None and verbose:
        kind = "time-mean" if file_info.is_time_mean else "snapshot"
        print(
            f"  Parsed filename: ice_sheet={file_info.ice_sheet}, "
            f"time={file_info.time_years:.2f} yr ({kind})"
        )

    # Override ice_sheet from filename when not set explicitly
    if file_info is not None and not nc_kwargs.get("ice_sheet"):
        nc_kwargs["ice_sheet"] = file_info.ice_sheet

    # Fall back to standard BISICLES grid origins when x0/y0 not supplied
    if (x0 is None or y0 is None) and epsg_code in UKESM_GRID_ORIGINS:
        origin = UKESM_GRID_ORIGINS[epsg_code]
        if x0 is None:
            x0 = origin["x0"]
        if y0 is None:
            y0 = origin["y0"]
        if verbose:
            print(f"  Using standard BISICLES grid origin: x0={x0}, y0={y0}")

    if verbose:
        print(f"Flattening {plot_file.name} onto level {level}...")

    xygrid = None
    if grid_spec is not None:
        print (f'output grid  x_min, x_max, nx, y_min, y_max, ny = {grid_spec}')
        x_min, x_max, nx, y_min, y_max, ny = grid_spec
        xygrid = np.linspace(x_min, x_max, nx), np.linspace(y_min, y_max, ny)

        
    flatten_data, _ = _flatten_plot_file(
        plot_file,
        exe_path=exe_path,
        level=level,
        x0=x0,
        y0=y0,
        xygrid=xygrid,
        keep_intermediate=keep_intermediate,
        intermediate_nc=intermediate_nc,
        verbose=verbose,
    )

    if verbose:
        bisicles_vars = list(flatten_data["variables"].keys())
        print(f"  Variables found: {bisicles_vars}")
        print(f"Writing per-variable CMIP7 NetCDF files to {output_dir}...")

    output_files = write_cmip7_per_variable_netcdfs(
        [(flatten_data, file_info)],
        output_dir,
        epsg_code=epsg_code,
        x0=x0,
        y0=y0,
        reference_year=reference_year,
        calendar=calendar,
        cmip7_only=cmip7_only,
        ice_density=ice_density,
        water_density=water_density,
        h_min=h_min,
        source_files=file_info.filename_pattern if file_info is not None else plot_file.name,
        grid_spec=grid_spec,
        **nc_kwargs,
    )

    if verbose:
        print(f"Done. Wrote {len(output_files)} variable file(s).")

    return output_files


def process_directory(
    directory,
    exe_path,
    output_dir=None,
    plot_pattern="plot.*.2d.hdf5",
    level=0,
    epsg_code=None,
    x0=None,
    y0=None,
    grid_spec=None,
    reference_year=1850,
    calendar="gregorian",
    cmip7_only=False,
    ice_density=918.0,
    water_density=1028.0,
    h_min=1.0,
    verbose=True,
    **nc_kwargs,
):
    """
    Flatten all BISICLES plot files in a directory and write one CMIP7 NetCDF
    per variable, with all timesteps on a single time axis.

    Each output file is named ``{cmip7_name}.nc`` and its ``time``
    dimension spans all input plot files in chronological order.

    Parameters
    ----------
    directory : str or Path
        Directory containing BISICLES plot HDF5 files.
    exe_path : str or Path
        Full path to the flatten executable.
    output_dir : str or Path, optional
        Directory to write output NetCDF files.  Defaults to the same directory
        as the input files.
    plot_pattern : str
        Glob pattern for plot files.
    level : int
        AMR level to flatten onto.
    epsg_code : int, optional
        EPSG code for the projection.
    x0 : float, optional
        X-coordinate of the lower-left corner of the domain (metres).
        Defaults to the standard BISICLES value for the given ``epsg_code``.
    y0 : float, optional
        Y-coordinate of the lower-left corner of the domain (metres).
        Defaults to the standard BISICLES value for the given ``epsg_code``.
    grid_spec : tuple (x_min, x_max, nx, y_min, y_max, ny)
        specify an output grid on different to the BISICLES grid but on the same
        projection . This will be a regular nx by ny grid,
        with x_min <= x <= x_max & similar for y
    reference_year : int
        Reference year for the CF time axis.
    calendar : str
        CF calendar name: ``"gregorian"`` (default) or ``"360_day"`` for
        UKESM-coupled runs.
    cmip7_only : bool
        If True, only write CMIP7-named variables.
    ice_density : float
        Ice density in kg m-3 for the flotation-based mask fallback (default 918.0).
    water_density : float
        Ocean water density in kg m-3 for the flotation-based mask fallback
        (default 1028.0).
    h_min : float
        Minimum ice thickness threshold in metres for the mask/iceFrac fallback
        (default 1.0 m).
    verbose : bool
        Print progress messages.
    **nc_kwargs
        Passed to :func:`write_cmip7_per_variable_netcdfs`.

    Returns
    -------
    list of Path
        Paths of output NetCDF files written (one per variable).
    """
    from .diagnostics import find_plot_files

    directory = Path(directory)
    if output_dir is None:
        output_dir = directory
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    # Fall back to standard BISICLES grid origins when x0/y0 not supplied
    if (x0 is None or y0 is None) and epsg_code in UKESM_GRID_ORIGINS:
        origin = UKESM_GRID_ORIGINS[epsg_code]
        if x0 is None:
            x0 = origin["x0"]
        if y0 is None:
            y0 = origin["y0"]
        if verbose:
            print(f"Using standard BISICLES grid origin: x0={x0}, y0={y0}")

    plot_files = find_plot_files(directory, pattern=plot_pattern)
    if verbose:
        print(f"Found {len(plot_files)} plot files in {directory}")

    # Collect flattened data from all plot files before writing
    all_data = []

    xygrid = None
    if grid_spec is not None:
        print (f'output grid  x_min, x_max, nx, y_min, y_max, ny = {grid_spec}')
        x_min, x_max, nx, y_min, y_max, ny = grid_spec
        xygrid = np.linspace(x_min, x_max, nx), np.linspace(y_min, y_max, ny)

    for i, pf in enumerate(plot_files):
        if verbose:
            print(f"  [{i+1}/{len(plot_files)}] Flattening {pf.name}...")

        # Parse filename for time info; override ice_sheet from first parseable file
        file_info = parse_bisicles_filename(pf)
        if file_info is not None and not nc_kwargs.get("ice_sheet"):
            nc_kwargs["ice_sheet"] = file_info.ice_sheet

        flatten_data, file_info = _flatten_plot_file(
            pf,
            exe_path=exe_path,
            level=level,
            x0=x0,
            y0=y0,
            xygrid=xygrid,
            verbose=False,
        )
        all_data.append((flatten_data, file_info))

    if verbose:
        print(f"Writing per-variable CMIP7 NetCDF files to {output_dir}...")

    output_files = write_cmip7_per_variable_netcdfs(
        all_data,
        output_dir,
        epsg_code=epsg_code,
        x0=x0,
        y0=y0,
        reference_year=reference_year,
        calendar=calendar,
        cmip7_only=cmip7_only,
        ice_density=ice_density,
        water_density=water_density,
        h_min=h_min,
        source_files=next(
            (fi.filename_pattern for _, fi in all_data if fi is not None),
            plot_files[0].name if plot_files else None,
        ),
        grid_spec=grid_spec,
        **nc_kwargs,
    )

    if verbose:
        print(f"Done. Wrote {len(output_files)} variable file(s).")

    return output_files
