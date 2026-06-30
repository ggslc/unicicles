"""
Wrapper for the BISICLES diagnostics file tool and post-processing utilities
for writing CF-compliant scalar timeseries NetCDF files.

The diagnostics tool computes integrated quantities (total ice volume, area,
mass balance fluxes, etc.) from a BISICLES plot HDF5 file and writes them as
CSV rows.  This module:
  1. Locates and calls the diagnostics executable.
  2. Parses the CSV output into Python data structures.
  3. Assembles a time series from multiple plot files and writes a CF/CMIP7-
     compliant NetCDF file with a time dimension.

CSV output format from the diagnostics tool:
  csvheader,filename,time,maskNo,region,quantity,unit,value
  csvdata,<file>,<time>,<maskNo>,<region>,<quantity>,<unit>,<value>

Regions: entire | grounded | floating | ice | nonice
Quantities (units):
  volume (m3), volumeAbove (m3), fracArea (m2), area (m2)
  SMB (m3/a), BMB (m3/a), dhdt (m3/a), calving (m3/a),
  discharge (m3/a), flxDivFile (m3/a), flxDivReconstr (m3/a)
"""

import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from .cmip7_vars import SCALAR_MAPPING, ICE_DENSITY, SECS_PER_YEAR
from .cf_utils import (
    CF_CONVENTIONS,
    FILL_VALUE,
    ISMIP7_FILL_VALUE,
    get_global_attributes,
    period_to_cmip_frequency,
    add_time_variable,
    add_crs_variable,
    add_time_bounds,
    _ismip7_drs_filename,
)
from .filename_parser import parse_bisicles_filename


# ---------------------------------------------------------------------------
# Running the diagnostics tool
# ---------------------------------------------------------------------------

def run_diagnostics(
    plot_file,
    exe_path,
    out_file=None,
    append=False,
    ice_density=918.0,
    water_density=1028.0,
    gravity=9.81,
    h_min=1.0,
    mask_file=None,
    mask_no_start=0,
    mask_no_end=0,
):
    """
    Run the BISICLES diagnostics tool on a single plot file.

    Parameters
    ----------
    plot_file : str or Path
        Path to the BISICLES plot HDF5 file.
    exe_path : str or Path
        Full path to the diagnostics executable.
    out_file : str or Path, optional
        Path for the CSV output.  A temporary file is created if not given.
    append : bool
        If True, append to an existing out_file rather than overwriting.
    ice_density : float
        Ice density in kg m-3 (default 918.0).
    water_density : float
        Ocean water density in kg m-3 (default 1028.0).
    gravity : float
        Gravitational acceleration in m s-2 (default 9.81).
    h_min : float
        Minimum ice thickness threshold in metres (default 1.0).
    mask_file : str or Path, optional
        Path to an HDF5 mask file for regional diagnostics.
    mask_no_start : int
        First mask region index (default 0 = whole domain).
    mask_no_end : int
        Last mask region index (default 0 = whole domain only).

    Returns
    -------
    str
        Path to the CSV output file.

    Raises
    ------
    RuntimeError
        If the diagnostics executable exits with a non-zero return code.
    """
    if exe_path is None:
        raise ValueError(
            "exe_path is required. Provide the full path to the diagnostics executable, e.g.:\n"
            "  --exe-path /path/to/filetools/diagnostics2d.Linux.64.g++.gfortran.OPT.ex"
        )
    p = Path(exe_path)
    if p.is_dir():
        raise ValueError(
            f"exe_path points to a directory, not an executable: {exe_path}\n"
            "Provide the full path to the executable file, e.g.:\n"
            "  /path/to/filetools/diagnostics2d.Linux.64.g++.gfortran.OPT.ex"
        )
    if not p.is_file():
        raise FileNotFoundError(f"Diagnostics executable not found: {exe_path}")
    exe = str(p)

    _tmp_created = out_file is None
    if _tmp_created:
        fd, out_file = tempfile.mkstemp(suffix=".csv", prefix="bisicles_diag_")
        os.close(fd)

    args = [
        exe,
        f"plot_file={plot_file}",
        f"out_file={out_file}",
        f"ice_density={ice_density}",
        f"water_density={water_density}",
        f"gravity={gravity}",
        f"h_min={h_min}",
        f"mask_no_start={mask_no_start}",
        f"mask_no_end={mask_no_end}",
    ]
    if append:
        args.append("-append")
    if mask_file is not None:
        args.append(f"mask_file={mask_file}")

    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"diagnostics tool failed (return code {result.returncode}):\n"
            f"stdout: {result.stdout}\n"
            f"stderr: {result.stderr}"
        )

    return out_file


# ---------------------------------------------------------------------------
# Parsing diagnostics CSV output
# ---------------------------------------------------------------------------

def parse_diagnostics_csv(csv_path):
    """
    Parse the CSV output produced by the BISICLES diagnostics tool.

    Parameters
    ----------
    csv_path : str or Path
        Path to the CSV file.

    Returns
    -------
    list of dict
        Each dict has keys:
          'filename', 'time' (float, years), 'maskNo' (int),
          'region' (str), 'quantity' (str), 'unit' (str), 'value' (float)
    """
    records = []
    with open(csv_path) as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith("csvdata,"):
                continue
            parts = line.split(",")
            if len(parts) != 8:
                continue
            records.append(
                {
                    "filename": parts[1],
                    "time": float(parts[2]),
                    "maskNo": int(parts[3]),
                    "region": parts[4],
                    "quantity": parts[5],
                    "unit": parts[6],
                    "value": float(parts[7]),
                }
            )
    return records


# ---------------------------------------------------------------------------
# Finding BISICLES plot files in a directory
# ---------------------------------------------------------------------------

def find_plot_files(directory, pattern="plot.*.2d.hdf5"):
    """
    Find and sort BISICLES plot files in a directory.

    Plot files follow the naming convention ``plot.XXXXXX.2d.hdf5`` where
    XXXXXX is a zero-padded step counter.  Sorting by name gives time order.

    Parameters
    ----------
    directory : str or Path
        Directory to search.
    pattern : str
        Glob pattern for plot file names (default ``plot.*.2d.hdf5``).

    Returns
    -------
    list of Path
        Sorted list of plot file paths.
    """
    directory = Path(directory)
    files = sorted(directory.glob(pattern))
    if not files:
        raise FileNotFoundError(
            f"No files matching '{pattern}' found in {directory}"
        )
    return files


# ---------------------------------------------------------------------------
# Writing scalar timeseries to CF-compliant NetCDF
# ---------------------------------------------------------------------------

def _build_timeseries(records, filename_to_info=None):
    """
    Organise diagnostic records into a dict of time-indexed arrays.

    Parameters
    ----------
    records : list of dict
        Parsed CSV records from :func:`parse_diagnostics_csv`.
    filename_to_info : dict, optional
        Maps plot-file basename -> :class:`~.filename_parser.BISICLESFileInfo`.
        When supplied, the simulation time for each record is taken from the
        corresponding ``BISICLESFileInfo.time_years`` value rather than from the
        CSV ``time`` field (which is always 0 in UKESM-coupled runs).

    Returns
    -------
    times : sorted list of float
        Nominal time coordinate values (fractional years).
    data  : dict mapping (cmip7_name, maskNo) -> {time: value, ...}
    meta  : dict mapping cmip7_name -> SCALAR_MAPPING entry
    time_bounds : dict mapping nominal_time -> (start_years, end_years) or None
        Populated for time-mean files; None entries for snapshot files.
    """
    times_set = set()
    data = {}
    meta = {}
    time_bounds = {}  # nominal_time -> (start_years, end_years)

    for r in records:
        key = (r["region"], r["quantity"])
        if key not in SCALAR_MAPPING:
            continue
        m = SCALAR_MAPPING[key]
        cname = m["cmip7_name"]
        mask_no = r["maskNo"]

        # Determine nominal time from filename metadata if available
        if filename_to_info is not None:
            import os
            basename = os.path.basename(r["filename"])
            info = filename_to_info.get(basename)
        else:
            info = None

        if info is not None:
            t = info.time_years
            if info.is_time_mean:
                time_bounds[t] = (info.start_time_years, info.end_time_years)
            else:
                time_bounds[t] = None
        else:
            t = r["time"]
            time_bounds[t] =  (int(t)-(1),int(t))#(int(t),int(t)+1)
        
        times_set.add(t)
        
        var_key = (cname, mask_no)
        if var_key not in data:
            data[var_key] = {}
            meta[cname] = m
        data[var_key][t] = r["value"] * m["conversion_factor"]

    print("time_set", sorted(times_set))

    return sorted(times_set), data, meta, time_bounds


def write_diagnostics_netcdf(
    records,
    output_dir,
    reference_year=1850,
    calendar="standard",
    ice_sheet="",
    institution="",
    source="BISICLES adaptive mesh refinement ice sheet model",
    experiment="",
    variant_label="",
    extra_attrs=None,
    filename_to_info=None,
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
):
    """
    Write one CF-compliant scalar timeseries NetCDF per variable from parsed
    diagnostic records.

    Each output file is named ``{cmip7_name}.nc`` (or
    ``{cmip7_name}_mask{N}.nc`` for regional-mask variants) inside
    *output_dir* and contains a single data variable with a ``time`` dimension
    spanning all timesteps.

    Parameters
    ----------
    records : list of dict
        Parsed diagnostic records as returned by :func:`parse_diagnostics_csv`.
    output_dir : str or Path
        Directory for output NetCDF files.  Created if it does not exist.
    reference_year : int
        Reference year for the CF time axis (default 1850).
    calendar : str
        CF calendar name: ``"gregorian"`` (default) or ``"360_day"`` for
        UKESM-coupled runs.
    ice_sheet : str
        Ice sheet identifier written to global attributes (e.g. 'GrIS').
    institution : str
        Institution name for global attributes.
    source : str
        Source model description.
    experiment : str
        Experiment name.
    variant_label : str
        CMIP variant label (e.g. 'r1i1p1f3').
    extra_attrs : dict, optional
        Additional global attributes to set.
    filename_to_info : dict, optional
        Maps plot-file basename -> :class:`~.filename_parser.BISICLESFileInfo`.
        When supplied, simulation times and time bounds are derived from the
        filenames rather than from the CSV ``time`` field (which is always 0
        in UKESM-coupled runs).  Time-mean files will have a ``time_bnds``
        variable to satisfy CF-1.12 requirements.
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
    from netCDF4 import Dataset  # imported here to keep the module importable without netCDF4

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    times, data, meta, time_bounds = _build_timeseries(records, filename_to_info=filename_to_info)
    if not times:
        raise ValueError(
            "No CMIP7-mapped quantities found in the diagnostic records. "
            "Check that the diagnostics tool produced output and that "
            "SCALAR_MAPPING covers the desired (region, quantity) pairs."
        )
    
    print (times, time_bounds)
    if not frequency and filename_to_info:
        first_info = next(iter(filename_to_info.values()))
        if first_info.period:
            frequency = period_to_cmip_frequency(first_info.period)

    # In ISMIP7 mode, populate group/model from DRS fields when not set explicitly
    _group = group or source_id
    _model = ism_id

    epsg_code = '3031' 
    if ice_sheet == 'AIS':
        epsg_code = '3031'
    elif ice_sheet == 'GrIS':
        epsg_code = '3413'
        
    crs_str=f'epsg:{epsg_code}'

    global_attrs = get_global_attributes(
        institution=institution,
        source=source,
        experiment=experiment,
        variant_label=variant_label,
        ice_sheet=ice_sheet,
        source_files=next(iter(filename_to_info.values())).filename_pattern if filename_to_info else None,
        frequency=frequency,
        conventions="CF-1.12" if ismip7_mode else CF_CONVENTIONS,
        model_id=model_id,
        member_id=member_id,
        crs=crs_str,
        group=_group if ismip7_mode else group,
        contact_name=contact_name,
        contact_email=contact_email,
    )
    if ismip7_mode and _model:
        global_attrs["model"] = _model
    if extra_attrs:
        global_attrs.update(extra_attrs)

    # ISMIP7 requires f4 precision and the netCDF4 default fill value;
    # UKESM/CMIP7 output uses f8 precision and the conventional 1e20 fill value.
    _fill_value = ISMIP7_FILL_VALUE if ismip7_mode else FILL_VALUE
    _time_dtype = "f4" if ismip7_mode else "f8"
    _var_dtype  = "f4" if ismip7_mode else "f8"

    time_arr = np.asarray(times)

    # Determine whether any timesteps are time-mean (need time_bnds)
    has_time_bounds = any(v is not None for v in time_bounds.values())

    if has_time_bounds:
        start_years = []
        end_years = []
        for t in times:
            bnds = time_bounds.get(t)
            if bnds is not None:
                start_years.append(bnds[0])
                end_years.append(bnds[1])
            else:
                # Snapshot: degenerate bounds equal to the time point
                start_years.append(t)
                end_years.append(t)

    output_files = []

    # Write one file per (variable, mask) combination
    for (cname, mask_no), time_dict in sorted(data.items()):
        m = meta[cname]

        # File and variable name
        print(m['ismip7_type'])
        
        var_name = cname if mask_no == 0 else f"{cname}_mask{mask_no}"
        if ismip7_mode:
            out_path = output_dir / _ismip7_drs_filename(
                cname, ice_sheet, source_id, ism_id, ism_member_id,
                esm_id, forcing_member_id, experiment, set_counter,
                time_arr, mask_no=mask_no,
            )
        else:
            out_path = output_dir / f"{var_name}.nc"

        values = np.array([time_dict.get(t, np.nan) for t in times])

        with Dataset(str(out_path), "w", format="NETCDF4") as ds:
            ds.setncatts(global_attrs)
            
            add_crs_variable(ds, epsg_code, x0=None, y0=None)
            
            ds.variable_id = cname
            ds.variable_name = cname
            ds.title = m["long_name"] if ismip7_mode else f"UniCiCles (BISICLES) output from UKESM: {m['long_name']}"

            # Time dimension (unlimited record dimension) and variable
            ds.createDimension("time", None)
            
            time_bodge = 0.5
            if m['ismip7_type'] == 'ST':
                time_bodge = 1.0
            add_time_variable(ds, time_arr+time_bodge, reference_year=reference_year,
                              calendar=calendar, dtype=_time_dtype)

            # Time bounds for time-mean data
            if has_time_bounds:
                add_time_bounds(ds, start_years, end_years,
                                reference_year=reference_year, calendar=calendar,
                                dtype=_time_dtype)

            var = ds.createVariable(
                var_name, _var_dtype, ("time",), fill_value=_fill_value
            )
            _cast = np.float32 if _var_dtype == "f4" else np.float64
            var[:] = np.where(np.isnan(values), _fill_value, values).astype(_cast)
            var.standard_name = m["standard_name"]
            var.long_name = m["long_name"]
            var.units = m["cmip7_units"]
            if "cell_methods" in m:
                if has_time_bounds:
                    var.cell_methods = f"time: mean {m['cell_methods']}"
                else:
                    var.cell_methods = m["cell_methods"]
            if "comment" in m:
                var.comment = m["comment"]
            if mask_no != 0:
                var.mask_number = mask_no
                var.mask_comment = (
                    f"Restricted to drainage-basin mask region {mask_no}."
                )

        output_files.append(out_path)

    return output_files


# ---------------------------------------------------------------------------
# High-level convenience functions
# ---------------------------------------------------------------------------

def _check_not_cf_mean(plot_file):
    """Raise ValueError if *plot_file* is a time-mean CF-output file.

    The BISICLES diagnostics tool operates on instantaneous plot files and
    does not produce meaningful results for time-averaged CF-output files
    (``plot.CF-*.hdf5``).  These should be processed with the flatten
    workflow instead.
    """
    info = parse_bisicles_filename(plot_file)
    if info is not None and info.is_time_mean:
        raise ValueError(
            f"The diagnostics tool cannot be used with time-mean CF-output "
            f"files.\n"
            f"  File: {Path(plot_file).name}\n"
            f"Time-mean files (plot.CF-*.hdf5) contain spatially gridded data "
            f"averaged over a coupling window and are not compatible with the "
            f"BISICLES diagnostics executable.\n"
            f"Use 'bisicles-flatten' (or process_plotfile / process_directory "
            f"from bisicles_cmip7_postproc.flatten) to process this file."
        )


def process_single_file(
    plot_file,
    output_dir,
    exe_path,
    ice_density=918.0,
    water_density=1028.0,
    gravity=9.81,
    h_min=1.0,
    mask_file=None,
    mask_no_start=0,
    mask_no_end=0,
    reference_year=1850,
    calendar="gregorian",
    **nc_kwargs,
):
    """
    Run diagnostics on a single plot file and write one CF timeseries NetCDF
    per variable into *output_dir*.

    Parameters
    ----------
    plot_file : str or Path
    output_dir : str or Path
        Directory for output NetCDF files.  Created if it does not exist.
    exe_path : str or Path
        Full path to the diagnostics executable.
    ice_density, water_density, gravity, h_min : float
    mask_file : str or Path, optional
    mask_no_start, mask_no_end : int
    reference_year : int
    calendar : str
        CF calendar name: ``"gregorian"`` (default) or ``"360_day"`` for
        UKESM-coupled runs.
    **nc_kwargs
        Passed to :func:`write_diagnostics_netcdf`
        (institution, source, experiment, variant_label, ice_sheet, extra_attrs).

    Returns
    -------
    list of Path
        Paths of output NetCDF files written.
    """
    plot_file = Path(plot_file)
    _check_not_cf_mean(plot_file)
    file_info = parse_bisicles_filename(plot_file)
    
    filename_to_info = {plot_file.name: file_info} if file_info is not None else None

    # Override ice_sheet from filename when not set explicitly
    if file_info is not None and not nc_kwargs.get("ice_sheet"):
        nc_kwargs["ice_sheet"] = file_info.ice_sheet

    csv_file = None
    try:
        csv_file = run_diagnostics(
            plot_file,
            exe_path=exe_path,
            ice_density=ice_density,
            water_density=water_density,
            gravity=gravity,
            h_min=h_min,
            mask_file=mask_file,
            mask_no_start=mask_no_start,
            mask_no_end=mask_no_end,
        )
        records = parse_diagnostics_csv(csv_file)
        return write_diagnostics_netcdf(
            records,
            output_dir,
            reference_year=reference_year,
            calendar=calendar,
            filename_to_info=filename_to_info,
            **nc_kwargs,
        )
    finally:
        if csv_file and Path(csv_file).exists():
            os.unlink(csv_file)


def process_directory(
    directory,
    output_dir,
    exe_path,
    plot_pattern="plot.*.2d.hdf5",
    ice_density=918.0,
    water_density=1028.0,
    gravity=9.81,
    h_min=1.0,
    mask_file=None,
    mask_no_start=0,
    mask_no_end=0,
    reference_year=1850,
    calendar="gregorian",
    verbose=True,
    **nc_kwargs,
):
    """
    Run diagnostics on all plot files in a directory and write one CF
    timeseries NetCDF per variable, with all timesteps on a single time axis.

    Parameters
    ----------
    directory : str or Path
        Directory containing BISICLES plot HDF5 files.
    output_dir : str or Path
        Directory for output NetCDF files.  Created if it does not exist.
    exe_path : str or Path
        Full path to the diagnostics executable.
    plot_pattern : str
        Glob pattern for plot files (default ``plot.*.2d.hdf5``).
    ice_density, water_density, gravity, h_min : float
        Physical constants passed to the diagnostics tool.
    mask_file : str or Path, optional
        Optional HDF5 mask file for regional diagnostics.
    mask_no_start, mask_no_end : int
        Mask region index range.
    reference_year : int
        Reference year for the CF time axis.
    calendar : str
        CF calendar name: ``"gregorian"`` (default) or ``"360_day"`` for
        UKESM-coupled runs.
    verbose : bool
        Print progress messages.
    **nc_kwargs
        Passed to :func:`write_diagnostics_netcdf`.

    Returns
    -------
    list of Path
        Paths of output NetCDF files written.
    """
    plot_files = find_plot_files(directory, pattern=plot_pattern)
    if verbose:
        print(f"S Found {len(plot_files)} plot files in {directory}")

    # Build filename -> BISICLESFileInfo mapping for UKESM-style filenames
    filename_to_info = {}
    for pf in plot_files:
        info = parse_bisicles_filename(pf)
  
        if info is not None:
            filename_to_info[pf.name] = info
            
    if not filename_to_info:
        filename_to_info = None  # Fall back to CSV time values

    # Override ice_sheet from first parseable filename if not set
    if filename_to_info and not nc_kwargs.get("ice_sheet"):
        first_info = next(iter(filename_to_info.values()))
        nc_kwargs["ice_sheet"] = first_info.ice_sheet

    csv_file = None
    try:
        fd, csv_file = tempfile.mkstemp(suffix=".csv", prefix="bisicles_diag_timeseries_")
        os.close(fd)

        for i, pf in enumerate(plot_files):
            _check_not_cf_mean(pf)
            if verbose:
                print(f"  [{i+1}/{len(plot_files)}] {pf.name}")
            run_diagnostics(
                pf,
                exe_path=exe_path,
                out_file=csv_file,
                append=(i > 0),
                ice_density=ice_density,
                water_density=water_density,
                gravity=gravity,
                h_min=h_min,
                mask_file=mask_file,
                mask_no_start=mask_no_start,
                mask_no_end=mask_no_end,
            )

        if verbose:
            print(f"Parsing diagnostics and writing per-variable files to {output_dir}")

        records = parse_diagnostics_csv(csv_file)
        output_files = write_diagnostics_netcdf(
            records,
            output_dir,
            reference_year=reference_year,
            calendar=calendar,
            filename_to_info=filename_to_info,
            **nc_kwargs,
        )

        if verbose:
            print(f"Done. Wrote {len(output_files)} variable file(s).")

        return output_files
    finally:
        if csv_file and Path(csv_file).exists():
            os.unlink(csv_file)
