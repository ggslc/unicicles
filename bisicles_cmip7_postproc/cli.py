"""
Command-line interface for bisicles_cmip7_postproc.

Provides three entry points:

  bike-cmip7-postproc-diagnostics
    Run the BISICLES diagnostics tool on a plot file or directory of plot files
    and write one CF-compliant scalar timeseries NetCDF per variable.

  bike-cmip7-postproc-flatten
    Flatten BISICLES plot HDF5 file(s) onto a uniform grid and write one
    CMIP7/CF-compliant 2D spatial NetCDF per variable.

  bike-cmip7-postproc-run
    Run either workflow from a YAML or JSON config file.  All options that can
    be passed on the command line can be specified in the config file instead.

Usage examples
--------------
Single plot file (diagnostics – one timestep per variable):

    bike-cmip7-postproc-diagnostics \\
        --input  /run/output/plot.000050.2d.hdf5 \\
        --output-dir /run/postproc/ \\
        --ice-sheet GrIS --institution "University of Bristol"

Full run directory (multi-year timeseries, one file per variable):

    bike-cmip7-postproc-diagnostics \\
        --input  /run/output/ \\
        --output-dir /run/postproc/ \\
        --ice-sheet GrIS --experiment historical --variant-label r1i1p1f3

Single plot file (2D spatial fields – one file per variable):

    bike-cmip7-postproc-flatten \\
        --input  /run/output/plot.000050.2d.hdf5 \\
        --output-dir /run/postproc/ \\
        --level 2 --epsg 3413 \\
        --ice-sheet GrIS --institution "University of Bristol"

Full run directory (multi-year 2D fields, one file per variable):

    bike-cmip7-postproc-flatten \\
        --input  /run/output/ \\
        --output-dir /run/postproc/ \\
        --level 2 --epsg 3413 \\
        --ice-sheet GrIS --experiment historical --variant-label r1i1p1f3
"""

import argparse
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Shared argument parsing helpers
# ---------------------------------------------------------------------------

def _add_metadata_args(parser):
    """Add CMIP7 metadata arguments shared by both subcommands."""
    grp = parser.add_argument_group("CF/CMIP7 metadata")
    grp.add_argument("--institution", default="", help="Institution name.")
    grp.add_argument(
        "--source",
        default="BISICLES adaptive mesh refinement ice sheet model",
        help="Source model description.",
    )
    grp.add_argument("--experiment", default="", help="Experiment identifier.")
    grp.add_argument(
        "--variant-label", default="", dest="variant_label",
        help="CMIP variant label (e.g. r1i1p1f3).",
    )
    grp.add_argument(
        "--ice-sheet", default="", dest="ice_sheet",
        help="Ice sheet identifier (e.g. GrIS, AIS).",
    )
    grp.add_argument(
        "--frequency", default="", dest="frequency",
        help=(
            "CMIP output frequency string (e.g. 'yr', 'mon', 'day'). "
            "Auto-detected from the filename period field when not set."
        ),
    )
    grp.add_argument(
        "--reference-year", type=int, default=1850, dest="reference_year",
        help="Reference year for the CF time axis (default: 1850).",
    )
    grp.add_argument(
        "--calendar", default="gregorian", dest="calendar",
        choices=["standard", "gregorian", "360_day"],
        help=(
            "CF calendar for the time axis (default: gregorian, 365.25 days/year). "
            "'standard' is accepted as a synonym for 'gregorian'. "
            "Use '360_day' for UKESM-coupled runs (360 days/year)."
        ),
    )


def _add_exe_args(parser, tool_name):
    """Add executable path argument."""
    parser.add_argument(
        "--exe-path", required=True, dest="exe_path",
        metavar="PATH",
        help=f"Full path to the BISICLES {tool_name} executable file.",
    )


# ---------------------------------------------------------------------------
# diagnostics subcommand
# ---------------------------------------------------------------------------

def _build_diagnostics_parser():
    p = argparse.ArgumentParser(
        prog="bike-cmip7-postproc-diagnostics",
        description=(
            "Run BISICLES diagnostics on a plot file or run directory and "
            "write a CF-compliant scalar timeseries NetCDF."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument(
        "--config", metavar="FILE", default=None,
        help=(
            "YAML or JSON config file.  All options can be set in the config "
            "file; any options given on the command line override the config."
        ),
    )
    p.add_argument(
        "--input", "-i", required=False, default=None, metavar="PATH",
        help="BISICLES plot HDF5 file or directory of plot files.",
    )
    p.add_argument(
        "--output-dir", "-o", required=False, default=None, metavar="DIR",
        dest="output_dir",
        help=(
            "Output directory for per-variable CF NetCDF files.  "
            "One file per diagnostic variable is written into this directory.  "
            "Defaults to the directory containing the input file(s)."
        ),
    )
    p.add_argument(
        "--plot-pattern", default="plot.*.2d.hdf5", dest="plot_pattern",
        help=(
            "Glob pattern for plot files when --input is a directory "
            "(default: 'plot.*.2d.hdf5')."
        ),
    )
    _add_exe_args(p, "diagnostics")

    phys = p.add_argument_group("physical constants")
    phys.add_argument("--ice-density", type=float, default=918.0, dest="ice_density",
                      help="Ice density in kg m-3 (default: 918.0).")
    phys.add_argument("--water-density", type=float, default=1028.0, dest="water_density",
                      help="Ocean water density in kg m-3 (default: 1028.0).")
    phys.add_argument("--gravity", type=float, default=9.81,
                      help="Gravitational acceleration in m s-2 (default: 9.81).")
    phys.add_argument("--h-min", type=float, default=1.0, dest="h_min",
                      help="Minimum ice thickness threshold in m (default: 1.0).")

    mask = p.add_argument_group("regional mask (optional)")
    mask.add_argument("--mask-file", default=None, dest="mask_file", metavar="HDF5",
                      help="HDF5 mask file for regional diagnostics.")
    mask.add_argument("--mask-no-start", type=int, default=0, dest="mask_no_start",
                      help="First mask region index (default: 0).")
    mask.add_argument("--mask-no-end", type=int, default=0, dest="mask_no_end",
                      help="Last mask region index (default: 0).")

    _add_metadata_args(p)
    p.add_argument("--quiet", "-q", action="store_true",
                   help="Suppress progress messages.")
    return p


def run_diagnostics_cli(args=None):
    """Entry point for the ``bike-cmip7-postproc-diagnostics`` command."""
    print("SLC diagnostics")
    # Pre-scan for --config before full parsing so we can set defaults
    import sys as _sys
    _raw = args if args is not None else _sys.argv[1:]
    parser = _build_diagnostics_parser()
    if "--config" in _raw:
        idx = _raw.index("--config")
        _cfg_path = _raw[idx + 1]
        from .config import load_config, _normalise_diag_cfg
        _cfg = load_config(_cfg_path)
        _cfg.pop("tool", None)
        _normalise_diag_cfg(_cfg)
        parser.set_defaults(**_cfg)
    ns = parser.parse_args(args)

    if ns.input is None:
        parser.error("--input is required (or set 'input' in a --config file)")

    from .diagnostics import process_single_file, process_directory

    nc_kwargs = dict(
        institution=ns.institution,
        source=ns.source,
        experiment=ns.experiment,
        variant_label=ns.variant_label,
        ice_sheet=ns.ice_sheet,
        frequency=ns.frequency,
    )

    input_path = Path(ns.input)

    # Default output directory to the location of the input file(s)
    output_dir = ns.output_dir if ns.output_dir is not None else (
        str(input_path) if input_path.is_dir() else str(input_path.parent)
    )

    if input_path.is_dir():
        process_directory(
            directory=input_path,
            output_dir=output_dir,
            exe_path=ns.exe_path,
            plot_pattern=ns.plot_pattern,
            ice_density=ns.ice_density,
            water_density=ns.water_density,
            gravity=ns.gravity,
            h_min=ns.h_min,
            mask_file=ns.mask_file,
            mask_no_start=ns.mask_no_start,
            mask_no_end=ns.mask_no_end,
            reference_year=ns.reference_year,
            calendar=ns.calendar,
            verbose=not ns.quiet,
            **nc_kwargs,
        )
    elif input_path.is_file():
        process_single_file(
            plot_file=input_path,
            output_dir=output_dir,
            exe_path=ns.exe_path,
            ice_density=ns.ice_density,
            water_density=ns.water_density,
            gravity=ns.gravity,
            h_min=ns.h_min,
            mask_file=ns.mask_file,
            mask_no_start=ns.mask_no_start,
            mask_no_end=ns.mask_no_end,
            reference_year=ns.reference_year,
            calendar=ns.calendar,
            **nc_kwargs,
        )
    else:
        print(f"Error: {ns.input} is not a file or directory.", file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# flatten subcommand
# ---------------------------------------------------------------------------

def _build_flatten_parser():
    p = argparse.ArgumentParser(
        prog="bike-cmip7-postproc-flatten",
        description=(
            "Flatten BISICLES plot HDF5 file(s) onto a uniform grid and write "
            "CMIP7/CF-compliant 2D spatial NetCDF file(s)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument(
        "--config", metavar="FILE", default=None,
        help=(
            "YAML or JSON config file.  All options can be set in the config "
            "file; any options given on the command line override the config."
        ),
    )
    p.add_argument(
        "--input", "-i", required=False, default=None, metavar="PATH",
        help="BISICLES plot HDF5 file or directory of plot files.",
    )
    p.add_argument(
        "--output-dir", "-o", metavar="DIR", default=None, dest="output_dir",
        help=(
            "Output directory for per-variable NetCDF files.  "
            "One file per CMIP7 variable is written into this directory, "
            "named {cmip7_name}.nc.  "
            "Defaults to the same directory as the input file(s)."
        ),
    )
    p.add_argument(
        "--plot-pattern", default="plot.*.2d.hdf5", dest="plot_pattern",
        help="Glob pattern for plot files in directory mode (default: 'plot.*.2d.hdf5').",
    )
    _add_exe_args(p, "flatten")
    p.add_argument(
        "--level", type=int, default=0,
        help=(
            "AMR level to flatten onto. 0 = coarsest grid; higher values give "
            "finer resolution. Must be >= 0. (default: 0)"
        ),
    )
    p.add_argument(
        "--epsg", type=int, default=None, dest="epsg_code",
        help=(
            "EPSG code for the projection (e.g. 3413 for GrIS, 3031 for AIS). "
            "Overrides the value read from the HDF5 file metadata."
        ),
    )
    p.add_argument(
        "--x0", type=float, default=None, dest="x0",
        help=(
            "X-coordinate of the lower-left corner of the domain in metres. "
            "Defaults to the standard BISICLES grid value for the chosen EPSG code "
            "(GrIS/3413: -654650.0; AIS/3031: -3072000.0)."
        ),
    )
    p.add_argument(
        "--y0", type=float, default=None, dest="y0",
        help=(
            "Y-coordinate of the lower-left corner of the domain in metres. "
            "Defaults to the standard BISICLES grid value for the chosen EPSG code "
            "(GrIS/3413: -3385950.0; AIS/3031: -3072000.0)."
        ),
    )
    p.add_argument(
        "--x_min", type=float, default=None, dest="x_min",
        help= ("x coordinate of the lower-left corner of the output domain"))           
    p.add_argument(
        "--x_max", type=float, default=None, dest="x_max",
        help= ("x coordinate of the upper-right corner of the output domain"))
    p.add_argument(
        "--y_min", type=float, default=None, dest="y_min",
        help= ("y coordinate of the lower-left corner of the output domain")) 
    p.add_argument(
        "--y_max", type=float, default=None, dest="y_max",
        help= ("x coordinate of the upper right corner of the output domain")) 
    p.add_argument(
        "--nx", type=int, default=None, dest="nx",
        help= ("ouput grid x dimension")) 
    p.add_argument(
        "--ny", type=int, default=None, dest="ny",
        help= ("ouput grid y dimension")) 

    p.add_argument(
        "--cmip7-only", action="store_true", dest="cmip7_only",
        help=(
            "Only write CMIP7-standard variables. By default, unmapped BISICLES "
            "variables are also written with their original names."
        ),
    )
    p.add_argument(
        "--keep-intermediate", action="store_true", dest="keep_intermediate",
        help="Keep the intermediate raw flatten NetCDF file (single-file mode only).",
    )

    geom = p.add_argument_group(
        "flotation mask fallback",
        description=(
            "Physical constants used when the plot file does not contain an "
            "explicit mask or iceFrac variable.  In that case the grounded/"
            "floating partition is derived from the flotation criterion using "
            "ice thickness and bed topography."
        ),
    )
    geom.add_argument(
        "--ice-density", type=float, default=918.0, dest="ice_density",
        help="Ice density in kg m-3 (default: 918.0).",
    )
    geom.add_argument(
        "--water-density", type=float, default=1028.0, dest="water_density",
        help="Ocean water density in kg m-3 (default: 1028.0).",
    )
    geom.add_argument(
        "--h-min", type=float, default=1.0, dest="h_min",
        help=(
            "Minimum ice thickness in metres below which a cell is treated as "
            "ice-free when deriving the mask or iceFrac (default: 1.0 m)."
        ),
    )

    _add_metadata_args(p)
    p.add_argument("--quiet", "-q", action="store_true",
                   help="Suppress progress messages.")
    return p


def run_flatten_cli(args=None):
    """Entry point for the ``bike-cmip7-postproc-flatten`` command."""
    print("SLC flatten")
    import sys as _sys
    _raw = args if args is not None else _sys.argv[1:]
    parser = _build_flatten_parser()
    if "--config" in _raw:
        idx = _raw.index("--config")
        _cfg_path = _raw[idx + 1]
        from .config import load_config, _normalise_flatten_cfg
        _cfg = load_config(_cfg_path)
        _cfg.pop("tool", None)
        _normalise_flatten_cfg(_cfg)
        parser.set_defaults(**_cfg)
    ns = parser.parse_args(args)

    if ns.input is None:
        parser.error("--input is required (or set 'input' in a --config file)")

    grid_spec = None
    grid_spec_test = (ns.x_min, ns.x_max, ns.nx, ns.y_min, ns.y_max, ns.ny)
    if any(grid_spec_test):
        if not all(grid_spec_test):
            parser.error("must specify all of --x_min, --x_max, --nx, --y_min, --y_max, --ny, or none")
        else:
            grid_spec = grid_spec_test

    from .flatten import process_plotfile, process_directory

    nc_kwargs = dict(
        institution=ns.institution,
        source=ns.source,
        experiment=ns.experiment,
        variant_label=ns.variant_label,
        ice_sheet=ns.ice_sheet,
        frequency=ns.frequency,
    )

    input_path = Path(ns.input)

    # Default output directory to the location of the input file(s)
    output_dir = ns.output_dir if ns.output_dir is not None else (
        str(input_path) if input_path.is_dir() else str(input_path.parent)
    )

    if input_path.is_dir():
        process_directory(
            directory=input_path,
            output_dir=output_dir,
            plot_pattern=ns.plot_pattern,
            exe_path=ns.exe_path,
            level=ns.level,
            epsg_code=ns.epsg_code,
            x0=ns.x0,
            y0=ns.y0,
            reference_year=ns.reference_year,
            calendar=ns.calendar,
            cmip7_only=ns.cmip7_only,
            ice_density=ns.ice_density,
            water_density=ns.water_density,
            h_min=ns.h_min,
            grid_spec=grid_spec,
            verbose=not ns.quiet,
            **nc_kwargs,
        )
    elif input_path.is_file():
        process_plotfile(
            plot_file=input_path,
            output_dir=output_dir,
            exe_path=ns.exe_path,
            level=ns.level,
            epsg_code=ns.epsg_code,
            x0=ns.x0,
            y0=ns.y0,
            reference_year=ns.reference_year,
            calendar=ns.calendar,
            cmip7_only=ns.cmip7_only,
            ice_density=ns.ice_density,
            water_density=ns.water_density,
            h_min=ns.h_min,
            keep_intermediate=ns.keep_intermediate,
            grid_spec=grid_spec,
            verbose=not ns.quiet,
            **nc_kwargs,
        )
    else:
        print(f"Error: {ns.input} is not a file or directory.", file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# bike-cmip7-postproc-run: config-file entry point
# ---------------------------------------------------------------------------

def _build_run_parser():
    p = argparse.ArgumentParser(
        prog="bike-cmip7-postproc-run",
        description=(
            "Run a BISICLES post-processing workflow from a YAML or JSON config "
            "file.  The config file must contain a 'tool' key set to either "
            "'flatten' or 'diagnostics'.  All other options correspond to the "
            "keyword arguments of the underlying Python functions."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example config file (YAML):\n\n"
            "    tool: flatten\n"
            "    input: /scratch/cx209c/run/output/\n"
            "    output_dir: /scratch/cx209c/postproc/\n"
            "    epsg: 3413\n"
            "    level: 2\n"
            "    calendar: 360_day\n"
            "    ice_density: 918.0\n"
            "    institution: University of Bristol\n"
            "    experiment: historical\n"
            "    variant_label: r1i1p1f3\n"
            "    ice_sheet: GrIS\n"
        ),
    )
    p.add_argument(
        "config", metavar="CONFIG_FILE",
        help="Path to a YAML or JSON config file.",
    )
    p.add_argument(
        "--set", nargs=2, action="append", metavar=("KEY", "VALUE"),
        default=[],
        help=(
            "Override a config value from the command line.  Can be given "
            "multiple times.  Values are auto-cast to int/float/bool where "
            "possible.  Example: --set calendar 360_day --set level 2"
        ),
    )
    return p


def _cast_override(value):
    """Try to cast a string override value to int, float, or bool."""
    if value.lower() == "true":
        return True
    if value.lower() == "false":
        return False
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    return value


def run_config_cli(args=None):
    """Entry point for the ``bike-cmip7-postproc-run`` command."""
    parser = _build_run_parser()
    ns = parser.parse_args(args)

    overrides = {k: _cast_override(v) for k, v in ns.set}

    from .config import run_from_config
    run_from_config(ns.config, overrides=overrides if overrides else None)


# ---------------------------------------------------------------------------
# ISMIP7 standalone: shared metadata args
# ---------------------------------------------------------------------------

def _add_ismip7_metadata_args(parser):
    """Add ISMIP7 DRS metadata arguments for the ISM submission filename and global attributes."""
    grp = parser.add_argument_group("ISMIP7 DRS metadata")
    grp.add_argument(
        "--source-id", default="", dest="source_id",
        help="Modelling group name for ISMIP7 DRS (e.g. 'BristolGlaciology'). [required]",
    )
    grp.add_argument(
        "--ism-id", default="BISICLES", dest="ism_id",
        help="ISM name and version for ISMIP7 DRS (default: 'BISICLES').",
    )
    grp.add_argument(
        "--ism-member-id", default="m001", dest="ism_member_id",
        help="ISM choice variant identifier for ISMIP7 DRS (default: 'm001').",
    )
    grp.add_argument(
        "--esm-id", default="standalone", dest="esm_id",
        help=(
            "CMIP ESM name used for forcing, for ISMIP7 DRS "
            "(e.g. 'CESM2-WACCM'; default: 'standalone' for idealised forcing)."
        ),
    )
    grp.add_argument(
        "--forcing-member-id", default="f001", dest="forcing_member_id",
        help="Forcing choice variant identifier for ISMIP7 DRS (default: 'f001').",
    )
    grp.add_argument(
        "--set-counter", default="C001", dest="set_counter",
        help="Set counter for ISMIP7 DRS (e.g. 'C001', 'E001'; default: 'C001').",
    )
    grp.add_argument(
        "--group", default="", dest="group",
        help="Modelling group name for ISMIP7 global attribute 'group'.",
    )
    grp.add_argument(
        "--contact-name", default="", dest="contact_name",
        help="Contact person name(s) for ISMIP7 mandatory global attribute.",
    )
    grp.add_argument(
        "--contact-email", default="", dest="contact_email",
        help="Contact person email(s) for ISMIP7 mandatory global attribute.",
    )


# ---------------------------------------------------------------------------
# bike-ismip7-postproc-diagnostics
# ---------------------------------------------------------------------------

def _build_ismip7_diagnostics_parser():
    p = _build_diagnostics_parser()
    p.prog = "bike-ismip7-postproc-diagnostics"
    p.description = (
        "Run BISICLES diagnostics on a standalone plot file or run directory and "
        "write ISMIP7 DRS-named CF-compliant scalar timeseries NetCDF files."
    )
    _add_ismip7_metadata_args(p)
    # ISMIP7 requires "standard" calendar; override the CMIP7 default of "gregorian"
    p.set_defaults(calendar="standard")
    return p


def run_ismip7_diagnostics_cli(args=None):
    """Entry point for the ``bike-ismip7-postproc-diagnostics`` command."""
    import sys as _sys
    _raw = args if args is not None else _sys.argv[1:]
    parser = _build_ismip7_diagnostics_parser()
    if "--config" in _raw:
        idx = _raw.index("--config")
        _cfg_path = _raw[idx + 1]
        from .config import load_config, _normalise_diag_cfg
        _cfg = load_config(_cfg_path)
        _cfg.pop("tool", None)
        _normalise_diag_cfg(_cfg)
        parser.set_defaults(**_cfg)
    ns = parser.parse_args(args)

    if ns.input is None:
        parser.error("--input is required (or set 'input' in a --config file)")

    from .diagnostics import process_single_file, process_directory

    nc_kwargs = dict(
        institution=ns.institution,
        source=ns.source,
        experiment=ns.experiment,
        variant_label=ns.variant_label,
        ice_sheet=ns.ice_sheet,
        frequency=ns.frequency or "yr",
        ismip7_mode=True,
        source_id=ns.source_id,
        ism_id=ns.ism_id,
        ism_member_id=ns.ism_member_id,
        esm_id=ns.esm_id,
        forcing_member_id=ns.forcing_member_id,
        set_counter=ns.set_counter,
        group=ns.group,
        contact_name=ns.contact_name,
        contact_email=ns.contact_email,
    )

    input_path = Path(ns.input)
    output_dir = ns.output_dir if ns.output_dir is not None else (
        str(input_path) if input_path.is_dir() else str(input_path.parent)
    )

    if input_path.is_dir():
        process_directory(
            directory=input_path,
            output_dir=output_dir,
            exe_path=ns.exe_path,
            plot_pattern=ns.plot_pattern,
            ice_density=ns.ice_density,
            water_density=ns.water_density,
            gravity=ns.gravity,
            h_min=ns.h_min,
            mask_file=ns.mask_file,
            mask_no_start=ns.mask_no_start,
            mask_no_end=ns.mask_no_end,
            reference_year=ns.reference_year,
            calendar=ns.calendar,
            verbose=not ns.quiet,
            **nc_kwargs,
        )
    elif input_path.is_file():
        process_single_file(
            plot_file=input_path,
            output_dir=output_dir,
            exe_path=ns.exe_path,
            ice_density=ns.ice_density,
            water_density=ns.water_density,
            gravity=ns.gravity,
            h_min=ns.h_min,
            mask_file=ns.mask_file,
            mask_no_start=ns.mask_no_start,
            mask_no_end=ns.mask_no_end,
            reference_year=ns.reference_year,
            calendar=ns.calendar,
            **nc_kwargs,
        )
    else:
        print(f"Error: {ns.input} is not a file or directory.", file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# bike-ismip7-postproc-flatten
# ---------------------------------------------------------------------------

def _build_ismip7_flatten_parser():
    p = _build_flatten_parser()
    p.prog = "bike-ismip7-postproc-flatten"
    p.description = (
        "Flatten standalone BISICLES plot HDF5 file(s) onto a uniform grid and "
        "write ISMIP7 DRS-named CF-compliant 2D spatial NetCDF files."
    )
    _add_ismip7_metadata_args(p)
    # ISMIP7 requires "standard" calendar; override the CMIP7 default of "gregorian"
    p.set_defaults(calendar="standard")
    return p


def run_ismip7_flatten_cli(args=None):
    """Entry point for the ``bike-ismip7-postproc-flatten`` command."""
    import sys as _sys
    _raw = args if args is not None else _sys.argv[1:]
    parser = _build_ismip7_flatten_parser()
    if "--config" in _raw:
        idx = _raw.index("--config")
        _cfg_path = _raw[idx + 1]
        from .config import load_config, _normalise_flatten_cfg
        _cfg = load_config(_cfg_path)
        _cfg.pop("tool", None)
        _normalise_flatten_cfg(_cfg)
        parser.set_defaults(**_cfg)
    ns = parser.parse_args(args)

    if ns.input is None:
        parser.error("--input is required (or set 'input' in a --config file)")

    from .flatten import process_plotfile, process_directory

    nc_kwargs = dict(
        institution=ns.institution,
        source=ns.source,
        experiment=ns.experiment,
        variant_label=ns.variant_label,
        ice_sheet=ns.ice_sheet,
        frequency=ns.frequency or "yr",
        ismip7_mode=True,
        source_id=ns.source_id,
        ism_id=ns.ism_id,
        ism_member_id=ns.ism_member_id,
        esm_id=ns.esm_id,
        forcing_member_id=ns.forcing_member_id,
        set_counter=ns.set_counter,
        group=ns.group,
        contact_name=ns.contact_name,
        contact_email=ns.contact_email,
    )

    input_path = Path(ns.input)
    output_dir = ns.output_dir if ns.output_dir is not None else (
        str(input_path) if input_path.is_dir() else str(input_path.parent)
    )

    
    grid_spec=None
    grid_spec_test = (ns.x_min, ns.x_max, ns.nx, ns.y_min, ns.y_max, ns.ny)
    if any(grid_spec_test):
       if not all(grid_spec_test):
           parser.error("must specify all of --x_min, --x_max, nx, --y_min, --y_max, ny, or none")
       else:
           grid_spec = grid_spec_test

           
    if input_path.is_dir():
        process_directory(
            directory=input_path,
            output_dir=output_dir,
            plot_pattern=ns.plot_pattern,
            exe_path=ns.exe_path,
            level=ns.level,
            epsg_code=ns.epsg_code,
            x0=ns.x0,
            y0=ns.y0,
            grid_spec=grid_spec,
            reference_year=ns.reference_year,
            calendar=ns.calendar,
            cmip7_only=ns.cmip7_only,
            ice_density=ns.ice_density,
            water_density=ns.water_density,
            h_min=ns.h_min,
            verbose=not ns.quiet,
            **nc_kwargs,
        )
    elif input_path.is_file():
        process_plotfile(
            plot_file=input_path,
            output_dir=output_dir,
            exe_path=ns.exe_path,
            level=ns.level,
            epsg_code=ns.epsg_code,
            x0=ns.x0,
            y0=ns.y0,
            grid_spec=grid_spec,
            reference_year=ns.reference_year,
            calendar=ns.calendar,
            cmip7_only=ns.cmip7_only,
            ice_density=ns.ice_density,
            water_density=ns.water_density,
            h_min=ns.h_min,
            keep_intermediate=ns.keep_intermediate,
            verbose=not ns.quiet,
            **nc_kwargs,
        )
    else:
        print(f"Error: {ns.input} is not a file or directory.", file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# bike-ismip7-postproc-run: config-file entry point
# ---------------------------------------------------------------------------

def run_ismip7_config_cli(args=None):
    """Entry point for the ``bike-ismip7-postproc-run`` command."""
    parser = _build_run_parser()
    parser.prog = "bike-ismip7-postproc-run"
    parser.description = (
        "Run a standalone BISICLES ISMIP7 post-processing workflow from a YAML "
        "or JSON config file.  The config file must contain a 'tool' key set to "
        "either 'flatten' or 'diagnostics'.  ISMIP7 DRS output naming is applied "
        "automatically; set model_id and member_id in the config or use --set."
    )
    ns = parser.parse_args(args)

    overrides = {k: _cast_override(v) for k, v in ns.set}
    overrides["ismip7_mode"] = True  # always force ISMIP7 mode from this entry point

    from .config import run_from_config
    run_from_config(ns.config, overrides=overrides)


# ---------------------------------------------------------------------------
# Script entry points
# ---------------------------------------------------------------------------

def main_diagnostics():
    print ("main_diag")
    run_diagnostics_cli()


def main_flatten():
    print ("main_flatten")
    run_flatten_cli()


def main_run():
    run_config_cli()


def main_ismip7_diagnostics():
    run_ismip7_diagnostics_cli()


def main_ismip7_flatten():
    run_ismip7_flatten_cli()


def main_ismip7_run():
    print ("main_ismip7_run")
    run_ismip7_config_cli()


if __name__ == "__main__":
    # Allow running as: python -m bisicles_cmip7_postproc.cli {diagnostics|flatten|run} ...
    if len(sys.argv) > 1 and sys.argv[1] == "diagnostics":
        sys.argv.pop(1)
        run_diagnostics_cli()
    elif len(sys.argv) > 1 and sys.argv[1] == "flatten":
        sys.argv.pop(1)
        run_flatten_cli()
    elif len(sys.argv) > 1 and sys.argv[1] == "run":
        sys.argv.pop(1)
        run_config_cli()
    else:
        print("Usage: python -m bisicles_cmip7_postproc.cli {diagnostics|flatten|run} [args]")
        sys.exit(1)
