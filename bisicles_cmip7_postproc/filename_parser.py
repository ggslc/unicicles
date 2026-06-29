"""
Parse BISICLES/UKESM plot filenames to extract simulation time and metadata.

When BISICLES is coupled to UKESM, the model's internal time counter is reset
to zero at the start of each annual coupling window.  The actual simulation
time and ice-sheet identity must therefore be derived from the filename rather
than from the HDF5 file contents.

Three filename conventions are supported:

1. Old format (snapshot)::

       bisicles_cx209c_18510101_plot-AIS.hdf5

   A snapshot plotfile produced at the *end* of running from 18500101 to
   18510101.  The date in the filename is the snapshot date.

2. New format (snapshot)::

       bisicles_dx030c_1y_18880101-18890101_plot-AIS.hdf5

   A snapshot plotfile at the *end* of the time window
   18880101-18890101.  The snapshot date is the *end* date.

3. New format CF-output (time mean)::

       bisicles_dx030c_1y_18880101-18890101_plot.CF-AIS.hdf5

   Data contain time-averaged values over the period 18880101-18890101.
   The nominal time coordinate is the midpoint of the period; time bounds
   span the full period.
"""

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Regex patterns for the three filename formats
# ---------------------------------------------------------------------------

# Group names: suite_id, date (YYYYMMDD)
_RE_OLD = re.compile(
    r"bisicles_(?P<suite_id>[a-z0-9]+)_"
    r"(?P<date>\d{8})_"
    r"plot-(?P<ice_sheet>[A-Z]+(?:[a-z]*[A-Z]*)*)\.hdf5$",
    re.IGNORECASE,
)

# Group names: suite_id, period, start_date, end_date, ice_sheet
_RE_NEW_SNAPSHOT = re.compile(
    r"bisicles_(?P<suite_id>[a-z0-9]+)_"
    r"(?P<period>[^_]+)_"
    r"(?P<start_date>\d{8})-(?P<end_date>\d{8})_"
    r"plot-(?P<ice_sheet>[A-Z]+(?:[a-z]*[A-Z]*)*)\.hdf5$",
    re.IGNORECASE,
)

# Group names: suite_id, period, start_date, end_date, ice_sheet
_RE_NEW_CFMEAN = re.compile(
    r"bisicles_(?P<suite_id>[a-z0-9]+)_"
    r"(?P<period>[^_]+)_"
    r"(?P<start_date>\d{8})-(?P<end_date>\d{8})_"
    r"plot\.CF-(?P<ice_sheet>[A-Z]+(?:[a-z]*[A-Z]*)*)\.hdf5$",
    re.IGNORECASE,
)

_RE_MATT= re.compile(
    r"plot-(?P<ice_sheet>[A-Z]+(?:[a-z]*[A-Z]*)*)\.hdf5$",
    re.IGNORECASE,
)

# ---------------------------------------------------------------------------
# Helper: convert YYYYMMDD integer to a fractional year
# ---------------------------------------------------------------------------

def _yyyymmdd_to_year(yyyymmdd: int) -> float:
    """
    Convert a date integer (YYYYMMDD) to a fractional year.

    Uses a simple 365.25-day year, consistent with the CF time axis
    convention used elsewhere in this package.

    Parameters
    ----------
    yyyymmdd : int
        Date as an 8-digit integer, e.g. 18510101.

    Returns
    -------
    float
        Fractional year, e.g. 1851.0 for 18510101.
    """
    year = yyyymmdd // 10000
    month = (yyyymmdd // 100) % 100
    day = yyyymmdd % 100

    # Days elapsed since Jan 1 of this year (1-based day-of-year)
    _days_in_month = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    doy = sum(_days_in_month[:month]) + day  # day-of-year (1 = Jan 1)
    return year + (doy - 1) / 365.25


# ---------------------------------------------------------------------------
# BISICLESFileInfo dataclass
# ---------------------------------------------------------------------------

@dataclass
class BISICLESFileInfo:
    """
    Metadata extracted from a BISICLES/UKESM plot filename.

    Attributes
    ----------
    suite_id : str
        Rose suite identifier (e.g. 'cx209c', 'dx030c').
    ice_sheet : str
        Ice sheet identifier (e.g. 'AIS', 'GrIS').
    start_date : int
        Start date of the coupling window as YYYYMMDD (0 for old format).
    end_date : int
        End date (snapshot date for formats 1 and 2; end of averaging
        period for format 3).
    is_time_mean : bool
        True if the file contains time-averaged data (format 3).
    period : str
        Period string from the filename (e.g. '1y') or '' for old format.
    """

    suite_id: str
    ice_sheet: str
    start_date: int
    end_date: int
    is_time_mean: bool
    period: str = ""

    # ------------------------------------------------------------------
    # Derived properties
    # ------------------------------------------------------------------

    @property
    def snapshot_date(self) -> int:
        """End date as YYYYMMDD (the date of the snapshot or end of mean)."""
        return self.end_date

    @property
    def start_time_years(self) -> float:
        """Start of the coupling window as a fractional year."""
        return _yyyymmdd_to_year(self.start_date) if self.start_date else _yyyymmdd_to_year(self.end_date)

    @property
    def end_time_years(self) -> float:
        """End of the coupling window (snapshot date) as a fractional year."""
        return _yyyymmdd_to_year(self.end_date)

    @property
    def time_years(self) -> float:
        """
        Nominal time coordinate in fractional years.

        For snapshot files this is the snapshot date (``end_time_years``).
        For time-mean files this is the midpoint of the averaging period.
        """
        if self.is_time_mean and self.start_date:
            return 0.5 * (self.start_time_years + self.end_time_years)
        return self.end_time_years

    @property
    def cell_methods_time(self) -> str:
        """
        CF ``cell_methods`` time component describing how the data were
        sampled.

        Returns ``'time: point'`` for snapshots and
        ``'time: mean'`` for time-averaged outputs.
        """
        return "time: mean" if self.is_time_mean else "time: point"

    @property
    def filename_pattern(self) -> str:
        """
        Filename pattern for this file series with date fields replaced by
        YYYYMMDD placeholders, e.g.
        ``bisicles_dx030c_1y_YYYYMMDD-YYYYMMDD_plot-AIS.hdf5``.
        """
        if not self.period:  # old format
            return f"bisicles_{self.suite_id}_YYYYMMDD_plot-{self.ice_sheet}.hdf5"
        elif self.is_time_mean:
            return (
                f"bisicles_{self.suite_id}_{self.period}_YYYYMMDD-YYYYMMDD"
                f"_plot.CF-{self.ice_sheet}.hdf5"
            )
        else:
            return (
                f"bisicles_{self.suite_id}_{self.period}_YYYYMMDD-YYYYMMDD"
                f"_plot-{self.ice_sheet}.hdf5"
            )


# ---------------------------------------------------------------------------
# Public parsing function
# ---------------------------------------------------------------------------

def parse_bisicles_filename(path) -> Optional[BISICLESFileInfo]:
    """
    Extract simulation metadata from a BISICLES/UKESM plot filename.

    Parameters
    ----------
    path : str or Path
        Full path or bare filename of the plot file.

    Returns
    -------
    BISICLESFileInfo or None
        Parsed metadata, or ``None`` if the filename does not match any
        known BISICLES/UKESM convention.

    Examples
    --------
    >>> info = parse_bisicles_filename(
    ...     "bisicles_cx209c_18510101_plot-AIS.hdf5")
    >>> info.ice_sheet
    'AIS'
    >>> info.time_years   # end of 1850-01-01 → 1851-01-01 window
    1851.0
    >>> info.is_time_mean
    False

    >>> info = parse_bisicles_filename(
    ...     "bisicles_dx030c_1y_18880101-18890101_plot.CF-AIS.hdf5")
    >>> info.is_time_mean
    True
    >>> round(info.time_years, 4)   # midpoint of 1888-1889
    1888.5
    """
    name = Path(path).name

    # Format 3: time-mean CF output (check before snapshot — .CF. is distinctive)
    m = _RE_NEW_CFMEAN.match(name)
    if m:
        return BISICLESFileInfo(
            suite_id=m.group("suite_id"),
            ice_sheet=m.group("ice_sheet").upper(),
            start_date=int(m.group("start_date")),
            end_date=int(m.group("end_date")),
            is_time_mean=True,
            period=m.group("period"),
        )

    # Format 2: new snapshot
    m = _RE_NEW_SNAPSHOT.match(name)
    if m:
        return BISICLESFileInfo(
            suite_id=m.group("suite_id"),
            ice_sheet=m.group("ice_sheet").upper(),
            start_date=int(m.group("start_date")),
            end_date=int(m.group("end_date")),
            is_time_mean=True,#False,
            period=m.group("period"),
        )

    # Format 1: old snapshot
    m = _RE_OLD.match(name)
    if m:
        return BISICLESFileInfo(
            suite_id=m.group("suite_id"),
            ice_sheet=m.group("ice_sheet").upper(),
            start_date=int(m.group("date"))-1,
            end_date=int(m.group("date")),
            is_time_mean=True,#False,
            period="",
        )

    return None
