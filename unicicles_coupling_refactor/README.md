# unicicles_coupling_refactor
Intermediate version control for the UKESM-UniCiCles coupling code while it's undergoing refactoring.

Overview of new cylc tasks (and associated top-level scripts/programs):
1. ocean_to_regional_ice_sheet (nemo*_to_regional_bisicles.py)
2. ice_sheet (unicicles)
3. regional_ice_sheet_calving_to_ocean (regional_bisicles_calving_to_nemo.py)
4. ice_sheets_to_ocean (bisicles_global_to_nemo.py)
5. regional_ice_sheet_topo_to_ancil (regional_unicicles_topo_to_cap.py)
6. regional_ancil_topo (capOrog)
7. ice_sheets_ancils_to_atmos (unicicles_cap_global_to_um.py)
- coupled happens here
8. average_ocean_for_ice_sheets (average_nemo_for_bisicles.py)
9. atmos_to_ice_sheets (um_to_unicicles.py)

Quick tour of directories inside active/:
- python3/ is the latest version of the coupling code we're working on.
- python3_legacy/ is the original code before refactoring.
- python2_legacy/ is the python2 version of the original code (needed for Met Office HPC).
- fortran/ and shell/ are not being touched for now.
