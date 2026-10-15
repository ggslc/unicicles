#!/usr/bin/env bash

set -e

echo "check for follow on dump"
if [ -e $PREICE_ICE_DUMP_HDF ]; then
  export BISICLES_STATE=$PREICE_ICE_DUMP_HDF
else
  echo "expected BISICLES restart $PREICE_ICE_DUMP_HDF, NOT FOUND"
  exit 1
fi

if [ -e $PREICE_ICE_DUMP_NC ]; then
  export GLINT_STATE=$PREICE_ICE_DUMP_NC
else
  echo "expected Glint restart $PREICE_ICE_DUMP_NC, NOT FOUND"
  exit 1
fi

echo "bisicles starting from $BISICLES_STATE"
echo "Glint starting from $GLINT_STATE"

echo "Build the Glint and BISICLES namelists"
#construct the INI namelist for bisicles from the rose-app.conf
if [ -e $BISICLES_CONFIG ]; then rm $BISICLES_CONFIG; fi
rose config --file=rose-app-run.conf bisicles_config | rose env-cat >> $BISICLES_CONFIG;
echo " " >> $BISICLES_CONFIG;

#construct the INI namelist for glimmer/glint from the rose-app.conf
if [ -e $GLINT_CONFIG ]; then rm $GLINT_CONFIG; fi
for section in "grid" "sigma" "GLINT climate" "projection" "options" \
               "external_dycore_options" "time" "parameters"         \
               "CF output" "CF input"; do
  echo ["$section"] >> $GLINT_CONFIG;
  rose config --file=rose-app-run.conf "$section" | rose env-cat >> $GLINT_CONFIG;
  echo " " >> $GLINT_CONFIG;
done

#construct the top.config namelist for glimmer/glint from the rose-app.conf
if [ -e top.config ]; then rm top.config; fi
for section in "GLINT" "GLINT instance"; do
  echo ["$section"] >> top.config;
  rose config --file=rose-app-run.conf "$section" | rose env-cat >> top.config;
  echo " " >> top.config;
done

#launch unicicles
echo "Run unicicles"
rose mpi-launch unicicles

echo "rename and copy the output needed for other stages"
#rename and copy the dump/output we need to $DATAM
LAST_PLOT_HDF=`ls -1r plot.*[0-9].2d.hdf5 | head -1`
LAST_DUMP_HDF=`ls -1r chk.*hdf5 | head -1`
echo "final BISICLES output are $LAST_PLOT_HDF, $LAST_DUMP_HDF"
#
mv $LAST_PLOT_HDF $POSTICE_ICE_PLOT_HDF
mv $LAST_DUMP_HDF $POSTICE_ICE_DUMP_HDF
cdo -seltimestep,-1 *dummy.hot.nc temp.nc
#cdo on archer2 works in a container that doesn't know much of the filesystem
mv temp.nc $POSTICE_ICE_DUMP_NC


#check for CF plot files
shopt -s nullglob
files=( plot.*CF.2d.hdf5 )
if (( ${#files[@]} )); then 
  LAST_PLOT_CF_HDF=`ls -1r plot.*CF.2d.hdf5 | head -1`
  echo "also found CF plot files, $LAST_PLOT_CF_HDF"
  mv $LAST_PLOT_CF_HDF ${POSTICE_ICE_PLOT_CF_HDF}
  #mv plot.*CF.2d.hdf5 ${UNICICLES_DATAM}
fi
