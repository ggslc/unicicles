#!/bin/ksh
#=============================================
# Ensure that unicicles is setup properly
#=============================================
# Set the environment variables
#---------------------------------------------
coreCoupleName="_${RUNID}c_${UNICICLES_PERIOD}_${PRE_UNICICLES_DATE}-${START_UNICICLES_DATE}_icecouple.nc"
atmosInitCoupleFile="atmos${coreCoupleName}"
oceanInitCoupleFile="nemo${coreCoupleName}"

coreCalvingName="bisicles_${RUNID}c_${UNICICLES_PERIOD}_${PRE_UNICICLES_DATE}-${START_UNICICLES_DATE}_calving"
#---------------------------------------------
# Display environment variables
#---------------------------------------------
echo "Unicicles data directory:               ${UNICICLES_DATAM}"
echo "Atmosphere coupling file for start of unicicles: ${UNICICLES_ATMOS_COUPLE_START}"
echo "Ocean coupling file for start of unicicles:      ${UNICICLES_OCEAN_COUPLE_START}"
#---------------------------------------------
# Make sure that data directory exists
#---------------------------------------------
if test -d "${UNICICLES_DATAM}"; then
  echo "Note: directory ${UNICICLES_DATAM} already exists"
else
  echo "Note: creating ${UNICICLES_DATAM} ..."
  mkdir -p ${UNICICLES_DATAM}
fi
if test -d "${UNI_OASIS_RMP}"; then
  echo "Note: directory ${UNI_OASIS_RMP} already exists"
else
  echo "Note: creating ${UNI_OASIS_RMP} ..."
  mkdir -p ${UNI_OASIS_RMP}
fi
cd ${UNICICLES_DATAM}
#---------------------------------------------
# Create symbolic link to the initial atmosphere data
#---------------------------------------------
if test -f "${UNICICLES_ATMOS_COUPLE_START}"; then
  if test -h ${atmosInitCoupleFile}; then
    echo "Removing symbolic link ${atmosInitCoupleFile}"
    rm ${atmosInitCoupleFile}
  fi
  echo "Creating symbolic link from ${UNICICLES_ATMOS_COUPLE_START} to ${atmosInitCoupleFile} ..."
  ln -s ${UNICICLES_ATMOS_COUPLE_START} ${atmosInitCoupleFile}
else
  echo "ERROR: failed to find unicicles atmosphere file ${UNICICLES_ATMOS_COUPLE_START}"
  exit 1
fi
#---------------------------------------------
# Create symbolic links to the initial ocean data
#---------------------------------------------

if test -f "${NEMO_CALVING_START}"; then
  if test -h ${UNICICLES_CALVING_IN}.nc; then
    echo "Removing symbolic link ${UNICICLES_CALVING_IN}.nc"
    rm ${UNICICLES_CALVING_IN}.nc
  fi
  echo "Creating symbolic link from ${NEMO_CALVING_START} to ${UNICICLES_CALVING_IN}.nc"
  ln -s ${NEMO_CALVING_START} ${UNICICLES_CALVING_IN}.nc
else
  echo "ERROR: failed to find initial NEMO calving file ${NEMO_CALVING_START}. Not OK"
  exit 1
fi

if $UNICICLES_AIS; then
  if test -f "${UNICICLES_OCEAN_COUPLE_START}"; then
    if test -h ${oceanInitCoupleFile}; then
      echo "Removing symbolic link ${oceanInitCoupleFile}"
      rm ${oceanInitCoupleFile}
    fi
    echo "Creating symbolic link from ${UNICICLES_OCEAN_COUPLE_START} to ${oceanInitCoupleFile} ..."
    ln -s ${UNICICLES_OCEAN_COUPLE_START} ${oceanInitCoupleFile}
  else
    echo "ERROR: failed to find unicicles ocean file ${UNICICLES_OCEAN_COUPLE_START}. Not OK if running Antarctica."
    exit 1
  fi

  if test -f "${UNICICLES_CALVING_START_AIS}"; then
    if test -h ${coreCalvingName}-AIS.hdf5; then
      echo "Removing symbolic link ${coreCalvingName}-AIS.hdf5"
      rm ${coreCalvingName}-AIS.hdf5
    fi
    echo "Creating symbolic link from ${UNICICLES_CALVING_START_AIS} to ${coreCalvingName}-AIS.hdf5 ..."
    ln -s ${UNICICLES_CALVING_START_AIS} ${coreCalvingName}-AIS.hdf5
  else
    echo "ERROR: failed to find unicicles calving file ${UNICICLES_CALVING_START_AIS}. Not OK if running Antarctica."
    exit 1
  fi

  if test -f "${NEMO_BATHY_START}"; then
    if test -h ${UNICICLES_BATHY_IN}.nc; then
      echo "Removing symbolic link ${UNICICLES_BATHY_IN}.nc"
      rm ${UNICICLES_BATHY_IN}.nc
    fi
    echo "Creating symbolic link from ${NEMO_BATHY_START} to ${UNICICLES_BATHY_IN}.nc"
    ln -s ${NEMO_BATHY_START} ${UNICICLES_BATHY_IN}.nc
  else
    echo "ERROR: failed to find initial NEMO bathymetry file ${NEMO_BATHY_START}. Not OK if running Antarctica."
    exit 1
  fi

fi

if $UNICICLES_GRIS; then
  if test -f "${UNICICLES_CALVING_START_GRIS}"; then
    if test -h ${coreCalvingName}-GrIS.hdf5; then
      echo "Removing symbolic link ${coreCalvingName}-GrIS.hdf5"
      rm ${coreCalvingName}-GrIS.hdf5
    fi
    echo "Creating symbolic link from ${UNICICLES_CALVING_START_GRIS} to ${coreCalvingName}-GrIS.hdf5 ..."
    ln -s ${UNICICLES_CALVING_START_GRIS} ${coreCalvingName}-GrIS.hdf5
  else
    echo "ERROR: failed to find unicicles calving file ${UNICICLES_CALVING_START_GRIS}. Not OK if running Greenland."
    exit 1
  fi
fi

cd ${UNI_OASIS_RMP}
#---------------------------------------------
# Create symbolic link to the initial OASIS maps
#---------------------------------------------
if test -d "${START_UNICICLES_DATE}"; then
  echo "Note: directory ${UNI_OASIS_RMP}/${START_UNICICLES_DATE} already exists"
else
  echo "Note: linking ${UNI_OASIS_RMP}/${START_UNICICLES_DATE} ..."
  ln -s ${UNI_OASIS_RMP_START} ${START_UNICICLES_DATE} 
fi
