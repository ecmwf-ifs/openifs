#!/usr/bin/env bash

# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

set -o pipefail
set +v

SCRIPTS_DIR="$( cd $( dirname "${BASH_SOURCE[0]}" ) && pwd -P )"
if [[ $(basename ${SCRIPTS_DIR}) == "bin" ]]; then
  SCRIPTS_DIR="$( cd $( dirname "${SCRIPTS_DIR}/$(readlink "${BASH_SOURCE[0]}")" ) && pwd -P )"
fi

source ${SCRIPTS_DIR}/ecwam_runtime.sh
source ${SCRIPTS_DIR}/ecwam_parse_commandline.sh
source ${SCRIPTS_DIR}/ecwam_helper_functions.sh

function cleanup() {
  shopt -s nullglob
  LOG_DIR=${RUN_DIR}/logs/create_bathymetry
  mkdir -p ${LOG_DIR}
  rm -rf ${LOG_DIR}/*
  mv ${WORK_DIR}/stdout.log ${LOG_DIR}/
  mv ${WORK_DIR}/stderr.log ${LOG_DIR}/

  unset -f echo
  echo
  echo "Log files have been copied to    ${LOG_DIR}"

  cd ${RUN_DIR}
  rm -rf ${WORK_DIR}

  shopt -u nullglob
}

#.############################################################################
#.
#. run_bathymetry generates unformatted bathymetry input files
#. for use by run_preproc 
#.
#. THE INPUT BATHYMETRY FILES CAN BE ANY OF 
#.  - ETOPO2: etopo2_2006apr.dat.gz
#.  - ETOPO1: ETOPO1_Ice_g_int.xyz.gz
#.
#.############################################################################

#. There are three possible ways to produce the bathymetry for WAM
#. It is controlled by llobstrct and wambathy.
#. If llobstrct is set to false then you will need to supply the mean bathymetry in the prescribed
#. format (see bathymetry_global_1.0). !!! This is the old way !!!.
#. You will NOT be using the unresolved bathymetry scheme that is now available
#. in this version of WAM.
#. If llobstrct is true (the default at ECMWF), then
#. The mean bathymetry and unresolved bathymetry attenuations coefficients will be created 
#. from the high resolution bathymetry data etopo2 (see data file etopo2_2006apr.dat.gz). 
#. (as obtained from the US NGDC web site. This is the version of April 2006).
#. Note that the input bathymetry is ETOP02, which means that it has a resolution of 2 arc minutes
#.
#. or with ETOPO1 (1 arc minute data set)
#. you will need to get it from NOAA web site:
#. https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/xyz/
#. 
# Bathymetry source (ETOPO1 or ETOPO2)
# For fine horizontal resolution (< 0.1 degree), use at least 36 directions !!!!!!!!
#
# The "reference_levels" input file is used to reset the reference level in the high resolution data set
# within specified areas (see create_wam_bathymetry). It is used to blank some land areas
# that would otherwise be interpreted as below sea level (hence sea) and sea areas that
# are outside the area of interest.

# MODEL SETUP:
##############

wamresol=$(read_config grid)
wamnfre=$(read_config frequencies)
wambathy=$(read_config bathymetry)

# Directory where the output from this job will be saved
########################################################
mkdir -p ${RUN_DIR}
cd ${RUN_DIR}

WORK_DIR=${RUN_DIR}/workdir_create_bathymetry_$(date "+%Y.%m.%d-%H.%M.%S")
mkdir -p ${WORK_DIR}
cd ${WORK_DIR}


#. PRODUCE BATHYMETRY TO BE USED BY WAM
#######################################

source ${SCRIPTS_DIR}/ecwam_configure.sh
ecwam_configure $wamresol $wamnfre $wambathy

if [[ $wambathy = "aqua" ]] ; then
  rm -f ${RUN_DIR}/wam_topo
  touch ${RUN_DIR}/wam_topo
  echo "\n\t aqua bathymetry: empty file ${RUN_DIR}/wam_topo created"
  exit 0
elif [[ $wambathy = "ETOPO1" ]] ; then
  REFERENCE_LEVELS=${ETOPO1_reference_levels}
elif [[ $wambathy = "ETOPO2" ]] ; then
  REFERENCE_LEVELS=${ETOPO2_reference_levels}
else
  echo "\n\n\t UNSUPPORTED value for bathymetry: ${wambathy}\n"
  abort 1
fi
${SCRIPTS_DIR}/ecwam_retrieve.sh ${REFERENCE_LEVELS} ${DATA_DIR}/${REFERENCE_LEVELS}
ln -sf ${DATA_DIR}/${REFERENCE_LEVELS} ${WORK_DIR}/reference_levels

cat > for_md5 <<EOF
${ecwam_bathymetry_version}
EOF

if [ -f grid_description ]; then
  cat grid_description >> for_md5
else
  cat >> for_md5 <<EOF
${xdella}
${amosop}
${amonop}
${amowep}
${amoeap}
${iper}
${irgg}
${fr1}
${wamnfre}
${ifre1}
EOF
fi

cat reference_levels >> for_md5
md5=$(cat for_md5 | md5sum | awk '{print $1}')
WAM_TOPO=v${ecwam_bathymetry_version}/bathymetry_${wamresol}_nfre${wamnfre}_${wambathy}_${md5}

${SCRIPTS_DIR}/ecwam_retrieve.sh data/bathymetry/${WAM_TOPO} ${DATA_DIR}/data/bathymetry/${WAM_TOPO} || {
  echo "Could not retrieve data/bathymetry/${WAM_TOPO}"
}

if [ -f ${DATA_DIR}/data/bathymetry/${WAM_TOPO} ]; then
  echo "\n\n\t File ${WAM_TOPO} has been found\n"
else
  echo "\n\n\t File ${DATA_DIR}/data/bathymetry/${WAM_TOPO} was not found\n"
  echo     "\t It needs to be computed from $wambathy data set\n"

  if [[ $wambathy = "ETOPO1" ]] ; then
    echo "\n\n\t Getting ETOPO1 data set\n"
    log ${SCRIPTS_DIR}/ecwam_retrieve.sh ${ETOPO1} ${DATA_DIR}/${ETOPO1}
    ln -sf ${DATA_DIR}/${ETOPO1} ETOPO1_Ice_g_int.xyz
    CREATE_WAM_BATHYMETRY_EXE=${CREATE_WAM_BATHYMETRY_ETOPO1}
  else
    echo "\n\n\t Getting ETOPO2 data set\n"
    log ${SCRIPTS_DIR}/ecwam_retrieve.sh ${ETOPO2} ${DATA_DIR}/${ETOPO2}
    ln -sf ${DATA_DIR}/${ETOPO2} etopo2_2006apr.dat
    CREATE_WAM_BATHYMETRY_EXE=${CREATE_WAM_BATHYMETRY_ETOPO2}
  fi

  assert_executable_is_available ${CREATE_WAM_BATHYMETRY_EXE} || abort 4

  #. input files to create_wam_bathymetry
  cat > input_to_wam_bathymetry <<EOF
#User input to create_wam_bathymetry
EOF
  if [[ $wambathy = "ETOPO1" ]] ; then
    cat >> input_to_wam_bathymetry <<EOF
#Grib input:
F
#Grib output:
F
EOF
  fi

  cat >> input_to_wam_bathymetry <<EOF
#Grid resolution (XDELLA):
${xdella}
#Domain boundaries S,N,W,E (AMOSOP AMONOP AMOWEP AMOEAP) :
${amosop} ${amonop} ${amowep} ${amoeap}
#PERIODIC DOMAIN (IPER=0 for NON-periodic domain, IPER=1 for periodic domain):
${iper}
#REGULAR LATLON GRID (IRGG=0) OR IRREGULAR GRID (IRGG=1):
${irgg}
#REFERENCE INITIAL FREQUENCY (FR1):
${fr1}
#TOTAL NUMBER OF FREQUENCIES (NFRE) AND INDEX OF FR1 (IFRE1)
${wamnfre} ${ifre1} 
#DO YOU WANT TO CREATE OUTPUT FILES TO PLOT WITH METVIEW AS GEOPOINTS
F
#OVER WHICH AREA (WEST-EAST-SOUTH-NORTH):
-6 36 30 46
#DO YOU ALSO WANT THE ORIGINAl BATHYMETRY ON THAT AREA AS A METVIEW FILE WITH GEOPOINTS
F
EOF

  trace_ls $(pwd)

  echo "\n\n\t Starting ${LAUNCH_SERIAL} $(which ${CREATE_WAM_BATHYMETRY_EXE}) with input: \n"
  log cat input_to_wam_bathymetry
  echo "\n+ ${LAUNCH_SERIAL} $(which ${CREATE_WAM_BATHYMETRY_EXE}) < input_to_wam_bathymetry\n"
  START=$(date +%s)
  log ${LAUNCH_SERIAL} ${CREATE_WAM_BATHYMETRY_EXE} < input_to_wam_bathymetry || {
    echo "${CREATE_WAM_BATHYMETRY_EXE} FAILED"
    cleanup
    exit 1
  }
  END=$(date +%s)
  DIFF=$(( $END - $START ))
  echo "\n\n\t Running ${LAUNCH_SERIAL} $(which ${CREATE_WAM_BATHYMETRY_EXE}) took $DIFF seconds\n"

  trace_ls $(pwd)

  cwamresol=${cwamresol:-$(python3 -c "print(\"{:05d}\".format(int(1000*${xdella})))")}

  if [[ ! -r wam_topo_${cwamresol} ]] ; then
    echo "\n\n\t File wam_topo_${cwamresol} does not exist\n\n"
    abort 9
  fi
  mkdir -p ${DATA_DIR}/data/bathymetry/$(dirname ${WAM_TOPO})
  mv wam_topo_${cwamresol} ${DATA_DIR}/data/bathymetry/${WAM_TOPO}

fi

(cd ${RUN_DIR} && ln -sf ${DATA_DIR}/data/bathymetry/${WAM_TOPO} wam_topo)
echo "\n\t Bathymetry is available in DATA_DIR with symlink in RUN_DIR:\n\n\t   ${RUN_DIR}/wam_topo -> ${DATA_DIR}/data/bathymetry/${WAM_TOPO}"

cleanup