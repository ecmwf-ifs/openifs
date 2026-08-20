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
set -e
set +v

SCRIPTS_DIR="$( cd $( dirname "${BASH_SOURCE[0]}" ) && pwd -P )"
if [[ $(basename ${SCRIPTS_DIR}) == "bin" ]]; then
  SCRIPTS_DIR="$( cd $( dirname "${SCRIPTS_DIR}/$(readlink "${BASH_SOURCE[0]}")" ) && pwd -P )"
fi

ECWAM_CONTEXT=preset
source ${SCRIPTS_DIR}/ecwam_runtime.sh
source ${SCRIPTS_DIR}/ecwam_parse_commandline.sh
source ${SCRIPTS_DIR}/ecwam_helper_functions.sh

function cleanup() {
  shopt -s nullglob
  LOG_DIR=${RUN_DIR}/logs/preset
  mkdir -p ${LOG_DIR}
  rm -rf ${LOG_DIR}/*
  mv ${WORK_DIR}/stdout.log ${LOG_DIR}/
  mv ${WORK_DIR}/stderr.log ${LOG_DIR}/
  [[ -f ${WORK_DIR}/meminfo.txt    ]] && mv ${WORK_DIR}/meminfo.txt    ${LOG_DIR}/
  [[ -f ${WORK_DIR}/gstats.xml     ]] && mv ${WORK_DIR}/gstats.xml     ${LOG_DIR}/
  [[ -f ${WORK_DIR}/wam.log        ]] && mv ${WORK_DIR}/wam.log        ${LOG_DIR}/
  [[ -f ${WORK_DIR}/wam.log.1      ]] && mv ${WORK_DIR}/wam.log.*      ${LOG_DIR}/
  [[ -f ${WORK_DIR}/drhook.prof.1  ]] && mv ${WORK_DIR}/drhook*        ${LOG_DIR}/
  [[ -f ${WORK_DIR}/gstats.xml     ]] && mv ${WORK_DIR}/gstats.xml     ${LOG_DIR}/
  [[ -f ${WORK_DIR}/gstats.1       ]] && mv ${WORK_DIR}/gstats*        ${LOG_DIR}/
  [[ -f ${WORK_DIR}/statistics.log ]] && mv ${WORK_DIR}/statistics.log ${LOG_DIR}/

  # Delete symbolic links
  find ${WORK_DIR} -maxdepth 1 -type l -delete

  mkdir -p ${RUN_DIR}/restart/
  for outfile in ${WORK_DIR}/LAW* ${WORK_DIR}/BLS* ${WORK_DIR}/SGS* ; do
    mv ${outfile} ${RUN_DIR}/restart/
  done

  unset -f echo
  echo
  echo "Log     files have been copied to    ${LOG_DIR}"
  echo "Restart files have been copied to    ${RUN_DIR}/restart/"

  cd ${RUN_DIR}
  rm -rf ${WORK_DIR}
  rm workdir

  shopt -u nullglob
}


#.####################################################################################
#.
#. run_preset is used to generate the FIRST start files (cold start) for run_wamodel
#. To continue a simulation you will need to keep the restart files from that run
#. and put them at the right location in a manner similar to what is done with the
#. start files produced by this job (see run_wamodel).
#.
# The initial spectra are obtained from the Jonswap parametrisation.
# The input namelist NALINE (see below) determines how
# the parameters for the Jonswap parametrisation are determined. When
# IOPTI is 1 or 2 then the parameters are specified using a fetch law and
# the local wind provided as input as well. When IOPTI is 0 then the
# parameters are as specified in the namelist (see ALFA, FM, GAMMA, SA, SB,
# THETA). The fetch used is also specified as input (FETCH). The difference
# between IOPTI=1 and IOPTI=2 is no longer enforced since there is a 
# requirement that the minimum wind speed is 1m/s.
# For IOPTI=1 or 2 you will need a wind field (see DATA=
#.

#. Type of input/output files !!!!
##################################
#.
#. The wave spectra and fields needed to specify the initial wind and surface stress
#  can either be in grib format or in binary forms.
#  For this reason, one has to specify what format will be used to create
#  the initial data with preset and one has also to make sure that the same
#  format is specified in the namelist for the input fields for a wave model run
#  (see run_wamodel). Usually the same format is also used for the output of the wave model
#  run as some of the output is generally the input the the next run.
#  THE DEFAULT USED HERE IS PURE BINARY.

#  So in summary the following setting should specified in the input namelist
#  for PRESET                  |       for WAMODEL 
#.                             |
#.                      using binary format
#. LGRIBOUT=F                  |   LGRIBIN=F and LGRIBOUT=F
#.                             |
#.                      using GRIB format
#  LGRIBOUT=T                  |   LGRIBIN=T  and LGRIBOUT=T
#.                             |
#.
#. The grib format is commonly used at ECMWF, however, it will require tools
#. to archive and retrieve grib data (mars at ECMWF). Therefore the example
#. given here is for binary re-start files.
#.

#.  All output will be directed to directory specified by $WORK_DIR
#.  Make sure to clean it of the outputs if you rerun

#.  In the case of binary file outputs,
#.  the restart files (stress and spectra) in
#.  LAWYYYYMMDDHHmmss_000000000000 and BLSYYYYMMDDHHmmss_000000000000 in directory specified by $RUN_DIR

#.  The GRIB spectra will be saved in SGSYYYYMMDDHHmmss_000000000000
#.  where YYYYMMDDHHmmss is the starting date


#. ############
#. USER INPUTS:
#. ###########~
#. The example provided here creates initial conditions based on the ECMWF analysis winds
#. from date given by begofrn=
#. The wave model will then be run with analysis winds until begoffo= , updated every 6 hrs
#. and then with hourly forecast winds until endofrn=
#. Adapt accordingly !!!
#. (if no forecast is run, then set begoffo=endofrn)

#. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#. The necessary wind fields, as provided in wind_filename=wind_fields_${begofrn}_${begoffo}_${endofrn}
#. has to be copied to location specified by DATA=
#. If you do not use input files for the winds in grib format, you will need to adapt the WAM code
#. (see readwind.F). 
#. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# format  YYYYMMDDHHmmss (year-month-day-hour-minutes-seconds)
# begofrn=20180815000000
# begoffo=20180815120000
# endofrn=20180816000000

begofrn=$(read_config begin                --format="%Y%m%d%H%M%S")
endofrn=$(read_config end                  --format="%Y%m%d%H%M%S")
begoffo=$(read_config forcings.at[1].begin --format="%Y%m%d%H%M%S" --default=${endofrn} )

# MODEL SETUP:
##############
wamnang=$(read_config directions)
wamnfre=$(read_config frequencies)

forcings_file=$(read_config forcings.file)

opti=1
fetch=50000.0
fmax=0.2000000
llunstr=F
lgribout=F

assert_executable_is_available ${PRESET} || abort 4

# Directory where the output from this job will be saved
########################################################

workdir preset
ln -sf ../config.yml config.yml

#. Verify that all input files produced are available
#.###################################################

find_preproc_files

ln -sf ${RUN_DIR}/wam_grid_tables ${WORK_DIR}/wam_grid_tables

echo "+ ${SCRIPTS_DIR}/ecwam_retrieve.sh ${forcings_file} ${DATA_DIR}/${forcings_file}"

log ${SCRIPTS_DIR}/ecwam_retrieve.sh ${forcings_file} ${DATA_DIR}/${forcings_file} || {
  echo "ERROR: Could not retrieve forcings ${forcings_file}"
  exit 4
}

ln -s ${DATA_DIR}/${forcings_file} sfcwindin

# ${SCRIPTS_DIR}/generate_mars_request_oper_analysis_forecast.py --grid=O48 --date=${begofrn} --steps=12 --target=ecwam_forcings.grib | mars
# mv ecwam_forcings.grib ${RUN_DIR}/ecwam_forcings.grib
# ln -sf /perm/nawd/ecwam-concept-bundle/data/ecwam_forcings.grib sfcwindin

#. Run preset to generate the initial wave data
#. #############################################

#. Specify the input namelist.
#. The description of the input namelist is described in program preset.
cat > PREINFO <<EOF
&NALINE
 HEADER    = " WAVE MODEL INITIALISATION "
 CPATH     = "${WORK_DIR}"
 CDATEA    = "${begofrn}"
 IOPTI     = ${opti}
 ITEST     = 0
 ITESTB    = 0
 ALFA      = 1.800000E-02
 FM        = ${fmax}
 GAMMA     = 3.000000
 SA        = 7.000000E-02
 SB        = 9.000000E-02
 THETA     = 0.0
 FETCH     = ${fetch}
 LLUNSTR   = ${llunstr}
 ! IDELWI    =  21600
 CLTUNIT   = "S"
 LGRIBOUT  =  ${lgribout}
 MARSTYPE  = "an"
 YCLASS    = "rd"
 YEXPVER   = "wave"
 NANG      = ${wamnang}
 NFRE      = 36
 NFRE_RED  = ${wamnfre}
/
EOF


if ${dryrun}; then
  echo "*******************************************************************************"
  echo "WAVE MODEL PRESET DRYRUN"
  echo "*******************************************************************************"
  echo "#1 Go to workdir and set environment\n"
  echo "  cd ${RUN_DIR}/workdir"
  echo
  echo "#2 Execute preset model\n"
  echo "  $(abs_path ${PRESET}) < PREINFO"
  echo 
  echo "#3 Mark as complete for downstream scripts"
  echo "  touch ${RUN_DIR}/preset.completed"
  echo
  echo "#4 Clean up (delete) workdir"
  echo "   ${RUN_DIR}/workdir/cleanup.sh"
  exit 0
fi

echo "*******************************************************************************"
echo "PRESET START"
echo "\n+ ${LAUNCH_SERIAL} $(abs_path ${PRESET}) < PREINFO\n"
echo "*******************************************************************************"
echo "PREINFO:"
echo "*******************************************************************************"
log cat PREINFO
echo "*******************************************************************************"

START=$(date +%s)
log ${LAUNCH_SERIAL} ${PRESET} < PREINFO || {
  sleep 1
  echo
  echo "*******************************************************************************"
  echo "PRESET FAILED EXECUTION"
  echo "*******************************************************************************"
  trace_ls $(pwd)
  cleanup
  exit 14
}
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "\n\n\t Running ${LAUNCH_SERIAL} $(which ${PRESET}) took $DIFF seconds\n"

trace_ls ${WORK_DIR}

cleanup

touch ${RUN_DIR}/preset.completed

echo "==============================================================================="
echo "ecwam_run_preset end"
echo "==============================================================================="
