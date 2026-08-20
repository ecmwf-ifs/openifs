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

ECWAM_CONTEXT=preproc
source ${SCRIPTS_DIR}/ecwam_runtime.sh
source ${SCRIPTS_DIR}/ecwam_parse_commandline.sh
source ${SCRIPTS_DIR}/ecwam_helper_functions.sh

function cleanup() {
  shopt -s nullglob
  LOG_DIR=${RUN_DIR}/logs/preproc
  mkdir -p ${LOG_DIR}
  rm -rf ${LOG_DIR}/*
  mv ${WORK_DIR}/stdout.log ${LOG_DIR}/
  mv ${WORK_DIR}/stderr.log ${LOG_DIR}/

  # Delete symbolic links
  find ${WORK_DIR} -maxdepth 1 -type l -delete

  mv ${WORK_DIR}/wam_grid_tables  ${RUN_DIR}/wam_grid_tables

  for ip in 0 1 2; do
    mv ${WORK_DIR}/wam_subgrid_${ip} ${RUN_DIR}/wam_subgrid_${ip}
  done
  
  # Does not seem to be used further
  # mv PARWAM ${RUN_DIR}/wam_parwam

  unset -f echo
  echo
  echo "Log   files have been copied to    ${LOG_DIR}"
  echo "Other files have been copied to    ${RUN_DIR}"

  shopt -u nullglob
}

#.############################################################################
#.
#. run_preproc generates unformatted input files (bathymetry, grid and tables)
#. for use by run_preset and run_wamodel 
#.
#. THE INPUT BATHYMETRY FILES :
#. bathymetry_global_1.0
#. OR
#. ETOPO2: etopo2_2006apr.dat.gz
#. supplied here
#. ETOPO1: ETOPO1_Ice_g_int.xyz.gz
#. you will need to get it from NOAA web site:
#. https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/xyz/
#.
#. !!! THEY MUST BE COPIED TO DIRECTORY SPECIFIED BY DATA= (see below) !!!
#.
#.############################################################################

# MODEL SETUP:
##############

wamresol=$(read_config grid)
wamnfre=$(read_config frequencies)
wambathy=$(read_config bathymetry)

assert_executable_is_available ${PREPROC} || abort 33

# Directory where the output from this job will be saved
########################################################
mkdir -p $RUN_DIR
cd $RUN_DIR

workdir preproc
ln -sf ../config.yml config.yml

source ${SCRIPTS_DIR}/ecwam_configure.sh
ecwam_configure $wamresol $wamnfre $wambathy

if [[ $wambathy == "aqua" ]]; then
  laqua=true
else
  laqua=false
fi

#. INPUT FOR PREPRPOC
#####################
# input namelists for preproc (used to specify domain and resolutions
#                              as well as manual bathymetric corrections)
#. Note that namelist NACORR is only working when llobstrct is false
#. otherwise correct bathymetry using the input file reference levels 
#. By default, the model will used the wave physics for Ardhuin et al. 2010
#. for wind input and dissipation. The previous physics can be used instead
#  by setting IPHYS=0 (see below)
cat > procin <<EOF
&NALINE
  CLINE=     " PREPROC INPUT "
  NFRE=      36,
  NFRE_RED=  ${wamnfre},
  FR1=       ${fr1},
  IFRE1=     ${ifre1},
  IRGG=      ${irgg},
  XDELLA=    ${xdella},
  XDELLO=    ${xdella},
  AMOSOP=    ${amosop},
  AMONOP=    ${amonop},
  AMOWEP=    ${amowep},
  AMOEAP=    ${amoeap},
  LAQUA=     ${laqua},
  LLOBSTRCT= ${llobstrct},
  LLUNSTR =  ${llunstr},
  IFORM=     1,
  ITEST=     0,
  ITESTB=    4,
  IBOUNC=    0,
  IBOUNF=    0, AMOSOC= 0.0, AMONOC= 0.0, AMOWEC= 0.0, AMOEAC= 0.0
  CLDOMAIN=  '${cldomain}'
/
EOF

if [[ ${laqua} = true ]]; then
  echo "Aqua planet (laqua=true), no need for bathymetry"
  touch wam_topo
else
# Make sure we have the correct bathymetry
    ${SCRIPTS_DIR}/ecwam_run_create_bathymetry.sh --run-dir=${RUN_DIR} || {
    echo "Could not run ecwam_run_create_bathymetry"
    exit 1
  }
  ln -sf ${RUN_DIR}/wam_topo wam_topo
fi

# run preproc:
##############
echo "\n\nls -ltrch $(pwd)\n"
ls -ltrch

if ${dryrun}; then
  echo "*******************************************************************************"
  echo "WAVE MODEL PREPROC DRYRUN"
  echo "*******************************************************************************"
  echo "#1 Go to workdir and set environment\n"
  echo "  cd ${RUN_DIR}/workdir"
  echo
  echo "#2 Execute preset model\n"
  echo "  $(abs_path ${PREPROC})"
  echo 
  echo "#3 Mark as complete for downstream scripts"
  echo "  touch ${RUN_DIR}/preproc.completed"
  echo
  echo "#4 Clean up (delete) workdir"
  echo "   ${RUN_DIR}/workdir/cleanup.sh"
  exit 0
fi

echo "\n\n\t Starting ${LAUNCH_SERIAL} $(which ${PREPROC}) with namelist input:\n\n"
log cat procin
echo
echo "\n+ ${LAUNCH_SERIAL} $(which ${PREPROC})\n"
START=$(date +%s)
log ${LAUNCH_SERIAL} $(which ${PREPROC}) || {
   echo "\n\n\t ${PREPROC} FAILED\n\n"
   sleep 1
   abort 14
  }
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "\n\n\t Running ${LAUNCH_SERIAL} $(which ${PREPROC}) took $DIFF seconds\n"

trace_ls $(pwd)

cleanup

touch ${RUN_DIR}/preproc.completed
