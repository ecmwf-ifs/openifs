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

ECWAM_CONTEXT=model
source ${SCRIPTS_DIR}/ecwam_runtime.sh
source ${SCRIPTS_DIR}/ecwam_parse_commandline.sh
source ${SCRIPTS_DIR}/ecwam_helper_functions.sh

assert_executable_is_available ${MODEL} || abort 4

begofrn=$(read_config begin               --format="%Y%m%d%H%M%S")
endofrn=$(read_config end                 --format="%Y%m%d%H%M%S")
begoffo=$(read_config forcing.at[1].begin --format="%Y%m%d%H%M%S" --default=${endofrn} )

find_preproc_files
find_preset_files ${begofrn}

lgribin=F
lgribout=F
lnocdin=T

output_fields=$(read_config output.fields.name[:] --default="")
output_statistics=$(read_config output.statistics.name[:] --default="")

output_parameters_index=""
[[ -z "${output_fields}" ]] || output_parameters_index=$(${SCRIPTS_DIR}/ecwam_parameters.py --request=index ${output_fields})

OUTPUT_FLAGS=""
if [ $(read_config output.fields.format --default=grib) = grib ] ; then
  for index in ${output_parameters_index}; do
    [[ -z "${OUTPUT_FLAGS}" ]] || OUTPUT_FLAGS+=$'\n  '
    OUTPUT_FLAGS+="GFLAG(${index})=T,"
  done
else
  for index in ${output_parameters_index}; do
    [[ -z "${OUTPUT_FLAGS}" ]] || OUTPUT_FLAGS+=$'\n  '
    OUTPUT_FLAGS+="FFLAG(${index})=T,"
  done
fi

output_parameters_index=""
[[ -z "${output_statistics}" ]] || output_parameters_index=$(${SCRIPTS_DIR}/ecwam_parameters.py --request=index ${output_statistics})
for index in ${output_parameters_index}; do
    [[ -z "${OUTPUT_FLAGS}" ]] || OUTPUT_FLAGS+=$'\n  '
    OUTPUT_FLAGS+="NFLAG(${index})=T,"
done

NAWI=""
for i in $(seq 1 $(read_config forcings.at.size)); do
  [[ -z "${NAWI}" ]] || NAWI+=$'\n'
  forcing_end=$(read_config forcings.at[${i}-1].end --format=%Y%m%d%H%M%S)
  forcing_step=$(read_config forcings.at[${i}-1].timestep --format=seconds)
  NAWI+="&NAWI IDWI=${forcing_step}, IDWO=${forcing_step}, CLWOUT=\"${forcing_end}\" /"
done

NAOS=""
for restart_time in $(read_config output.restart.at[:].time --format=%Y%m%d%H%M%S --default=" "); do
  [[ -z "${NAOS}" ]] || NAOS+=$'\n'
  NAOS+="&NAOS CLSOUT=\"${restart_time}\" /"
done

# MODEL SETUP:
##############
wamnang=$(read_config directions)
wamnfre=$(read_config frequencies)

nproma=$(read_config nproma --default=24)

# read timesteps
phys_tstp=$(read_config physics.timestep --format=seconds --default=900)
adv_base_tstp=$(read_config advection.timestep --format=seconds --default=900)
adv_fast_tstp=$(read_config advection.fast_waves.timestep --format=seconds --default=$adv_base_tstp)
ifrelfmax=$(read_config advection.fast_waves.max_frequency --default=0)

# verify timesteps
if [ $(( $adv_base_tstp%$adv_fast_tstp )) -ne 0 ] ; then
   echo "ERROR: Base advection timestep should be a multiple of fast-wave advection timestep"
   exit 4
fi
if [ $(( $phys_tstp%$adv_base_tstp )) -ne 0 ] ; then
   echo "ERROR: Physics timestep should be a multiple of base advection timestep"
   exit 4
fi

if [[ $(read_config forcings.sea_ice --default=True) == "True" ]] ; then
  licerun=T
else
  licerun=F
fi

ppfreq=$(read_config output.fields.at[0].timestep --format=hours --default=02:00)

forcings_file=$(read_config forcings.file)

function cleanup() {
  shopt -s nullglob
  LOG_DIR=${RUN_DIR}/logs/model
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
  [[ -f ${WORK_DIR}/gstats         ]] && mv ${WORK_DIR}/gstats         ${LOG_DIR}/
  [[ -f ${WORK_DIR}/statistics.log ]] && mv ${WORK_DIR}/statistics.log ${LOG_DIR}/

  # Delete symbolic links
  find ${WORK_DIR} -maxdepth 1 -type l -delete

  mkdir -p ${RUN_DIR}/output/
  for outfile in ${WORK_DIR}/MPP*; do
    mv ${outfile} ${RUN_DIR}/output/
  done

  mkdir -p ${RUN_DIR}/restart/
  for outfile in ${WORK_DIR}/LAW* ${WORK_DIR}/BLS*; do
    mv ${outfile} ${RUN_DIR}/restart/
  done

  mv ${WORK_DIR}/wam_namelist ${RUN_DIR}/

  unset -f echo
  echo
  echo "Log     files have been copied to    ${LOG_DIR}"
  echo "Output  files have been copied to    ${RUN_DIR}/output/"
  echo "Restart files have been copied to    ${RUN_DIR}/restart/"

  cd ${RUN_DIR}
  rm -rf ${WORK_DIR}
  rm workdir # symlink

  shopt -u nullglob
}

workdir model
ln -sf ../config.yml config.yml
echo "*******************************************************************************"
log ${ECWAM_PROJECT_NAME} --info
echo "*******************************************************************************"

# General grid information and model constants
##############################################

for preproc_file in wam_grid_tables wam_subgrid_0 wam_subgrid_1 wam_subgrid_2; do
  if [[ -r ${RUN_DIR}/${preproc_file} ]] ; then
    ln -s  ${RUN_DIR}/${preproc_file} ${preproc_file}
  else
    echo "\n\n\t\tPREPROC file ${preproc_file} was not found in ${RUN_DIR}\n\n"
    trace_ls ${RUN_DIR}
    exit 1
  fi
done

# Initial stress and spectrum
#############################

for preset_file in LAW${begofrn}_000000000000 BLS${begofrn}_000000000000; do
  if [[ -r ${RUN_DIR}/restart/${preset_file} ]] ; then
    ln -sf ${RUN_DIR}/restart/${preset_file} ${preset_file}
  else
    echo "\n\n\tPRESET FILE ${preset_file} was not found in ${RUN_DIR}/restart"
    trace_ls ${RUN_DIR}/restart
    exit 1
  fi
done

# Forcing
#########

log ${SCRIPTS_DIR}/ecwam_retrieve.sh ${forcings_file} ${DATA_DIR}/${forcings_file} || {
  echo "ERROR: Could not retrieve forcings ${forcings_file}"
  exit 4
}

ln -s ${DATA_DIR}/${forcings_file} sfcwindin

# NAMELIST
##########

cat > wam_namelist << EOF
&NALINE
  NANG                  = ${wamnang},
  NFRE                  = 36,
  NFRE_RED              = ${wamnfre},
  CLHEADER              = " WAVE MODEL ",
  CBPLTDT               = "${begofrn}",
  CEPLTDT               = "${endofrn}",
  CDATEF                = "${begoffo}",
  DELPRO_LF             = ${adv_fast_tstp},
  IFRELFMAX             = ${ifrelfmax},
  IDELPRO               = ${adv_base_tstp},
  IDELT                 = ${phys_tstp},
  IDELINT               = ${ppfreq},
  IREST                 = 1,
  LFDBIOOUT             = F,
  LFDB                  = F,
  IPHYS                 = 1,
  ISHALLO               = 0,
  ISNONLIN              = 0,
  LBIWBK                = T,
  LLCAPCHNK             = T,
  LLGCBZ0               = F,
  LLNORMAGAM            = F,
  IPROPAGS              = 2,
  LSUBGRID              = F,
  IREFRA                = 0,
  LICERUN               = ${licerun},
  LMASKICE              = T,
  LWAMRSETCI            = T,
  NGRIB_VERSION         = 2,
  LL_GRID_SIMPLE_MATRIX = F,
  YCLASS                = "rd",
  YEXPVER               = "wave",
  ISTREAM               = 1045,
  CPATH                 = "${WORK_DIR}",
  LGRIBIN               = ${lgribin},
  LGRIBOUT              = ${lgribout},
  LNOCDIN               = ${lnocdin},
  NPROMA_WAM            = ${nproma},
  LL1D                  = F,
  LFRSTFLD              = T,
  IDELRES               = 0,  ! regular output for restart spectra, ignored if NAOS sections exists
  LRSTPARALR            = F,
  LRSTPARALW            = F,
  LSECONDORDER          = F,
  LWVFLX_SNL            = F,
  LLNORMWAMOUT          = T,
  LLNORMWAMOUT_GLOBAL   = T,
  CNORMWAMOUT_FILE      = "statistics.log",
  ${OUTPUT_FLAGS}
/
${NAWI}
${NAOS}
EOF

trace_ls $(pwd)

echo "*******************************************************************************"
echo "NAMELIST INPUT: wam_namelist"
echo "*******************************************************************************"
log cat wam_namelist
echo

precision=$(${ECWAM_PROJECT_NAME} --precision)
if ${dryrun}; then
  echo "*******************************************************************************"
  echo "WAVE MODEL DRYRUN"
  echo "*******************************************************************************"
  echo "#1 Go to workdir\n"
  echo "  cd ${RUN_DIR}/workdir"
  echo
  echo "#2 Execute wave model\n"
  echo "  $(abs_path ${MODEL})"
  echo 
  echo "#3 Validate results\n"
  echo "  ${SCRIPTS_DIR}/ecwam_validation.py config.yml statistics.log --section=validation.${precision}_precision"
  echo
  echo "#4 Clean up workdir when ready"
  echo "  ${RUN_DIR}/workdir/cleanup.sh"
  exit 0
fi

echo "*******************************************************************************"
echo "WAVE MODEL START"
echo "\n+ ${LAUNCH_PARALLEL} $(which ${MODEL})\n"
echo "*******************************************************************************"


START=$(date +%s)
log ${LAUNCH_PARALLEL} $(which ${MODEL}) || {
  sleep 1
  echo
  echo "*******************************************************************************"
  echo "WAVE MODEL FAILED EXECUTION"
  echo "*******************************************************************************"
  trace_ls $(pwd)
  cleanup
  exit 14
}

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "\n\n\t Running ${LAUNCH_PARALLEL} $(which ${MODEL}) took $DIFF seconds\n"

trace_ls $(pwd)

if [ "$(read_config validation.${precision}_precision --default=NOTFOUND)" != NOTFOUND ]; then
  echo "\n\n"
  echo "*******************************************************************************"
  echo "Validation"
  echo "*******************************************************************************"
  mkdir -p ${RUN_DIR}/logs/model
  cp statistics.log ${RUN_DIR}/logs/model/statistics.log
  echo " + ${SCRIPTS_DIR}/ecwam_validation.py ${RUN_DIR}/config.yml ${RUN_DIR}/logs/model/statistics.log --section=validation.${precision}_precision"
  log ${SCRIPTS_DIR}/ecwam_validation.py ${RUN_DIR}/config.yml ${RUN_DIR}/logs/model/statistics.log --section=validation.${precision}_precision || {
    cleanup
    exit 1
  }
fi

cleanup

touch ${RUN_DIR}/model.completed
