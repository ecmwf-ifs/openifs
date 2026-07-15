# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

# HELPER FUNCTIONS:
###################

SCRIPTS_DIR="$( cd $( dirname "${BASH_SOURCE[0]}" ) && pwd -P )"

function log () {
  # Send stdout/stderr both to file and terminal
  ("$@" | tee -a stdout.log) 3>&1 1>&2 2>&3 | tee -a stderr.log
}

function echo () {
  log builtin echo -e "$@"
}

function assert_executable_is_available {
  command -v $1 >/dev/null 2>&1 || {
    echo "\n\n\t \"$1\" was not found in PATH\n"
    echo "PATH: $PATH\n"
  }
  command -v $1 >/dev/null 2>&1
}

function trace_ls() {
  echo "\n\n+ ls -ltrch ${1}\n"
  log ls -ltrch ${1}
}

function read_config() {
  ${SCRIPTS_DIR}/ecwam_read_config.py ${RUN_DIR}/config.yml "$@"
}

function workdir() {
  cd ${RUN_DIR}
  WORK_DIR=${RUN_DIR}/workdir_${1}_$(date "+%Y.%m.%d-%H.%M.%S")
  mkdir -p ${WORK_DIR}
  [[ -e workdir ]] && rm workdir
  ln -sf ${WORK_DIR} workdir
  ls -l ${RUN_DIR}
  [[ -f ${RUN_DIR}/stdout.log ]] && mv ${RUN_DIR}/stdout.log ${WORK_DIR}/stdout.log
  [[ -f ${RUN_DIR}/stderr.log ]] && mv ${RUN_DIR}/stderr.log ${WORK_DIR}/stderr.log
  cd ${WORK_DIR}

  cat > cleanup.sh <<EOF
#!/usr/bin/env bash
RUN_DIR=${RUN_DIR}
WORK_DIR=${WORK_DIR}
$(declare -f cleanup)
cleanup
EOF
  chmod 755 cleanup.sh
}

function abs_path {
  builtin echo "$(cd $(dirname $(which ${1})) && pwd)/${1}"
}

function find_preproc_files() {
  SEARCH_LOCATION=${ECWAM_PREPROC_RUN_DIR:-${RUN_DIR}}
  found=true
  files=(wam_grid_tables wam_subgrid_0 wam_subgrid_1 wam_subgrid_2)
  for file in "${files[@]}"; do
    if [[ ! ( -r ${SEARCH_LOCATION}/${file} ) ]] ; then
      found=false
    fi
  done
  if $found ; then
    for file in "${files[@]}"; do
      if [[ ! ( -r ${RUN_DIR}/${file} ) ]] ; then
        ln -sf ${SEARCH_LOCATION}/${file} ${RUN_DIR}/${file}
      fi
    done
  else
    builtin echo
    builtin echo "ERROR: Grid files not found in ${SEARCH_LOCATION}:"
    for file in "${files[@]}"; do
      if [[ ! ( -r ${SEARCH_LOCATION}/${file} ) ]] ; then
        builtin echo "         - ${file}"
      fi
    done
    builtin echo
    builtin echo "       To generate, run:"
    builtin echo
    builtin echo "           $(abs_path ${ECWAM_PROJECT_NAME}-run-preproc)"
    exit 1
  fi
  builtin echo "preproc files have been found in ${SEARCH_LOCATION}"
}

function find_preset_files() {
  date=${1}
  SEARCH_LOCATION=${ECWAM_PRESET_RUN_DIR:-${RUN_DIR}}
  found=true
  files=(restart/LAW${date}_000000000000 restart/BLS${date}_000000000000)
  for file in "${files[@]}"; do
    if [[ ! ( -r ${SEARCH_LOCATION}/${file} ) ]] ; then
      found=false
    fi
  done
  if $found ; then
    mkdir -p ${RUN_DIR}/restart
    for file in "${files[@]}"; do
      if [[ ! ( -r ${RUN_DIR}/${file} ) ]] ; then
        echo "        ln -sf ${SEARCH_LOCATION}/${file} ${RUN_DIR}/${file}"
        ln -sf ${SEARCH_LOCATION}/${file} ${RUN_DIR}/${file}
      fi
    done
  else
    builtin echo
    builtin echo "ERROR: Initial condition files not found in ${SEARCH_LOCATION}:"
    for file in "${files[@]}"; do
      if [[ ! ( -r ${SEARCH_LOCATION}/${file} ) ]] ; then
        builtin echo "         - ${file}"
      fi
    done
      builtin echo
      builtin echo "       To generate, run:"
      builtin echo
      builtin echo "           $(abs_path ${ECWAM_PROJECT_NAME}-run-preset)"
    exit 1
  fi
  builtin echo "preset files have been found in ${SEARCH_LOCATION}"
}
