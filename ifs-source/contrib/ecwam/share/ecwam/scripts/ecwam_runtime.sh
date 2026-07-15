# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

SCRIPTS_DIR="$( cd $( dirname "${BASH_SOURCE[0]}" ) && pwd -P )"

# Version of topography result, to be increased when results for same run would be different.
# e.g. when a bug is fixed.
export ecwam_bathymetry_version=1

ECWAM_CACHE_PATH_DEFAULT=${HOME}/cache/ecwam
[[ ${HPCPERM} ]] && ECWAM_CACHE_PATH_DEFAULT=${HPCPERM}/cache/ecwam

export ECWAM_CACHE_PATH=${ECWAM_CACHE_PATH:-${ECWAM_CACHE_PATH_DEFAULT}}
export ecwam_ROOT=${ecwam_ROOT:-${SCRIPTS_DIR}/../../..}
export RUN_DIR=${ECWAM_RUN_DIR:-$(pwd)}
export DATA_DIR=${ECWAM_DATA_DIR:-${ECWAM_CACHE_PATH}}
export DOWNLOAD_DIR=${DOWNLOAD_DIR:-$(pwd)}

# Names of executables
ECWAM_PROJECT_NAME=$(basename $(cd ${SCRIPTS_DIR}/.. && pwd -P) )
export CREATE_WAM_BATHYMETRY_ETOPO1=${ECWAM_PROJECT_NAME}-create_wam_bathymetry_ETOPO1
export CREATE_WAM_BATHYMETRY_ETOPO2=${ECWAM_PROJECT_NAME}-create_wam_bathymetry
export PREPROC=${PREPROC:-${ECWAM_PROJECT_NAME}-preproc}
export PRESET=${PRESET:-${ECWAM_PROJECT_NAME}-preset}
export MODEL=${CHIEF:-${ECWAM_PROJECT_NAME}-chief}

export ETOPO1=data/bathymetry/ETOPO1/ETOPO1_Ice_g_int.xyz
export ETOPO2=data/bathymetry/ETOPO2/etopo2_2006apr.dat
export ETOPO1_reference_levels=data/bathymetry/ETOPO1/reference_levels.v2
export ETOPO2_reference_levels=data/bathymetry/ETOPO2/reference_levels.v1

# Executable runtime environment
export LAUNCH=${LAUNCH:-}
if [[ ${LAUNCH} == "" && $(uname) == "Darwin" && -x "$(command -v mpirun)" ]]; then
  # See README: Avoids SIGILL in debug builds during MPI_INIT
  export LAUNCH="ecwam-launch"
fi

export LAUNCH_PARALLEL=${LAUNCH_PARALLEL:-${LAUNCH}}
export LAUNCH_SERIAL=${LAUNCH_SERIAL:-${LAUNCH}}

export MBX_SIZE=${MBX_SIZE:-150000000}
export MPL_MBX_SIZE=${MPL_MBX_SIZE:-${MBX_SIZE}}

[[ ! (-d ${ecwam_ROOT}/bin) ]] || export PATH=${ecwam_ROOT}/bin:${ecwam_ROOT}/share/${ECWAM_PROJECT_NAME}/scripts:$PATH
