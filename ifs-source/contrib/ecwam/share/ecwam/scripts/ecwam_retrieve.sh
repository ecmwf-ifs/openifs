#!/usr/bin/env bash

# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

set -ea
set -o pipefail

FILE_PATH=$1
FILE=$(basename ${FILE_PATH})
TARGET_FILE_PATH=${2:-$(pwd)/${FILE_PATH}}
if [[ "${TARGET_FILE_PATH}" == */ ]] ; then
  TARGET_FILE_PATH=${TARGET_FILE_PATH}/${FILE}
fi
TARGET_FILE=$(basename ${TARGET_FILE_PATH})
LOCATION=$(dirname ${TARGET_FILE_PATH})

SCRIPTS_DIR="$( cd $( dirname "${BASH_SOURCE[0]}" ) && pwd -P )"

source ${SCRIPTS_DIR}/ecwam_runtime.sh

################################################################################################
# 1. Return if FILE is already available

if [ -f "${TARGET_FILE_PATH}" ]; then
  echo "File is available in LOCATION: ${TARGET_FILE_PATH}"
  exit
fi

################################################################################################
# 2. Find FILE in search locations on disk

declare -a SEARCH_LOCATIONS

# Add search location in ECWAM_CACHE_PATH, a writeable directory
if [[ ! -z ${ECWAM_CACHE_PATH} ]] ; then
  SEARCH_LOCATIONS+=(${ECWAM_CACHE_PATH})
fi

# Add search locations from ':' separated list "ECWAM_CACHE_PATH"
if [[ ! -z ${ECWAM_DATA_PATH} ]] ; then
  for d in $(echo ${ECWAM_DATA_PATH} | tr ":" "\n"); do
    SEARCH_LOCATIONS+=(${d})
  done
fi

SEARCH_LOCATIONS+=(${ecwam_ROOT}/share/ecwam)
SEARCH_LOCATIONS+=(${PWD})

for SEARCH_LOCATION in ${SEARCH_LOCATIONS[*]}; do
  if [ -f "${SEARCH_LOCATION}/${FILE_PATH}" ]; then
    echo "File was found in: ${SEARCH_LOCATION}/${FILE_PATH}"
    mkdir -p ${LOCATION}
    ln -sf ${SEARCH_LOCATION}/${FILE_PATH} ${TARGET_FILE_PATH}
    echo "Symbolic link is available in LOCATION: ${TARGET_FILE_PATH}"
    exit
  fi
  if [ -f "${SEARCH_LOCATION}/${FILE}" ]; then
    echo "File was found in: ${SEARCH_LOCATION}/${FILE}"
    mkdir -p ${LOCATION}
    ln -sf ${SEARCH_LOCATION}/${FILE} ${TARGET_FILE_PATH}
    echo "Symbolic link is available in LOCATION: ${TARGET_FILE_PATH}"
    exit
  fi
done

################################################################################################
# 3. FILE is not found in any path hierarchy. Download from get.ecmwf.int

ECWAM_URL=https://get.ecmwf.int/repository/ecwam

function command_exists () {
    type "$1" &> /dev/null ;
}

function trace () {
  echo "+ $@"
  eval "$@"
}


function url_exists () {
  curl -o/dev/null -sfI "$1"
}

function decompress () {
  case ${1} in
    *.gz)  echo "+ gzip  -ckd ${1} > ${2}.decompress"; gzip  -ckd ${1} > ${2}.decompress ;;
    *.bz2) echo "+ bzip2 -ckd ${1} > ${2}.decompress"; bzip2 -ckd ${1} > ${2}.decompress ;;
    *.zst) echo "+ zstd  -ckd ${1} > ${2}.decompress"; zstd  -ckd ${1} > ${2}.decompress ;;
    *)     trace cp ${1} ${2}.decompress ;;
  esac
  trace mv ${2}.decompress ${2}
}

function url() {
  if   ( command_exists zstd  ) && ( url_exists ${ECWAM_URL}/${1}.zst ); then
    echo ${ECWAM_URL}/${1}.zst
  elif ( command_exists bzip2 ) && ( url_exists ${ECWAM_URL}/${1}.bz2 ); then
    echo ${ECWAM_URL}/${1}.bz2
  elif ( command_exists gzip  ) && ( url_exists ${ECWAM_URL}/${1}.gz  ); then
    echo ${ECWAM_URL}/${1}.gz
  elif ( url_exists ${ECWAM_URL}/${1} ); then
    echo ${ECWAM_URL}/${1}
  else
    echo "URL_NOT_FOUND"
  fi
}

if [ ! -f "${DOWNLOAD_DIR}/${FILE}" ]; then
  mkdir -p ${DOWNLOAD_DIR} && cd ${DOWNLOAD_DIR}

  ### From here ...
  # This is not ideal, just saves a download of compressed ETOPO1/ETOPO2 dataset on HPC.
  # Better would be that the uncompressed bathymetry is present in a directory part of ECWAM_DATA_PATH
  if [ -f ~rdx/data/wave/bathymetry_data/${FILE}.gz ]; then
    FILE_RDX_PATH=~rdx/data/wave/bathymetry_data/${FILE}.gz
    echo "File was found in ${FILE_RDX_PATH}"
    FILE_COMPRESSED=$(basename ${FILE_RDX_PATH})
    ln -sf ${FILE_RDX_PATH} ${FILE_COMPRESSED}
  ### ... until here
  else
    URL=$(url ${FILE_PATH})
    url_exists ${URL} || {
      echo "${FILE_PATH} or ${FILE} were not found in any search location."
      echo "A suitable download URL is also not available."
      echo "Tried: ${ECWAM_URL}/${FILE_PATH}"
      exit 1
    }
    FILE_COMPRESSED=$(basename ${URL})
  fi
  if [ ! -f ${FILE_COMPRESSED} ]; then
    echo "Downloading from URL: ${URL}"
    trace curl -L -C - "${URL}" --output ${FILE_COMPRESSED}.download
    trace mv ${FILE_COMPRESSED}.download ${FILE_COMPRESSED}
    echo "Download completed"
  fi
  if [ "${FILE}" != "${FILE_COMPRESSED}" ]; then
    echo "Decompressing ${FILE_COMPRESSED}"
    decompress ${FILE_COMPRESSED} ${FILE}
    echo "Decompression completed"
    rm ${FILE_COMPRESSED}
  fi
fi

cmp -s ${DOWNLOAD_DIR}/${FILE} ${TARGET_FILE_PATH} || {
  echo "Moving file to LOCATION"
  mkdir -p ${LOCATION}
  trace mv ${DOWNLOAD_DIR}/${FILE} ${TARGET_FILE_PATH}
}

echo "File is available in LOCATION: ${TARGET_FILE_PATH}"
