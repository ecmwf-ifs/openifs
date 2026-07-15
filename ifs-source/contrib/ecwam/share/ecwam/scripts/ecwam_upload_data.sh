#!/usr/bin/env bash

# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

# Following https://confluence.ecmwf.int/display/CIT/Management+of+Nexus+repos

SCRIPTS_DIR="$( cd $( dirname "${BASH_SOURCE[0]}" ) && pwd -P )"

source ${SCRIPTS_DIR}/ecwam_runtime.sh
export NEXUS_SYSTEM=DOWNLOAD_INSTANCE # get.ecmwf.int
export NEXUS_REPOSITORY=ecwam
export ECMWF_USERNAME=${ECMWF_USERNAME:-$USER}
export ECMWF_PASSWORD=${ECMWF_PASSWORD:-$(read -sp 'ECMWF_PASSWORD: ' passvar; echo $passvar)}

# Install nexus-cli
command -v nexus-cli >/dev/null 2>&1 || {
  if [ ! -d $HOME/.venv_nexus ] ; then
    python3 -m venv $HOME/.venv_nexus
    source $HOME/.venv_nexus/bin/activate
    pip3 install nexus-toolkit -U -i https://get.ecmwf.int/repository/pypi-all/simple
  fi
  source $HOME/.venv_nexus/bin/activate
}


# Utlity functions

function check_uploaded {
  nexus-cli search --query=${1} >/dev/null 2>&1
}

function upload {
  yes | nexus-cli upload --local-path=${1} --remote-path=${location}
}

function compress () {
  case ${2} in
    *.gz)   echo "    + gzip  -9 -c ${1}  > ${2}";  time gzip  -9 -c $(readlink -e ${1})  > ${2} ;;
    *.bz2)  echo "    + bzip2 -9 -c ${1}  > ${2}";  time bzip2 -9 -c $(readlink -e ${1})  > ${2} ;;
    *.zst)  echo "    + zstd -19 -f ${1} -o ${2}";  time zstd -19 -f $(readlink -e ${1}) -o ${2} ;;
  esac
}

# Upload bathymetry files from ${ECWAM_CACHE_PATH}/data/bathymetry/v${ecwam_bathymetry_version}

location=data/bathymetry/v${ecwam_bathymetry_version}
cd ${ECWAM_CACHE_PATH}/${location}
COMPRESSION_DIR=${ECWAM_CACHE_PATH}/${location}
mkdir -p ${COMPRESSION_DIR}
shopt -s extglob
for f in bathymetry!(*.bz2|*.zst|*.gz) ; do
  echo "Processing ${ECWAM_CACHE_PATH}/${location}/$f"

  for ext in ".zst" ".gz" ".bz2"; do
    fc=$f$ext
    check_uploaded ${location}/${fc} && echo "  - ${location}/${fc} already uploaded" || {
      if [ ! -f ${COMPRESSION_DIR}/${fc} ]; then
        echo "  - Compressing ${f} -> ${COMPRESSION_DIR}/${fc}"
        compress ${f} ${COMPRESSION_DIR}/${fc}
      fi
      echo "  - Uploading ${fc} -> ${location}/${fc}"
      upload ${COMPRESSION_DIR}/${fc}
    }
  done
done

echo
echo "Finished uploading bathymetry"
echo "Compressed files can be deleted from ${COMPRESSION_DIR}"


# Upload forcing files from ${ECWAM_CACHE_PATH}/data/forcings

echo
location=data/forcings
cd ${ECWAM_CACHE_PATH}/${location}
COMPRESSION_DIR=${ECWAM_CACHE_PATH}/${location}
mkdir -p ${COMPRESSION_DIR}
shopt -s extglob
for f in *.grib ; do
  echo "Processing ${ECWAM_CACHE_PATH}/${location}/$f"

  # Upload grib file without additional compresion
  check_uploaded ${location}/${f} && echo "  - ${location}/${f} already uploaded" || {
      echo "  - Uploading ${f} -> ${location}/${f}"
      upload $PWD/${f}
  }

  # Upload grib file with additional compression
  for ext in ".zst"; do
    fc=$f$ext
    check_uploaded ${location}/${fc} && echo "  - ${location}/${fc} already uploaded" || {
      echo "need to upload ${location}/${fc}"
      if [ ! -f ${COMPRESSION_DIR}/${fc} ]; then
        echo "  - Compressing ${f} -> ${COMPRESSION_DIR}/${fc}"
        compress ${f} ${COMPRESSION_DIR}/${fc}
      fi
      echo "  - Uploading ${fc} -> ${location}/${fc}"
      upload ${COMPRESSION_DIR}/${fc}
    }
  done
done

echo
echo "Finished uploading fordings"
echo "Compressed files can be deleted from ${COMPRESSION_DIR}"
