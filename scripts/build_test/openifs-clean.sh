#!/usr/bin/env bash
#
# (C) Copyright 2011- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

#  This script will remove the source, build, ecbundle and *.out or *.log 
#  from the openifs directory 
#
# Check that OIFS_HOME is defined and we are inside OIFS_HOME
set -eu
set -o pipefail

check_user_input() 
{
  # used in check_arguments
  while true; do
    printf "%s (y/n): " "$1"
    read -r choice
    case "$choice" in
      [Yy])
        return 0
        ;;
      [Nn])
        return 1
        ;;
      *)
      echo "Please answer 'y' or 'n'."
        ;;
    esac
  done
}

force=false
while getopts "f" OPTION ; do

  case "${OPTION}" in 
    f) force=true ;;
    *) usage_info "incorrect command line argument passed" ;;
  esac
done

echo "===================================================================================" 
echo "[WARNING]: This script removes source, build, ecbundle and *.out or *.log from OpenIFS home" 
if [[ "$force" == "false" ]]; then
  if check_user_input "[INPUT] Do you want to proceed with a full clean: Y/n"; then
    echo "[INFO]: User input is yes, so full clean requested - CONTINUE" 
  else 
    echo "[INFO]: User check is no, so cancel full clean - EXIT " 
    exit 1
  fi
else 
  echo "[INFO]: -f - No checks - CONTINUE"
fi
echo "==================================================================================="

if [ -n "${OIFS_HOME-}" ] ;  then 
  echo "[INFO]: OIFS_HOME path exists and is $OIFS_HOME"
  if [ -d "$OIFS_HOME" ]; then
    echo "[INFO]: cd $OIFS_HOME"
    cd "$OIFS_HOME" || { echo "[ERROR]: $OIFS_HOME does not exist - EXIT" ; exit 1; }
  fi
else
  echo "[ERROR]: OIFS_HOME variable is unset or empty." 
  echo "         Please edit and source oifs-config.edit_me.sh"
  echo "         to set OIFS_HOME to the OpenIFS directory." 
  echo "         EXITING"

  exit 1
fi

echo "[INFO]: BEGIN full clean - rm build, source, ecbundle and *.out or *.log" 
if ! rm -rf build source ecbundle *.log *.out ; then 
  echo "[ERROR]: Failed to remove Build and/or source - EXITING" 
  exit 1
else
  echo "[INFO]: Removed build and source successful" 
  echo "[INFO]: END full clean - rm build, source, ecbundle and *.out or *.log"   
fi 