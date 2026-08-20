#!/usr/bin/env bash

filename=ICMGGepc8INIUA

nbo3=$(grib_ls -w shortName=o3 ${filename} | tail -n 1 | cut -f 1 -d ' ')
if [[ "${nbo3}" -eq 0 ]] ; then
  grib_copy -w shortName=q ${filename}  ${filename}_q
  grib_set -d 1e-6 -s shortName=o3 ${filename}_q ${filename}_o3
  newname=$(echo ${filename} | sed -e "s/INIUA$/IMIUA/")
  cat ${filename}_o3 ${filename} > ${newname}
  rm -f ${filename}_q ${filename}_o3
fi

