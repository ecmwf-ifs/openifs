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
#===============================================================================
# Build & run script for IFS Single Column Model
#------------------------------------------------
# Modifications history:
#  Mar 2020  R. Forbes  Heavily modified for 46r1
#  Oct 2021  R.Forbes  Updated for 47r1/47r3
#  Oct 2023  R.Forbes  Updated for 48r1
#
# To run the SCM:
#
# Default with no arguments:
#  run.scm
#  will assume expt_name="ref", namelist="namelist", scm_in="scm_in.nc"
#  (all must be in the current scm-run directory)
#
# With options:
#  run.scm [-c case_name] [-d duration] [-L nlevs] [-n namelist] [-o outdir]
#          [-s scm_in] [-t timestep] [-v] [-x expt_name] [-X exec]
#   -c case_name name of the case study used for namelist and output directory
#   -d duration  duration of run in hours (default 72h)
#   -L nlevs     number of levels (default 137)
#   -n namelist  namelist file (default namelist.case_study if case_study defined,
#                 otherwise ./namelist in current directory
#   -o outdir    working/output subdirectory (default current directory)
#   -s scm_in    SCM input/forcing netcdf file (defaults ./scm_in.nc in current directory)
#   -t timestep  timestep in seconds (default 450s)
#   -v           verbose turns on set -x (default off)
#   -x expt_name shortname to identify experiment
#   -X exec      path and filename of the executable
#
#===============================================================================

while getopts c:d:L:n:o:s:t:v:x:X: name
do
  case $name in
  c)
    CASE_NAME="$OPTARG"
    ;;
  d)
    DURATION="$OPTARG"
    ;;
  L)
    NFLEVG="$OPTARG"
    ;;
  n)
    NAMELIST="$OPTARG"
    ;;
  o)
    SCM_OUT="$OPTARG"
    ;;
  s)
    SCM_IN="$OPTARG"
    USE_FORCING_DIR=1
    ;;
  t)
    TIMESTEPOPT="$OPTARG"
    ;;
  v)
    VERBOSE=1
    ;;
  x)
    EXPT_NAME="$OPTARG"
    ;;
  X)
    EXEC="$OPTARG"
    ;;
  ?)
   echo "Usage: run.scm [-c case_name] [-d duration] [-L nlevs] [-n namelist] [-o outdir] "
   echo "               [-s scm_in] [-t timestep] [-v] [-x expt_name] [-X exec]"
   echo "-c case_name name of the case study used for namelist and output directory"
   echo "-d duration  duration of run in hours (default 72h)"
   echo "-L nlevs     number of levels (default 137)"
   echo "-n namelist  namelist file (default namelist.case_study if case_study defined,"
   echo "              otherwise ./namelist in current directory"
   echo "-o outdir    working/output subdirectory (default current directory)"
   echo "-s scm_in    SCM input/forcing netcdf file (defaults ./scm_in.nc in current directory)"
   echo "-t timestep  timestep in seconds (default 900s)"
   echo "-v           verbose turns on "set -ex" (default off)"
   echo "-x expt_name shortname to identify experiment"
   echo "-X exec      path and filename of the executable"
   exit 2;;
  esac
done

set -e
[ $VERBOSE ] && set -ex

#module load scm

# Increase ulimit otherwise executable fails
ulimit -s unlimited

########################################################
# Define variables
########################################################
export EC_MEMINFO=0
export LGRIB_API=0

# Set scm run directory (current directory)
RUNDIR=$SCM_RUNDIR

# Set scm project directory
#cd ../
PROJDIR=$SCM_PROJDIR
#cd ../
VERSIONDIR=$SCM_VERSIONDIR

cd "$RUNDIR"

# Set scm case input
FORCINGDIR=${VERSIONDIR}/scm-cases
NAMELISTDIR=${VERSIONDIR}/scm-cases

# Directory for ancillary input data (radiation)
IFSDATADIR=$VERSIONDIR/scm-ancillary
IFSDATALIST=`ls $IFSDATADIR`

# Define defaults for arguments
WORKDIR=${WORKDIR:-$RUNDIR}
NFLEVG=${NFLEVG:-137}       # model level: 19,31,40,50,60,91,137

# Check for subcase
CASE1_NAME=`echo ${CASE_NAME} | cut -d':' -f1`
CASE2_NAME=`echo ${CASE_NAME} | cut -d':' -f2`

# Set input SCM forcing
SCM_IN=${SCM_IN:-$RUNDIR/scm_in.nc}
if [ "${USE_FORCING_DIR}" ] ; then
  SCM_IN=${FORCINGDIR}/${CASE1_NAME}/${SCM_IN}
fi

TIMESTEP=${TIMESTEPOPT:-900}
DURATION=${DURATION:-'h72'}
EXPT_NAME=${EXPT_NAME:-ref}
      [ $VERBOSE ] && echo 0 $NAMELIST ${CASE2_NAME}

if [ ! ${NAMELIST} ] ; then
  if [ ${CASE2_NAME} ] ; then
    NAMELIST=$RUNDIR/namelist.${CASE2_NAME}
      [ $VERBOSE ] && echo 1 $NAMELIST
    if [ ! -f "${NAMELIST}" ] ; then
      NAMELIST=${NAMELISTDIR}/${CASE1_NAME}/namelist.${CASE2_NAME}
      [ $VERBOSE ] && echo 2 $NAMELIST
    fi
  else
    NAMELIST=$RUNDIR/namelist
  fi
      [ $VERBOSE ] && echo 3 $NAMELIST ${CASE2_NAME}
fi

SCMOUT_NAME=
if [ ${CASE_NAME} ] ; then SCMOUT_NAME=${SCMOUT_NAME}_${CASE2_NAME} ; fi
if [ ${EXPT_NAME} ] ; then SCMOUT_NAME=${SCMOUT_NAME}_${EXPT_NAME} ; fi
if [ ${TIMESTEPOPT} ]  ; then SCMOUT_NAME=${SCMOUT_NAME}_${TIMESTEP}s ; fi
WORKDIR=${RUNDIR}/scmout${SCMOUT_NAME}
export PLOT_NAME=${SCMOUT_NAME}

MASTER=$EXEC


echo
echo "*************************************************************************"
echo "SCM VERSION DIRECTORY" ${VERSIONDIR}
echo "SCM PROJECT DIRECTORY" ${PROJDIR}
echo "SCM RUN DIRECTORY    " ${RUNDIR}
echo "SCM WORKING DIRECTORY" ${WORKDIR}
echo "SCM EXPT NAME        " ${EXPT_NAME}
echo "USING EXECUTABLE     " ${MASTER}
echo "USING NAMELIST       " ${NAMELIST}
echo "USING IFSDATA        " ${IFSDATADIR}
echo "USING SCM INPUT      " ${SCM_IN}
echo "USING TIMESTEP       " ${TIMESTEP}
echo "NUMBER OF LEVELS     " ${NFLEVG}
echo "OUTPUT IN            " ${WORKDIR}
echo "*************************************************************************"
echo

###########################################
# Make working directory if it doesnt exist
###########################################
mkdir -p $WORKDIR

########################################
# Set up input files
########################################

# Copy IFS data files to working directory
cd ${IFSDATADIR}
cp ${NAMELIST} ${WORKDIR}/
NAMELISTFILE=`basename ${NAMELIST}`
cp ${IFSDATALIST} ${WORKDIR}/

# Change to working directory
cd $WORKDIR
chmod u+rw ${IFSDATALIST}

# Link forcing data
ln -sf $SCM_IN scm_in.nc

# Link vtable to appropriate vertical level info file
# /home/rd/rdx/data/ifs/vtable_L$NFLEVG
rm -f vtable
ln -s $VERSIONDIR/scm-ancillary/vtable_L$NFLEVG vtable
# Add vtable to end of namelist
cat vtable >> ${NAMELISTFILE}

# Change number of vertical levels in namelist file
# and timestep and duration of run
rm -f fort.4

# macOS (Darwin) based on FreeBSD, sed expects arg after -i option
SYS=$(uname -s)
if [ "${SYS}" = "Darwin" ]; then
   I_OPT='-i.old'
else
   I_OPT='-i'
fi

sed ${I_OPT} "s/TSTEP.*/TSTEP=${TIMESTEP},/" ${NAMELISTFILE}
sed ${I_OPT} "s/CSTOP.*/CSTOP='${DURATION}',/" ${NAMELISTFILE}
sed ${I_OPT} "s/NFLEVG.*/NFLEVG=${NFLEVG},/" ${NAMELISTFILE}
sed ${I_OPT} "s/NCEXTR.*/NCEXTR=${NFLEVG},/" ${NAMELISTFILE}
mv ${NAMELISTFILE} fort.4
cat ${VERSIONDIR}/scm-ancillary/namelist_extra >> fort.4

########################################
# Run SCM executable
########################################
export DR_HOOK=1

$MASTER

# Remove IFS data files from working directory
rm -f $IFSDATALIST

# Change fort. output files to meaningful filenames
mv fort.4 ${NAMELISTFILE}_input

exit
