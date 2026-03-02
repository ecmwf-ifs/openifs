# (C) Copyright 2011- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

#
#
#   oifs-config.edit_me.sh
#
#
#   This script sets the environment for OpenIFS 48r1
#
#
#   Read this script using the command:
#
#   source ./oifs-config.edit_me.sh
#
#
#--- set machine specific settings -----------------------------

export OIFS_HOST="local"
export OIFS_PLATFORM="local"

#--- set principal OIFS variables ------------------------------

export OIFS_CYCLE=48r1
export OIFS_CLIMATE="climate.v020"

#---Base code assumes openifs and openifs-expt are installed
#---in $HOME. Either these can be changed by the user------------
export OIFS_HOME="${HOME}/openifs"

#---Central location for non-openifs bundle source code, which
#---is required for the build and run. This only needs to be set
#---if not using the standard download of source to OIFS_HOME.
export OIFS_CENTRAL_SRC="${HOME}/openifs-bundle-src"

#---It is recommended that the openifs-expt and oifs_data dir
#---exist in a location designed for permanent storage-----------
export OIFS_EXPT="${HOME}/openifs-expt"
export OIFS_DATA_DIR="${OIFS_HOME}/openifs-data"

#---Set the path for the arch directory. Depending on system,i.e.,
#---all libs are installed on the sytem, this is not required,
#---so set to an empty string OIFS_ARCH=""
export OIFS_ARCH="./arch/${OIFS_HOST}/${OIFS_PLATFORM}"

#---Set the path for the directory that contains bin and share,
#---both of which are produced by the OpenIFS build. As standard,
#---this will be build, but can be install if -i option is used
#---in the build process.
export OIFS_BLD_PARENT="${OIFS_HOME}/build"

#---Path to the executable for 3d global model. This is the
#---default path for the exe, produced by openifs-test.sh.
#---SP means single precision. To run double precision change
#---SP to DP
export OIFS_EXEC="${OIFS_BLD_PARENT}/bin/ifsMASTER.SP"

#---Default assumed paths, only change if you know what you are doing
#---Path to the build script openifs-test.sh
export OIFS_TEST="${OIFS_HOME}/scripts/build_test"
#---Path to log for openifs-test.sh  script
export OIFS_LOGFILE="${OIFS_HOME}/openifs-test.log"
#---Path to dir containing scripts to run OpenIFS experiment
export OIFS_RUN_SCRIPT="${OIFS_HOME}/scripts/exp_3d"
#---Path to OpenIFS docker scripts and yaml config for docker
export OIFS_DOCKER="${OIFS_HOME}/scripts/docker"

alias oenv="env -0 | sort -z | tr '\0' '\n' | grep -a OIFS_"

echo -e "\nOpenIFS environment variables are:"
echo "------------------------------------------------------"
env -0 | sort -z | tr '\0' '\n' | grep -a OIFS_
echo

#---Path to the executable for the SCM. This is the
#---default path for the exe, produced by openifs-test.sh.
#---DP means double precision. To run single precision change
#---DP to SP
export SCM_EXEC="${OIFS_BLD_PARENT}/bin/MASTER_scm.SP"

#---Default assumed paths, only change if you know what you are doing
export SCM_TEST="${OIFS_HOME}/scripts/scm"
export SCM_VERSIONDIR="${OIFS_EXPT}/scm_openifs/48r1"
export SCM_PROJDIR="${SCM_VERSIONDIR}/scm-projects"
export SCM_RUNDIR="${SCM_PROJDIR}/ref48r1"
export SCM_LOGFILE="${SCM_RUNDIR}/scm_run_log.txt"

alias scm_env="env -0 | sort -z | tr '\0' '\n' | grep -a SCM_"

echo -e "\nSCM environment variables are:"
echo "------------------------------------------------------"
env -0 | sort -z | tr '\0' '\n' | grep -a SCM_
echo
