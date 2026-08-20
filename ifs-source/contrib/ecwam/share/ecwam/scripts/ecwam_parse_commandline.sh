# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

function usage() {
  echo "USAGE:"
  echo "    $(basename ${0}) --help"
  echo "    $(basename ${0}) --completed [--run-dir RUN_DIR]"
  echo "    $(basename ${0}) [--run-dir RUN_DIR] [--config CONFIG] [--launch LAUNCH_CMD] [--cache CACHE_PATH] [--halt]"
if [[ ${ECWAM_CONTEXT:-unset} == "model" ]]; then
  echo "                     [--np NTASKS] [--nt NTHREADS]"
fi
}

function help() {
  usage
  echo
  echo "Optional arguments"
  echo "    -h, --help                 Show this help message and exit"
  echo "    --completed                Instead of running, return 'true'/'false' whether execution is completed."
  echo "                               This checks for a file ${ECWAM_CONTEXT}.completed in the run-dir"
  echo "                               In order to reset completion status to 'false', simply remove that file."
  echo "    -c, --config CONFIG        YAML configuration file. If not provided,"
  echo "                               use ${RUN_DIR}/config.yml"
  echo "    -r, --run-dir RUN_DIR      Run directory. Overrides \${ECWAM_RUN_DIR}, see below."
  echo "    -l, --launch LAUNCH        Used as prefix to launch execution, e.g. \"mpirun -np <NTASKS>\", or \"ddt\""
  echo "                               This overrides defaults set by 'LAUNCH' environment variable"
  echo "    --cache CACHE_PATH         Path where downloaded and computed data will be stored"
  echo "                               Default: ${ECWAM_CACHE_PATH}"
  echo "    --halt                     Flag to stop before launching the binary execution, e.g. to allow for manually"
  echo "                               launching, or checking of generated inputs."
if [[ ${ECWAM_CONTEXT:-unset} == "model" ]]; then
  echo "    -np NTASKS                 Number of MPI tasks (ignored when LAUNCH is set, either via command line or environment)"
  echo "    -nt NTHREADS               Number of OpenMP threads per MPI task. This exports OMP_NUM_THREADS=<NTRHEADS>"
fi
  echo ""
  echo "    "
  echo
  echo "Advanced options via environment variables"
  echo
if [[ "${ECWAM_CONTEXT:-unset}" == "model" ]] ; then
  echo "    DR_HOOK=1                   # Enables drhook signal handling backtraces"
  echo "    DR_HOOK_OPT=PROF            # Enables drhook profiling, also requires DR_HOOK=1"
  echo "    GSTATS=1                    # Enables gstats timer output"
fi
  echo "    LAUNCH=<COMMAND>            # Used as prefix to launch execution, e.g. \"mpirun -np <NTASKS>\", or \"ddt\""
  echo "    OMP_NUM_THREADS=<NTHREADS>  # OpenMP threads"
  echo "    ECWAM_RUN_DIR               # Default run directory. If not set, current directory is assumed."
  echo "    ECWAM_CACHE_PATH            # Path where downloaded and computed data will be stored."
  echo "    ECWAM_DATA_PATH             # ':'-separated list of search locations for downloaded and computed data."
  echo "                                # Note that \$ECWAM_CACHE_PATH is foremost also a search location, taking precedence."
}

config=""
dryrun=false
do_check_completed=false
run_dir=""
cache=""
nthreads=0
ntasks=0
launch_cmd=""

while test $# -gt 0; do

    # Split --option=value in $opt="--option" and $val="value"

    opt=""
    val=""
    without_equal_sign=true
    case "$1" in
    --*=*)
      opt=`echo "$1" | sed 's/=.*//'`
      val=`echo "$1" | sed 's/--[_a-zA-Z0-9-]*=//'`
      without_equal_sign=false
      ;;
    --*)
      opt=$1
      val=$2
      ;;
    -*)
      opt=$1
      val=$2
      ;;
    *)
      break
      ;;
    esac

    # echo "debug_opt $opt $val $without_equal_sign"

    # Parse options
    case "$opt" in
      --help|-h)
        help; exit 0
        ;;
      --completed)
        do_check_completed=true
        ;;
      --halt)
        dryrun=true
        ;;
      --config|-c)
        if [[ ${val::1} == "-" || ${val} == "" ]] ; then
         echo "Argument for --config|-c CONFIG is missing"
         echo; usage; exit 1
        fi
        config=$val
        [[ $without_equal_sign == true ]] && shift
        ;;
      --run-dir|-r)
        if [[ ${val::1} == "-" || ${val} == "" ]] ; then
         echo "Argument for --run-dir|-r RUN_DIR is missing"
         echo; usage; exit 1
        fi
        run_dir=${val}
        [[ $without_equal_sign == true ]] && shift
        ;;
      --launch)
        if [[ ${val::1} == "-" || ${val} == "" ]] ; then
         echo "Argument for --launch|-l LAUNCH_CMD is missing"
         echo; usage; exit 1
        fi
        launch_cmd=${val}
        [[ $without_equal_sign == true ]] && shift
        ;;
      --cache)
        if [[ ${val::1} == "-" || ${val} == "" ]] ; then
         echo "Argument for --cache CACHE_PATH is missing"
         echo; usage; exit 1
        fi
        cache=${val}
        [[ $without_equal_sign == true ]] && shift
        ;;
      -np)
        if [[ ${val::1} == "-" || ${val} == "" ]] ; then
         echo "Argument for --np NTASKS is missing"
         echo; usage; exit 1
        fi
        ntasks=${val}
        [[ $without_equal_sign == true ]] && shift
        ;;
      -nt)
        if [[ ${val::1} == "-" || ${val} == "" ]] ; then
         echo "Argument for --nt NTHREADS is missing"
         echo; usage; exit 1
        fi
        nthreads=${val}
        [[ $without_equal_sign == true ]] && shift
        ;;
      *)
        echo "Unknown option: $opt"; echo; usage; exit 1
        ;;
    esac
    [[ $# -gt 0 ]] && shift
done

if [[ ${run_dir} != "" ]]; then
  if [[ -d ${run_dir} ]]; then
    if [[ ${do_check_completed} == true ]]; then
      echo "false";
      exit 0
    fi
  fi
  RUN_DIR=${run_dir}
fi

if [[ ${do_check_completed} == true ]]; then
  if [[ -f ${RUN_DIR}/${ECWAM_CONTEXT}.completed ]]; then
    echo "true"
  else
    echo "false"
  fi
  exit 0
fi

if [[ ${launch_cmd} != "" ]]; then
  if [[ ${ECWAM_CONTEXT:-unset} == "model" ]]; then
    LAUNCH_PARALLEL=${launch_cmd}
  else
    LAUNCH_SERIAL=${launch_cmd}
  fi
fi

if [[ ${nthreads} != 0 ]]; then
  export OMP_NUM_THREADS=${nthreads}
fi
if [[ ${ECWAM_CONTEXT:-unset} == "model" ]]; then
  if [[ ${ntasks} != 0 ]]; then
    if [[ "$LAUNCH_PARALLEL" == "" ]]; then
      export LAUNCH_PARALLEL=ecwam-launch
    fi
    if [[ "$LAUNCH_PARALLEL" =~ .*"ecwam-launch".* ]]; then
      export ECWAM_LAUNCH_NTASKS=${ntasks};
    else
      if [[ ${ntasks} != 1 ]]; then
        echo "WARNING: ignoring argument '-np ${ntasks}'"
      fi
    fi
  fi
fi

# Create RUN_DIR if not exists, and make absolute path
mkdir -p ${RUN_DIR}
RUN_DIR="$( cd ${RUN_DIR} && pwd -P )"

if [[ ${ECWAM_CONTEXT:-unset} != "unset" ]]; then
echo "==============================================================================="
echo "${ECWAM_PROJECT_NAME}-run-${ECWAM_CONTEXT}"
echo "==============================================================================="
fi

if [[ ${cache} != "" ]]; then
  ECWAM_CACHE_PATH=${cache}
fi

if [[ ${config} == "" ]]; then
  if [[ -f ${RUN_DIR}/config.yml ]]; then
    config=${RUN_DIR}/config.yml
  else
    config=${PWD}/config.yml
  fi
fi


if [[ ! -r ${config} ]]; then
  echo "Configuration file [${config}] was not found"
  echo
  usage
  exit 1
fi

if [ -z ${DO_ONCE+x} ]; then
  echo "*******************************************************************************"
  echo "RUN_DIR:  ${RUN_DIR}"
  echo "DATA_DIR: ${DATA_DIR}"
  echo "*******************************************************************************"
  echo "CONFIGURATION: ${config}"
  echo "*******************************************************************************"
  cat ${config}
  echo "*******************************************************************************"
  export DO_ONCE=true

  # Store symbolic link of configuration in RUN_DIR
  #################################################
  if [ -f ${RUN_DIR}/config.yml ]; then
    cmp --silent -- "${RUN_DIR}/config.yml" "${config}" || {
      echo "ERROR: Content of ${config} is different from content of already present ${RUN_DIR}/config.yml"
      echo "       - To use existing file, specify --config=${RUN_DIR}/config.yml"
      echo "       - To continue with new configuration, try again after deleting file ${RUN_DIR}/config.yml"
      exit 1
    }
  else
    config_dirname=$( cd $( dirname "${config}" ) && pwd -P )
    config_basename=$(basename ${config})
    [ -d ${RUN_DIR} ] || {
      echo "+ mkdir -p ${RUN_DIR}"
      mkdir -p ${RUN_DIR}
    }
    echo "+ ln -sf ${config_dirname}/${config_basename} ${RUN_DIR}/config.yml"
    ln -sf ${config_dirname}/${config_basename} ${RUN_DIR}/config.yml
  fi

fi
