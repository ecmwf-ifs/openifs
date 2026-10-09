<p align="center">
<img width="300" height="90" alt="OpenIFS logo" src="docs/img/openifs-logo.png" />
</p>


# ECMWF OpenIFS

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

This repository contains code and scripts needed to build and run OpenIFS and the OpenIFS Single-Column Model.

## Contact

Contact information for OpenIFS Support is available on the OpenIFS home page: https://openifs.ecmwf.int/wiki. Support is given on a best-effort basis by the developers.

In addition to https://openifs.ecmwf.int/wiki, the [OpenIFS User Forums](https://forum.ecmwf.int/) are available to post support questions. These are monitored by the OpenIFS support team as well as members of the OpenIFS user community.

## Licence

[Apache License 2.0](LICENSE). In applying this licence, ECMWF does not waive the privileges and immunities granted to it by virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.

## Contributing

Contributions to OpenIFS are welcome. In order to do so, please create a pull request with your contribution and sign the contributor licence agreement (CLA).

## Supported Platforms

* Linux

Other UNIX-like operating systems, e.g. macOS, may also work out of the box, as long as the correct dependencies are installed.

## Prerequisites

The minimum software packages required to run OpenIFS on Linux (and UNIX-like operating systems) are the following:

* git
* cmake
* openmpi
* python3 python3-ruamel.yaml python3-yaml python3-venv
* libomp-dev
* libboost-dev libboost-date-time-dev libboost-filesystem-dev libboost-serialization-dev libboost-program-options-dev
* netcdf-bin libnetcdf-dev libnetcdff-dev
* libatlas-base-dev
* liblapack-dev
* libeigen3-dev
* bison
* flex

> Note: OpenIFS, as with the IFS, is constantly tested with a wide range of compilers, e.g. GNU/gcc, Intel and Cray. Even with this testing, we cannot and do not guarantee that all release branches will be compatible with all compiler versions.

## Installing and Building OpenIFS

> Note: It is recommended to build and install OpenIFS only from **supported release** branches rather than from unsupported branches or the main branch.

### Clone OpenIFS

OpenIFS is available directly from this repository and it can be extracted by either cloning or downloading the package:

* Extract just the release branch using a shallow clone that targets a specific release branch, e.g.
  * `git clone --depth 1 --branch openifs-lts/CY48R1.1 --single-branch https://github.com/ecmwf-ifs/openifs.git openifs-48r1.1`
* Extract a tagged release (shallow clone):
  * `git clone --depth 1 --branch openifs-48r1.1_20260825 --single-branch https://github.com/ecmwf-ifs/openifs.git openifs-48r1.1_20260825`

> Note: cloning a tag will result in a detached HEAD. If you plan to make commits, create a branch at that tag after cloning:

```bash
cd TAG
git switch -c my-branch-at-TAG
```

where TAG is the tag in the repository.

### Building OpenIFS

In this section, the build and test process is defined assuming the prerequisites exist and OpenIFS has been extracted.

Section [Docker install of OpenIFS](#docker-install-of-openifs) describes how to automate this process using a Docker container.

#### Set up the platform configuration file

The OpenIFS model requires a number of global environment variables to be set for both installation and runs. These environment variables are defined and set in the `oifs-config.edit_me.sh` file, which can be found at the top level of your extracted OpenIFS package.

The most important environment variable in `oifs-config.edit_me.sh` is `OIFS_HOME`, which is required by both model build and run scripts. For a description of other variables, please refer to [OpenIFS-env-vars](docs/oifs_env_vars.md).

Once edited, the platform configuration file is loaded using the following command:

```bash
source /path/to/file/location/oifs-config.edit_me.sh
```

For example, if you extracted OpenIFS into `$HOME/openifs`, the platform file would be loaded using

```bash
source $HOME/openifs/oifs-config.edit_me.sh
```

#### OpenIFS build

The build and initial test of OpenIFS is controlled by the script `$OIFS_HOME/scripts/build_test/openifs-test.sh`.

Once the platform configuration file, `oifs-config.edit_me.sh`, has been sourced (see previous section), OpenIFS can be built using the following commands:

```bash
cd $OIFS_HOME
$OIFS_TEST/openifs-test.sh -cb
```

where:

`$OIFS_TEST` is defined in the platform configuration file (`oifs-config.edit_me.sh`) as `$OIFS_HOME/scripts/build_test`.

* `-c` creates a `source` directory in `$OIFS_HOME`, which is used to collect all the sources defined in the `bundle.yml`, in preparation for the build
* `-b` builds `source`. This step creates the directory `build` in `$OIFS_HOME`, which is used to build and store the OpenIFS and SCM executables.

For more details about `openifs-test.sh` and the available options, please refer to [OpenIFS-build-options](docs/oifs_build_options.md).

### Test OpenIFS build

Once executables are successfully built, they can be tested using the following command:

```bash
cd $OIFS_HOME
$OIFS_TEST/openifs-test.sh -t
```

where

* `-t` invokes the testing simulations, which are coarse resolution T21 tests, comprising 21 3-D NWP tests with and without chemistry and 1 SCM test (based on TWP-ICE).

> Note:
> OpenIFS build and test can be run together using `$OIFS_TEST/openifs-test.sh -cbt`.
> The defaults in `oifs-config.edit_me.sh` set the host and site to `local`, which assumes that all the dependencies are installed and available locally.
> If running on an HPC, it is probably necessary to use an arch file. For example, if running on ECMWF HPC either set `OIFS_HOST` and `OIFS_PLATFORM` in `oifs-config.edit_me.sh` as follows:

```bash
export OIFS_HOST="ecmwf"
export OIFS_PLATFORM="hpc2020"
```

> or use `$OIFS_TEST/openifs-test.sh -cbt --arch=./arch/ecmwf/hpc2020/gnu` or `$OIFS_TEST/openifs-test.sh -cbt --arch=./arch/ecmwf/hpc2020/intel`, depending on whether you want to use the Intel or GNU compiler.

If everything has worked correctly with the build of OpenIFS, then all tests should have passed and `openifs-test.sh` returns the following

```bash
[INFO]: Good news - ctest has passed
        openifs is ready for experiment and SCM testing
----------------------------------------------------------------
END ifstest on OpenIFS build
```

A 100% pass with `$OIFS_TEST/openifs-test.sh -cbt` shows that the low-resolution (T21) ifs-test cases can run to completion on the chosen system. These tests do not check bit comparability with known good output. If this is a requirement, e.g., if a user makes a code change and needs to test whether the code has led to unexpected behaviour in the code, then please refer to [OpenIFS-test-options](docs/oifs_test_options.md).

### Docker install of OpenIFS

The previous section, [Installing and Building OpenIFS](#installing-and-building-openifs), describes the pre-requisites and build process for OpenIFS on a generic Linux-based system.

[create-oifs-docker.py](scripts/bootstrap/docker/create-oifs-docker.py) and associated scripts and configuration automate the process described in section [Installing and Building OpenIFS](#installing-and-building-openifs), by creating a Docker container, installing OpenIFS and dependencies and then building OpenIFS and running the tests.

* Please go to [OpenIFS Docker Builder](scripts/bootstrap/docker/README.md) for details about the Docker install.

[create-oifs-docker.py](scripts/bootstrap/docker/create-oifs-docker.py) and the resulting Docker development have been tested on macOS but can be applied to other systems, as long as Docker is installed and the appropriate Python dependencies are available.

## Install the static input data files for OpenIFS

OpenIFS requires **static input data** to run forecast experiments. Many of these static data files **are specific to** the respective **model cycle**. Please see the **release notes** for the respective model cycle to find the corresponding **climate version**.

* This static data must be located in `$OIFS_DATA_DIR`.
* As a minimum, you will require the packages `rtables.tar.gz` and `ifsdata.tar.gz` and **at least one** of the packages for a specific horizontal grid resolution (e.g. `48r1_climate.v020_159.tar.gz` for T159 at model cycle 48r1 using climate version v020).

Create the required directory structure, then download and install the static files by following these steps:

```bash
# Replace this path with the actual installation path.
source "/path/to/installation/oifs-config.edit_me.sh"

# Step 1: Create the directory structure for static files.
mkdir -p \
  "$OIFS_DATA_DIR/ifsdata" \
  "$OIFS_DATA_DIR/rtables" \
  "$OIFS_DATA_DIR/$OIFS_CLIMATE"

# Step 2: Download and extract radiation table files.
cd "$OIFS_DATA_DIR/rtables"
wget "https://openifs.ecmwf.int/data/ifsdata/$OIFS_CYCLE/rtables/rtables.tar.gz"
tar -xvzf rtables.tar.gz

# Step 3: Download and extract resolution-independent data files.
cd "$OIFS_DATA_DIR/ifsdata"
wget "https://openifs.ecmwf.int/data/ifsdata/$OIFS_CYCLE/ifsdata/ifsdata.tar.gz"
tar -xvzf ifsdata.tar.gz

# Step 4: Download and extract resolution-dependent data files (here for T159).
cd "$OIFS_DATA_DIR/$OIFS_CLIMATE"
archive="${OIFS_CYCLE}_${OIFS_CLIMATE}_159.tar.gz"

wget "https://openifs.ecmwf.int/data/ifsdata/$OIFS_CYCLE/$OIFS_CLIMATE/$archive"
tar -xvzf "$archive"
```

You should repeat the final Step 4 above for all additional grid resolutions that you intend to use. Browse all available grid resolutions here: `https://openifs.ecmwf.int/data/ifsdata/${OIFS_CYCLE}/${OIFS_CLIMATE}/`

## Run a standard OpenIFS 3-D NWP experiment

The following is a step-by-step guide on how to run a global forecast experiment. A more detailed description of these steps can be found in [How to run global forecast experiments](docs/oifs_howto_run_experiments.md).

### Set up the experiment directory

An **example forecast experiment** has been prepared. You will need to download the **experiment data pack** matching your model cycle. We recommend selecting the N80 (T159) model grid.

* All experiments have a four-character alphanumeric identifier which we will refer to as `<experiment-id>` in the documentation below.

* A listing of the available experiments, including their experiment ID and download information, is shown in the table found in document [List of OpenIFS Experiments](docs/oifs_experiment_list.md). You should **download the experiment data pack** for your model cycle.

* Set variable `OIFS_EXPT` in `oifs-config.edit_me.sh` to point to a suitable location path for your model experiment and extract the downloaded experiment data pack to this location.

* Copy the model run scripts and experiment configuration file into the experiment directory.

To carry out the above steps, run the following commands in your terminal:

```bash
# Replace this path with the actual installation path:
source "/path/to/installation/oifs-config.edit_me.sh"

cd $OIFS_EXPT

# Replace the download URL in the next command with the link to your experiment pack,
# ensuring it ends with the correct experiment-id tarball:
wget <download-url>/<experiment-id>.tar.gz

# Extract the experiment data pack and change into the experiment directory;
# replace <experiment-id> with the name of your experiment:
tar -xvzf <experiment-id>.tar.gz
cd $OIFS_EXPT/<experiment-id>/2016092500

# Copy required scripts and configuration files from your OpenIFS installation:
cp $OIFS_RUN_SCRIPT/exp-config.h .
cp $OIFS_RUN_SCRIPT/oifs-run .
cp $OIFS_RUN_SCRIPT/run-oifs.ecmwf-hpc2020.job .
```

Edit the experiment configuration file `exp-config.h` and
* replace the OIFS_EXPID value with your experiment's identifier
* replace the OIFS_RES value with your experiment's resolution (here: 159 for the N80 grid)
* replace the OIFS_GRIDTYPE value with your experiment's type (here: l for the N80 grid)
* adjust the number of MPI tasks and OpenMP threads to values that are suitable for your computing system.

**exp-config.h:**
```bash
#--- required variables for this experiment:
 
# this is specific for each experiment:
OIFS_EXPID="<experiment-id>"  #  your experiment ID
OIFS_RES="159"          #  the spectral grid resolution (here: T159)
OIFS_GRIDTYPE="l"       #  the grid type, either 'l' for linear reduced grid, 
                        #  or 'o' for cubic octahedral grid

# note: use of the batch job script will overwrite these values:
OIFS_NPROC=8            #  the number of MPI tasks
OIFS_NTHREAD=4          #  the number of OpenMP threads

# postprocessing is optional but recommended:
OIFS_PPROC=false        #  enable/disable postprocessing of model output after the model run
OUTPUT_ROOT=$(pwd)      #  folder where pproc output is created (only used if 
                        #  OIFS_PPROC=true). In this example an output folder is 
                        #  created in the experiment directory.

LFORCE=true             #  overwrite existing symbolic links in the experiment directory
LAUNCH=""               #  the platform specific run command for the MPI environment
                        #  (e.g. "mpirun", "srun", etc). Setting this variable overwrites 
                        #  any platform-specific default run commands defined in oifs-run
 
#--- optional variables that can be set for this experiment:
 
#OIFS_NAMELIST='my-fort.4'               # custom atmospheric model namelist file
#OIFS_EXEC="<custom-path>/ifsMASTER.DP"  # model executable to be used for this experiment
```

### Running the experiment

Now the model run can be started. Depending on the available hardware, the experiment can either be run interactively or as a batch job.

#### Running interactively

If your computing system is capable of running programs in a hybrid parallel configuration (MPI/OpenMP), then you can run the model interactively.

* To run the experiment interactively, execute the `oifs-run` script from the command line in your terminal.

```bash
# run interactively:
source </path/to/installation>/oifs-config.edit_me.sh
cd $OIFS_EXPT/<experiment-id>/2016092500
./oifs-run
```

#### Running a batch job

This method is the preferred way to run OpenIFS, as it is more efficient and it allows more flexibility in using the available hardware resources. 

* The job wrapper script `run-oifs.ecmwf-hpc2020.job` has been designed for the ECMWF hpc2020 HPC and might need adjusting for your local computing system.
* Edit this file and adjust the header lines for the batch scheduler as required.
* Edit the variable `PLATFORM_CFG` to point to your `oifs-config.edit_me.sh` file.

Once you have made the appropriate changes, the job can be submitted:

```bash
# run as slurm batch job:
source </path/to/installation>/oifs-config.edit_me.sh
cd $OIFS_EXPT/<experiment-id>/2016092500
sbatch ./run-oifs.ecmwf-hpc2020.job
```

>NOTE: The job wrapper script will read the `exp-config.h` file and use the selected values. The exceptions are `LAUNCH`, which is set to "srun" for batch jobs, and `OIFS_NPROC` & `OIFS_NTHREAD` for which values from the batch job headers are used. The job wrapper script modifies the `exp-config.h` file accordingly prior to calling the `oifs-run` script.


## Run a standard OpenIFS SCM case

Since the release of OpenIFS 48r1, the Single-Column Model (SCM) has been included in the source code and is built by default when OpenIFS is built. In this section we present an overview of how to set up and run the SCM.

### Setting up and building the SCM

As with all OpenIFS operations, the SCM depends on environment variables defined in `oifs-config.edit_me.sh`, e.g. for model cycle 48r1 see below:

```bash
#---Path to the executable for the SCM. This is the
#---default path for the exe, produced by openifs-test.sh.
#---SP means single precision. To run double precision change
#---SP to DP
export SCM_EXEC="${OIFS_BLD_PARENT}/bin/MASTER_scm.SP"

#---Default assumed paths, only change if you know what you are doing
export SCM_TEST="${OIFS_HOME}/scripts/scm"
export SCM_VERSIONDIR="${OIFS_EXPT}/scm_openifs/${OIFS_CYCLE}"
export SCM_PROJDIR="${SCM_VERSIONDIR}/scm-projects"
export SCM_RUNDIR="${SCM_PROJDIR}/ref${OIFS_CYCLE}"
export SCM_LOGFILE="${SCM_RUNDIR}/scm_run_log.txt"
```

SCM environment variables depend on the `OIFS_HOME` and `OIFS_EXPT`, which are also defined by sourcing `oifs-config.edit_me.sh`.

Before attempting to run the SCM, please follow the instructions in section [Set up the platform configuration file](#set-up-the-platform-configuration-file).

### SCM standard test-case package

The standard test-case package consists of three test cases, each representative of different cloudy regimes:

* DYCOMS - marine stratocumulus case
* BOMEX - trade-wind cumulus case
* TWPICE - a multi-day deep convective case

For ease of use with the standard OpenIFS environment variables, we recommend that the SCM test-case package is installed in `$OIFS_EXPT`. For example, in the template `oifs-config.edit_me.sh`, `$OIFS_EXPT=${HOME}/openifs-expt`. In this scenario, the directory `scm_openifs` needs to be in `$OIFS_EXPT` or `${HOME}/openifs-expt/`.

We show below the installation steps:

```bash
# Replace this path with the actual installation path:
source "/path/to/installation/oifs-config.edit_me.sh"

cd $OIFS_EXPT
wget https://openifs.ecmwf.int/data/scm/${OIFS_CYCLE}/scm_openifs_${OIFS_CYCLE}.tar.gz
tar -xvzf scm_openifs_${OIFS_CYCLE}.tar.gz
```

> Note: The extracted SCM package is small, ~45 MB, and the data produced by a standard individual SCM simulation is also small. However, if a user is planning to perform many simulations and store the data, which is often the case, the disk space usage can become large. If this is the plan, then a user may need to consider installing the SCM test-case package on a larger disk than `$HOME`.

### Run the SCM

Once the SCM test-case package installation has been completed, the SCM is run using the `callscm` script, which is a wrapper for the main `run.scm`. Both scripts can be found in `$SCM_TEST`, which is set in the `oifs-config.edit_me.sh` file to `${OIFS_HOME}/scripts/scm`.

`callscm` includes default settings, which are the three cases, with a 450 s timestep and an experiment name of ref-oifs-scm. To run with these settings, enter the following

```bash
cd $OIFS_HOME
$SCM_TEST/callscm
```

> Note: If running on the ECMWF HPC, the MPI environment needs to be loaded to avoid runtime MPI errors and the SCM failing when run with `callscm`. Use the following to load the environment

```bash
# If OpenIFS and SCM built with Intel compiler use
module load prgenv/intel
module load intel-mpi
# if OpenIFS and SCM built with GNU compiler use
module load prgenv/gnu
module load gcc/11.2.0
module load openmpi/4.1.1.1
```

`callscm` (with defaults, i.e. no arguments) will run the DYCOMS, BOMEX and TWPICE cases with the SCM and create an output directory in `$SCM_RUNDIR/scmout_DYCOMS_ref-oifs-scm_450s`, which contains the diagnostic output from the SCM. In addition, the file `scm_run_log.txt` will be created in `$SCM_RUNDIR`. This file contains the print output from the SCM, which is useful for checking all the sources and paths for a simulation.

#### `callscm` command-line options

Some of the `callscm` defaults can be changed through command-line options, e.g.

```
callscm -h -c <case_name or list of case_names> -t <timestep or list of timesteps>
        -x <expt_name>
where:
-h is help which returns basic usage options and exits
-c case_name or list of case_names (space delimited) of the case study
   used for namelist and output directory. Default list is
   "DYCOMS BOMEX TWPICE"
-t timestep or list of timesteps in seconds. The default is 450s. An
   example of a list is "1800 900 300"
-x expt_name shortname to identify experiment. Default is ref-oifs-scm
```

For example, if a user wanted to run the BOMEX case with timesteps of 1800 s and 900 s and an experiment name of "bomex_test", they would enter the following

```bash
$SCM_TEST/callscm -c BOMEX -t "1800 900" -x "bomex_test"
```

This command results in the following output directories: `$SCM_RUNDIR/scmout_BOMEX_bomex_test_1800s` and `$SCM_RUNDIR/scmout_BOMEX_bomex_test_900s`.
