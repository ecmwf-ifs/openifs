
# ECWMF OpenIFS

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

This repository contains code and scripts need to build and run the OpenIFS and OpenIFS Single-Column Model.

## Contact

Contact information for OpenIFS Support is available on the OpenIFS home page: https://openifs.ecmwf.int/wiki. Support is given on a best-effort basis by the developers.

In addtion to https://openifs.ecmwf.int/wiki, the [OpenIFS User Forums](https://forum.ecmwf.int/) are available to post support questions. These are monitored by the OpenIFS support team as well as members of the OpenIFS user community.

## Licence

License: [Apache License 2.0](LICENSE) In applying this licence, ECMWF does not waive the privileges and immunities granted to it by virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.

## Contributing

Contributions to OpenIFS are welcome. In order to do so, please create a pull request with your contribution and sign the contributors license agreement (CLA).

## Supported Platforms

* Linux 

Other UNIX-like operating systems, e.g. Mac OS, may work too out of the box, as long as the correct dependencies are installed.

## Pre-requisites

The minimum software packages required to run OpenIFS on Linux (and UNIX-like operating systems) is the following

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

> Note: OpenIFS, as with the IFS, is constantly tested with a wide range of compilers, e.g. gnu/gcc, intel and cray. Even with this testing, we cannot and do not guarantee all release branches will be compatible with all compiler versions.

## Installing and Building OpenIFS

OpenIFS is available direct from this repository and it can be extracted by either cloning or downloading the package using 

* Extract the entire OpenIFS repository by either executing the following command in the directory where you want to extract OpenIFS 
  * `git clone https://github.com/ecmwf-ifs/openifs.git`
* Extract just the release branch using a shallow clone that targets a specific release branch, e.g.
  * `git clone --depth 1 --branch release/openifs-48r1 --single-branch https://github.com/ecmwf-ifs/openifs.git openifs-48r1`
* Extract a tagged release (shallow clone):
  * `git clone --depth 1 --branch TAG --single-branch https://github.com/ecmwf-ifs/openifs.git openifs-TAG`

> Note: cloning a tag will result in a detached HEAD. If you plan to make commits, create a branch at that tag after cloning:

```bash
cd openifs-TAG
git switch -c my-branch-at-TAG
```

### Building OpenIFS

#### Set up the platform configuration file

The OpenIFS model requires a number of Linux global environment variables to be set for both installation and runs. These environment variables are defined and set in the `oifs-config.edit_me.sh` file, which can be found in the top-level of your extracted OpenIFS package.

The most important environment variable in `oifs-config.edit_me.sh` is `OIFS_HOME`, which is required by both model build and run scripts. For description of other variables please refer to [OpenIFS-env-vars](docs/oifs_env_vars.md).

Once edited the platform configuration file is loaded using the following command: 

```bash
source /path/to/file/location/oifs-config.edit_me.sh
```

For example, if you extracted OpenIFS into `$HOME/openifs`, the platform file would be loaded using 

```bash
source $HOME/openifs/oifs-config.edit_me.sh
```

#### OpenIFS build

The build of OpenIFS and optional running of initial tests, which broadly test the build, is controlled by the build script `$OIFS_HOME/scripts/build_test/openifs-test.sh`. To run the build process and the tests use the following commands (assumes the platform configuration file has been sourced) :

```bash
cd $OIFS_HOME
$OIFS_TEST/openifs-test.sh -cb
```

where:

`$OIFS_TEST` is defined in the platform configuration file (`oifs-config.edit_me.sh`) as `$OIFS_HOME/scripts/build_test`.

* `-c`  creates `source` directory in `$OIFS_HOME`, which is used to collects all the sources defined in the `bundle.yml`, in preparation for the build
* `-b`  builds the source. This step creates the directory `build` in `$OIFS_HOME`, which is used to build and store the OpenIFS and SCM executables.
  
For more details about `openifs-test.sh` and the available options please refer to [OpenIFS-build-options](docs/oifs_build_options.md).

### Test OpenIFS build

Once executables are successfully built, they can be tested using the following command 

```bash
cd $OIFS_HOME
$OIFS_TEST/openifs-test.sh -t
```

where

* `-t` invokes the testing simulations, which are coarse resolutions t21 tests, comprising of 21 3-D NWP tests with and without chemistry and 1 SCM test (based on TWP-ICE).

> Note: OpenIFS build and test can be run together using `$OIFS_TEST/openifs-test.sh -cbt`

> Note: The defaults in `oifs-config.edit_me.sh` set the host and site to `local`, which assumes that all the dependencies are installed and defined. This works in the docker install. 
> If running on an HPC, it is probably necessary to use an arch file. For example, if running on ECMWF HPC either set `OIFS_HOST` and `OIFS_PLATFORM` in `oifs-config.edit_me.sh`, i.e.
  
```bash
export OIFS_HOST="ecmwf"
export OIFS_PLATFORM="hpc2020"
```  

> or use `$OIFS_TEST/openifs-test.sh -cbt --arch=./arch/ecmwf/hpc2020/gnu` or `$OIFS_TEST/openifs-test.sh -cbt --arch=./arch/ecmwf/hpc2020/intel`, depending if you want to use the Intel or gnu compiler.

If everything has worked correctly with the build of OpenIFS, then all tests should have passed and the `openifs-test.sh` returns the following

```bash
[INFO]: Good news - ctest has passed
        openifs is ready for experiment and SCM testing
----------------------------------------------------------------
END ifstest on OpenIFS build
```

> NOTE: 100% pass with `$OIFS_TEST/openifs-test.sh -cbt` shows that the low resolution (t21) ifs-test cases can run to completion on the chosen system. These tests do not check bit comparibility with known good output. If this is a requirement, e.g., if a user makes a code change and needs to test whether the code has led to unexpected behaviour in the code, then please refer to [OpenIFS-test-options](docs/oifs_test_options.md).

## Run a standard OpenIFS 3-D NWP experiment

### Set up the experiment directory

An **example forecast experiment** has been prepared for OpenIFS 48r1. The experiment ID is `ab2a` and you can download the tarball from here: https://openifs.ecmwf.int/data/experiments/48r1/2016-09-25_Karl/ab2a.tar.gz

* Set variable `OIFS_EXPT` in `oifs-config.edit_me.sh` to point to a suitable location for your model experiments.
* Extract the example forecast experiment `ab2a.tar.gz` to this location. We will refer to directory `$OIFS_EXPT/ab2a/2016092500/` as the **experiment directory**.
> NOTE: The `OIFS_EXPT` location for model experiments will often be in a different location from the model installation `OIFS_HOME`. In general, you will require more disk space for experiments, depending on the model grid resolution, the duration of the forecast experiment and the output frequency of model results.

```
cd $OIFS_EXPT
wget https://openifs.ecmwf.int/data/experiments/48r1/2016-09-25_Karl/ab2a.tar.gz
tar -xvzf ab2a.tar.gz
```
* Ensure the **namelist files** for the atmopsheric model (fort.4) and for the wave model (wam_namelist) are found in the experiment directory. Backup copies should remain in the `ecmwf` subfolder:

```
source oifs-config.edit_me.sh  # always do this first
cd $OIFS_EXPT/ab2a/2016092500
cp ./ecmwf/fort.4 .
cp ./ecmwf/wam_namelist .
```
* Copy the model run scripts and the experiment configuration file into the experiment directory.
    - `oifs-run`: This is a generic **run script** which executes the binary model program file.
    - `exp-config.h`: This is the **experiment configuration file** that determines settings for your experiment; this file is read by oifs-run.
    - `run-oifs.ecmwf-hpc2020.job`: This is a **wrapper script** that calls oifs-run and submits a non-interactive job to the batch scheduler on the ECMWF hpc2020; to use this script you may need to adjust it for your HPC system.

```
cd $OIFS_EXPT/ab2a/2016092500
cp $OIFS_HOME/scripts/exp_3d/oifs-run .
cp $OIFS_HOME/scripts/exp_3d/exp-config.h .
cp $OIFS_HOME/scripts/exp_3d/run-oifs.ecmwf-hpc2020.job .
```

### Customise the experiment parameters

Here we describe you can customise the way how OpenIFS will run the forecast experiment.

#### Customising the namelist:

* You can edit the atmospheric model namelist file `fort.4`. It contains Fortran namelists which control model settings and switches.
* An important switch is `CSTOP` in namelist `NAMRIP`. Set this to the desired duration of the forecast experiment.
* Experiment `ab2a` can be run for up to 144 hours (6 days) by setting `CSTOP='h144'`.

#### Customising the experiment configuration file: 

* You should customise the `exp-config.h` file which determines the settings used for this experiment.
* The oifs-run script will read the settings from this file.
* Alternatively, the settings can be passed to the oifs-run script via command line parameters, which takes precedence over the exp-config.h settings. 
> NOTE: The preference should be to **always set up an exp-config.h** for each experiment. If no exp-config.h file is found in the experiment directory (and no command line parameters are provided when calling oifs-run), then oifs-run will revert to its own default values which are not appropriate.

**exp-config.h:**
```
#--- required variables for this experiment:
 
# this is specific for each experiment:
OIFS_EXPID="ab2a"       # your experiment ID
OIFS_RES="255"          # the spectral grid resolution (here: T255)
OIFS_GRIDTYPE="l"       # the grid type, either 'l' for linear reduced grid, 

# use of the batch job script will overwrite these values:
OIFS_NPROC=8            # the number of MPI tasks
OIFS_NTHREAD=4          # the number of OpenMP threads

# postprocessing is optional but recommended:
OIFS_PPROC=true         # enable postprocessing of model output after the model run
OUTPUT_ROOT=$(pwd)      # folder where pproc output is created (only used if 
                        #   OIFS_PPROC=true). In this example an output folder is 
                        #   created in the experiment directory.

LFORCE=true             # overwrite existing symbolic links in the experiment directory
LAUNCH=""               # the platform specific run command for the MPI environment
                        #   (e.g. "mpirun", "srun", etc).
 
#--- optional variables that can be set for this experiment:
 
#OIFS_NAMELIST='my-fort.4'               # custom atmospheric model namelist file
#OIFS_EXEC="<custom-path>/ifsMASTER.DP"  # model executable to be used for this experiment
```

### Running the experiment

After all edits to the namelists (`fort.4`) and to the experiment configuration file (`exp-config.h`) have been completed the model run can be started.

Depending on the available hardware experiments can either be run interactively or as a batch job.

#### Running a batch job:

This method is the preferred way to run OpenIFS, as it is more efficient and it allows more flexibility in using the available hardware resources. 

* A job wrapper script that is suitable for the locally available batch scheduler needs to be used to call `oifs-run`.
* We include an example job wrapper script `run-oifs.ecmwf-hpc2020.job` in `$OIFS_HOME/scripts/exp3d`, which is suitable for the ECMWF hpc2020 HPC. This system uses the SLURM batch job scheduler.
    - As described above, this script is copied to the experiment directory because it needs to be located here, to run an experiment.
* `run-oifs.ecmwf-hpc2020.job` needs to be edited with the following essential and optional changes
    - Intially run-oifs.ecmwf-hpc2020.job sets the PLATFORM_CFG variable as follows:
    ```
    # set OpenIFS platform environment:
    PLATFORM_CFG="/path/to/your/config/oifs-config.edit_me.sh"
    ```
    - It is important to change "/path/to/your/config/oifs-config.edit_me.sh" to the actual path for the oifs-config.edit_me.sh, e.g., "$HOME/openifs/oifs-config.edit_me.sh"
    - You will need to adjust the batch scheduler header lines as required for your local system.
    - The `run-oifs.ecmwf-hpc2020.job` script will update the entries in `exp-config.h` with the correct settings dermined by the batch job headers.
    - For information, the LAUNCH command for batch job submission is set to "srun" without any further options, because all required parallel environment settings are provided through the SLURM script headers.

Once you have made the appropriate changes the job can be submitted:

```
# run as slurm batch job:
source oifs-config.edit_me.sh
cd $OIFS_EXPT/ab2a/2016092500
sbatch ./run-oifs.ecmwf-hpc2020.job
```

>NOTE: The job wrapper script will read the exp-config.h file and adopt the selected values. The exceptions are `LAUNCH`, which is set to "srun" for batch jobs, and `OIFS_NPROC` & `OIFS_NTHREAD` for which values from the batch job headers are used. The job wrapper script modifies the `exp-config.h` file accordingly prior to calling the `oifs-run` script.

#### Running interactively:

This example is shown for the ECMWF hpc2020, where running the model interactively **should be fine for lower grid resolutions** up to T255L91. 

* In order to run the experiment interactively, execute the `oifs-run` script from the command line in your terminal.
* If no command line parameters are provided with the `oifs-run` command, then the values from the `exp-config.h` will be used.
* In `exp-config.h` set `OIFS_NPROC=8` and `OIFS_NTHREAD=4`.
* In `exp-config.h` the `LAUNCH` variable should remain empty, i.e. `LAUNCH=""` and no `--runcmd` parameter should be provided in the command line.

The `oifs-run` script will in this case use its default launch parameters:  `srun -c${OIFS_NPROC} --mem=64GB --time=60`  which will work fine with `OIFS_NPROC=8` for experiment `ab2a`. 

```
# run interactively:
cd $OIFS_EXPT/ab2a/2016092500
./oifs-run
```

### Postprocessing the model output

Postprocessing creates a unique output folder and  groups all model output fields and diagnostics into individual GRIB files with ascending forecast time step. Also, a copy of the atmospheric model namelist file `fort.4`, as well as the `ifs.stat` and `NODE.01_001` log files are moved into the output folder.

Postprocessing can be done either directly by `oifs-run` or in a separate step by running the `run-pproc.ecmwf-hpc2020.job` script.

#### Postprocessing with oifs-run:

If in the `exp-config.h` file the `OIFS_PPROC` variable has been set to `true` (or if the `--pproc` command line parameter was used) then the model output in the experiment directory is further processed after completing the model run.

* In this case the script will generate a folder called `output_YYYMMDD_HHMMSS`, with YYYYMMDD being the current date and HHMMSS the current time. 
* This avoids accidental modification or overwriting of any previous results when the model experiment is repeated.
* For convenience a symbolic link output is set to the most recently generated model output. If the model run is repeated and a new `output_YYYMMDD_HHMMSS` folder is generated, the symbolic link will be updated to point to the most recent output folder.
* The variable `OUTPUT_ROOT` in `exp-config.h` determines where this ouput folder will be created. The default location is inside the experiment directory, but when assigning another path to `OUTPUT_ROOT` this could be created elsewhere.

#### Postprocessing with a separate script:

The postprocessing task is carried out by a single processor. During this time the batch job's computing resources will remain allocated for `oifs-run` which can be inefficient. 

As an alternative we include the script `$OIFS_HOME/scripts/exp3d/run-pproc.ecmwf-hpc2020.job` which completes the postprocessing in isolation after the model forecast experiment has ended. This script can be submitted as a serial job.

If you want to use this script follow these steps:

* In `exp-config.h` set `OIFS_PPROC=false` and run the forecast experiment to completion.
* Copy the run-pproc script to the `OIFS_EXPT` location: `cp $OIFS_HOME/scripts/exp3d/run-pproc.ecmwf-hpc2020.job $OIFS_EXPT`
* Edit the run-pproc script and update the variable `PLATFORM_CFG` as required for your installation.
* The script allows postprocessing of multiple experiments in sequence. Ensure the experiment IDs of all experiments you wish to process are included in the `EXP_LIST` heredoc. 
>NOTE: In our worked example you need the line entry `ab2a`.
* Submit the script with the `sbatch` command (or as appropriate for your system) as a serial job.


## Run a standard OpenIFS SCM case
