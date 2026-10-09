# Summary of OpenIFS environment variables

Building and running OpenIFS depends on a set of global environment variables. These variables are set in `oifs-config.edit_me.sh`, which can be found at the top level of your OpenIFS package. Here we describe some of the variables and their purpose.

* `OIFS_HOME` is the most important environment variable, which describes the location of the OpenIFS model installation and is generally the path where the git repository was extracted.  For example, if you cloned the repository as `openifs` inside your `$HOME` directory then you should set
  * `export OIFS_HOME=$HOME/openifs/`
* `OIFS_CYCLE` describes the OpenIFS model cycle which is named after the corresponding parent IFS cycle (e.g. 48r1).
* `OIFS_CLIMATE` describes the version of the static climate input data that can be used for the given cycle (e.g. `climate.v020`).
* `OIFS_EXPT` path to the location of the OpenIFS experiments. This is a requirement for the OpenIFS 3D model and the SCM experiments.
* `OIFS_ARCH` - if available for a system, this variable describes the location of the arch directory,  which provides specific information about the system and compiler, e.g. `$OIFS_HOME/arch/ecmwf/hpc2020/gnu`
  * Such a directory is not always required, i.e., if a system has all the appropriate libraries installed. If this is the case `OIFS_ARCH` can be set to an empty string, i.e., `OIFS_ARCH=""`
* `OIFS_DATA_DIR` - describes the location of climatological input files that are required to run OpenIFS. These have been installed on the ECMWF HPC in a central and accessible location and the information is organised by model cycle.
  * If you do not have access to the ECMWF hpc2020 file system, or if you wish to install the climatological input files in a local directory of your choice, then you can download the required data for your model cycle from this site: [OpenIFS download: ifsdata](https://openifs.ecmwf.int/data/ifsdata/).
    * As a minimum you will require the packages `${OIFS_CYCLE}/rtables/rtables.tar.gz` and `${OIFS_CYCLE}/ifsdata/ifsdata.tar.gz`, where `${OIFS_CYCLE}` matches the environment variable.
    * You will also need to download the package for each of your selected horizontal grid resolutions from `${OIFS_CYCLE}/${OIFS_CLIMATE}/`. For example, for a T255 grid at model cycle 48r1 you will at the very least need the package `48r1_climate.v020_255.tar.gz`. Installing the packages for all supported horizontal grids will require a lot of disk space and is therefore not recommended.
    * Download and extract the files in these tarballs into the same directory in your chosen location, and set the variable `OIFS_DATA_DIR` to the path of this directory.
* `OIFS_EXEC` describes the location and file name for the OpenIFS 3D executable, e.g. `$OIFS_HOME/build/bin/ifsMASTER.SP`, which is the single-precision executable for the 3D model.

> `$OIFS_HOME/oifs-config.edit_me.sh` also contains environment variables that relate to the Single-Column Model (SCM).