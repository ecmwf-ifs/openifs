# Summary of OpenIFS environment variables

Building and running OpenIFS depends on a set of global environment variables. These variables are set in `oifs-config.edit_me.sh`, which can be found in the top-level of your OpenIFS package. Here we describe some of the variables and their purpose

* `OIFS_HOME` is the most important environment variable, which describes the location of the OpenIFS model installation and is general the path where the git repository was extracted.  For example, if you named clone the repository `openifs` and cloned it into your `$HOME` directory then you should set
  * `export OIFS_HOME=$HOME/openifs/`
* `OIFS_CYCLE` describes the model input/experiment data for a given cycle (e.g. 48r1) for which this model configuration can be used.
* `OIFS_EXPT`  path to the location of the openifs experiments. This is a requirement for the OpenIFS 3D model and the SCM experiments.
* `OIFS_ARCH` - if available for a system, this variable describes the location of the arch directory,  which provides specific information about the system and compiler, e.g. `$OIFS_HOME/arch/ecmwf/hpc2020/gnu`
  * Such a directory is not always required, i.e., if a system has all the appropriate libraries installed. If this is the case `OIFS_ARCH` can be set to an empty string, i.e., `OIFS_ARCH=""`
* `OIFS_DATA_DIR` - describes the location of climatological input files that are required to run OpenIFS. These have been installed on the ECMWF HPC in a central and accessible location and the information is organised by model cycle.
  * If you do not have access to the ECMWF hpc2020 file system, or if you wish to install the climatological input files in a local directory of your choice, then you can download the required data from this site: [OpenIFS data: ifsdata](https://sites.ecmwf.int/openifs/openifs-data/ifsdata/).
    * As a minimum you will require the packages `ifsdata_rtables_<OIFS_CYCLE>.tgz` and `ifsdata_climatology_<OIFS_CYCLE>.tgz`, where `OIFS_CYCLE` matches the environment variable.
    * You will also need to download the package for your selected horizontal grid resolution, which in the case of our worked 48r1 case, for example, is for a T255 grid, and therefore you will at the very least need the package `ifsdata_48r1_climate.v020_255.tgz`. Installing the packages for all supported horizontal grids will require a lot of disk space and is therefore not recommended.
    * Download and extract the files in these tarballs into the same directory in your chosen location, and the filepath of this directory should then be set to variable `OIFS_DATA_DIR`.
* `OIFS_EXEC` describes the location and file name for the OpenIFS 3D executable, e.g. `$OIFS_HOME/build/bin/ifsMASTER.DP` , which is a double-precision executable for the 3D model.

> `$OIFS_HOME/oifs-config.edit_me.sh` also contains environment variables that relate to the Single-Column Model (SCM).