# OpenIFS build options (more detail)

The OpenIFS build system is based on a combination of [ecbuild](https://github.com/ecmwf/ecbuild), which is an extension of the widely-used [CMake metabuilder](https://cmake.org/), and [ecbundle](https://github.com/ecmwf/ecbundle) software packages, which are both developed by ECMWF and used in the ECMWF IFS model.

In a similar fashion to the IFS, OpenIFS employs a software bundle through the `bundle.yml`, which defines the software packages and versions that are required for the build, i.e., ifs-source and ifs-test version plus all the associated software libraries, such as eccodes, multio, atlas etc. 

The build of the OpenIFS through invoking `ecbuild` and `ecbundle` and using the defined software bundle in the bundle.yml is controlled by the `openifs-test.sh` script, which is a derivative of the `ifs-git-tools` script, `ifstest`.

A standard build of OpenIFS can be achieved with the following

```bash
source /path/to/file/location/oifs-config.edit_me.sh
#
# e.g if OpenIFS is stored a directory openifs in a home directory:
# source $HOME/openifs/oifs-config.edit_me.sh
#
cd $OIFS_HOME
$OIFS_TEST/openifs-test.sh -cbt
```

where

`$OIFS_TEST` is defined in the platform configuration file (`oifs-config.edit_me.sh`) as `$OIFS_HOME/scripts/build_test` 

and the command line options are

* `-c`  creates `source` directory in `$OIFS_HOME`, which is used to collect all the sources defined in the `bundle.yml`, in preparation for the build
* `-b`  builds the source. This step creates the directory `build` in `$OIFS_HOME`, which is used to build and store the OpenIFS and SCM executables.
* `-t` will run the ifs-test t21 tests, which comprise 21 coarse resolution (t21) 3-D NWP tests with and without chemistry and 1 SCM test (based on TWP-ICE)

This command and the associated options create and store 3-D OpenIFS double and single precision master executables, `ifsMASTER.DP` and `ifsMASTER.SP`, respectively. These executables are stored in `$OIFS_BLD_PARENT/bin`, which by default points to `$OIFS_HOME/build/bin`. Further, the location and name of the executable for OpenIFS experiments is defined as `$OIFS_EXEC`, which by default is the single precision executable, i.e., `${OIFS_BLD_PARENT}/bin/ifsMASTER.SP`.

The standard build command also creates and stores the double and single precision Single Column Model (SCM) executables (`MASTER_scm.DP` and `MASTER_scm.SP`, respectively), which are also located in `$OIFS_HOME/build/bin` (`$OIFS_BLD_PARENT/bin`). The name and location of the SCM executable is defined by `$SCM_EXEC`, which, by default, is the single precision executable, i.e., `${OIFS_BLD_PARENT}/bin/MASTER_scm.SP`.

## Additional openifs-test build options

In addition to the `-cbt` options, the following additional options exist

* `-j` allows the user to define the number of threads used to create and build openifs, e.g.
  
    ```bash
    $OIFS_TEST/openifs-test.sh -cbt -j 16
    # This will build OpenIFS using 16 threads
    ```

    * If `-j` is not present as a command line argument the number of threads will default to:
      * 64 threads for ATOS
      * 8 threads for all other systems
  * If `-j` is present, then the command line number of threads takes precedence over the defaults for any system
    * The default of 8 threads can be increased for larger systems. Such an increase may speed up the build. For advice on this please liaise with the system administrator.
    * It is important to note that 8 threads can be problematic for older and lower spec systems. For example, it will not work with a 4 core intel m3. In this case reduce the number of threads to 2 for example.

* `-u` signals the `openifs-test.sh` to create (`-c`) the source directory using symbolic links that point to the central location for the non-OpenIFS sources, e.g. `eccodes`, `multio`, `fckit`, `eckit`, `metkit`, `fdb5`, `atlas` and `fcm`
* `-i` uses the `install.sh`, which is created by `ecbuild`, to install the OpenIFS in an install directory.

`-u` and `-i` are useful options for significantly reducing the file numbers stored on a system following a successful installation and build of OpenIFS. For more information, please refer to [Minimal OpenIFS install](oifs_howto_setup_central_source.md).

## ecbundle options

### Default build options

There are 5 default build options that are defined directly in `openifs-test.sh`. Of the 5, 3 options cannot be negated by command line options, while 2 can be negated by selecting the opposite option from the command line.

The following build options cannot be negated

* `--openifs-only` 
  * ensures code for OpenIFS is included in the build
* `--with-scmec` 
  * builds the SCM executable alongside the OpenIFS 3D executable
* `--arch=$OIFS_ARCH`
  * defines the path to the arch directory defined in $OIFS_HOME/oifs-config.edit_me.sh

The following default build options, can be negated from the command line by selecting the opposite option 

* `--with-single-precision`
  * Ensures OpenIFS is built in single precision, as well as the default double precision
* `--init-snan`
  * initialises uninitialised variables to nan

### Command line ecbundle options

In addition to the default options, there is a range of command-line options defined in `$OIFS_HOME/bundle.yml`, which can be invoked when `openifs-test.sh -b` is executed.

The main available options that most users may want to use are

* `--clean`
  * force a clean build by removing the `$OIFS_HOME/build` directory. 
  * This is particularly useful when changing cmake files or adding new functionality, but it does mean the build takes longer. 
  * If there is no --clean and there is a build, then the build will be incremental (and quicker). 
* `--arch=<add path of arch file>` 
  * Allows a user to over-ride `$OIFS_ARCH` loaded while sourcing `oifs-config.edit_me.sh`. 
  * The expected path structure for `<add path of arch file>`  is `./arch/<site>/<platform>/<compiler>/<compiler_version>`. For example `--arch=./arch/ecmwf/hpc2020/gnu` , where the site is ecmwf, the platform is hpc2020 and the compiler is gnu. The compiler version is not stated because the path contains a default link.
* `--without-single-precision` or `--without-double-precision`
  * By default both single and double precision executables are built. These options negate the default option by excluding double or single precision build, which speeds up the build
* `--build-type=DEBUG`
  * Builds the executable in debug mode with bounds checking, some trapping and no compiler optimisation
* `--no-init-snan`
  * negates the default `init-snan` option. It is recommended that this should only be used for debugging issues.