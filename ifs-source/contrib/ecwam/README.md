ecWAM
*****

The ECMWF wave model ecWAM

Introduction
============

The ECMWF Ocean Wave Model (ecWAM) describes the development and evolution of wind generated surface waves and their height, direction and period.
ecWAM is solely concerned with ocean wave forecasting and does not model the ocean itself: dynamical modelling of the ocean can be done by an ocean model such as NEMO.

- ecWAM may be used as a standalone tool that can produce a wave forecast driven by external forcings provided via GRIB files.
- Alternatively it can be used in a coupled mode where it provides feedback and receives forcings from
  * the atmospheric forecast model IFS.
  * the dynamic ocean model NEMO.

For more information, please go to https://confluence.ecmwf.int/display/FUG/2.2+Ocean+Wave+Model+-+ECWAM

License
=======

ecWAM is distributed under the Apache License Version 2.0.
See `LICENSE` file for details.

Installing ecWAM
================

Supported Platforms
-------------------

- Linux
- Apple MacOS

Other UNIX-like operating systems may work too out of the box.

Requirements
------------
- Fortran and C compiler, and optionally C++ compiler
- CMake (see https://cmake.org)
- ecbuild (see https://github.com/ecmwf/ecbuild)
- fiat (see https://github.com/ecmwf-ifs/fiat)
- eccodes (see https://github.com/ecmwf/eccodes)
- field_api (see https://github.com/ecmwf-ifs/field_api)

Further optional dependencies:
- MPI Fortran libraries
- multio (see https://github.com/ecmwf/multio)
- ocean model (e.g. NEMO or FESOM)
- fypp (see https://github.com/aradi/fypp)

Some driver scripts to run tests and validate results rely on availability of:
- md5sum (part of GNU Coreutils; on MacOS, install with `brew install coreutils`)
- Python with pyyaml package

Generating the derived-type data structures (read end of next section) requires:
- Python with pyyaml package
- fypp

Building ecWAM
--------------

Environment variables

    $ export ecbuild_ROOT=<path-to-ecbuild>
    $ export MPI_HOME=<path-to-MPI>
    $ export fiat_ROOT=<path-to-fiat>
    $ export eccodes_ROOT=<path-to-eccodes>
    $ export field_api_ROOT=<path-to-field_api>
    $ export CC=<path-to-C-compiler>
    $ export FC=<path-to-Fortran-compiler>
    $ export CXX=<path-to-C++-compiler>

If you want to pre-download or install extra data files, run in the source-directory before CMake (re)configuration:

    $ share/ecwam/data/populate.sh

You must compile ecWAM out-of-source, so create a build-directory

    $ mkdir build && cd build

Configuration of the build happens through standard CMake

    $ cmake <path-to-source>

Extra options can be added to the `cmake` command to control the build:

 - `-DCMAKE_BUILD_TYPE=<Debug|RelWithDebInfo|Release|Bit>` default=RelWithDebInfo (typically `-O2 -g`)
 - `-DENABLE_TESTS=<ON|OFF>`
 - `-DENABLE_MPI=<ON|OFF>`
 - `-DENABLE_OMP=<ON|OFF>`
 - `-DCMAKE_INSTALL_PREFIX=<install-prefix>`

More options to control compilation flags, only when defaults are not sufficient

 - `-DOpenMP_Fortran_FLAGS=<flags>`
 - `-DCMAKE_Fortran_FLAGS=<fortran-flags>`
 - `-DCMAKE_C_FLAGS=<c-flags>`

Once this has finished successfully, run `make` and `make install`.

An informational tool `ecwam [--help] [--info] [--version] [--git]` is available upon compilation
and can be used the to verify compilation options and version information of ecWAM.

Optionally, tests can be run to check succesful compilation, when the feature TESTS is enabled (`-DENABLE_TESTS=ON`, default ON)

    $ ctest

### Generate derived-types data structures
The derived-types storing grid-point data in ecWam can be configured in `src/ecwam/yowfield_mod_config.yaml`. If the fypp preprocessor and Python + pyyaml are found, then the configuration file is used to expand the accompanying `src/ecwam/yowfield_mod.fypp` into Fortran derived-type objects. The glue-code required to turn the derived-types members into FIELD API objects is also generated. If fypp and pyyaml are not found, then the existing `src/ecwam/yowfield_mod.F90` is used. Generation of the derived-type data structures can be disabled by passing the following build time option: `-DENABLE_GEN_DERIV_TYPES=OFF`.

## Build using ecWAM bundle

Another way of building ecWAM is to use the bundle definition included in `package/bundle`:

    $ ./package/bundle/ecwam-bundle create  --bundle package/bundle/bundle.yml # Checks out dependency packages
    $ ./package/bundle/ecwam-bundle build [--build-type=<build-type>] [--arch=<path-to-arch>] [--option]

The bundle also facilitates setting environment variables and compiler flags relevant to certain architectures by specifying the corresponding arch file at the build step. For example, to build on the ECMWF Atos system using Intel compilers and the hpcx-openmpi `MPI` library:

`--arch=package/bundle/arch/ecmwf/hpc2020/intel/2021.4.0/hpcx-openmpi/2.9.0`

The following options can also be configured during the bundle build step:
 - `--without-mpi` - Disable MPI
 - `--without-omp` - Disable OpenMP
 - `--single-precision` - Build single-precision variant of ecWAM

 Finally, additional `CMake` options can also be set during the bundle build step:

`--cmake-"OPTION=<arg>"`

Running ecWAM
=============

Following are instructions to run ecWAM as a standalone wave forecasting tool.

A YAML configuration file is required. See [tests](tests) directory for examples.

To run, use the commands listed in following steps, located in the build or install `bin` directory.
Each of following commands can be inquired with the `--help` argument.
If `--run-dir` argument is not specified, or `ECWAM_RUN_DIR` environment variable is not set,
then the current directory is used.
If `--config` argument is not specified, it is assumed that a file called `config.yml`
is present in the `<run-dir>`.


1) Create bathymetry and grid tables

```shell
ecwam-run-preproc --run-dir=<run-dir> --config=<path-to-config.yml>
```

This command generates bathymetry data files as specified by configuration options.
As bathymetry data files are large and require heavy computations they are
cached for later use in a directory which can be chosen with the `--cache` argument, or
`ECWAM_CACHE_PATH` environment variable.
By default the cache path will be `$HOME/cache/ecwam` unless on the ECMWF HPC it is in
`$HPCPERM/cache/ecwam`.
Bathymetry data files can also be searched for in a hierarchy of cache-like directories
specified with the `ECWAM_DATA_PATH` variable containing a ':'-separated list of paths
(like `$PATH`). If not found, they are attempted to be downloaded from URL
https://get.ecmwf.int/repository/ecwam. If still not available, they will be computed.
The cache path will then be populated with computed, or downloaded data,
or with symbolic links to found data in the `ECWAM_DATA_PATH`s.

Grid tables are always computed and never cached. THey are placed in the `<run-dir>`

2) Create initial conditions

```shell
ecwam-run-preset --run-dir=<run-dir> --config=<path-to-config.yml>
```

As a result files, binary files of the form `<run-dir>/restart/BLS*` and `<run-dir>/restart/LAW*` are created.
They contain all initial conditions required for the wave model run from "cold start".
This command requires surface wind and sea-ice-cover input, at initial simulation time, provided in GRIB format.
The configuration file must specify this file. For several benchmark or test cases,
they are retrieved in similar fashion as the bathymetry files (see above).

This package also contains some scripts to generate MARS requests to retrieve data from the ECWMF operational
forecast or from the ERA5 reanalysis data set. This is useful to generate new tests or for longer runs.

3) Run wave model

```shell
ecwam-run-model --run-dir=<run-dir> --config=<path-to-config.yml>
```

With initial conditions, forcings, and grid tables in place we can run the actual wave model.
The advection and physics time step needs to be configured via the configuration file in accordance
the grid resolution.
The configuration file offers options to output GRIB output fields at regular time intervals, or
binary restart files similar to the initial condition files generated in step 2.
After the run, the output files will be in `<run-dir>/output/` and log files will be in `<run-dir>/logs/model`.
One log file called `statistics.log` contains computed norms which can be used to validate results.
Such validation will occur automatically when the configuration file contains a `validation` section.
See [tests](tests) directory for example configuration files.

Running with MPI
----------------

Above commands by default run without MPI.
To use MPI, or apply a custom command prefixed to the binary execution,
there are following options:

- Use argument `--launch="ecwam-launch -np <NTASKS> -nt <NTHREADS>"`.
  `ecwam-launch` is an internal "smart" launcher that chooses a good launcher
  depending on availability and the used platform. It will also set
  `export OMP_NUM_THREADS=<NTHREADS>` for you. The unit-tests are automatically
  configured to use this when invoked via `ctest`.

- Use arguments `-np <NTASKS> -nt <NTHREADS>`. This is equivalent to the above, and
  internally uses the `ecwam-launch` launcher.

- Use any other custom command, e.g. `--launch="srun -n <NTASKS> -c <NTHREADS>"` for full control.
  Note that this does *not* automatically export the `OMP_NUM_THREADS` variable.

Note that only `ecwam-run-model` currently supports MPI.

Known issues
============

1) On macOS arm64 with gfortran 12.2, and Open MPI 4.1.4, and with compilation
   with flag `-ffpe-trap=overflow`, the execution of `ecwam-preproc` and `ecwam-chief`
   needs to be launched with `mpirun -np 1`, even for serial runs in order to avoid
   a floating point exception during during call to `MPI_INIT`.
   The flag `-ffpe-trap=overflow` is set e.g. for `Debug` build type.
   Floating point exceptions on arm64 manifest as a `SIGILL`.


Reporting Bugs
==============

Please report bugs using a [GitHub issue](https://github.com/ecmwf-ifs/ecwam/issues).
Support is given on a best-effort basis by package developers.

Contributing
============

Contributions to ecWAM are welcome.
In order to do so, please open a [GitHub issue](https://github.com/ecmwf-ifs/ecwam/issues) where
a feature request or bug can be discussed.
Then create a [pull request](https://github.com/ecmwf-ifs/ecwam/pulls) with your contribution.
All contributors to the pull request need to sign the
[contributors license agreement (CLA)](http://claassistant.ecmwf.int/ecmwf-ifs/ecwam).
