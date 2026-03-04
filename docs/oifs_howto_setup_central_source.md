# OpenIFS build : Setup and Use a central location for source

## Introduction

With the release of OpenIFS 48r1, the build system for OpenIFS moved from [Flexible Configuration Management (FCM) system](https://metomi.github.io/fcm/doc/user_guide/introduction.html) to [ecbuild](https://github.com/ecmwf/ecbuild) system. This change in build system, along with some code changes[^1], enables much closer alignment of OpenIFS with IFS and, importantly the ability to run and develop OpenIFS within the IFS framework. [ecbuild](https://github.com/ecmwf/ecbuild) is CMake-based build system that consists of collection of CMake macros and functions[^2], the "ease the managing of software build systems".

ecbuild, just like CMake, uses out-of-source builds. In the context of OpenIFS (and IFS), the build is based on a defined software bundle for a given cycle, in which versions of the separate required software, e.g. eccodes, openifs-source, eckit, are defined. The bonus of the bundle system, is that it is ensures that a user is building consistent libraries for a given cycle, which significantly simplifies the initial set-up of OpenIFS on a user system. The downside of the bundle build is that each new build will, by design, clone the repositories or link to local sources, mainly IFS/OpenIFS source,  of the defined bundle software. The cloning of repositories creates the `source` directory and is part of the `openifs-test` create step. The software packages in `source` will be built in the `build` directory. Both `source` and `build` are created by `ecbuild`, if they are not present.

The process of cloning to `source` and building in `build` results in a lot of files (\<50000) for each IFS branch or OpenIFS package that is built. This can cause problems for systems that have low file quotas. This how-to presents details about the standard build and then some extra build options to reduce the number of files per build.

## Develop branch

<!-- # TODO: git.ecmwf.int references need fixing -->

The functionality described below can be found in the [OpenIFS-48r1 develop branch](https://git.ecmwf.int/projects/OIFS/repos/openifs-48r1/browse?at=refs%2Fheads%2Fdevelop). This branch is the latest tagged release of OpenIFS-48r1, with the merge of [Pull-Request for OIFS-605](https://git.ecmwf.int/projects/OIFS/repos/openifs-48r1/commits/f594a0c8269b24c4565869d8f1f7a3a767402dc7#scripts%2Fbuild_test%2Fopenifs-test.sh).

## Standard build

Since the OpenIFS 48r1, the build is controlled using the [`openifs-test.sh`](https://git.ecmwf.int/projects/OIFS/repos/openifs-48r1/browse/scripts/build_test/openifs-test.sh) script, which is derivative of the ifs-git-tools, ifstest script.

Details about how to perform a standard build with `openifs-test.sh` can be found in [OpenIFS 48r1 User Guide - 3. Build OpenIFS](https://confluence.ecmwf.int/display/OIFS/Getting+started#Gettingstarted-BuildOpenIFS).

In addition to these instructions, `openifs-test.sh` has recently been expanded to include an install step, using `-i` option, i.e.,

`$OIFS_TEST/openifs-test.sh -cbti`

where

- `-c` - creates the source repository and clones all the required software defined in the `bundle.yml`
- `-b` - performs a standard build
- `-t` - run standard tests to check they are passing
- `-i` - NEW - uses the `install.sh`, which is created by ecbuild, to install the OpenIFS in an install directory.

The standard build retains (and uses) the source, build and more recently, install directory, which means that subsequent builds are a lot faster, because `ecbuild` can do incremental builds. However, source and build result in > 50000 files for a build, which, can be problematic for some systems.

## Setup a central location for non-OpenIFS sources

An alternative, reduced file build option is to create a central location for the bundled software packages that most users will never change. For example, with OpenIFS 48r1 these packages are  `ecbuild`, `eccodes`, `multio`, `fckit`, `eckit`, `metkit`, `fdb5`, `atlas` and `fcm`. A central source location can be set-up running

- [openifs-48r1/scripts/build_test/setup_central_source.sh](https://git.ecmwf.int/projects/OIFS/repos/openifs-48r1/browse/scripts/build_test/setup_central_source.sh?at=refs%2Fheads%2Ffeature%2FOIFS-605-implement-install-script)
  - At the time of writing this script is on a feature branch, awaiting merge.

To use `setup_central_source.sh` do the following:

- Set up the OpenIFS environment with `source /path/to/<openifs_package_dir>/oifs-config.edit_me.sh`, e.g. `source ~/openifs-48r1/oifs-config.edit_me.sh`
  - As well as the normal changes to oifs-config.edit_me.sh, it is important to set the following
    - `OIFS_CENTRAL_SRC` - this is the path for the top-level central location, default is `$HOME/openifs-bundle-src`. This may need to be changed depending on your system.
    - `OIFS_CYCLE` - A directory names OIFS_CYCLE will be created in OIFS_CENTRAL_SRC, e.g., for openifs-48r1, OIFS_CYCLE is 48r1, so by default the directory for the central location will be `$HOME/openifs-bundle-src/48r1`
- Execute `setup_central_source.sh` with
  - `./scripts/build_test/setup_central_source.sh` or `$OIFS_TEST/setup_central_source.sh`
  > Note: This script will download the non-OpenIFS sources to `$OIFS_CENTRAL_SRC/$CYCLE`.

Once complete, the central source will be installed and the build of OpenIFS can start.

> Note: If setting up a central package location for multiple users, `setup_central_source.sh` only needs to be run by one user, probably sysadmin and the central location has to be readable by all OpenIFS users.

<details>
<summary>Click to expand for more details about how `setup_central_source.sh` works.</summary>

- create the directory `$OIFS_CENTRAL_SRC/$OIFS_CYCLE`
- Execute `$OIFS_TEST/openifs-test.sh -c`
  - This will create the standard `$OIFS_HOME/source` directory.
- Once `$OIFS_TEST/openifs-test.sh -c` is complete, all the directories in `$OIFS_HOME/source` that are not symbolic links are copied to ``$OIFS_CENTRAL_SRC/$OIFS_CYCLE``
  - At 48r1, the packages copied are  `ecbuild`, `eccodes`, `multio`, `fckit`, `eckit`, `metkit`, `fdb5`, `atlas` and `fcm`
  - `setup_central_source.sh` checks that the directories copied are then same as a saved list. If not, then a warning will be logged. It is worth paying attention to the log output, particularly when changing the cycle.
- Finally, `$OIFS_HOME/source` is removed.

</details>
</br>

## Build OpenIFS using the central source

To Build OpenIFS using the central location for the non-OpenIFS sources, first use the following command

`$OIFS_TEST/openifs-test.sh -c -u`

where `-u` is a new argument that, rather than downloading/cloning the non-OpenIFS sources, signals the openifs-test to create (`-c`) the source directory using symbolic links that point to the central location for the non-OpenIFS sources, i.e., the `-u` invokes the following function,

``` bash
set_ifs_bundle_dir () {

  mkdir -p ${OIFS_HOME}/source
  #
  CENTRAL_SOURCE=${OIFS_CENTRAL_SRC}/${OIFS_CYCLE}
  #
  export IFS_BUNDLE_ECBUILD_DIR=$CENTRAL_SOURCE/ecbuild
  export IFS_BUNDLE_ECCODES_DIR=$CENTRAL_SOURCE/eccodes
  export IFS_BUNDLE_MULTIO_DIR=$CENTRAL_SOURCE/multio
  export IFS_BUNDLE_FCKIT_DIR=$CENTRAL_SOURCE/fckit
  export IFS_BUNDLE_ECKIT_DIR=$CENTRAL_SOURCE/eckit
  export IFS_BUNDLE_METKIT_DIR=$CENTRAL_SOURCE/metkit
  export IFS_BUNDLE_FDB5_DIR=$CENTRAL_SOURCE/fdb5
  export IFS_BUNDLE_ATLAS_DIR=$CENTRAL_SOURCE/atlas
  export IFS_BUNDLE_FCM_DIR=$CENTRAL_SOURCE/fcm
  #
}
```

where `OIFS_CENTRAL_SRC` is either set in `oifs-config.edit_me.sh` or defaults to `$HOME/openifs-bundle-src`

This will create a source directory that includes a symbolic links pointing to your central location for non-OpenIFS sources, plus the standard links for ifs-source, ifs-test, ifs_sp, ifs-dp

Once created then execute the normal OpenIFS build command:

`$OIFS_TEST/openifs-test.sh -<cu>bti --<plus any command line options>`

> Note: `-cu` is only required if this has not been run already and the source directory needs to be created or updated

Option `-i` is not always needed but has to be used if there is an intention to reduce the file count by removing the `build` directory following a successful install.

Once an install directory has been successfully produced using `-i`, then the `build` directory can be removed from a user directory, which reduces the file count from ~26000 to ~300.

> Note: Once the build directory has been removed, incremental builds are no longer an option, so any rebuild has to be a full build, which will increase the file count

[^1]: The required code changes consist of changes to `.cmake` files and the inclusion of preprocessor statements in the IFS fortran code, which permit OpenIFS to be built directly from the IFS using a build option `--openifs-only`. These changes are in a maintained branch for IFS CY48R1 and have been merged into the IFS main at IFS CY49R1.

[^2]: See <https://ecbuild.readthedocs.io/en/latest/> for list of and details about the various macros and functions.
