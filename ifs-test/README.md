# Overview

This is the IFS QA Test Suite.

# Usage

To run all available tests, use `ctest` in the `ifs-test` build folder to 
launch them. The tests can also be launched with more control, by running the
`ifs-run` or `ifs-run.py` scripts in the corresponding test folders.

# Environment variables

A couple of environment variables can be used to control the test procedure:

* `IFS_TEST_BITIDENTICAL` 
  * If set to 1/init, reference results are written to the test folders.
  * If set to 2/check, the tests compare the results to the stored reference
    results.
  * If not set or 0, any result validation is disabled.
* `IFS_TEST_RUN_LARGE` is used to enable (set to 1) or disable (0) all larger 
  tests (currently >= tl159).
* `IFS_TEST_TOLERANCE` controls the tolerance when validating results (only used
  if `IFS_TEST_BITIDENTICAL=2`. If the value is <=0, results are validated by
  checking for bit-identicality. Otherwise, the given tolerance is used when
  comparing the results.
* `IFS_TEST_LAUNCHER_FLAGS` are additional flags that are passed to the job
  scheduler when starting a test.
* `IFS_TEST_LEGACY`: If set to any other value than 0, this activates the 
  legacy mode where results are compared using the old ifs-grep-norms.pl
  script when checking for bit-identicality.
  Not compatible with IFS_TEST_RUN_LARGE=1 or IFS_TEST_TOLERANCE!=0.
* Aditional options can be passed to the ifsbench-based tests by specifying the
  `IFSBENCH_<PARAM_NAME>` variable, e.g. `IFSBENCH_ARCH=lumi_g`.

# CMake flags

* `IFS_TEST_TIMEOUT`: Specify a default timeout for the tests (set to negative
  value to disable timeouts).
* `IFS_TEST_CUSTOM_DATA_SOURCE`: Specify a directory on the filesystem which
  will be used as the primary data source for fetching input data.
* `IFS_TEST_FORCE_FETCH=ON/OFF`: Enable/disable fetching of external input data
  if input data is already available. Default: On. 

# Data management

Binary files are stored externally (usually in ECPDS) and are automatically pulled 
in by CMake. 

## Fetching files

Files are automatically fetched from ECPDS in CMake, which uses the 
`bin/storage.py fetch` command to place the files in `build_dir/.cache`. The
fetch will move the files to the directories that are specified in `storage.yaml`
(`target_path` for files, `extract_path` for the content of archives).  

## Caching files

The files from ECPDS can also be cached locally. This can be done by 
running the `bin/boostrap` script. This installs a Python virtual
environment with all necessary dependencies into `.boostrap/venv` and calls
`bin/storage.py cache` to put all files that are specified in `storage.yaml`
into the `.boostrap/cache` directory. By using `bin/storage.py cache` directly,
the data may also be put in a different place. 
The `bin/bootstrap` or `bin/storage.py cache` commands put the ECPDS data into
the cache directory, keeping the original directory structure - the files won't
be moved to `target_path/extract_path`. 
 
Please note that the bootstrap script (and `bin/storage.py`) require at least
Python 3.8.

To use the cached data when building `ifs-test`, set the 
`IFS_TEST_CUSTOM_DATA_SOURCE` CMake flag to the cache directory
(`source_dir/.bootstrap/cache` by default).


## Adding new files

Adding a new file that should be fetched from ECPDS can be done in these steps:

1. Create the new file.
2. Add the file to the `files` section in `storage.yaml`. `source_path`
   specifies the relative location in ECPDS, `target_path` the location where
   it will be placed after being fetched from ECPDS. If your file is an
   archive, you can also specify `extract_path` which is the directory where
   the extracted files will be placed.
3. Add the actual file to ECPDS (see _Adding files to ECPDS_). If you do not have
   write access to ECPDS, please specify the files that should be uploaded to
   ECPDS in your pull request so that an IFS maintainer can upload it.


## Updating archives

If you have to update the content of an archive, the following workflow can be
used:
1. Build `ifs-test`, either via `git ifstest` or `ifs-bundle`.
2. Check `storage.yaml` for the `extract_path` that belongs to the archive that
   you want to updated.
3. The `build_dir/ifs-test/.cache` directory contains the files that were
   fetched from ECPDS, as well as the extracted archives. 
   Update the files in `build_dir/ifs-test/.cache/<extract_dir>`.
4. Re-run the tests to check that your updated files are working as intended.
5. Return to the ifs-test source directory and update the `source_path` in 
   `storage.yaml` for the corresponding archive. All archives have names like 
   `YYYY_MM_DD_some_name` so update the date component to today's date.
6. Repack the data by using `bin/storage.py` and put the resulting archives
   into `output_dir`:
   ```
   python3 bin/storage.py pack storage.yaml build_dir/ifs-test/.cache output_dir
   ```
   You can also use the `--match` flag to limit the packing to certain
   archives.
   ```
   python3 bin/storage.py pack --match ifsdata storage.yaml build_dir/ifs-test/.cache output_dir
   ```
    
   Please note: This requires a `python3` installation that includes the
   `click` and `yaml` modules. If you have run the `bin/bootstrap` script, you
   can use `ifs-test-source/.bootstrap/venv/bin/python3`.
7. Add the archive in `output_dir` to ECPDS (see _Adding files to ECPDS_). 
   If you do not have write access to ECPDS, please specify the files that
   should be uploaded to ECPDS in your pull request so that an IFS maintainer
   can upload it.


## Adding files to ECPDS

If new binary files are needed, these have to be uploaded to ECPDS. On
hpc2020 this is done like this:
```
module load mspds
mspds -echost aux -destination ifs_test_inidata -lifetime 4000d -source <filename> -target <cycle>/<type>/<date>_<filename>

# Example
mspds -echost aux -destination ifs_test_inidata -lifetime 4000d -source 2024-02-20-ifsdata.tar.gz -target 49r1/shared/2024_02_20_ifsdata.tar.gz
``` 
Specifying the lifetime is important as this determines how long the data will reside in ECPDS (and how long it can be downloaded via wget)!

Write access to the `ifs_test_inidata` ECPDS bucket is limited to a number of
IFS maintainers. If you lack the write permissions, please state which files
must be uploaded to ECPDS in your pull request.

# Documentation

To build the Ford documentation for this repository:

1. ```cd doc```
1. ```rm -rf ./html``` (if this output directory already exists, the build fails)
1. ```./make_ford.sh```
1. ```firefox ./html/index.html```
