# ifsbench based tests

Some of the tests in ifs-test are based on the `ifsbench` (Python) package. For
these tests, the input data is stored in some central location and the
functionality in the `ifsbench` package is used to build an IFS-compatible run
directory.


## Experiment specification

Each such test is described by a `.yaml` file with the following content:

* `exp`: 
  * `local_path`: The path to the `.yaml` file that was created by `ifs-bench.py
    unpack-experiment`. This `.yaml` file essentially contains the paths to the
    actual data.

* `env` (optional): This block may contain arbitrary key/value pairs that are loaded
  as environment variables when running the actual IFS executable.

* `run_opts` (optional): The parameters given in this block are used as the 
  default values for the ifsbench run options. The most common values are 
  * `nproc`: The number of processes that are used.
  * `nthread`: The number of threads per process.
  * `tasks_per_node`: The number of MPI processes per node.

* `relative_paths` (optional): In some cases, some paths in the experiment must 
  be altered in order to work. `relative_paths` is a yaml list which entries consist of three 
  components which are then converted internally to a 
  `ifsbench.SpecialRelativePath` object.

## Running the tests

The tests can be run in three different ways:
* By using the ifs-run.py runner.
* By using CTest.
* By calling git ifstest.

### Using the ifs-run runner

The central script here is the `ifs_run.py` script that is placed in the 
`build_dir/bin` directory when `ifs-test` is built. In the simplest case, 
one only has to pass the `from_yaml` command and the experiment file to this 
script to run it:
```
# Run the tco399 test.
ifs-run.py from_yaml tco399.yml
```
On the Bologna machine, it may be necessary to specify that the tests should be
executed on the compute nodes, instead of the general purpose nodes by adding
'-q np' to the launch options (otherwise, the necessary number of cores can't
be allocated):
```
# Run the tco399 test.
ifs-run.py from_yaml tco399.yml --launch-options='-q np'
```

Various options can be passed to this script, for example the number of MPI 
processes or the directory where the tests are run:
```
# Run in a given directory (temporary directory otherwise).
ifs-run.py from_yaml tco399.yml --nproc=256 --run-dir=/hpcperm/myhome/tmp_tco399
```

One can also use this script to generate reference results for this test...
```
# Verify results against a given reference.
ifs-run.py from_yaml tco399.yml –r /some/reference/path --update-reference
```

... and use this reference data to verify the results later on
```
# Verify results against a given reference.
ifs-run.py from_yaml tco399.yml –r /some/reference/path --validate
```

### Using CTest or git ifstest

As `git ifstest` calls CTest internally, these two approaches are very similar. 
In the case of CTest, go to the `ifs-test` build directory and simply run 
`ctest`. For `git ifstest`, go to the `ifs-source` main directory and run 
`IGT_TEST_LAUNCHER='salloc -n 128 -q np --mem=20GB' git ifstest -t`.
(Currently, `git ifstest` only reserves 8 cores for the tests which is not
sufficient for the larger tests. Therefore we have to launch `git ifstest` with
a sufficient number of tasks).

Several environment variables control the testing: 
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


## Technical details

The `ifsbench` based tests are controlled by the following scripts:

### ifs-run

The `share/ifs-test/ifs-run.py` script is configured by CMake and 
* stores a number of configuration values (cycle number, used IFS source
  directory, build directory, etc.).
* holds a list of `ifsbench`-based Python runscripts that can be used as
  separate commands using this script.

### from_yaml

The `tests/ifsbench/from_yaml` file is the logical core of the `ifsbench`-based
tests. It parses the previously specified experiment `.yaml` file and actually
runs it. 

### ifsbench_runner

The `share/ifs_test/ifsbench_runner.sh` script is a simple wrapper around
`ifs-run.py`. It does two things:
* Checks the `IFS_TEST_BITIDENTICAL` environment variable and adjusts the
  arguments to `ifs-run.py` accordingly.
* Skips a test if a test was set to be large (`IS_LARGE` in
  `cmake/add_ifsbench_test.cmake`) but `IFS_TEST_RUN_LARGE` is not 1. 
