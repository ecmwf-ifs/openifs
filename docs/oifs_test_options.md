# OpenIFS test options (more detail)

An OpenIFS build can be tested, initially, using the `-t` option with `openifs-test.sh`, e.g. 

```bash
$OIFS_TEST/openifs-test.sh -t
```

`-t` on its own will always run 22 cases, i.e., 21 coarse resolution (t21) 3-D NWP tests with and without chemistry and 1 SCM test (based on TWP-ICE) and will not test bit comparison. Hence, there are also some extra options to enable individual or group tests and bit comparison testing. 

## Running individual or groups of tests

With OpenIFS 48r1, the openifs-test suite includes 22 tests and can take up to 10 mins to run (depending on system). With future cycles the number of tests will increase as the test coverage is improved, so the run time for `openifs-test.sh -t` will increase. 

To speed this up, either a group of tests or one test can be run using the command line option `–R`, e.g.

```bash
$OIFS_TEST/openifs-test.sh -t -R "test_fc”
```

which will run all the tests with `test_fc` in their name

Under the hood the `-R` option, uses a regular expression, so to run an individual case, append the `test_name` with `$`, e.g.

```bash
$OIFS_TEST/openifs-test.sh -t -R "test_fc$”
```

Running a group or individual test can speed up testing, particularly if you are frequently running `openifs-test –t`  during your development workflow. It is important to note, however, that while it is useful and efficient to run individual or groups of tests, it is recommended that all tests are run intermittently when making code changes, so that unexpected/unintended impacts of changes are captured.

## Bit identical testing

Bit identical testing over multiple tests can be a very useful tool to ensure that any science/code changes only impact the results that they are expected to impact, while other results are unaffected.

By default, `openifs-test.sh` as described only checks that an installation and/or code change builds and runs without any major runtime errors; it does not perform bit-comparison testing against a reference set of known-good output (KGO).

In order to perform a bit comparison testing it is first necessary to create the KGO from a "known-good" source of the installation of the code. The creation of the KGO should be done prior to any code changes. KGO can be created by running the test with the environment variable `IFS_TEST_BITIDENTICAL=init`:

```bash
IFS_TEST_BITIDENTICAL=init $OIFS_TEST/openifs-test.sh -t
```

Note that only -t  is used, so the above assumes that OpenIFS has already been built.

When each test that supports bit identical testing runs, the above will filter out all the norms from the NODE file and write the data to a `SAVED_NORMS` file. These `SAVED_NORMS` represent the KGO to which other tests can be compared.

Once changes have been made to the source code, you can recompile and check that the new norms are identical to the previous ones in `SAVED_NORMS`, by setting `IFS_TEST_BITIDENTICAL=check`:

```bash
IFS_TEST_BITIDENTICAL=check $OIFS_TEST/openifs-test.sh -t
```

This time, each supported test will compare the norms in the NODE file with the previously created reference in `SAVED_NORMS`. If they differ then the test will fail. The tests will also fail if there is no reference file because the tests haven't been run in `init` mode.

It can be useful to do this incrementally as you make a series of changes for a single feature that you expect to be bit-identical, to catch any point at which bit-reproducibility is accidentally lost.

It is important to note that running with `--clean` will wipe any existing `SAVED_NORMS`; also note that `--build-type=BIT` and `--build-type=DEBUG` will produce different norms.

Finally, at present there's currently no built-in mechanism for automatically testing that two different tests produce identical results (e.g. with different namelist settings, or different numbers of MPI tasks or OpenMP threads).
