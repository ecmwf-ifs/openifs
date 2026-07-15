#! @IFSBENCH_Python3_EXECUTABLE@
# (C) Copyright 2011- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction


"""
Run script for ifsbench-based test cases.

This runner takes an experiment yaml file as input which contains the paths
to the experiment data as well as some other default values (see for example
tests/cy48/README.md).
To just run a test, simply pass the experiment yaml file:
```
python3 ifsbench-test.py /path/to/experiment
```

This runs the test in a temporary folder without any additional checks. To
run the test in a given folder (which also keeps all files in place), use
the --run-dir flag:
```
python3 ifsbench-test.py /path/to/experiment --run-dir=/some/folder
```

The default number of processes/threads/processes per node can be set in the
experiment yaml file. They can also be overridden by the corresponding command
line arguments:
```
python3 ifsbench-test.py /path/to/experiment --nproc=256 --nthread=2
```

To use reference data, the path to this reference data must be provided using
the `-r` flag. Then one can use either the `--update-reference` flag to
create reference data at the path given by `-r` or `--validate` to compare
against the reference data:
```
# Create reference data.
python3 ifsbench-test.py /path/to/experiment -r /reference/path --update-reference

# Compare against reference data.
python3 ifsbench-test.py /path/to/experiment -r /reference/path --validate
```
"""

import shutil
import sys
import tempfile
from contextlib import nullcontext
from pathlib import Path

import numpy
import yaml

import click

from ifsbench import (
    reference_options, run_options, ExperimentFiles,
    ExperimentFilesBenchmark, DrHook, IFS, SpecialRelativePath,
    IFSNamelist, cli, RunRecord, debug, error
)

try:
    from .config import ifstest_config
except ImportError:
    ifstest_config = {}

# Obtain IFS source and build dirs or provide default values
_SOURCE_DEFAULT = ifstest_config.get('sourcedir', Path(__file__).resolve().parent.parent.parent)
_BUILD_DEFAULT = ifstest_config.get('builddir', Path.cwd().parent.parent.parent)
_CYCLE_DEFAULT = ifstest_config.get('cycle', 'default')
_PREC_DEFAULT = ifstest_config.get('precision', 'dp')

def _parse_yaml(ctx, param, experiment):
    """
    Callback function for click. It opens the experiment yaml file, extracts
    the parameters in run_opts and sets them as the default values for the
    arguments in the @run_options wrapper.
    """
    with open(experiment, 'r') as f:
        yaml_data = yaml.safe_load(f)

        if 'run_opts' in yaml_data:
            # Add the defaults to the ctx.default_map map (which may be None).
            if ctx.default_map is None:
                ctx.default_map = dict(yaml_data['run_opts'])
            else:
                ctx.default_map.update(dict(yaml_data['run_opts']))


    # We have to return the actual value of the 'experiment' argument,
    # otherwise its value would be None in the main function.
    return experiment

def _compare_records(record, ref_record, tolerance):
    """
    Check that the results in two given RunRecords objects are equal up to a
    given tolerance. If this is the case, return True, otherwise False.
    """
    diff_spectral = record.spectral_norms - ref_record.spectral_norms
    diff_grid = record.gridpoint_norms - ref_record.gridpoint_norms

    success = True

    for frame in (diff_spectral, diff_grid):
        # Do the comparison column wise - this makes it easier to print debug
        # information.
        for column, series in frame.items():

            # Convert a pandas Series (essentially a column of a table) to a
            # numpy array in order to use all the available numpy tools.
            diff = numpy.array(series)

            diff_norm = numpy.linalg.norm(diff, numpy.inf)
            # Use the maximum norm to essentially check that |value| < tol
            # holds for all entries in the column.
            if diff_norm > tolerance:
                error(f"Tolerance is not satisfied in column {column} "
                      f"({diff_norm} > {tolerance}).")
                error(f"Difference is {diff}.")
                success = False
            else:
                debug(f"Tolerance is satisfied in column {column} "
                      f"({diff_norm} <= {tolerance}).")
                debug(f"Difference is {diff}.")

    return success

# Some click-magic is going on here... click will call the callback function
# that is specified in the 'experiment' argument, extract the default run
# options from this experiment file and use them as the default values for
# the argument handling inside the `run_options` wrapper.
@cli.command('from_yaml', context_settings={"auto_envvar_prefix": "IFSBENCH"})
@run_options
@reference_options
@click.argument(
    'experiment',
    type=click.Path(exists=True),
    callback=_parse_yaml,
)
@click.option('--run-dir', type=click.Path(), default=None,
              help='Run directory for the tests (temporary directory by default)')
@click.option('--tasks-per-node', default=1,
              help='Number of MPI processes per node')
@click.option('--drhook/--no-drhook', default=True,
              help='Turn DrHook profiling on or off')
@click.option('--cycle', default=_CYCLE_DEFAULT, type=str,
              help='The IFS cycle that is used')
@click.option('--build-dir', type=click.Path(), default=_BUILD_DEFAULT,
              help='Path to the IFS build directory')
@click.option('--source-dir', type=click.Path(), default=_SOURCE_DEFAULT,
              help='Path to the IFS source directory')
@click.option('--single-precision/--double-precision', type=bool,
              default=(_PREC_DEFAULT == 'sp'),
              help='Use single-precision/double-precision')
@click.option('--tolerance', type=float, default=0.0,
              help="If larger than zero, run the comparison with the given tolerance.")
@click.option('--custom-entry', multiple=True, default=[],
              help="Set namelist entry. Format KEY:VALUE or KEY. Nested keys "
              "should be separated by /. If only KEY is given, this entry is deleted.")
@click.option('--gpu/--no-gpu', default=False,
              help='Use GPUs on/off')
def from_yaml(runopts, refopts, experiment, run_dir, tasks_per_node, drhook,
              cycle, build_dir, source_dir, single_precision, tolerance,
              custom_entry, gpu):
    """
    Execute a ifsbench-based test.
    """
    source_dir = source_dir.resolve()
    build_dir = build_dir.resolve()

    # Additional environment variables that are used during the test execution.
    env = {}
    paths = []

    # List of namelist entries that will be modified.
    namelist_mod = {}

    # Load the yaml file that holds the test setup.
    with open(experiment, 'r') as f:
        exp_data = yaml.safe_load(f)

        # Load the path to the inidata directory or the inidata ifsbench YAML
        # file.
        exp_path = exp_data['exp']['local_path']
        exp_path = Path(exp_path)

        # If a directory is given, look inside it for the ifsbench file summary
        # YAML file.
        if exp_path.is_dir():
            # Find the experiment yaml file in this directory.
            # Match everything that ends with .yml or .yaml.
            exp_path = list(exp_path.glob('*.y*ml'))

            if len(exp_path) != 1:
                raise RuntimeError(f'Found {len(exp_path)} experiment configurations, instead of one!')
            else:
                exp_path = exp_path[0]

        if 'env' in exp_data:
            env = dict(exp_data['env'])

        if 'relative_paths' in exp_data:
            # Convert the entries in the relative_paths list to
            # ifsbench.SpecialRelativePath objects.
            for entry in exp_data['relative_paths']:
                inp = entry['in']
                out = entry['out']
                match_str = entry['match'].strip().lower()

                if match_str == 'left':
                    match = SpecialRelativePath.NameMatch.LEFT_ALIGNED
                elif match_str == 'right':
                    match = SpecialRelativePath.NameMatch.RIGHT_ALIGNED
                elif match_str == 'exact':
                    match = SpecialRelativePath.NameMatch.EXACT
                else:
                    match = SpecialRelativePath.NameMatch.FREE

                paths.append(SpecialRelativePath.from_filename(
                    inp, out, match=match))

        if 'namelist' in exp_data:
            namelist_mod = dict(exp_data['namelist'])

    # Loop over the custom namelist entries and split them into key/value
    # pairs.
    for e in custom_entry:
        tmp = e.split(':', 1)

        if len(tmp) == 1:
            # Only key is given - this entry should be deleted.
            namelist_mod[tmp[0]] = None
        elif len(tmp) == 2:
            # The values arrive as strings and should be converted to the
            # proper datatype before being added to the namelist. Here we just
            # use the yaml parsing ability to convert a string to the
            # appropriate data type.
            value = yaml.safe_load(tmp[1])
            namelist_mod[tmp[0]] = value


    class IFSBenchTest(ExperimentFilesBenchmark):
        """
        Class for the IFSBench test cases. It adds a few default path
        conversions and also includes test-specific path conversions that
        are read from the experiment yaml file.
        """

        # Some default conversions of relative paths.
        special_paths = [
            SpecialRelativePath.from_filename(
                r'wam_grid_tables|wam_subgrid_\d', r'\g<match>',
                match=SpecialRelativePath.NameMatch.LEFT_ALIGNED),

            # We have to modify the wam_namelist later on, therefore we rename
            # the original wam_namelist here. Later, we will create a modified
            # namelist and call it `wam_namelist`.
            SpecialRelativePath.from_filename(
                r'wam_namelist', r'wam_namelist_template',
                match=SpecialRelativePath.NameMatch.EXACT),
            SpecialRelativePath.from_filename(
                r'wam_', r'\g<post>', match=SpecialRelativePath.NameMatch.LEFT_ALIGNED),
            SpecialRelativePath.from_filename(
                r'rtablel_\d+', r'ifs/\g<name>', match=SpecialRelativePath.NameMatch.EXACT),
            SpecialRelativePath.from_dirname(
                'ifsdata', r'ifsdata\g<child>', match=SpecialRelativePath.NameMatch.EXACT),

            # Rename the fort.4 namelist to something else as fort.4 is
            # overwritten by the tests.
            SpecialRelativePath.from_filename(
                'fort.4', r'namelist_template', match=SpecialRelativePath.NameMatch.EXACT),
        ]

        # Add the test-specific path conversions.
        special_paths += paths

    if run_dir is None:
        # If no run_dir is given, create a temporary directory for these runs.
        run_dir_manager = tempfile.TemporaryDirectory(dir=Path.cwd())
    else:
        # If a run directory is given, create a dummy context around run_dir.
        # Here, we use the absolute path (resolve) as using relative paths may
        # cause issues when we pass the path to IFS (as it may run in a
        # different directory.
        run_dir_manager = nullcontext(Path(run_dir).resolve())

    precision = 'sp' if single_precision else 'dp'

    # Define the default IFS build and source directory
    ifs = IFS.create_cycle(
        cycle=cycle, prec=precision,
        sourcedir=source_dir, builddir=build_dir
    )

    # Enable/disable DrHook.
    drhook = DrHook.PROF if drhook else DrHook.OFF

    # Use a context manager for the run directory. This will clean up all
    # stuff, if we used a temporary directory.
    with run_dir_manager as run_dir:
        run_dir = Path(run_dir).resolve()
        untar_dir = Path(run_dir)/'untar'
        untar_dir.mkdir(parents=True, exist_ok=True)

        # Create benchmark to perform path resolution and sanity checking
        exp_files = ExperimentFiles.from_tarball(
            summary_file=exp_path, 
            input_dir=exp_path.parent,
            output_dir=untar_dir, 
            with_ifsdata=True,
            verify_checksum=True
        )

        benchmark = IFSBenchTest.from_experiment_files(exp_files=exp_files,
            ifs=ifs, rundir=run_dir)
        benchmark.check_input()  # Ensure all input files are found

        # The wam_namelist file contains a CPATH variable that by default points
        # to the folder of the original experiment directory. As this directory
        # may not exist anymore (or is not readable), we change it to the run
        # directory.
        wam_nlst = IFSNamelist(template=run_dir/'wam_namelist_template')
        wam_nlst['NALINE']['CPATH'] = str(run_dir)
        wam_nlst.write(run_dir/'wam_namelist', force=True)

        # Take the namelist template, modify some entries (if requested) and
        # save it as namelist_in.
        namelist = IFSNamelist(template=run_dir/'namelist_template')

        # Update the namelist entries, using the data in namelist_mode.
        # If the given value is None, delete the given entry.

        for key, value in namelist_mod.items():
            # Split the key into its different components/levels (or whatever
            # it is called in nested Fortran namelists).
            subkeys = key.split('/')

            tmp = namelist
            for subkey in subkeys[:-1]:
                if subkey not in tmp:
                    tmp[subkey] = {}

                tmp = tmp[subkey]

            if value is None:
                if subkeys[-1] in tmp:
                    del tmp[subkeys[-1]]
            else:
                tmp[subkeys[-1]] = value


        namelist.write(run_dir/'namelist_in')

        gpus_per_task = 1 if gpu else 0

        # Execute benchmark with given options and read RunRecord object
        record = benchmark.run(
            namelist=run_dir/'namelist_in', arch=runopts.arch,
            env=env, drhook=drhook,
            nproc=runopts.nproc, tasks_per_node=tasks_per_node, nthread=runopts.nthread,
            nproc_io=runopts.nproc_io, fclen=runopts.forecast_length,
            launch_cmd=runopts.launch_cmd, launch_user_options=runopts.launch_options,
            logfile=run_dir/'ifs.log', gpus_per_task=gpus_per_task
        )

        # Do some reference handling if a reference path was given.
        if refopts.path is not None:
            # Read user-given reference record or load default one
            reference = Path(refopts.path).resolve()

            if refopts.update:
                # Store record of this run as the new reference
                record.write(filepath=reference, comment=refopts.comment)

            if refopts.validate:
                if tolerance <= 0:
                    record.validate(refpath=reference, exit_on_error=True)
                else:
                    ref_record = RunRecord.from_json(reference)
                    success = _compare_records(record, ref_record, tolerance)
                    if not success:
                        sys.exit(1)
                    else:
                        sys.exit(0)


@cli.command('save_reference', context_settings={"auto_envvar_prefix": "IFSBENCH"})
@click.argument('node_file', type=click.Path(exists=True))
@click.argument('output_file', type=click.Path())
def save_reference(node_file, output_file):
    """
    Take a NODE file and convert it to an ifsbench-style JSON result file.
    """
    record = RunRecord.from_run(node_file)
    record.write(output_file)


@cli.command('validate_reference', context_settings={"auto_envvar_prefix": "IFSBENCH"})

@click.argument('node_file', type=click.Path(exists=True))
@click.argument('ref_file', type=click.Path(exists=True))
@click.option('--tolerance', type=float, default=0.0,
    help="If larger than zero, run the comparison with the given tolerance.")
def validate_reference(node_file, ref_file, tolerance):
    """
    Compare the results that are stored in a NODE file against an ifsbench-style
    JSON result file.

    :param tolerance: If <= 0, check the results for bit-identicality.
                      If > 0, use some tolerance.
    """

    node_file = Path(node_file).resolve()
    ref_file = Path(ref_file).resolve()

    record = RunRecord.from_run(node_file)

    if tolerance <= 0.0:
        record.validate(ref_file, exit_on_error=True)
    else:
        record = RunRecord.from_run(node_file)
        ref_record = RunRecord.from_json(ref_file)

        success = _compare_records(record, ref_record, tolerance)
        if not success:
            sys.exit(1)
        else:
            sys.exit(0)

if __name__ == "__main__":
    cli()
