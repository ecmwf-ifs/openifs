#! @IFSBENCH_Python3_EXECUTABLE@
# (C) Copyright 2011- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

from pathlib import Path
import click

from ifsbench import (
    cli, reference_options, run_options,
    Benchmark, DrHook, IFS, symlink_data,
    error
)

try:
    from .config import ifstest_config
except ImportError:
    ifstest_config = {}

# Obtain IFS source and build dirs or provide default values
sourcedir = ifstest_config.get('sourcedir', Path(__file__).resolve().parent.parent.parent)
builddir = ifstest_config.get('builddir', Path.cwd().parent.parent.parent)
cachedir = ifstest_config.get('cachedir', builddir/'.cache')
testdir = sourcedir/'tests/t21'
cycle = ifstest_config.get('cycle', 'default')
prec = ifstest_config.get('precision', 'dp')
nml_template = ifstest_config.get('nml_template', sourcedir/'share/ifs_test/namelist')

# Define the default IFS build and source directory
ifs = IFS.create_cycle(cycle=cycle, prec=prec, sourcedir=sourcedir, builddir=builddir, nml_template=nml_template)

# TODO: Placeholder until we have a more "pyhtonic" way of composing
# preset sets of env flags
default_env = {
    'OMP_STACKSIZE': '64M',
    'CRAYBLAS_AUTOTUNING_OFF': '1',
    'ATP_ENABLED': '1',
    'MPICH_RANK_REORDER_DISPLAY': '1',
    'GRIB_API_IO_BUFFER_SIZE': '4194304',
    'GRIB_API_WRITE_ON_FAIL': '0',
    'GRIB_API_LARGE_CONSTANT_FIELDS': '1',
    'FDB_DEBUG': 'no',
    'FDB_IO': 'aio',
    'TRACEBACK_LEVEL': '2',
    'MPICH_MAX_THREAD_SAFETY': 'multiple',
}


@cli.group()
def t21():
    """
    Benchmark suite to test various model configurations at T21 resolution.
    """


@t21.command('fc', context_settings={"auto_envvar_prefix": "IFSBENCH"})
@run_options
@reference_options
@click.option('--drhook/--no-drhook', default=False, help='Turn DrHook profiling on or off')
def t21_fc(runopts, refopts, drhook):
    """
    Run a T21 forecast benchmark and validate against reference.
    """

    class T21FC(Benchmark):
        """
        Definition of a T21 forecast benchmark
        """

        input_files = ['ICMSHepc8INIT', 'ICMGGepc8INIT', 'ICMGGepc8INIUA']

        namelist = testdir/'t21_fc.nml'

        reference_path = testdir/'t21_fc_reference'

    # Set sub-directory to run benchmark in
    rundir = Path.cwd()/'t21_fc'

    # Enable/disable DrHook based on user flag
    drhook = DrHook.PROF if drhook else DrHook.OFF

    # Create benchmark to perform path resolution and sanity checking
    benchmark = T21FC.from_files(
        ifs=ifs, srcdir=cachedir/'t21/inidata', rundir=rundir, ifsdata=cachedir/'ifsdata'
    )
    benchmark.check_input()  # Ensure all input files are found

    # Execute benchmark with given options and read RunRecord object
    record = benchmark.run(namelist=benchmark.namelist, arch=runopts.arch,
                           nproc=runopts.nproc, nthread=runopts.nthread,
                           env=default_env.copy(), drhook=drhook,
                           fclen=runopts.forecast_length, launch_cmd=runopts.launch_cmd,
                           launch_user_options=runopts.launch_options)

    # Read user-given reference record or load default one
    reference = refopts.path if refopts.path is not None else benchmark.reference_path

    if refopts.update:
        # Store record of this run as the new reference
        record.write(filepath=reference, comment=refopts.comment)

    if refopts.validate:
        # Validate run against reference data and exit(-1) if records don't match
        record.validate(refpath=reference, exit_on_error=True)


@t21.command('compo-fc', context_settings={"auto_envvar_prefix": "IFSBENCH"})
@run_options
@reference_options
@click.option('--drhook/--no-drhook', default=False, help='Turn DrHook profiling on or off')
def t21_compo_fc(runopts, refopts, drhook):
    """
    Run a T21 "Compo" forecast and validate against reference.
    """

    class T21CompoFC(Benchmark):
        """
        Definition of a T21 "Compo" forecast benchmark
        """
        expid = 'hmec'

        input_files = [f'ICMSH{expid}INIT', f'ICMGG{expid}INIT',
                       f'ICMGG{expid}INIUA', f'ICMCL{expid}INIT',
                       'tropo_look_up_cbmhybrid.dat', 'OMI.data.extraterrest',
                       'aerosol_reduce.dat', 'uars_ratio.txt', 'haloe_ch4clim.dat']

        namelist = testdir/'t21_compo_fc.nml'

        reference_path = testdir/'t21_compo_fc_reference'

    # Set sub-directory to run benchmark in
    rundir = Path.cwd()/'t21_compo_fc'

    # Enable/disable DrHook based on user flag
    drhook = DrHook.PROF if drhook else DrHook.OFF

    # Create benchmark to perform path resolution and sanity checking
    srcdir = [cachedir/'t21/inidata', cachedir/'compo_data']

    benchmark = T21CompoFC.from_files(ifs=ifs, srcdir=srcdir, rundir=rundir, ifsdata=cachedir/'ifsdata')
    benchmark.check_input()  # Ensure all input files are found

    # Add a mid-level symlink for chemistry input data (I know!)
    symlink_data(cachedir/'compo_data/21_full', rundir/'21_full')

    # Execute benchmark with given options and read RunRecord object
    record = benchmark.run(namelist=benchmark.namelist, arch=runopts.arch,
                           nproc=runopts.nproc, nthread=runopts.nthread,
                           env=default_env.copy(), drhook=drhook,
                           fclen=runopts.forecast_length, launch_cmd=runopts.launch_cmd,
                           launch_user_options=runopts.launch_options)

    # Read user-given reference record or load default one
    reference = refopts.path if refopts.path is not None else benchmark.reference_path

    if refopts.update:
        # Store record of this run as the new reference
        record.write(filepath=reference, comment=refopts.comment)

    if refopts.validate:
        # Validate run against reference data and exit(-1) if records don't match
        record.validate(refpath=reference, exit_on_error=True)


if __name__ == "__main__":
    cli()
