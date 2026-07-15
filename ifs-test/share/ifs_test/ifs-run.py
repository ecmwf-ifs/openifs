#!@IFSBENCH_Python3_EXECUTABLE@
# (C) Copyright 2011- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction


"""
Run script for ifs-test

This head script allows to execute any of the tests defined in ifs-test.
It is populated with configuration values from the CMake build (such as
IFS version, precision, source and build directories etc.).
"""

from pathlib import Path
from importlib import import_module

from ifsbench import warning, cli

# Set some default values for tests
from ifsbench_runner.config import ifstest_config

ifstest_config['sourcedir'] = Path("@PROJECT_SOURCE_DIR@")
ifstest_config['builddir'] = Path("@IFS_BIN_DIR@").parent
ifstest_config['cachedir'] = Path("@DATA_CACHE_DIR@")
ifstest_config['sharedir'] = Path("@SHARE_SOURCE_DIR@")
ifstest_config['cycle'] = f'cy{"@IFS_CYCLE@".lower()}'
ifstest_config['precision'] = "@IFS_PRECISION@".lower()
ifstest_config['nml_template'] = Path("@IFS_NAMELIST@")
ifstest_config['libblackdir'] = Path("@IFS_BLACKLIST_LIBRARY_PATH@")

# Collect all available ifs_run test groups
ifs_tests = "@IFSBENCH_RUNNER_LIST@".split(";")
for test in ifs_tests:
    try:
        import_module("ifsbench_runner." + test)
    except ImportError:
        warning(f'Failed to import ifsbench_runner.{test}')

if __name__ == '__main__':
    cli(auto_envvar_prefix='IFSBENCH')
