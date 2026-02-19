# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )
find_package( ecbuild 3.3 REQUIRED )
project( fiat_test_install VERSION 0.0.0 LANGUAGES Fortran )

find_package( fiat REQUIRED )

ecbuild_add_executable( TARGET  main_dp
                        SOURCES main.F90
                        LIBS    fiat parkind_dp )

ecbuild_add_executable( TARGET  main_sp
                        SOURCES main.F90
                        LIBS    fiat parkind_sp )
