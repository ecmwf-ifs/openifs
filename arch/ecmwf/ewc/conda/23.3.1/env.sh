# (C) Copyright 2011- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

# Empty environment
CC=/opt/conda/envs/ecmwf-lab/bin/mpicc
CXX=/opt/conda/envs/ecmwf-lab/bin/mpicxx
LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}:/opt/conda/envs/ecmwf-lab/lib
