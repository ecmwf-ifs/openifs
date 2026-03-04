# (C) Copyright 2011- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

#
#
#   exp-config.h
#
#   OpenIFS experiment configuration file
#
#
#   This file is required in the experiment directory for oifs-run.
#   It needs to be read using command:
#
#     source ./exp-config.h
#
#

#--- required variables for this experiment:

OIFS_EXPID="<expid>"    # your experiment ID
OIFS_RES="<res>"        # spectral resolution
OIFS_GRIDTYPE="l"       # 'l'=linear reduced, 'o'=cubic octahedral
OIFS_NPROC=8            # no of MPI tasks
OIFS_NTHREAD=4          # no of OpenMP threads
OIFS_PPROC=false        # enable postprocessing of model output
OUTPUT_ROOT=$(pwd)      # where output folder for pproc is created
LFORCE=true             # overwrite option
LAUNCH=""               # overwrite platform specific run command

#--- optional variables that can be set for this experiment:
#
#OIFS_NAMELIST='my-fort.4'               # custom atmospheric model namelist file
#OIFS_EXEC="<new-location>/ifsMASTER.DP" # model exec to be used for this experiment
#
