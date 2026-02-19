# (C) Copyright 2011- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

# Source me to get the correct configure/build/run environment

# Store tracing and disable (module is *way* too verbose)
{ tracing_=${-//[^x]/}; set +x; } 2>/dev/null

module_load() {
  echo "+ module load $*"
  module load $*
}
module_unload() {
  echo "+ module unload $*"
  module unload $*
}
module_purge() {
  echo "+ module purge"
  module purge
}

# Unload all modules to be certain
[[ ${IFS_RUNTIME_ENV:-unset} == "unset" ]] && module_purge

# Load modules

module_load prgenv/intel
module_load intel/2021.4.0
module_load intel-mpi/2021.4.0
module_load intel-mkl/19.0.5

# Don't load these modules if env.sh is used as part of the IFS runtime environment - only the modules above are required
if [[ ${IFS_RUNTIME_ENV:-unset} == "unset" ]]; then
  module_load python3/3.8.8-01
  module_load fftw/3.3.9
  module_load netcdf4/4.7.4
  module_load hdf5/1.10.6
  module_load eigen/3.3.7 
  module_load cmake/3.20.2
  module_load ninja/1.10.0
  module_load fcm/2019.05.0
  module_load aec/1.0.4
fi

# Setting required for bit reproducibility with Intel MKL:
export MKL_CBWR=AUTO,STRICT

# Record the RPATH in the executable
export LD_RUN_PATH=$LD_LIBRARY_PATH

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

