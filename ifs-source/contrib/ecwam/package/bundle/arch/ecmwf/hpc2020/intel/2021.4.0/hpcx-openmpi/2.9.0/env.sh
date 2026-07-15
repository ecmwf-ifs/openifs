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
module_purge

# Load modules
module_load prgenv/intel
module_load intel/2021.4.0
module_load hpcx-openmpi/2.9.0
module_load intel-mkl/19.0.5
module_load fftw/3.3.9
module_load netcdf4/4.7.4
module_load hdf5/1.10.6
module_load eigen/3.3.7 
module_load cmake/3.25.2
module_load ninja/1.10.0
module_load fcm/2019.05.0
module_load aec/1.0.4

# Setting required for bit reproducibility with Intel MKL:
export MKL_CBWR=AUTO,STRICT

# Record the RPATH in the executable
export LD_RUN_PATH=$LD_LIBRARY_PATH

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

