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

module_purge
module_load gcc/12.3.0
module_load openmpi/4.1.5-gcc
module_load hdf5/1.14.1-2-gcc-openmpi
module_load pnetcdf/1.12.3-gcc-openmpi
module_load netcdf/c-4.9.2_fortran-4.6.1_cxx4-4.3.1_hdf5-1.14.1-2_pnetcdf-1.12.3-gcc-openmpi
module_load cmake/3.29.2
module_load aec/1.1.2-gcc
module_load mkl/2024.1
module_load python/3.12.1-gcc
module_load ucx/1.16.0-gcc
module_load fftw/3.3.10-gcc-ompi

# Setting required for bit reproducibility with Intel MKL:
export MKL_CBWR=AUTO,STRICT

# Record the RPATH in the executable
export LD_RUN_PATH=$LD_LIBRARY_PATH

# Undo stack size limitation enforced by Python module (prevent segfault during
# runtime)
ulimit -s unlimited

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null
