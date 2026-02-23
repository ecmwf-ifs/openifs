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
module_load intel/2023.2.0
module_load impi/2021.10.0
module_load hdf5/1.14.1-2
module_load pnetcdf/1.12.3
module_load netcdf/2023-06-14
module_load fftw/3.3.10
module_load mkl/2024.1
module_load ucx/1.16.0
module_load aec/1.1.2
module_load cmake/3.29.2
module_load python/3.12.1

export TBBMALLOC_DIR="/apps/GPP/ONEAPI/2023.2.0/tbb/2021.10.0/lib/intel64/gcc4.8"
export TBBROOT="/apps/GPP/ONEAPI/2023.2.0/tbb/2021.10.0"
export TBB_MALLOC_USE_HUGE_PAGES=1
export TBB_MALLOC_SET_HUGE_SIZE_THRESHOLD=0
export I_MPI_FABRICS="shm:ofi"
export I_MPI_OFI_PROVIDER="verbs"
export FI_PROVIDER="verbs"
export I_MPI_PLATFORM="spr"
export UCX_TLS="rc,sm,self"   # or rc,self if no shared memory

# Setting required for bit reproducibility with Intel MKL:
export MKL_CBWR=AUTO,STRICT

# Record the RPATH in the executable
export LD_RUN_PATH=$LD_LIBRARY_PATH

# Undo stack size limitation enforced by Python module (prevent segfault during
# runtime)
ulimit -s unlimited

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null
