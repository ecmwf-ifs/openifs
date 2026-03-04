LAUNCHER="mpirun"
LAUNCHER_MEM_FLAG=""
LAUNCHER_NPROC_FLAG="-np"
LAUNCHER_NTHREAD_FLAG=""
LAUNCHER_OTHER_FLAGS=""

if [[ "${ECPLATFORM:-"unset"}" == "hpc2020" ]] ; then
    LAUNCHER="srun"
    LAUNCHER_MEM_FLAG="--mem-per-cpu=1G"
    LAUNCHER_NPROC_FLAG="-n"
    LAUNCHER_NTHREAD_FLAG="--cpus-per-task"
    LAUNCHER_OTHER_FLAGS="--gres=ssdtmp:0"
    export OMP_PROC_BIND=true
    export OMP_PLACES=threads
fi
