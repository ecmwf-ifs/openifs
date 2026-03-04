#!/bin/ksh
trap 'set +x; echo "Received signal, aborting ..."; exit 1' 1 2 3 15 
#
# Run the 3 test examples provided in the Dr.Hook distribution
#
# These examples have been tested on ibm_power4 platform
# where memory profiling and Mflop-counters are definitely available.
#
# Change the FC and LIBS first if you want to run it on your platform.
#
# Usage : drhook_ex.ksh [$ARCH]
#
#
# Please note, that these example rely on serial MPI-lib
# i.e. run only using one MPI-task
#
# In order to run with multiple processors, refer to your local
# MPI documentation on how to include MPI-lib proper and how to run MPI-
# applications on your platform. Also, do not forget to remove references
# to lib mpi_serial !
#

set -veu

if [[ $# -eq 1 ]] ; then
  export ARCH=$1
else
  export ARCH=${ARCH:=ibm_power4}
fi

has_mpi=0
NPES=1
ECYOD=""

if [[ "$ARCH" = ibm_power* ]] ; then
#-- IBM Power4
  FC="xlf90_r -g -qsmp=omp -qextname -qsuffix=cpp=F90 -qfree=F90 -I."
  CC="xlc_r -g -qsmp=noauto -I."
#-- The following assumes that drhook.c has been compiled with -DHPM
#-- Requires odbdummy because of ec_bind()
  LIBS="-L. -ldrhook -lmpi_serial -L/usr/pmapi/lib -lpmapi -lessl -lodbdummy"
#-- If not compiled with -DHPM, then supply pm_* dummies
#  LIBS="-L. -ldrhook -lmpi_serial -lodbdummy -lessl"
# Note: -lessl i.e. ESSL-library is needed due to jsort.F
elif [[ "$ARCH" = rs6000 ]] ; then
#-- IBM Power3 (at ECMWF)
  FC="xlf90_r -g -qsmp=omp -qextname -qsuffix=cpp=F90 -qfree=F90 -I."
  CC="xlc_r -g -qsmp=noauto -I."
  LIBS="-L. -ldrhook -lmpi_serial -lessl"
elif [[ "$ARCH" = linux ]] ; then
#-- Linux with Portland Group F90 (for IA32)
  FC="pgf90 -mp -I. -g -Ktrap=fp -O1 -fPIC"
  CC="gcc -g -m32 -ansi -Werror"
  LIBS="-L. -ldrhook -lmpi_serial -Wl,-export-dynamic"
elif [[ "$ARCH" = @(linuxg95|cygwin) ]] ; then
#-- Linux with GNU g95
  FC="g95 -I. -g -fno-second-underscore"
  CC="gcc -g -m32 -ansi -Werror"
  LIBS="-L. -ldrhook -lmpi_serial"
elif [[ "$ARCH" = ia64 ]] ; then
#-- Linux with Intel F90
  FC="ifort -openmp -cpp -I. -g -fpe0 -O1 -assume byterecl -assume cc_omp -traceback \
            -assume cc_omp -assume underscore -lowercase -nomixed_str_len_arg -fPIC"
  CC="gcc -g -ansi -Werror"
  LIBS="-L. -ldrhook -lmpi_serial -Wl,-export-dynamic"
elif [[ "$ARCH" = amd64 ]] ; then
  MPICH_ROOT=${MPICH_ROOT:=/not/available}
  if [[ -d $MPICH_ROOT && -x $MPICH_ROOT/bin/mpirun ]] ; then
    #-- Linux with Portland Group F90 -- for AMD64 with MPICH
    MPFC="$MPICH_ROOT/bin/mpif90 -f90=pgf90 -mp -I. -g -Ktrap=fp -O1 -fPIC -tp amd64"
    MPLIBS="-L. -ldrhook -L$MPICH_ROOT/lib -lmpichf90 -lmpich -ldl -Wl,-export-dynamic"
    NPES=2
    mpirun="./mpirun.linux -np $NPES"
    has_mpi=1
  fi
  #-- Linux with Portland Group F90 -- for AMD64
  FC="pgf90 -mp -I. -g -Ktrap=fp -O1 -fPIC -tp amd64"
  LIBS="-L. -ldrhook -lmpi_serial -Wl,-export-dynamic"
  CC="gcc -g -m64 -ansi -Werror"
elif [[ "$ARCH" = @(*alpha) ]] ; then
#-- ECMWF Dec Alpha (nowadays an outdated toy machine)
  FC="f90 -convert big_endian -O0 -assume byterecl"
  CC="cc -std -g -O0"
  LIBS="-L. -ldrhook -lmpi_serial"
elif [[ "$ARCH" = i86pc ]] ; then
#-- ECMWF Intel/Solaris
  FC="f90 -xfilebyteorder=big8:%all -I. -M."
  CC=${CC:="gcc -g -ansi -Werror"}
  I86PC_GCC=${I86PC_GCC:="-L/opt/sfw/lib/gcc-lib/i386-pc-solaris2.9/2.95.3 -lgcc"}
  #suse 10.0: I86PC_GCC=${I86PC_GCC:="-L/usr/lib/gcc/i586-suse-linux/4.0.2 -lgcc"}
  LIBS="-L. -ldrhook -lmpi_serial $I86PC_GCC"
  export DBXDEBUGGER=${DBXDEBUGGER:=1}
elif [[ "$ARCH" = "necsx" ]] ; then
  FC="./run_fe sxf90 -DNECSX -V -USX -Cvopt -Wf,-pvctl fullmsg vwork=stack,-L fmtlist transform"
  CC="./run_fe sxcc -DNECSX" 
  LIBS="-L. -ldrhook -lmpi_serial -lgen -lcpp -lpthread"
  export DBXDEBUGGER=${DBXDEBUGGER:=1}
elif [[ "$ARCH" = @(cray_xt*) ]] ; then
  FC="ftn -g"
  CC="cc -g"
  LIBS="-L. -ldrhook -ldl -Wl,-export-dynamic"
  #-- yod stuff
  MPFC="$FC"
  MPLIBS="$LIBS"
  NPES=2
  mpirun="yod -sz $NPES"
  has_mpi=1
  ECYOD="yod -sz 1 -SN"
else
  echo "***Error: Unrecognized ARCH=$ARCH"
  echo "          Please edit this file $0"
  exit 1
fi

if [[ -f ./ODB_FCLIBS ]] ; then
  ODB_FCLIBS=$(cat ./ODB_FCLIBS)
else
  ODB_FCLIBS=""
fi

#-- In case of shareable linkage ...

LD_LIBRARY_PATH=${LD_LIBRARY_PATH:=.}
if [[ "$LD_LIBRARY_PATH" != "." ]] ; then
  LD_LIBRARY_PATH=".:${LD_LIBRARY_PATH}"
fi
export LD_LIBRARY_PATH

#-- Initialize (outside this script) the following to 0/1 to prevent/allow running GNU-debugger
# export GNUDEBUGGER=${GNUDEBUGGER:=0}
#-- Initialize (outside this script) the following to 0/1 to prevent/allow running DBX-debugger
# export DBXDEBUGGER=${DBXDEBUGGER:=0}

set -xv

#-------------------------------------------------------------------------
# Example#0 : Show wordsizes in F90 & C of this arch (these are MPI/serial)
#-------------------------------------------------------------------------

$ECYOD ./sizeof.x || :

$ECYOD ./kind.x || :

#-------------------------------------------------------------------------
# Example#1 : Generate Dr.Hook, run and fail in divide by zero
#-------------------------------------------------------------------------

./insert_drhook drhook_ex1.F90
$FC -I. _drhook_ex1.F90 -o drhook_ex1.x $LIBS
env DR_HOOK=1 $ECYOD ./drhook_ex1.x || :

#-------------------------------------------------------------------------
# Example#2 : With the previous case fixed experience with watch point
#             Testing also mpl_arg_mod
#-------------------------------------------------------------------------

./insert_drhook drhook_ex2.F90
$FC -I. _drhook_ex2.F90 -o drhook_ex2.x $LIBS
env DR_HOOK=1 $ECYOD ./drhook_ex2.x 1000 || :

#-------------------------------------------------------------------------
# Example#3 : Experience with different profilings. Also silence Dr.Hook!
#-------------------------------------------------------------------------

./insert_drhook drhook_ex3.F90
$FC -I. _drhook_ex3.F90 -o drhook_ex3.x $LIBS
export DR_HOOK=1

export DR_HOOK_SILENT=1
export DR_HOOK_SHOW_PROCESS_OPTIONS=0

#-- Wall clock profile
env DR_HOOK_OPT=wallprof $ECYOD ./drhook_ex3.x || :
mv drhook.prof.1 drhook.prof.1.wallprof
cat drhook.prof.1.wallprof

#-- CPU-time profile
env DR_HOOK_OPT=cpuprof $ECYOD ./drhook_ex3.x || :
mv drhook.prof.1 drhook.prof.1.cpuprof
cat drhook.prof.1.cpuprof

#-- Mflop-counter profile
env DR_HOOK_OPT=hpmprof $ECYOD ./drhook_ex3.x || :
mv drhook.prof.1 drhook.prof.1.hpmprof
cat drhook.prof.1.hpmprof

#-- Wall clock & memory usage with timeline
env DR_HOOK_OPT=wall,mem DR_HOOK_TIMELINE=-1 $ECYOD ./drhook_ex3.x > drhook.timeline 2>&1 || :
ls -l drhook.timeline

#-- Memory profile (only)
env DR_HOOK_OPT=memprof $ECYOD ./drhook_ex3.x || :
mv drhook.prof.1-mem drhook.prof.1.memprof
cat drhook.prof.1.memprof

unset DR_HOOK
unset DR_HOOK_SILENT
unset DR_HOOK_SHOW_PROCESS_OPTIONS

#-------------------------------------------------------------------------
# Example#4 : Check whether your Dr.Hook works with C-main program
#-------------------------------------------------------------------------

rc4=0
$CC -c drhook_ex4.c || {
  set +xv
  echo "***Error: Your Dr.Hook C-main program test did not compile"
  rc4=1
}

if [[ $rc4 -eq 0 ]] ; then
  $CC drhook_ex4.o -o ./drhook_ex4.x $LIBS $ODB_FCLIBS || {
    set +xv
    echo "***Error: Your Dr.Hook C-main program test did not link"
    rc4=2
  }
fi

if [[ $rc4 -eq 0 ]] ; then
  #-- Test 4a) : Trigger SIGFPE and fail
  env DR_HOOK=1 $ECYOD ./drhook_ex4.x SIGFPE "ARG TEST" 'test arg' `pwd` || :
fi

if [[ $rc4 -eq 0 ]] ; then
  #-- Test 4b) : Do not fail and produce wall clock profile
  export DR_HOOK_SILENT=1
  export DR_HOOK_SHOW_PROCESS_OPTIONS=0
  export DR_HOOK_PROFILE=cdrhook.prof.%d
  env DR_HOOK=1 DR_HOOK_OPT=wallprof $ECYOD ./drhook_ex4.x || :
  cat cdrhook.prof.1
  unset DR_HOOK_PROFILE
  unset DR_HOOK_SILENT
  unset DR_HOOK_SHOW_PROCESS_OPTIONS
fi

set -xv

#-------------------------------------------------------------------------
# Example#5 : Check MPI implementation, if applicable
#-------------------------------------------------------------------------

rc5=0

export DR_HOOK=1
export DR_HOOK_SILENT=1
export DR_HOOK_SHOW_PROCESS_OPTIONS=0

if [[ $has_mpi -eq 1 ]] ; then
  # "*** Testing parallel MPI-implementation ***"

  ./insert_drhook drhook_ex5.F90
  $MPFC -I. _drhook_ex5.F90 -o drhook_ex5.x $MPLIBS || rc5=$?

  if [[ $rc5 -eq 0 ]] ; then
    #-- Wall clock profile
    env DR_HOOK_OPT=wallprof $mpirun ./drhook_ex5.x one two three || rc5=$?

    if [[ ! -f drhook.prof.1 ]] ; then
    #-- Something went wrong ... reverting to serial version
      rc5=1
      has_mpi=0
    fi

    if [[ $rc5 -eq 0 ]] ; then
      typeset np=0
      while [[ $np -lt $NPES ]]
      do
        ((np+=1))
        mv drhook.prof.$np drhook.prof.$np.wallprof
        cat drhook.prof.$np.wallprof
      done
    fi
  fi

fi

if [[ $has_mpi -eq 0 ]] ; then
  rc5=0
  # "*** Testing MPI-implementation with a dummy MPI-serial version ***"

  ./insert_drhook drhook_ex5.F90
  $FC -I. _drhook_ex5.F90 -o drhook_ex5.x $LIBS || rc5=$?

  if [[ $rc5 -eq 0 ]] ; then
    #-- Wall clock profile
    env DR_HOOK_OPT=wallprof $ECYOD ./drhook_ex5.x one two three || rc5=$?

    if [[ $rc5 -eq 0 ]] ; then
      mv drhook.prof.1 drhook.prof.1.wallprof
      cat drhook.prof.1.wallprof
    fi
  fi

fi

unset DR_HOOK_PROFILE
unset DR_HOOK_SILENT
unset DR_HOOK_SHOW_PROCESS_OPTIONS

#-------------------------------------------------------------------------

((rc=rc4+rc5))

exit $rc
