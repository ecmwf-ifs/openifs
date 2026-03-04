! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

! These functions are to be used within drhook C methods, to avoid having OMP pragmas there.

subroutine drhook_run_omp_parallel_ipfstr(NTIDS, FUNC, CDSTR)
! Usage:
! ------
! void func( const char* string, long strlength ) { ... }
! extern void drhook_run_omp_parallel_ipfstr_(const int *,  void (*func)(const char *, long),
!                                             const char *, /*hidden*/ long);
! drhook_run_omp_parallel_ipfstr_(&ntids,func,string,strlen(string));
use, intrinsic :: iso_c_binding, only : c_char, c_int
implicit none
INTEGER(KIND=C_INT), INTENT(IN) :: NTIDS
EXTERNAL :: FUNC
CHARACTER(LEN=*,KIND=C_CHAR), INTENT(IN) :: CDSTR
!$OMP PARALLEL NUM_THREADS(NTIDS)
CALL FUNC(CDSTR)
!$OMP END PARALLEL
end subroutine drhook_run_omp_parallel_ipfstr

subroutine drhook_run_omp_parallel_ipfipipipdpstr(NTIDS, FUNC, KTIDS, TARGET_OMPTID, TARGET_SIG, START_TIME, CDSTR)
use, intrinsic :: iso_c_binding, only : c_char, c_int, c_double
implicit none
INTEGER(KIND=C_INT), INTENT(IN) :: NTIDS, KTIDS, TARGET_OMPTID, TARGET_SIG
REAL(KIND=C_DOUBLE), INTENT(IN) :: START_TIME
CHARACTER(LEN=*,KIND=C_CHAR), INTENT(IN) :: CDSTR
EXTERNAL :: FUNC
!$OMP PARALLEL NUM_THREADS(NTIDS)
CALL FUNC(KTIDS, TARGET_OMPTID, TARGET_SIG, START_TIME, CDSTR)
!$OMP END PARALLEL
end subroutine drhook_run_omp_parallel_ipfipipipdpstr

subroutine drhook_run_omp_parallel_get_cycles(NTIDS, NCYCLES)
use, intrinsic :: iso_c_binding, only : c_int, c_long_long
use ec_parkind, only : JPIM, JPIB
implicit none
INTEGER(KIND=C_INT), INTENT(IN) :: NTIDS
INTEGER(KIND=C_LONG_LONG), INTENT(INOUT) :: NCYCLES(0:NTIDS-1)
INTEGER(KIND=JPIM) :: IOMPTID
INTEGER(KIND=JPIM) OMP_GET_THREAD_NUM
INTEGER(KIND=C_LONG_LONG), EXTERNAL :: ec_get_cycles ! from ec_get_cycles.c
INTEGER(KIND=C_LONG_LONG) :: ICYCLES
#ifdef _OPENMP
EXTERNAL OMP_GET_THREAD_NUM
#else
OMP_GET_THREAD_NUM() = 0
#endif
!-- Obtain per OpenMP-thread CPU-cycles increment since last call
!$OMP PARALLEL NUM_THREADS(NTIDS) PRIVATE(IOMPTID,ICYCLES) SHARED(NCYCLES)
IOMPTID = OMP_GET_THREAD_NUM()
ICYCLES = ec_get_cycles()
NCYCLES(IOMPTID) = ICYCLES - NCYCLES(IOMPTID)
!$OMP END PARALLEL
end subroutine drhook_run_omp_parallel_get_cycles
