! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE ECSORT_MIX
USE EC_PARKIND   , ONLY : JPIM     ,JPIB,    JPRM    ,JPRD
USE YOMHOOK    , ONLY : LHOOK, DR_HOOK, JPHOOK
USE OML_MOD     , ONLY : OML_MAX_THREADS, OML_MY_THREAD, OML_IN_PARALLEL
USE STRHANDLER_MOD, ONLY : TOUPPER
USE MPL_MODULE , ONLY : MPL_MYRANK, MPL_NPROC

!$
!..   Author: Sami Saarinen, ECMWF, 10/02/98
!     Fixes : Sami Saarinen, ECMWF, 08/11/99 : Sub-arrays go now correctly (look for addrdiff)
!                                              Genuine real(4) sort "re-habilitated"
!                                              sizeof_int, _real4 & _real8 HARDCODED !
!             Sami Saarinen, ECMWF, 11/10/00 : REAL*4 version included (REAL_M)
!             Sami Saarinen, ECMWF, 28/11/03 : Calls to DR_HOOK added manually (on top of CY28)
!             Sami Saarinen, ECMWF, 18/02/05 : 64-bit integer sorting introduced (for CY30)
!             Sami Saarinen, ECMWF, 22/02/05 : Using genuine 64-bit rsort64() => one-pass through data
!             Sami Saarinen, ECMWF, 06/07/05 : "current_method" made OpenMP-thread aware (for max. # of threads = NTHRDS)
!             Sami Saarinen, ECMWF, 07/07/05 : Quick-sort method finally arrived (and applicable to multikeys, too)
!                                              Quick-sort the default for scalar machines ("non-VPP"), VPPs is radix-sort
!             Sami Saarinen, ECMWF, 03/07/07 : Quick-sort method uses stable approach (was not guaranteed so before)
!             Sami Saarinen, ECMWF, 15/10/07 : Subroutines put into a common file ../include/ecsort_shared.h and
!                                              preprocessed from there
!             Sami Saarinen, ECMWF, 15/10/07 : NTHRDS increased from 32 to 64
!             Sami Saarinen, ECMWF, 16/10/07 : The default sorting method can be overriden via export EC_SORTING_METHOD=[<number>|<string>]
!             Sami Saarinen, ECMWF, 30/10/07 : Support for CountingSort added as part of QuickSort speedup
!             Sami Saarinen, ECMWF, 31/10/07 : CALL SORTING_METHOD() now prints the prevailing method from EC_SORTING_METHOD
!             Sami Saarinen, ECMWF, 01/11/07 : CountingSort implemented as independent method
!             Sami Saarinen, ECMWF, 06/11/07 : index_adj added as an optional argument to init_index ; new routine adjust_index()
!             Sami Saarinen, ECMWF, 07/11/07 : threshold length "nomp" (see below). Override with EC_SORTING_NOMP.
!             Sami Saarinen, ECMWF, 12/11/07 : Gnome-sort -- the easiest sort on Earth (and very slow for large arrays)
!             Sami Saarinen, ECMWF, 15/11/07 : OpenMP-sorting still under development. Do NOT override the EC_SORTING_NOMP yet.
!             Sami Saarinen, ECMWF, 05/12/07 : When export EC_SORTING_INFO=0, then no info messages are printed from CALL sorting_method()
!             Sami Saarinen, ECMWF, 20/12/07 : export EC_SORTING_INFO=0, rather than =1 is now the default --> less hassling output


IMPLICIT NONE
SAVE
PRIVATE

INTEGER(KIND=JPIM), PARAMETER :: NTHRDS = 64 ! ***Note: A hardcoded max number of threads !!!

INTEGER(KIND=JPIM), PARAMETER :: SIZEOF_INT4  = 4
INTEGER(KIND=JPIM), PARAMETER :: SIZEOF_INT8  = 8
INTEGER(KIND=JPIM), PARAMETER :: SIZEOF_REAL4 = 4
INTEGER(KIND=JPIM), PARAMETER :: SIZEOF_REAL8 = 8

INTEGER(KIND=JPIM), PARAMETER :: MIN_METHOD = 1
INTEGER(KIND=JPIM), PARAMETER :: MAX_METHOD = 5

INTEGER(KIND=JPIM), PARAMETER :: RADIXSORT_METHOD = 1
INTEGER(KIND=JPIM), PARAMETER :: HEAPSORT_METHOD  = 2
INTEGER(KIND=JPIM), PARAMETER :: QUICKSORT_METHOD = 3
INTEGER(KIND=JPIM), PARAMETER :: COUNTINGSORT_METHOD = 4
INTEGER(KIND=JPIM), PARAMETER :: GNOMESORT_METHOD = 5

CHARACTER(LEN=12), PARAMETER :: METHOD_NAME(MIN_METHOD:MAX_METHOD) = &
     & (/&
     &   'RADIXSORT   ' &
     &  ,'HEAPSORT    ' &
     &  ,'QUICKSORT   ' &
     &  ,'COUNTINGSORT' &
     &  ,'GNOMESORT   ' &
     &  /)

!-- Select such method for default_method, which also works for multikey sorts
!   Vector machines should choose radixsort_method, others quicksort_method (oh, sorry, countingsort_method!)
!
!   Note: Occasionally radixsort_method may be faster on non-vector machines, too
INTEGER(KIND=JPIM) :: DEFAULT_METHOD = COUNTINGSORT_METHOD
INTEGER(KIND=JPIM) :: CURRENT_METHOD(NTHRDS) = COUNTINGSORT_METHOD

!-- A threshold length after which OpenMP in sorting, merging, copying may kick in.
!   Override with EC_SORTING_NOMP. Detected while initializing/calling SORTING_METHOD
!   Non-positive values (<= 0) indicate that OpenMP will NOT be attempted at all
INTEGER(KIND=JPIM) :: NOMP = -1

!-- EC_SORTING_INFO = MPL-proc id [1..$NPES] to print the info message when CALL sorting_method()
!                   : 0 (no print; the default)
!                   : 1 (print on MPL-task one)
!                   : > 1 and <= $NPES (some other MPL-task than 1 prints)
INTEGER(KIND=JPIM) :: NSINFO = 0

INTERFACE KEYSORT
MODULE PROCEDURE &
     &INT4_KEYSORT_1D, INT4_KEYSORT_2D, &
     &INT8_KEYSORT_1D, INT8_KEYSORT_2D, &
     &REAL8_KEYSORT_1D, REAL8_KEYSORT_2D, &
     &REAL4_KEYSORT_1D, REAL4_KEYSORT_2D
END INTERFACE

INTERFACE SORTING_METHOD
MODULE PROCEDURE INT_SORTING_METHOD, STR_SORTING_METHOD
END INTERFACE

PUBLIC :: KEYSORT
PUBLIC :: INIT_INDEX, GET_RANK, ADJUST_INDEX
PUBLIC :: SORTING_METHOD
PUBLIC :: RADIXSORT_METHOD, HEAPSORT_METHOD, QUICKSORT_METHOD
PUBLIC :: COUNTINGSORT_METHOD, GNOMESORT_METHOD

CONTAINS

!----------------------------
!--   Public subroutines   --
!----------------------------

SUBROUTINE INT_SORTING_METHOD(INEW, IOLD)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)  :: INEW
INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: IOLD
INTEGER(KIND=JPIM) :: ITMP, IMYPROC, INPES
INTEGER(KIND=JPIM) :: ITID
CHARACTER(LEN=20) CLENV
LOGICAL, SAVE :: LLFIRST = .TRUE.
LOGICAL LLOMP, LLHOOK_OK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-- This maybe called from the very first call of DR_HOOK_UTIL ...
LLHOOK_OK = LHOOK .AND. (PRESENT(INEW) .OR. PRESENT(IOLD))
IF (LLHOOK_OK) CALL DR_HOOK('ECSORT_MIX:INT_SORTING_METHOD',0,ZHOOK_HANDLE)
ITID = OML_MY_THREAD()
LLOMP = OML_IN_PARALLEL()
IF (PRESENT(IOLD)) IOLD = CURRENT_METHOD(ITID)
ITMP = -1
IF (PRESENT(INEW)) THEN
  ITMP = INEW
ELSE IF (.NOT.LLOMP) THEN ! Override the default method (only if outside the OpenMP)
  ITMP = -1 ! no change
  IF (LLFIRST) THEN ! Do once per execution only
    INPES = MPL_NPROC()
    CALL GET_ENVIRONMENT_VARIABLE('EC_SORTING_INFO',CLENV) ! ../support/env.c
    IF (CLENV /= ' ') THEN
      ITMP = NSINFO
      READ(CLENV,'(i20)',END=89,ERR=89) ITMP
      GOTO 88
89    CONTINUE
      ITMP = NSINFO ! no change
88    CONTINUE
      IF (ITMP <= 0) THEN
        NSINFO = 0
      ELSE IF (ITMP >= 1 .AND. ITMP <= INPES) THEN
        NSINFO = ITMP
      ENDIF
    ENDIF
    IMYPROC = MPL_MYRANK()
    CALL GET_ENVIRONMENT_VARIABLE('EC_SORTING_METHOD',CLENV) ! ../support/env.c
    IF (CLENV /= ' ') THEN
       IF (IMYPROC == NSINFO) WRITE(0,'(a)')'<EC_SORTING_METHOD='//TRIM(CLENV)
       CALL TOUPPER(CLENV)
       SELECT CASE (CLENV)
       CASE ('RADIX', 'RADIXSORT')
          ITMP = RADIXSORT_METHOD
       CASE ('HEAP', 'HEAPSORT')
          ITMP = HEAPSORT_METHOD
       CASE ('QUICK', 'QUICKSORT', 'QSORT')
          ITMP = QUICKSORT_METHOD
       CASE ('COUNT', 'COUNTINGSORT', 'COUNTSORT')
          ITMP = COUNTINGSORT_METHOD
       CASE ('GNOME', 'GNOMESORT')
          ITMP = GNOMESORT_METHOD
       CASE ('DEFAULT', 'DEF')
          ITMP = -1 ! no change
       CASE DEFAULT
          READ(CLENV,'(i20)',END=99,ERR=99) ITMP
          GOTO 98
99        CONTINUE
          ITMP = -1 ! no change
98        CONTINUE
       END SELECT
    ENDIF

    IF (ITMP < MIN_METHOD .OR. ITMP > MAX_METHOD) ITMP = DEFAULT_METHOD
    DEFAULT_METHOD = ITMP
    CURRENT_METHOD(:) = ITMP
    IF (IMYPROC == NSINFO) &
         & WRITE(0,'(a,i1,a)')'>EC_SORTING_METHOD=',DEFAULT_METHOD,&
         & ' # '//METHOD_NAME(DEFAULT_METHOD)

    CALL GET_ENVIRONMENT_VARIABLE('EC_SORTING_NOMP',CLENV)
    IF (CLENV /= ' ') THEN
       IF (IMYPROC == NSINFO) WRITE(0,'(a)')'<EC_SORTING_NOMP='//TRIM(CLENV)
       READ(CLENV,'(i20)',END=199,ERR=199) ITMP
       GOTO 198
199    CONTINUE
       ITMP = NOMP ! no change
198    CONTINUE
       NOMP = ITMP
    ENDIF

    IF (IMYPROC == NSINFO) THEN
      WRITE(CLENV,'(i20)') NOMP
      WRITE(0,'(a)')'>EC_SORTING_NOMP='//TRIM(ADJUSTL(CLENV))
    ENDIF

    ITMP = DEFAULT_METHOD
    LLFIRST = .FALSE.
  ENDIF
ENDIF

IF (ITMP < MIN_METHOD .OR. ITMP > MAX_METHOD) ITMP = DEFAULT_METHOD
IF (LLOMP) THEN ! Only this thread sees the change
  CURRENT_METHOD(ITID) = ITMP
ELSE ! All threads see the change
  CURRENT_METHOD(:) = ITMP
ENDIF
IF (LLHOOK_OK) CALL DR_HOOK('ECSORT_MIX:INT_SORTING_METHOD',1,ZHOOK_HANDLE)
END SUBROUTINE INT_SORTING_METHOD


SUBROUTINE STR_SORTING_METHOD(CDNEW, IOLD)
CHARACTER(LEN=*), INTENT(IN) :: CDNEW
INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: IOLD
CHARACTER(LEN=LEN(CDNEW)) CLNEW
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ECSORT_MIX:STR_SORTING_METHOD',0,ZHOOK_HANDLE)
CLNEW = CDNEW
CALL TOUPPER(CLNEW)
SELECT CASE (CLNEW)
CASE ('RADIX', 'RADIXSORT')
  CALL SORTING_METHOD(RADIXSORT_METHOD, IOLD)
CASE ('HEAP', 'HEAPSORT')
  CALL SORTING_METHOD(HEAPSORT_METHOD, IOLD)
CASE ('QUICK', 'QUICKSORT', 'QSORT')
  CALL SORTING_METHOD(QUICKSORT_METHOD, IOLD)
CASE ('COUNT', 'COUNTINGSORT', 'COUNTSORT')
  CALL SORTING_METHOD(COUNTINGSORT_METHOD, IOLD)
CASE ('GNOME', 'GNOMESORT')
  CALL SORTING_METHOD(GNOMESORT_METHOD, IOLD)
CASE ('DEFAULT', 'DEF')
  CALL SORTING_METHOD(DEFAULT_METHOD, IOLD)
CASE DEFAULT
  CALL SORTING_METHOD(DEFAULT_METHOD, IOLD)
END SELECT
IF (LHOOK) CALL DR_HOOK('ECSORT_MIX:STR_SORTING_METHOD',1,ZHOOK_HANDLE)
END SUBROUTINE STR_SORTING_METHOD


SUBROUTINE INIT_INDEX(INDEX, INDEX_ADJ)
INTEGER(KIND=JPIM), INTENT(OUT):: INDEX(:)
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: INDEX_ADJ
INTEGER(KIND=JPIM) :: I, N
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ECSORT_MIX:INIT_INDEX',0,ZHOOK_HANDLE)
N = SIZE(INDEX)
IF (PRESENT(INDEX_ADJ)) THEN
  DO I=1,N
    INDEX(I) = I + INDEX_ADJ
  ENDDO
ELSE
  DO I=1,N
    INDEX(I) = I
  ENDDO
ENDIF
IF (LHOOK) CALL DR_HOOK('ECSORT_MIX:INIT_INDEX',1,ZHOOK_HANDLE)
END SUBROUTINE INIT_INDEX


SUBROUTINE ADJUST_INDEX(INDEX, INDEX_ADJ)
INTEGER(KIND=JPIM), INTENT(INOUT):: INDEX(:)
INTEGER(KIND=JPIM), INTENT(IN) :: INDEX_ADJ
INTEGER(KIND=JPIM) :: I, N
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ECSORT_MIX:ADJUST_INDEX',0,ZHOOK_HANDLE)
IF (INDEX_ADJ /= 0) THEN
  N = SIZE(INDEX)
  DO I=1,N
    INDEX(I) = INDEX(I) + INDEX_ADJ
  ENDDO
ENDIF
IF (LHOOK) CALL DR_HOOK('ECSORT_MIX:ADJUST_INDEX',1,ZHOOK_HANDLE)
END SUBROUTINE ADJUST_INDEX


SUBROUTINE GET_RANK(INDEX, RANK, INDEX_ADJ)
INTEGER(KIND=JPIM), INTENT(IN) :: INDEX(:)
INTEGER(KIND=JPIM), INTENT(OUT):: RANK(:)
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: INDEX_ADJ
INTEGER(KIND=JPIM) :: I, N
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ECSORT_MIX:GET_RANK',0,ZHOOK_HANDLE)
N = MIN(SIZE(INDEX),SIZE(RANK))
IF (PRESENT(INDEX_ADJ)) THEN
  DO I=1,N
    RANK(INDEX(I)+INDEX_ADJ) = I
  ENDDO
ELSE
  DO I=1,N
    RANK(INDEX(I)) = I
  ENDDO
ENDIF
IF (LHOOK) CALL DR_HOOK('ECSORT_MIX:GET_RANK',1,ZHOOK_HANDLE)
END SUBROUTINE GET_RANK


#undef INT_VERSION
#undef REAL_VERSION

!-- Create version for INTEGER(KIND=JPIM)
#define INT_VERSION   4
#define REAL_VERSION  0
#include "ecsort_shared.h"
#undef INT_VERSION
#undef REAL_VERSION

!-- Create version for INTEGER(KIND=JPIB)
#define INT_VERSION   8
#define REAL_VERSION  0
#include "ecsort_shared.h"
#undef INT_VERSION
#undef REAL_VERSION

!-- Create version for REAL(KIND=JPRM)
#define INT_VERSION   0
#define REAL_VERSION  4
#include "ecsort_shared.h"
#undef INT_VERSION
#undef REAL_VERSION

!-- Create version for REAL(KIND=JPRD)
#define INT_VERSION   0
#define REAL_VERSION  8
#include "ecsort_shared.h"
#undef INT_VERSION
#undef REAL_VERSION

END MODULE ECSORT_MIX
