! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE OML_MOD

!-- the following system specific omp_lib-module is not always available (e.g. pgf90)
!! use omp_lib

USE EC_PARKIND  ,ONLY : JPIM, JPIB
USE EC_LUN, ONLY : NULERR

!**SS/18-Feb-2005
!--Dr.Hook references removed, because these locks may also be
!  called from within drhook.c itself !! 
!--Also, there could be considerable & unjustified overhead
!  when using Dr.Hook in such a low level

!**SS/15-Dec-2005
!--The size of lock-variables are now OMP_LOCK_KIND as of in OMP_LIB,
!  and OMP_LOCK_KIND is aliased to OML_LOCK_KIND
!  OMP_LOCK_KIND is usually 4 in 32-bit addressing mode
!                           8 in 64-bit addressing mode
!--M_OML_LOCK changed to M_EVENT and kept as 32-bit int
!--OML_FUNCT changed to OML_TEST_EVENT
!--M_LOCK initialized to -1
!--M_EVENT initialized to 0
!--Added intent(s)
!--Support for omp_lib (but not always available)
!--Locks can now also be set/unset OUTSIDE the parallel regions
!--Added routine OML_TEST_LOCK (attempts to set lock, but if *un*successful, does NOT  block)
!--Buffer-zone for M_LOCK; now a vector of 2 elements in case problems/inconsistencies with OMP_LOCK_KIND 4/8

!**SS/22-Feb-2006
!--Locking routines are doing nothing unless OMP_GET_MAX_THREADS() > 1
!  This is to avoid unacceptable deadlocks/timeouts with signal handlers when
!  the only thread receives signal while inside locked region
!--Affected routines: OML_TEST_LOCK()  --> always receives .TRUE.
!                     OML_SET_LOCK()   --> sets nothing
!                     OML_UNSET_LOCK() --> unsets nothing
!                     OML_INIT_LOCK()  --> inits nothing

!**SS/11-Sep-2006
!--Added OML_DEBUG feature

!**REK/18-Jul-2007
!--Protected OML_DESTROY_LOCK

!**REK/07-Sep-2007
!--Add OMP FLUSH feature

!**SS/05-Dec-2007
!--Added routine OML_NUM_THREADS([optional_new_number_of_threads])
! 1) To adjust [reduce] the number of threads working in concert
!    Accepts only # of threads between 1 and the max # of threads (i.e. from export OMP_NUM_THREADS=<value>)
! 2) Returns the previous active number of threads
! 3) Can be called from outside the OpenMP-parallel region only

!**SS/14-Dec-2007
!--The routine OML_NUM_THREADS() now also accepts character string (= environment variable)
!  as the sole argumentoz
!--You could now set effective number of threads (<= $OMP_NUM_THREADS) to the value of
!  particular environment variable; f.ex.:
!  export OML_MSGPASS_OBSDATA_READ=8 and call to OML_NUM_THREADS('OML_MSGPASS_OBSDATA_READ')
!  would set the effective no. of threads to (max) 8 when reading obs. wiz msgpass_obsdata

!**SS/09-May-2008
!-- OML_NUM_THREADS() did not work as expected since I misinterpreted the meaning of 
!   the OpenMP-function OMP_GET_NUM_THREADS()
!-- With two PRIVATE [to this module] variables the bug will get sorted out
!   + a new routine OML_INIT() was added (to be called from MPL_INIT or so)

!**FV/27-May-2009
!-- OML_GET_NUM_THREADS()

!**FV/09-Oct-2018
!-- OML_INCR_COUNTER to increment and check before launching synchronization

!**Willem Deconinck/02-Feb-2022
!-- Deprecate OML_GET_NUM_THREADS() and remove overloads as they were not used
!   Use OML_GET_MAX_THREADS instead.
!-- Initialize OML_ABSMAX_THREADS within OML_INIT and abort if OML_INIT was forgotten
!   before accessing OML_ABSMAX_THREADS

!**Sami Saarinen/19-Feb-2022
!-- Final clean-ups & fixes for this module
!:: Rename (and bring back) OML_NUM_THREADS(<value>) as OML_SET_NUM_THREADS(<value>)
!:: Bring back OML_GET_NUM_THREADS() == the *actual* number of threads WHEN IN a parallel region (otherwise == 1)
!:: The OML_GET_MAX_THREADS() is the *currently* available max threads, which could go up|down
!   if OML|OMP_SET_NUM_THREADS(<value>) was called explicitly (NOT recommended -- messes up thread affinity)
!:: Got rid of the confusing OML_ABSMAX_THREADS() -- instead relying solely on OML_MAX_THREADS()
!   Now N_OML_MAX_THREADS holds the max allowed value (upon OML_INIT())
!   this variable never changes and *should* be the high water mark for number of threads ever per run
!   equivalent to the value of "export OMP_NUM_THREADS=<value>" or the system shell default; e.g. that of "nproc --all"
!   Caveat & warning: User may still increase beyond the N_OML_MAX_THREADS by explicitly calling OMP_SET_NUM_THREADS(<value>) !!!
!:: OML_GET_MAX_THREADS() is now an alias to OML_MAX_THREADS()
!:: OML_INIT()'s initialization of N_OML_MAX_THREADS only done if OUTSIDE the parallel region -- otherwise task aborted

IMPLICIT NONE

SAVE

PRIVATE

LOGICAL :: OML_DEBUG = .FALSE.

INTERFACE OML_SET_NUM_THREADS
   MODULE PROCEDURE &
        & OML_SET_NUM_THREADS_INT, &
        & OML_SET_NUM_THREADS_STR
END INTERFACE OML_SET_NUM_THREADS

INTERFACE OML_GET_MAX_THREADS
   MODULE PROCEDURE &
        & OML_MAX_THREADS
END INTERFACE OML_GET_MAX_THREADS

PUBLIC OML_WAIT_EVENT, OML_SET_EVENT, OML_INCR_EVENT, &
   &   OML_MY_THREAD,  OML_MAX_THREADS , OML_GET_MAX_THREADS, OML_OMP, &
   &   OML_IN_PARALLEL, OML_TEST_EVENT, OML_INCR_COUNTER, &
   &   OML_UNSET_LOCK, OML_INIT_LOCK, OML_SET_LOCK, OML_DESTROY_LOCK, &
   &   OML_LOCK_KIND, OML_TEST_LOCK, OML_DEBUG, OML_SET_NUM_THREADS, &
   &   OML_INIT, OML_GET_NUM_THREADS

!-- The following should normally be 4 in 32-bit addressing mode
!                                    8 in 64-bit addressing mode
! Since system specific omp_lib-module is not always available (e.g. pgf90)
! we hardcode OML_LOCK_KIND to JPIB (usually 8) for now
!!INTEGER(KIND=JPIM), PARAMETER :: OML_LOCK_KIND = OMP_LOCK_KIND
INTEGER(KIND=JPIM), PARAMETER :: OML_LOCK_KIND = JPIB

!-- Note: Still JPIM !!
INTEGER(KIND=JPIM) :: M_EVENT = 0

!-- Note: OML_LOCK_KIND, not JPIM !!
INTEGER(KIND=OML_LOCK_KIND) :: M_LOCK(2) = (/-1, -1/)

!-- The two PRIVATE [to this module] variables
INTEGER(KIND=JPIM) :: N_OML_MAX_THREADS = -1

!-- OMP function declarations (that require knowledge of their type) now in one place:

#include "abor1.intfb.h"

CONTAINS

SUBROUTINE OML_INIT()
CHARACTER(LEN=*), PARAMETER :: CLMSG = 'Fatal coding error : Cannot CALL OML_INIT() INSIDE the OpenMP parallel region'
!$ INTEGER(KIND=JPIM) :: OMP_GET_MAX_THREADS
!$ LOGICAL :: OMP_IN_PARALLEL
IF (N_OML_MAX_THREADS == -1) THEN
   IF (OML_IN_PARALLEL()) THEN
      WRITE(NULERR,'(1X,A)') CLMSG
      CALL ABOR1(CLMSG)
      RETURN
   ENDIF
   N_OML_MAX_THREADS = 1
!$ N_OML_MAX_THREADS = OMP_GET_MAX_THREADS()
ENDIF
END SUBROUTINE OML_INIT

FUNCTION OML_OMP()
LOGICAL :: OML_OMP
   OML_OMP=.FALSE.
!$ OML_OMP=.TRUE.
END FUNCTION OML_OMP

FUNCTION OML_IN_PARALLEL()
LOGICAL :: OML_IN_PARALLEL
!$ INTEGER(KIND=JPIM) :: OMP_GET_NUM_THREADS
!$ LOGICAL :: OMP_IN_PARALLEL
   OML_IN_PARALLEL=.FALSE.
!$ OML_IN_PARALLEL=((OMP_GET_NUM_THREADS() > 1).AND.OMP_IN_PARALLEL())
END FUNCTION OML_IN_PARALLEL

FUNCTION OML_TEST_LOCK(MYLOCK)
LOGICAL :: OML_TEST_LOCK
INTEGER(KIND=OML_LOCK_KIND),INTENT(INOUT),OPTIONAL :: MYLOCK
!$ INTEGER(KIND=JPIM) :: OMP_GET_MAX_THREADS
!$ LOGICAL :: OMP_TEST_LOCK
OML_TEST_LOCK = .TRUE.
!$ IF (OMP_GET_MAX_THREADS() > 1) THEN
!$   IF (PRESENT(MYLOCK)) THEN
!$     OML_TEST_LOCK = OMP_TEST_LOCK(MYLOCK)
!$   ELSE
!$     OML_TEST_LOCK = OMP_TEST_LOCK(M_LOCK(1))
!$   ENDIF
!$ ENDIF
END FUNCTION OML_TEST_LOCK

SUBROUTINE OML_UNSET_LOCK(MYLOCK)
INTEGER(KIND=OML_LOCK_KIND),INTENT(INOUT),OPTIONAL :: MYLOCK
!$ INTEGER(KIND=JPIM) :: OMP_GET_MAX_THREADS
!$ IF (OMP_GET_MAX_THREADS() > 1) THEN
!$   IF (PRESENT(MYLOCK)) THEN
!$     CALL OMP_UNSET_LOCK(MYLOCK)
!$   ELSE
!$     CALL OMP_UNSET_LOCK(M_LOCK(1))
!$   ENDIF
!$ ENDIF
END SUBROUTINE OML_UNSET_LOCK

SUBROUTINE OML_SET_LOCK(MYLOCK)
INTEGER(KIND=OML_LOCK_KIND),INTENT(INOUT),OPTIONAL :: MYLOCK
!$ INTEGER(KIND=JPIM) :: OMP_GET_MAX_THREADS
!$ IF (OMP_GET_MAX_THREADS() > 1) THEN
!$   IF (PRESENT(MYLOCK)) THEN
!$     CALL OMP_SET_LOCK(MYLOCK)
!$   ELSE
!$     CALL OMP_SET_LOCK(M_LOCK(1))
!$   ENDIF
!$ ENDIF
END SUBROUTINE OML_SET_LOCK

SUBROUTINE OML_INIT_LOCK(MYLOCK)
INTEGER(KIND=OML_LOCK_KIND),INTENT(INOUT),OPTIONAL :: MYLOCK
!$ INTEGER(KIND=JPIM) :: OMP_GET_MAX_THREADS
!$ IF (OMP_GET_MAX_THREADS() > 1) THEN
!$   IF (PRESENT(MYLOCK)) THEN
!$     CALL OMP_INIT_LOCK(MYLOCK)
!$   ELSE
!$     CALL OMP_INIT_LOCK(M_LOCK(1))
!$   ENDIF
!$ ENDIF
END SUBROUTINE OML_INIT_LOCK

SUBROUTINE OML_DESTROY_LOCK(MYLOCK)
INTEGER(KIND=OML_LOCK_KIND),INTENT(INOUT),OPTIONAL :: MYLOCK
!$ INTEGER(KIND=JPIM) :: OMP_GET_MAX_THREADS
!$ IF (OMP_GET_MAX_THREADS() > 1) THEN
!$   IF (PRESENT(MYLOCK)) THEN
!$     CALL OMP_DESTROY_LOCK(MYLOCK)
!$   ELSE
!$     CALL OMP_DESTROY_LOCK(M_LOCK(1))
!$   ENDIF
!$ ENDIF
END SUBROUTINE OML_DESTROY_LOCK

FUNCTION OML_TEST_EVENT(K,MYEVENT)
LOGICAL :: OML_TEST_EVENT
INTEGER(KIND=JPIM),INTENT(IN) :: K,MYEVENT
!$OMP FLUSH
IF (K.EQ.MYEVENT) THEN
 OML_TEST_EVENT =.TRUE.
ELSE
 OML_TEST_EVENT=.FALSE.
ENDIF
END FUNCTION OML_TEST_EVENT

SUBROUTINE OML_WAIT_EVENT(K,MYEVENT)
INTEGER(KIND=JPIM),INTENT(IN) :: K
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: MYEVENT
IF (PRESENT(MYEVENT)) THEN
  DO
    IF (OML_TEST_EVENT(K,MYEVENT)) EXIT
  ENDDO
ELSE
  DO
    IF (OML_TEST_EVENT(K,M_EVENT)) EXIT
  ENDDO
ENDIF
END SUBROUTINE OML_WAIT_EVENT

SUBROUTINE OML_SET_EVENT(K,MYEVENT)
INTEGER(KIND=JPIM),INTENT(IN) :: K
INTEGER(KIND=JPIM),INTENT(OUT),OPTIONAL :: MYEVENT
IF (PRESENT(MYEVENT)) THEN
  MYEVENT=K
ELSE
  M_EVENT=K
ENDIF
END SUBROUTINE OML_SET_EVENT

SUBROUTINE OML_INCR_EVENT(K,MYEVENT)
INTEGER(KIND=JPIM) :: K
INTEGER(KIND=JPIM),INTENT(INOUT),OPTIONAL :: MYEVENT
!$OMP FLUSH
IF (PRESENT(MYEVENT)) THEN
!$OMP ATOMIC
  MYEVENT=MYEVENT+K
ELSE
!$OMP ATOMIC
  M_EVENT=M_EVENT+K
ENDIF
!$OMP FLUSH
END SUBROUTINE OML_INCR_EVENT

SUBROUTINE OML_INCR_COUNTER(K,KMAX)
INTEGER(KIND=JPIM),INTENT(INOUT) :: K
INTEGER(KIND=JPIM),INTENT(IN) :: KMAX
! Increment by 1
K=K+1
! Security check
IF (K > KMAX)  THEN
   WRITE(NULERR,'("OML_INCR_COUNTER: ILOCK > SIZE(KLOCK), ILOCK=",I0," SIZE(KLOCK)=",I0)') K,KMAX
   CALL ABOR1('PLEASE INCREASE JP_LOCKS IN CALL_SL_AD')
   RETURN
ENDIF
END SUBROUTINE OML_INCR_COUNTER


FUNCTION OML_MY_THREAD()
INTEGER(KIND=JPIM) :: OML_MY_THREAD
!$ INTEGER(KIND=JPIM) :: OMP_GET_THREAD_NUM
   OML_MY_THREAD = 1
!$ OML_MY_THREAD = OMP_GET_THREAD_NUM() + 1
END FUNCTION OML_MY_THREAD

FUNCTION OML_MAX_THREADS() ! alias OML_GET_MAX_THREADS()
INTEGER(KIND=JPIM) :: OML_MAX_THREADS
!$ INTEGER(KIND=JPIM) :: OMP_GET_MAX_THREADS
IF (N_OML_MAX_THREADS == -1) CALL OML_INIT() ! Harmless as usually called OUTSIDE the parallel region (checked)
   OML_MAX_THREADS = 1
!$ OML_MAX_THREADS = OMP_GET_MAX_THREADS()
END FUNCTION OML_MAX_THREADS

FUNCTION OML_GET_NUM_THREADS()
INTEGER(KIND=JPIM) :: OML_GET_NUM_THREADS
!$ INTEGER(KIND=JPIM) :: OMP_GET_NUM_THREADS
   OML_GET_NUM_THREADS = 1
!$ OML_GET_NUM_THREADS = OMP_GET_NUM_THREADS()
END FUNCTION OML_GET_NUM_THREADS

FUNCTION OML_SET_NUM_THREADS_INT(KOMP_SET_THREADS)
INTEGER(KIND=JPIM) :: OML_SET_NUM_THREADS_INT
INTEGER(KIND=JPIM),INTENT(IN) :: KOMP_SET_THREADS
!$ INTEGER(KIND=JPIM) :: OMP_GET_MAX_THREADS
!$ LOGICAL :: OMP_IN_PARALLEL
IF (N_OML_MAX_THREADS == -1) CALL OML_INIT() ! Harmless as usually called OUTSIDE the parallel region (checked)
   OML_SET_NUM_THREADS_INT = 1
!$ OML_SET_NUM_THREADS_INT = OMP_GET_MAX_THREADS()
!$ IF (KOMP_SET_THREADS /= OML_SET_NUM_THREADS_INT) THEN
!$   IF (KOMP_SET_THREADS >= 1 .AND. KOMP_SET_THREADS <= N_OML_MAX_THREADS) THEN
!- This is the absolute max no. of threads allowed   --> ^^^^^^^^^^^^^^^^^ <--
!$     IF (.NOT.OMP_IN_PARALLEL()) THEN ! Change *only* if called from OUTSIDE the OpenMP-parallel region
!$       CALL OMP_SET_NUM_THREADS(KOMP_SET_THREADS) ! Warning: could mess up your thread affinity !!
!$     ENDIF
!$   ENDIF
!$ ENDIF
END FUNCTION OML_SET_NUM_THREADS_INT

FUNCTION OML_SET_NUM_THREADS_STR(CD_ENV)
INTEGER(KIND=JPIM) :: OML_SET_NUM_THREADS_STR
CHARACTER(LEN=*),INTENT(IN) :: CD_ENV
!$ character(len=20) CLvalue
!$ INTEGER(KIND=JPIM) :: ivalue
!$ INTEGER(KIND=JPIM) :: OMP_GET_MAX_THREADS
   OML_SET_NUM_THREADS_STR = 1
!$ OML_SET_NUM_THREADS_STR = OMP_GET_MAX_THREADS()
!$ IF (LEN(CD_ENV) > 0) THEN
!$   CALL GET_ENVIRONMENT_VARIABLE(CD_ENV,CLvalue)
!$   IF (CLvalue /= ' ') THEN
!$     READ(CLvalue,'(i20)',end=99,err=99) ivalue
!$     OML_SET_NUM_THREADS_STR = OML_SET_NUM_THREADS_INT(ivalue)
!$   ENDIF
!$ 99 continue
!$ ENDIF
END FUNCTION OML_SET_NUM_THREADS_STR

!================================================================================================================================
! C bindings ( Signatures in oml.h must match )
!================================================================================================================================

SUBROUTINE OML_SET_DEBUG_BINDC(KONOFF) BIND(C,name="oml_set_debug")
USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
IMPLICIT NONE
INTEGER(KIND=C_INT), VALUE, INTENT(IN) :: KONOFF
IF (KONOFF == 0) THEN 
      OML_DEBUG = .FALSE.
ELSE
      OML_DEBUG = .TRUE.
ENDIF
END SUBROUTINE

FUNCTION OML_GET_DEBUG_BINDC() BIND(C,name="oml_get_debug") RESULT(KRET)
USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
IMPLICIT NONE
INTEGER(KIND=C_INT) :: KRET
KRET = 0
IF (OML_DEBUG) KRET = 1
END FUNCTION

SUBROUTINE OML_INIT_LOCKID_WITH_NAME_BINDC(KMYLOCK,CDLOCKNAME) BIND(C,NAME="oml_init_lockid_with_name")
USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INTPTR_T, C_CHAR
IMPLICIT NONE
INTEGER(KIND=OML_LOCK_KIND), INTENT(INOUT) :: KMYLOCK
CHARACTER(KIND=C_CHAR,LEN=1), INTENT(IN)   :: CDLOCKNAME(*)
INTEGER(KIND=C_INTPTR_T), EXTERNAL :: LOC_ADDR
CALL OML_INIT_LOCK(KMYLOCK)
IF (OML_DEBUG) WRITE(0,'(1x,a,2i20)') &
      & 'oml_init_lockid_with_name "'//from_c_str(CDLOCKNAME)//'" :',KMYLOCK,LOC_ADDR(KMYLOCK)

CONTAINS

      function from_c_str(s) result(string)
      use, intrinsic :: iso_c_binding, only : c_char, c_null_char
      character(kind=c_char,len=1), intent(in) :: s(*)
      character(len=:), allocatable :: string
      integer i, nchars
      i = 1
      do
          if (s(i) == c_null_char) exit
          i = i + 1
      enddo
      nchars = i - 1  ! Exclude null character from Fortran string
      allocate( character(len=(nchars)) :: string )
      do i=1,nchars
          string(i:i) = s(i)
      enddo
      end function
END SUBROUTINE

SUBROUTINE OML_INIT_LOCKID_BINDC(KMYLOCK) BIND(C,NAME="oml_init_lockid")
USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INTPTR_T
IMPLICIT NONE
INTEGER(KIND=OML_LOCK_KIND), INTENT(INOUT) :: KMYLOCK
INTEGER(KIND=C_INTPTR_T), EXTERNAL :: LOC_ADDR
CALL OML_INIT_LOCK(KMYLOCK)
IF (OML_DEBUG) WRITE(0,'(1x,2i20)') &
      & 'oml_init_lockid :',KMYLOCK,LOC_ADDR(KMYLOCK)
END SUBROUTINE

SUBROUTINE OML_INIT_LOCK_BINDC() BIND(C,NAME="oml_init_lock")
IMPLICIT NONE
CALL OML_INIT_LOCK()
END SUBROUTINE

FUNCTION OML_TEST_LOCKID_BINDC(KMYLOCK) BIND(C,NAME="oml_test_lockid") RESULT(KISSET)
USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
IMPLICIT NONE
INTEGER(KIND=C_INT) :: KISSET
INTEGER(KIND=OML_LOCK_KIND), INTENT(INOUT) :: KMYLOCK
KISSET = 1
IF (.NOT.OML_TEST_LOCK(KMYLOCK)) KISSET = 0
END FUNCTION

FUNCTION OML_TEST_LOCK_BINDC() BIND(C,NAME="oml_test_lock") RESULT(KISSET)
USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
IMPLICIT NONE
INTEGER(KIND=C_INT) :: KISSET
KISSET = 1
IF (.NOT.OML_TEST_LOCK()) KISSET = 0
END FUNCTION

SUBROUTINE OML_SET_LOCKID_BINDC(KMYLOCK) BIND(C,NAME="oml_set_lockid")
USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INTPTR_T
USE EC_PARKIND, ONLY : JPRD
IMPLICIT NONE
INTEGER(KIND=OML_LOCK_KIND), INTENT(INOUT) :: KMYLOCK
INTEGER(KIND=C_INTPTR_T), EXTERNAL :: LOC_ADDR
REAL(KIND=JPRD), EXTERNAL :: UTIL_WALLTIME
IF (OML_DEBUG) WRITE(0,'(1x,f20.6,1x,i3,a,2i20)') &
      & UTIL_WALLTIME(),OML_MY_THREAD(),': oml_set_lockid >>',KMYLOCK,LOC_ADDR(KMYLOCK)
CALL OML_SET_LOCK(KMYLOCK)
IF (OML_DEBUG) WRITE(0,'(1x,f20.6,1x,i3,a,2i20)') &
      & UTIL_WALLTIME(),OML_MY_THREAD(),': oml_set_lockid <<',KMYLOCK,LOC_ADDR(KMYLOCK)
END SUBROUTINE

SUBROUTINE OML_SET_LOCK_BINDC() BIND(C,NAME="oml_set_lock")
IMPLICIT NONE
CALL OML_SET_LOCK()
END SUBROUTINE

SUBROUTINE OML_UNSET_LOCK_BINDC() BIND(C,NAME="oml_unset_lock")
IMPLICIT NONE
CALL OML_UNSET_LOCK()
END SUBROUTINE

SUBROUTINE OML_UNSET_LOCKID_BINDC(KMYLOCK) BIND(C,NAME="oml_unset_lockid")
USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INTPTR_T
USE EC_PARKIND, ONLY : JPRD
IMPLICIT NONE
INTEGER(KIND=OML_LOCK_KIND), INTENT(INOUT) :: KMYLOCK
INTEGER(KIND=C_INTPTR_T), EXTERNAL :: LOC_ADDR
REAL(KIND=JPRD), EXTERNAL :: UTIL_WALLTIME
IF (OML_DEBUG) WRITE(0,'(1x,f20.6,1x,i3,a,2i20)') &
      & UTIL_WALLTIME(),OML_MY_THREAD(),': oml_unset_lockid >>',KMYLOCK,LOC_ADDR(KMYLOCK)
CALL OML_UNSET_LOCK(KMYLOCK)
IF (OML_DEBUG) WRITE(0,'(1x,f20.6,1x,i3,a,2i20)') &
      & UTIL_WALLTIME(),OML_MY_THREAD(),': oml_unset_lockid <<',KMYLOCK,LOC_ADDR(KMYLOCK)
END SUBROUTINE

FUNCTION OML_IN_PARALLEL_BINDC() BIND(C,NAME="oml_in_parallel") RESULT(KISPAR_REGION)
USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
IMPLICIT NONE
INTEGER(KIND=C_INT) :: KISPAR_REGION
KISPAR_REGION = 0
IF (OML_IN_PARALLEL()) KISPAR_REGION = 1
END FUNCTION

FUNCTION OML_GET_MAX_THREADS_BINDC() BIND(C,NAME="oml_get_max_threads") RESULT(KTIDS)
USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
IMPLICIT NONE
INTEGER(KIND=C_INT) :: KTIDS
KTIDS = OML_GET_MAX_THREADS()
END FUNCTION

FUNCTION OML_GET_NUM_THREADS_BINDC() BIND(C,NAME="oml_get_num_threads") RESULT(KTIDS)
USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
IMPLICIT NONE
INTEGER(KIND=C_INT) :: KTIDS
KTIDS = OML_GET_NUM_THREADS()
END FUNCTION

FUNCTION OML_MY_THREAD_BINDC() BIND(C,NAME="oml_my_thread") RESULT(KMYTID)
USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
IMPLICIT NONE
INTEGER(KIND=C_INT) :: KMYTID
KMYTID = OML_MY_THREAD()
END FUNCTION

END MODULE OML_MOD
