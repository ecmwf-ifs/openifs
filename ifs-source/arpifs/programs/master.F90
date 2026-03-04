! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

PROGRAM MASTER

USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : STDOUT=>OUTPUT_UNIT
#ifdef WITH_FCKIT
USE FCKIT_MODULE, ONLY : FCKIT_MAIN, FCKIT_LOG, FCKIT_EXCEPTION, FCKIT_EXCEPTION_HANDLER
#endif
USE PARKIND1   ,ONLY : JPIM, JPIB, JPRB, JPRD
USE YOMHOOK    ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE MPL_MODULE ,ONLY : MPL_INIT, MPL_END
USE OML_MOD    ,ONLY : OML_MAX_THREADS
USE OMP_LIB    ,ONLY : OMP_SET_NUM_THREADS
USE STACK_MIX  ,ONLY : INIT_STACK
USE YOMIO_SERV ,ONLY : IO_SERV_C001
USE YOMFP_SERV ,ONLY : FP_SERV_C001
USE MPL_MPIF, ONLY : MPI_COMM_WORLD, MPI_THREAD_MULTIPLE, MPI_THREAD_SINGLE, MPI_WTIME
USE MPL_MODULE, ONLY : LMPLUSERCOMM, MPLUSERCOMM, LTHSAFEMPI, LINITMPI_VIA_MPL
#if defined(WITH_OASIS) || defined(WITH_NEMO)
USE COUPLING
#endif

IMPLICIT NONE

! P. Marguinaud : 01-Jan-2001 : Support for IO server
! R. El Khatib  : 22-Mar-2011 : conditional support for IBM high perf monitoring
! P. Marguinaud : 10-Oct-2013 : Change IO server init & exit routines

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
LOGICAL :: LLHOOK_SAVE
INTEGER (KIND=JPIM) :: ICOMM_IFS_PLUS_IFS_IO_SERV
INTEGER (KIND=JPIM) :: ICOMM_IFS_PLUS_NEMO_IO_SERV

LOGICAL :: LLINIT, LLSERV, LLCOUPACTIVE, LLMPI1, LLSTARTUPCOST
INTEGER :: IERR, ICOL, INUM, IREQUIRED,IPROVIDED,IME

LOGICAL :: LLNEMOIO, LLNEMOIOSERVER
CHARACTER(LEN=512) :: CLIOSERV_LOGFILE
#ifdef WITH_FCKIT
PROCEDURE(FCKIT_EXCEPTION_HANDLER), POINTER :: FUNPTR
#endif
INTEGER(KIND=JPIM) :: IER
CHARACTER(LEN=20) :: CL_MPI_EPOCH ! 10-digits in seconds + dot + 9-digits nanosecs
REAL(KIND=JPRD) :: ZMPI_INIT(2)

#include "abor1.intfb.h"
#include "cnt0.intfb.h"
!#include "couplo4_inimpi.intfb.h"
!#include "couplo4_endmpi.intfb.h"
#include "io_serv_init_part1.intfb.h"
#include "io_serv_init_part2.intfb.h"
#include "io_serv_exit.intfb.h"
#include "io_serv_run.intfb.h"
#include "fp_serv_init_part1.intfb.h"
#include "fp_serv_init_part2.intfb.h"
#include "fp_serv_exit.intfb.h"
#include "ininemoio.intfb.h"
#include "ininemoio2.intfb.h"
#include "ec_meminfo.intfb.h"

LLHOOK_SAVE = LHOOK
LHOOK = .FALSE.
LLNEMOIOSERVER = .FALSE.
LLNEMOIO = .FALSE.

IME = -1
LLSTARTUPCOST = .FALSE. ! If true, then display MPI startup cost (only ever to happen on the global master task IME == 0)
ZMPI_INIT(:) = 0

! XIOS and MPI initialization
! #ifdef 1
! #ifdef WITH_XIOS
! CALL SUXIOS_INI
! #elif WITH_CPLNG
! CALL CPLNG_INIT
! #endif
! #endif

! OASIS3 or OASIS4 interface must be initialized before any DR_HOOK call.

#if defined(WITH_OASIS)
CALL CPL_INIT(LLCOUPACTIVE)
#endif

#ifdef WITH_NEMO
#ifndef MPI1
LLMPI1 = .FALSE.
IREQUIRED = MPI_THREAD_MULTIPLE
IPROVIDED = MPI_THREAD_SINGLE
CALL ININEMOIO (LLNEMOIO,LLNEMOIOSERVER,IREQUIRED,IPROVIDED,LLMPI1)
LTHSAFEMPI = (IPROVIDED >= IREQUIRED)
#else
LLMPI1 = .TRUE.
IREQUIRED = -1
IPROVIDED = -1
CALL ININEMOIO (LLNEMOIO,LLNEMOIOSERVER,IREQUIRED,IPROVIDED,LLMPI1)
LTHSAFEMPI = .FALSE.
#endif
#endif

CALL MPI_INITIALIZED (LLINIT, IERR)
IF (IERR /= 0) CALL ABOR1 ('MASTER: MPI_INITIALIZED FAILED')

IF (.NOT. LLINIT) THEN
#ifndef MPI1
  LLMPI1 = .FALSE.
  IREQUIRED = MPI_THREAD_MULTIPLE
  IPROVIDED = MPI_THREAD_SINGLE
  CALL MPI_INIT_THREAD(IREQUIRED,IPROVIDED,IERR)
  IF (IERR /= 0) CALL ABOR1 ('MASTER: MPI_INIT_THREAD FAILED')
  LTHSAFEMPI = (IPROVIDED >= IREQUIRED)
#else
  LLMPI1 = .TRUE.
  IREQUIRED = -1
  IPROVIDED = -1
  CALL MPI_INIT(IERR)
  IF (IERR /= 0) CALL ABOR1 ('MASTER: MPI_INIT FAILED')
  LTHSAFEMPI = .FALSE.
#endif
  LINITMPI_VIA_MPL = .TRUE. ! To re-instate ec_meminfo-call from within ec_mpi_finalize @ mpl_end()
ENDIF

IF (LINITMPI_VIA_MPL.OR.LLNEMOIO) THEN

  CALL MPI_Comm_rank(MPI_COMM_WORLD, IME, IERR)

  ! MPI_Init* overhead -- reference time from MPI_EPOCH variable (obtained via command date +%s.%N)
  CALL GET_ENVIRONMENT_VARIABLE('MPI_EPOCH',CL_MPI_EPOCH)
  IF (CL_MPI_EPOCH /= ' ') THEN
     CALL MPI_Barrier(MPI_COMM_WORLD, IERR)
     IF (IME == 0) THEN
        ZMPI_INIT(2) = MPI_Wtime() ! *not* from ifsaux/support/env.c
        !READ(CL_MPI_EPOCH,'(f20.0)',err=999,end=999) ZMPI_INIT(1)
        READ(CL_MPI_EPOCH,'(f20.0)') ZMPI_INIT(1)
        LLSTARTUPCOST = .TRUE.
!999     CONTINUE
     ENDIF
  ENDIF

  ! Print out thread safety etc. messages -- must use MPI_Comm_rank since MPL not initialized just yet
  IF (IME == 0) THEN
     WRITE(0,'(1X,A,4(1X,I0),2(1X,L1))') &
          & 'MAIN: IREQUIRED, MPI_THREAD_MULTIPLE, MPI_THREAD_SINGLE, IPROVIDED, LTHSAFEMPI, LLMPI1 =',&
          &        IREQUIRED, MPI_THREAD_MULTIPLE, MPI_THREAD_SINGLE, IPROVIDED, LTHSAFEMPI, LLMPI1
  ENDIF
ENDIF

! Every program needs to be initialised
#ifdef WITH_FCKIT
IF(.NOT.LLNEMOIOSERVER) THEN
  CALL FCKIT_MAIN%INITIALISE()

#ifndef ONT_REGISTER_FCKIT_HANDLER
  ! Register ABOR1 as fckit's exception handler
  FUNPTR => ABOR1_EXCEPTION_HANDLER
  CALL FCKIT_EXCEPTION%SET_HANDLER( FUNPTR )
#endif
ENDIF
#endif

IF (LMPLUSERCOMM) THEN
  ICOMM_IFS_PLUS_IFS_IO_SERV = MPLUSERCOMM
ELSE
  ICOMM_IFS_PLUS_IFS_IO_SERV = MPI_COMM_WORLD
ENDIF





IO_SERV_C001%LIO_SERVER=.FALSE.


! Setup default for eckit/fckit enabled logging
!  - Worker with proc==1 writes to STDOUT, others disabled
!  - IO-servers write each to a server specific file IOSERV.<IOPROC>
! In further routines this may be ammended to write to NULOUT (see SU0YOMA)
#ifdef WITH_FCKIT
IF(.NOT.LLNEMOIOSERVER) THEN
  IF( FCKIT_MAIN%TASKID()+1 == 1 ) THEN
    CALL FCKIT_LOG%SET_FORTRAN_UNIT(UNIT=STDOUT,STYLE=FCKIT_LOG%PREFIX)
  ELSEIF( IO_SERV_C001%LIO_SERVER ) THEN
    WRITE(CLIOSERV_LOGFILE,'(A,I0)') "IOSERV.",IO_SERV_C001%MYPROC_IO
    CALL FCKIT_LOG%SET_FILE(TRIM(CLIOSERV_LOGFILE),STYLE=FCKIT_LOG%PREFIX)
  ELSE
    CALL FCKIT_LOG%RESET()
  ENDIF
ENDIF
#endif





FP_SERV_C001%LFP_SERVER=.FALSE.


LLSERV = IO_SERV_C001%LIO_SERVER .OR. FP_SERV_C001%LFP_SERVER

! Create a communicator for IFS compute tasks and NEMO/IO server

IF (LLSERV) THEN
  ICOL = 2
  INUM = 0
ELSE
  ICOL = 1
  CALL MPI_COMM_RANK (MPI_COMM_WORLD, INUM, IERR)
ENDIF
CALL MPI_COMM_SPLIT (MPI_COMM_WORLD, ICOL, INUM, ICOMM_IFS_PLUS_NEMO_IO_SERV, IERR)

IF (.NOT.LLSERV .AND. LLNEMOIO) THEN
#ifdef WITH_NEMO
  CALL ININEMOIO2 (ICOMM_IFS_PLUS_NEMO_IO_SERV)
#endif
ENDIF

IF (.NOT.LLNEMOIOSERVER) THEN
IF (LLSTARTUPCOST) THEN
   WRITE(0,'(1X,A,F30.6,A)') 'MAIN: MPI startup cost = ',ZMPI_INIT(2) - ZMPI_INIT(1),' secs'
   !WRITE(0,'(1X,A,3F30.9,A)') 'MAIN: MPI startup cost = ',ZMPI_INIT(1),ZMPI_INIT(2),ZMPI_INIT(2) - ZMPI_INIT(1),' secs'
ENDIF

LHOOK = LLHOOK_SAVE

IF (IO_SERV_C001%LIO_SERVER) THEN

  ! Set number of threads for IO servers
  ! Can't use OML functions here as OML_INIT has not yet been called, but needs to be before MPL_INIT
  CALL OMP_SET_NUM_THREADS(IO_SERV_C001%NTHREAD_IO)

  CALL EC_MEMINFO(-1,"master:io-server",ICOMM_IFS_PLUS_IFS_IO_SERV,KBARR=1,KIOTASK=1,KCALL=0)
  CALL IO_SERV_RUN(IO_SERV_C001)

ELSE

  CALL MPL_INIT() ! Don't rely on DR_HOOK or anything implicit to setup MPL_INIT

  IF (LHOOK) CALL DR_HOOK('MASTER',0,ZHOOK_HANDLE)
  CALL INIT_STACK(1)

  CALL EC_MEMINFO(-1,"master:computation",ICOMM_IFS_PLUS_IFS_IO_SERV,KBARR=1,KIOTASK=0,KCALL=0)
  CALL CNT0(LDCOUPACTIVE=LLCOUPACTIVE,KCOMM=ICOMM_IFS_PLUS_NEMO_IO_SERV)

  IF (LHOOK) CALL DR_HOOK('MASTER',1,ZHOOK_HANDLE)

ENDIF

!  Every program needs to be finalised (optional, but good practise)
#if defined(WITH_FCKIT)
CALL FCKIT_MAIN%FINALISE()
#endif

! OASIS3 or OASIS4 interface must be finalised after last DR_HOOK call.
#if defined(WITH_OASIS)
CALL CPL_FINALIZE(LLCOUPACTIVE)
#endif

!  CALL COUPLO4_ENDMPI

CALL MPL_END(KERROR=IER) ! Does not fail

#ifdef WITH_NEMO
IF (LLNEMOIO) CALL ENDNEMOIO()
#endif

!#ifdef 1
! XIOS and MPI finalization
! #ifdef WITH_XIOS
! CALL SUXIOS_FIN
! #elif WITH_CPLNG
! CALL CPLNG_FINALIZE
! #endif
!#endif

ENDIF

END PROGRAM MASTER
