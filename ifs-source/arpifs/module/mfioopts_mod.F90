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

MODULE MFIOOPTS_MOD

!**** *MFIOOPTS_MOD*  - A structure holding IO options

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE*
!      Original : 11-09-2012

!     Modifications :
!     P. Marguinaud : 10-10-2013 : Remove unused options
!     P. Marguinaud : 10-10-2014 : Add LROOO (read fields in the same order as
!                                  in the input file)

USE PARKIND1, ONLY : JPRB, JPIM
USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE MFIOOPTS
  LOGICAL             :: LDOCP         = .FALSE. ! Compress data
  LOGICAL             :: LDOGF         = .FALSE. ! Gather fields on #1
  LOGICAL             :: LDOWR         = .FALSE. ! Write data
  LOGICAL             :: LDOW1         = .FALSE. ! #1 writes data
  LOGICAL             :: LDOWN         = .FALSE. ! All procs write to their own file (1 file/forecast term)
  LOGICAL             :: LDOCW         = .FALSE. ! All procs write to the same file
  LOGICAL             :: LUCFA         = .FALSE. ! Use C layer for writing FA files
  LOGICAL             :: LUSE_IO_SERV  = .FALSE. ! Use IO server
  LOGICAL             :: LFACPT_OPENMP = .FALSE. ! Write with FA & OpenMP
  LOGICAL             :: LFADPT_OPENMP = .FALSE. ! Read with FA & OpenMP
  LOGICAL             :: LRECV_NOORDER = .FALSE. ! MPI #1 receives messages with order
  LOGICAL             :: LFACRW        = .FALSE. ! LFI read/write in C
  INTEGER (KIND=JPIM) :: NIOROOT       = 1_JPIM  ! Number of IO roots
  INTEGER (KIND=JPIM) :: NFACTM        = 0_JPIM  ! "Facteur multiplicatif"
  LOGICAL             :: LFILN         = .FALSE. ! Files are numbered by MPI rank
  LOGICAL             :: LROOO         = .FALSE. ! Adapt to read files with fields in the wrong order
END TYPE MFIOOPTS

TYPE MFIOFLAG
  LOGICAL             :: LLWRSFX       = .FALSE.
  LOGICAL             :: LLWRSPECA     = .FALSE.
  LOGICAL             :: LLWRGRIDA     = .FALSE.
  LOGICAL             :: LLWRGRIDUA    = .FALSE.
  LOGICAL             :: LLWRXFU       = .FALSE.
  LOGICAL             :: LLWRFU        = .FALSE.
END TYPE MFIOFLAG

PRIVATE

PUBLIC :: MFIOOPTS, MFIOOPTS_GETOPTS, &
        & MFIOFLAG, MFIOOPTS_GETFLAG

SAVE

CONTAINS

SUBROUTINE MFIOOPTS_GETFLAG(YDCFU,YDXFU,YGFL,YDPHY,YDFLAG,CDCONF,KSTEP)

USE YOMCT0, ONLY : NCONF, LELAM, LFBDAP, LOBSC1
USE YOMCT1, ONLY : LWRSPEC
USE YOMCFU, ONLY : TCFU
USE YOMPHY, ONLY : TPHY
USE YOMXFU, ONLY : TXFU
USE YOM_YGFL, ONLY : TYPE_GFLD

TYPE(TCFU)          ,INTENT(INOUT):: YDCFU
TYPE(TPHY)          ,INTENT(INOUT):: YDPHY
TYPE(TXFU)          ,INTENT(INOUT):: YDXFU
TYPE(TYPE_GFLD)     ,INTENT(INOUT):: YGFL
TYPE (MFIOFLAG),     INTENT (OUT) :: YDFLAG
CHARACTER (LEN=1),   INTENT (IN)  :: CDCONF
INTEGER (KIND=JPIM), INTENT (IN)  :: KSTEP

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK ('MFIOOPTS_MOD:MFIOOPTS_GETFLAG',0,ZHOOK_HANDLE)

YDFLAG%LLWRSFX = (CDCONF == 'x')

YDFLAG%LLWRSPECA = ((CDCONF == 'A').OR.(CDCONF == 'f').OR.(CDCONF == 'F').OR.(CDCONF == 'D')&
                 & .OR.(CDCONF == 'Q')).AND.LWRSPEC

YDFLAG%LLWRGRIDA = ((YDPHY%LMPHYS.AND.(CDCONF == 'A'.OR.CDCONF == 'f'.OR.CDCONF == 'F'.OR.CDCONF == 'Q')&
                 & .OR.LOBSC1.OR.NCONF == 131).OR.&
                 & (LELAM.AND.NCONF /= 601)).AND.(.NOT.(CDCONF == 'x'))

YDFLAG%LLWRGRIDUA = (YGFL%NUMGPFLDS > 0).AND.(.NOT.(CDCONF == 'x'))

YDFLAG%LLWRXFU = YDXFU%LXFU.AND.LFBDAP.AND.(CDCONF == 'A'.OR.CDCONF == 'F'.OR.CDCONF == 'f'.OR.CDCONF == 'Q') 

YDFLAG%LLWRFU = YDCFU%LCUMFU.AND.LFBDAP.AND.(CDCONF == 'A'.OR.CDCONF == 'F'.OR.CDCONF == 'f'.OR.CDCONF == 'Q').AND.(KSTEP > 0) 

IF (LHOOK) CALL DR_HOOK ('MFIOOPTS_MOD:MFIOOPTS_GETFLAG',1,ZHOOK_HANDLE)

END SUBROUTINE MFIOOPTS_GETFLAG

SUBROUTINE MFIOOPTS_GETOPTS (YDOPTS, YDIOS)

! Reckon IO options from NDISTIO and other parameters

USE YOMMP0, ONLY : NDISTIO
USE YOMIO_SERV, ONLY : IO_SERV

TYPE (MFIOOPTS), INTENT (INOUT)        :: YDOPTS
TYPE (IO_SERV),  INTENT (IN), OPTIONAL :: YDIOS

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK ('MFIOOPTS_MOD:MFIOOPTS_GETOPTS',0,ZHOOK_HANDLE)

IF (PRESENT (YDIOS)) THEN
  YDOPTS%LUSE_IO_SERV = YDIOS%LIO_SERV_WR
ELSE
  YDOPTS%LUSE_IO_SERV = .FALSE.
ENDIF

YDOPTS%LDOCP = (.NOT. YDOPTS%LUSE_IO_SERV) .AND. (NDISTIO(4) == 0) 
YDOPTS%LDOGF = (.NOT. YDOPTS%LUSE_IO_SERV) .AND. (NDISTIO(3) == 0) .AND. YDOPTS%LDOCP .AND. (NDISTIO(1) == 0)
YDOPTS%LDOWR = (.NOT. YDOPTS%LUSE_IO_SERV) .AND. (NDISTIO(2) == 0) .AND. YDOPTS%LDOCP .AND. (NDISTIO(1) /= 3)
YDOPTS%LDOCW = (.NOT. YDOPTS%LUSE_IO_SERV) .AND. (NDISTIO(1) == 3) .AND. YDOPTS%LDOCP

YDOPTS%LDOW1 = (.NOT. YDOPTS%LUSE_IO_SERV) .AND. (NDISTIO(1) == 0)
YDOPTS%LDOWN = (.NOT. YDOPTS%LUSE_IO_SERV) .AND. (NDISTIO(1) == 1) .OR.  YDOPTS%LDOCW
YDOPTS%LUCFA = (.NOT. YDOPTS%LUSE_IO_SERV) .AND. (NDISTIO(8) /= 0)

YDOPTS%LFACPT_OPENMP = (.NOT. YDOPTS%LUSE_IO_SERV) .AND. (NDISTIO(6) /= 0)
YDOPTS%LFADPT_OPENMP = (.NOT. YDOPTS%LUSE_IO_SERV) .AND. (NDISTIO(12) /= 0)
YDOPTS%LRECV_NOORDER = (.NOT. YDOPTS%LUSE_IO_SERV) .AND. (NDISTIO(7) /= 0)

YDOPTS%NIOROOT = MAX (NDISTIO(11), 1)
YDOPTS%LFACRW  = NDISTIO (9) /= 0
YDOPTS%NFACTM  = NDISTIO (10)

YDOPTS%LFILN = (YDOPTS%LDOWN .AND. .NOT. YDOPTS%LDOCW) .OR. (YDOPTS%NIOROOT > 1)
YDOPTS%LROOO = NDISTIO (13) /= 0


IF (LHOOK) CALL DR_HOOK ('MFIOOPTS_MOD:MFIOOPTS_GETOPTS',1,ZHOOK_HANDLE)

END SUBROUTINE MFIOOPTS_GETOPTS

END MODULE MFIOOPTS_MOD

