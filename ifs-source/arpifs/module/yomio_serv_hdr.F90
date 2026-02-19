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

MODULE YOMIO_SERV_HDR

!**** *YOMIO_SERV_HDR*  - Defines field header for io_serv messages

!     Author. 
!     ------- 
!      Philippe Marguinaud  *METEO-FRANCE*
!      Original : 01-01-2011

!      Modifications :
!      P.Marguinaud : 11-09-2012 : Refactor using IOFLDDESC_MOD
!      P.Marguinaud : 10-10-2013 : Freeze & thaw routines
!      P.Marguinaud : 10-10-2014 : Add NIO_SERV_HDR_TYPE_WFL (write field)
!                                  and NIO_SERV_HDR_TYPE_RFL (read field)

USE PARKIND1, ONLY : JPIM, JPRB
USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM), PARAMETER :: NIO_SERV_HDR_TYPE_NON = 0, &
                               & NIO_SERV_HDR_TYPE_WFL = 1, &  ! Write field
                               & NIO_SERV_HDR_TYPE_RFL = 2, &  ! Read field
                               & NIO_SERV_HDR_TYPE_STO = 3, &  ! Stop
                               & NIO_SERV_HDR_TYPE_CLO = 4     ! Close file

INTEGER(KIND=JPIM), PARAMETER :: NIO_SERV_HDR_IDOM_FPA = 99, & ! All Fullpos domains (full)
                               & NIO_SERV_HDR_IDOM_MDL = 0 , & ! Model domain (full)
                               & NIO_SERV_HDR_IDOM_MDP = 1 , & ! Partial model domain
                               & NIO_SERV_HDR_IDOM_WV  = 2     !  Wave model domain (full)

!                                                      12345678
CHARACTER(LEN=8),   PARAMETER :: CIO_SERV_HDR_MAGIC = "IO_SERV_"

TYPE IO_SERV_HDR1
  SEQUENCE
  CHARACTER(LEN=8)   :: CMAGIC = CIO_SERV_HDR_MAGIC !  Magic string
! Field specific
! General
  INTEGER(KIND=JPIM) :: NTYPE = -1_JPIM             !  NIO_SERV_HDR_TYPE_WFL (write field) 
                                                    !  NIO_SERV_HDR_TYPE_RFL (read field)
                                                    !  NIO_SERV_HDR_TYPE_STO (stop)
                                                    !  NIO_SERV_HDR_TYPE_CLO (close)
  INTEGER(KIND=JPIM) :: ISTEP   = 0_JPIM            !  Model step
  INTEGER(KIND=JPIM) :: IDATEF (22) = &             !  Date
& (/ 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, &
&    0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM /)
  INTEGER(KIND=JPIM) :: IDOM    = 0_JPIM            !  Domain, 0=model, 99=Fullpos domains
  INTEGER(KIND=JPIM) :: IHH     = 0_JPIM            !  Hour   (for output file)
  INTEGER(KIND=JPIM) :: IMM     = 0_JPIM            !  Minute (for output file)
  INTEGER(KIND=JPIM) :: ISS     = 0_JPIM            !  Second (for output file)
  INTEGER(KIND=JPIM) :: MYPROC  = 0_JPIM            !  Origin task rank
  INTEGER(KIND=JPIM) :: IPAD    = 0_JPIM
END TYPE IO_SERV_HDR1

TYPE IO_SERV_HDR2
  SEQUENCE
  CHARACTER(LEN=8)   :: CMAGIC    = CIO_SERV_HDR_MAGIC !  Magic string
  INTEGER(KIND=JPIM) :: NFLDS     =  0_JPIM            !  Number of fields to be send
  INTEGER(KIND=JPIM) :: NSIZEG    =  0_JPIM            !  Global field size
  INTEGER(KIND=JPIM) :: NSIZEL    =  0_JPIM            !  Local field size
  INTEGER(KIND=JPIM) :: MYSETW    =  0_JPIM            !  YOMMP
  INTEGER(KIND=JPIM) :: MYSETV    =  0_JPIM            !  YOMMP
  INTEGER(KIND=JPIM) :: MYPROC    =  0_JPIM            !  YOMMP
  INTEGER(KIND=JPIM) :: MYSETA    =  0_JPIM            !  YOMMP
  INTEGER(KIND=JPIM) :: MYSETB    =  0_JPIM            !  YOMMP
  INTEGER(KIND=JPIM) :: IHASFLDSC =  0_JPIM            !  Field headers present
  INTEGER(KIND=JPIM) :: ISSPEC    =  0_JPIM            !  Spectral field
  INTEGER(KIND=JPIM) :: ITAG      =  0_JPIM            !  Tag
  INTEGER(KIND=JPIM) :: NPETOT    =  0_JPIM            !  Number of processing elements involved
                                                       !  to create this field
  INTEGER(KIND=JPIM) :: NPELOC    =  0_JPIM            !  Rank of processing element
  INTEGER(KIND=JPIM) :: NPESET    =  0_JPIM            !  Set of processing element
END TYPE IO_SERV_HDR2

TYPE IO_SERV_HDRG
  SEQUENCE
  TYPE (IO_SERV_HDR1) :: YH1
  TYPE (IO_SERV_HDR2) :: YH2
  INTEGER(KIND=JPIM)  :: IPAD      =  0_JPIM        !  Skip words after header
  INTEGER(KIND=JPIM)  :: ISPARE    =  0_JPIM
END TYPE IO_SERV_HDRG

INTEGER(KIND=JPIM) :: NIO_SERV_HDRF_SIZE = 0, &
                    & NIO_SERV_HDRG_SIZE = 0, &
                    & NIO_SERV_HDR1_SIZE = 0, &
                    & NIO_SERV_HDR2_SIZE = 0

SAVE

CONTAINS

#define LHOOK .FALSE.
REAL (KIND=JPRB) FUNCTION IO_SERV_FREEZE_HDRG (YDHG)
DIMENSION :: IO_SERV_FREEZE_HDRG (NIO_SERV_HDRG_SIZE)
TYPE (IO_SERV_HDRG), INTENT (IN) :: YDHG
REAL (KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK ('YOMIO_SERV_HDR:IO_SERV_FREEZE_HDRG',0,ZHOOK_HANDLE)
IO_SERV_FREEZE_HDRG = TRANSFER (YDHG, IO_SERV_FREEZE_HDRG)
IF (LHOOK) CALL DR_HOOK ('YOMIO_SERV_HDR:IO_SERV_FREEZE_HDRG',1,ZHOOK_HANDLE)
END FUNCTION IO_SERV_FREEZE_HDRG

TYPE (IO_SERV_HDRG) FUNCTION IO_SERV_THAW_HDRG (P)
REAL (KIND=JPRB), INTENT (IN) :: P (NIO_SERV_HDRG_SIZE)
REAL (KIND=JPHOOK) :: ZHOOK_HANDLE
#include "abor1.intfb.h"
IF (LHOOK) CALL DR_HOOK ('YOMIO_SERV_HDR:IO_SERV_THAW_HDRG',0,ZHOOK_HANDLE)
IO_SERV_THAW_HDRG = TRANSFER (P, IO_SERV_THAW_HDRG)
IF (IO_SERV_THAW_HDRG%YH1%CMAGIC /= CIO_SERV_HDR_MAGIC) &
  & CALL ABOR1 ('YOMIO_SERV_HDR:IO_SERV_THAW_HDRG: DETECTED CORRUPTED IO_SERV_HDRG')
IF (LHOOK) CALL DR_HOOK ('YOMIO_SERV_HDR:IO_SERV_THAW_HDRG',1,ZHOOK_HANDLE)
END FUNCTION IO_SERV_THAW_HDRG

REAL (KIND=JPRB) FUNCTION IO_SERV_FREEZE_HDR1 (YDH1)
DIMENSION :: IO_SERV_FREEZE_HDR1 (NIO_SERV_HDR1_SIZE)
TYPE (IO_SERV_HDR1), INTENT (IN) :: YDH1
REAL (KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK ('YOMIO_SERV_HDR:IO_SERV_FREEZE_HDR1',0,ZHOOK_HANDLE)
IO_SERV_FREEZE_HDR1 = TRANSFER (YDH1, IO_SERV_FREEZE_HDR1)
IF (LHOOK) CALL DR_HOOK ('YOMIO_SERV_HDR:IO_SERV_FREEZE_HDR1',1,ZHOOK_HANDLE)
END FUNCTION IO_SERV_FREEZE_HDR1

TYPE (IO_SERV_HDR1) FUNCTION IO_SERV_THAW_HDR1 (P)
REAL (KIND=JPRB), INTENT (IN) :: P (NIO_SERV_HDR1_SIZE)
REAL (KIND=JPHOOK) :: ZHOOK_HANDLE
#include "abor1.intfb.h"
IF (LHOOK) CALL DR_HOOK ('YOMIO_SERV_HDR:IO_SERV_THAW_HDR1',0,ZHOOK_HANDLE)
IO_SERV_THAW_HDR1 = TRANSFER (P, IO_SERV_THAW_HDR1)
IF (IO_SERV_THAW_HDR1%CMAGIC /= CIO_SERV_HDR_MAGIC) &
  & CALL ABOR1 ('YOMIO_SERV_HDR:IO_SERV_THAW_HDR1: DETECTED CORRUPTED IO_SERV_HDR1')
IF (LHOOK) CALL DR_HOOK ('YOMIO_SERV_HDR:IO_SERV_THAW_HDR1',1,ZHOOK_HANDLE)
END FUNCTION IO_SERV_THAW_HDR1
#undef LHOOK

END MODULE YOMIO_SERV_HDR

