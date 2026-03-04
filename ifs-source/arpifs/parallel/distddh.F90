! (C) Copyright 1989- Meteo-France.

SUBROUTINE DISTDDH(PCOMM,LDFLAG)

!**** *DISTDDH* Communicate and calculate global max.
!               Distance of single points for sumddh

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *DISTDDH(...)

!        Explicit arguments : PCOMM       - things to communicate
!        -------------------- LDFLAG      - Flag processor that owns
!                                          global naximum 

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      MPP Group *ECMWF*
!      Original: 95-10-01 

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : NPRCIDS, MYPROC, NPROC
USE YOMTAG   , ONLY : MTAGDDH1
USE YOMLUN   , ONLY : NULERR
USE MPL_MODULE

IMPLICIT NONE

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCOMM(7) 
LOGICAL           ,INTENT(OUT)   :: LDFLAG 
REAL(KIND=JPRB) ::ZCOMMIN(7),ZCOMMOUT(8)

INTEGER(KIND=JPIM) :: IBFLEN, IMSGLEN, ITAG1, ITAG2, JROC, IMAXPROC

REAL(KIND=JPRB) :: ZDMAX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('DISTDDH',0,ZHOOK_HANDLE)
IBFLEN=SIZE(ZCOMMIN)
LDFLAG=.FALSE.

ITAG1=MTAGDDH1

IF (MYPROC == 1) THEN
  ZDMAX=PCOMM(1)
  ZCOMMOUT(1:7)=PCOMM(:)
  ZCOMMOUT(8)=NPRCIDS(1)
  IMAXPROC=NPRCIDS(1)
  DO JROC=2,NPROC
    CALL MPL_RECV(ZCOMMIN(1:IBFLEN),KSOURCE=NPRCIDS(JROC),KTAG=ITAG1, &
     & KOUNT=IMSGLEN,CDSTRING='DISTDDH:')  
    IF (IMSGLEN /= IBFLEN) THEN
      WRITE(NULERR,*)' MASTERRECV: IMSGLEN,IBFLEN ',IMSGLEN,IBFLEN
      CALL ABOR1('DISTDDH: RECEIVED'//' MESSAGE LENGTH OF WRONG SIZE')
    ENDIF

    IF(ZCOMMIN(1) > ZDMAX) THEN
      ZDMAX=ZCOMMIN(1)
      ZCOMMOUT(1:7)=ZCOMMIN(:)
!       record in ZCOMMOUT(8) the PE number which owns the point which is
!       closest to the DDH point.
!       If more than one PE has the same distance to a DDH point,
!       only the first such processor is selected.
      ZCOMMOUT(8)=NPRCIDS(JROC)
      IMAXPROC=NPRCIDS(JROC)
    ENDIF
  ENDDO

ELSE
  ZCOMMIN(:)=PCOMM(:)
  CALL MPL_SEND(ZCOMMIN(1:IBFLEN),KDEST=NPRCIDS(1),KTAG=ITAG1, &
   & CDSTRING='DISTDDH:')  
ENDIF

! The global max is now available on the root node. Now it will be
! broadcasted back.

IBFLEN=SIZE(ZCOMMOUT)
ITAG2=MTAGDDH1+1
IF (MYPROC == 1) THEN
  DO JROC=2,NPROC
    CALL MPL_SEND(ZCOMMOUT(1:IBFLEN),KDEST=NPRCIDS(JROC),KTAG=ITAG2, &
     & CDSTRING='DISTDDH:')  
  ENDDO
ELSE
  CALL MPL_RECV(ZCOMMOUT(1:IBFLEN),KSOURCE=NPRCIDS(1),KTAG=ITAG2, &
   & KOUNT=IMSGLEN,CDSTRING='DISTDDH:')  
  IF (IMSGLEN /= IBFLEN) CALL ABOR1('DISTDDH : RECEIVED'//&
   & ' MESSAGE LENGTH OF WRONG SIZE')  
  ZDMAX=ZCOMMOUT(1)
  IMAXPROC=ZCOMMOUT(8)
ENDIF

PCOMM(:)=ZCOMMOUT(1:7)

!      LDFLAG = true for the processor that owns the point
!              false otherwise

IF(IMAXPROC == MYPROC) LDFLAG=.TRUE.

IF (LHOOK) CALL DR_HOOK('DISTDDH',1,ZHOOK_HANDLE)
END SUBROUTINE DISTDDH
