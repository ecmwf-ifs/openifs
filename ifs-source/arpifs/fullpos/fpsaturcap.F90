! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPSATURCAP(YDAFN,YDFPGEO,YDGEOMETRY,KAPTR,PABUF)

!**** *FPSATURCAP*  - Processing previously done in FPSPECFITG

!     PURPOSE.
!     --------
!        Cap humidity to saturation

!**   INTERFACE.
!     ----------
!       *CALL* *FPSATURCAP*

!        EXPLICIT ARGUMENTS
!        --------------------
!     PABUF    : post-processed fields 

!        IMPLICIT ARGUMENTS
!        --------------------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      Tomas Wilhelmsson *ECMWF*
!      ORIGINAL : 2014-06-10

!     MODIFICATIONS.
!     --------------
!      R. El Khatib 08-Jul-2022 Contribution to the encapsulation of YOMCST and YOETHF
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARFPOS  , ONLY : JPOSDYN
USE YOMFPGEO , ONLY : TFPGEO
USE YOMAFN   , ONLY : TAFN
USE YOMCST   , ONLY : YDCST=>YRCST ! allows use of included functions. REK.
USE YOETHF   , ONLY : YDTHF=>YRTHF ! allows use of included functions. REK.

!     ------------------------------------------------------------------
IMPLICIT NONE

TYPE(TAFN),         INTENT(IN)    :: YDAFN
TYPE(TFPGEO),       INTENT(IN)    :: YDFPGEO
TYPE(GEOMETRY),     INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN)    :: KAPTR(JPOSDYN)
REAL(KIND=JPRB),    INTENT(INOUT) :: PABUF(:,:,:)

INTEGER(KIND=JPIM) :: JBLOC, JL, IST, IEND, JI

REAL(KIND=JPRB) :: ZQS ! Saturation humidity at level
REAL(KIND=JPRB) :: ZPRESH(YDFPGEO%NFPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPRESF(YDFPGEO%NFPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gphpre.intfb.h"

#include "fcttre.func.h"

!      -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPSATURCAP',0,ZHOOK_HANDLE)
ASSOCIATE(TFP=>YDAFN%TFP, YDDIMV=>YDGEOMETRY%YRDIMV)
ASSOCIATE(NFPBLOCS=>YDFPGEO%NFPBLOCS, NFPROMA=>YDFPGEO%NFPROMA, NFPEND=>YDFPGEO%NFPEND, &
 & RETV=>YDCST%RETV, RLVTT=>YDCST%RLVTT, RLSTT=>YDCST%RLSTT, RTT=>YDCST%RTT, &
 & R2ES=>YDTHF%R2ES, R3LES=>YDTHF%R3LES, R3IES=>YDTHF%R3IES, R4LES=>YDTHF%R4LES, &
 & R4IES=>YDTHF%R4IES, R5LES=>YDTHF%R5LES, R5IES=>YDTHF%R5IES, R5ALVCP=>YDTHF%R5ALVCP, &
 & R5ALSCP=>YDTHF%R5ALSCP, RALVDCP=>YDTHF%RALVDCP, RALSDCP=>YDTHF%RALSDCP, &
 & RTWAT=>YDTHF%RTWAT, RTICE=>YDTHF%RTICE, RTICECU=>YDTHF%RTICECU, RTWAT_RTICE_R=>YDTHF%RTWAT_RTICE_R, &
 & RTWAT_RTICECU_R=>YDTHF%RTWAT_RTICECU_R, &
 & NFLEVG=>YDDIMV%NFLEVG)

!      -----------------------------------------------------------

! Cap specific humidity to saturation 
DO JBLOC=1,NFPBLOCS
  IST =1
  IEND=NFPEND(JBLOC)
    
  DO JI=IST,IEND
    ZPRESH(JI,NFLEVG) = PABUF(JI,KAPTR(TFP%LNSP%ICOD),JBLOC)
  ENDDO
  CALL GPHPRE(NFPROMA,NFLEVG,IST,IEND,YDGEOMETRY%YRVAB,YDGEOMETRY%YRCVER,ZPRESH,PRESF=ZPRESF)

  DO JL=1,NFLEVG
    DO JI=IST,IEND
      ZQS=FOEEWM(PABUF(JI,KAPTR(TFP%T%ICOD)+JL-1,JBLOC))/ZPRESF(JI,JL)
      ZQS=MIN(0.5_JPRB,ZQS)
      ZQS=ZQS/(1.0_JPRB-RETV*ZQS)
      PABUF(JI,KAPTR(TFP%Q%ICOD)+JL-1,JBLOC)=MAX(0.0_JPRB,MIN(ZQS, &
        & PABUF(JI,KAPTR(TFP%Q%ICOD)+JL-1,JBLOC)))
    ENDDO
  ENDDO
ENDDO
    
!      -----------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FPSATURCAP',1,ZHOOK_HANDLE)
!      -----------------------------------------------------------

END SUBROUTINE FPSATURCAP
