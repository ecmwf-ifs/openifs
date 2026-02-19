SUBROUTINE SRTM_INIT(NWVCONTINUUM)

!-- read in the basic coefficients to configure RRTM_SW
!- creates module YOESRTWN with BG, NSPA, NSPB, WAVENUM1, WAVENUM2,
!  DELWAVE, PREF, PREFLOG, TREF

USE PARKIND1  ,ONLY : JPIM , JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

USE PARSRTM  , ONLY : JPG, JPSW
USE YOESRTM  , ONLY : NGN
USE YOESRTWN , ONLY : NG, NGM, WT, NGC, RWGT, WTSM
!USE YOESRTWN , ONLY : NG, NGM, WT, NGC, NGN, RWGT, WTSM
!USE YOMLUN   , ONLY : NULOUT

IMPLICIT NONE

! Water vapour continuum model (0=default MT_CKD2.5, 1=CAVIAR)
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NWVCONTINUUM

! Local variables
INTEGER(KIND=JPIM) :: IGC, IGCSM, IBND, IG, IND, IPR, IPRSM
REAL(KIND=JPRB)    :: ZWTSUM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!#include "susrtmcf.intfb.h"
#include "susrtm.intfb.h"
#include "srtm_kgb16.intfb.h"
#include "srtm_kgb17.intfb.h"
#include "srtm_kgb18.intfb.h"
#include "srtm_kgb19.intfb.h"
#include "srtm_kgb20.intfb.h"
#include "srtm_kgb21.intfb.h"
#include "srtm_kgb22.intfb.h"
#include "srtm_kgb23.intfb.h"
#include "srtm_kgb24.intfb.h"
#include "srtm_kgb25.intfb.h"
#include "srtm_kgb26.intfb.h"
#include "srtm_kgb27.intfb.h"
#include "srtm_kgb28.intfb.h"
#include "srtm_kgb29.intfb.h"
!#include "susrtop.intfb.h"

#include "modify_wv_continuum.intfb.h"

#include "srtm_cmbgb16.intfb.h"
#include "srtm_cmbgb17.intfb.h"
#include "srtm_cmbgb18.intfb.h"
#include "srtm_cmbgb19.intfb.h"
#include "srtm_cmbgb20.intfb.h"
#include "srtm_cmbgb21.intfb.h"
#include "srtm_cmbgb22.intfb.h"
#include "srtm_cmbgb23.intfb.h"
#include "srtm_cmbgb24.intfb.h"
#include "srtm_cmbgb25.intfb.h"
#include "srtm_cmbgb26.intfb.h"
#include "srtm_cmbgb27.intfb.h"
#include "srtm_cmbgb28.intfb.h"
#include "srtm_cmbgb29.intfb.h"

IF (LHOOK) CALL DR_HOOK('SRTM_INIT',0,ZHOOK_HANDLE)

!CALL SUSRTMCF
CALL SUSRTM

!-- read in the molecular absorption coefficients

CALL SRTM_KGB16
CALL SRTM_KGB17
CALL SRTM_KGB18
CALL SRTM_KGB19
CALL SRTM_KGB20
CALL SRTM_KGB21
CALL SRTM_KGB22
CALL SRTM_KGB23
CALL SRTM_KGB24
CALL SRTM_KGB25
CALL SRTM_KGB26
CALL SRTM_KGB27
CALL SRTM_KGB28
CALL SRTM_KGB29

IF (PRESENT(NWVCONTINUUM)) THEN
  ! Modify the shortwave water vapour continuum, if requested
  CALL MODIFY_WV_CONTINUUM(NWVCONTINUUM)
END IF

!-- read in the cloud optical properties
!- creates module YOESRTOP with EXTLIQ1, SSALIQ1, ASYLIQ1, 
!  EXTICE3, SSAICE3, ASYICE3, FDLICE3  

!-- RRTM_SW cloud optical properties are not used
!   SRTM_CLDPROP is not called
!   no need to call SUSRTOP

!CALL SUSRTOP ( -1 )


!Mike Iacono 20050804
!-- Perform g-point reduction from 16 per band (224 total points) to
!-- a band dependent number (112 total points) for all absorption
!-- coefficient input data and Planck fraction input data.
!-- Compute relative weighting for new g-point combinations.

IGCSM = 0
!WRITE(NULOUT,9001) JPSW,JPG
9001 FORMAT(1X,'srtm_init JPSW=',I3,' JPG=',I3)
DO IBND = 1,JPSW
  IPRSM = 0
!  WRITE(NULOUT,9002) IBND,NGC(IBND)
9002 FORMAT(1X,'srtm_init NGC(',I3,')=',I3)
  IF (NGC(IBND) < JPG) THEN
    DO IGC = 1,NGC(IBND)
      IGCSM = IGCSM + 1
      ZWTSUM = 0.
!      WRITE(NULOUT,9003) IGC,IGCSM,NGN(IGCSM)
9003  FORMAT(1X,'srtm_init IGC=',I3,' NGN(',I3,')=',I3)
      DO IPR = 1, NGN(IGCSM)
        IPRSM = IPRSM + 1
!        WRITE(NULOUT,9004) IPR,IPRSM,WT(IPRSM)
9004    FORMAT(1X,'srtm_init IPR=',I3,' WT(',I3,')=',E14.7)
        ZWTSUM = ZWTSUM + WT(IPRSM)
      ENDDO
!      WRITE(NULOUT,9005) IGC,ZWTSUM
9005  FORMAT(1X,'srtm_init WTSM(',I3,')=',E14.7)
      WTSM(IGC) = ZWTSUM
    ENDDO

!    WRITE(NULOUT,9006) IBND+15,NG(IBND+15)
9006 FORMAT(1X,'srtm_init NG(',I3,')=',I3)
    DO IG = 1,NG(IBND+15)
      IND = (IBND-1)*JPG + IG
!!      WRITE(NULOUT,9007) IND,NGM(IND), IG,WT(IG), WTSM(NGM(IND)), IND,RWGT(IND)
9007 FORMAT(1X,'srtm_init NGM(',I3,')=',I3,' WT(',I3,')=',E13.7,' WTSM=',E13.7,' RWGT(',I3,')=',E13.7)
      RWGT(IND) = WT(IG)/WTSM(NGM(IND))
!      WRITE(NULOUT,9007) IND,NGM(IND),IG,WT(IG),WTSM(NGM(IND)),IND,RWGT(IND)
    ENDDO
  ELSE
    DO IG = 1,NG(IBND+15)
      IGCSM = IGCSM + 1
      IND = (IBND-1)*JPG + IG
      RWGT(IND) = 1.0
    ENDDO
  ENDIF
ENDDO

CALL SRTM_CMBGB16
CALL SRTM_CMBGB17
CALL SRTM_CMBGB18
CALL SRTM_CMBGB19
CALL SRTM_CMBGB20
CALL SRTM_CMBGB21
CALL SRTM_CMBGB22
CALL SRTM_CMBGB23
CALL SRTM_CMBGB24
CALL SRTM_CMBGB25
CALL SRTM_CMBGB26
CALL SRTM_CMBGB27
CALL SRTM_CMBGB28
CALL SRTM_CMBGB29

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_INIT',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_INIT

