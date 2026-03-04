! (C) Copyright 1989- Meteo-France.

SUBROUTINE SCAN2M_HPOS(YDRQPHY,YDRQCLI,YDNAMFPSCI,YDAFN,LDCTLCLIM,KSTGLO,KBL,KEND,LDFPOSHOR,YDGEOMETRY,YDSURF,YDCFU,YDXFU, &
 & YDMODEL,PCFUBUF,PXFUBUF,PWSXI,PWDXI,PBUF)

!****-------------------------------------------------------------------
!**** *SCAN2M_HPOS* - Fullpos grid-point space computations
!****-------------------------------------------------------------------
!     Purpose.   post-processing in grid-point space
!     --------   

!**   Interface.
!     ----------
!        *CALL* *SCAN2M_HPOS (..)

!        Explicit arguments :
!        ------------------  
!          LDCTLCLIM - control input clim against model

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   See includes below.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Ryad El Khatib *Meteo-France* 

! Modifications
! -------------
!   original 18-Jul-2012 from SCAN2M
!   R. El Khatib  04-Dec-2012 Fix bounds checking issues created by recent "cleanings"
!   G. Balsamo 14-Jan-2014 add lake prognostics SP_SL
!   R. El Khatib 13-Dec-2012 Fullpos buffers reshaping
!   T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!   R. El Khatib 27-Jul-2016 interpolations over C+I+E
!   R. El Khatib 20-Sep-2016 Split HPOS
!   E.Dutra/G.Arduini Jan 2018: SP_SG 4 Dimensions, snow multi-layer 
!-----------------------------------------------------------------------------

USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMCFU             , ONLY : TCFU
USE YOMXFU             , ONLY : TXFU
USE TYPE_MODEL         , ONLY : MODEL
USE YOMAFN             , ONLY : TAFN
USE YOMFP4L            , ONLY : TRQFP
USE YOMFPC             , ONLY : TNAMFPSCI

!-----------------------------------------------------------------------------

IMPLICIT NONE

TYPE (TRQFP), INTENT(IN) :: YDRQPHY
TYPE (TRQFP), INTENT(IN) :: YDRQCLI
TYPE(TNAMFPSCI)   ,INTENT(IN) :: YDNAMFPSCI
TYPE(TAFN)        ,INTENT(IN) :: YDAFN
LOGICAL           ,INTENT(IN) :: LDCTLCLIM
INTEGER(KIND=JPIM),INTENT(IN) :: KSTGLO
INTEGER(KIND=JPIM),INTENT(IN) :: KBL
INTEGER(KIND=JPIM),INTENT(IN) :: KEND
LOGICAL           ,INTENT(IN) :: LDFPOSHOR
TYPE(GEOMETRY)    ,INTENT(IN) :: YDGEOMETRY
TYPE(TSURF)       ,INTENT(IN) :: YDSURF
TYPE(TCFU)        ,INTENT(IN) :: YDCFU
TYPE(TXFU)        ,INTENT(IN) :: YDXFU
TYPE(MODEL)       ,INTENT(IN) :: YDMODEL
REAL(KIND=JPRB)   ,INTENT(IN) :: PCFUBUF(YDGEOMETRY%YRDIM%NPROMA,YDCFU%NFDCFU,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(IN) :: PXFUBUF(YDGEOMETRY%YRDIM%NPROMA,YDXFU%NFDXFU,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(IN) :: PWSXI
REAL(KIND=JPRB)   ,INTENT(IN) :: PWDXI
REAL(KIND=JPRB)   ,INTENT(OUT) :: PBUF(YDGEOMETRY%YRDIM%NPROMA,YDRQPHY%NFIELDG)


INTEGER(KIND=JPIM) :: JFLD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "checkinclim.intfb.h"
#include "hpos.intfb.h"
#include "hpos_cfu.intfb.h"
#include "hpos_xfu.intfb.h"

!-----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SCAN2M_HPOS',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, &
 & YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, YDPHY1=>YDMODEL%YRML_PHY_MF%YRPHY1)
ASSOCIATE(NFIELDG=>YDRQPHY%NFIELDG,GFP_PHYDS=>YDAFN%GFP_PHYDS,NPROMA=>YDDIM%NPROMA, NGPTOT=>YDGEM%NGPTOT, NGPBLKS=>YDDIM%NGPBLKS, &
 & NFPCLI=>YDNAMFPSCI%NFPCLI, NFLEVG=>YDDIMV%NFLEVG, YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, &
 & SD_SM=>YDSURF%SD_SM, SD_VA=>YDSURF%SD_VA, SD_VC=>YDSURF%SD_VC, SD_VD=>YDSURF%SD_VD, SD_VF=>YDSURF%SD_VF, SD_VP=>YDSURF%SD_VP, &
 & SD_VV=>YDSURF%SD_VV, SD_VX=>YDSURF%SD_VX, SD_WS=>YDSURF%SD_WS, SD_X2=>YDSURF%SD_X2, SP_CI=>YDSURF%SP_CI, SP_RR=>YDSURF%SP_RR, &
 & SD_OC=>YDSURF%SD_OC, &
 & SP_SB=>YDSURF%SP_SB, SP_SG=>YDSURF%SP_SG, SP_SL=>YDSURF%SP_SL)


IF (LDCTLCLIM) THEN
  CALL CHECKINCLIM(YDNAMFPSCI,NPROMA,YDSURF,KEND,YDGEOMETRY%YROROG(KBL)%OROG,SD_VX(:,:,KBL),SD_VV(:,:,KBL))
ENDIF
CALL HPOS(YDRQPHY,YDRQCLI,YDNAMFPSCI,YDAFN,LDFPOSHOR,YDGEOMETRY,YDSURF,YDMODEL,KEND,KSTGLO,&
 & SP_SB(:,:,:,KBL),SP_SG(:,:,:,KBL),SP_SL(:,:,KBL),SP_RR(:,:,KBL),SP_CI(:,:,KBL),&
 & SD_WS(:,:,KBL),SD_VD(:,:,KBL),SD_VX(:,:,KBL),SD_VF(:,:,KBL),SD_VV(:,:,KBL),&
 & SD_VP(:,:,KBL),SD_VA(:,:,KBL),SD_VC(:,:,KBL),SD_X2(:,:,KBL),SD_SM(:,:,:,KBL),&
 & SD_OC(:,:,KBL), &
 & PWSXI,PWDXI,PBUF)

IF (YDCFU%NFDCFU > 0) THEN
  DO JFLD=1,NFIELDG
    CALL HPOS_CFU(YDAFN,LDFPOSHOR,YDCFU,NPROMA,NFLEVG,KEND,YDRQPHY%ICOD(JFLD),NFPCLI,PCFUBUF(:,:,KBL),PBUF(:,JFLD))
  ENDDO
ENDIF

IF (YDXFU%NFDXFU > 0) THEN
  DO JFLD=1,NFIELDG
    CALL HPOS_XFU(YDAFN,LDFPOSHOR,YDXFU,YDPHY,YDPHY1,NPROMA,NFLEVG,YDSURF,KEND,YDRQPHY%ICOD(JFLD),NFPCLI, &
     & YDGEOMETRY%YRGSGEOM(KBL)%GM,PWSXI,SP_RR(:,:,KBL),SD_VX(:,:,KBL),SD_VF(:,:,KBL),SD_VV(:,:,KBL), &
     & PXFUBUF(:,:,KBL),PBUF(:,JFLD))
  ENDDO
ENDIF

! ---------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SCAN2M_HPOS',1,ZHOOK_HANDLE)
END SUBROUTINE SCAN2M_HPOS
