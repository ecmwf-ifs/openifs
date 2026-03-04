! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE UVRADI_LAYER(YDSURF, &
  & YDMODEL,KDIM, PAUX, AUXL, GEMSL, SURFL, STATE_T0, STATE_TMP, PDIAG, PSURF, UVP1, CHEM)


!**** *UVRADI_LAYER* - Layer routine calling time stepping of surface scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! PAUX     : Derived variables for general auxiliary quantities
! AUXL     : Derived variables for local auxiliary quantities
! GEMSL    : Derived variables for local GEMS quantities
! SURFL    : Derived variables for local surface quantities
! state_t0 : Derived variable for actual time model state 
! state_tmp: Derived variable for updated  model state 
! PDIAG    : Derived variable for diagnostics quantities
! PSURF    : Derived variables for general surface quantities

!     ==== Input/output ====
! PGFL     : GFL structure


!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!      Original : 2012-12-04  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!      K. Yessad (July 2014): Move some variables.
!      ABozzo Sep 2016: included switch to pass O3 from chemistry 
!-----------------------------------------------------------------------

USE PARKIND1           , ONLY : JPIM, JPRB
USE TYPE_MODEL         , ONLY : MODEL
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE VARIABLE_MODULE    , ONLY : VARIABLE_3D
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMPHYDER          , ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
 &                              SURF_AND_MORE_TYPE, SURF_AND_MORE_LOCAL_TYPE, AUX_DIAG_LOCAL_TYPE, &
 &                              GEMS_LOCAL_TYPE, AUX_DIAG_TYPE
USE YOMCT3             , ONLY : NSTEP
USE YOMCT0             , ONLY : LIFSMIN, LIFSTRAJ
USE YOECLDP            , ONLY : NCLDQI, NCLDQL

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT (IN)   :: AUXL
TYPE (GEMS_LOCAL_TYPE)         , INTENT (IN)   :: GEMSL
TYPE (SURF_AND_MORE_LOCAL_TYPE), INTENT (IN)   :: SURFL
TYPE (STATE_TYPE),TARGET       , INTENT (IN)   :: STATE_T0
TYPE (STATE_TYPE),TARGET       , INTENT (IN)   :: STATE_TMP
TYPE (AUX_DIAG_TYPE)           , INTENT (IN)   :: PDIAG
TYPE (SURF_AND_MORE_TYPE)      , INTENT (INOUT):: PSURF
TYPE(VARIABLE_3D)              , INTENT(INOUT) :: UVP1
TYPE(VARIABLE_3D)              , INTENT(INOUT) :: CHEM(:)
! REAL(KIND=JPRB)                , INTENT(INOUT) :: PGFL(KDIM%KLON,KDIM%KLEV,YDMODEL%YRML_GCONF%YGFL%NDIM)
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JCAER,JCHEM,IO3_CHEM
INTEGER(KIND=JPIM) :: ITIME, IZT, JL, JK

REAL(KIND=JPRB) :: ZRMUZ
REAL(KIND=JPRB) :: ZUVC(KDIM%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NUV), ZUVT(KDIM%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NUV)

REAL(KIND=JPRB) :: ZO3_MMR(KDIM%KLON,KDIM%KLEV)

TYPE (STATE_TYPE), POINTER :: STATE
REAL(KIND=JPRB), POINTER :: ZPRES(:,:), ZPRESF(:,:), ZDELP(:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "uvradi.intfb.h"
#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('UVRADI_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YDEUVRAD=>YDMODEL%YRML_PHY_RAD%YREUVRAD,YDRIP=>YDMODEL%YRML_GCONF%YRRIP, &
 & YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, &
 & YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2,YGFL=>YDMODEL%YRML_GCONF%YGFL,YDDYNA=>YDMODEL%YRML_DYN%YRDYNA)
ASSOCIATE(NACTAERO=>YGFL%NACTAERO, NDIM=>YGFL%NDIM, YUVP=>YGFL%YUVP, &
 & LUVPOUT=>YGFL%LUVPOUT,LERADIMPL=>YDEPHY%LERADIMPL, &
 & NAER=>YDERAD%NAER, NSW=>YDERAD%NSW, NUV=>YDERAD%NUV, &
 & LUVDBG=>YDEUVRAD%LUVDBG, NUVTIM=>YDEUVRAD%NUVTIM, RMUZUV=>YDEUVRAD%RMUZUV, &
 & NSTART=>YDRIP%NSTART, YSD_VD=>YDSURF%YSD_VD, &
 & YSD_VF=>YDSURF%YSD_VF, YSP_RR=>YDSURF%YSP_RR,TSPHY=>YDPHY2%TSPHY, &
 & LO3_CHEM_UV=>YDEUVRAD%LO3_CHEM_UV, YO3=>YGFL%YO3)
!     ------------------------------------------------------------------

!*         0.  Do the job :-)   

! Define state and pressure
IF(LERADIMPL) THEN
  STATE=>STATE_TMP
  ZPRES=>PAUX%PRS1
  ZPRESF=>PAUX%PRSF1
  ZDELP=>PAUX%PDELP1
ELSE
  STATE=>STATE_T0
  ZPRES=>PAUX%PAPRS
  ZPRESF=>PAUX%PAPRSF
  ZDELP=>PAUX%PDELP
ENDIF

IF (NAER == 1) THEN
  JCAER=6
ELSE
  JCAER=0
ENDIF

ZRMUZ=MAX(0._JPRB, MAXVAL(PAUX%PMU0(KDIM%KIDIA:KDIM%KFDIA)))

ITIME=NINT(TSPHY)
IF (YDDYNA%LTWOTL) THEN
  IZT=NINT(TSPHY*(REAL(NSTEP,JPRB)+0.5_JPRB))
ELSE
  IZT=ITIME*NSTEP
ENDIF

IO3_CHEM=-999
IF (LO3_CHEM_UV) THEN  
! find prognostic ozone in chemistry array
  DO JCHEM=1,YGFL%NCHEM
    IF (TRIM(YGFL%YCHEM(JCHEM)%CNAME) == 'O3' )  IO3_CHEM = JCHEM
  ENDDO

  IF (IO3_CHEM > 0  ) THEN 
    DO JK=1,KDIM%KLEV
     DO JL=KDIM%KIDIA,KDIM%KFDIA
        ZO3_MMR(JL,JK)=CHEM(IO3_CHEM)%PH9(JL,JK)
     ENDDO
    ENDDO
  ELSE
    CALL ABOR1('UVRADI_LAYER: chemical ozone is not defined')
  ENDIF
ELSE
  IF (YO3%LGP) THEN  
    ZO3_MMR=STATE%O3
  ELSE
    CALL ABOR1('UVRADI_LAYER: IFS ozone is not defined')
  ENDIF
ENDIF

IF (ZRMUZ > RMUZUV .AND. IZT <= NUVTIM ) THEN

  LUVDBG=.FALSE.
  IF (NSTEP == NSTART) THEN
    LUVDBG=.TRUE.
  ENDIF

  CALL UVRADI &
    & ( YDMODEL, KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON , KDIM%KLEV , JCAER, NSW   , NUV  , &
    &   GEMSL%ZAEROP(:,:,1:NACTAERO),SURFL%ZALBD, SURFL%ZALBP, ZPRES, ZPRESF, AUXL%ZCCNL, AUXL%ZCCNO , STATE%A   , &
    &   PAUX%PGELAM,PAUX%PCLON, PAUX%PSLON, ZDELP, PAUX%PGEMU, PAUX%PMU0 , ZO3_MMR , STATE%Q , PDIAG%PQSAT, &
    &   STATE%CLD(:,:,NCLDQI)   , STATE%CLD(:,:,NCLDQL)   , PSURF%PSD_VF(:,YSD_VF%YLSM%MP) , STATE%T , &
    &   PSURF%PSP_RR(:,YSP_RR%YT%MP9)  , &
    &   ZUVC , ZUVT ,PSURF%PSD_VD(:,YSD_VD%YUVBEDCS%MP),PSURF%PSD_VD(:,YSD_VD%YUVBED%MP), &
    &   UVP1)

ELSEIF(LUVPOUT .AND. NACTAERO > 0 .AND. .NOT.LIFSMIN .AND. .NOT. LIFSTRAJ) THEN

  UVP1%P(KDIM%KIDIA:KDIM%KFDIA,1) = PAUX%PMU0(KDIM%KIDIA:KDIM%KFDIA)
  UVP1%P(KDIM%KIDIA:KDIM%KFDIA,2:2*NUV+6) = 0._JPRB

ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('UVRADI_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE UVRADI_LAYER
