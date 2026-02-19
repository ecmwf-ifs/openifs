! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_PHY3_LAYER(YDSURF, &
 ! Input quantities
  & YDMODEL,KDIM, PAUX, STATE, SURFL, AUXL, PDIAG, PCHEM2AER, &
 ! Input/Output quantities
  & PGFL, PSURF, FLUX, GEMSL, PRAD)

!**** *AER_PHY3_LAYER* - Layer routine calling the last part of prognostics aerosol scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! PAUX     : Derived variables for general auxiliary quantities
! state    : Derived variable for model state
! SURFL    : Derived variable for local surface quantities.
! AUXL     : Derived variables for local auxiliary quantities
! PDIAG    : Derived variables for diagnostics quantities
! PCHEM2AER  (KLON,KLEV,6)     : TENDENCY OF Selected TRACERS because specific
!            reactions (kg/kg s-1)

!     ==== Input/output ====
! PGFL     : GFL arrays
! PSURF    : Derived variables for general surface quantities
! FLUX     : Derived variable for fluxes
! GEMSL    : Derived variables for local GEMS quantities


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
!      Original : 2012-12-05  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!     JJMorcrette 20131001 MACC-related diagnostics & DMS-related variables
!     E. Dutra/G.Arduini, Jan 2018: snow multi-layer: pre-compute total snow mass

!-----------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1  ,ONLY : JPIM,    JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST    ,ONLY : RA, RPI

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
   & AUX_DIAG_LOCAL_TYPE, SURF_AND_MORE_TYPE, FLUX_TYPE,  &
   & GEMS_LOCAL_TYPE, SURF_AND_MORE_LOCAL_TYPE, AUX_RAD_TYPE, AUX_DIAG_TYPE
USE YOECLDP   , ONLY : NCLDQL,  NCLDQI,  NCLDQR, NCLDQS
USE YOE_AERODIAG, ONLY : NPAERAOT, NPAERLISI_VAR, NPAERLISI_WVL, JPAERO_WVL_AOD
USE TM5_CHEM_MODULE    , ONLY : NCHEM2AER

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (SURF_AND_MORE_LOCAL_TYPE), INTENT (IN)   :: SURFL
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT (IN)   :: AUXL
TYPE (AUX_DIAG_TYPE)           , INTENT (IN)   :: PDIAG
REAL(KIND=JPRB)                , INTENT (IN)   :: PCHEM2AER(KDIM%KLON,KDIM%KLEV,NCHEM2AER)
REAL(KIND=JPRB)                , INTENT(INOUT) :: PGFL(KDIM%KLON,KDIM%KLEV,YDMODEL%YRML_GCONF%YGFL%NDIM)
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (FLUX_TYPE)               , INTENT(INOUT) :: FLUX
TYPE (GEMS_LOCAL_TYPE)         , INTENT(INOUT) :: GEMSL
TYPE (AUX_RAD_TYPE)            , INTENT(INOUT) :: PRAD
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JAER, JVAR, JWVL, JGFL, JK, JL
REAL(KIND=JPRB) :: ZTH(KDIM%KLON,KDIM%KLEV+1),ZSNM(KDIM%KLON)
REAL(KIND=JPRB) :: ZAERO_WVL_DIAG(KDIM%KLON,YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG,YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG_TYPES)
REAL(KIND=JPRB),ALLOCATABLE :: ZAREA(:)

INTEGER(KIND=JPIM) :: IIFSCHEM ! for looping over wavelengths wrt RADAER code
INTEGER(KIND=JPIM) :: IWETOX_ON, ISO2WETOXBYO3, IDRYOX_IN_AER, IWETOX_IN_AER
!options set in ukca_chem_ifs and passed to ukca_aero_step for chemistry
REAL(KIND=JPRB) :: ZDELSO2_DRY_OH(KDIM%KLON,KDIM%KLEV),ZDELSO2_WET_O3(KDIM%KLON,KDIM%KLEV)
REAL(KIND=JPRB) :: ZDELSO2_WET_H2O2(KDIM%KLON,KDIM%KLEV)

!FIXME: proper 3d fields and arbitrary diagnostic "sets" shouldn't share the same dimension,
!  and this gets messy when KLEV changes!
REAL(KIND=JPRB) :: ZGLOMAP_DIAGS(KDIM%KLON,MAX(KDIM%KLEV,60),60)

INTEGER(KIND=JPIM) :: IO3_CHEM
REAL(KIND=JPRB) :: ZO3(KDIM%KLON,KDIM%KLEV)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "aer_phy3.intfb.h"




#include "aer_diagglomap.intfb.h"
#include "abor1.intfb.h"


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('AER_PHY3_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YGFL=>YDMODEL%YRML_GCONF%YGFL,YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, &
 & YDCOMPO=>YDMODEL%YRML_CHEM%YRCOMPO,YDCHEM=>YDMODEL%YRML_CHEM%YRCHEM)
ASSOCIATE(TSPHY=>YDPHY2%TSPHY, &
 & YSP_RR=>YDSURF%YSP_RR, YSP_SG=>YDSURF%YSP_SG, YSD_VF=>YDSURF%YSD_VF, &
 & YSP_SB=>YDSURF%YSP_SB, YSD_VD=>YDSURF%YSD_VD, NDIM=>YGFL%NDIM, &
 & LAERAOT=>YGFL%LAERAOT, LAERLISI=>YGFL%LAERLISI, &
 & NAERO_WVL_DIAG=>YGFL%NAERO_WVL_DIAG, &
 & NAERO_WVL_DIAG_TYPES=>YGFL%NAERO_WVL_DIAG_TYPES, &
 & NACTAERO=>YGFL%NACTAERO, &
 & LAERCHEM=>YGFL%LAERCHEM, &
 & AERO_SCHEME=>YDCOMPO%AERO_SCHEME, &
 & YO3=>YGFL%YO3, &
 & NCHEM=>YGFL%NCHEM, YCHEM=>YGFL%YCHEM, LCHEM_O3RAD=>YDCHEM%LCHEM_O3RAD)
!     ------------------------------------------------------------------

!*         0.     Preliminary computation

DO JK=2,KDIM%KLEV
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    ZTH(JL,JK)=(STATE%T(JL,JK-1)*PAUX%PAPRSF(JL,JK-1)&
      & *(PAUX%PAPRSF(JL,JK)-PAUX%PAPRS(JL,JK))&
      & +STATE%T(JL,JK)*PAUX%PAPRSF(JL,JK)*(PAUX%PAPRS(JL,JK)-PAUX%PAPRSF(JL,JK-1)))&
      & *(1.0_JPRB/(PAUX%PAPRS(JL,JK)*(PAUX%PAPRSF(JL,JK)-PAUX%PAPRSF(JL,JK-1))))
  ENDDO
ENDDO
DO JL=KDIM%KIDIA,KDIM%KFDIA
  ZTH(JL,1)=STATE%T(JL,1)-PAUX%PAPRSF(JL,1)*(STATE%T(JL,1)-ZTH(JL,2))&
    & /(PAUX%PAPRSF(JL,1)-PAUX%PAPRS(JL,2))
  ZTH(JL,KDIM%KLEV+1)=PSURF%PSP_RR(JL,YSP_RR%YT%MP9)
ENDDO

ZSNM = SUM(PSURF%PSP_SG(:,:,YSP_SG%YF%MP9),DIM=2) ! pre-compute total snow mass

! Select ozone to use in aerosol diagnostics
IF (NCHEM > 0 .AND. LCHEM_O3RAD) THEN
  IO3_CHEM=-999
  DO JGFL=1,YGFL%NCHEM
    IF (TRIM(YGFL%YCHEM(JGFL)%CNAME) == 'O3' )  IO3_CHEM = JGFL
  ENDDO
  IF (IO3_CHEM <= 0) THEN
    CALL ABOR1('AER_PHY3_LAYER: LCHEM_O3RAD set, but no chemistry O3 found')
  ENDIF
  ZO3(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) = PGFL(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,YGFL%YCHEM(IO3_CHEM)%MP)
ELSEIF (.NOT. YO3%LGP) THEN
  CALL ABOR1('AER_PHY3_LAYER: LCHEM_O3RAD unset, but no standalone gridpoint O3 found')
ELSE
  ZO3(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) = STATE%O3(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)
ENDIF

SELECT CASE (TRIM(AERO_SCHEME))

  CASE ("glomap")
    CALL ABOR1("OIFS - glomap should never be called from OIFS, EXIT")


  CASE ("aer")

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL AER_PHY3

    CALL AER_PHY3 &
       &( YDMODEL, &
       &  KDIM%KIDIA  , KDIM%KFDIA  , KDIM%KLON   , KDIM%KTDIA  , KDIM%KLEV , KDIM%KTILES , &
       &  KDIM%KFLDX  , KDIM%KLEVX, &
       &  GEMSL%ITRAC , GEMSL%IAERO ,  GEMSL%ICHEM, &
       &  PAUX%PRS1   , PAUX%PRSF1  , GEMSL%ZAEROP, GEMSL%ZCAERO, GEMSL%ZCEN , PAUX%PAPHIF, &
       &  FLUX%PFPLCL , FLUX%PFPLCN , FLUX%PFPLSL , FLUX%PFPLSN , PAUX%PGELAT, PAUX%PGELAM, &
       &  STATE%A , STATE%CLD(:,:,NCLDQI), STATE%CLD(:,:,NCLDQL), STATE%CLD(:,:,NCLDQR), &
       &  STATE%CLD(:,:,NCLDQS), PDIAG%PCOVPTOT, PDIAG%ZLU  ,&
       &  FLUX%PFCCQL,FLUX%PFCCQN,FLUX%PFCSQL,FLUX%PFCSQN,&
       &  ZO3          , STATE%Q      , STATE%T      , ZTH          , GEMSL%ZTENC     , GEMSL%ZCFLX  , &
       &  GEMSL%ZAERDDP, GEMSL%ZAERSDM, GEMSL%ZAERSRC, GEMSL%ZAERWS , GEMSL%ZAERGUST  , GEMSL%ZAERUST, GEMSL%ZAERMAP, &
       &  GEMSL%ZCLAERS, GEMSL%ZPRAERS, PCHEM2AER, &
       &  SURFL%ZALBD  , GEMSL%ZTAUAER, SURFL%ZFRTI  , PSURF%PSD_VF(:,YSD_VF%YLSM%MP)  , &
       &  ZSNM , AUXL%ZWND  , PSURF%PSP_SB(:,1,YSP_SB%YQ%MP9) , &
       &  GEMSL%ZAERFLX, GEMSL%ZAERLIF, &
       &  FLUX%PAERODDF, PSURF%PSD_VF(:,YSD_VF%YFCA1%MP),PSURF%PSD_VF(:,YSD_VF%YFCA2%MP),&
       &  TSPHY  , PGFL   , &
       &  PSURF%PSD_VD(:,YSD_VD%YODSS%MP), PSURF%PSD_VD(:,YSD_VD%YODDU%MP),&
       &  PSURF%PSD_VD(:,YSD_VD%YODOM%MP), PSURF%PSD_VD(:,YSD_VD%YODBC%MP), PSURF%PSD_VD(:,YSD_VD%YODSU%MP),&
       &  PSURF%PSD_VD(:,YSD_VD%YODNI%MP),PSURF%PSD_VD(:,YSD_VD%YODAM%MP), PSURF%PSD_VD(:,YSD_VD%YODSOA%MP),&
       &  PSURF%PSD_VD(:,YSD_VD%YODVFA%MP),PSURF%PSD_VD(:,YSD_VD%YODVSU%MP),PSURF%PSD_VD(:,YSD_VD%YODTOACC%MP),&
       &  PSURF%PSD_VD(:,YSD_VD%YAEPM1%MP),PSURF%PSD_VD(:,YSD_VD%YAEPM25%MP),PSURF%PSD_VD(:,YSD_VD%YAEPM10%MP),&
       &  GEMSL%ZAERAOT, GEMSL%ZAERLISI, PSURF%PSD_XA, &
       &  ZAERO_WVL_DIAG &
       &)


  CASE DEFAULT

    CALL ABOR1(" NO AEROSOL SCHEME "//TRIM(AERO_SCHEME) )

END SELECT

DO JAER=1,6
  DO JK=1,KDIM%KLEV
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      PRAD%PTAUAER(JL,JK,JAER)=GEMSL%ZTAUAER(JL,JK,JAER)
    ENDDO
  ENDDO
ENDDO

IF (LAERAOT) THEN
  DO JVAR=1,NPAERAOT
    DO JK=1,KDIM%KLEV
      DO JL=KDIM%KIDIA,KDIM%KFDIA
        PGFL(JL,JK,YGFL%YAERAOT(JVAR)%MP) = GEMSL%ZAERAOT(JL,JK,JVAR)
      ENDDO
    ENDDO
  ENDDO
ENDIF

IF (LAERLISI) THEN
  JGFL=0
  DO JVAR=1,NPAERLISI_VAR
    DO JWVL=1,NPAERLISI_WVL
      JGFL=JGFL+1
      DO JK=1,KDIM%KLEV
        DO JL=KDIM%KIDIA,KDIM%KFDIA
          PGFL(JL,JK,YGFL%YAERLISI(JGFL)%MP) = GEMSL%ZAERLISI(JL,JK,JWVL,JVAR)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF

DO JVAR=1,NAERO_WVL_DIAG_TYPES
  DO JWVL=1,NAERO_WVL_DIAG
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      PSURF%PSD_VD(JL,YSD_VD%YAERO_WVL_DIAG(JWVL,JVAR)%MP) = ZAERO_WVL_DIAG(JL,JWVL,JVAR)
    ENDDO
  ENDDO
ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_PHY3_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE AER_PHY3_LAYER
