! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AERINI_LAYER(YDGEOMETRY,YDSURF,&
  & YDMODEL,KDIM,PAUX,STATE,R,S,PSURF,SURFL,GEMSL)

!**** *AERINI_LAYER* - Layer routine calling time stepping of first part
!                      of prognostic aerosol computations


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

!     ==== Input/output ====
! PSURF    : Derived variables for general surface quantities
! SURFL    : Derived variables for local surface quantities
! GEMSL    : Derived variable for local GEMS quantities.



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
!      Original : 2012-12-03  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!     JJMorcrette 20131001 DMS-related variables and diagnostics
!     SRemy       20150126 Add arguments for biomass burning injection heights
!     SRemy       20160426 SOA from CO
!     SRemy       20170421 Altitude of volcanoes for SO2 emissions
!     E. Dutra/G.Arduini, Jan 2018: snow multi-layer: pre-compute total snow mass
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE VARIABLE_MODULE, ONLY : VARIABLE_3D
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMPHYDER, ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE,&
   & SURF_AND_MORE_TYPE, SURF_AND_MORE_LOCAL_TYPE, GEMS_LOCAL_TYPE
USE YOMCT3   , ONLY : NSTEP
USE YOMCST    ,ONLY : RA, RPI





!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)                 ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)                    ,INTENT(INOUT) :: YDSURF
TYPE(MODEL)                    ,INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)          ,INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                ,INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              ,INTENT (IN)   :: STATE
TYPE(VARIABLE_3D)              ,INTENT (IN)   :: R, S
TYPE (SURF_AND_MORE_TYPE)      ,INTENT(INOUT) :: PSURF
TYPE (SURF_AND_MORE_LOCAL_TYPE),INTENT(INOUT) :: SURFL
TYPE (GEMS_LOCAL_TYPE)         ,INTENT(INOUT) :: GEMSL
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JL, ISNH3_C,ISHNO3_C,JK,JT
! names begin with "IX" to avoid confusion with the species indices used in the chemical schemes
INTEGER(KIND=JPIM), PARAMETER :: IXNH3 = 1_JPIM
INTEGER(KIND=JPIM), PARAMETER :: IXHNO3 = 2_JPIM
! "default" species index values (for species not used in comparisons)
INTEGER(KIND=JPIM), PARAMETER :: IXNOTUSED = -1_JPIM
INTEGER(KIND=JPIM), DIMENSION(YDMODEL%YRML_GCONF%YGFL%NCHEM) :: IXCHEM
! IFS-GLOMAP SPECIFIC
REAL(KIND=JPRB),ALLOCATABLE :: ZGLO_MASS_EMS(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZGLO_NUM_EMS(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZAREA(:)
REAL(KIND=JPRB)             :: ZNO3(KDIM%KLON,KDIM%KLEV)
REAL(KIND=JPRB)             :: ZHNO3(KDIM%KLON,KDIM%KLEV)
REAL(KIND=JPRB)             :: ZNH3(KDIM%KLON,KDIM%KLEV)
REAL(KIND=JPRB)             :: ZTHNO3(KDIM%KLON,KDIM%KLEV)
REAL(KIND=JPRB)             :: ZTNH3(KDIM%KLON,KDIM%KLEV)
REAL(KIND=JPRB)             :: ZNH4(KDIM%KLON,KDIM%KLEV)
REAL(KIND=JPRB)             :: ZSO4(KDIM%KLON,KDIM%KLEV)



REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZSNM(KDIM%KLON)
!-----------------------------------------------------------------------

#include "aer_wind.intfb.h"
#include "aer_phy2.intfb.h"
#include "abor1.intfb.h"




!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('AERINI_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
  & YDMP=>YDGEOMETRY%YRMP, YGFL=>YDMODEL%YRML_GCONF%YGFL,YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, &
  & YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD, YDCOMPO=>YDMODEL%YRML_CHEM%YRCOMPO, &
  & YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM, &
  & YDRIP=>YDMODEL%YRML_GCONF%YRRIP)
ASSOCIATE(NACTAERO=>YGFL%NACTAERO, &
 & LECSRAD=>YDERAD%LECSRAD, NCSRADF=>YDERAD%NCSRADF, NSW=>YDERAD%NSW, &
 & NSTART=>YDRIP%NSTART, &
 & LAERNITRATE => YDCOMPO%LAERNITRATE, &
 & YCHEM=>YGFL%YCHEM, &
 & AERO_SCHEME=>YDCOMPO%AERO_SCHEME, &
 & YSD_VD=>YDSURF%YSD_VD, YSD_VF=>YDSURF%YSD_VF, YSP_RR=>YDSURF%YSP_RR, &
 & YSP_SB=>YDSURF%YSP_SB, YSP_SG=>YDSURF%YSP_SG, TSPHY=>YDPHY2%TSPHY)
!     ------------------------------------------------------------------

IXCHEM(:) = IXNOTUSED
DO JT=1,YGFL%NCHEM
  SELECT CASE (YCHEM(JT)%CNAME)
    CASE ("NH3")
      IXCHEM(JT) = IXNH3
    CASE ("HNO3")
      IXCHEM(JT) = IXHNO3
  END SELECT
ENDDO

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL AER_WIND AND AER_PHY2


CALL AER_WIND&
 &( KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KTILES,&
 &STATE%U(:,KDIM%KLEV),STATE%V(:,KDIM%KLEV),STATE%T(:,KDIM%KLEV),STATE%Q(:,KDIM%KLEV),&
 &PAUX%PAPRS(:,KDIM%KLEV),PAUX%PGEOM1(:,KDIM%KLEV),&
 &PSURF%PUSTRTI,PSURF%PVSTRTI,PSURF%PAHFSTI,PSURF%PEVAPTI,SURFL%ZFRTI,GEMSL%ZAZ0M,&
 &  GEMSL%ZAERWS, GEMSL%ZAERGUST, GEMSL%ZAERUST )

ZSNM = SUM(PSURF%PSP_SG(:,:,YSP_SG%YF%MP9),DIM=2) ! pre-compute total snow mass

SELECT CASE (TRIM(AERO_SCHEME))

  CASE ("glomap")
    CALL ABOR1("OIFS - glomap should never be called from OIFS, EXIT")


  CASE ("aer")

    CALL AER_PHY2 &
 &( YDGEOMETRY, YDMODEL, KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON , KDIM%KTDIA, KDIM%KLEV, KDIM%KFLDX , KDIM%KLEVX,&
 &  KDIM%KTILES, KDIM%KSTGLO, GEMSL%ITRAC, GEMSL%IAERO, NSW,&
 &  PAUX%PRS1 , PAUX%PRSF1, PAUX%PAPHI, STATE%T , PAUX%PVERVEL, GEMSL%ZCEN  , PAUX%PGEOMH,&
 &  PSURF%PSD_VD(:,YSD_VD%YALB%MP), SURFL%ZALBD, PSURF%PSD_VF(:,YSD_VF%YALUVD%MP),&
 &  PSURF%PSD_VF(:,YSD_VF%YAERDEP%MP),PSURF%PSD_VF(:,YSD_VF%YAERLTS%MP),PSURF%PSD_VF(:,YSD_VF%YAERSCC%MP),&
 &  GEMSL%ZAERWS, GEMSL%ZAERGUST, GEMSL%ZAERUST,&
 &  PSURF%PSD_VF(:,YSD_VF%YURBF%MP), R%P(:,KDIM%KLEV), S%P(:,KDIM%KLEV), &
 &  PSURF%PSD_VF(:,YSD_VF%YSO2DD%MP),&
 &  PSURF%PSD_VF(:,YSD_VF%YDSF%MP), &
 &  PSURF%PSD_VF(:,YSD_VF%YDSZ%MP), &
 &  PSURF%PSD_VF(:,YSD_VF%YDMSO%MP),&
 &  PSURF%PSD_VF(:,YSD_VF%YCI%MP),&
 &  SURFL%ZFRTI, PSURF%PSD_VF(:,YSD_VF%YLSM%MP), PSURF%PSD_VF(:,YSD_VF%YCLK%MP), STATE%Q  ,PSURF%PSD_VF(:,YSD_VF%YSST%MP), &
 &  ZSNM , PSURF%PSP_RR(:,YSP_RR%YT%MP9)  , PAUX%PGELAM, PAUX%PGELAT, PAUX%PGEMU, SURFL%ZHSDFOR,&
 &  STATE%U(:,KDIM%KLEV) , STATE%V(:,KDIM%KLEV) , PSURF%PSP_SB(:,1,YSP_SB%YQ%MP9)   , TSPHY, GEMSL%ZAZ0M,&
 &  GEMSL%ICHEM,&
 &  PSURF%ITVL, PSURF%ITVH, PSURF%PCVL, PSURF%PCVH, PSURF%PLAIL, PSURF%PLAIH ,&
 &  GEMSL%ZCFLX, GEMSL%ZTENC,GEMSL%ZDDVLC,&
 &  GEMSL%ZDMSO, GEMSL%ZLDAY, GEMSL%ZLISS, GEMSL%ZSO2, GEMSL%ZTDMS,&
 &  GEMSL%ZAERDDP, GEMSL%ZAERSDM, GEMSL%ZAERSRC, GEMSL%ZAERMAP, GEMSL%ZAERFLX, GEMSL%ZAERLIF,&
 &  GEMSL%ZDMSI, GEMSL%ZODMS, PSURF%PSD_XA &
 &)

    IF (LECSRAD .AND. NCSRADF == 1 .AND. NACTAERO >= 12) THEN
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    PSURF%PSD_XA(JL,15,14)=SURFL%ZHSDFOR(JL)                        ! stand.dev. form orography
    PSURF%PSD_XA(JL,16,14)=SURFL%ZFRTI(JL,4)                        ! tile fraction
    PSURF%PSD_XA(JL,17,14)=SURFL%ZFRTI(JL,8)                        ! tile fraction
    PSURF%PSD_XA(JL,18,14)=0._JPRB
    PSURF%PSD_XA(JL,19,14)=PSURF%PSD_VD(JL,YSD_VD%YALB%MP)          ! background albedo
    PSURF%PSD_XA(JL,20,14)=PSURF%PSD_VF(JL,YSD_VF%YALUVD%MP)        ! diffuse UVis albedo
    PSURF%PSD_XA(JL,21,14)=SURFL%ZALBD(JL,1)                        ! diffuse albedo 1st interval
    PSURF%PSD_XA(JL,22,14)=SUM(PSURF%PSP_SG(JL,1:KDIM%KLEVSN,YSP_SG%YF%MP9),DIM=1)          ! snow SNS
! - DMS-related
    PSURF%PSD_XA(JL,23,14)=GEMSL%ZDMSO(JL)
    PSURF%PSD_XA(JL,24,14)=GEMSL%ZLDAY(JL)
    PSURF%PSD_XA(JL,25,14)=GEMSL%ZSO2(JL)
    PSURF%PSD_XA(JL,26,14)=GEMSL%ZAERWS(JL)
    PSURF%PSD_XA(JL,27,14)=GEMSL%ZLISS(JL)
    PSURF%PSD_XA(JL,28,14)=GEMSL%ZTDMS(JL)
    PSURF%PSD_XA(JL,29,14)=GEMSL%ZODMS(JL)
    PSURF%PSD_XA(JL,30,14)=GEMSL%ZDMSI(JL)
  ENDDO
   ENDIF

  CASE DEFAULT

    CALL ABOR1(" NO AEROSOL SCHEME "//TRIM(AERO_SCHEME) )

END SELECT


!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AERINI_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE AERINI_LAYER
