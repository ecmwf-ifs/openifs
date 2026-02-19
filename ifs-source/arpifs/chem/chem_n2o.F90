! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

 SUBROUTINE CHEM_N2O &
 &    (YDML_GCONF,YDEPHY,YDML_CHEM,YDPHY2,KSTEP, KIDIA  , KFDIA , KLON, KLEV , KVCLIS,  &
 &     PTSTEP ,PDELP, PRS1, PRSF1, PGEOH, PTP, PKOZO,PPTROPO, &
 &     PCSZA, PGELAT, PGELAM, &
 &     PGEMU,   PCEN , PTENC1, POUT) 

!**   DESCRIPTION 
!     ----------
!
!   routine for C-IFS-N2O stratospheric chemistry ,
!           derived from BASCOE chemistry
!
!
!
!**   INTERFACE.
!     ----------
!          *CHEM_N2O* IS CALLED FROM *CHEM_MAIN*.

! INPUTS:
! -------
! KSTEP : Time step 
! KIDIA :  Start of Array  
! KFDIA :  End  of Array 
! KLON  :  Length of Arrays 
! KLEV  :  Number of Levels
! KVCLIS                      : Number Cariolle chemistry coefficinets 
! PKOZO(KLON,KLEV,KVCLIS)     : PHOTOCHEMICAL COEFFICIENTS COMPUTED FROM A 2D PHOTOCHEMICAL MODEL (KVCLIS=8)!
! PPTROPO  (KLON)              : Tropopause pressure           (Pa)
! PTSTEP:  Time step in seconds 
! PDELP(KLON,KLEV)            : PRESSURE DELTA in PRESSURE UNITES      (Pa)
! PRS1(KLON,0:KLEV)           : HALF-LEVEL PRESSURE           (Pa)
! PRSF1(KLON,KLEV)            : FULL-LEVEL PRESSURE           (Pa)
! PTP     (KLON,KLEV)         :  TEMPERATURE                  (K)
! PCSZA(KLON)                 : COS of Solar Zenit Angle
! PGELAM(KLON)                : LONGITUDE (RADIANS)
! PGELAT(KLON)                : LATITUDE (RADIANS) 
! PGEMU(KLON)                 : SINE OF LATITUDE
! PCEN(KLON,KLEV,NCHEM)       : CONCENTRATION OF TRACERS           (kg/kg)
!
! OUTPUTS:
! -------
! PTENC1  (KLON,KLEV,NCHEM)     : TENDENCY OF N2O BECAUSE OF CHEMISTRY (kg/kg s-1), no update
! POUT (KLON,KLEV,5)            : additional output
!
! LOCAL:
! -------
!
! ZCVM0(KLON,NCHEM)       : initial volume ratios OF TRACERS           (molec/cm3)
! ZCVM (KLON,NCHEM)       : final   volume ratios OF TRACERS           (molec/cm3)
!
!
!     AUTHOR.
!     -------
!        VINCENT HUIJNEN    *KNMI*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2015-05-01



USE MODEL_CHEM_MOD , ONLY : MODEL_CHEM_TYPE
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOEPHY   , ONLY : TEPHY
USE YOMPHY2  , ONLY : TPHY2
USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
! NCHEM : number of chemical species
! YCHEM : Data structure with Chemistry meta data 
! USE YOMLUN   , ONLY : NULOUT ,NULERR
USE YOMCST   , ONLY : RMD, RG , RPI , RNAVO
USE BASCOE_J_MODULE, ONLY :  J_N2O, J_O3_O1D


IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

TYPE(TEPHY)       ,INTENT(INOUT):: YDEPHY
TYPE(MODEL_CHEM_TYPE),INTENT(INOUT):: YDML_CHEM
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(TPHY2)       ,INTENT(INOUT):: YDPHY2
INTEGER(KIND=JPIM),INTENT(IN) :: KSTEP, KIDIA , KFDIA , KLON , KLEV, KVCLIS
REAL(KIND=JPRB),INTENT(IN)    :: PTSTEP
REAL(KIND=JPRB),INTENT(IN)    :: PDELP(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PRSF1(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PRS1(KLON,0:KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PGEOH(KLON,0:KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PTP(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PPTROPO(KLON)                   
REAL(KIND=JPRB),INTENT(IN)    :: PKOZO(KLON,KLEV,KVCLIS)
REAL(KIND=JPRB),INTENT(OUT)   :: PTENC1(KLON,KLEV,YDML_GCONF%YGFL%NCHEM)
REAL(KIND=JPRB),INTENT(IN)    :: PCEN(KLON,KLEV,YDML_GCONF%YGFL%NCHEM) 
REAL(KIND=JPRB),INTENT(IN)    :: PCSZA(KLON)                  
REAL(KIND=JPRB),INTENT(IN)    :: PGELAT(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PGELAM(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PGEMU(KLON)
REAL(KIND=JPRB),INTENT(OUT)   :: POUT(KLON,KLEV,5)  


!*       0.5   LOCAL VARIABLES
!              ---------------

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

! * Lat /Lon time
! REAL(KIND=JPRB) , DIMENSION(KLON)   :: ZLAT
! REAL(KIND=JPRB) , DIMENSION(KLON)   :: ZLON

! * counters
INTEGER(KIND=JPIM) :: JK, JL, JT, JLEV

! * chemical data 
REAL(KIND=JPRB) , DIMENSION(KLON,YDML_GCONF%YGFL%NCHEM)     :: ZCVM
REAL(KIND=JPRB) , DIMENSION(KLON,YDML_GCONF%YGFL%NCHEM)     :: ZCVM0 
REAL(KIND=JPRB) , DIMENSION(KLON)                :: ZDENS 
REAL(KIND=JPRB)                                  :: ZAIRDM


! * TRACER data
INTEGER(KIND=JPIM) :: IO3,IN2O

!* Selection of photorates
INTEGER(KIND=JPIM) :: J_DISS(2)

! * Photolysis data: cloud/ozone info, actinic fluxes, photolysis rates
REAL(KIND=JPRB)                          :: ZCOLO3_DU(KLON,KLEV)


REAL(KIND=JPRB)                           :: ZCST  ! molec/cm2 -> DU
REAL(KIND=JPRB)                           :: ZVMRO3,ZSZA, ZHGT, ZO1D_MMR
REAL(KIND=JPRB)                           :: ZAJVAL(KLON,2),ZO1D_VM(KLON)
REAL(KIND=JPRB),PARAMETER                 :: ZSMALL_VMR=1.0E-30                        


! Surface boundary conditions (BASCOE)
REAL(KIND=JPRB) , DIMENSION(KLON)         :: ZTENBC
REAL(KIND=JPRB)                           :: ZBCVAL

! Variables used to compute altitude
REAL(KIND=JPRB)       :: ZPSURF_STD,ZSURF_H, ZTHKNESS,ZHGT_BASCOE(KLON,KLEV)

! * O3 tendency from Cariolle chemistry 
REAL(KIND=JPRB) , DIMENSION(KLON,KLEV)   :: ZTENO3COR

! ------------------------------------------------------------------
#include "fcttim.func.h"
!-------------------------------------------------------------------
! #include "abor1.intfb.h"
#include "n2o_j_interp.intfb.h"
#include "n2o_o1d_interp.intfb.h"
#include "n2o_solver.intfb.h"
#include "o3chem.intfb.h"

!-----------------------------------------------------------------------
! chemistry scheme name - this will later also come from external input
IF (LHOOK) CALL DR_HOOK('CHEM_N2O',0,ZHOOK_HANDLE )
!-----------------------------------------------------------------------
! chemistry scheme name - this will later also come from external input
ASSOCIATE(LCHEM_JOUT=>YDML_CHEM%YRCHEM%LCHEM_JOUT, NCHEM=>YDML_GCONF%YGFL%NCHEM, YCHEM=>YDML_GCONF%YGFL%YCHEM )
!-----------------------------------------------------------------------


! Lat / Lon
!DO JL=KIDIA,KFDIA
!  ZLAT(JL)=(180.0_JPRB/RPI)*PGELAT(JL)
!  ZLON(JL)=(180.0_JPRB/RPI)*PGELAM(JL)
!ENDDO

! Initialize tendencies

PTENC1(:,:,:) = 0.0

! Initialize extra output
POUT(:,:,:) = 0.0

! Initialize requested photorates
J_DISS=(/J_N2O,J_O3_O1D/)

! Find tracer indices
  IO3=-1
  DO JT=1,NCHEM
    IF (TRIM (YCHEM(JT)%CNAME) == 'O3' ) IO3=JT 
  ENDDO 
  IN2O=-1
  DO JT=1,NCHEM
    IF (TRIM (YCHEM(JT)%CNAME) == 'N2O' ) IN2O=JT 
  ENDDO 


! ZRGI=1.0_JPRB/RG

! BASCOE variant to compute overhead column [DU]
!-----------------------------------------------------------------------
! Calculate overhead ozone columns *at* levels:  Zcst * SUM( vmr * delta_p )
! where delta_p is between levels and vmr are at mid-levels
!-----------------------------------------------------------------------
ZCST  = (1.0E-4/RG)*( RNAVO / (1.0E-3* RMD) ) * 1.0E3 / 2.687E19  ! molec/cm2 -> DU
DO JL=KIDIA,KFDIA
  ! convert to mixing ratio
  ZVMRO3=MAX(PCEN(JL,1,IO3) / YCHEM(IO3)%RMOLMASS *RMD ,0._JPRB) 
  ! convert to DU
  ZCOLO3_DU(JL,1)  = ZCST * ZVMRO3 * ( PRSF1(JL,1) - 0._JPRB )
ENDDO

DO JLEV=2,KLEV
  DO JL=KIDIA,KFDIA
    ! convert to mixing ratio
    ZVMRO3=MAX(0.5_JPRB*(PCEN(JL,JLEV,IO3)+PCEN(JL,JLEV-1,IO3)) / YCHEM(IO3)%RMOLMASS *RMD ,0._JPRB) 
    ! convert to DU
    ZCOLO3_DU(JL,JLEV)=ZCOLO3_DU(JL,JLEV-1)+ZCST * ZVMRO3* ( PRSF1(JL,JLEV) - PRSF1(JL,JLEV-1))
  ENDDO
ENDDO


! BASCOE way to compute model altitude, first surface model level:
ZPSURF_STD=101325._JPRB ! std p at surf (Pa)
DO JL=KIDIA,KFDIA

  IF( PRS1(JL,KLEV) < ZPSURF_STD ) THEN
    ZSURF_H = 7._JPRB*LOG( ZPSURF_STD / PRS1(JL,KLEV) )
  ELSE
    ZSURF_H=0.0_JPRB
  ENDIF
  ZTHKNESS = PTP(JL,KLEV)*287./9.806*LOG(PRS1(JL,KLEV)/PRSF1(JL,KLEV))
  ZHGT_BASCOE(JL,KLEV) = ZSURF_H + 1.E-3*ZTHKNESS
ENDDO

DO JLEV=KLEV-1,1,-1
  DO JL=KIDIA,KFDIA
     ZTHKNESS=0.5*(PTP(JL,JLEV+1)+PTP(JL,JLEV))*287./9.806*&
     &                 LOG(PRSF1(JL,JLEV+1)/PRSF1(JL,JLEV))
     ZHGT_BASCOE(JL,JLEV)=ZHGT_BASCOE(JL,JLEV+1)+1E-3*ZTHKNESS
  ENDDO
ENDDO



!2.0 loop over levels for solving chemistry
DO JK=1,KLEV
 
! switch on chemistry from 400hPa onwards:
IF ( PRSF1(KIDIA,JK) < 40000_JPRB ) THEN

!2.1 Loop over lon/lat
  DO JL=KIDIA,KFDIA
     
!   ! switch on chemistry from 400hPa onwards:
!   IF ( PRSF1(JL,JK) < 40000_JPRB ) THEN

    ZSZA = ACOS(PCSZA(JL))*180_JPRB/RPI

    IF( ZSZA >= 96._JPRB ) THEN
      ZAJVAL(JL,:) = 0._JPRB
      ZO1D_MMR = 0._JPRB
    ELSE
      !VH ZHGT = PGEOH(JL,JK-1)  * ZRGI *1e-3 ! height in km
      ! ZHGT = 0.5*(PGEOH(JL,JK-1) + PGEOH(JL,JK))  * ZRGI *1e-3 ! height in km
      ! 
      ZHGT = ZHGT_BASCOE(JL,JK)

      ! Set maximimum height to ~110 km altitude 
      ZHGT = MIN(ZHGT,109.9_JPRB)

      ! ZCOLO3_DU: o3 overhead column in DU, converted from kg/m2
      ! compute photolysis rate ZAJVAL for selected reaction 
      CALL N2O_J_INTERP( ZSZA, ZHGT  , ZCOLO3_DU(JL,JK), J_DISS, ZAJVAL(JL,:) )
      !*  Estimate corresponding O1D field from lookup-table

      CALL N2O_O1D_INTERP(YDML_GCONF%YGFL,ZAJVAL(JL,:),PRSF1(JL,JK),ZO1D_MMR)
    ENDIF



       
!*  Air density (molec/cm3) 
    ZDENS(JL) = 7.24291E16_JPRB*PRSF1(JL,JK)/PTP(JL,JK) 
!*  Air density mutiplied with RMD (dry air molar mass) for efficiency  
    ZAIRDM = ZDENS(JL) * RMD


!*  convert tracer concentrations from kg/kg to molec/cm3
    DO JT=1,NCHEM
!*     assure positivity for initial concentrations
      ZCVM0(JL,JT) = MAX(PCEN(JL,JK,JT) / YCHEM(JT)%RMOLMASS *ZAIRDM ,ZSMALL_VMR) 
      ZCVM(JL,JT)  = ZCVM0(JL,JT)   
    ENDDO
    ZO1D_VM(JL) = ZO1D_MMR / 16. * ZAIRDM

    IF (LCHEM_JOUT) THEN
      POUT(JL,JK,3) =  ZO1D_MMR
    ENDIF

!VH Outside this loop...  ENDIF ! Pres < 400 hPa
  ENDDO ! Loop over JL


!* Fixed concentrations - check ZDENS / loop!
!   ZO2 = 0.209*ZDENS    ! O2 number density
!   ZN2 = 0.781*ZDENS    ! N2 number density


! ----------------------------------------------------------------------
!  Compute heterogeneous reaction rates (ZRHET)
! ----------------------------------------------------------------------

    ! Output something which one finds interesting
    IF (LCHEM_JOUT) THEN
      DO JL=KIDIA,KFDIA
        POUT(JL,JK,2) =  ZAJVAL(JL,2)
      ENDDO
    ENDIF
   

! ----------------------------------------------------------------------
! Call simple N2O solver... 
! ----------------------------------------------------------------------

    CALL N2O_SOLVER(YDML_GCONF%YGFL,KIDIA,KFDIA,KLON,IN2O,PTSTEP,ZAJVAL,ZCVM0,ZO1D_VM,ZCVM,PRSF1(:,JK),PTP)



!5.0 convert N2O concentration tendencies to mass mixing ratio
   JT=IN2O 
   DO JL=KIDIA,KFDIA
!*  Air density mutiplied with RMD (dry air molar mass) for efficiency  
     ZAIRDM = ZDENS(JL) * RMD
     PTENC1(JL,JK,JT) =   (ZCVM(JL,JT)-ZCVM0(JL,JT)) * YCHEM(JT)%RMOLMASS /(ZAIRDM*PTSTEP)
   ENDDO

ENDIF ! IF P<400 hPa

ENDDO ! loop over levels

! ----------------------------------------------------------------------
! Call O3 solver... 
! ----------------------------------------------------------------------

! for O3
 CALL O3CHEM(YDML_GCONF%YRRIP,YDEPHY,YDML_CHEM%YROZO,YDPHY2,KIDIA,KFDIA,KLON,1,KLEV,KVCLIS,PGEMU,PCSZA,PRS1,PRSF1,PKOZO, &
   & PDELP,PTP,PCEN(:,:,IO3),ZTENO3COR)

! 6.1 replace tendecies for ozone above ZPMAXO3CAR with Cariolle tendency.
! The stratosphere is defined at zonal mean ozone levels exceeding 150 ppb
! based on a climatology. This boundary can approximately be described
! by the function P = 230-148*( cos(lat) )^4 (hPa)
! stratospheric Cariolle chemistry at (P-20) hPa is applied.

 DO JK=1,KLEV
   DO JL=KIDIA,KFDIA
!     ZPMAXO3CAR=PPTROPO(JL) 
!     IF (PRSF1(JL,JK)<= ZPMAXO3CAR)  THEN
!        ! completely stratospheric region: P < ZPMIN
        PTENC1(JL,JK,IO3) =  ZTENO3COR(JL,JK) 
!     ENDIF
   ENDDO
 ENDDO    

! 
! Very simple boundary condition for N2O
! convert to mass ratio [kg/kg]
ZBCVAL = 3.22E-7* YCHEM(IN2O)%RMOLMASS / RMD 
DO JL=KIDIA,KFDIA
  ZTENBC(JL) = ZBCVAL - PCEN(JL,KLEV,IN2O)
  PTENC1(JL,KLEV,IN2O) =  ZTENBC(JL)  / PTSTEP
  ! store this special LBC budget contribution in POUT(:,x,x) 
  POUT(JL,1,1)=PTSTEP *PTENC1(JL,KLEV,IN2O)*PDELP(JL,KLEV) / RG 
ENDDO




END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CHEM_N2O',1,ZHOOK_HANDLE )
END SUBROUTINE CHEM_N2O





