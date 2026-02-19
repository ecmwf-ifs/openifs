! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE UPD_GHGCLIM( YDERAD, YDERDI, KINDAT, KMINUT, KSTEP, KUPGHG )

! PURPOSE
! -------
! Update the greenhouse gas climatologies by interpolating from the
! monthly averages to the current model time
!
! Before cycle 47 the routine SU_GHGCLIM both loaded (via calling many
! subroutines) and updated the gas climatologies; this routine
! includes a cleaned version of just the update part.  Loading of the
! data from file is now performed in the YOMCLIM module.
!
! INTERFACE
! ---------
! UPD_GHGCLIM is called from UPDTIM
!
!
! AUTHORS
! -------
! R. Hogan, ECMWF 2019-01-25
!
! MODIFICATIONS
! -------------
!

USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOERAD    ,ONLY : TERAD
USE YOMCLIM   ,ONLY : YGHGCLIM ! Climatology from NetCDF file
USE YOERADGHG ,ONLY : YRADGHG  ! Climatology interpolated to current model time
USE YOERDI    ,ONLY : TERDI
! Molar masses
USE YOMCST    ,ONLY : RMD, RMNO2, RMO3

IMPLICIT NONE

TYPE(TERAD)       ,INTENT(INOUT) :: YDERAD
TYPE(TERDI)       ,INTENT(INOUT) :: YDERDI
INTEGER(KIND=JPIM),INTENT(IN)    :: KINDAT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KMINUT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP
INTEGER(KIND=JPIM),INTENT(INOUT) :: KUPGHG

! Integrated number of days through the year at the end of each month
REAL(KIND=JPRB), PARAMETER :: PP_DAY(*) = &
 & [ 31.00_JPRB,  59.25_JPRB,  90.25_JPRB, 120.25_JPRB, 151.25_JPRB, 181.25_JPRB, &
 &  212.25_JPRB, 243.25_JPRB, 273.25_JPRB, 304.25_JPRB, 334.25_JPRB, 365.25_JPRB]
! Number of minutes since the beginning of the year at the mid-point
! of each month
REAL(KIND=JPRB), PARAMETER :: PP_MIN(*) = &
 & [ 22320._JPRB,  64980._JPRB, 107640._JPRB, 151560._JPRB, 195480._JPRB, 239400._JPRB, &
 &  283320._JPRB, 327960._JPRB, 371880._JPRB, 415800._JPRB, 459720._JPRB, 503640._JPRB]

! Day and month
INTEGER(KIND=JPIM) :: IDY, IMN

! Number of minutes through the year
REAL(KIND=JPRB) :: ZMIN

! Indices of months between which we interpolate
INTEGER(KIND=JPIM) :: IM1, IM2, JM
! Month weight
REAL(KIND=JPRB)    :: ZTIMI

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

! -----------------------------------------------------------------------
#include "abor1.intfb.h"
#include "fcttim.func.h"

IF (LHOOK) CALL DR_HOOK('UPD_GHGCLIM',0,ZHOOK_HANDLE)

ASSOCIATE(NGHGRAD=>YDERAD%NGHGRAD, &
 & RCO2=>YDERDI%RCARDI, RCFC11=>YDERDI%RCFC11, RCFC12=>YDERDI%RCFC12, &
 & RCH4=>YDERDI%RCH4, RN2O=>YDERDI%RN2O, RHCFC22=>YDERDI%RCFC22, RCCL4=>YDERDI%RCCL4)
! -----------------------------------------------------------------------

!                 -----------------------------------
!          1.     TIME INDEX WITHIN THE CLIMATOLOGIES
!                 -----------------------------------
IDY=NDD(KINDAT)-1
IMN=NMM(KINDAT)
IF (IMN == 1) THEN
  ZMIN = REAL(IDY,KIND(ZMIN))*1440._JPRB + KMINUT
ELSEIF (IMN == 2) THEN
  IF(IDY == 28) IDY=IDY-1
  ! A DAY IN FEB. IS 28.25*24*60/28=1452.8571min LONG.
  ZMIN = 44640._JPRB+REAL(IDY,KIND(ZMIN))*1452.8571_JPRB+KMINUT
ELSE
  ZMIN = (PP_DAY(IMN-1)+REAL(IDY,KIND(ZMIN)))*1440._JPRB+KMINUT
ENDIF
! 525960=MINUTES IN A SIDERAL YEAR (365.25d)
ZMIN = MOD(ZMIN, 525960._JPRB)

IM1=0
IM2=0
IF (ZMIN <= PP_MIN(1)) THEN
  IM1=12
  IM2=1
  ZTIMI=(PP_MIN(1)-ZMIN)/44640._JPRB
ELSEIF(ZMIN > PP_MIN(12)) THEN
  IM1=12
  IM2=1
  ZTIMI=(548280._JPRB-ZMIN)/44640._JPRB
! 548280.=(365.25d + 15.5d)*24*60
ELSE
  DO JM=1,11
    IF (ZMIN > PP_MIN(JM) .AND. ZMIN <= PP_MIN(JM+1)) THEN
      IM1=JM
      IM2=JM+1
      ZTIMI=(ZMIN-PP_MIN(IM2))/(PP_MIN(IM1)-PP_MIN(IM2))
    ENDIF
  ENDDO
  IF (IM1 == 0 .OR. IM2 == 0) THEN
    CALL ABOR1('Problem with time interpolation in UPD_GHGCLIM!')
  ENDIF
ENDIF

!                 -----------------------------------
!          2.     INITIALIZE CURRENT-TIME CONCENTRATIONS
!                 -----------------------------------

! This is now done in suecrad.F90

!                 -----------------------------------
!          3.     INTERPOLATE CONCENTRATIONS IN TIME
!                 -----------------------------------
!
!  Field is a/ transformed from vmr to mmr (kg/kg)
!           b/ potentially modified to account for the ongoing increase in GHG
!           c/ multiplied by the pressure difference across the layer
!
! NGHGRAD index for defining which trace gas is taken from the 2D (lat-height) climatologies
! = 1   CO2
! = 2   CH4
! = 3   N2O
! = 4   NO2
! = 5   CFC11
! = 6   CFC12
! = 7   CFC22
! = 8   CCl4
! = 11  CO2+CH4
! = 12  CO2+CH4+N2O
! = 15  CO2+CH4+N2O+CFC11
! = 16  CO2+CH4+N2O+CFC11+CFC12
! = 17  CO2+CH4+N2O+CFC11+CFC12+O3
! = 18  CO2+CH4+N2O+CFC11+CFC12+O3+CFC22
! = 19  CO2+CH4+N2O+CFC11+CFC12+O3+CFC22+CCl4
! = 20  CO2+CH4+N2O+CFC11+CFC12+CFC22+CCl4+NO2

! The fifth argument in the INTERP_GAS routine provides the scaling,
! the ratio of the target mass mixing ratio to the reference volume
! mixing ratio: the climatic variations in CO2, CH4, N2O, CFC11,
! CFC12, HCFC22 and CCl4 are accounted for in UPDRGAS, but not NO2 or
! O3, so for the latter we enter the ratio of the molar mass of the
! gas to the molar mass of dry air.

IF (NGHGRAD == 1 .OR. NGHGRAD >= 10) THEN
  CALL INTERP_GAS(IM1, IM2, ZTIMI, YRADGHG%PRESSURE_HL,    RCO2/YGHGCLIM%SURF_MEAN_CO2, &
       &  YGHGCLIM%VMR_CO2,   YRADGHG%MASS_CO2)
ENDIF

IF (NGHGRAD == 2 .OR. NGHGRAD >= 11) THEN
  CALL INTERP_GAS(IM1, IM2, ZTIMI, YRADGHG%PRESSURE_HL,    RCH4/YGHGCLIM%SURF_MEAN_CH4, &
       &  YGHGCLIM%VMR_CH4,   YRADGHG%MASS_CH4)
ENDIF

IF (NGHGRAD == 3 .OR. NGHGRAD >= 12) THEN
  CALL INTERP_GAS(IM1, IM2, ZTIMI, YRADGHG%PRESSURE_HL,    RN2O/YGHGCLIM%SURF_MEAN_N2O, &
       &  YGHGCLIM%VMR_N2O,   YRADGHG%MASS_N2O)
ENDIF

IF (NGHGRAD == 4 .OR. NGHGRAD >= 20) THEN
  CALL INTERP_GAS(IM1, IM2, ZTIMI, YRADGHG%PRESSURE_HL,    RMNO2/RMD, &
       &  YGHGCLIM%VMR_NO2,   YRADGHG%MASS_NO2)
ENDIF

IF (NGHGRAD == 5 .OR. NGHGRAD >= 15) THEN
  CALL INTERP_GAS(IM1, IM2, ZTIMI, YRADGHG%PRESSURE_HL,    RCFC11/YGHGCLIM%SURF_MEAN_CFC11, &
       &  YGHGCLIM%VMR_CFC11, YRADGHG%MASS_CFC11)
ENDIF

IF (NGHGRAD == 6 .OR. NGHGRAD >= 16) THEN
  CALL INTERP_GAS(IM1, IM2, ZTIMI, YRADGHG%PRESSURE_HL,    RCFC12/YGHGCLIM%SURF_MEAN_CFC12, &
       &  YGHGCLIM%VMR_CFC12, YRADGHG%MASS_CFC12)
ENDIF

IF (NGHGRAD == 7 .OR. NGHGRAD == 17 .OR. NGHGRAD >= 21) THEN
  CALL INTERP_GAS(IM1, IM2, ZTIMI, YRADGHG%PRESSURE_HL,    RMO3/RMD, &
       &  YGHGCLIM%VMR_O3,    YRADGHG%MASS_O3)
ENDIF

IF (NGHGRAD == 8 .OR. NGHGRAD >= 18) THEN
  CALL INTERP_GAS(IM1, IM2, ZTIMI, YRADGHG%PRESSURE_HL,    RHCFC22/YGHGCLIM%SURF_MEAN_HCFC22, &
       &  YGHGCLIM%VMR_HCFC22,YRADGHG%MASS_HCFC22)
ENDIF

IF (NGHGRAD == 9 .OR. NGHGRAD >= 19) THEN
  CALL INTERP_GAS(IM1, IM2, ZTIMI, YRADGHG%PRESSURE_HL,    RCCL4/YGHGCLIM%SURF_MEAN_CCL4, &
       &  YGHGCLIM%VMR_CCL4,  YRADGHG%MASS_CCL4)
ENDIF


END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('UPD_GHGCLIM',1,ZHOOK_HANDLE)

CONTAINS

! Interpolate gas layer masses in time from volume mixing ratio
! (unitless) to layer mass (kg m-2)
SUBROUTINE INTERP_GAS(KM1, KM2, PTIMI, PPRESSURE_HL, PSCALING, PGAS_VMR, PGAS_MASS)

  IMPLICIT NONE

  ! Indices of the months between which we are interpolating
  INTEGER(KIND=JPIM), INTENT(IN) :: KM1, KM2
  ! Time interpolation weight
  REAL(KIND=JPRB)   , INTENT(IN) :: PTIMI
  ! Scaling to allow for climatic change in concentrations, in units
  ! of the ratio of the target mass mixing ratio to the reference
  ! volume mixing ratio; if the target mmr is to be the same as the
  ! reference then enter the ratio of molar mass of the gas to the
  ! molar mass of dry air.
  REAL(KIND=JPRB)   , INTENT(IN) :: PSCALING
  ! Half-level pressure (Pa), *indexing from 0*
  REAL(KIND=JPRB)   , INTENT(IN) :: PPRESSURE_HL(0:)
  ! Gas volume mixing ratio climatology (unitless), dimensioned
  ! (latitude, pressure, month)
  REAL(KIND=JPRB)   , INTENT(IN) :: PGAS_VMR(:,:,:)
  ! Output gas layer mass (kg m-2) at the current model time,
  ! dimensioned (latitude, pressure)
  REAL(KIND=JPRB)   , INTENT(OUT):: PGAS_MASS(:,:)

  ! Loop index and maximum index for levels
  INTEGER(KIND=JPIM) :: JLEV, ILEV

  REAL(KIND=JPRB)    :: ZSCALING

  ILEV = UBOUND(PPRESSURE_HL,1)

  DO JLEV = 1,ILEV
    ZSCALING = (PPRESSURE_HL(JLEV)-PPRESSURE_HL(JLEV-1)) * PSCALING
    PGAS_MASS(:,JLEV) = ZSCALING * (PGAS_VMR(:,JLEV,IM2) &
         &                + PTIMI*(PGAS_VMR(:,JLEV,IM1)-PGAS_VMR(:,JLEV,IM2)))
  ENDDO

END SUBROUTINE INTERP_GAS


END SUBROUTINE UPD_GHGCLIM
