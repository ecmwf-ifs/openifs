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

SUBROUTINE UPDRGAS(YDDYNA,YDERAD,YDERDI,YDRIP,PSOLINC)

!**** *UPDRGAS*

!     PURPOSE.
!     --------

!     Updates the concentrations of the well mixed gases for the
!     radiative calculations, and also the solar forcing. This routine
!     is only called if YDERAD%LHGHG is TRUE.

!**   INTERFACE.
!     ----------

!     CALL UPDRGAS     from *UPDTIM*

!        EXPLICIT ARGUMENTS :
!        --------------------

!     AUTHORS.
!     --------

!       JJMorcrette, ECMWF, 00/01/18

!       JJMorcrette 050119 revised time variations of uniformly mixed gases
!       T Stockdale 2010-04-27  CMIP5 solar forcing data
!       H Hersbach  2011-04-01  CMIP5-recommended RCP greenhouse gases
!       K. Yessad (July 2014): Move some variables.
!       C. Roberts/R. Senan 2017-01-11: CMIP6 forcings 
!       R. Hogan    2017-11-29  INTENT(IN) where appropriate
!       S. Massart   19-Feb-2019 Solar constant optimisation
!       R. Hogan    2019-03-11  Read GHG/TSI timeseries from NetCDF file
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RDAY,  RI0, RMD, RMCO2, RMCH4, RMN2O, RMNO2, RMCFC11, RMCFC12, RMHCFC22, RMCCL4
USE YOMLUN   , ONLY : NULOUT 
USE YOMRIP0  , ONLY : NINDAT
USE YOMDYNA  , ONLY : TDYNA
USE YOMRIP   , ONLY : TRIP
USE YOERAD   , ONLY : TERAD
USE YOERDI   , ONLY : TERDI
USE YOMGHGTIMESERIES,ONLY   : YGHGTIMESERIES   ! GHG multi-annual timeseries
USE YOMSOLARIRRADIANCE,ONLY : YSOLARIRRADIANCE ! TSI multi-annual timeseries

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDYNA),     INTENT(IN)           :: YDDYNA
TYPE(TERAD),     INTENT(IN)           :: YDERAD ! Configuration information
TYPE(TERDI),     INTENT(INOUT)        :: YDERDI ! Output gas concentrations
TYPE(TRIP),      INTENT(IN)           :: YDRIP  ! Time information
REAL(KIND=JPRB), INTENT(IN), OPTIONAL :: PSOLINC ! Solar contant

REAL(KIND=JPRB) :: ZCO2    , ZCH4    , ZN2O    , ZNO2    , ZCFC11    , ZCFC12,     ZHCFC22,     ZCCL4
REAL(KIND=JPRB) :: ZCO2RMWG, ZCH4RMWG, ZN2ORMWG, ZNO2RMWG, ZCFC11RMWG, ZCFC12RMWG, ZHCFC22RMWG, ZCCL4RMWG

! Year as a real number
REAL(KIND=JPRB) :: ZYEAR

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('UPDRGAS',0,ZHOOK_HANDLE)
ASSOCIATE(NHINCSOL=>YDERAD%NHINCSOL, NSCEN=>YDERAD%NSCEN, &
 & RCARDI=>YDERDI%RCARDI, RCFC11=>YDERDI%RCFC11, RCFC12=>YDERDI%RCFC12, &
 & RCH4=>YDERDI%RCH4, RN2O=>YDERDI%RN2O, RNO2=>YDERDI%RNO2, &
 & RHCFC22=>YDERDI%RCFC22, RCCL4=>YDERDI%RCCL4, &
 & RSOLINC=>YDERDI%RSOLINC, YDECMIP=>YDRIP%YRECMIP)
!     ------------------------------------------------------------------

! Get time in decimal years
IF (YDECMIP%NCMIPFIXYR > 0) THEN
  ! Fix time at the centre of the requested year
  ZYEAR = REAL(YDECMIP%NCMIPFIXYR,JPRB) + 0.5_JPRB
ELSE
  ! Get current time
  ZYEAR = GET_YEAR(YDRIP)
ENDIF


!*         1     SURFACE GREENHOUSE GAS CONCENTRATIONS FROM FILE
!                -----------------------------------------------

! Ratio of molar masses with that of dry air to convert volume to mass
! mixing ratio
ZCO2RMWG   = RMCO2   / RMD
ZCH4RMWG   = RMCH4   / RMD
ZN2ORMWG   = RMN2O   / RMD
ZNO2RMWG   = RMNO2   / RMD
ZCFC11RMWG = RMCFC11 / RMD
ZCFC12RMWG = RMCFC12 / RMD
ZHCFC22RMWG= RMHCFC22/ RMD
ZCCL4RMWG  = RMCCL4  / RMD

IF (.NOT. ASSOCIATED(YGHGTIMESERIES)) THEN
  CALL ABOR1('UPDRGAS: Must call SUECRAD first so that YGHGTIMESERIES is allocated')
ENDIF

! Interpolate the values pre-read from NetCDF file to current forecast time
CALL YGHGTIMESERIES%GET(ZYEAR, ZCO2, ZCH4, ZN2O, ZCFC11, ZCFC12, ZHCFC22, ZCCL4)

ZNO2 = 500.E-13_JPRB ! Note that NO2 is not represented in RRTM so it is just a placeholder
  
WRITE(NULOUT,'(a,f0.3,a)') 'UPDRGAS: Surface greenhouse gas concentrations for decimal year ', ZYEAR, ':'
WRITE(NULOUT,'(a,f0.2,a,f0.2,a,f0.2,a,f0.2,a,f0.2,a,f0.2,a,f0.2,a)') 'UPDRGAS:   CO2 = ',ZCO2*1.E6, ' ppmv, CH4 = ',ZCH4*1.E9, &
     &  ' ppbv, N2O = ',ZN2O*1.E9,' ppbv, CFC11 = ',ZCFC11*1.E12,' pptv, CFC12 = ',ZCFC12*1.E12, ' pptv, HCFC22 = ', &
     &  ZHCFC22*1.E12, ' pptv, CCl4 = ', ZCCL4*1.E12, ' pptv'

! Convert to mass mixing ratio
RCARDI  = ZCO2    * ZCO2RMWG
RCH4    = ZCH4    * ZCH4RMWG
RN2O    = ZN2O    * ZN2ORMWG
RCFC11  = ZCFC11  * ZCFC11RMWG
RCFC12  = ZCFC12  * ZCFC12RMWG
RHCFC22 = ZHCFC22 * ZHCFC22RMWG
RCCL4   = ZCCL4   * ZCCL4RMWG
RNO2    = ZNO2    * ZNO2RMWG


!*         2     SOLAR IRRADIANCE
!                ----------------

IF (PRESENT(PSOLINC)) THEN 
  RSOLINC=PSOLINC
  WRITE(NULOUT,*)'UPDRGAS, OVERWRITE RSOLINC WITH ', PSOLINC
ELSE
  IF (NHINCSOL == 0) THEN
    RSOLINC=RI0
  ELSE
    IF (.NOT. ASSOCIATED(YSOLARIRRADIANCE)) THEN
      CALL ABOR1('UPDRGAS: Must call SUECRAD first so that YSOLARIRRADIANCE is allocated')
    ENDIF
    CALL YSOLARIRRADIANCE%GET(ZYEAR, RSOLINC)
    WRITE(NULOUT,'(a,f0.3,a)') 'UPDRGAS:   Total Solar Irradiance: ', RSOLINC, ' W m-2'
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('UPDRGAS',1,ZHOOK_HANDLE)


CONTAINS

  ! ------------------------------------------------------------------
  ! Get the current year as a real number, accurate to the nearest day
  FUNCTION GET_YEAR(YDRIP)

    USE YOMRIP0  , ONLY : NINDAT
    USE YOMCT2   , ONLY : NSTAR2
    USE YOMCT0   , ONLY : LNF 

    ! For time functions NCCAA and NMM
#include "fcttim.func.h"

    TYPE(TRIP),  INTENT(IN)    :: YDRIP  ! Time information
    REAL(KIND=JPRB)            :: GET_YEAR

    ! Integrated number of days through the year at the end of each month
    REAL(KIND=JPRB), PARAMETER :: PP_DAY(*) = &
         & [0.0_JPRB,     31.00_JPRB,  59.25_JPRB,  90.25_JPRB, 120.25_JPRB, 151.25_JPRB, &
         &  181.25_JPRB, 212.25_JPRB, 243.25_JPRB, 273.25_JPRB, 304.25_JPRB, 334.25_JPRB]

    INTEGER(KIND=JPIM) :: ISTADD ! Days since start of forecast
    INTEGER(KIND=JPIM) :: ITIMESTEP, IZT
    INTEGER(KIND=JPIM) :: IYEAR, IMONTH, IDAY

    ! Days since start of year
    REAL(KIND=JPRB) :: ZDAY

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('UPDRGAS:GET_YEAR',0,ZHOOK_HANDLE)

    IF (.NOT.LNF.AND.YDRIP%NSTADD == 0) THEN
      ! IN CASE OF RESTART:
      ITIMESTEP=NINT(YDRIP%TSTEP)
      IF (YDDYNA%LTWOTL) THEN
        IZT=NINT(YDRIP%TSTEP*(REAL(NSTAR2,JPRB)+0.5_JPRB))
      ELSE
        IZT=ITIMESTEP*NSTAR2
      ENDIF
      ISTADD=IZT/NINT(RDAY)
    ELSE
      ISTADD=YDRIP%NSTADD
    ENDIF

    ! Obtain integer year, month and day of start of forecast
    IYEAR  = NCCAA(NINDAT)
    IMONTH = NMM(NINDAT)
    IDAY   = NDD(NINDAT)

    ZDAY = REAL(IDAY,JPRB) - 1.0_JPRB + PP_DAY(IMONTH) + REAL(ISTADD,JPRB)

    GET_YEAR = REAL(IYEAR,JPRB) + ZDAY / 365.25_JPRB

    IF (LHOOK) CALL DR_HOOK('UPDRGAS:GET_YEAR',1,ZHOOK_HANDLE)

  END FUNCTION GET_YEAR

END SUBROUTINE UPDRGAS
