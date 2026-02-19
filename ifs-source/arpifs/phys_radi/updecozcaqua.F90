! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#ifdef RS6K
@PROCESS NOOPTIMIZE
#endif
!pgi$r opt=0 
SUBROUTINE UPDECOZCAQUA(YDERAD,YDERDI,YDEOZOC,KINDAT,KMINUT)

!**** *UPDECOZCAQUA* - DEFINES A ZONALLY SYMMETRIC LAT-HEIGHT OZONE DISTRIBUTION
!                 SYMMETRIZED ABOUT THE EQUATOR, CORRESPONDING TO THE ANNUAL MEAN
!                 CLIMATOLOGY USED IN AMIP II.

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *UPDECOZCAQUA* FROM *UPDTIM*

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
!     ==== OUTPUTS ===
! ROZT    :                : AMOUNT OF OZONE (KG/KG) 

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        Separated out from SUECOZCAQUA : SUEC does the initial setup, this is called at appropriate time steps from updtim

!     AUTHOR.
!     -------
!     O. Marsden  - January 2018

!     MODIFICATIONS.
!     --------------
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB,  JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RPI
USE YOERDI   , ONLY : TERDI
USE YOEOZOC  , ONLY : TEOZOC
USE YOERAD   , ONLY : TERAD

IMPLICIT NONE

!     -----------------------------------------------------------------

!*       0.1   ARGUMENTS.
!              ----------

TYPE(TERAD)       ,INTENT(IN)    :: YDERAD
TYPE(TERDI)       ,INTENT(IN)    :: YDERDI
TYPE(TEOZOC)      ,INTENT(INOUT) :: YDEOZOC
INTEGER(KIND=JPIM),INTENT(IN)    :: KINDAT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KMINUT 

!*       0.2   LOCAL ARRAYS.
!              -------------

REAL(KIND=JPRB), ALLOCATABLE, SAVE :: ZYTIME(:), ZMDAY(:), ZLAT(:)

INTEGER(KIND=JPIM) :: IDY, IM, IM1, IM2, IMN, JK, JL

REAL(KIND=JPRB) :: ZTIMI, ZXTIME, ZINCH, ZINCL, ZINCROZ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

#include "fcttim.func.h"      
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('UPDECOZCAQUA',0,ZHOOK_HANDLE)
ASSOCIATE(NPERTOZ=>YDERAD%NPERTOZ, RPERTOZ=>YDERAD%RPERTOZ, &
 & RO3=>YDERDI%RO3)

!! Ugly trick to keep this stuff local, but to avoid doing it at each time step. 
!! Should be improvable once IFS time-handling is cleaned up. O Marsden
IF ( .NOT. ALLOCATED(ZMDAY) ) THEN

ALLOCATE(ZMDAY(12))
ALLOCATE(ZYTIME(12))
ALLOCATE(ZLAT(64))

ZMDAY = (/&
 & 31._JPRB,    59.25_JPRB,  90.25_JPRB, 120.25_JPRB, 151.25_JPRB, 181.25_JPRB&
 & ,  212.25_JPRB, 243.25_JPRB, 273.25_JPRB, 304.25_JPRB, 334.25_JPRB, 365.25_JPRB&
 & /)  
ZYTIME= (/&
 & 22320._JPRB,  64980._JPRB, 107640._JPRB, 151560._JPRB, 195480._JPRB, 239400._JPRB&
 & , 283320._JPRB, 327960._JPRB, 371880._JPRB, 415800._JPRB, 459720._JPRB, 503640._JPRB&
 & /)

ZLAT = (/ -8.78638E+01_JPRB, -8.50965E+01_JPRB, -8.23129E+01_JPRB, -7.95256E+01_JPRB, -7.67369E+01_JPRB&
    &,-7.39475E+01_JPRB, -7.11577E+01_JPRB, -6.83678E+01_JPRB, -6.55776E+01_JPRB, -6.27873E+01_JPRB&
    &,-5.99970E+01_JPRB, -5.72066E+01_JPRB, -5.44162E+01_JPRB, -5.16257E+01_JPRB, -4.88352E+01_JPRB&
    &,-4.60447E+01_JPRB, -4.32542E+01_JPRB, -4.04636E+01_JPRB, -3.76731E+01_JPRB, -3.48825E+01_JPRB&
    &,-3.20919E+01_JPRB, -2.93014E+01_JPRB, -2.65108E+01_JPRB, -2.37202E+01_JPRB, -2.09296E+01_JPRB&
    &,-1.81390E+01_JPRB, -1.53484E+01_JPRB, -1.25578E+01_JPRB, -9.76710E+00_JPRB, -6.97650E+00_JPRB&
    &,-4.18590E+00_JPRB, -1.39530E+00_JPRB,  1.39530E+00_JPRB,  4.18590E+00_JPRB,  6.97650E+00_JPRB&
    &, 9.76710E+00_JPRB,  1.25578E+01_JPRB,  1.53484E+01_JPRB,  1.81390E+01_JPRB,  2.09296E+01_JPRB&
    &, 2.37202E+01_JPRB,  2.65108E+01_JPRB,  2.93014E+01_JPRB,  3.20919E+01_JPRB,  3.48825E+01_JPRB&
    &, 3.76731E+01_JPRB,  4.04636E+01_JPRB,  4.32542E+01_JPRB,  4.60447E+01_JPRB,  4.88352E+01_JPRB&
    &, 5.16257E+01_JPRB,  5.44162E+01_JPRB,  5.72066E+01_JPRB,  5.99970E+01_JPRB,  6.27873E+01_JPRB&
    &, 6.55776E+01_JPRB,  6.83678E+01_JPRB,  7.11577E+01_JPRB,  7.39475E+01_JPRB,  7.67369E+01_JPRB&
    &, 7.95256E+01_JPRB,  8.23129E+01_JPRB,  8.50965E+01_JPRB,  8.78638E+01_JPRB /)

ENDIF

!     ------------------------------------------------------------------

!*         1.     TIME INDEX WITHIN OZONE CLIMATOLOGY
!                 -----------------------------------

IDY=NDD(KINDAT)-1
IMN=NMM(KINDAT)
IF (IMN == 1) THEN
  ZXTIME=REAL(IDY,KIND(ZXTIME))*1440._JPRB + KMINUT
ELSEIF (IMN == 2) THEN
  IF(IDY == 28) IDY=IDY-1
! A DAY IN FEB. IS 28.25*24*60/28=1452.8571min LONG.
  ZXTIME=44640._JPRB+REAL(IDY,KIND(ZXTIME))*1452.8571_JPRB+KMINUT
ELSE
  ZXTIME=(ZMDAY(IMN-1)+REAL(IDY,KIND(ZXTIME)))*1440._JPRB+KMINUT
ENDIF
! 525960=MINUTES IN A SIDERAL YEAR (365.25d)
ZXTIME=MOD(ZXTIME,525960._JPRB)

IM1=0
IM2=0
IF (ZXTIME <= ZYTIME(1)) THEN
  IM1=12
  IM2=1
  ZTIMI=(ZYTIME(1)-ZXTIME)/44640._JPRB
ELSEIF(ZXTIME > ZYTIME(12)) THEN
  IM1=12
  IM2=1
  ZTIMI=(548280._JPRB-ZXTIME)/44640._JPRB
! 548280.=(365.25d + 15.5d)*24*60
ELSE
  DO IM=1,11
    IF (ZXTIME > ZYTIME(IM) .AND. ZXTIME <= ZYTIME(IM+1)) THEN
      IM1=IM
      IM2=IM+1
      ZTIMI=(ZXTIME-ZYTIME(IM2))/(ZYTIME(IM1)-ZYTIME(IM2))
    ENDIF
  ENDDO
  IF (IM1 == 0.OR. IM2 == 0 ) THEN
    CALL ABOR1('Problem with time interpolation in UPDECOzc!')
  ENDIF
ENDIF

!*         2.0    TIME INTERPOLATED FIELD
!                 -----------------------

!*( Field is also transformed in kg/kg! )

DO JL=1,64
  DO JK=1,59
    YDEOZOC%ROZT(JL,JK)=RO3 * YDEOZOC%ZOZCLAQUA(JL,JK)
  ENDDO
  YDEOZOC%ROZT(JL, 0)=RO3 * 5.21105E-02_JPRB
  YDEOZOC%ROZT(JL,60)=YDEOZOC%ROZT(JL,59)
ENDDO
ZINCROZ=RPERTOZ
IF (RPERTOZ /= 0.0_JPRB) THEN
  DO JL=6,14
    ZINCL=SIN( FLOAT(JL-6)/4._JPRB*RPI/2._JPRB )
!-- perturbation of Ozone climatology above tropopause (10-100 hPa)
    IF (NPERTOZ == 0) THEN
      DO JK=23,29
        ZINCH=SIN( FLOAT(JK-23)/3.0_JPRB*RPI/2.0_JPRB )
        YDEOZOC%ROZT(JL,JK)=YDEOZOC%ROZT(JL,JK) * (1.0_JPRB + ZINCROZ*ZINCL*ZINCH)
      ENDDO
!-- perturbation of Ozone climatology below tropopause (50-500 hPa)
    ELSEIF (NPERTOZ == 1) THEN
      DO JK=27,33
        ZINCH=SIN( FLOAT(JK-27)/3.0_JPRB*RPI/2.0_JPRB )
        YDEOZOC%ROZT(JL,JK)=YDEOZOC%ROZT(JL,JK) * (1.0_JPRB + ZINCROZ*ZINCL*ZINCH)
      ENDDO
    ENDIF
  ENDDO
ENDIF


!     VECTOR OF LATITUDES FOR OZONE CLIMATOLOGY:

DO JL=1,64
  YDEOZOC%RSINC(JL)=SIN(ZLAT(JL)*RPI/180._JPRB)
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('UPDECOZCAQUA',1,ZHOOK_HANDLE)
END SUBROUTINE UPDECOZCAQUA
