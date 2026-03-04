! (C) Copyright 1987- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
SUBROUTINE SUCST(KULOUT,KDAT,KSSS,KPRINTLEV)

!**** *SUCST * - Routine to initialize the constants of the model.

!     Purpose.
!     --------
!           Initialize and print the common YOMCST + initialize
!         date and time of YOMRIP.

!**   Interface.
!     ----------
!        *CALL* *SUCST (..)

!        Explicit arguments :
!        --------------------

!        KULOUT  - logical unit for the output
!        KDAT    - date in the form AAAAMMDD
!        KSSS    - number of seconds in the day
!        KPRINTLEV - printing level

!        Implicit arguments :
!        --------------------
!        COMMON YOMCST
!        COMMON YOMRIP

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-10-15
!        Additions : 90-07-30 (J.-F. Geleyn)
!                    91-11-15 (M. Deque)
!                    96-08-12 M.Hamrud - Reduce printing
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        A. Agusti-Panareda 2020-11-17 Add CO2 molar mass
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRM, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RPI      ,RCLUM    ,RHPLA    ,RKBOL    ,&
 & RNAVO    ,RDAY     ,REA      ,REPSM    ,RSIYEA   ,&
 & RSIDAY   ,ROMEGA   ,RA       ,RG       ,R1SA     ,&
 & RSIGMA   ,RI0      ,R        ,RMD      ,RMV      ,&
 & RMO3     ,RD       ,RV       ,RCPD     ,RCPV     ,&
 & RCVD     ,RCVV     ,RKAPPA   ,RETV     ,RCW      ,&
 & RCS      ,RLVTT    ,RLSTT    ,RLVZER   ,RLSZER   ,&
 & RLMLT    ,RTT      ,RATM     ,RDT      ,RESTT    ,&
 & RALPW    ,RBETW    ,RGAMW    ,RALPS    ,RBETS    ,&
 & RGAMS    ,RALPD    ,RBETD    ,RGAMD    ,RMCO2    ,&
 & RANCO2, RCO2REFYEAR , RNCO2YEARS
USE YOMRIP   , ONLY : RTIMST   ,RTIMTR

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDAT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSSS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPRINTLEV 
INTEGER(KIND=JPIM) :: IA, ID, IDAT, IM, ISSS, J

REAL(KIND=JPRD) :: ZDE, ZET, ZJU, ZRS, ZRSREL, ZTETA, ZTI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "fctast.h"
#include "fcttrm.h"
#include "fcttim.h"
!      -----------------------------------------------------------------

!*       1.    DEFINE FUNDAMENTAL CONSTANTS.
!              -----------------------------

IF (LHOOK) CALL DR_HOOK('SUCST',0,ZHOOK_HANDLE)
RPI=2.0_JPRB*ASIN(1.0_JPRB)
RCLUM=299792458._JPRB
RHPLA=6.6260755E-34_JPRB
RKBOL=1.380658E-23_JPRB
RNAVO=6.0221367E+23_JPRB

!     ------------------------------------------------------------------

!*       2.    DEFINE ASTRONOMICAL CONSTANTS.
!              ------------------------------

RDAY=86400._JPRB
REA=149597870000._JPRB
REPSM=0.409093_JPRB

RSIYEA=365.25_JPRB*RDAY*2.0_JPRB*RPI/6.283076_JPRB
RSIDAY=RDAY/(1.0_JPRB+RDAY/RSIYEA)
ROMEGA=2.0_JPRB*RPI/RSIDAY

IDAT=KDAT
ISSS=KSSS
ID=NDD(IDAT)
IM=NMM(IDAT)
IA=NCCAA(IDAT)
ZJU=RJUDAT(IA,IM,ID)
ZTI=RTIME(IA,IM,ID,ISSS)
RTIMST=ZTI
RTIMTR=ZTI
ZTETA=RTETA(ZTI)
ZRS=RRS(ZTETA)
ZDE=RDS(ZTETA)
ZET=RET(ZTETA)
ZRSREL=ZRS/REA

!     ------------------------------------------------------------------

!*       3.    DEFINE GEOIDE.
!              --------------

RG=9.80665_JPRB
RA=6371229._JPRB
R1SA=REAL(1.0_JPRB/REAL(RA,KIND(1.0_JPRB)),KIND(R1SA))

!     ------------------------------------------------------------------

!*       4.    DEFINE RADIATION CONSTANTS.
!              ---------------------------
RSIGMA=2.0_JPRB * RPI**5 * RKBOL**4 /(15._JPRB* RCLUM**2 * RHPLA**3)
 ! Special case for single precision just activate in this case to keep bit-identical results 
IF ( JPRB == JPRM ) THEN
  RSIGMA=2.0_JPRB * RPI**5 * (RKBOL/RHPLA)**3 * RKBOL /(15._JPRB* RCLUM**2)
ENDIF
  
RI0=1370._JPRB

!     ------------------------------------------------------------------

!*       5.    DEFINE THERMODYNAMIC CONSTANTS, GAS PHASE.
!              ------------------------------------------

R=RNAVO*RKBOL
RMD=28.9644_JPRB
RMCO2=44.0095_JPRB
RMV=18.0153_JPRB
RMO3=47.9942_JPRB
RD=1000._JPRB*R/RMD
RV=1000._JPRB*R/RMV
RCPD=3.5_JPRB*RD
RCVD=RCPD-RD
RCPV=4._JPRB *RV
RCVV=RCPV-RV
RKAPPA=RD/RCPD
RETV=RV/RD-1.0_JPRB



!     ------------------------------------------------------------------

!*       6.    DEFINE THERMODYNAMIC CONSTANTS, LIQUID PHASE.
!              ---------------------------------------------

RCW=4218._JPRB

!     ------------------------------------------------------------------

!*       7.    DEFINE THERMODYNAMIC CONSTANTS, SOLID PHASE.
!              --------------------------------------------

RCS=2106._JPRB

!     ------------------------------------------------------------------

!*       8.    DEFINE THERMODYNAMIC CONSTANTS, TRANSITION OF PHASE.
!              ----------------------------------------------------

RTT=273.16_JPRB
RDT=11.82_JPRB
RLVTT=2.5008E+6_JPRB
RLSTT=2.8345E+6_JPRB
RLVZER=RLVTT+RTT*(RCW-RCPV)
RLSZER=RLSTT+RTT*(RCS-RCPV)
RLMLT=RLSTT-RLVTT
RATM=100000._JPRB

!     ------------------------------------------------------------------

!*       9.    SATURATED VAPOUR PRESSURE.
!              --------------------------

RESTT=611.14_JPRB
RGAMW=(RCW-RCPV)/RV
RBETW=RLVTT/RV+RGAMW*RTT
RALPW=LOG(RESTT)+RBETW/RTT+RGAMW*LOG(RTT)
RGAMS=(RCS-RCPV)/RV
RBETS=RLSTT/RV+RGAMS*RTT
RALPS=LOG(RESTT)+RBETS/RTT+RGAMS*LOG(RTT)
RGAMD=RGAMS-RGAMW
RBETD=RBETS-RBETW
RALPD=RALPS-RALPW

!     ------------------------------------------------------------------

!*       10.    GLOBAL ATMOSPHERIC CO2 FROM CMIP6 (1950-2030) 
!              taken from timeseries file used in radiation (cy47r3)
!              to be used in dtforc as forcing in offline model 
!              for photosynthesis
!              --------------------------------------------------------

RCO2REFYEAR=1950_JPIM
RANCO2=(/312.821_JPRB, 313.014_JPRB, 313.342_JPRB, 313.730_JPRB, 314.095_JPRB, &
       & 314.415_JPRB, 314.698_JPRB, 314.992_JPRB, 315.345_JPRB, 315.807_JPRB, &
       & 316.625_JPRB, 317.299_JPRB, 318.044_JPRB, 318.650_JPRB, 319.333_JPRB, &
       & 319.816_JPRB, 320.880_JPRB, 321.480_JPRB, 322.389_JPRB, 323.251_JPRB, &
       & 324.783_JPRB, 325.400_JPRB, 327.349_JPRB, 329.909_JPRB, 330.756_JPRB, &
       & 330.827_JPRB, 331.545_JPRB, 333.353_JPRB, 335.010_JPRB, 336.605_JPRB, &
       & 338.705_JPRB, 340.059_JPRB, 340.644_JPRB, 342.266_JPRB, 344.008_JPRB, &
       & 345.459_JPRB, 346.903_JPRB, 348.775_JPRB, 351.276_JPRB, 352.894_JPRB, &
       & 354.073_JPRB, 355.353_JPRB, 356.229_JPRB, 356.925_JPRB, 358.254_JPRB, &
       & 360.239_JPRB, 362.005_JPRB, 363.252_JPRB, 365.933_JPRB, 367.845_JPRB, &
       & 369.125_JPRB, 370.673_JPRB, 372.835_JPRB, 375.411_JPRB, 376.987_JPRB, &
       & 378.907_JPRB, 381.010_JPRB, 382.603_JPRB, 384.739_JPRB, 386.280_JPRB, &
       & 388.717_JPRB, 390.944_JPRB, 393.016_JPRB, 395.725_JPRB, 397.547_JPRB, &
       & 399.948_JPRB, 403.117_JPRB, 405.813_JPRB, 408.832_JPRB, 411.946_JPRB, &
       & 415.155_JPRB, 418.461_JPRB, 421.843_JPRB, 425.283_JPRB, 428.784_JPRB, &
       & 432.347_JPRB, 435.975_JPRB, 439.670_JPRB, 443.4319_JPRB, 447.2629_JPRB, 451.1639_JPRB/)


!*      10.    PRINTS

IF (KPRINTLEV >= 1) THEN
  WRITE(KULOUT,'(''0*** Constants of the ICM   ***'')')
  WRITE(KULOUT,'('' *** Fundamental constants ***'')')
  WRITE(KULOUT,'(''           PI = '',E13.7,'' -'')')RPI
  WRITE(KULOUT,'(''            c = '',E13.7,''m s-1'')')RCLUM
  WRITE(KULOUT,'(''            h = '',E13.7,''J s'')')RHPLA
  WRITE(KULOUT,'(''            K = '',E13.7,''J K-1'')')RKBOL
  WRITE(KULOUT,'(''            N = '',E13.7,''mol-1'')')RNAVO
  WRITE(KULOUT,'('' *** Astronomical constants ***'')')
  WRITE(KULOUT,'(''          day = '',E13.7,'' s'')')RDAY
  WRITE(KULOUT,'('' half g. axis = '',E13.7,'' m'')')REA
  WRITE(KULOUT,'('' mean anomaly = '',E13.7,'' -'')')REPSM
  WRITE(KULOUT,'('' sideral year = '',E13.7,'' s'')')RSIYEA
  WRITE(KULOUT,'(''  sideral day = '',E13.7,'' s'')')RSIDAY
  WRITE(KULOUT,'(''        omega = '',E13.7,'' s-1'')')ROMEGA

  WRITE(KULOUT,'('' The initial date of the run is :'')')
  WRITE(KULOUT,'(1X,I8,1X,I5,5X,I4,1X,I2,1X,I2)')IDAT,ISSS,IA,IM,ID
  WRITE(KULOUT,'('' The Julian date is : '',F11.2)') ZJU
  WRITE(KULOUT,'('' Time of the model  : '',F15.2,'' s'')')ZTI
  WRITE(KULOUT,'('' Distance Earth-Sun : '',E13.7,'' m'')')ZRS
  WRITE(KULOUT,'('' Relative Dist. E-S : '',E13.7,'' m'')')ZRSREL
  WRITE(KULOUT,'('' Declination        : '',F12.5)') ZDE
  WRITE(KULOUT,'('' Eq. of time        : '',F12.5,'' s'')')ZET
  WRITE(KULOUT,'('' ***         Geoide         ***'')')
  WRITE(KULOUT,'(''      Gravity = '',E13.7,'' m s-2'')')RG
  WRITE(KULOUT,'('' Earth radius = '',E13.7,'' m'')')RA
  WRITE(KULOUT,'('' Inverse E.R. = '',E13.7,'' m'')')R1SA
  WRITE(KULOUT,'('' ***        Radiation       ***'')')
  WRITE(KULOUT,'('' Stefan-Bol.  = '',E13.7,'' W m-2 K-4'')')  RSIGMA
  WRITE(KULOUT,'('' Solar const. = '',E13.7,'' W m-2'')')RI0
  WRITE(KULOUT,'('' *** Thermodynamic, gas     ***'')')
  WRITE(KULOUT,'('' Perfect gas  = '',e13.7)') R
  WRITE(KULOUT,'('' Dry air mass = '',e13.7)') RMD
  WRITE(KULOUT,'('' Vapour  mass = '',e13.7)') RMV
  WRITE(KULOUT,'('' Ozone   mass = '',e13.7)') RMO3
  WRITE(KULOUT,'('' Dry air cst. = '',e13.7)') RD
  WRITE(KULOUT,'('' Vapour  cst. = '',e13.7)') RV
  WRITE(KULOUT,'(''         Cpd  = '',e13.7)') RCPD
  WRITE(KULOUT,'(''         Cvd  = '',e13.7)') RCVD
  WRITE(KULOUT,'(''         Cpv  = '',e13.7)') RCPV
  WRITE(KULOUT,'(''         Cvv  = '',e13.7)') RCVV
  WRITE(KULOUT,'(''      Rd/Cpd  = '',e13.7)') RKAPPA
  WRITE(KULOUT,'(''     Rv/Rd-1  = '',e13.7)') RETV
  WRITE(KULOUT,'('' *** Thermodynamic, liquid  ***'')')
  WRITE(KULOUT,'(''         Cw   = '',E13.7)') RCW
  WRITE(KULOUT,'('' *** thermodynamic, solid   ***'')')
  WRITE(KULOUT,'(''         Cs   = '',E13.7)') RCS
  WRITE(KULOUT,'('' *** Thermodynamic, trans.  ***'')')
  WRITE(KULOUT,'('' Fusion point  = '',E13.7)') RTT
  WRITE(KULOUT,'('' RTT-Tx(ew-ei) = '',E13.7)') RDT
  WRITE(KULOUT,'(''        RLvTt  = '',E13.7)') RLVTT
  WRITE(KULOUT,'(''        RLsTt  = '',E13.7)') RLSTT
  WRITE(KULOUT,'(''        RLv0   = '',E13.7)') RLVZER
  WRITE(KULOUT,'(''        RLs0   = '',E13.7)') RLSZER
  WRITE(KULOUT,'(''        RLMlt  = '',E13.7)') RLMLT
  WRITE(KULOUT,'('' Normal press. = '',E13.7)') RATM
  WRITE(KULOUT,'('' Latent heat :  '')')
  WRITE(KULOUT,'(10(1X,E10.4))') (10._JPRB*J,J=-4,4)
  WRITE(KULOUT,'(10(1X,E10.4))') (RLV(RTT+10._JPRB*J),J=-4,4)
  WRITE(KULOUT,'(10(1X,E10.4))') (RLS(RTT+10._JPRB*J),J=-4,4)
  WRITE(KULOUT,'('' *** Thermodynamic, satur.  ***'')')
  WRITE(KULOUT,'('' Fusion point = '',E13.7)') RTT
  WRITE(KULOUT,'(''      es(Tt)  = '',e13.7)') RESTT
  WRITE(KULOUT,'('' es(T) :  '')')
  WRITE(KULOUT,'(10(1X,E10.4))') (10._JPRB*J,J=-4,4)
  WRITE(KULOUT,'(10(1X,E10.4))') (ESW(RTT+10._JPRB*J),J=-4,4)
  WRITE(KULOUT,'(10(1X,E10.4))') (ESS(RTT+10._JPRB*J),J=-4,4)
  WRITE(KULOUT,'(10(1X,E10.4))') (ES (RTT+10._JPRB*J),J=-4,4)
  WRITE(KULOUT,'('' CMPI6 atmospheric global mean CO2 :  '')')
  WRITE(KULOUT,'(20(1X,E10.4))') (RANCO2(J),J=60,79)
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUCST',1,ZHOOK_HANDLE)
END SUBROUTINE SUCST

