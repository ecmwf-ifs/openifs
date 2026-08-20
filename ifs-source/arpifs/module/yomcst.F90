! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOMCST

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

TYPE TCST
!*    Common of physical constants
!     You will find the meanings in the annex 1 of the documentation
  
! A1.0 Fundamental constants
! * RPI          : number Pi
! * RCLUM        : light velocity
! * RHPLA        : Planck constant
! * RKBOL        : Bolzmann constant
! * RNAVO        : Avogadro number
  REAL(KIND=JPRB) :: RPI
  REAL(KIND=JPRB) :: RCLUM
  REAL(KIND=JPRB) :: RHPLA
  REAL(KIND=JPRB) :: RKBOL
  REAL(KIND=JPRB) :: RNAVO
  
! A1.1 Astronomical constants
! * RDAY         : duration of the solar day
! * RDAYI        : invariant time unit of 86400s
! * RHOUR        : duration of the solar hour 
! * REA          : astronomical unit (mean distance Earth-sun)
! * REPSM        : polar axis tilting angle
! * RSIYEA       : duration of the sideral year
! * RSIDAY       : duration of the sideral day
! * ROMEGA       : angular velocity of the Earth rotation
  REAL(KIND=JPRB) :: RDAY
  REAL(KIND=JPRB) :: RDAYI
  REAL(KIND=JPRB) :: RHOUR
  REAL(KIND=JPRB) :: REA
  REAL(KIND=JPRB) :: REPSM
  REAL(KIND=JPRB) :: RSIYEA
  REAL(KIND=JPRB) :: RSIDAY
  REAL(KIND=JPRB) :: ROMEGA
  
! A1.2 Geoide
! * RA           : Earth radius
! * RG           : gravity constant
! * R1SA         : 1/RA
! * DEG2RAD      : 180./RPI
  REAL(KIND=JPRB) :: RA
  REAL(KIND=JPRB) :: RG
  REAL(KIND=JPRB) :: R1SA
  REAL(KIND=JPRB) :: DEG2RAD
  
! A1.3 Radiation
! * RSIGMA       : Stefan-Bolzman constant
! * RI0          : solar constant
  REAL(KIND=JPRB) :: RSIGMA
  REAL(KIND=JPRB) :: RI0
  
! A1.4 Thermodynamic gas phase
! * R            : perfect gas constant
! * RMD          : dry air molar mass
! * RMV          : vapour water molar mass
! * RMO3         : ozone molar mass
! * RD           : R_dry (dry air constant)
! * RV           : R_vap (vapour water constant)
! * RCPD         : Cp_dry (dry air calorific capacity at constant pressure)
! * RCPV         : Cp_vap (vapour calorific capacity at constant pressure)
! * RCVD         : Cv_dry (dry air calorific capacity at constant volume)
! * RCVV         : Cv_vap (vapour calorific capacity at constant volume)
! * RKAPPA       : Kappa = R_dry/Cp_dry
! * RETV         : R_vap/R_dry - 1
! * RMCO2        : CO2 (carbon dioxide) molar mass
! * RMCH4        : CH4 (methane) molar mass
! * RMN2O        : N2O molar mass
! * RMCO         : CO (carbon monoxide) molar mass
! * RMHCHO       : HCHO molar mass
! * RMNO2        : NO2 (nitrogen dioxide) molar mass
! * RMSO2        : SO2 (sulfur dioxide) molar mass
! * RMSO4        : SO4 (sulphate) molar mass
! * RMCFC11      : CFC11 molar mass
! * RMCFC12      : CFC12 molar mass
! * RMHCFC12     : HCFC22 molar mass
! * RMCCL4       : CCl4 molar mass
  REAL(KIND=JPRB) :: R
  REAL(KIND=JPRB) :: RMD
  REAL(KIND=JPRB) :: RMV
  REAL(KIND=JPRB) :: RMO3
  REAL(KIND=JPRB) :: RD
  REAL(KIND=JPRB) :: RV
  REAL(KIND=JPRB) :: RCPD
  REAL(KIND=JPRB) :: RCPV
  REAL(KIND=JPRB) :: RCVD
  REAL(KIND=JPRB) :: RCVV
  REAL(KIND=JPRB) :: RKAPPA
  REAL(KIND=JPRB) :: RETV
  REAL(KIND=JPRB) :: RMCO2
  REAL(KIND=JPRB) :: RMCH4
  REAL(KIND=JPRB) :: RMN2O
  REAL(KIND=JPRB) :: RMCO
  REAL(KIND=JPRB) :: RMHCHO
  REAL(KIND=JPRB) :: RMNO2
  REAL(KIND=JPRB) :: RMSO2
  REAL(KIND=JPRB) :: RMSO4
  REAL(KIND=JPRB) :: RMCFC11
  REAL(KIND=JPRB) :: RMCFC12
  REAL(KIND=JPRB) :: RMHCFC22
  REAL(KIND=JPRB) :: RMCCL4
  
! A1.5,6 Thermodynamic liquid,solid phases
! * RCW          : Cw (calorific capacity of liquid water)
! * RCS          : Cs (calorific capacity of solid water)
  REAL(KIND=JPRB) :: RCW
  REAL(KIND=JPRB) :: RCS
  
! A1.7 Thermodynamic transition of phase
! * RATM         : pre_n = "normal" pressure
! * RTT          : Tt = temperature of water fusion at "pre_n"
! * RLVTT        : RLvTt = vaporisation latent heat at T=Tt
! * RLSTT        : RLsTt = sublimation latent heat at T=Tt
! * RLVZER       : RLv0 = vaporisation latent heat at T=0K
! * RLSZER       : RLs0 = sublimation latent heat at T=0K
! * RLMLT        : RLMlt = melting latent heat at T=Tt
! * RDT          : Tt - Tx(ew-ei)
  REAL(KIND=JPRB) :: RATM
  REAL(KIND=JPRB) :: RTT
  REAL(KIND=JPRB) :: RLVTT
  REAL(KIND=JPRB) :: RLSTT
  REAL(KIND=JPRB) :: RLVZER
  REAL(KIND=JPRB) :: RLSZER
  REAL(KIND=JPRB) :: RLMLT
  REAL(KIND=JPRB) :: RDT
  
! A1.8 Curve of saturation
! * RESTT        : es(Tt) = saturation vapour tension at T=Tt
! * RGAMW        : Rgamw = (Cw-Cp_vap)/R_vap
! * RBETW        : Rbetw = RLvTt/R_vap + Rgamw*Tt
! * RALPW        : Ralpw = log(es(Tt)) + Rbetw/Tt + Rgamw*log(Tt)
! * RGAMS        : Rgams = (Cs-Cp_vap)/R_vap
! * RBETS        : Rbets = RLsTt/R_vap + Rgams*Tt
! * RALPS        : Ralps = log(es(Tt)) + Rbets/Tt + Rgams*log(Tt)
! * RALPD        : Ralpd = Ralps - Ralpw
! * RBETD        : Rbetd = Rbets - Rbetw
! * RGAMD        : Rgamd = Rgams - Rgamw
  REAL(KIND=JPRB) :: RESTT
  REAL(KIND=JPRB) :: RGAMW
  REAL(KIND=JPRB) :: RBETW
  REAL(KIND=JPRB) :: RALPW
  REAL(KIND=JPRB) :: RGAMS
  REAL(KIND=JPRB) :: RBETS
  REAL(KIND=JPRB) :: RALPS
  REAL(KIND=JPRB) :: RALPD
  REAL(KIND=JPRB) :: RBETD
  REAL(KIND=JPRB) :: RGAMD
END TYPE TCST

!*    Common of physical constants
!     You will find the meanings in the annex 1 of the documentation

! A1.0 Fundamental constants
! * RPI          : number Pi
! * RCLUM        : light velocity
! * RHPLA        : Planck constant
! * RKBOL        : Bolzmann constant
! * RNAVO        : Avogadro number
REAL(KIND=JPRB),PROTECTED :: RPI
REAL(KIND=JPRB),PROTECTED :: RCLUM
REAL(KIND=JPRB),PROTECTED :: RHPLA
REAL(KIND=JPRB),PROTECTED :: RKBOL
REAL(KIND=JPRB),PROTECTED :: RNAVO

! A1.1 Astronomical constants
! * RDAY         : duration of the solar day
! * RDAYI        : invariant time unit of 86400s
! * RHOUR        : duration of the solar hour 
! * REA          : astronomical unit (mean distance Earth-sun)
! * REPSM        : polar axis tilting angle
! * RSIYEA       : duration of the sideral year
! * RSIDAY       : duration of the sideral day
! * ROMEGA       : angular velocity of the Earth rotation
REAL(KIND=JPRB),PROTECTED :: RDAY
REAL(KIND=JPRB),PROTECTED :: RDAYI
REAL(KIND=JPRB),PROTECTED :: RHOUR
REAL(KIND=JPRB),PROTECTED :: REA
REAL(KIND=JPRB),PROTECTED :: REPSM
REAL(KIND=JPRB),PROTECTED :: RSIYEA
REAL(KIND=JPRB),PROTECTED :: RSIDAY
REAL(KIND=JPRB),PROTECTED :: ROMEGA

! A1.2 Geoid
! * RA           : Planetary radius
! * RA_EARTH     : Earth radius
! * RG           : Gravity constant
! * R1SA         : 1/RA
! * DEG2RAD      : 180./RPI
REAL(KIND=JPRB), PROTECTED :: RA
REAL(KIND=JPRB), PROTECTED :: RA_EARTH
REAL(KIND=JPRB), PROTECTED :: RG
REAL(KIND=JPRB), PROTECTED :: R1SA
REAL(KIND=JPRB), PROTECTED :: DEG2RAD
REAL(KIND=JPRB), PROTECTED :: RAD2DEG

! A1.3 Radiation
! * RSIGMA       : Stefan-Bolzman constant
! * RI0          : solar constant
REAL(KIND=JPRB),PROTECTED :: RSIGMA
REAL(KIND=JPRB),PROTECTED :: RI0

! A1.4 Thermodynamic gas phase
! * R            : perfect gas constant
! * RMD          : dry air molar mass
! * RMV          : vapour water molar mass
! * RMO3         : ozone molar mass
! * RD           : R_dry (dry air constant)
! * RV           : R_vap (vapour water constant)
! * RCPD         : Cp_dry (dry air calorific capacity at constant pressure)
! * RCPV         : Cp_vap (vapour calorific capacity at constant pressure)
! * RCVD         : Cv_dry (dry air calorific capacity at constant volume)
! * RCVV         : Cv_vap (vapour calorific capacity at constant volume)
! * RKAPPA       : Kappa = R_dry/Cp_dry
! * RETV         : R_vap/R_dry - 1
! * RMCO2        : CO2 (carbon dioxide) molar mass
! * RMCH4        : CH4 (methane) molar mass
! * RMN2O        : N2O molar mass
! * RMCO         : CO (carbon monoxide) molar mass
! * RMHCHO       : HCHO molar mass
! * RMNO2        : NO2 (nitrogen dioxide) molar mass
! * RMSO2        : SO2 (sulfur dioxide) molar mass
! * RMSO4        : SO4 (sulphate) molar mass
! * RMCFC11      : CFC11 molar mass
! * RMCFC12      : CFC12 molar mass
! * RMHCFC12     : HCFC22 molar mass
! * RMCCL4       : CCl4 molar mass
REAL(KIND=JPRB),PROTECTED :: R
REAL(KIND=JPRB),PROTECTED :: RMD
REAL(KIND=JPRB),PROTECTED :: RMV
REAL(KIND=JPRB),PROTECTED :: RMO3
REAL(KIND=JPRB),PROTECTED :: RD
REAL(KIND=JPRB),PROTECTED :: RV
REAL(KIND=JPRB),PROTECTED :: RCPD
REAL(KIND=JPRB),PROTECTED :: RCPV
REAL(KIND=JPRB),PROTECTED :: RCVD
REAL(KIND=JPRB),PROTECTED :: RCVV
REAL(KIND=JPRB),PROTECTED :: RKAPPA
REAL(KIND=JPRB),PROTECTED :: RETV
REAL(KIND=JPRB),PROTECTED :: RMCO2
REAL(KIND=JPRB),PROTECTED :: RMCH4
REAL(KIND=JPRB),PROTECTED :: RMN2O
REAL(KIND=JPRB),PROTECTED :: RMCO
REAL(KIND=JPRB),PROTECTED :: RMHCHO
REAL(KIND=JPRB),PROTECTED :: RMNO2
REAL(KIND=JPRB),PROTECTED :: RMSO2
REAL(KIND=JPRB),PROTECTED :: RMSO4
REAL(KIND=JPRB),PROTECTED :: RMCFC11
REAL(KIND=JPRB),PROTECTED :: RMCFC12
REAL(KIND=JPRB),PROTECTED :: RMHCFC22
REAL(KIND=JPRB),PROTECTED :: RMCCL4

! A1.5,6 Thermodynamic liquid,solid phases
! * RCW          : Cw (calorific capacity of liquid water)
! * RCS          : Cs (calorific capacity of solid water)
REAL(KIND=JPRB),PROTECTED :: RCW
REAL(KIND=JPRB),PROTECTED :: RCS

! A1.7 Thermodynamic transition of phase
! * RATM         : pre_n = "normal" pressure
! * RTT          : Tt = temperature of water fusion at "pre_n"
! * RLVTT        : RLvTt = vaporisation latent heat at T=Tt
! * RLSTT        : RLsTt = sublimation latent heat at T=Tt
! * RLVZER       : RLv0 = vaporisation latent heat at T=0K
! * RLSZER       : RLs0 = sublimation latent heat at T=0K
! * RLMLT        : RLMlt = melting latent heat at T=Tt
! * RDT          : Tt - Tx(ew-ei)
REAL(KIND=JPRB),PROTECTED :: RATM
REAL(KIND=JPRB),PROTECTED :: RTT
REAL(KIND=JPRB),PROTECTED :: RLVTT
REAL(KIND=JPRB),PROTECTED :: RLSTT
REAL(KIND=JPRB),PROTECTED :: RLVZER
REAL(KIND=JPRB),PROTECTED :: RLSZER
REAL(KIND=JPRB),PROTECTED :: RLMLT
REAL(KIND=JPRB),PROTECTED :: RDT

! A1.8 Curve of saturation
! * RESTT        : es(Tt) = saturation vapour tension at T=Tt
! * RGAMW        : Rgamw = (Cw-Cp_vap)/R_vap
! * RBETW        : Rbetw = RLvTt/R_vap + Rgamw*Tt
! * RALPW        : Ralpw = log(es(Tt)) + Rbetw/Tt + Rgamw*log(Tt)
! * RGAMS        : Rgams = (Cs-Cp_vap)/R_vap
! * RBETS        : Rbets = RLsTt/R_vap + Rgams*Tt
! * RALPS        : Ralps = log(es(Tt)) + Rbets/Tt + Rgams*log(Tt)
! * RALPD        : Ralpd = Ralps - Ralpw
! * RBETD        : Rbetd = Rbets - Rbetw
! * RGAMD        : Rgamd = Rgams - Rgamw
REAL(KIND=JPRB),PROTECTED :: RESTT
REAL(KIND=JPRB),PROTECTED :: RGAMW
REAL(KIND=JPRB),PROTECTED :: RBETW
REAL(KIND=JPRB),PROTECTED :: RALPW
REAL(KIND=JPRB),PROTECTED :: RGAMS
REAL(KIND=JPRB),PROTECTED :: RBETS
REAL(KIND=JPRB),PROTECTED :: RALPS
REAL(KIND=JPRB),PROTECTED :: RALPD
REAL(KIND=JPRB),PROTECTED :: RBETD
REAL(KIND=JPRB),PROTECTED :: RGAMD

! NaN value
CHARACTER(LEN=8), PARAMETER :: CSNAN = &
  & CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//CHAR(244)//CHAR(127)
REAL(KIND=JPRB),PROTECTED :: RSNAN

LOGICAL,PROTECTED :: L_HAS_BEEN_SETUP=.FALSE.

TYPE (TCST), TARGET :: YRCST


CONTAINS 

SUBROUTINE TCST_INIT(YDCST)

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   JPHOOK, DR_HOOK

IMPLICIT NONE

TYPE(TCST), INTENT(IN) :: YDCST

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMCST:TCST_INIT',0,ZHOOK_HANDLE)

  RPI       =   YDCST%RPI      
  RCLUM     =   YDCST%RCLUM    
  RHPLA     =   YDCST%RHPLA    
  RKBOL     =   YDCST%RKBOL    
  RNAVO     =   YDCST%RNAVO    
  RDAY      =   YDCST%RDAY     
  RDAYI     =   YDCST%RDAYI    
  RHOUR     =   YDCST%RHOUR    
  REA       =   YDCST%REA      
  REPSM     =   YDCST%REPSM    
  RSIYEA    =   YDCST%RSIYEA   
  RSIDAY    =   YDCST%RSIDAY   
  ROMEGA    =   YDCST%ROMEGA   
  RA        =   YDCST%RA       
  RG        =   YDCST%RG       
  R1SA      =   YDCST%R1SA     
  DEG2RAD   =   YDCST%DEG2RAD  
  RSIGMA    =   YDCST%RSIGMA   
  RI0       =   YDCST%RI0      
  R         =   YDCST%R        
  RMD       =   YDCST%RMD      
  RMV       =   YDCST%RMV      
  RMO3      =   YDCST%RMO3     
  RD        =   YDCST%RD       
  RV        =   YDCST%RV       
  RCPD      =   YDCST%RCPD     
  RCPV      =   YDCST%RCPV     
  RCVD      =   YDCST%RCVD     
  RCVV      =   YDCST%RCVV     
  RKAPPA    =   YDCST%RKAPPA   
  RETV      =   YDCST%RETV     
  RMCO2     =   YDCST%RMCO2    
  RMCH4     =   YDCST%RMCH4    
  RMN2O     =   YDCST%RMN2O    
  RMCO      =   YDCST%RMCO     
  RMHCHO    =   YDCST%RMHCHO   
  RMNO2     =   YDCST%RMNO2    
  RMSO2     =   YDCST%RMSO2    
  RMSO4     =   YDCST%RMSO4    
  RMCFC11   =   YDCST%RMCFC11  
  RMCFC12   =   YDCST%RMCFC12  
  RMHCFC22  =   YDCST%RMHCFC22 
  RMCCL4    =   YDCST%RMCCL4   
  RCW       =   YDCST%RCW      
  RCS       =   YDCST%RCS      
  RATM      =   YDCST%RATM     
  RTT       =   YDCST%RTT      
  RLVTT     =   YDCST%RLVTT    
  RLSTT     =   YDCST%RLSTT    
  RLVZER    =   YDCST%RLVZER   
  RLSZER    =   YDCST%RLSZER   
  RLMLT     =   YDCST%RLMLT    
  RDT       =   YDCST%RDT      
  RESTT     =   YDCST%RESTT    
  RGAMW     =   YDCST%RGAMW    
  RBETW     =   YDCST%RBETW    
  RALPW     =   YDCST%RALPW    
  RGAMS     =   YDCST%RGAMS    
  RBETS     =   YDCST%RBETS    
  RALPS     =   YDCST%RALPS    
  RALPD     =   YDCST%RALPD    
  RBETD     =   YDCST%RBETD    
  RGAMD     =   YDCST%RGAMD    

IF (LHOOK) CALL DR_HOOK('YOMCST:TCST_INIT',1,ZHOOK_HANDLE)

END SUBROUTINE TCST_INIT

SUBROUTINE SETUP_CONSTANTS(KULOUT,KPRINTLEV)

!**** *SETUP_CONSTANTS * - Routine to initialize the constants of the model.

!     Purpose.
!     --------
!           Initialize and print the common YOMCST + initialize
!         date and time of YOMRIP0.

!**   Interface.
!     ----------
!        *CALL* *SETUP_CONSTANTS (..)

!        Explicit arguments :
!        --------------------

!        KULOUT  - logical unit for the output
!        KPRINTLEV - printing level

!        Implicit arguments :
!        --------------------

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
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      A.Alias   : 07-29-11 read values from NAMSCEN
!      A.Voldoire: 09-08    YOMSCEN removed as NAMSCEN has been modified
!      A.Voldoire: 11-03    LASTRF to prevent any drift in insolation
!      N. Semane+P.Bechtold   add RDAYI/RHOUR
!      K. Yessad (July 2014): Move some variables.
!      F. Vana  05-Mar-2015  Support for single precision
!      M.Hamrud    : from SUCST
!      R. Hogan 25-Jan-2019  Added molar masses of CFC11, CFC12, HCFC22 and CCl4
!      R. El Khatib 08-Jul-2022 Contribution to the encapsulation of YOMCST and YOETHF
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULNAM
USE YOMDYNCORE,ONLY : RPLRADI  ,RCORIOI   ,LAPE    ,RPLRG

!      -----------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPRINTLEV 

!      -----------------------------------------------------------------

INTEGER(KIND=JPIM) :: J
TYPE(TCST) :: YDCST ! allows use of included functions. Sorry Doctor. REK.

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      -----------------------------------------------------------------

NAMELIST/NAMSCEN/ RI0

#include "posnam.intfb.h"
#include "abor1.intfb.h"

#include "fcttrm.func.h"

!      -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('YOMCST:SETUP_CONSTANTS',0,ZHOOK_HANDLE)
ASSOCIATE(RPI=>YDCST%RPI, RCLUM=>YDCST%RCLUM, RHPLA=>YDCST%RHPLA, RKBOL=>YDCST%RKBOL, RNAVO=>YDCST%RNAVO, RDAY=>YDCST%RDAY, &
 & RDAYI=>YDCST%RDAYI, RHOUR=>YDCST%RHOUR, REA=>YDCST%REA, REPSM=>YDCST%REPSM, RSIYEA=>YDCST%RSIYEA, RSIDAY=>YDCST%RSIDAY, &
 & ROMEGA=>YDCST%ROMEGA, RA=>YDCST%RA, RG=>YDCST%RG, R1SA=>YDCST%R1SA, DEG2RAD=>YDCST%DEG2RAD, RSIGMA=>YDCST%RSIGMA, &
 & RI0=>YDCST%RI0, R=>YDCST%R, RMD=>YDCST%RMD, RMV=>YDCST%RMV, RMO3=>YDCST%RMO3, RD=>YDCST%RD, RV=>YDCST%RV, RCPD=>YDCST%RCPD, &
 & RCPV=>YDCST%RCPV, RCVD=>YDCST%RCVD, RCVV=>YDCST%RCVV, RKAPPA=>YDCST%RKAPPA, RETV=>YDCST%RETV, RMCO2=>YDCST%RMCO2, &
 & RMCH4=>YDCST%RMCH4, RMN2O=>YDCST%RMN2O, RMCO=>YDCST%RMCO, RMHCHO=>YDCST%RMHCHO, RMNO2=>YDCST%RMNO2, RMSO2=>YDCST%RMSO2, &
 & RMSO4=>YDCST%RMSO4, RMCFC11=>YDCST%RMCFC11, RMCFC12=>YDCST%RMCFC12, RMHCFC22=>YDCST%RMHCFC22, RMCCL4=>YDCST%RMCCL4, &
 & RCW=>YDCST%RCW, RCS=>YDCST%RCS, RATM=>YDCST%RATM, RTT=>YDCST%RTT, RLVTT=>YDCST%RLVTT, RLSTT=>YDCST%RLSTT, RLVZER=>YDCST%RLVZER, &
 & RLSZER=>YDCST%RLSZER, RLMLT=>YDCST%RLMLT, RDT=>YDCST%RDT, RESTT=>YDCST%RESTT, RGAMW=>YDCST%RGAMW, RBETW=>YDCST%RBETW, &
 & RALPW=>YDCST%RALPW, RGAMS=>YDCST%RGAMS, RBETS=>YDCST%RBETS, RALPS=>YDCST%RALPS, RALPD=>YDCST%RALPD, RBETD=>YDCST%RBETD, &
 & RGAMD=>YDCST%RGAMD)
!      -----------------------------------------------------------------

IF(L_HAS_BEEN_SETUP) CALL ABOR1('YOMCST:SETUP_CONSTANTS - ALREADY CALLED ONCE')
!*       1.    DEFINE FUNDAMENTAL CONSTANTS.
!              -----------------------------

RPI=2.0_JPRB*ASIN(1.0_JPRB)
RCLUM=299792458._JPRB
RHPLA=6.6260755E-34_JPRB
RKBOL=1.380658E-23_JPRB
RNAVO=6.0221367E+23_JPRB

!     ------------------------------------------------------------------

!*       2.    DEFINE ASTRONOMICAL CONSTANTS.
!              ------------------------------

RDAYI=86400._JPRB
RDAY=86400._JPRB
RHOUR=3600._JPRB
IF (RCORIOI /= 0.0_JPRB) THEN
RDAY=86400._JPRB*RCORIOI
RHOUR=3600._JPRB*RCORIOI
ENDIF 
REA=149597870000._JPRB
IF( LAPE ) THEN
  ! aqua-planet special, obliquity == 0.0
  REPSM=0.0_JPRB
ELSE
  REPSM=0.409093_JPRB
ENDIF

RSIYEA=365.25_JPRB*RDAY*2.0_JPRB*RPI/6.283076_JPRB
RSIDAY=RDAY/(1.0_JPRB+RDAY/RSIYEA)
IF (RCORIOI /= 0.0_JPRB) THEN
ROMEGA=2.0_JPRB*RPI/RSIDAY
ELSE
ROMEGA=0.0_JPRB
ENDIF

!     ------------------------------------------------------------------

!*       3.    DEFINE GEOIDE.
!              --------------

RA_EARTH=6371229._JPRB
RA=RA_EARTH*RPLRADI
IF( LAPE ) THEN
  RG=9.79764_JPRB*RPLRG
! intercomparison RA=6371000.
ELSE
  RG=9.80665_JPRB
ENDIF
R1SA=REAL(1.0_JPRB/REAL(RA,KIND(1.0_JPRB)),KIND(R1SA))
DEG2RAD=RPI/180.0_JPRB
RAD2DEG=180._JPRB/RPI
!     ------------------------------------------------------------------

!*       4.    DEFINE RADIATION CONSTANTS.
!              ---------------------------

!RSIGMA=2.0_JPRB * RPI**5 * RKBOL**4 /(15._JPRB* RCLUM**2 * RHPLA**3)
RSIGMA=2.0_JPRB * RPI**5 * (RKBOL/RHPLA)**3 * RKBOL /(15._JPRB* RCLUM**2) 
IF( LAPE ) THEN
  RI0=1365._JPRB
ELSE
!  RI0=1370._JPRB ! pre-CY35R3 value
  RI0=1366._JPRB
ENDIF

CALL POSNAM(NULNAM,'NAMSCEN')
READ       (NULNAM, NAMSCEN)

!     ------------------------------------------------------------------

!*       5.    DEFINE THERMODYNAMIC CONSTANTS, GAS PHASE.
!              ------------------------------------------

R=RNAVO*RKBOL
RMD=28.9644_JPRB
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
RMCO2=44.0095_JPRB
RMCH4=16.04_JPRB
RMN2O=44.013_JPRB
RMCO=28.01_JPRB
RMHCHO=30.03_JPRB
RMNO2=46.01_JPRB
RMSO2=64.056_JPRB
RMSO4=96.052_JPRB
RMCFC11=137.3686_JPRB
RMCFC12=120.914_JPRB
RMHCFC22=86.469_JPRB
RMCCL4=153.823_JPRB
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

  WRITE(KULOUT,'('' ***         Geoid          ***'')')
  WRITE(KULOUT,'(''      Gravity = '',E13.7,'' m s-2'')')RG
  WRITE(KULOUT,'('' Earth radius = '',E13.7,'' m'')')RA_EARTH
  WRITE(KULOUT,'('' Plan. radius = '',E13.7,'' m'')')RA
  WRITE(KULOUT,'('' Inverse P.R. = '',E13.7,'' m'')')R1SA
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
ENDIF

! NAN value
RSNAN = TRANSFER (CSNAN, RSNAN)
L_HAS_BEEN_SETUP = .TRUE.

! Initialize module structure (it would be better to pass YRCST in argument to SETUP_CONSTANTS) :
YRCST = YDCST

! Initialize module individual variables from module structure :
END ASSOCIATE
CALL TCST_INIT(YDCST)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('YOMCST:SETUP_CONSTANTS',1,ZHOOK_HANDLE)
END SUBROUTINE SETUP_CONSTANTS


SUBROUTINE RVEQUALSRD
USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
IMPLICIT NONE
!very dirty code used at MF with the switch LDRYTL. See cva1.F90
!Need to be checked and hopefully removed soon...
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('YOMCST:RVEQUALSRD',0,ZHOOK_HANDLE)
RV=RD
YRCST%RV=YRCST%RD
IF (LHOOK) CALL DR_HOOK('YOMCST:RVEQUALSRD',1,ZHOOK_HANDLE)
END SUBROUTINE RVEQUALSRD
!    ------------------------------------------------------------------
END MODULE YOMCST
