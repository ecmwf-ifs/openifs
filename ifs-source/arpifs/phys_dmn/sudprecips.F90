! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUDPRECIPS(YDPHY,YDRIP,KULOUT)

!**** *SUDPRECIPS*   - Initialize structure YDPRECIPS controlling
!                  constants

!     Purpose.
!     --------
!           Initialize YDPRECIPS, the structure that contains the parameters
!           for the diagnostic of precipitation type.

!**   Interface.
!     ----------
!        *CALL* *SUDPRECIPS(...)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        STRUCTURE YDPRECIPS

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Documentation ARPEGE

!     Author.
!     -------
!      I.Etchevers .
!      Original : 2018-07-12

!     Modifications.
!     --------------
!      Y. Seity : 2018-07-17 add RDHAIL* 
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK,    JPHOOK
USE YOMLUN   , ONLY : NULNAM
USE YOMPHY   , ONLY : TPHY
USE YOMRIP   , ONLY : TRIP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TPHY) , INTENT(INOUT), TARGET :: YDPHY
TYPE(TRIP) , INTENT(INOUT), TARGET :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

REAL(KIND=JPRB), POINTER ::  HDPRECIPS
REAL(KIND=JPRB), POINTER ::  HDCLWC
REAL(KIND=JPRB), POINTER ::  RDHAIL1
REAL(KIND=JPRB), POINTER ::  RDHAIL2
REAL(KIND=JPRB), POINTER ::  RDSEUIL1
REAL(KIND=JPRB), POINTER ::  RDSEUIL2
REAL(KIND=JPRB), POINTER ::  RDSEUIL3
REAL(KIND=JPRB), POINTER ::  RDSEUIL4
REAL(KIND=JPRB), POINTER ::  RDSEUIL5
REAL(KIND=JPRB), POINTER ::  RDCLWC
REAL(KIND=JPRB), POINTER ::  RPRECSEUIL
REAL(KIND=JPRB), POINTER ::  RHTOP
REAL(KIND=JPRB), POINTER ::  RTPW
REAL(KIND=JPRB), POINTER ::  RAWARM
REAL(KIND=JPRB), POINTER ::  RACOLD
INTEGER(KIND=JPIM), POINTER ::  NDTPREC
INTEGER(KIND=JPIM), POINTER ::  NDTPREC2
INTEGER(KIND=JPIM), POINTER ::  NDTPRECCUR
INTEGER(KIND=JPIM), POINTER ::  NDTPRECCUR2

#include "posnam.intfb.h"
#include "namdprecips.nam.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUDPRECIPS',0,ZHOOK_HANDLE)
!Associate for variables not in the include namelists nor allocated in the
!routine 
ASSOCIATE(YDPRECIPS=>YDPHY%YRDPRECIPS,NDPRECPERIOD=>YDPHY%NDPRECPERIOD,NDPRECPERIOD2=>YDPHY%NDPRECPERIOD2,&
                & TSTEP=>YDRIP%TSTEP)
!     ------------------------------------------------------------------
!include namelists variables, or variables allocated in the routine
HDPRECIPS  => YDPRECIPS%HDPRECIPS
HDCLWC     => YDPRECIPS%HDCLWC
RDHAIL1    => YDPRECIPS%RDHAIL1
RDHAIL2    => YDPRECIPS%RDHAIL2
RDSEUIL1   => YDPRECIPS%RDSEUIL1
RDSEUIL2   => YDPRECIPS%RDSEUIL2
RDSEUIL3   => YDPRECIPS%RDSEUIL3
RDSEUIL4   => YDPRECIPS%RDSEUIL4
RDSEUIL5   => YDPRECIPS%RDSEUIL5
RDCLWC     => YDPRECIPS%RDCLWC
RPRECSEUIL => YDPRECIPS%RPRECSEUIL
RHTOP      => YDPRECIPS%RHTOP
RTPW       => YDPRECIPS%RTPW
RAWARM     => YDPRECIPS%RAWARM
RACOLD     => YDPRECIPS%RACOLD
NDTPREC    => YDPRECIPS%NDTPREC
NDTPREC2   => YDPRECIPS%NDTPREC2
NDTPRECCUR => YDPRECIPS%NDTPRECCUR
NDTPRECCUR2=> YDPRECIPS%NDTPRECCUR2

!*       1.    Set default values.
!              -------------------

!        1.1 Set implicit default values

HDPRECIPS=2.0_JPRB
HDCLWC=100._JPRB
RDHAIL1=8._JPRB
RDHAIL2=16._JPRB
RDSEUIL1=0.001_JPRB
RDSEUIL2=0.2_JPRB
RDSEUIL3=0.6_JPRB
RDSEUIL4=0.8_JPRB
RDSEUIL5=0.0003_JPRB !in mm/s = 1mm/h rain/drizzle
RDCLWC=0.00035_JPRB
RPRECSEUIL=0.00003_JPRB !in mm/s = 0.1mm/h
RHTOP=4000._JPRB
RTPW=273.15_JPRB
RACOLD=4500._JPRB
RAWARM=400._JPRB
NDTPREC=INT(NDPRECPERIOD/TSTEP)
NDTPREC2=INT(NDPRECPERIOD2/TSTEP)
NDTPRECCUR=1_JPIM
NDTPRECCUR2=1_JPIM

!*       2.    Modify default values.
!              ----------------------

CALL POSNAM(NULNAM,'NAMDPRECIPS')
READ(NULNAM,NAMDPRECIPS)

!     ------------------------------------------------------------------

!*       3.    Print final values.
!              -------------------

WRITE(UNIT=KULOUT,FMT='('' STRUCTURE YDPRECIPS '')')
WRITE(UNIT=KULOUT,FMT='('' HDPRECIPS = '',E10.4,'' HDCLWC = '',E10.4 &
 & ,'' RDHAIL1 = '',E10.4,'' RDHAIL2 = '',E10.4 &
 & ,'' RDSEUIL1 = '',E10.4,'' RDSEUIL2 = '',E10.4 &
 & ,'' RDSEUIL3 = '',E10.4,'' RDSEUIL4 = '',E10.4 &
 & ,'' RDSEUIL5 = '',E10.4 &
 & ,'' RDCLWC = '',E10.4,'' RPRECSEUIL = '',E10.4&
 & ,'' RHTOP = '',E10.4,'' RTPW = '',E10.4&
 & ,'' RAWARM = '',E10.4,'' RACOLD = '',E10.4&
 & ,'' NDTPREC = '',I8,'' NDTPRECCUR  = '',I8&
 & )')&
 & HDPRECIPS,HDCLWC,RDHAIL1,RDHAIL2,RDSEUIL1,RDSEUIL2,RDSEUIL3,RDSEUIL4,&
 & RDSEUIL5,RDCLWC,RPRECSEUIL,RHTOP,RTPW,RAWARM,RACOLD,NDTPREC,NDTPRECCUR

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUDPRECIPS',1,ZHOOK_HANDLE)
END SUBROUTINE SUDPRECIPS
