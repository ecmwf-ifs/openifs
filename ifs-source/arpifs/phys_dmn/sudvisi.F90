! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUDVISI(YDPHY,YDARPHY,KULOUT)

!**** *SUDVISI*   - Initialize structure YDVISI controlling
!                  constants

!     Purpose.
!     --------
!           Initialize YDVISI, the structure that contains the parameters
!           for the diagnostic of precipitation type.

!**   Interface.
!     ----------
!        *CALL* *SUDVISI(...)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        STRUCTURE YDVISI

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
!      Original : 2019-05-16

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK,    JPHOOK
USE YOMLUN   , ONLY : NULNAM
USE YOMPHY   , ONLY : TPHY
USE YOMARPHY , ONLY : TARPHY

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TPHY) , INTENT(INOUT), TARGET :: YDPHY
TYPE(TARPHY), INTENT(INOUT), TARGET :: YDARPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

REAL(KIND=JPRB) , POINTER :: HVISI
REAL(KIND=JPRB) , POINTER :: COEF_CM1
REAL(KIND=JPRB) , POINTER :: COEF_CM2
REAL(KIND=JPRB) , POINTER :: COEF_CM3
REAL(KIND=JPRB) , POINTER :: COEF_CM4
REAL(KIND=JPRB) , POINTER :: COEF_RM1
REAL(KIND=JPRB) , POINTER :: COEF_RM2
REAL(KIND=JPRB) , POINTER :: COEF_IM1
REAL(KIND=JPRB) , POINTER :: COEF_IM2
REAL(KIND=JPRB) , POINTER :: COEF_SM1
REAL(KIND=JPRB) , POINTER :: COEF_SM2
REAL(KIND=JPRB) , POINTER :: COEF_GM1
REAL(KIND=JPRB) , POINTER :: COEF_GM2

#include "posnam.intfb.h"
#include "namdvisi.nam.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUDVISI',0,ZHOOK_HANDLE)
!Associate for variables not in the include namelists nor allocated in the
!routine 
ASSOCIATE(YDVISI=>YDPHY%YRDVISI, LMPA=>YDARPHY%LMPA)
!     ------------------------------------------------------------------
!include namelists variables, or variables allocated in the routine

HVISI   => YDVISI%HVISI
COEF_CM1 => YDVISI%COEF_CM1
COEF_CM2 => YDVISI%COEF_CM2
COEF_CM3 => YDVISI%COEF_CM3
COEF_CM4 => YDVISI%COEF_CM4
COEF_RM1 => YDVISI%COEF_RM1
COEF_RM2 => YDVISI%COEF_RM2
COEF_IM1 => YDVISI%COEF_IM1
COEF_IM2 => YDVISI%COEF_IM2
COEF_SM1 => YDVISI%COEF_SM1
COEF_SM2 => YDVISI%COEF_SM2
COEF_GM1 => YDVISI%COEF_GM1
COEF_GM2 => YDVISI%COEF_GM2


!*       1.    Set default values.
!              -------------------

!        1.1 Set implicit default values

! Visibilities coefficients
IF (LMPA) THEN  !AROME
   HVISI=5.0_JPRB     
   COEF_CM1=0.07649_JPRB
   COEF_CM2=0.15602_JPRB
   COEF_CM3=0.01937_JPRB
   COEF_CM4=0.92246_JPRB
   COEF_RM1=2.5_JPRB
   COEF_RM2=0.75_JPRB
   COEF_IM1=163.9_JPRB
   COEF_IM2=1_JPRB
   COEF_SM1=10.4_JPRB
   COEF_SM2=0.78_JPRB
   COEF_GM1=2.4_JPRB
   COEF_GM2=0.78_JPRB
ELSE         !ARPEGE
   HVISI= 10.0_JPRB      
   COEF_CM1=0.04412_JPRB
   COEF_CM2=-0.10459_JPRB
   COEF_CM3=-0.0044_JPRB
   COEF_CM4=0.16917_JPRB
   COEF_RM1=2.5_JPRB
   COEF_RM2=0.75_JPRB
   COEF_IM1=163.9_JPRB
   COEF_IM2=1_JPRB
   COEF_SM1=10.4_JPRB
   COEF_SM2=0.78_JPRB
   COEF_GM1=2.4_JPRB
   COEF_GM2=0.78_JPRB
ENDIF



!*       2.    Modify default values.
!              ----------------------

CALL POSNAM(NULNAM,'NAMDVISI')
READ(NULNAM,NAMDVISI)

!     ------------------------------------------------------------------

!*       3.    Print final values.
!              -------------------

WRITE(UNIT=KULOUT,FMT='('' STRUCTURE YDVISI '')')
WRITE(UNIT=KULOUT,FMT='('' HDVISI = '',E10.4,'' COEF_CM1 = '',E10.4 &
 & ,'' COEF_CM2 = '',E10.4,'' COEF_CM3 = '',E10.4 &
 & ,'' COEF_CM4 = '',E10.4,'' COEF_RM1 = '',E10.4 &
 & ,'' COEF_RM2 = '',E10.4,'' COEF_IM1 = '',E10.4 &
 & ,'' COEF_IM2 = '',E10.4,'' COEF_SM1 = '',E10.4 &
 & ,'' COEF_SM2 = '',E10.4,'' COEF_GM1 = '',E10.4 &
 & ,'' COEF_GM2 = '',E10.4 &
 & )')&
 & HVISI,COEF_CM1,COEF_CM2,COEF_CM3,COEF_CM4,COEF_RM1,COEF_RM2,COEF_IM1,COEF_IM2,&
 & COEF_SM1,COEF_SM2,COEF_GM1,COEF_GM2

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUDVISI',1,ZHOOK_HANDLE)
END SUBROUTINE SUDVISI
