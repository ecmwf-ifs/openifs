! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUCLD(YDECLD,KLEV,PETA)

!**** *SUCLD*   - INITIALIZE COMMON YOECLD CONTROLLING *CLOUD*

!     PURPOSE.
!     --------
!           INITIALIZE YOECLD

!**   INTERFACE.
!     ----------
!        CALL *SUCLD* FROM *SUPHEC*
!              -----        ------

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOECLD

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE
!     "INTEGRATED FORECASTING SYSTEM"

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-12-15
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        R.Forbes      15-Oct-2007 Added RDECORR_CF,_CW
!        P. Marguinaud 04-Oct-2016 Port to single precision

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOECLD   , ONLY : TECLD

IMPLICIT NONE

TYPE(TECLD)       ,INTENT(INOUT) :: YDECLD
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PETA(KLEV) 
!      ----------------------------------------------------------------

INTEGER(KIND=JPIM) :: JK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------

!*       1.    SET VALUES
!              ----------

IF (LHOOK) CALL DR_HOOK('SUCLD',0,ZHOOK_HANDLE)
ASSOCIATE(CETA=>YDECLD%CETA, LOMEGA=>YDECLD%LOMEGA, RANVA=>YDECLD%RANVA, &
 & RANVB=>YDECLD%RANVB, RANVH=>YDECLD%RANVH, RCCA=>YDECLD%RCCA, &
 & RCCB=>YDECLD%RCCB, RCCC=>YDECLD%RCCC, RCFCT=>YDECLD%RCFCT, &
 & RCLWMR=>YDECLD%RCLWMR, RCSCAL=>YDECLD%RCSCAL, RDECORR_CF=>YDECLD%RDECORR_CF, &
 & RDECORR_CW=>YDECLD%RDECORR_CW, REPSCR=>YDECLD%REPSCR, REPSEC=>YDECLD%REPSEC, &
 & RETAHB=>YDECLD%RETAHB, RETAMB=>YDECLD%RETAMB, RGAMMAS=>YDECLD%RGAMMAS, &
 & RLOIA=>YDECLD%RLOIA, RLOIB=>YDECLD%RLOIB, RLOIC=>YDECLD%RLOIC, &
 & RLOID=>YDECLD%RLOID, RLONIA=>YDECLD%RLONIA, RLONIB=>YDECLD%RLONIB, &
 & RRHH=>YDECLD%RRHH, RRHL=>YDECLD%RRHL)
RANVA  = 2._JPRB
RANVB  = 0.3_JPRB
RANVH  = 0.4_JPRB
RCCA   = 0.125_JPRB
RCCB   = 1.5_JPRB
RCCC   = 0.8_JPRB
RCFCT  = 0.400_JPRB
RCSCAL = 1.0E+11_JPRB

! Decorrelation depth scale for generalised cloud overlap (km)
RDECORR_CF = 2.0_JPRB    ! cloud fraction    
RDECORR_CW = 1.0_JPRB    ! condensate

RETAHB = 0.45_JPRB
RETAMB = 0.80_JPRB

RLOIA  = 1.0E+02_JPRB
RLOIB  =-10.00_JPRB
RLOIC  =-0.9_JPRB
RLOID  = 5.0_JPRB

RLONIA = -0.1_JPRB
RLONIB = -10.0_JPRB

RRHH   = 0.9_JPRB
RRHL   = 0.70_JPRB

RGAMMAS= 0.05_JPRB
RCLWMR = 1.E-04_JPRB
LOMEGA =.TRUE.

IF (JPRD == JPRB) THEN
  REPSEC = 1.0E-12_JPRB
ELSE
  REPSEC = 1.0E-6_JPRB
ENDIF

REPSCR = 1.0E-12_JPRB

DO JK=1,KLEV
  CETA(JK)=PETA(JK)
ENDDO

!     -----------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUCLD',1,ZHOOK_HANDLE)
END SUBROUTINE SUCLD
