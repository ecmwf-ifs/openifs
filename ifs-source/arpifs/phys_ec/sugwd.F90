! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUGWD(YDSTA,YDEGWD,YDEPHY,KULOUT,KLEV,PVAH,PVBH)

!**** *SUGWD* INITIALIZE COMMON YOEGWD CONTROLLING GRAVITY WAVE DRAG

!     PURPOSE.
!     --------
!           INITIALIZE YOEGWD, THE COMMON THAT CONTROLS THE
!           GRAVITY WAVE DRAG PARAMETRIZATION.

!**   INTERFACE.
!     ----------
!        CALL *SUGWD* FROM *SUPHEC*
!              -----        ------

!        EXPLICIT ARGUMENTS :
!        --------------------
!        KULOUT      : LOGICAL UNIT FOR THE OUTPUT
!        PVAH,PVBH   : VERTICAL COORDINATE TABLE
!        KLEV        : NUMBER OF MODEL LEVELS

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOEGWD

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        MARTIN MILLER             *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 90-01-01       ALSO : 95-01-20
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P.Bechtold    14-Jan-2009 precompute levels for RF
!        P.Bechtold    04-Oct-2010 simplify precomp of NGWDLIM
!        N.Semane+P.Bechtold    04-10-2010 replace 3600s by RHOUR for small planet
!        I. Sandu      15-03-2013  adjustment of GKWAKE parameter
!        R. El Khatib  10-Aug-2011 More proper abort check when EC physics is inactive
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE YOMSTA   , ONLY : TSTA
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOEPHY   , ONLY : TEPHY
USE YOEGWD   , ONLY : TEGWD
USE YOMLUN   , ONLY : NULNAM
USE YOMCST   , ONLY : RHOUR

IMPLICIT NONE

TYPE(TSTA)        ,INTENT(IN)   :: YDSTA
TYPE(TEGWD),TARGET,INTENT(INOUT):: YDEGWD
TYPE(TEPHY)       ,INTENT(INOUT):: YDEPHY
INTEGER(KIND=JPIM),INTENT(IN)   :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)   :: KULOUT 
REAL(KIND=JPRB)   ,INTENT(IN)   :: PVAH(KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)   :: PVBH(KLEV+1) 
!      ----------------------------------------------------------------

INTEGER(KIND=JPIM) :: JK

REAL(KIND=JPRB) :: ZPM1R(KLEV), ZPR, ZSIGT, ZPLIM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


REAL(KIND=JPRB), POINTER :: GTENLIM
REAL(KIND=JPRB), POINTER :: GKWAKE 
REAL(KIND=JPRB), POINTER :: GKDRAG
LOGICAL, POINTER :: LRDIFF_STRATO, LDIAG_STRATO
INTEGER, POINTER :: NDIFF_STRATO

#include "namgwd.nam.h"

!#include "abor1.intfb.h"
#include "posnam.intfb.h"


!*       1.    SET THE VALUES OF THE PARAMETERS
!              --------------------------------

IF (LHOOK) CALL DR_HOOK('SUGWD',0,ZHOOK_HANDLE)
ASSOCIATE(GFRCRIT=>YDEGWD%GFRCRIT,  &
 & GRCRIT=>YDEGWD%GRCRIT, GRFPLM=>YDEGWD%GRFPLM, GSSEC=>YDEGWD%GSSEC, &
 & GTSEC=>YDEGWD%GTSEC, GVSEC=>YDEGWD%GVSEC, NGWDLIM=>YDEGWD%NGWDLIM, &
 & NGWDTOP=>YDEGWD%NGWDTOP, NKTOPG=>YDEGWD%NKTOPG, &
 & LEPHYS=>YDEPHY%LEPHYS, &
 & STPRE=>YDSTA%STPRE)
! Associate pointers for variables in namelist
GTENLIM       => YDEGWD%GTENLIM
LRDIFF_STRATO => YDEGWD%LRDIFF_STRATO
LDIAG_STRATO  => YDEGWD%LDIAG_STRATO
NDIFF_STRATO  => YDEGWD%NDIFF_STRATO
GKDRAG        => YDEGWD%GKDRAG
GKWAKE        => YDEGWD%GKWAKE


ZSIGT=0.94_JPRB
ZPR  =80000._JPRB

!  As a security measure when running with few levels,
!  we force NKTOPG=KLEV to make sure it is defined

NKTOPG=KLEV

DO JK=KLEV,1,-1
  ZPM1R(JK)=0.5_JPRB*(PVAH(JK)+PVAH(JK+1)+ZPR*(PVBH(JK)+PVBH(JK+1)))
  IF((ZPM1R(JK)/ZPR) >= ZSIGT)THEN
    NKTOPG=JK
  ENDIF
ENDDO

GRFPLM=9.9_JPRB
NGWDTOP=1
DO JK=1,KLEV
  IF(ZPM1R(JK)<GRFPLM)THEN
    NGWDTOP=JK
  ENDIF
ENDDO

GKDRAG =0.3_JPRB
GKWAKE =1._JPRB
!  Revised gwd parameter values
GKDRAG =0.15_JPRB
!
!GKWAKE =2._JPRB
!38r2
!GKWAKE =1.3_JPRB

GKWAKE =3._JPRB

GRCRIT =0.25_JPRB
GFRCRIT=0.40_JPRB

!      ----------------------------------------------------------------

!*       2.    SET VALUES OF SECURITY PARAMETERS
!              ---------------------------------

GVSEC=0.10_JPRB
GSSEC=1.E-12_JPRB

GTSEC=1.E-07_JPRB

!      ----------------------------------------------------------------

!*       3.    SET VALUES OF PARAMETERS FOR LIMITING
!*             THE WIND TENDENCIES IN STRATOSPHERE AND MESOSPHERE

! WIND TENDENCIES LIMITED ABOVE 50hPa
IF(KLEV > 19)THEN
  ZPLIM=5000._JPRB
ELSE
  ZPLIM=1000._JPRB
ENDIF

NGWDLIM=0
DO JK=KLEV,1,-1
  IF(STPRE(JK) >= ZPLIM) THEN
    NGWDLIM=JK
  ENDIF
ENDDO


! * Reset NGWDLIM in order to have NGWDLIM+1 < 1 if ECMWF physics is not called,
!   because the following test on PVBH(NGWDLIM+1) must NEVER generate an
!   abor1 when ECMWF physics is not called.
IF( .NOT.LEPHYS ) THEN
  NGWDLIM=0
ENDIF

WRITE(KULOUT,*) ' SUGWD: NGWDLIM=',NGWDLIM

! WIND TENDENCIES LIMITED TO LESS OR EQUAL TO 20.m/s per hour
!  Recommend a value of at least 80.m/s per hour
GTENLIM=20.0_JPRB/RHOUR
GTENLIM=80.0_JPRB/RHOUR

! Reduced vertical diffusion in stratosphere
LRDIFF_STRATO=.FALSE.
NDIFF_STRATO=5
! Diagnostics: stratosphere wave fluxes
LDIAG_STRATO=.FALSE.

CALL POSNAM(NULNAM,'NAMGWD')
READ (NULNAM,NAMGWD)

WRITE(KULOUT,*) ' SUGWD: GTENLIM=',GTENLIM, ' SUGWD: LRDIFF_STRATO=',LRDIFF_STRATO
WRITE(KULOUT,*) ' SUGWD: NDIFF_STRATO=',NDIFF_STRATO, 'SUGWD: GKWAKE=',GKWAKE, 'SUGWD: GKDRAG=',GKDRAG

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUGWD',1,ZHOOK_HANDLE)
END SUBROUTINE SUGWD
