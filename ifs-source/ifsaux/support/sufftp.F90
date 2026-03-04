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

SUBROUTINE SUFFTP(KDLON,KFFTP0,LDODD)

!**** *SUFFTP*  - Initialize possible numbers for FFT's

!     Purpose.
!     --------
!           Initialize possible numbers for FFT's

!**   Interface.
!     ----------
!        *CALL* *SUFFTP(KDLON,KFFTP0)

!        Explicit arguments :
!        --------------------
!        LDODD : .TRUE. if odd numbers are allowed
!        KFFTP0(KDLON) 0 ---> possible for FFT
!                      1 ---> impossible for FFT

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    none.
!     ----------

!     Reference.
!     ----------
!        ARPEGE documentation

!     Author.
!     -------
!        Philippe Courtier  *DMN*

!     Modifications.
!     --------------
!        Creation of cycle 4
!        R. El Khatib 15-May-2013 Option to allow or not odd numbers
!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KDLON
INTEGER(KIND=JPIM), INTENT(OUT) :: KFFTP0(KDLON)
LOGICAL, INTENT(IN) :: LDODD

INTEGER(KIND=JPIM) :: I2MAX, I3MAX, I5MAX, IMIN, IN, J, J2, J3, J5


REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('REDUCED_GRID',0,ZHOOK_HANDLE)

!*       1.    Initialize up to IMIN-1
!              --------------------

IMIN=16
IF(KDLON <= IMIN-1)THEN
  CALL ABOR1(' ERROR IN SUFFTP0 ')
ENDIF

DO J=1,IMIN-1
  KFFTP0(J)=1
ENDDO

!     ------------------------------------------------------------------

!*       2.    Initialize from 16.
!              -------------------

I2MAX=INT(LOG(REAL(KDLON,JPRB))/LOG(2._JPRB))+1
I3MAX=INT(LOG(REAL(KDLON,JPRB))/LOG(3._JPRB))+1
I5MAX=INT(LOG(REAL(KDLON,JPRB))/LOG(5._JPRB))+1

DO J=IMIN,KDLON

  IF (.NOT.LDODD) THEN
    IF(MOD(J,2) /= 0)THEN
      KFFTP0(J)=1
      GOTO 202
    ENDIF
  ENDIF

  IN=J
  DO J2=1,I2MAX
    IF(MOD(IN,2) == 0)THEN
      IN=IN/2
    ENDIF
  ENDDO

  DO J3=1,I3MAX
    IF(MOD(IN,3) == 0)THEN
      IN=IN/3
    ENDIF
  ENDDO

  DO J5=1,I5MAX
    IF(MOD(IN,5) == 0)THEN
      IN=IN/5
    ENDIF
  ENDDO

  IF(IN /= 1)THEN
    KFFTP0(J)=1
  ELSE
    KFFTP0(J)=0
  ENDIF
  202 CONTINUE
ENDDO

! Memo :
! this is old stuff and is not used anymore, unfortunately some resolutions like T799 have been created like this !!
!  The possibility of 625 or 1250 points in the T426, TL639 grids
!   should be eliminated due to the radiation code (memo from Lars)      
!IF(KDLON >= 640) THEN
!  KFFTP0(625)=1
!  IF(KDLON >= 1250) THEN
!    KFFTP0(1250)=1
!  ENDIF
!ENDIF

IF (LHOOK) CALL DR_HOOK('REDUCED_GRID',1,ZHOOK_HANDLE)

END SUBROUTINE SUFFTP
