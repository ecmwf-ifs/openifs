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

SUBROUTINE SUCSLINT(LDSLHD,LDCOMAD,LDVSPLIP,LDHV,LDQM,LDQMH,LDQM3D,LDQML3D,&
  & LDINTLIN,LDVWENO,CDSLINT)

!**** *SUCSLINT*   - Set kind of semi-Lagrangian interpolation
!                    of GFL fields.

!     Purpose.
!     --------
!           Set the kind of SL interpolation for a given
!           GFL field.

!**   Interface.
!     ----------
!        *CALL* *SUCSLINT(...)

!        Explicit arguments :
!        --------------------
!       Input:
!        LDSLHD    -  key to activate SLHD. 
!        LDCOMAD   -  key to activate COMAD. 
!        LDVSPLIP  -  key to use vertical monotonic spline for the
!                     high order SL interpolation.
!        LDHV      -  key to use vertical hermite cubic interpolation
!                     for the high order SL interpolation.
!        LDQM      -  key to use quasi-monotonic interpolator for
!                     the high order SL interpolation.
!        LDQMH     -  key to use horizontally quasi-monotonic interpolator
!                     for the high order SL interpolation.
!        LDQM3D    -  key to use quasi-monotonic interpolator directly applied
!                     in 3-d for the high order SL interpolation.
!        LDINTLIN  -  key to use linear interpolation
!        LDVWENO   -  key to use vertical quintic interpolation of WENO family

!       Output:
!        CDSLINT   -  the name of the resulting SL interpolation.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Called by SUDYN

!     Reference.
!     ----------
!        Documentation about SL advection

!     Author.
!     -------
!        Filip VANA and Karim YESSAD
!        Meteo-France

! Modifications
! -------------
!   Original : NOVEMBER 2004
!   28-Aug-2007 F. Vana: removing LRSPLINE attribute
!   30-Jun-2008 J. Masek   New SLHD interpolators.
!   30-Mar-2011 S. Malardel  add LDINTLIN
!   S. Malardel (Nov 2013): COMAD weights for SL interpolations
!   Feb 2014    M. Diamantakis add LAITQM3D
!   Feb 2016    M. Diamantakis add LAITQML3D
!   20-Feb-2019 F. Vana: Vertical quintic interpolation of WENO type
!   Nov-2019 F. Vana:  new interpolation type LAITRWENOQM3
!------------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL,INTENT(IN)            :: LDSLHD
LOGICAL,INTENT(IN)            :: LDCOMAD
LOGICAL,INTENT(IN)            :: LDVSPLIP
LOGICAL,INTENT(IN)            :: LDHV
LOGICAL,INTENT(IN)            :: LDQM
LOGICAL,INTENT(IN)            :: LDQMH
LOGICAL,INTENT(IN)            :: LDQM3D
LOGICAL,INTENT(IN)            :: LDQML3D
LOGICAL,INTENT(IN)            :: LDINTLIN
LOGICAL,INTENT(IN)            :: LDVWENO
CHARACTER(LEN=12),INTENT(OUT) :: CDSLINT

!     ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUCSLINT',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.   SET GFL KIND OF SL INTERPOLATION
!             --------------------------------

IF(LDSLHD) THEN

!* SLHD interpolations

  IF(LDQM) THEN
    CDSLINT='LAITQM(SLD) '
  ELSEIF(LDQMH) THEN
    CDSLINT='LAITQMH(SLD)'
  ELSE
    CDSLINT='LAITRI(SLD) '
  ENDIF

ELSE

!*  non SLHD interpolations

  IF(LDVSPLIP) THEN
      CDSLINT='LAITVSPCQM  '
  ELSEIF(LDINTLIN) THEN
    IF(LDCOMAD) THEN
      CDSLINT='LIN(MAD)    '
    ELSE
      CDSLINT='LIN         '
    ENDIF
  ELSEIF(LDVWENO) THEN
    IF(LDQM) THEN
      CDSLINT='LAITRWENOQM '
    ELSEIF(LDQM3D) THEN
      CDSLINT='LAITRWENOQM3'
    ELSEIF(LDQMH) THEN
      CDSLINT='LAITRWENOQMH'
    ELSE
      CDSLINT='LAITRWENO   '
    ENDIF
  ELSEIF(LDHV) THEN
    IF(LDQM) THEN
      CDSLINT='LAIHVTQM    '
    ELSEIF(LDQMH) THEN
      CDSLINT='LAIHVTQMH   '
    ELSE
      CDSLINT='LAIHVT      '
    ENDIF
  ELSE
    IF(LDQM3D) THEN
      IF(LDCOMAD) THEN
        CDSLINT='LAITQM3(MAD)'
      ELSE
        CDSLINT='LAITQM3D    '
      ENDIF
    ELSEIF(LDQML3D) THEN
      IF(LDCOMAD) THEN
        CDSLINT='LAITQML(MAD)'
      ELSE
        CDSLINT='LAITQML3D   '
      ENDIF
    ELSEIF(LDQM) THEN
      IF(LDCOMAD) THEN
        CDSLINT='LAITQM(MAD) '
      ELSE
        CDSLINT='LAITQM      '
      ENDIF
    ELSEIF(LDQMH) THEN
      IF(LDCOMAD) THEN
        CDSLINT='LAITQMH(MAD)'
      ELSE
        CDSLINT='LAITQMH     '
      ENDIF
    ELSE
      IF(LDCOMAD) THEN
        CDSLINT='LAITRI(MAD) '
      ELSE
        CDSLINT='LAITRI      '
      ENDIF
    ENDIF
  ENDIF

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUCSLINT',1,ZHOOK_HANDLE)
END SUBROUTINE SUCSLINT
