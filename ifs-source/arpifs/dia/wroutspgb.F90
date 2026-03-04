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

SUBROUTINE WROUTSPGB(YDDIM,PSPFLD,KSHMAX,KGRIBIO,CDLEVTYPE,KRESOL,CDFNSH,PTSTEP)

!**** *WROUTSPGB* _ GRIB codes and writes out spectral fields

!     Purpose.
!     --------
!     Write out spectral fields in GRIB

!**   Interface.
!     ----------
!        *CALL* *WROUTSPGB*(...)

!        Explicit arguments :    
!        --------------------

!        Implicit arguments :      The state variables of the model
!        --------------------

!     Method.
!     -------
!        See documentation
!        - spectral part would only work if input-resol = output resol
!          for vertical and horizontal resolution

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud  *ECMWF*
!        Original : 01-12-14 Adapted from WR****

!     Modifications.
!     --------------
!       P.Towers : 02-11-13 Fixes for fewer writers than gatherers
!       P.Towers : 03-04-23 Added ISETFIELDCOUNTFDB logic
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!        G. Radnoti:   25-Aug-2004 treatment of tstep<0 for LTWINC
!        M.Hamrud      01-Dec-2005 Generalized IO scheme
!       R. El Khatib : 23-Oct-2008 use suofname
!       R. El Khatib : 20-Aug-2012 optional argument KRESOL and consequences
!       G. Carver :    22-May-2013 fixed write mode for non-FDB I/O
!       K. Yessad (July 2014): Move some variables.
!       R. El Khatib : 17-Mar-2016 CFPATH passed to suofname
!     ------------------------------------------------------------------

USE YOMDIM   , ONLY : TDIM
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : NOUTTYPE
USE IOSTREAM_MIX , ONLY : SETUP_IOSTREAM, SETUP_IOREQUEST, IO_PUT,&
 & CLOSE_IOSTREAM, TYPE_IOSTREAM , TYPE_IOREQUEST, Y_IOSTREAM_FDB,&
 & CLOSE_IOREQUEST

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIM)        , INTENT(IN) :: YDDIM
INTEGER(KIND=JPIM), INTENT(IN) :: KSHMAX 
REAL(KIND=JPRB)   , INTENT(IN) :: PSPFLD(:,:) 
INTEGER(KIND=JPIM), INTENT(IN) :: KGRIBIO(4,KSHMAX) 
CHARACTER(LEN=1)  , INTENT(IN) :: CDLEVTYPE 
INTEGER(KIND=JPIM), INTENT(IN) :: KRESOL
CHARACTER(LEN=210), INTENT(IN) :: CDFNSH
REAL(KIND=JPRB)   , INTENT(IN), OPTIONAL :: PTSTEP

!     ------------------------------------------------------------------

CHARACTER :: CLEVT*2

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

TYPE(TYPE_IOSTREAM) :: YL_IOSTREAM
TYPE(TYPE_IOREQUEST) :: YL_IOREQUEST

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WROUTSPGB',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

IF(CDLEVTYPE == 'm') THEN
  CLEVT='ML'
ELSEIF(CDLEVTYPE == 'p') THEN
  CLEVT='PL'
ELSEIF(CDLEVTYPE == 'v') THEN
  CLEVT='PV'
ELSEIF(CDLEVTYPE == 't') THEN
  CLEVT='TH'
ELSE
  CALL ABOR1('WROUTSPGB:UNKNOWN LEVEL TYPE')
ENDIF

CALL SETUP_IOREQUEST(YL_IOREQUEST,'SPECTRAL_FIELDS',LDGRIB=.TRUE.,&
 & KGRIB2D=KGRIBIO(1,:),KLEVS2D=KGRIBIO(2,:),KBSET2D=KGRIBIO(3,:),&
 & CDLEVTYPE=CLEVT,KRESOL=KRESOL,PTSTEP=PTSTEP)

IF(NOUTTYPE == 2) THEN
  CALL IO_PUT(Y_IOSTREAM_FDB,YL_IOREQUEST,PR2=PSPFLD)
ELSE
  CALL SETUP_IOSTREAM(YL_IOSTREAM,'CIO',TRIM(CDFNSH),CDMODE='a',KIOMASTER=1)
  CALL IO_PUT(YL_IOSTREAM,YL_IOREQUEST,PR2=PSPFLD)
  CALL CLOSE_IOSTREAM(YL_IOSTREAM)
ENDIF
 
CALL CLOSE_IOREQUEST(YL_IOREQUEST)

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WROUTSPGB',1,ZHOOK_HANDLE)
END SUBROUTINE WROUTSPGB
