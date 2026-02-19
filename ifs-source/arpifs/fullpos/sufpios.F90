! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPIOS(KFPGRIB,KFPSURFEX,CDFPDIR,CDFPDOM,CDFPFN,CDFPCLIFNAME,CDFPSFXFNAME,YDNAMFPIOS)

!**** *SUFPIOS* - SET UP FULLPOS I/O SCHEME

!     Purpose.   TO SET UP COMMMON BLOCK YOMFPIOS WHICH CONTAINS   -
!     --------   PARAMETERS FOR USING THE MIO PACKAGE ON WORK
!                FILES. OPENS WORK FILES.
!                SET UP CONTROL ARRAYS AND LENGHTS OF BUFFERS
!                IF NO WORKFILES

!**   Interface.
!     ----------
!        *CALL* *SUFPIOS*

!        Explicit arguments :
!        --------------------
!           KFPGRIB      : level of GRIB encoding
!           KFPSURFEX    : Surfex usage for interoperability ISBA => Surfex 
!           CDFPDIR      : path or prefix for the output files
!           CDFPDOM      : array of names of the output domains
!           CDFPFN       : array of partial output filenames (filenames without extensions)
!           CDFPCLIFNAME : array of filename of climatology file on target geometry
!           CDFPSFXFNAME : array of filename of surfex climatology file on target geometries

!        Implicit arguments :
!        --------------------
!         See modules above.

!     Method.
!     -------

!     Externals.
!     ----------
!       SUFPSC2B.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 94-04-08

!     Modifications.
!     --------------
!      R. El Khatib : 01-08-07 Pruning options
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib : 09-Dec-2015 NFPWRITE
!      R. El Khatib : 09-Dec-2015 NFPADDING
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE PARFPOS  , ONLY : JPOSDOM
USE YOMLUN   , ONLY : NULOUT   ,NULNAM
USE YOMCT0   , ONLY : CNMEXP, LARPEGEF
USE YOMOPH0  , ONLY : CFNCLIMOUT, CFPEXTSFX
USE YOMFPIOS , ONLY : TNAMFPIOS

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KFPGRIB
INTEGER(KIND=JPIM), INTENT(IN) :: KFPSURFEX
CHARACTER(LEN=*), INTENT(IN) :: CDFPDIR
CHARACTER(LEN=*), INTENT(IN) :: CDFPDOM(:)
CHARACTER(LEN=*), INTENT(OUT) :: CDFPFN(:)
CHARACTER(LEN=*), INTENT(OUT) :: CDFPCLIFNAME(:)
CHARACTER(LEN=*), INTENT(OUT) :: CDFPSFXFNAME(:)
TYPE(TNAMFPIOS), TARGET, INTENT(OUT) :: YDNAMFPIOS

! CDFPFN      : partial output filenames (filenames without extensions)
! CFPCLIFNAME : filename of climatology file on target geometry
! CFPSFXFNAME : filename of surfex climatology file on target geometries
CHARACTER(LEN=180) :: CFPFN      (JPOSDOM)
CHARACTER(LEN=180) :: CFPCLIFNAME(JPOSDOM)
CHARACTER(LEN=180) :: CFPSFXFNAME(JPOSDOM)

INTEGER(KIND=JPIM), POINTER :: NFPWRITE, NFPDIGITS, NFPXFLD

! namphmse should not be read : what we need are variables specific to fullpos I/Os.
! For now namphmse is read for continuity with the older namelists
LOGICAL, POINTER :: LFTZERO, LPGDFWR, LHISFWR
REAL (KIND=JPRB), POINTER ::  XZSEPS
INTEGER(KIND=JPIM), POINTER :: NSURFEXCTL


INTEGER(KIND=JPIM) :: J, ISTAT, ILASTCHAR, IFPDOM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "posnam.intfb.h"
#include "posname.intfb.h"

#include "namfpios.nam.h"
#include "namphmse.nam.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUFPIOS',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

WRITE(NULOUT,'('' == Full-Pos : setup I/O handling == '')')

!*       0. POINTERS TO THE NAMELIST STRUCTURE
!           ----------------------------------

NFPXFLD=>YDNAMFPIOS%NFPXFLD
NFPWRITE=>YDNAMFPIOS%NFPWRITE
NFPDIGITS=>YDNAMFPIOS%NFPDIGITS
LPGDFWR=>YDNAMFPIOS%LFPPGDFWR
LHISFWR=>YDNAMFPIOS%LFPHISFWR

! already read by a namelist
YDNAMFPIOS%NFPGRIB=KFPGRIB

ILASTCHAR=LEN_TRIM(CFNCLIMOUT)

!*       1.    I/O PARAMETERS UNDER USER CONTROL.
!              ----------------------------------

NFPXFLD=-999
NFPWRITE=1
IF (LARPEGEF) THEN
  NFPDIGITS=4
ELSE
  NFPDIGITS=6
ENDIF

IFPDOM=SIZE(CDFPDOM)
DO J=1, IFPDOM
  CFPFN(J)=TRIM(CDFPDIR)//CNMEXP(1:4)//CDFPDOM(J)
  ! traditional clim filename on target geometry
  CFPCLIFNAME(J)=TRIM(CFNCLIMOUT)//TRIM(CDFPDOM(J))
  ! Surfex clim filename on target geometry
  CFPSFXFNAME(J)=CFNCLIMOUT(1:ILASTCHAR-1)//TRIM(CFPEXTSFX)//CFNCLIMOUT(ILASTCHAR:ILASTCHAR)//TRIM(CDFPDOM(J))
ENDDO


!*       2.   READ NAMELIST
!             -------------

CALL POSNAM(NULNAM,'NAMFPIOS')
READ(NULNAM,NAMFPIOS)

IF (KFPSURFEX == 1) THEN
  LPGDFWR=.FALSE.
  LHISFWR=.TRUE.
  ! posname for transition because I wish we could move later these variables to yomfpios
  CALL POSNAME(NULNAM,'NAMPHMSE',ISTAT)
  IF (ISTAT == 0) THEN
    READ (NULNAM, NAMPHMSE)
  ENDIF
ENDIF

!*       3. SAVE FILENAMES
!           --------------

DO J=1,IFPDOM
  CDFPFN(J)=TRIM(CFPFN(J))
  CDFPCLIFNAME(J)=TRIM(CFPCLIFNAME(J))
  CDFPSFXFNAME(J)=TRIM(CFPSFXFNAME(J))
ENDDO

IF (LHOOK) CALL DR_HOOK('SUFPIOS',1,ZHOOK_HANDLE)

END SUBROUTINE SUFPIOS
