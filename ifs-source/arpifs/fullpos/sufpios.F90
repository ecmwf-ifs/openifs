! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPIOS(KFPGRIB,KFPSURFEX,PTSTEP,KSTOP,CDFPDIR,CDFPDOM,CDFPFN,CDFPCLIFNAME,CDFPSFXFNAME,CDMODEL,LDEXTERN,YDNAMFPIOS)

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
!           PTSTEP       : forecasting model time step
!           KSTOP        : forecasting model number of time steps
!           CDFPDIR      : path or prefix for the output files
!           CDFPDOM      : array of names of the output domains
!           CDFPFN       : array of partial output filenames (filenames without extensions)
!           CDFPCLIFNAME : array of filename of climatology file on target geometry
!           CDFPSFXFNAME : array of filename of surfex climatology file on target geometries
!           CDMODEL      : FA model name
!           LDEXTERN     : .TRUE. to write fields in a separate GRIB2 file

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
!      R. El Khatib : 23-Jun-2021 NFRFPDI, NFPDITS
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE PARFPOS  , ONLY : JPOSDOM
USE YOMLUN   , ONLY : NULOUT   ,NULNAM
USE YOMCT0   , ONLY : CNMEXP, LARPEGEF, NCONF, LECMWF
USE YOMFA    , ONLY : CMODEL, LEXTERN
USE YOMOPH0  , ONLY : CFNCLIMOUT, CFPEXTSFX
USE YOMFPIOS , ONLY : TNAMFPIOS

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KFPGRIB
INTEGER(KIND=JPIM), INTENT(IN) :: KFPSURFEX
REAL(KIND=JPRB),    INTENT(IN) :: PTSTEP
INTEGER(KIND=JPIM), INTENT(IN) :: KSTOP
CHARACTER(LEN=*), INTENT(IN) :: CDFPDIR
CHARACTER(LEN=*), INTENT(IN) :: CDFPDOM(:)
CHARACTER(LEN=*), INTENT(OUT) :: CDFPFN(:)
CHARACTER(LEN=*), INTENT(OUT) :: CDFPCLIFNAME(:)
CHARACTER(LEN=*), INTENT(OUT) :: CDFPSFXFNAME(:)
CHARACTER(LEN=*), INTENT(OUT) :: CDMODEL(:)
LOGICAL,          INTENT(OUT) :: LDEXTERN(:)
TYPE(TNAMFPIOS), TARGET, INTENT(OUT) :: YDNAMFPIOS

! CDFPFN, CFPFN : partial output filenames (filenames without extensions)
! CFPCLIFNAME : filename of climatology file on target geometry
! CFPSFXFNAME : filename of surfex climatology file on target geometries
CHARACTER(LEN=180) :: CFPFN      (JPOSDOM) ! sorry Doctor, the variable are in namelist, I won't change it.
CHARACTER(LEN=180) :: CFPCLIFNAME(JPOSDOM)
CHARACTER(LEN=180) :: CFPSFXFNAME(JPOSDOM)
CHARACTER(LEN=64)  :: CFAMODEL(JPOSDOM)
LOGICAL            :: LFAEXTERN(JPOSDOM)

INTEGER(KIND=JPIM), POINTER :: NFPWRITE, NFPDIGITS, NFPXFLD

! namphmse should not be read : what we need are variables specific to fullpos I/Os.
! For now namphmse is read for continuity with the older namelists
LOGICAL, POINTER :: LFTZERO, LPGDFWR, LHISFWR
REAL (KIND=JPRB), POINTER ::  XZSEPS
INTEGER(KIND=JPIM), POINTER :: NSURFEXCTL
INTEGER(KIND=JPIM), POINTER :: NFPDITS(:), NFPDITSMIN(:), NFRFPDI


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
NFPDITS=>YDNAMFPIOS%NFPDITS
NFPDITSMIN=>YDNAMFPIOS%NFPDITSMIN
NFRFPDI=>YDNAMFPIOS%NFRFPDI
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
  CFAMODEL(J)=CMODEL
  LFAEXTERN(J)=LEXTERN
ENDDO

IF (NCONF==1) THEN
  NFPDITS(:)=0
  IF (LECMWF) THEN
    NFRFPDI=1
  ELSE
    ! reduce norms frequency to 6 hours for computational time savings (unless the forecast is less, or no forecast)
    NFRFPDI=MIN(MAX(1,KSTOP),NINT(6._JPRB*3600._JPRB/PTSTEP))
  ENDIF
ENDIF

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

!*       3. SAVE FILENAMES & GRIB-SPECIFIC HANDLING PARAMETERS
!           --------------------------------------------------

DO J=1,IFPDOM
  CDFPFN(J)=TRIM(CFPFN(J))
  CDFPCLIFNAME(J)=TRIM(CFPCLIFNAME(J))
  CDFPSFXFNAME(J)=TRIM(CFPSFXFNAME(J))
  CDMODEL(J)=TRIM(CFAMODEL(J))
  LDEXTERN(J)=LFAEXTERN(J)
ENDDO

IF (LHOOK) CALL DR_HOOK('SUFPIOS',1,ZHOOK_HANDLE)

END SUBROUTINE SUFPIOS
