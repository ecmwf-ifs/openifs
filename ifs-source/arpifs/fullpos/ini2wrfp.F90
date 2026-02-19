! (C) Copyright 1989- Meteo-France.

SUBROUTINE INI2WRFP(YDRQSP,YDRQGP,YDRQPHY,CDCONF,LDPADDING,PUNDEF,LDYN,LDCOSP,KFIELDS,KFPDOM,KVGRIB,YDTFP_DYNDS,YDGFP_PHYDS,YDFLDSC)

!**** *INI2WRFP*  - PREPARE TO WRITE OUT THE HORIZONTALLY POST-PROCESSED 
!                   FIELDS TO ARPEGE/ALADIN FILES -FIELD-DEPENDANT VARIABLES

!     PURPOSE.
!     --------
!           To compute variables which depend of the fields themselves, 
!           prior to writing out to post-processed files

!**   INTERFACE.
!     ----------
!       *CALL* *INI2WRFP(...)

!        EXPLICIT ARGUMENTS
!        --------------------
!        CDCONF  : configuration of work
!        LDPADDING : padding domains wider than the model one
!        PUNDEF  : value for gridpoints outside the input domain
!        LDYN    : .TRUE. if dynamic fields ; else physical fields
!        LDCOSP  : .TRUE. if spectral fields ; else gridpoint fields
!        KFIELDS : number of fields
!        KFPDOM  : number of subdomains
!        KVGRIB  : level of GRIB encoding
!        YDTFP_DYNDS : overall fullpos dynamic fields descriptors
!        YDGFP_PHYDS : overall fullpos physical fields descriptors
!        YDFLDSC : Field descriptors

!        IMPLICIT ARGUMENTS
!        ------------------
!         See modules above

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 98-12-17 from wrhfpsm/dm and wrsfpsm/dm

!     MODIFICATIONS
!     --------------
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      R. El Khatib : 03-05-27 Enable Pressure level below 1000 hPa despite *FA*
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib 17-Jul-2013 FABEC post-processing
!      P.Marguinaud : 10-10-2013 : Add extra arguments
!      P.Marguinaud : 01-10-2013 : Make LDMASK optional
!      R. El Khatib 20-Nov-2015 : simplify configuration of work
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE FULLPOS_MIX, ONLY : FULLPOS_TYPE
USE TYPE_FPDSPHYS, ONLY : FPDSPHY
USE YOMFP4L      , ONLY : TRQFP
USE IOFLDDESC_MOD, ONLY : IOFLDDESC

IMPLICIT NONE

TYPE (TRQFP),  INTENT(IN) :: YDRQSP
TYPE (TRQFP),  INTENT(IN) :: YDRQGP
TYPE (TRQFP),  INTENT(IN) :: YDRQPHY
CHARACTER(LEN=1)  ,INTENT(IN)    :: CDCONF
LOGICAL           ,INTENT(IN)    :: LDPADDING
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUNDEF
LOGICAL           ,INTENT(IN)    :: LDYN
LOGICAL           ,INTENT(IN)    :: LDCOSP
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPDOM
INTEGER(KIND=JPIM),INTENT(IN)    :: KVGRIB
TYPE(FULLPOS_TYPE),INTENT(IN)    :: YDTFP_DYNDS(:)
TYPE(FPDSPHY),     INTENT(IN)    :: YDGFP_PHYDS(:)
TYPE (IOFLDDESC)  ,INTENT(INOUT) :: YDFLDSC(KFIELDS)

!     ZUNIT   : unit for level value that will be written out on files

!     IFLDA  : absolute field pointer 
!     JFL    : local field pointer

!     YDRQSP%ICOD : field pointers for concatenated spectral arrays
!     YDRQSP%ZLEV : level pointers for concatenated spectral arrays

REAL(KIND=JPRB) :: ZUNIT

INTEGER(KIND=JPIM) :: ID, IFLDA, JFL, JD, IIL
INTEGER(KIND=JPIM) :: IMASK(KFIELDS,KFPDOM)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
!     ------------------------------------------------------------------

!*       1. PREPARATIONS
!           ------------

IF (LHOOK) CALL DR_HOOK('INI2WRFP',0,ZHOOK_HANDLE)

YDFLDSC(:)%LSPEC = LDCOSP

YDFLDSC(:)%NGRIBL = KVGRIB
YDFLDSC(:)%ICPLB  = -1
YDFLDSC(:)%ICPLS  = -1

IF (LDYN) THEN
  YDFLDSC(:)%CPREF=CDCONF
  IF (CDCONF == 'V') THEN
!   Potential vorticity levels are written in deci-PVU :
    ZUNIT=10000000._JPRB
  ELSE
    ZUNIT=1.0_JPRB
  ENDIF
ENDIF

!*       2. SETUP PARAMETERS FOR PACKING AND FIELD NAME CONSTRUCTION
!           --------------------------------------------------------

IMASK(:,:)=0
DO JFL=1,KFIELDS
  IF (CDCONF=='I') THEN
!   Not ready for upper air physical post-processing
    DO JD=1,YDRQPHY%IDOM(JFL)
      ID=YDRQPHY%IDMP(JD,JFL)
      IMASK(JFL,ID)=1
    ENDDO
    IFLDA=YDRQPHY%ICOD(JFL)
    YDFLDSC(JFL)%JBITS=YDGFP_PHYDS(IFLDA)%IBITS
    YDFLDSC(JFL)%ILEVG=0
    YDFLDSC(JFL)%CPREF=YDGFP_PHYDS(IFLDA)%CLNAME(1:4)
    YDFLDSC(JFL)%CSUFF=YDGFP_PHYDS(IFLDA)%CLNAME(5:)
  ELSEIF(LDCOSP) THEN
    IFLDA=YDRQSP%ICOD(JFL)
    YDFLDSC(JFL)%JBITS=YDTFP_DYNDS(IFLDA)%IBITS
    IMASK(:,:)=1
    IF (YDRQSP%LLSURF(JFL)) THEN
      YDFLDSC(JFL)%CPREF=YDTFP_DYNDS(IFLDA)%CLNAME(1:4)
      YDFLDSC(JFL)%ILEVG=0
      YDFLDSC(JFL)%CSUFF=YDTFP_DYNDS(IFLDA)%CLNAME(5:16)
    ELSE
      YDFLDSC(JFL)%ILEVG=NINT(YDRQSP%ZLEV(JFL)*ZUNIT)
      YDFLDSC(JFL)%CSUFF=YDTFP_DYNDS(IFLDA)%CLNAME(1:12)
      SELECT CASE (CDCONF)
      CASE ('P')
!         This enables pressure level below 1000 hPa in spite of *FA* limit.
!         level value is then written on 4 digits instead of 3
        IF (YDFLDSC(JFL)%ILEVG > 100000) THEN
          YDFLDSC(JFL)%CPREF(2:)=' '
          WRITE(YDFLDSC(JFL)%CPREF(2:7),'(I6)') YDFLDSC(JFL)%ILEVG
        ENDIF
      CASE ('K')
!         Iso- surface can be compute either from the bottom (B) or from the top (T)
        YDFLDSC(JFL)%CPREF(2:)=' '
        YDFLDSC(JFL)%CPREF(2:2)='B'
        IF (YDFLDSC(JFL)%ILEVG < 0) THEN
          YDFLDSC(JFL)%CPREF(2:2)='T'
          YDFLDSC(JFL)%ILEVG=-NINT(YDRQSP%ZLEV(JFL)*ZUNIT)
        ENDIF
      END SELECT
    ENDIF
  ELSE
    DO JD=1,YDRQGP%IDOM(JFL)
      ID=YDRQGP%IDMP(JD,JFL)
      IMASK(JFL,ID)=1
    ENDDO
    IFLDA=YDRQGP%ICOD(JFL)
    YDFLDSC(JFL)%JBITS=YDTFP_DYNDS(IFLDA)%IBITS
    IF (YDRQGP%LLSURF(JFL)) THEN
      YDFLDSC(JFL)%CPREF=YDTFP_DYNDS(IFLDA)%CLNAME(1:4)
      YDFLDSC(JFL)%ILEVG=0
      YDFLDSC(JFL)%CSUFF=YDTFP_DYNDS(IFLDA)%CLNAME(5:16)
    ELSE
      YDFLDSC(JFL)%ILEVG=NINT(YDRQGP%ZLEV(JFL)*ZUNIT)
      YDFLDSC(JFL)%CSUFF=YDTFP_DYNDS(IFLDA)%CLNAME(1:12)
      SELECT CASE (CDCONF)
      CASE ('P')
!         This enables pressure level below 1000 hPa in spite of *FA* limit.
!         level value is then written on 4 digits instead of 3
        IF (YDFLDSC(JFL)%ILEVG > 100000) THEN
          YDFLDSC(JFL)%CPREF(2:)=' '
          WRITE(YDFLDSC(JFL)%CPREF(2:7),'(I6)') YDFLDSC(JFL)%ILEVG
        ENDIF
      CASE ('K')
!         Iso- surface can be computed either from the bottom (B) or from the top (T)
        YDFLDSC(JFL)%CPREF(2:)=' '
        YDFLDSC(JFL)%CPREF(2:2)='B'
        IF (YDFLDSC(JFL)%ILEVG < 0) THEN
          YDFLDSC(JFL)%CPREF(2:2)='T'
          YDFLDSC(JFL)%ILEVG=-NINT(YDRQGP%ZLEV(JFL)*ZUNIT)
        ENDIF
      END SELECT
    ENDIF
! Hollow fields
    DO IIL = 1, SIZE (YDTFP_DYNDS(IFLDA)%ILEVHOLI)
      IF ((YDFLDSC (JFL)%ILEVG >= YDTFP_DYNDS(IFLDA)%ILEVHOLI (IIL)) &
  & .AND. (YDFLDSC (JFL)%ILEVG <= YDTFP_DYNDS(IFLDA)%ILEVHOLF (IIL))) THEN
        YDFLDSC (JFL)%NGRIBL = 4_JPIM
        YDFLDSC (JFL)%ICPLS = YDTFP_DYNDS(IFLDA)%ICPLSIZE
        YDFLDSC (JFL)%ICPLB = YDTFP_DYNDS(IFLDA)%ICPLBITS
      ENDIF
    ENDDO
  ENDIF
ENDDO

WHERE (YDFLDSC%ICPLB < 0) 
  YDFLDSC%ICPLB = YDFLDSC%JBITS
ENDWHERE

! Fullpos mask
IF (KFPDOM > KIND(YDFLDSC(:)%IFPMASK)*8) CALL ABOR1('INI2WRFP : KFPDOM TOO BIG FOR KIND(IFPMASK) !')
YDFLDSC%IFPMASK = 0
DO JFL = 1, KFIELDS
  DO JD = KFPDOM, 1, -1
    IF (IMASK(JFL,JD) /= 0) THEN
      YDFLDSC(JFL)%IFPMASK = IBSET (YDFLDSC (JFL)%IFPMASK, JD-1)
    ENDIF
  ENDDO
ENDDO

IF (LDPADDING) THEN
  YDFLDSC (:)%XUNDF = PUNDEF
  YDFLDSC (:)%LUNDF = .TRUE.
ENDIF

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('INI2WRFP',1,ZHOOK_HANDLE)
END SUBROUTINE INI2WRFP
