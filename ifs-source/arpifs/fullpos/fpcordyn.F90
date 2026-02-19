! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPCORDYN(YDQTYPE,YDNAMFPSCI,YDTFP,KFIELDS,KGPST,KGPEND,PROW,KFPROMA,CDCONF)

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMAFN   , ONLY : ALL_FULLPOS_TYPES
USE YOMFPC   , ONLY : TNAMFPSCI
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN


!**** *FPCORDYN*  - FULL-POS corrector for dynamic fields after 
!                   horizontal interpolations

!     PURPOSE.
!     --------
!        To correct dynamic fields after the horizontal interpolations that 
!        may create overshoots, or to recover the requested fields.

!**   INTERFACE.
!     ----------
!       *CALL* *FPCORDYN*

!        EXPLICIT ARGUMENTS
!        --------------------

!        INPUT:
!         KFIELDS  : number of fields in row
!         KGPST    : first output point in row.
!         KGPEND   : last output point in row.
!         KFPROMA  : length of the output row.
!         CDCONF   : configuration of work
!           'A' = gridpoint fields
!           'F' = spectrally fitted gridpoint fields
!           'B' = all dynamical fields (default)
!        INPUT/OUTPUT:
!         PROW     : fields rows to correct.

!        IMPLICIT ARGUMENTS
!        --------------------
!          See #include below. 

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!      Calls FPHOR12.

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        See documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 98-09-02 from FPINTDYN

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 01-08-07 Pruning options
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib : 04-09-29 Fix on o3
!      R. El Khatib 28-Jul-2016 Recode LWIDER_DOM
!      Y. Bouteloup: 11-Jan-2017 Add deep convective variables
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(IN) :: YDQTYPE
TYPE (TNAMFPSCI),  INTENT(IN) :: YDNAMFPSCI
TYPE(ALL_FULLPOS_TYPES), INTENT(IN) :: YDTFP
INTEGER(KIND=JPIM),      INTENT(IN) :: KFIELDS 
INTEGER(KIND=JPIM),      INTENT(IN) :: KFPROMA 
INTEGER(KIND=JPIM),      INTENT(IN) :: KGPST 
INTEGER(KIND=JPIM),      INTENT(IN) :: KGPEND 
REAL(KIND=JPRB),         INTENT(INOUT) :: PROW(KFPROMA,KFIELDS) 
CHARACTER(LEN=1),        INTENT(IN), OPTIONAL :: CDCONF 

REAL(KIND=JPRB) :: ZDUM(KFPROMA), ZMIN(KFPROMA), ZMAX(KFPROMA)
INTEGER(KIND=JPIM) :: ICOR, JFLD, JI, IJ, IPTR, INC, JL, IFLD
CHARACTER(LEN=1) :: CLCONF 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
#include "fphor12.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPCORDYN',0,ZHOOK_HANDLE)
ASSOCIATE(LFPLOSP=>YDNAMFPSCI%LFPLOSP, LFPRH100=>YDNAMFPSCI%LFPRH100)

IF (PRESENT(CDCONF)) THEN
  CLCONF=CDCONF
ELSE
  CLCONF='B'
ENDIF

IPTR=1
DO JFLD=1,YDQTYPE%NFPOSDYN
  IJ=YDQTYPE%NFPTRDYN(JFLD)
  IF ((YDQTYPE%ISF(IJ) == 0 .AND. CLCONF=='A').OR.(YDQTYPE%ISF(IJ) >= 1 .AND. CLCONF=='F').OR.CLCONF=='B') THEN

    INC=YDQTYPE%ILEV(IJ)

    IF (IJ==YDTFP%HU%ICOD) THEN
      ZMIN(KGPST:KGPEND)=0.0_JPRB
      IF (LFPRH100) THEN
!     0. < Relative moisture < 100.
      ZMAX(KGPST:KGPEND)=100.0_JPRB
    ELSE
!     0. < Relative moisture < 1.
      ZMAX(KGPST:KGPEND)=1.0_JPRB
    ENDIF
    ICOR=6
  ELSEIF(IJ==YDTFP%O3MX%ICOD) THEN
!   Ozone > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%CPF%ICOD .OR. IJ==YDTFP%SPF%ICOD) THEN
!   precipitation fluxes > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%L%ICOD) THEN
!   Cloud water > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%I%ICOD) THEN
!   Ice > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%CLF%ICOD) THEN
!   0. < Cloud fraction < 1.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ZMAX(KGPST:KGPEND)=1.0_JPRB
    ICOR=6
  ELSEIF(IJ==YDTFP%SN%ICOD) THEN
!   Snow > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%RR%ICOD) THEN
!   Rain > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%GR%ICOD) THEN
!   Graupel > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%HL%ICOD) THEN
!   Haill > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%LCONV%ICOD) THEN
!   Lconv > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%ICONV%ICOD) THEN
!   Iconv > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%RCONV%ICOD) THEN
!   Rconv > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%SCONV%ICOD) THEN
!   Sconv > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%UAL%ICOD) THEN
!   Ual > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%DAL%ICOD) THEN
!   Dal > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%UOM%ICOD) THEN
!   Uom > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%DOM%ICOD) THEN
!   Dom > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%TKE%ICOD) THEN
!   TKE > 1E-06
    ZMIN(KGPST:KGPEND)=0.000001_JPRB
    ICOR=4
  ELSEIF(IJ==YDTFP%EDR%ICOD) THEN
!     EDR > 1E-06
      ZMIN(KGPST:KGPEND)=0.000001_JPRB
      ICOR=4
  ELSEIF(IJ==YDTFP%SP%ICOD) THEN
!   After interpolation of log(Ps) for better consistency with PHIs,
!   Recover Ps or keep Log(Ps) :
    IF (.NOT.LFPLOSP) THEN
      DO JL=IPTR,IPTR+INC-1
        DO JI=KGPST,KGPEND
          PROW(JI,JL)=EXP(PROW(JI,JL))
        ENDDO
      ENDDO
    ENDIF
    ICOR=0
  ELSEIF(IJ==YDTFP%PJET%ICOD) THEN
!   Jet pressure > 0.
    ZMIN(KGPST:KGPEND)=0.0_JPRB
    ICOR=4
  ELSE
    ICOR=0
  ENDIF
  IF (ICOR /= 0) THEN
    DO JL=IPTR,IPTR+INC-1
      IFLD=JL
      CALL FPHOR12(IFLD,KGPST,KGPEND,PROW,KFIELDS,KFPROMA,ICOR,ZDUM,ZDUM, &
       & ZMIN,ZMAX,ZDUM,ZDUM)  
    ENDDO
  ENDIF
  IPTR=IPTR+INC

  ENDIF
ENDDO
IF (IPTR-1 /= KFIELDS) CALL ABOR1('FPCORDYN : INTERNAL ERROR')

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FPCORDYN',1,ZHOOK_HANDLE)
END SUBROUTINE FPCORDYN

