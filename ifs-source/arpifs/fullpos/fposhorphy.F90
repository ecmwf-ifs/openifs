! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPOSHORPHY(YDRQPHY,YDCLIMO,YDNAMFPINT,YDNAMFPSCI,YDAFN,YDFPGEOMETRY,YDFPSUW,YDFPWSTD,YDFPSTRUCT,YDFPGEO_DEP, &
 & KFLDIN,PBUF1,KCOD,PFP,YDRQAUX)

!**** *FPOSHORPHY*  - HORIZONTAL POST-PROCESSING

!     PURPOSE.
!     --------
!        PERFORM THE HORIZONTAL INTERPOLATIONS FOR PHYSICAL FIELDS

!        Computations are DM-local if distributed memory.

!**   INTERFACE.
!     ----------
!       *CALL* *FPOSHORPHY*

!        EXPLICIT ARGUMENTS
!        --------------------

!     All dummy arguments are input ones.

!     PBUF1  : buffer 1 which contain the fields to interpolate
!     PFP  : interpolated fields

!        IMPLICIT ARGUMENTS
!        ------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        See documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 01-08-07 Pruning options
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!      R. El Khatib  20-May-2005 NFPWIDE moved to YOMWFPDS
!      K. Yessad: 28-02-2007 DM-environment optimisations in FULL-POS.
!      G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      R. El Khatib  24-Jul-2012 LFPDISTRIB replaced by NFPDISTRIB
!      R. El Khatib 13-Dec-2012 Fullpos buffers reshaping
!      P.Marguinaud 01-Oct-2014 PREP fields interpolations
!      R. El Khatib 04-Aug-2016 protection for Boyd biperiodicization
!     ------------------------------------------------------------------

USE PARKIND1     , ONLY : JPIM     ,JPRB
USE YOMHOOK      , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMFPGEOMETRY, ONLY : TFPGEOMETRY
USE YOMFPGEO     , ONLY : TFPGEO
USE TYPE_FPOSBUF , ONLY : FPOSBUF
USE YOMFPC       , ONLY : TNAMFPINT, TNAMFPSCI
USE EINT_MOD     , ONLY : SL_STRUCT
USE YOMWFPB      , ONLY : TFPWSTD, TFPSUW
USE YOMAFN       , ONLY : TAFN
USE YOMFP4L, ONLY : TRQFP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TRQFP), INTENT(IN) :: YDRQPHY
TYPE (FPOSBUF)   , INTENT(IN) :: YDCLIMO
TYPE (TNAMFPINT),  INTENT(IN) :: YDNAMFPINT
TYPE (TNAMFPSCI),  INTENT(IN) :: YDNAMFPSCI
TYPE(TAFN),        INTENT(IN) :: YDAFN
TYPE (TFPGEOMETRY),INTENT(IN) :: YDFPGEOMETRY
TYPE (TFPSUW),     INTENT(IN) :: YDFPSUW
TYPE (TFPWSTD),    INTENT(IN) :: YDFPWSTD
TYPE(SL_STRUCT),   INTENT(IN) :: YDFPSTRUCT
TYPE (TFPGEO),     INTENT(IN) :: YDFPGEO_DEP
INTEGER(KIND=JPIM),INTENT(IN) :: KFLDIN
REAL(KIND=JPRB),   INTENT(IN) :: PBUF1(YDFPSTRUCT%NASLB1*KFLDIN) 
INTEGER(KIND=JPIM), INTENT(IN) :: KCOD(YDRQPHY%NFIELDG)
REAL(KIND=JPRB),   INTENT(OUT) :: PFP(YDFPGEO_DEP%NFPROMA,YDRQPHY%NFIELDG,YDFPGEO_DEP%NFPBLOCS)
TYPE (TRQFP), INTENT(IN), OPTIONAL :: YDRQAUX

!     ------------------------------------------------------------------

!     Handling horizontal interpolations of physical fields : 
!     LLCLI=.T. if the field is overwritten by the climatology
!     LLNIL=.T. if the field cannot be "directly" interpolated
!     LLINT=.T. if the field is to be interpolated
!     IHINTS : number of times a field has been computed

LOGICAL :: LLCLI(YDRQPHY%NFIELDG), LLNIL(YDRQPHY%NFIELDG), LLINT(YDRQPHY%NFIELDG)
INTEGER(KIND=JPIM) :: IHINTS(YDRQPHY%NFIELDG)

INTEGER(KIND=JPIM) :: IEND, IST, JBLOC, JFLD, IERR  , IOFF, IAUX

LOGICAL :: LLINTER  ! .TRUE. if actual interpolation on physical fields

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "fpcliphy.intfb.h"
#include "fpintphy.intfb.h"
#include "fpnilphy.intfb.h"
#include "fpsampl.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPOSHORPHY',0,ZHOOK_HANDLE)
ASSOCIATE(NFIELDG=>YDRQPHY%NFIELDG, YDRQCLI=>YDCLIMO%YRQPHY, &
 & NFPINPHY=>YDNAMFPINT%NFPINPHY, NFPSFXINT=>YDNAMFPINT%NFPSFXINT, &
 & NFPBOYD=>YDNAMFPSCI%NFPBOYD, NFPCLI=>YDNAMFPSCI%NFPCLI, NFPMASK=>YDNAMFPSCI%NFPMASK, &
 & GFP=>YDAFN%GFP, LFPOSHOR=>YDFPGEOMETRY%LFPOSHOR, NSLWIDE=>YDFPSTRUCT%NSLWIDE, &
 & NAFPB1=>YDFPSTRUCT%NASLB1, NFPBLOCS_DEP=>YDFPGEO_DEP%NFPBLOCS, NFPROMA_DEP=>YDFPGEO_DEP%NFPROMA, &
 & NFPEND_DEP=>YDFPGEO_DEP%NFPEND,NFPNUMD_DEP=>YDFPGEO_DEP%NFPNUMD)

!     ------------------------------------------------------------------

CALL GSTATS(1911,0)

!     ------------------------------------------------------------------

!*       1. PREPARATIONS
!           ------------

LLINT(:)=.FALSE.

! * Compute LLCLI:
LLCLI(:)=.FALSE.
IF (NFPCLI > 0) THEN
  CALL FPCLIPHY(YDRQCLI,GFP,NFIELDG,KCOD,IST,IEND,NFPROMA_DEP,LDCLI=LLCLI)
ENDIF

! * Compute LLNIL:
LLNIL(:)=.FALSE.
IF (LFPOSHOR) THEN
  LLINTER=NFPINPHY > 0
  IAUX=0
  CALL FPNILPHY(YDRQCLI,YDAFN,NFIELDG,KCOD,IST,IEND,NFPROMA_DEP,IAUX,LLINTER,YDRQAUX=YDRQAUX,LDNIL=LLNIL)
ENDIF

!     ------------------------------------------------------------------

!*       2. CALCULATIONS DOING HORIZONTAL INTERPOLATIONS OR SAMPLING
!           --------------------------------------------------------

! These calculations must be done on the departure geometry DM-environment.

!$OMP PARALLEL PRIVATE(JBLOC,IST,IEND,LLINT,JFLD,IOFF)
!$OMP DO SCHEDULE(DYNAMIC,1)
DO JBLOC=1,NFPBLOCS_DEP

  IST =1
  IEND=NFPEND_DEP(JBLOC)
  IOFF=(JBLOC-1)*NFPROMA_DEP

    IF (NFPBOYD /= 0) THEN
      ! Initialize non-interpolated fields to avoid computation on NaN values during biperiodicization : 
      DO JFLD=1,NFIELDG
        IF (LLCLI(JFLD).OR.LLNIL(JFLD)) THEN
          PFP(IST:IEND,JFLD,JBLOC)=0._JPRB
        ENDIF
      ENDDO
    ENDIF

    ! Compute pronostic fields : 
    IF (LFPOSHOR) THEN
      ! Interpolate fields :  
      CALL FPINTPHY(YDRQPHY,YDNAMFPINT,YDAFN,YDFPSUW,YDFPWSTD,NFPMASK,NAFPB1,NSLWIDE,NFIELDG,IST,IEND, &
       & PFP(:,:,JBLOC),PBUF1,NFPROMA_DEP,KFLDIN,LLCLI,LLNIL,LLINT,JBLOC)
    ELSE
      ! Sampling : 
      CALL FPSAMPL(YDNAMFPINT,YDFPWSTD,NAFPB1,NFIELDG,IST,IEND,NFPROMA_DEP,KFLDIN,PBUF1,PFP(:,:,JBLOC),&
       & LLCLI,LLINT,JBLOC)
    ENDIF


ENDDO
!$OMP END DO
!$OMP END PARALLEL

!     ------------------------------------------------------------------

!*       3. SOME CHECKINGS ABOUT LLCLI, LLNIL, LLINT
!           ----------------------------------------

! * Check that all fields have been computed once and only once:
  IERR=0
  DO JFLD=1,NFIELDG
    IHINTS(JFLD)=0
    IF (LLCLI(JFLD)) IHINTS(JFLD)=IHINTS(JFLD)+1
    IF (LLNIL(JFLD)) IHINTS(JFLD)=IHINTS(JFLD)+1
    IF (LLINT(JFLD)) IHINTS(JFLD)=IHINTS(JFLD)+1
    IF (IHINTS(JFLD) > 1) IERR=1
  ENDDO
  IF (IERR == 1) CALL ABOR1('FPOSHORPHY : INTERNAL ERROR ON IHINTS')

CALL GSTATS(1911,1)

!-----------------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FPOSHORPHY',1,ZHOOK_HANDLE)
END SUBROUTINE FPOSHORPHY
