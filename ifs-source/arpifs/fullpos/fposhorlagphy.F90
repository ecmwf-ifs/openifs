! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPOSHORLAGPHY(YDRQPHY,YDCLIMO,PWSXI,PWDXI,YDNAMFPSCI,YDAFN,YDFPGEO,YDEPHY,YDPHY,YDPHY1,LDINTERPOL,PFP,LDFPOSHOR, &
 & KWIC,YDAUX)

!**** *FPOSHORLAGPHY*  - HORIZONTAL POST-PROCESSING - Lagged part for physical fields

!     PURPOSE.
!     --------
!        PERFORM THE CORRECTIONS AFTER HORIZONTAL INTERPOLATIONS

!        Computations are DM-local if distributed memory.

!**   INTERFACE.
!     ----------
!       *CALL* *FPOSHORLAGPHY*

!        EXPLICIT ARGUMENTS
!        --------------------

!     PFP    : post-processed fields 
!     LDFPOSHOR : .TRUE. if horizontal interpolations mechanism
!     LDINTERPOL : .TRUE. if actual interpolation -not the nearest point is taken
!     KWIC   : isolated lake/island indicator
!     YDAUX  : auxilary fields  structure

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
!      ORIGINAL : 24-09-2012 from FPOSHOR

!     MODIFICATIONS.
!     --------------
!      R. El Khatib 13-Dec-2012 Fullpos buffers reshaping, memory savings
!      R. El Khatib 27-Sep-2013 Cleaning
!      R. El Khatib 28-Jul-2016 Recode LWIDER_DOM
!     ------------------------------------------------------------------

USE YOMPHY1  , ONLY : TPHY1
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMFPGEO , ONLY : TFPGEO
USE YOMFPC   , ONLY : TNAMFPSCI
USE TYPE_FPOSBUF, ONLY : FPOSBUF
USE YOMPHY   , ONLY : TPHY
USE YOEPHY   , ONLY : TEPHY
USE YOMAFN   , ONLY : TAFN
USE YOMFP4L, ONLY : TRQFP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TRQFP),  INTENT(IN) :: YDRQPHY
TYPE (FPOSBUF),    INTENT(IN) :: YDCLIMO
REAL(KIND=JPRB)   ,INTENT(IN) :: PWSXI
REAL(KIND=JPRB)   ,INTENT(IN) :: PWDXI
TYPE (TNAMFPSCI),  INTENT(IN) :: YDNAMFPSCI
TYPE(TAFN),        INTENT(IN) :: YDAFN
TYPE(TFPGEO),      INTENT(IN) :: YDFPGEO
TYPE(TEPHY)       ,INTENT(IN) :: YDEPHY
TYPE(TPHY)        ,INTENT(IN) :: YDPHY
TYPE(TPHY1)       ,INTENT(IN) :: YDPHY1
LOGICAL           ,INTENT(IN) :: LDINTERPOL
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFP(YDFPGEO%NFPROMA,YDRQPHY%NFIELDG,YDFPGEO%NFPBLOCS)
LOGICAL           ,INTENT(IN)    :: LDFPOSHOR
INTEGER(KIND=JPIM),INTENT(IN), OPTIONAL :: KWIC(:,:)
TYPE (FPOSBUF),  INTENT(IN), OPTIONAL :: YDAUX

!     ------------------------------------------------------------------

!     * output climatology and geometry:
REAL(KIND=JPRB) :: Z1SGM2(YDFPGEO%NFPROMA),ZAUX(0,0)

INTEGER(KIND=JPIM) :: IEND, IST, JBLOC, JFP
TYPE (TRQFP) :: YLRQAUX

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "fpcliphy.intfb.h"
#include "fpcorphy.intfb.h"
#include "fpgeophy.intfb.h"
#include "fpnilphy.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPOSHORLAGPHY',0,ZHOOK_HANDLE)
ASSOCIATE(NFIELDG=>YDRQPHY%NFIELDG, YDRQCLI=>YDCLIMO%YRQPHY, RFPBUF=>YDCLIMO%FPBUF, NFPVSO=>YDCLIMO%YRQPHY%NFIELDG, &
 & NFPCLI=>YDNAMFPSCI%NFPCLI, GFP=>YDAFN%GFP, GFP_PHYDS=>YDAFN%GFP_PHYDS, &
 & NFPBLOCS=>YDFPGEO%NFPBLOCS, NFPROMA=>YDFPGEO%NFPROMA, NFPEND=>YDFPGEO%NFPEND, &
 & RFPGM=>YDFPGEO%RFPGM, RFPNORX=>YDFPGEO%RFPNORX, RFPNORY=>YDFPGEO%RFPNORY)

!     ------------------------------------------------------------------

CALL GSTATS(1911,0)

!     ------------------------------------------------------------------

!*       1. CORRECTION OF INTERPOLATED/BIPERIODICIZED DATA 
!           ----------------------------------------------

YLRQAUX%NFIELDG = 0

DO JBLOC=1,NFPBLOCS

  IST =1
  IEND=NFPEND(JBLOC)

  ! Overwrite constant fields with climatology : 
  IF (NFPCLI > 0) THEN
    CALL FPCLIPHY(YDRQCLI,GFP,NFIELDG,YDRQPHY%ICOD,IST,IEND,NFPROMA,PCLI=RFPBUF(:,:,JBLOC),PROW=PFP(:,:,JBLOC))
  ENDIF

  ! Compute pronostic fields : 
  IF (LDFPOSHOR) THEN
    ! Compute fields that cannot be interpolated : 
    IF (PRESENT(YDAUX)) THEN
      CALL FPNILPHY(YDRQCLI,YDAFN,NFIELDG,YDRQPHY%ICOD,IST,IEND,NFPROMA,YDAUX%YRQPHY%NFIELDG,LDINTERPOL, &
       & YDRQAUX=YDAUX%YRQPHY,PAUX=YDAUX%FPBUF(:,:,JBLOC),PROW=PFP(:,:,JBLOC))
    ELSE
      CALL FPNILPHY(YDRQCLI,YDAFN,NFIELDG,YDRQPHY%ICOD,IST,IEND,NFPROMA,YLRQAUX%NFIELDG,LDINTERPOL,PROW=PFP(:,:,JBLOC))
    ENDIF
  ENDIF

  ! Corrections after interpolations or climatology usage : 
  IF ((LDFPOSHOR.OR.NFPCLI > 0) .AND. (YDEPHY%LEPHYS.OR.YDPHY%LMPHYS)) THEN
    IF (PRESENT(YDAUX)) THEN
      IF (PRESENT(KWIC)) THEN
        CALL FPCORPHY(YDRQCLI,YDNAMFPSCI,YDAFN,YDEPHY,YDPHY,YDPHY1,IST,IEND,NFPROMA,LDFPOSHOR,NFIELDG,YDRQPHY%ICOD, &
         & PFP(:,:,JBLOC),YDAUX%YRQPHY%NFIELDG,YDAUX%FPBUF(:,:,JBLOC),YDAUX%YRQPHY,NFPVSO,RFPBUF(:,:,JBLOC),PWSXI,PWDXI, &
         & KWIC=KWIC(:,JBLOC))
      ELSE
        CALL FPCORPHY(YDRQCLI,YDNAMFPSCI,YDAFN,YDEPHY,YDPHY,YDPHY1,IST,IEND,NFPROMA,LDFPOSHOR,NFIELDG,YDRQPHY%ICOD, &
         & PFP(:,:,JBLOC),YDAUX%YRQPHY%NFIELDG,YDAUX%FPBUF(:,:,JBLOC),YDAUX%YRQPHY,NFPVSO,RFPBUF(:,:,JBLOC),PWSXI,PWDXI)
      ENDIF
    ELSE
      IF (PRESENT(KWIC)) THEN
        CALL FPCORPHY(YDRQCLI,YDNAMFPSCI,YDAFN,YDEPHY,YDPHY,YDPHY1,IST,IEND,NFPROMA,LDFPOSHOR,NFIELDG,YDRQPHY%ICOD, &
         & PFP(:,:,JBLOC),YLRQAUX%NFIELDG,ZAUX,YLRQAUX,NFPVSO,RFPBUF(:,:,JBLOC),PWSXI,PWDXI,KWIC=KWIC(:,JBLOC))
      ELSE
        CALL FPCORPHY(YDRQCLI,YDNAMFPSCI,YDAFN,YDEPHY,YDPHY,YDPHY1,IST,IEND,NFPROMA,LDFPOSHOR,NFIELDG,YDRQPHY%ICOD, &
         & PFP(:,:,JBLOC),YLRQAUX%NFIELDG,ZAUX,YLRQAUX,NFPVSO,RFPBUF(:,:,JBLOC),PWSXI,PWDXI)
      ENDIF
    ENDIF
  ENDIF

  ! Go to output geometry
  IF (LDFPOSHOR) THEN
    DO JFP=IST,IEND
      Z1SGM2(JFP)=1.0_JPRB/(RFPGM(JFP,JBLOC)**2)
    ENDDO
    CALL FPGEOPHY(YDRQPHY,IST,IEND,NFPROMA,RFPNORX(:,JBLOC),RFPNORY(:,JBLOC),Z1SGM2,NFIELDG,GFP_PHYDS(:)%IORDR,PFP(:,:,JBLOC))
  ENDIF

ENDDO

CALL GSTATS(1911,1)

!-----------------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FPOSHORLAGPHY',1,ZHOOK_HANDLE)
END SUBROUTINE FPOSHORLAGPHY
