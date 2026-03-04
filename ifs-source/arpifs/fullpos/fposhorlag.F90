! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPOSHORLAG(YDQTYPE,YDNAMFPSCI,YDAFN,YDFPGEO,KFLDOUT,PFP)

!**** *FPOSHORLAG*  - HORIZONTAL POST-PROCESSING - Lagged part for dynamical fields

!     PURPOSE.
!     --------
!        PERFORM THE CORRECTIONS AFTER HORIZONTAL INTERPOLATIONS

!        Computations are DM-local if distributed memory.

!**   INTERFACE.
!     ----------
!       *CALL* *FPOSHORLAG*

!        EXPLICIT ARGUMENTS
!        --------------------
!          PFP    : post-processed fields 

!     IMPLICIT ARGUMENTS
!     ------------------

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

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMFPGEO , ONLY : TFPGEO
USE YOMAFN   , ONLY : TAFN
USE YOMFPC   , ONLY : TNAMFPSCI
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(IN) :: YDQTYPE
TYPE(TNAMFPSCI)   , INTENT(IN) :: YDNAMFPSCI
TYPE(TAFN)        , INTENT(IN) :: YDAFN
TYPE(TFPGEO)      , INTENT(IN) :: YDFPGEO
INTEGER(KIND=JPIM), INTENT(IN) :: KFLDOUT
REAL(KIND=JPRB)   , INTENT(INOUT) :: PFP(YDFPGEO%NFPROMA,KFLDOUT,YDFPGEO%NFPBLOCS)

!     ------------------------------------------------------------------

! output geometry:
REAL(KIND=JPRB) :: Z1SGM2(YDFPGEO%NFPROMA)
INTEGER(KIND=JPIM) :: IEND, IST, JBLOC, JFP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "fpcordyn.intfb.h"
#include "fpgeo.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPOSHORLAG',0,ZHOOK_HANDLE)
ASSOCIATE(TFP=>YDAFN%TFP, TFP_DYNDS=>YDAFN%TFP_DYNDS, &
 & NFPBLOCS=>YDFPGEO%NFPBLOCS, NFPROMA=>YDFPGEO%NFPROMA, NFPEND=>YDFPGEO%NFPEND, &
 & RFPGM=>YDFPGEO%RFPGM, RFPNORX=>YDFPGEO%RFPNORX, RFPNORY=>YDFPGEO%RFPNORY)

!*       1. CORRECTION OF INTERPOLATED/BIPERIODICIZED DATA 
!           ----------------------------------------------

CALL GSTATS(1911,0)

DO JBLOC=1,NFPBLOCS
  IST =1
  IEND=NFPEND(JBLOC)
  DO JFP=IST,IEND
    Z1SGM2(JFP)=1.0_JPRB/(RFPGM(JFP,JBLOC)**2)
  ENDDO
  ! Corrections afer interpolations : 
  CALL FPCORDYN(YDQTYPE,YDNAMFPSCI,TFP,KFLDOUT,IST,IEND,PFP(:,:,JBLOC),NFPROMA)
  ! Go to output geometry : 
  CALL FPGEO(YDQTYPE,IST,IEND,NFPROMA,RFPNORX(:,JBLOC),RFPNORY(:,JBLOC),Z1SGM2,KFLDOUT,TFP_DYNDS(:)%IORDR,PFP(:,:,JBLOC))
ENDDO

CALL GSTATS(1911,1)

!-----------------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FPOSHORLAG',1,ZHOOK_HANDLE)
END SUBROUTINE FPOSHORLAG
