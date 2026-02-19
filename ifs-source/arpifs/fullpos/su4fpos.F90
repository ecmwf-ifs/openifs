! (C) Copyright 1989- Meteo-France.

SUBROUTINE SU4FPOS(YDNAMFPSCI,YDAFN,YDFPCNT,KFLEVG,KLEVSNOW,LDPHYS,KTIME,KFPDOM,CDFPDOM,YDFPVAB,YGFL,YDEAERATM, &
 & LDLAGGED,LDOCEDELAY,LDOCE,YDFPFIELDS,YDRQPHY,YDRQDYN)

!**** *SU4FPOS*  - INITIALIZE FULL-POS level 4

!     PURPOSE.
!     --------
!        INITIALIZE COMMON BLOCKS YOMFP4 AND PTRFP4, WHICH CONTAIN THE REQUESTS
!           FOR POST-PROCESSING AT A GIVEN TIME STEP.

!**   INTERFACE.
!     ----------
!       *CALL* *SU4FPOS*

!        EXPLICIT ARGUMENTS
!        --------------------
!        KTIME  : time in model
!        KFPDOM : number of subdomains
!        CDFPDOM: names of the subdomains
!        LDLAGGED : Lagged call
!        LDOCEDELAY : Post-processing physics fields at the end of the step
!        LDOCE   : Lagged call for ocean pp fields

!        IMPLICIT ARGUMENTS
!        --------------------
!        NONE.

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib 17-Jul-2013 FABEC post-processing
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE PARKIND1  , ONLY :  JPIM     ,JPRB
USE YOMHOOK   , ONLY :  LHOOK    ,DR_HOOK, JPHOOK
USE YOMLUN    , ONLY :  NULOUT
USE YOMFPCNT  , ONLY : TFPCNT
USE TYPE_FPFIELDS, ONLY : TFPFIELDS
USE YOM4FPOS  , ONLY : TRQFPDYN
USE YOMCT0    , ONLY : LECMWF
USE YOMVERT   , ONLY : TVAB
USE YOMAFN    , ONLY : TAFN
USE YOMFPC    , ONLY : TNAMFPL, TNAMFPSCI
USE YOMFP4L   , ONLY : TRQFP
USE YOM_YGFL  , ONLY : TYPE_GFLD
USE YOEAERATM , ONLY : TEAERATM

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TNAMFPSCI),   INTENT(IN)  :: YDNAMFPSCI
TYPE(TAFN),        INTENT(IN)  :: YDAFN
TYPE(TFPCNT),      INTENT(IN)  :: YDFPCNT
INTEGER(KIND=JPIM),INTENT(IN)  :: KFLEVG
LOGICAL,           INTENT(IN)  :: LDPHYS
INTEGER(KIND=JPIM),INTENT(IN)  :: KTIME
INTEGER(KIND=JPIM),INTENT(IN)  :: KFPDOM
CHARACTER(LEN=*),  INTENT(IN)  :: CDFPDOM(KFPDOM)
TYPE(TVAB),        INTENT(IN)  :: YDFPVAB
TYPE(TYPE_GFLD),   INTENT(IN)  :: YGFL
TYPE(TEAERATM),    INTENT(IN)  :: YDEAERATM
LOGICAL,           INTENT(IN)  :: LDLAGGED
LOGICAL,           INTENT(IN)  :: LDOCEDELAY
LOGICAL,           INTENT(IN)  :: LDOCE
TYPE (TFPFIELDS),  INTENT(OUT) :: YDFPFIELDS
TYPE(TRQFP),       INTENT(OUT) :: YDRQPHY
TYPE(TRQFPDYN),    INTENT(OUT) :: YDRQDYN
INTEGER(KIND=JPIM),INTENT(IN)  :: KLEVSNOW



!     ------------------------------------------------------------------

CHARACTER(LEN=120) :: CLYFILE
LOGICAL :: LLNOPPFIL
TYPE(TNAMFPL) :: YLNAMFPL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "ppreq.intfb.h"
#include "sufpc.intfb.h"
#include "sufpfields.intfb.h"
#include "sufpdyn.intfb.h"
#include "sufpphy.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SU4FPOS',0,ZHOOK_HANDLE)
ASSOCIATE(LFPCNT=>YDFPCNT%LFPCNT, LFPNAMELIST=>YDFPCNT%LFPNAMELIST, CNAM=>YDFPCNT%CNAM, &
 & TFP_DYNDS=>YDAFN%TFP_DYNDS)
!     ------------------------------------------------------------------

!*       1. INITIALIZE REQUESTS ON PHYSICS AND FLUXES
!           -----------------------------------------

WRITE(NULOUT,'(''--- Set up Fullpos fields list'')')
IF (LECMWF) THEN
  !*    Initialize Fullpos fields request
  IF (LFPCNT) THEN
    CALL PPREQ(KTIME,LLNOPPFIL,CLYFILE)
    CALL SUFPC(YDNAMFPL=YLNAMFPL,LDPRINT=.FALSE.,CDPATH=CLYFILE,LDNAMELIST=LFPNAMELIST)
  ELSEIF (CNAM == ' ') THEN
    CALL SUFPC(YDNAMFPL=YLNAMFPL,LDPRINT=.FALSE.)
  ELSE
    CALL SUFPC(YDNAMFPL=YLNAMFPL,LDPRINT=.FALSE.,CDPATH=CNAM)
  ENDIF
  CALL SUFPFIELDS(YLNAMFPL,KFLEVG,KLEVSNOW,YDFPFIELDS,YDFPVAB,LDPRINT=.TRUE.)
  LLNOPPFIL=.TRUE. ! pp file was read as an alternative standard namelist and fullpos has been reconstructed
ELSE
  !*    Initialize Fullpos fields request
  IF (CNAM == ' ') THEN
    CALL SUFPC(YDNAMFPL=YLNAMFPL,LDPRINT=.FALSE.)
  ELSE
    CALL SUFPC(YDNAMFPL=YLNAMFPL,LDPRINT=.FALSE.,CDPATH=CNAM)
  ENDIF
  CALL SUFPFIELDS(YLNAMFPL,KFLEVG,KLEVSNOW,YDFPFIELDS,YDFPVAB,LDPRINT=.TRUE.)
  IF (LFPCNT) THEN
    CALL PPREQ(KTIME,LLNOPPFIL,CLYFILE) ! a selection file is read as a complement to the standard namelist
  ELSE
    LLNOPPFIL=.TRUE.
  ENDIF
ENDIF


!*       1.2 PHYSICAL FIELDS

CALL SUFPPHY(YDAFN,LDPHYS,YDFPFIELDS,KFPDOM,CDFPDOM,LLNOPPFIL,CLYFILE,YGFL,YDEAERATM,LDLAGGED,LDOCEDELAY,LDOCE,YDRQPHY)

!*       2. INITIALIZE REQUESTS ON DYNAMICS
!           -------------------------------

CALL SUFPDYN(YDNAMFPSCI,YDFPFIELDS,KFPDOM,CDFPDOM,TFP_DYNDS(:)%IGRIB,TFP_DYNDS(:)%LLSRF,TFP_DYNDS(:)%CLNAME, &
 & LLNOPPFIL,CLYFILE,YDRQDYN)

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SU4FPOS',1,ZHOOK_HANDLE)
END SUBROUTINE SU4FPOS
