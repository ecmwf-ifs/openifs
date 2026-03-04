! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPOSHOR(YDQTYPE,YDNAMFPINT,YDTFP_DYNDS,YDFPGEOMETRY,YDFPWSTD,YDFPSTRUCT,YDFPGEO_DEP,KFLDIN,PBUF1,KFLDOUT,PFP)

!**** *FPOSHOR*  - HORIZONTAL POST-PROCESSING FOR DYNAMICAL FIELDS

!     PURPOSE.
!     --------
!        PERFORM THE HORIZONTAL INTERPOLATIONS AND THE PBL DISPLACEMENT

!        Computations are DM-local if distributed memory.

!**   INTERFACE.
!     ----------
!       *CALL* *FPOSHOR*

!        EXPLICIT ARGUMENTS
!        --------------------

!     All dummy arguments are input ones.

!     KFLDIN : number of fields in input buffer
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
USE EINT_MOD     , ONLY : SL_STRUCT
USE YOMWFPB      , ONLY : TFPWSTD
USE FULLPOS_MIX  , ONLY : FULLPOS_TYPE
USE YOMFPC       , ONLY : TNAMFPINT
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(IN) :: YDQTYPE
TYPE(TNAMFPINT),   INTENT(IN) :: YDNAMFPINT
TYPE(FULLPOS_TYPE),INTENT(IN) :: YDTFP_DYNDS(:)
TYPE (TFPGEOMETRY),INTENT(IN) :: YDFPGEOMETRY
TYPE (TFPWSTD),    INTENT(IN) :: YDFPWSTD
TYPE(SL_STRUCT),   INTENT(IN) :: YDFPSTRUCT
TYPE (TFPGEO),     INTENT(IN) :: YDFPGEO_DEP
INTEGER(KIND=JPIM),INTENT(IN) :: KFLDIN
REAL(KIND=JPRB),   INTENT(IN) :: PBUF1(YDFPSTRUCT%NASLB1*KFLDIN) 
INTEGER(KIND=JPIM),INTENT(IN) :: KFLDOUT
REAL(KIND=JPRB),   INTENT(OUT) :: PFP(YDFPGEO_DEP%NFPROMA,KFLDOUT,YDFPGEO_DEP%NFPBLOCS)

!     ------------------------------------------------------------------

!     Handling horizontal interpolations of physical fields : 
!     LLCLI=.T. if the field is overwritten by the climatology
!     LLINT=.T. if the field is to be interpolated

LOGICAL :: LLCLI(KFLDOUT), LLINT(KFLDOUT)

INTEGER(KIND=JPIM) :: IEND, IST, JBLOC, JFLD, IOFF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "fpintdyn.intfb.h"
#include "fpsampl.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPOSHOR',0,ZHOOK_HANDLE)
ASSOCIATE(LFPOSHOR=>YDFPGEOMETRY%LFPOSHOR, &
 & NAFPB1=>YDFPSTRUCT%NASLB1,NFPBLOCS_DEP=>YDFPGEO_DEP%NFPBLOCS, NFPROMA_DEP=>YDFPGEO_DEP%NFPROMA, &
 & NFPEND_DEP=>YDFPGEO_DEP%NFPEND,NFPNUMD_DEP=>YDFPGEO_DEP%NFPNUMD)

!     ------------------------------------------------------------------

CALL GSTATS(1911,0)

!     ------------------------------------------------------------------

!*       1. PREPARATIONS
!           ------------

LLINT(:)=.FALSE.
LLCLI(:)=.FALSE.

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

    IF (LFPOSHOR) THEN
      ! Interpolate fields :  
      CALL FPINTDYN(YDQTYPE,YDNAMFPINT,YDFPWSTD,NAFPB1,KFLDOUT,IST,IEND,YDTFP_DYNDS(:)%INTER,PBUF1, &
       & NFPROMA_DEP,KFLDIN,JBLOC,NFPNUMD_DEP(IOFF+IST:IOFF+IEND),PFP(:,:,JBLOC))
    ELSE
      ! Sampling on all fields : 
      CALL FPSAMPL(YDNAMFPINT,YDFPWSTD,NAFPB1,KFLDOUT,IST,IEND,NFPROMA_DEP,KFLDIN,PBUF1,PFP(:,:,JBLOC),&
       & LLCLI,LLINT,JBLOC)
    ENDIF

ENDDO
!$OMP END DO
!$OMP END PARALLEL

CALL GSTATS(1911,1)

!-----------------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FPOSHOR',1,ZHOOK_HANDLE)
END SUBROUTINE FPOSHOR
