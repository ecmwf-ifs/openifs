! (C) Copyright 1989- Meteo-France.

SUBROUTINE CPCLIMI(YDCLIMO,YDNAMFPINT,YDNAMFPSCI,YDAFN,KFPXFLD,YDFPGEOMETRY,YDFPWSTD,YDFPSTRUCT,YDGEOMETRY,YDSURF,YDMODEL,&
 & KFLD,KOD,KORDER,PHALO,PFPOUT)

!**** *CPCLIMI*  - Interpolate climatology fields

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!       *CALL* *CPCLIMI*

!        EXPLICIT ARGUMENTS
!        --------------------

!     INPUT:
!     ------
!        KFPXFLD : maximum number of fields to extract at a time
!        KFLD         : number of gridpoint fields in output gridpoint buffer 
!        KOD          : internal fields code
!        KORDER       : order of derivative of the fields

!     OUTPUT:
!     -------
!        PHALO        : interpolation buffer
!        PFPOUT       : output gridpoint buffer

!        IMPLICIT ARGUMENTS
!        --------------------
!        See lower #include.

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!         Called by SUFPWFPBUF.

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!      K. YESSAD (MARCH 1997) after routine CPCLIMO.

!     MODIFICATIONS.
!     --------------
!      G.Mozdzynski 02-10-01: support for radiation on-demand comms
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!      R. El Khatib  20-May-2005 NFPWIDE moved to YOMWFPDS
!      R. El Khatib  17-Oct-2008 bf for ALADIN fullpos with equal distribution
!                             of C+I+E (RDISTR_E=1.)
!      K. Yessad Dec 2008: merge SLCOMM+SLCOMM1 -> SLCOMM.
!      K. Yessad Dec 2008: merge the different (E)SLEXTPOL.. -> (E)SLEXTPOL.
!      G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!      G. Mozdzynski (May 2012): further cleaning
!      R. El Khatib  04-Dec-2012 Fix bounds checking issues created by recent "cleanings"
!      R. El Khatib 27-Sep-2013 Boyd periodization in Fullpos-2 (differentiation
!      between interpolation grid and output grid)
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      R. El Khatib 27-Jul-2016 interpolations over C+I+E
!      R. El Khatib 09-Sep-2016 Cleaning + temporary option for Boyd
!      E.Dutra/G.Arduini Jan 2018: change of HPOS dimension of PSP_SG, snow multi-layer
!     ------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN             , ONLY : NULOUT
USE YOMCT0             , ONLY : LELAM
USE YOMMP0             , ONLY : NPROC
USE YOMFPC             , ONLY : TNAMFPSCI, TNAMFPINT, LALLOFP
USE YOMWFPB            , ONLY : TFPWSTD, TFPSUW
USE EINT_MOD           , ONLY : SL_STRUCT
USE YOMAFN             , ONLY : TAFN
USE YOMFP4L      , ONLY : TRQFP
USE TYPE_FPRQPHYS      , ONLY : ALLOCATE_FPRQPHY
USE YOMFPGEOMETRY      , ONLY : TFPGEOMETRY, LFPDISTRIB
USE TYPE_FPOSBUF       , ONLY : FPOSBUF

!-----------------------------------------------------------------------------

IMPLICIT NONE

TYPE(FPOSBUF)     ,INTENT(IN) :: YDCLIMO
TYPE(TNAMFPINT)   ,INTENT(IN) :: YDNAMFPINT
TYPE(TNAMFPSCI)   ,INTENT(IN) :: YDNAMFPSCI
TYPE(TAFN)        ,INTENT(IN) :: YDAFN
INTEGER(KIND=JPIM),INTENT(IN) :: KFPXFLD
TYPE (TFPGEOMETRY),INTENT(IN) :: YDFPGEOMETRY
TYPE (TFPWSTD)    ,INTENT(IN) :: YDFPWSTD
TYPE(SL_STRUCT)   ,INTENT(IN) :: YDFPSTRUCT
TYPE(GEOMETRY)    ,INTENT(IN) :: YDGEOMETRY
TYPE(TSURF)       ,INTENT(IN) :: YDSURF
TYPE(MODEL)       ,INTENT(IN) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN) :: KFLD
INTEGER(KIND=JPIM),INTENT(IN) :: KOD (KFLD)
INTEGER(KIND=JPIM),INTENT(IN) :: KORDER(KFLD)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PHALO(YDFPSTRUCT%NASLB1,KFLD) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PFPOUT(YDFPGEOMETRY%YFPGEO_DEP%NFPROMA,KFLD,YDFPGEOMETRY%YFPGEO_DEP%NFPBLOCS)
!-----------------------------------------------------------------------------

LOGICAL :: LLINC, LLCTLCLIM, LLINTERPOL
INTEGER(KIND=JPIM) :: JFLD, JKGLO, IBL, JROF, IST, IEND, JBLOC, IRFPOS, IOFF
INTEGER(KIND=JPIM) :: IDUMARR(2)

REAL(KIND=JPRB) :: ZCORE(YDGEOMETRY%YRDIM%NPROMA,KFLD), ZPARITY(KFLD)
REAL(KIND=JPRB) :: ZROW(YDFPGEOMETRY%YFPGEO%NFPROMA,KFLD,YDFPGEOMETRY%YFPGEO%NFPBLOCS)
REAL(KIND=JPRB) :: ZAUX(0,0), ZCLI(0,0:0)
REAL(KIND=JPRB) :: ZWSXI, ZWDXI

TYPE (TFPSUW) :: YLFPSUW
TYPE (TNAMFPSCI) :: YLNAMFPSCI, YLNAMFPSCI2
TYPE (TRQFP) :: YLRQCLI, YLRQAUX

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------------

#include "hpos.intfb.h"
#include "fptensor.intfb.h"
#include "fphalo.intfb.h"
#include "slcomm.intfb.h"
#include "slextpol.intfb.h"
#include "eslextpol.intfb.h"
#include "fposhorphy.intfb.h"
#include "fptrdtoa.intfb.h"
#include "fptratod.intfb.h"
#include "fposhorlagphy.intfb.h"
#include "fpcorphy.intfb.h"

!-----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPCLIMI',0,ZHOOK_HANDLE)
ASSOCIATE(GFP_PHYDS=>YDAFN%GFP_PHYDS,LFPOSHOR=>YDFPGEOMETRY%LFPOSHOR, YDFPUSERGEO=>YDFPGEOMETRY%YFPUSERGEO,&
 & YDFPGEO_DEP=>YDFPGEOMETRY%YFPGEO_DEP, YDFPGEO=>YDFPGEOMETRY%YFPGEO, YDFPGIND=>YDFPGEOMETRY%YFPGIND, &
 & YDDIM=>YDGEOMETRY%YRDIM,YDGEM=>YDGEOMETRY%YRGEM, &
 & SD_SM=>YDSURF%SD_SM, SD_VA=>YDSURF%SD_VA, SD_VC=>YDSURF%SD_VC, SD_VD=>YDSURF%SD_VD, SD_VF=>YDSURF%SD_VF, SD_VP=>YDSURF%SD_VP, &
 & SD_VV=>YDSURF%SD_VV, SD_VX=>YDSURF%SD_VX, SD_WS=>YDSURF%SD_WS, SD_X2=>YDSURF%SD_X2, SP_CI=>YDSURF%SP_CI, SP_RR=>YDSURF%SP_RR, &
 & SP_SB=>YDSURF%SP_SB, SP_SG=>YDSURF%SP_SG, SP_SL=>YDSURF%SP_SL, SD_OC=>YDSURF%SD_OC)
ASSOCIATE(NFPBOYD=>YDNAMFPSCI%NFPBOYD, NFPCLI=>YDNAMFPSCI%NFPCLI, NFPINPHY=>YDNAMFPINT%NFPINPHY, &
 & NPROMA=>YDDIM%NPROMA, NGPTOT=>YDGEM%NGPTOT, &
 & NFPBLOCS_DEP=>YDFPGEO_DEP%NFPBLOCS, NFPROMA_DEP=>YDFPGEO_DEP%NFPROMA, &
 & NFPEND_DEP=>YDFPGEO_DEP%NFPEND, NFPRGPL_DEP=>YDFPGEO_DEP%NFPRGPL, &
 & NFPBLOCS=>YDFPGEO%NFPBLOCS, NFPROMA=>YDFPGEO%NFPROMA, NFPEND=>YDFPGEO%NFPEND)
!-----------------------------------------------------------------------------

!*       4.1 Setup request

WRITE(NULOUT,'('' CPCLIMI CALLED '')')

CALL ALLOCATE_FPRQPHY(GFP_PHYDS,YLRQCLI,SIZE(YDFPUSERGEO),KFLD,KOD,NULOUT,LALLOFP)
YLRQAUX%NFIELDG=0

!*       4.2 Get halo and interpolated/output field

! Lake option is disabled because the lake/island index is not yet computed (it needs the target land-sea mask itself !)

! Disable climatology for interpolations since we are computing climatology-like data (correct ??) and at least for Boyd extension
! corrections:
YLNAMFPSCI=YDNAMFPSCI
YLNAMFPSCI%NFPCLI=0
YLNAMFPSCI%NFPLAKE=0
LLCTLCLIM=.FALSE.

! For corrections on target grid climatology may naturally be used (we want to recover the clim data if available)
YLNAMFPSCI2=YDNAMFPSCI
YLNAMFPSCI2%NFPLAKE=0

ZPARITY(:)=1.0_JPRB

ZWSXI=HUGE(1._JPRB)
ZWDXI=HUGE(1._JPRB)

CALL GSTATS(1434,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,IEND,IBL,JROF,JFLD,IOFF,ZCORE)
DO JKGLO=1,NGPTOT,NPROMA

  IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
  IBL=(JKGLO-1)/NPROMA+1

  CALL HPOS(YLRQCLI,YDCLIMO%YRQPHY,YLNAMFPSCI,YDAFN,LFPOSHOR,YDGEOMETRY,YDSURF,YDMODEL,IEND,JKGLO,&
   & SP_SB(:,:,:,IBL),SP_SG(:,:,:,IBL),SP_SL(:,:,IBL),SP_RR(:,:,IBL),SP_CI(:,:,IBL),&
   & SD_WS(:,:,IBL),SD_VD(:,:,IBL),SD_VX(:,:,IBL),SD_VF(:,:,IBL),SD_VV(:,:,IBL),&
   & SD_VP(:,:,IBL),SD_VA(:,:,IBL),SD_VC(:,:,IBL),SD_X2(:,:,IBL),SD_SM(:,:,:,IBL),&
   & SD_OC(:,:,IBL), &
   & ZWSXI,ZWDXI,ZCORE)

  IF (LFPOSHOR) THEN
    CALL FPTENSOR(1,IEND,NPROMA,KFLD,KORDER,YDGEOMETRY%YRGSGEOM(IBL)%GM,ZCORE,ZPARITY)
  ENDIF
  CALL FPHALO(YDFPSTRUCT,NPROMA,KFLD,JKGLO,IEND,KFLD,ZCORE,PHALO)
ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1434,1)


IF (YDFPSTRUCT%NSLWIDE > 0) THEN

!*          COMMUNICATION BETWEEN PROCESSORS TO MAKE HALOS

  IF (NPROC > 1) THEN
    LLINC=.FALSE.
    CALL SLCOMM(YDFPSTRUCT,IDUMARR,KFLD,LLINC,0,PHALO)
  ENDIF

!*          POLAR OR E-ZONE EXTENSIONS

  IF (LELAM) THEN
    CALL ESLEXTPOL(YDGEOMETRY,YDFPSTRUCT,KFLD,IDUMARR,1,PHALO)
  ELSE
    CALL SLEXTPOL(YDGEOMETRY%YRDIM,YDFPSTRUCT,KFLD,IDUMARR,1,ZPARITY,PHALO)
  ENDIF

ENDIF



LLINTERPOL=NFPINPHY > 0

IF (LFPDISTRIB(YDFPGEOMETRY)) THEN

  CALL FPOSHORPHY(YLRQCLI,YDCLIMO,YDNAMFPINT,YLNAMFPSCI,YDAFN,YDFPGEOMETRY,YLFPSUW,YDFPWSTD,YDFPSTRUCT,YDFPGEO_DEP,KFLD,PHALO, &
   & YLRQCLI%ICOD,PFPOUT)
  CALL FPTRDTOA(KFPXFLD,YDFPGEO_DEP,YDFPGIND,YDFPGEO,KFLD,PFPOUT,ZROW)
  ! Correction over the arrival geometry
  CALL FPOSHORLAGPHY(YLRQCLI,YDCLIMO,ZWSXI,ZWDXI,YLNAMFPSCI2,YDAFN,YDFPGEO,YDMODEL%YRML_PHY_EC%YREPHY,YDMODEL%YRML_PHY_MF%YRPHY, &
   & YDMODEL%YRML_PHY_MF%YRPHY1,LLINTERPOL,ZROW,LFPOSHOR)
  IF (NFPBOYD /=0) THEN
    IRFPOS=0
    ! Correction over the (larger) departure geometry
    DO JBLOC=1,NFPBLOCS_DEP
      IST =1
      IEND=NFPEND_DEP(JBLOC)
      CALL FPCORPHY(YDCLIMO%YRQPHY,YLNAMFPSCI,YDAFN,YDMODEL%YRML_PHY_EC%YREPHY,YDMODEL%YRML_PHY_MF%YRPHY, &
       & YDMODEL%YRML_PHY_MF%YRPHY1,IST,IEND,NFPROMA_DEP,LFPOSHOR,KFLD,YLRQCLI%ICOD,PFPOUT(:,:,JBLOC),YLRQAUX%NFIELDG,ZAUX, &
       & YLRQAUX,IRFPOS,ZCLI(:,:),ZWSXI,ZWDXI)
      !   Though the call to fpcorphy should be theorically correct, 
      !   setting to zero seems to lead to smoother fields in the extension zone
      !   - as it was coded before - but it could have been 1. as well ...
      !   this option is not documented yet. REK
      IF (NFPBOYD==1) THEN
        PFPOUT(IST:IEND,:,JBLOC)=0._JPRB
      ELSEIF (NFPBOYD==-1) THEN
        PFPOUT(IST:IEND,:,JBLOC)=1._JPRB
      ENDIF
    ENDDO
  ENDIF
  ! Return to departure geometry (=> Combine corrections for Boyd)
  CALL FPTRATOD(KFPXFLD,YDFPGEO_DEP,YDFPGIND,YDFPGEO,KFLD,ZROW,PFPOUT)

ELSE

  CALL FPOSHORPHY(YLRQCLI,YDCLIMO,YDNAMFPINT,YLNAMFPSCI,YDAFN,YDFPGEOMETRY,YLFPSUW,YDFPWSTD,YDFPSTRUCT,YDFPGEO_DEP,KFLD,PHALO, &
   & YLRQCLI%ICOD,PFPOUT)
  CALL FPOSHORLAGPHY(YLRQCLI,YDCLIMO,ZWSXI,ZWDXI,YLNAMFPSCI2,YDAFN,YDFPGEO,YDMODEL%YRML_PHY_EC%YREPHY,YDMODEL%YRML_PHY_MF%YRPHY, &
   & YDMODEL%YRML_PHY_MF%YRPHY1,LLINTERPOL,PFPOUT,LFPOSHOR)

ENDIF

! ---------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPCLIMI',1,ZHOOK_HANDLE)
END SUBROUTINE CPCLIMI
