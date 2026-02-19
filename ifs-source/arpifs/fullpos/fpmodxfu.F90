! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPMODXFU(YDGFP,YDXFU,KFPCOD,KFLD,KFLDPTR)

!**** *FPMODXFU*  - Preliminary modification of instantaneous fluxes 
!                   for Post-Processing
!     PURPOSE.
!     --------
!        To prepare input fields before horizontal interpolations

!**   INTERFACE.
!     ----------
!       *CALL* *FPMODXFU*

!        EXPLICIT ARGUMENTS
!        --------------------
!        YDXFU     - fluxes structure
!        KFPCOD    - Fullpos fields code
!        KFLD      - number of fields to modify
!        KFLDPTR   - array of field pointers for extraction

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
!        RYAD EL KHATIB *METEO-FRANCE*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 96-08-06
!        R. El Khatib : 01-08-07 Pruning options
!    R. El Khatib : 03-04-17 Fullpos improvemnts
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        Y.Seity      15-Jan-2007 add Graupel and Hail
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMAFN   , ONLY : ALL_FPDSPHY_TYPES

USE YOMXFU   , ONLY : TXFU

IMPLICIT NONE

TYPE(ALL_FPDSPHY_TYPES),  INTENT(IN) :: YDGFP
TYPE(TXFU)        , INTENT(IN)  :: YDXFU
INTEGER(KIND=JPIM), INTENT(IN)  :: KFPCOD(:)
INTEGER(KIND=JPIM), INTENT(OUT) :: KFLD
INTEGER(KIND=JPIM), INTENT(OUT) :: KFLDPTR(YDXFU%NFDXFU)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!*    1. FIND FIELDS TO MODIFY 
!        ---------------------

IF (LHOOK) CALL DR_HOOK('FPMODXFU',0,ZHOOK_HANDLE)
ASSOCIATE(MXPLCL=>YDXFU%YXFUPT%MXPLCL, MXPLSL=>YDXFU%YXFUPT%MXPLSL, MXPLCN=>YDXFU%YXFUPT%MXPLCN, &
 & MXPLSN=>YDXFU%YXFUPT%MXPLSN, MXPLCG=>YDXFU%YXFUPT%MXPLCG, MXPLSG=>YDXFU%YXFUPT%MXPLSG, &
 & MXPLCH=>YDXFU%YXFUPT%MXPLCH, MXPLSH=>YDXFU%YXFUPT%MXPLSH)

KFLD=0
IF (ANY(KFPCOD(:)==YDGFP%XLSP%ICOD)) THEN   ! Large scale rain
  KFLD=KFLD+1
  KFLDPTR(KFLD)=MXPLSL
ENDIF
IF (ANY(KFPCOD(:)==YDGFP%XCP%ICOD)) THEN    ! Convective rain
  KFLD=KFLD+1
  KFLDPTR(KFLD)=MXPLCL
ENDIF
IF (ANY(KFPCOD(:)==YDGFP%XLSS%ICOD)) THEN   ! Large Scale Snow fall
  KFLD=KFLD+1
  KFLDPTR(KFLD)=MXPLSN
ENDIF
IF (ANY(KFPCOD(:)==YDGFP%XCSF%ICOD)) THEN   ! Convective  Snow fall
  KFLD=KFLD+1
  KFLDPTR(KFLD)=MXPLCN
ENDIF
IF (ANY(KFPCOD(:)==YDGFP%XLSG%ICOD)) THEN   ! Large Scale Graupel fall
  KFLD=KFLD+1
  KFLDPTR(KFLD)=MXPLSG
ENDIF
IF (ANY(KFPCOD(:)==YDGFP%XCSG%ICOD)) THEN   ! Convective  Graupel fall
  KFLD=KFLD+1
  KFLDPTR(KFLD)=MXPLCG
ENDIF
IF (ANY(KFPCOD(:)==YDGFP%XLSH%ICOD)) THEN   ! Large Scale Hail fall
  KFLD=KFLD+1
  KFLDPTR(KFLD)=MXPLSH
ENDIF
IF (ANY(KFPCOD(:)==YDGFP%XCSH%ICOD)) THEN   ! Convective  Hail fall
  KFLD=KFLD+1
  KFLDPTR(KFLD)=MXPLCH
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FPMODXFU',1,ZHOOK_HANDLE)
END SUBROUTINE FPMODXFU
