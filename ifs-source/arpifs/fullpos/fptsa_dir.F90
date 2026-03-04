! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPTSA_DIR(YDQTYPE,YDTFP,KRESOL,PSCA,PVOR,PDIV,PUMEAN,PVMEAN,PSPBFP,PUMEANFP,PVMEANFP)

!**** *FPTSA_DIR * - Full-Pos Transfer Spectral Array - Direct

!     Purpose
!     -------
!       To transfer spectral arrays out of transforms to the post-processing 
!       spectral array.

!**   Interface.  CALL FPTSA_DIR(...)
!     ---------- 

!     Explicit arguments : 
!     --------------------
!        KRESOL : resolution indicator
!        PSCA   : input scalar fields
!        PVOR   : input vorticity fields
!        PDIV   : input divergence fields
!        PUMEAN : input U - mean wind (for LAM)
!        PVMEAN : input V - mean wind (for LAM)
!        PSPBFP   : output spectral fields
!        PUMEANFP : output U - mean wind (for LAM)
!        PVMEANFP : output V - mean wind (for LAM)

!     Externals.
!     ----------
!       None.

!     Reference.
!     ----------
!        Fullpos technical & users guide.

!     Author.
!     -------
!        Ryad El Khatib *METEO-FRANCE*

!     Modifications.
!     --------------
!        Original : 03-02-11 From TRANSDIR_FP
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        R. El Khatib  20-Apr-2007 Optimisation.
!        R. El Khatib 03-Aug-2012 Data arrays in arguments + fpuv2kp (from spos)
!        R. El Khatib 22-Dec-2016 Optimisation
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : MYSETV
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN
USE YOMAFN   , ONLY : ALL_FULLPOS_TYPES

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(IN) :: YDQTYPE
TYPE(ALL_FULLPOS_TYPES), INTENT(IN) :: YDTFP
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KRESOL
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSCA(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVOR(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIV(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)    :: PUMEAN(:) 
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)    :: PVMEAN(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)             :: PSPBFP(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT)    :: PUMEANFP(:) 
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT)    :: PVMEANFP(:) 

INTEGER(KIND=JPIM) :: J, JF, INC, JLEV, IJ, JS, ISPEC2, IFPSCA ,IFPVEC ,IFPUVMN, IFPSPB
INTEGER(KIND=JPIM) :: IFLD(SIZE(PSPBFP,1)), IARR(SIZE(PSPBFP,1))

INTEGER(KIND=JPIM) :: IVOR     ! last pointer of vorticity fields in PVOR
INTEGER(KIND=JPIM) :: IDIV     ! last pointer of divergence fields in PDIV
INTEGER(KIND=JPIM) :: ISCA     ! last pointer of scalar fields in PSCA
INTEGER(KIND=JPIM) :: ISPB     ! last pointer of divergence fields in PSPBFP
INTEGER(KIND=JPIM) :: IMNU     ! last pointer of mean U
INTEGER(KIND=JPIM) :: IMNV     ! last pointer of mean V

LOGICAL :: LLETRANS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

#include "updtrans.intfb.h"
#include "fpuv2kp.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPTSA_DIR',0,ZHOOK_HANDLE)

!      1. Select model

CALL UPDTRANS(KRESOL,LLETRANS)
IF (LLETRANS .AND. &
  & (.NOT.(PRESENT(PUMEAN)) .OR. .NOT.(PRESENT(PVMEAN)) &
  & .OR. .NOT.(PRESENT(PUMEANFP)) .OR. .NOT.(PRESENT(PVMEANFP))) &
  & ) THEN
  CALL ABOR1('FPTSA_DIR : MEAN WIND IS MISSING')
ENDIF

!      2. Set dimensions

IFPSPB=SIZE(PSPBFP,1)
IFPSCA=SIZE(PSCA,1)
IFPVEC=SIZE(PVOR,1)
ISPEC2=SIZE(PSPBFP,2)

IF (SIZE(PDIV,1) /= IFPVEC) THEN
  CALL ABOR1('FPTSA_DIR : INTERNAL ERROR ON SIZE(PVOR,PDIV)') 
ENDIF
IF (LLETRANS) THEN
  IF (SIZE(PUMEAN) /= IFPVEC) THEN
    CALL ABOR1('FPTSA_DIR : INTERNAL ERROR ON SIZE(PVOR,PUMEAN)') 
  ENDIF
  IF (SIZE(PVMEAN) /= IFPVEC) THEN
    CALL ABOR1('FPTSA_DIR : INTERNAL ERROR ON SIZE(PVOR,PVMEAN)') 
  ENDIF
  IFPUVMN=SIZE(PUMEANFP)
  IF (SIZE(PVMEANFP) /= IFPUVMN) THEN
    CALL ABOR1('FPTSA_DIR : INTERNAL ERROR ON SIZE(PUMEANFP,PVMEANFP)') 
  ENDIF
ENDIF
IF (ISPEC2 > 0) THEN
  IF (SIZE(PSCA,2) /= ISPEC2) THEN
    CALL ABOR1('FPTSA_DIR : INTERNAL ERROR ON SIZE(PSPBFP,PSCA)') 
  ENDIF
  IF (SIZE(PVOR,2) /= ISPEC2) THEN
    CALL ABOR1('FPTSA_DIR : INTERNAL ERROR ON SIZE(PSPBFP,PVOR)') 
  ENDIF
  IF (SIZE(PDIV,2) /= ISPEC2) THEN
    CALL ABOR1('FPTSA_DIR : INTERNAL ERROR ON SIZE(PSPBFP,PDIV)') 
  ENDIF
ENDIF


!      3. Set fields correspondences

IVOR=1
IDIV=1
ISCA=1
ISPB=1
IMNU=1
IMNV=1
DO J=1,YDQTYPE%NFPOSDYN
  IJ=YDQTYPE%NFPTRDYN(J)
  IF (YDQTYPE%ISF(IJ) >= 1) THEN
    SELECT CASE (YDQTYPE%ILED(IJ))
    CASE (0)
!       Scalar
      DO JLEV=1,YDQTYPE%ILEV(IJ)
        IF (YDQTYPE%ISET(JLEV,IJ)==MYSETV) THEN
          INC=YDQTYPE%ISPD(JLEV,IJ)
          DO JF=1,INC
            IFLD(ISPB)=ISCA
            IARR(ISPB)=0
            ISPB=ISPB+1
          ENDDO
          ISCA=ISCA+1
        ENDIF
      ENDDO
    CASE (1,2)
!       U or Vor
      DO JLEV=1,YDQTYPE%ILEV(IJ)
        IF (YDQTYPE%ISET(JLEV,IJ)==MYSETV) THEN
          INC=YDQTYPE%ISPD(JLEV,IJ)
          DO JF=1,INC
            IFLD(ISPB)=IVOR
            IARR(ISPB)=1
            ISPB=ISPB+1
          ENDDO
!           Mean U wind if U and lam
          IF (LLETRANS .AND. YDQTYPE%ILED(IJ) /= 2) THEN
            PUMEANFP(IMNU:IMNU+INC-1)=PUMEAN(IVOR)
            IMNU=IMNU+INC
          ENDIF
          IVOR=IVOR+1
!           Skip divergence if Vor :
          IF (YDQTYPE%ILED(IJ) == 2) IDIV=IDIV+1
        ENDIF
      ENDDO
    CASE (:-1,3)
!       V or Div
      DO JLEV=1,YDQTYPE%ILEV(IJ)
        IF (YDQTYPE%ISET(JLEV,IJ)==MYSETV) THEN
          INC=YDQTYPE%ISPD(JLEV,IJ)
          DO JF=1,INC
            IFLD(ISPB)=IDIV
            IARR(ISPB)=-1
            ISPB=ISPB+1
          ENDDO 
!           Mean V wind if V and lam
          IF (LLETRANS .AND. YDQTYPE%ILED(IJ) /= 3) THEN
            PVMEANFP(IMNV:IMNV+INC-1)=PVMEAN(IDIV)
            IMNV=IMNV+INC 
          ENDIF
          IDIV=IDIV+1
!           Skip vorticity if Div :
          IF (YDQTYPE%ILED(IJ) == 3) IVOR=IVOR+1
        ENDIF
      ENDDO
    END SELECT
  ENDIF
ENDDO
IF (ISCA-1 /= IFPSCA) CALL ABOR1('FPTSA_DIR : INTERNAL ERROR ISCA')
IF (IDIV+IVOR-2 /= 2*IFPVEC) CALL ABOR1('FPTSA_DIR : INTERNAL ERROR IDIV+IVOR')
IF (LLETRANS) THEN
  IF (IMNU+IMNV-2 /= 2*IFPUVMN) CALL ABOR1 &
   & ('FPTSA_DIR : INTERNAL ERROR IMNU/IMNV')  
ENDIF

!      4. Fill spectral post-processing arrays

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(J,JS)
DO JS=1,ISPEC2
!DIR$ NOINTERCHANGE
  DO J=1,IFPSPB
    IF (IARR(J)==0) THEN
      PSPBFP(J,JS)=PSCA(IFLD(J),JS)
    ELSEIF (IARR(J)==1) THEN
      PSPBFP(J,JS)=PVOR(IFLD(J),JS)
    ELSEIF (IARR(J)==-1) THEN
      PSPBFP(J,JS)=PDIV(IFLD(J),JS)
    ENDIF 
  ENDDO
ENDDO
!$OMP END PARALLEL DO

!      5. Change U & V to khi & Psi

CALL FPUV2KP(YDQTYPE,YDTFP,KRESOL,PSPBFP)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPTSA_DIR',1,ZHOOK_HANDLE)

END SUBROUTINE FPTSA_DIR
