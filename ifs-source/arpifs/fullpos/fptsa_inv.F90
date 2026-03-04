! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPTSA_INV(YDQTYPE,KRESOL,PSPBFP,PUMEANFP,PVMEANFP,PSCA,PVOR,PDIV,PUMEAN,PVMEAN)

!**** *FPTSA_INV * - Full-Pos Transfer Spectral Array - Inverse

!     Purpose
!     -------
!       To transfer the post-processing spectral array into the spectral
!       transforms arrays.

!**   Interface.  CALL FPTSA_INV(...)
!     ---------- 

!     Explicit arguments : 
!     --------------------
!        KRESOL : Resolution indicator
!        PSPBFP   : input spectral fields
!        PUMEANFP : input U - mean wind (for LAM)
!        PVMEANFP : input V - mean wind (for LAM)
!        PSCA   : output scalar fields
!        PVOR   : output vorticity fields
!        PDIV   : output divergence fields
!        PUMEAN : output U - mean wind (for LAM)
!        PVMEAN : output V - mean wind (for LAM)

!     Externals.
!     ----------
!        None.

!     Reference.
!     ----------
!        Fullpos technical & users guide.

!     Author.
!     -------
!        Ryad El Khatib *METEO-FRANCE*

!     Modifications.
!     --------------
!        Original : 03-02-11 From TRANSINV_FP
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        R. El Khatib  20-Apr-2007 Optimisation.
!        R. El Khatib 03-Aug-2012 Data arrays in arguments
!        R. El Khatib 22-Dec-2016 Optimisation
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : MYSETV
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(IN) :: YDQTYPE
INTEGER(KIND=JPIM) ,INTENT(IN)          :: KRESOL
REAL(KIND=JPRB) ,INTENT(IN)             :: PSPBFP(:,:) 
REAL(KIND=JPRB) ,OPTIONAL,INTENT(IN)    :: PUMEANFP(:) 
REAL(KIND=JPRB) ,OPTIONAL,INTENT(IN)    :: PVMEANFP(:)
REAL(KIND=JPRB) ,INTENT(OUT)   :: PSCA(:,:) 
REAL(KIND=JPRB) ,INTENT(OUT)   :: PVOR(:,:) 
REAL(KIND=JPRB) ,INTENT(OUT)   :: PDIV(:,:) 
REAL(KIND=JPRB) ,OPTIONAL,INTENT(OUT) :: PUMEAN(:) 
REAL(KIND=JPRB) ,OPTIONAL,INTENT(OUT) :: PVMEAN(:) 
INTEGER(KIND=JPIM) :: J, INC, IJ, JLEV, JS, JF, IFPSPB, ISPEC2, IFPSCA ,IFPVEC, IFPUVMN
INTEGER(KIND=JPIM) :: IFLDVOR(SIZE(PVOR,1)), IFLDDIV(SIZE(PVOR,1))
INTEGER(KIND=JPIM) :: IARR(SIZE(PSPBFP,1)), IFLDSCA(SIZE(PSCA,1))

INTEGER(KIND=JPIM) :: IVOR     ! last pointer of vorticity fields in PVOR
INTEGER(KIND=JPIM) :: IDIV     ! last pointer of divergence fields in PDIV
INTEGER(KIND=JPIM) :: ISCA     ! last pointer of scalar fields in PSCA
INTEGER(KIND=JPIM) :: ISPB     ! last pointer of scalar fields in PSPBFP
INTEGER(KIND=JPIM) :: IMNU     ! last pointer of mean U
INTEGER(KIND=JPIM) :: IMNV     ! last pointer of mean V

LOGICAL :: LLETRANS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
#include "updtrans.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPTSA_INV',0,ZHOOK_HANDLE)

!      1. Select model

CALL UPDTRANS(KRESOL,LLETRANS)
IF (LLETRANS .AND. &
  & (.NOT.(PRESENT(PUMEAN)) .OR. .NOT.(PRESENT(PVMEAN)) &
  & .OR. .NOT.(PRESENT(PUMEANFP)) .OR. .NOT.(PRESENT(PVMEANFP))) &
  & ) THEN
  CALL ABOR1('FPTSA_INV : MEAN WIND IS MISSING')
ENDIF

!      2. Set dimensions

IFPSPB=SIZE(PSPBFP,1)
IFPSCA=SIZE(PSCA,1)
IFPVEC=SIZE(PVOR,1)
ISPEC2=SIZE(PSPBFP,2)

IF (SIZE(PDIV,1) /= IFPVEC) THEN
  CALL ABOR1('FPTSA_INV : INTERNAL ERROR ON SIZE(PVOR,PDIV)') 
ENDIF
IF (LLETRANS) THEN
  IF (SIZE(PUMEAN) /= IFPVEC) THEN
    CALL ABOR1('FPTSA_INV : INTERNAL ERROR ON SIZE(PVOR,PUMEAN)') 
  ENDIF
  IF (SIZE(PVMEAN) /= IFPVEC) THEN
    CALL ABOR1('FPTSA_INV : INTERNAL ERROR ON SIZE(PVOR,PVMEAN)') 
  ENDIF
  IFPUVMN=SIZE(PUMEANFP)
  IF (SIZE(PVMEANFP) /= IFPUVMN) THEN
    CALL ABOR1('FPTSA_INV : INTERNAL ERROR ON SIZE(PUMEANFP,PVMEANFP)') 
  ENDIF
ENDIF
IF (ISPEC2 > 0) THEN
  IF (SIZE(PSCA,2) /= ISPEC2) THEN
    CALL ABOR1('FPTSA_INV : INTERNAL ERROR ON SIZE(PSPBFP,PSCA)') 
  ENDIF
  IF (SIZE(PVOR,2) /= ISPEC2) THEN
    CALL ABOR1('FPTSA_INV : INTERNAL ERROR ON SIZE(PSPBFP,PVOR)') 
  ENDIF
  IF (SIZE(PDIV,2) /= ISPEC2) THEN
    CALL ABOR1('FPTSA_INV : INTERNAL ERROR ON SIZE(PSPBFP,PDIV)') 
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
    CASE (1)
!       U  => put Vor
      DO JLEV=1, YDQTYPE%ILEV(IJ)
        IF (YDQTYPE%ISET(JLEV,IJ)==MYSETV) THEN
          INC=YDQTYPE%ISPD(JLEV,IJ)
          DO JF=0,INC-1
            IFLDVOR(IVOR+JF)=ISPB+JF
            IARR(ISPB+JF)=1
          ENDDO
          ISPB=ISPB+INC
          IF (LLETRANS) THEN
            PUMEAN(IVOR:IVOR+INC-1)=PUMEANFP(IMNU:IMNU+INC-1)
            IMNU=IMNU+INC
          ENDIF
          IVOR=IVOR+INC
        ENDIF
      ENDDO
    CASE (:-1)
!       V => put Div
      DO JLEV=1, YDQTYPE%ILEV(IJ)
        IF (YDQTYPE%ISET(JLEV,IJ)==MYSETV) THEN
          INC=YDQTYPE%ISPD(JLEV,IJ)
          DO JF=0,INC-1
            IFLDDIV(IDIV+JF)=ISPB+JF
            IARR(ISPB+JF)=-1
          ENDDO
          ISPB=ISPB+INC
          IF (LLETRANS) THEN
            PVMEAN(IDIV:IDIV+INC-1)=PVMEANFP(IMNV:IMNV+INC-1)
            IMNV=IMNV+INC
          ENDIF
          IDIV=IDIV+INC
        ENDIF
      ENDDO
    CASE (0,2,3)
!       "Scalar"
      DO JLEV=1, YDQTYPE%ILEV(IJ)
        IF (YDQTYPE%ISET(JLEV,IJ)==MYSETV) THEN
          INC=YDQTYPE%ISPD(JLEV,IJ)
          DO JF=0,INC-1
            IFLDSCA(ISCA+JF)=ISPB+JF
            IARR(ISPB+JF)=0
          ENDDO 
          ISPB=ISPB+INC
          ISCA=ISCA+INC
        ENDIF
      ENDDO
    END SELECT
  ENDIF
ENDDO
IF (ISCA-1 /= IFPSCA) CALL ABOR1('FPTSA_INV : INTERNAL ERROR ISCA')
IF (IDIV+IVOR-2 /= 2*IFPVEC) CALL ABOR1('FPTSA_INV : INTERNAL ERROR IDIV+IVOR')
IF (ISPB-1 /= IFPSCA+2*IFPVEC) CALL ABOR1('FPTSA_INV : INTERNAL ERROR ISPB')
IF (LLETRANS) THEN
  IF (IMNU+IMNV-2 /= 2*IFPUVMN) CALL ABOR1 &
   & ('FPTSA_INV : INTERNAL ERROR IMNU/IMNV')  
ENDIF

!      4. Fill spectral post-processing arrays

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JF,JS)
DO JS=1,ISPEC2
!DIR$ NOINTERCHANGE
  DO JF=1,IFPVEC
    IF (IARR(IFLDVOR(JF))==1) THEN
      PVOR(JF,JS)=PSPBFP(IFLDVOR(JF),JS)
    ENDIF
  ENDDO
!DIR$ NOINTERCHANGE
  DO JF=1,IFPVEC
    IF (IARR(IFLDDIV(JF))==-1) THEN
      PDIV(JF,JS)=PSPBFP(IFLDDIV(JF),JS)
    ENDIF
  ENDDO
!DIR$ NOINTERCHANGE
  DO JF=1,IFPSCA
    IF (IARR(IFLDSCA(JF))==0) THEN
      PSCA(JF,JS)=PSPBFP(IFLDSCA(JF),JS)
    ENDIF
  ENDDO
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPTSA_INV',1,ZHOOK_HANDLE)

END SUBROUTINE FPTSA_INV
