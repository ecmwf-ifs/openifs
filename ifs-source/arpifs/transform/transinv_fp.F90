! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE TRANSINV_FP(YDQTYPE,KRESOL,PSCA,PVOR,PDIV,PUMEAN,PVMEAN,PGP)

!**** *TRANSINV_FP * - Inverse transforms for model

!     Purpose.  Perform inverse transform (spectral to gridpoint)
!     --------

!**   Interface.  CALL TRANSINV_FP(...)
!     ---------- 

!     Explicit arguments : 
!     --------------------
!        KRESOL : spectral transforms resolution indicator
!        PSCA   : input scalar fields
!        PVOR   : input vorticity fields
!        PDIV   : input divergence fields
!        PUMEAN : input U - mean wind (for LAM)
!        PVMEAN : input V - mean wind (for LAM)
!        PGP    : output gridpoint fields

!        Called by TRANSINVH

!     Externals.
!     ----------
!     INV_TRANS - inverse transform (TRANS library)

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-10-25
!    R. El Khatib : 01-03-07 Fix
!    R. El Khatib : 01-03-28 optimisation
!    R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!    R. El Khatib : 03-04-17 Fullpos improvments
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!    R. El Khatib : 03-Aug-2012 Data arrays in arguments + merge with etransdir_fp
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : NPROC
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(IN)    :: YDQTYPE
INTEGER(KIND=JPIM),INTENT(IN)       :: KRESOL
REAL(KIND=JPRB),INTENT(IN)          :: PSCA(:,:) 
REAL(KIND=JPRB),INTENT(IN)          :: PVOR(:,:) 
REAL(KIND=JPRB),INTENT(IN)          :: PDIV(:,:) 
REAL(KIND=JPRB),INTENT(IN),OPTIONAL :: PUMEAN(:) 
REAL(KIND=JPRB),INTENT(IN),OPTIONAL :: PVMEAN(:)
REAL(KIND=JPRB),INTENT(OUT)         :: PGP(:,:,:)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IVSETSC(YDQTYPE%NFPISCAG), IVSETVD(YDQTYPE%NFPIVECG)
INTEGER(KIND=JPIM) :: IFPSCA, IFPVEC, IPROMA, ISPEC2, IPRTRW, IFPTRV
LOGICAL :: LLETRANS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "inv_trans.h"
#include "einv_trans.h"
#include "trans_inq.h"
#include "etrans_inq.h"

#include "sufpvset_inv.intfb.h"
#include "updtrans.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TRANSINV_FP',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!    1. Setups

CALL UPDTRANS(KRESOL,LLETRANS)
IF (LLETRANS) THEN
  CALL ETRANS_INQ(KRESOL=KRESOL,KSPEC2=ISPEC2,KPRTRW=IPRTRW)
ELSE
  CALL TRANS_INQ(KRESOL=KRESOL,KSPEC2=ISPEC2,KPRTRW=IPRTRW)
ENDIF
IF (LLETRANS .AND. YDQTYPE%NFPIVECG > 0) THEN
  IF (.NOT.(PRESENT(PUMEAN)) .AND. .NOT.(PRESENT(PVMEAN))) THEN
    CALL ABOR1('TRANSINV_FP : MEAN WIND IS MISSING')
  ENDIF
ENDIF
IFPSCA=SIZE(PSCA,1)
IFPVEC=SIZE(PVOR,1)
IPROMA=(SIZE(PGP,1))

IF (SIZE(PDIV,1) /= IFPVEC) THEN
  CALL ABOR1('TRANSINV_FP : INTERNAL ERROR ON SIZE(PVOR,PDIV)') 
ENDIF
IF (LLETRANS) THEN
  IF (SIZE(PUMEAN) /= IFPVEC) THEN
    CALL ABOR1('TRANSINV_FP : INTERNAL ERROR ON SIZE(PVOR,PUMEAN)') 
  ENDIF
  IF (SIZE(PVMEAN) /= IFPVEC) THEN
    CALL ABOR1('TRANSINV_FP : INTERNAL ERROR ON SIZE(PVOR,PVMEAN)') 
  ENDIF
ENDIF
IF (ISPEC2 > 0) THEN
  IF (SIZE(PSCA,2) /= ISPEC2) THEN
    CALL ABOR1('TRANSINV_FP : INTERNAL ERROR ON SIZE(ISPEC2,PSCA)') 
  ENDIF
  IF (SIZE(PVOR,2) /= ISPEC2) THEN
    CALL ABOR1('TRANSINV_FP : INTERNAL ERROR ON SIZE(ISPEC2,PVOR)') 
  ENDIF
  IF (SIZE(PDIV,2) /= ISPEC2) THEN
    CALL ABOR1('TRANSINV_FP : INTERNAL ERROR ON SIZE(ISPEC2,PDIV)') 
  ENDIF
ENDIF
IFPTRV=NPROC/IPRTRW
CALL SUFPVSET_INV(YDQTYPE,IFPTRV,IVSETSC,IVSETVD)

!    2. Spectral transforms

IF (YDQTYPE%NFPISCAG > 0 .AND. YDQTYPE%NFPIVECG > 0) THEN
  IF (LLETRANS) THEN
    CALL EINV_TRANS(PSPVOR=PVOR(:,:),PSPDIV=PDIV(:,:),PSPSCALAR=PSCA(:,:),&
     & PMEANU=PUMEAN(:),PMEANV=PVMEAN(:),&
     & KRESOL=KRESOL,&
     & KPROMA=IPROMA,KVSETUV=IVSETVD(:),KVSETSC=IVSETSC(:),PGP=PGP(:,:,:))
  ELSE
    CALL INV_TRANS(PSPVOR=PVOR(:,:),PSPDIV=PDIV(:,:),PSPSCALAR=PSCA(:,:),&
     & KRESOL=KRESOL,&
     & KPROMA=IPROMA,KVSETUV=IVSETVD(:),KVSETSC=IVSETSC(:),PGP=PGP(:,:,:))  
  ENDIF
ELSEIF (YDQTYPE%NFPIVECG > 0 .AND. YDQTYPE%NFPISCAG == 0) THEN
  IF (LLETRANS) THEN
    CALL EINV_TRANS(PSPVOR=PVOR(:,:),PSPDIV=PDIV(:,:),KPROMA=IPROMA,&
     & PMEANU=PUMEAN(:),PMEANV=PVMEAN(:),&
     & KRESOL=KRESOL,KVSETUV=IVSETVD(:),PGP=PGP(:,:,:))
  ELSE
    CALL INV_TRANS(PSPVOR=PVOR(:,:),PSPDIV=PDIV(:,:),KPROMA=IPROMA,&
     & KRESOL=KRESOL,KVSETUV=IVSETVD(:),PGP=PGP(:,:,:))  
  ENDIF
ELSEIF (YDQTYPE%NFPISCAG > 0 .AND. YDQTYPE%NFPIVECG == 0) THEN
  IF (LLETRANS) THEN
    CALL EINV_TRANS(PSPSCALAR=PSCA(:,:),KPROMA=IPROMA,KVSETSC=IVSETSC(:),&
     & KRESOL=KRESOL,PGP=PGP(:,:,:))  
  ELSE
    CALL INV_TRANS(PSPSCALAR=PSCA(:,:),KPROMA=IPROMA,KVSETSC=IVSETSC(:),&
     & KRESOL=KRESOL,PGP=PGP(:,:,:))  
  ENDIF
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TRANSINV_FP',1,ZHOOK_HANDLE)

END SUBROUTINE TRANSINV_FP
