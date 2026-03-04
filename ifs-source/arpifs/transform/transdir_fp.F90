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

SUBROUTINE TRANSDIR_FP(YDQTYPE,KRESOL,PGP,PSCA,PVOR,PDIV,PUMEAN,PVMEAN)

!**** *TRANSDIR_FP * - Direct transforms for Full-Pos

!     Purpose.  Perform direct transform (gridpoint to spectral)
!     --------

!**   Interface.  CALL TRANSDIR_FP(...)
!     ---------- 

!     Explicit arguments : 
!     --------------------
!        KRESOL : spectral transforms resolution indicator
!        PGP    : input gridpoint fields
!        PSCA   : output scalar fields
!        PVOR   : output vorticity fields
!        PDIV   : output divergence fields
!        PUMEAN : output U - mean wind (for LAM)
!        PVMEAN : output V - mean wind (for LAM)

!        Called by TRANSDIRH

!     Externals.
!     ----------
!     DIR_TRANS - inverse transform (TRANS library)

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-10-25
!    R. El Khatib : 01-03-28 More appropriate ordering
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

TYPE (TYPE_FPRQDYN),  INTENT(IN)     :: YDQTYPE
INTEGER(KIND=JPIM),INTENT(IN)        :: KRESOL
REAL(KIND=JPRB),INTENT(IN)           :: PGP(:,:,:)
REAL(KIND=JPRB),INTENT(OUT)          :: PSCA(:,:) 
REAL(KIND=JPRB),INTENT(OUT)          :: PVOR(:,:) 
REAL(KIND=JPRB),INTENT(OUT)          :: PDIV(:,:) 
REAL(KIND=JPRB),INTENT(OUT),OPTIONAL :: PUMEAN(:) 
REAL(KIND=JPRB),INTENT(OUT),OPTIONAL :: PVMEAN(:)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IVSETSC(YDQTYPE%NFPSCAG), IVSETVD(YDQTYPE%NFPVECG)
INTEGER(KIND=JPIM) :: IFPSCA, IFPVEC, IPROMA, ISPEC2, IPRTRW, IFPTRV
LOGICAL :: LLETRANS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "dir_trans.h"
#include "edir_trans.h"
#include "trans_inq.h"
#include "etrans_inq.h"

#include "sufpvset_dir.intfb.h"
#include "updtrans.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TRANSDIR_FP',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!      1. Set dimensions

CALL UPDTRANS(KRESOL,LLETRANS)
IF (LLETRANS) THEN
  CALL ETRANS_INQ(KRESOL=KRESOL,KSPEC2=ISPEC2,KPRTRW=IPRTRW)
ELSE
  CALL TRANS_INQ(KRESOL=KRESOL,KSPEC2=ISPEC2,KPRTRW=IPRTRW)
ENDIF
IF (LLETRANS .AND. YDQTYPE%NFPVECG > 0) THEN
  IF (.NOT.(PRESENT(PUMEAN)) .AND. .NOT.(PRESENT(PVMEAN))) THEN
    CALL ABOR1('TRANSDIR_FP : MEAN WIND IS MISSING')
  ENDIF
ENDIF
IFPSCA=SIZE(PSCA,1)
IFPVEC=SIZE(PVOR,1)
IPROMA=(SIZE(PGP,1))

IF (SIZE(PDIV,1) /= IFPVEC) THEN
  CALL ABOR1('TRANSDIR_FP : INTERNAL ERROR ON SIZE(PVOR,PDIV)') 
ENDIF
IF (LLETRANS) THEN
  IF (SIZE(PUMEAN) /= IFPVEC) THEN
    CALL ABOR1('TRANSDIR_FP : INTERNAL ERROR ON SIZE(PVOR,PUMEAN)') 
  ENDIF
  IF (SIZE(PVMEAN) /= IFPVEC) THEN
    CALL ABOR1('TRANSDIR_FP : INTERNAL ERROR ON SIZE(PVOR,PVMEAN)') 
  ENDIF
ENDIF
IF (ISPEC2 > 0) THEN
  IF (SIZE(PSCA,2) /= ISPEC2) THEN
    CALL ABOR1('TRANSDIR_FP : INTERNAL ERROR ON SIZE(ISPEC2,PSCA)') 
  ENDIF
  IF (SIZE(PVOR,2) /= ISPEC2) THEN
    CALL ABOR1('TRANSDIR_FP : INTERNAL ERROR ON SIZE(ISPEC2,PVOR)') 
  ENDIF
  IF (SIZE(PDIV,2) /= ISPEC2) THEN
    CALL ABOR1('TRANSDIR_FP : INTERNAL ERROR ON SIZE(ISPEC2,PDIV)') 
  ENDIF
ENDIF

IFPTRV=NPROC/IPRTRW
CALL SUFPVSET_DIR(YDQTYPE,IFPTRV,IVSETSC,IVSETVD)

!    2. Spectral transforms

IF (YDQTYPE%NFPSCAG > 0 .AND. YDQTYPE%NFPVECG > 0) THEN
  IF (LLETRANS) THEN
    CALL EDIR_TRANS(PSPVOR=PVOR(:,:),PSPDIV=PDIV(:,:),PSPSCALAR=PSCA(:,:),&
     & PMEANU=PUMEAN(:),PMEANV=PVMEAN(:),&
     & KRESOL=KRESOL,&
     & KPROMA=IPROMA,KVSETUV=IVSETVD(:),KVSETSC=IVSETSC(:),PGP=PGP(:,:,:))
  ELSE
    CALL DIR_TRANS(PSPVOR=PVOR(:,:),PSPDIV=PDIV(:,:),PSPSCALAR=PSCA(:,:),&
     & KRESOL=KRESOL,&
     & KPROMA=IPROMA,KVSETUV=IVSETVD(:),KVSETSC=IVSETSC(:),PGP=PGP(:,:,:))  
  ENDIF
ELSEIF (YDQTYPE%NFPSCAG > 0 ) THEN
  IF (LLETRANS) THEN
    CALL EDIR_TRANS(PSPSCALAR=PSCA(:,:),KPROMA=IPROMA,KVSETSC=IVSETSC(:),&
     & KRESOL=KRESOL,PGP=PGP(:,:,:))  
  ELSE
    CALL DIR_TRANS(PSPSCALAR=PSCA(:,:),KPROMA=IPROMA,KVSETSC=IVSETSC(:),&
     & KRESOL=KRESOL,PGP=PGP(:,:,:))  
  ENDIF
ELSEIF (YDQTYPE%NFPVECG > 0) THEN
  IF (LLETRANS) THEN
    CALL EDIR_TRANS(PSPVOR=PVOR(:,:),PSPDIV=PDIV(:,:),KPROMA=IPROMA,&
     & PMEANU=PUMEAN(:),PMEANV=PVMEAN(:),&  
     & KRESOL=KRESOL, KVSETUV=IVSETVD(:),PGP=PGP(:,:,:))
  ELSE
    CALL DIR_TRANS(PSPVOR=PVOR(:,:),PSPDIV=PDIV(:,:),KPROMA=IPROMA,&
     & KRESOL=KRESOL, KVSETUV=IVSETVD(:),PGP=PGP(:,:,:))
  ENDIF
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TRANSDIR_FP',1,ZHOOK_HANDLE)
END SUBROUTINE TRANSDIR_FP
