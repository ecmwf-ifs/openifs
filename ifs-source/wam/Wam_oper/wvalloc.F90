! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      SUBROUTINE WVALLOC

! ----------------------------------------------------------------------

!**** *WVALLOC* - WAVE MODEL ARRAY ALLOCATION

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO

      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWMEAN  , ONLY : INTFLDS
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSHAL  , ONLY : WVPRPT
      USE YOWSPEC  , ONLY : FF_NOW   ,FL1
      USE YOWWIND  , ONLY : FF_NEXT

      USE YOWNEMOFLDS , ONLY : WAM2NEMO, NEMO2WAM

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: NPR, MAXLEN

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
 
      IF (LHOOK) CALL DR_HOOK('WVALLOC',0,ZHOOK_HANDLE)

!     1.  ALLOCATE NECESSARY ARRAYS
!         -------------------------

      IF (.NOT.ALLOCATED(WVPRPT)) ALLOCATE(WVPRPT(NPROMA_WAM, NFRE, NCHNK))

      IF (.NOT.ALLOCATED(FF_NOW)) ALLOCATE(FF_NOW(NPROMA_WAM, NCHNK))

      IF (.NOT.ALLOCATED(FL1)) THEN
        ALLOCATE(FL1(NPROMA_WAM, NANG, NFRE, NCHNK))
        FL1(:,:,:,:) = 0.0_JWRB
      ENDIF



      IF (.NOT.ALLOCATED(INTFLDS)) THEN 
        ALLOCATE(INTFLDS(NPROMA_WAM, NCHNK))
        INTFLDS(:,:)%PHIEPS  = 0.0_JWRB
        INTFLDS(:,:)%PHIAW   = 0.0_JWRB
        INTFLDS(:,:)%TAUOC   = 0.0_JWRB
        INTFLDS(:,:)%STRNMS  = 0.0_JWRB
        INTFLDS(:,:)%ALTWH   = ZMISS 
        INTFLDS(:,:)%CALTWH  = ZMISS 
        INTFLDS(:,:)%RALTCOR = ZMISS 
      ENDIF


      IF (.NOT.ALLOCATED(FF_NEXT)) ALLOCATE(FF_NEXT(NPROMA_WAM, NCHNK))


      IF (.NOT.ALLOCATED(WAM2NEMO)) THEN
        ALLOCATE(WAM2NEMO(NPROMA_WAM, NCHNK))
        WAM2NEMO(:,:)%NSWH = 0.0_JWRO
        WAM2NEMO(:,:)%NMWP = 0.0_JWRO
        WAM2NEMO(:,:)%NPHIEPS = 0.0_JWRO
        WAM2NEMO(:,:)%NEMOPHIF = 0.0_JWRO
        WAM2NEMO(:,:)%NTAUOC = 0.0_JWRO
        WAM2NEMO(:,:)%NEMOTAUX = 0.0_JWRO
        WAM2NEMO(:,:)%NEMOTAUY = 0.0_JWRO
        WAM2NEMO(:,:)%NEMOUSTOKES = 0.0_JWRO
        WAM2NEMO(:,:)%NEMOVSTOKES = 0.0_JWRO
        WAM2NEMO(:,:)%NEMOWSWAVE = 0.0_JWRO
        WAM2NEMO(:,:)%NEMOSTRN = 0.0_JWRO
      ENDIF

      IF (.NOT.ALLOCATED(NEMO2WAM)) THEN
        ALLOCATE(NEMO2WAM(NPROMA_WAM, NCHNK))
        NEMO2WAM(:,:)%NEMOSST = 0.0_JWRO
        NEMO2WAM(:,:)%NEMOCICOVER = 0.0_JWRO
        NEMO2WAM(:,:)%NEMOCITHICK = 0.0_JWRO
        NEMO2WAM(:,:)%NEMOUCUR = 0.0_JWRO
        NEMO2WAM(:,:)%NEMOVCUR = 0.0_JWRO
      ENDIF

      IF (LHOOK) CALL DR_HOOK('WVALLOC',1,ZHOOK_HANDLE)
 
      END SUBROUTINE WVALLOC
