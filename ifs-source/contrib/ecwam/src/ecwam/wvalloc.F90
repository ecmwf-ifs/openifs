! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

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

      INTEGER(KIND=JWIM) :: NPR, MAXLEN, ICHNK

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
 
      IF (LHOOK) CALL DR_HOOK('WVALLOC',0,ZHOOK_HANDLE)

!     1.  ALLOCATE NECESSARY ARRAYS
!         -------------------------

      IF (.NOT.ALLOCATED(WVPRPT%WAVNUM))THEN
         CALL WVPRPT%ALLOC(NPROMA_WAM, NFRE, NCHNK)
      ENDIF

      IF (.NOT.ALLOCATED(FF_NOW%UWND)) THEN
         CALL FF_NOW%ALLOC(NPROMA_WAM, NCHNK)
      ENDIF

      IF (.NOT.ALLOCATED(FL1)) THEN
        ALLOCATE(FL1(NPROMA_WAM, NANG, NFRE, NCHNK))
        FL1(:,:,:,:) = 0.0_JWRB
      ENDIF



      IF (.NOT.ALLOCATED(INTFLDS%PHIEPS)) THEN
        CALL INTFLDS%ALLOC(NPROMA_WAM, NCHNK)
        DO ICHNK=1,NCHNK
          INTFLDS%PHIEPS(:, ICHNK)  = 0.0_JWRB
          INTFLDS%PHIAW(:, ICHNK)   = 0.0_JWRB
          INTFLDS%TAUOC(:, ICHNK)   = 0.0_JWRB
          INTFLDS%STRNMS(:, ICHNK)  = 0.0_JWRB
          INTFLDS%ALTWH(:, ICHNK)   = ZMISS
          INTFLDS%CALTWH(:, ICHNK)  = ZMISS
          INTFLDS%RALTCOR(:, ICHNK) = ZMISS
        ENDDO
      ENDIF


      IF (.NOT.ALLOCATED(FF_NEXT%UWND)) THEN
         CALL FF_NEXT%ALLOC(NPROMA_WAM,NCHNK)
      ENDIF


      IF (.NOT.ALLOCATED(WAM2NEMO%NSWH)) THEN
        CALL WAM2NEMO%ALLOC(NPROMA_WAM, NCHNK)
        DO ICHNK=1,NCHNK
           WAM2NEMO%NSWH(:,ICHNK) = 0.0_JWRO
           WAM2NEMO%NMWP(:,ICHNK) = 0.0_JWRO
           WAM2NEMO%NPHIEPS(:,ICHNK) = 0.0_JWRO
           WAM2NEMO%NEMOPHIF(:,ICHNK) = 0.0_JWRO
           WAM2NEMO%NTAUOC(:,ICHNK) = 0.0_JWRO
           WAM2NEMO%NEMOTAUX(:,ICHNK) = 0.0_JWRO
           WAM2NEMO%NEMOTAUY(:,ICHNK) = 0.0_JWRO
           WAM2NEMO%NEMOUSTOKES(:,ICHNK) = 0.0_JWRO
           WAM2NEMO%NEMOVSTOKES(:,ICHNK) = 0.0_JWRO
           WAM2NEMO%NEMOWSWAVE(:,ICHNK) = 0.0_JWRO
           WAM2NEMO%NEMOSTRN(:,ICHNK) = 0.0_JWRO
        ENDDO
      ENDIF

      IF (.NOT.ALLOCATED(NEMO2WAM%NEMOSST)) THEN
        CALL NEMO2WAM%ALLOC(NPROMA_WAM, NCHNK)
        DO ICHNK=1,NCHNK
           NEMO2WAM%NEMOSST(:,ICHNK) = 0.0_JWRO
           NEMO2WAM%NEMOCICOVER(:,ICHNK) = 0.0_JWRO
           NEMO2WAM%NEMOCITHICK(:,ICHNK) = 0.0_JWRO
           NEMO2WAM%NEMOUCUR(:,ICHNK) = 0.0_JWRO
           NEMO2WAM%NEMOVCUR(:,ICHNK) = 0.0_JWRO
        ENDDO
      ENDIF

      IF (LHOOK) CALL DR_HOOK('WVALLOC',1,ZHOOK_HANDLE)
 
      END SUBROUTINE WVALLOC
