! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      SUBROUTINE WDFLUXES (KIJS, KIJL,                 &
     &                     MIJ,                        &
     &                     FL1, XLLWS,                 &
     &                     WVPRPT,                     &
     &                     WVENVI, FF_NOW,             &
     &                     INTFLDS, WAM2NEMO)

! ----------------------------------------------------------------------

!**** *WDFLUXES* - IF NEEDED, IT EVALUATES SINPUT AND SDISSIP
!****              FOR THE CALCULATION OF OUTPUT WAVE DEPENDANT FLUXES.
!****              THIS IS NORMALLY DONE IN *IMPLSCH* BUT AT INITILISATION,
!****              THE SOURCES TERMS ARE NOT YET EVALUATED.
!****              ALSO DETERMINE THE SURFACE STOKES DRIFT ANd STRIN IN THE ICE.

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *WDFLUXES (KIJS, KIJL,
!    &                    MIJ,
!    &                    FL1, XLLWS,
!    &                    WVPRPT,
!    &                    WVENVI, FF_NOW, INTFLDS, WAM2NEMO)

!          *KIJS*   - INDEX OF FIRST GRIDPOINT.
!          *KIJL*   - INDEX OF LAST GRIDPOINT.
!          *FL1*    - SPECTRUM(INPUT).
!          *XLLWS*  - TOTAL WINDSEA MASK FROM INPUT SOURCE TERM
!          *WVPRPT* - WAVE PROPERTIES FIELDS
!          *WVENVI* - WAVE ENVIRONMENT
!          *FF_NOW* - FORCING FIELDS AT CURRENT TIME.
!          *INTFLDS*-  INTEGRATED/DERIVED PARAMETERS
!          *WAM2NEMO* FIELDS PASSED FROM WAM TO NEMO

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
     &                         INTGT_PARAM_FIELDS, WAVE2OCEAN

      USE YOWCOUP  , ONLY : LWFLUX   ,LWVFLX_SNL, LWNEMOCOUSTRN
      USE YOWCOUT  , ONLY : LWFLUXOUT 
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : WSEMEAN_MIN, ROWATERM1

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "femeanws.intfb.h"
#include "fkmean.intfb.h"
#include "sdissip.intfb.h"
#include "sinflx.intfb.h"
#include "snonlin.intfb.h"
#include "stokestrn.intfb.h"
#include "wnfluxes.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM),  DIMENSION(KIJS:KIJL), INTENT(OUT) :: MIJ

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FL1
      TYPE(FREQUENCY), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WVPRPT
      TYPE(ENVIRONMENT), DIMENSION(KIJS:KIJL), INTENT(IN) :: WVENVI
      TYPE(FORCING_FIELDS), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: FF_NOW
      TYPE(INTGT_PARAM_FIELDS), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: INTFLDS
      TYPE(WAVE2OCEAN), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: WAM2NEMO
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, K, M
      INTEGER(KIND=JWIM) :: ICALL, NCALL

      REAL(KIND=JWRB) :: TAU, XN, PHIDIAG, TAUO
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TAUW_LOC  ! TAUW should not be updated do use a local array
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TAUWDIR_LOC  ! TAUW should not be updated do use a local array
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: RAORW
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EMEAN, FMEAN, HALP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EMEANWS, FMEANWS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: F1MEAN, AKMEAN, XKMEAN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: PHIWA
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: COSWDIF, SINWDIF2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: FLM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: RHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: FLD, SL, SPOS

      LOGICAL :: LCFLX
      LOGICAL :: LUPDTUS

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WDFLUXES',0,ZHOOK_HANDLE)

ASSOCIATE(DEPTH => WVENVI%DEPTH, &
 &        INDEP => WVENVI%INDEP, &
 &        WAVNUM => WVPRPT%WAVNUM, &
 &        CINV => WVPRPT%CINV, &
 &        XK2CG => WVPRPT%XK2CG, &
 &        STOKFAC => WVPRPT%STOKFAC, &
 &        WSWAVE => FF_NOW%WSWAVE, &
 &        WDWAVE => FF_NOW%WDWAVE, &
 &        UFRIC => FF_NOW%UFRIC, &
 &        Z0M => FF_NOW%Z0M, &
 &        Z0B => FF_NOW%Z0B, &
 &        CHRNCK => FF_NOW%CHRNCK, &
 &        AIRD => FF_NOW%AIRD, &
 &        WSTAR => FF_NOW%WSTAR, &
 &        CICOVER => FF_NOW%CICOVER, &
 &        WSEMEAN => INTFLDS%WSEMEAN, &
 &        WSFMEAN => INTFLDS%WSFMEAN)

!*    1. INITIALISATION.
!        ---------------

      LCFLX=LWFLUX.OR.LWFLUXOUT
! ----------------------------------------------------------------------

!*    1.2 COMPUTATION OF RELEVANT SOURCE FUNCTIONS.
!         -----------------------------------------

      CALL FKMEAN(KIJS, KIJL, FL1, WAVNUM,                         &
     &            EMEAN, FMEAN, F1MEAN, AKMEAN, XKMEAN)

      TAUW_LOC(:) = 0.0_JWRB
      TAUWDIR_LOC(:) = WDWAVE(:)

      LUPDTUS = .FALSE.
      FMEANWS(:) = FMEAN(:)
      FLM(:,:) = 0.0_JWRB

      DO IJ=KIJS,KIJL
        RAORW(IJ) = MAX(AIRD(IJ), 1.0_JWRB) * ROWATERM1
      ENDDO

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          COSWDIF(IJ,K) = COS(TH(K)-WDWAVE(IJ))
          SINWDIF2(IJ,K) = SIN(TH(K)-WDWAVE(IJ))**2
        ENDDO
      ENDDO

      NCALL = 1
      ICALL = 1
      CALL SINFLX (ICALL, NCALL, KIJS, KIJL,                              &
     &             LUPDTUS,                                               &
     &             FL1,                                                   &
     &             WAVNUM, CINV, XK2CG,                                   &
     &             WSWAVE, WDWAVE, AIRD, RAORW, WSTAR, CICOVER,           &
     &             COSWDIF, SINWDIF2,                                     &
     &             FMEAN, HALP, FMEANWS,                                  &
     &             FLM,                                                   &
     &             UFRIC, TAUW_LOC, TAUWDIR_LOC, Z0M, Z0B, CHRNCK, PHIWA, &
     &             FLD, SL, SPOS,                                         &
     &             MIJ, RHOWGDFTH, XLLWS)

      IF (LCFLX) THEN

        CALL SDISSIP (KIJS, KIJL, FL1 ,FLD, SL,    &
     &                INDEP, WAVNUM, XK2CG,        &
     &                EMEAN, F1MEAN, XKMEAN,       &
     &                UFRIC, COSWDIF, RAORW) 

        IF (.NOT. LWVFLX_SNL) THEN
          CALL WNFLUXES (KIJS, KIJL,                        &
     &                   MIJ, RHOWGDFTH,                    &
     &                   CINV,                              &
     &                   SL, CICOVER,                       &
     &                   PHIWA,                             &
     &                   EMEAN, F1MEAN, WSWAVE, WDWAVE,     &
     &                   UFRIC, AIRD, INTFLDS, WAM2NEMO,    &
     &                  .FALSE.)
        ENDIF

        CALL SNONLIN (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN)

        IF (LWVFLX_SNL) THEN
          CALL WNFLUXES (KIJS, KIJL,                        &
     &                   MIJ, RHOWGDFTH,                    &
     &                   CINV,                              &
     &                   SL, CICOVER,                       &
     &                   PHIWA,                             &
     &                   EMEAN, F1MEAN, WSWAVE, WDWAVE,     &
     &                   UFRIC, AIRD, INTFLDS, WAM2NEMO,    &
     &                  .FALSE.)
        ENDIF

        IF (LWFLUX) THEN
         CALL FEMEANWS(KIJS, KIJL, FL1, XLLWS, EMEANWS, FMEANWS)

          DO IJ=KIJS,KIJL
            IF (EMEANWS(IJ) < WSEMEAN_MIN) THEN
              WSEMEAN(IJ) = WSEMEAN_MIN 
              WSFMEAN(IJ) = 2._JWRB*FR(NFRE)
            ELSE
              WSEMEAN(IJ) = EMEANWS(IJ)
              WSFMEAN(IJ) = FMEANWS(IJ) 
            ENDIF
          ENDDO
        ENDIF

        CALL STOKESTRN(KIJS, KIJL, FL1, WAVNUM, STOKFAC, DEPTH, FF_NOW, INTFLDS, WAM2NEMO)

      ENDIF
! ----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('WDFLUXES',1,ZHOOK_HANDLE)

END SUBROUTINE WDFLUXES

