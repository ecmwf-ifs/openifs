! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE OUTSTEP0 (WVENVI, WVPRPT, FF_NOW, INTFLDS,  &
 &                   WAM2NEMO, NEMO2WAM, FL1)

! ----------------------------------------------------------------------

!**** *OUTSTEP0* - OUTPUTS INITIAL CONDITION AND/OR FORECAST STEP 0

!**   INTERFACE.
!     ----------

!     *CALL* *OUTSTEP0 (WVENVI, WVPRPT, FF_NOW, INTFLDS,
!    &                  WAM2NEMO, NEMO2WAM, FL1)
!        *WVENVI*    WAVE ENVIRONMENT FIELDS
!        *WVPRPT*    WAVE PROPERTIES FIELDS
!        *FF_NOW*    FORCING FIELDS AT CURRENT TIME.
!        *INTFLDS*   INTEGRATED/DERIVED PARAMETERS
!        *WAM2NEMO*  WAVE FIELDS PASSED TO NEMO
!        *NEMO2WAM*  FIELDS FRON OCEAN MODEL to WAM

! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
     &                         INTGT_PARAM_FIELDS, WAVE2OCEAN, OCEAN2WAVE

      USE YOWCOUT  , ONLY : JPPFLAG  ,FFLAG    ,GFLAG    ,                        & 
     &                      IRCD     ,IRU10    , IRALTHS  ,IRALTHSC  ,IRALTRC ,   &
     &                      NGOUT    ,NIPRMOUT , LFDB     ,LRSTST0  ,LWAMANOUT
      USE YOWFRED  , ONLY : TH
      USE YOWGRIBHD, ONLY : LGRHDIFS 
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE
      USE YOWMESPAS, ONLY : LGRIBOUT ,LNOCDIN  ,LWAVEWIND 
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,LLUNSTR
      USE YOWSTAT  , ONLY : CDTPRO   ,CDTINTT  ,IREST   , MARSTYPE ,    &
     &                      LLSOURCE , LANAONLY ,LFRSTFLD
      USE YOWTEST  , ONLY : IU06
      USE YOWTEXT  , ONLY : LRESTARTED
#ifdef WAM_HAVE_UNWAM
      USE UNWAM    , ONLY : EXCHANGE_FOR_FL1
#endif

      USE WAM_MULTIO_MOD, ONLY : WAM_MULTIO_FLUSH
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE EC_LUN   , ONLY : NULERR
      
! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "mpcrtbl.intfb.h"
#include "outwint.intfb.h"
#include "outwpsp.intfb.h"
#include "outbs.intfb.h"
#include "outint.intfb.h"
#include "outspec.intfb.h"
#include "setice.intfb.h"
#include "wdfluxes.intfb.h"

      TYPE(ENVIRONMENT), INTENT(INOUT)                                         :: WVENVI
      TYPE(FREQUENCY), INTENT(INOUT)                                           :: WVPRPT
      TYPE(FORCING_FIELDS), INTENT(INOUT)                                      :: FF_NOW
      TYPE(INTGT_PARAM_FIELDS), INTENT(INOUT)                                  :: INTFLDS
      TYPE(WAVE2OCEAN), INTENT(INOUT)                                          :: WAM2NEMO
      TYPE(OCEAN2WAVE), INTENT(IN)                                             :: NEMO2WAM
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1


      INTEGER(KIND=JWIM) :: IJ, K, ICHNK, KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(NPROMA_WAM, NCHNK) :: MIJ

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG) :: COSWDIF
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:,:) :: BOUTST0
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK) :: XLLWS

      CHARACTER(LEN= 2) :: MARSTYPEBAK
      CHARACTER(LEN=14) :: CDTINTTBAK

      LOGICAL :: LLFLUSH
      LOGICAL, DIMENSION(JPPFLAG) :: FFLAGBAK, GFLAGBAK

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OUTSTEP0',0,ZHOOK_HANDLE)


      LLFLUSH = .FALSE.
      LRSTST0=.FALSE.

      CDTINTTBAK=CDTINTT
      CDTINTT=CDTPRO


!     1.0  WAVES FLUXES AT INITIAL TIME
!          ----------------------------
      IF (LLSOURCE) THEN

!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK)
        DO ICHNK = 1, NCHNK
          CALL WDFLUXES (1, NPROMA_WAM,                          &
     &                   MIJ(:,ICHNK),                           &
     &                   FL1(:,:,:,ICHNK), XLLWS(:,:,:,ICHNK),   &
     &                   WVPRPT%WAVNUM(:,:,ICHNK), WVPRPT%CINV(:,:,ICHNK),  WVPRPT%CGROUP(:,:,ICHNK), &
     &                   WVPRPT%XK2CG(:,:,ICHNK), WVPRPT%STOKFAC(:,:,ICHNK), &
     &                   WVENVI%DEPTH(:,ICHNK), &
     &                   FF_NOW%WSWAVE(:,ICHNK), FF_NOW%WDWAVE(:,ICHNK), &
     &                   FF_NOW%AIRD(:,ICHNK), FF_NOW%WSTAR(:,ICHNK), &
     &                   FF_NOW%USTRA(:,ICHNK), FF_NOW%VSTRA(:,ICHNK), &
     &                   FF_NOW%CICOVER(:,ICHNK), &
     &                   FF_NOW%UFRIC(:,ICHNK), FF_NOW%Z0M(:,ICHNK), &
     &                   FF_NOW%Z0B(:,ICHNK), FF_NOW%CHRNCK(:,ICHNK), FF_NOW%CITHICK(:,ICHNK), &
     &                   INTFLDS%WSEMEAN(:,ICHNK), INTFLDS%WSFMEAN(:,ICHNK), &
     &                   INTFLDS%USTOKES(:,ICHNK), INTFLDS%VSTOKES(:,ICHNK), INTFLDS%STRNMS(:,ICHNK), &
     &                   INTFLDS%TAUXD(:,ICHNK), INTFLDS%TAUYD(:,ICHNK), INTFLDS%TAUOCXD(:,ICHNK), &
     &                   INTFLDS%TAUOCYD(:,ICHNK), INTFLDS%TAUOC(:,ICHNK), INTFLDS%PHIOCD(:,ICHNK), &
     &                   INTFLDS%PHIEPS(:,ICHNK), INTFLDS%PHIAW(:,ICHNK), &
     &                   WAM2NEMO%NEMOUSTOKES(:,ICHNK), WAM2NEMO%NEMOVSTOKES(:,ICHNK), WAM2NEMO%NEMOSTRN(:,ICHNK), &
     &                   WAM2NEMO%NPHIEPS(:,ICHNK), WAM2NEMO%NTAUOC(:,ICHNK), WAM2NEMO%NSWH(:,ICHNK), &
     &                   WAM2NEMO%NMWP(:,ICHNK), WAM2NEMO%NEMOTAUX(:,ICHNK), WAM2NEMO%NEMOTAUY(:,ICHNK), &
     &                   WAM2NEMO%NEMOWSWAVE(:,ICHNK), WAM2NEMO%NEMOPHIF(:,ICHNK))
        ENDDO
!$OMP   END PARALLEL DO

        IF (LLUNSTR) THEN
!!!!! this will need to be adapted to use FL1 
!!! is it still needed ?
        WRITE(NULERR,*) '!!! ********************************* !!'
        WRITE(NULERR,*) '!!! in outstep0. Not yet ready !!!' 
        WRITE(NULERR,*) '!!! ********************************* !!'
        CALL ABORT1
!!!!           CALL EXCHANGE_FOR_FL1(FL1)
        ENDIF

      ELSE
        MIJ(:,:) = NFRE
        XLLWS(:,:,:,:) = 0.0_JWRB
      ENDIF


!     1.2 SET SPECTRA ON ICE POINTS TO ZERO
!         ---------------------------------
      IF (LICERUN .AND. LMASKICE .AND. LLSOURCE) THEN

        CALL GSTATS(1439,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK, K, COSWDIF)
        DO ICHNK = 1, NCHNK
          DO K = 1, NANG
            COSWDIF(:,K) = COS(TH(K)-FF_NOW%WDWAVE(:,ICHNK))
          ENDDO

          CALL SETICE(1, NPROMA_WAM, FL1(:,:,:,ICHNK) ,                &
     &                FF_NOW%CICOVER(:,ICHNK), COSWDIF)

        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1439,1)

      ENDIF

!     2.0 OUTPUT POINT SPECTRA (not usually used at ECMWF)
!         -----------------------------------------------
      IF (NGOUT > 0 ) CALL OUTWPSP (FL1, FF_NOW)


!*    3.0 SAVE INITIAL INTEGRATED FIELDS (if needed)
!         ------------------------------
      IF ((MARSTYPE == 'cf' .OR. MARSTYPE == 'pf' .OR.                 &
     &     MARSTYPE == 'fc' .OR. MARSTYPE == '4v' .OR.                 &
     &     LANAONLY         .OR. LFRSTFLD             )                &
     &    .AND. LWAMANOUT) THEN

        IF (LFRSTFLD) THEN
          FFLAGBAK=FFLAG
          GFLAGBAK=GFLAG
          MARSTYPEBAK=MARSTYPE
          IF (MARSTYPE == 'fg'.OR.MARSTYPE == 'an') THEN
            IF (.NOT.LNOCDIN) THEN
              FFLAG(IRCD)=.FALSE.
              GFLAG(IRCD)=.FALSE.
            ENDIF
            IF (LWAVEWIND) THEN
              FFLAG(IRU10)=.FALSE.
              GFLAG(IRU10)=.FALSE.
            ENDIF
            GFLAG(IRALTHS)=.FALSE.
            GFLAG(IRALTHSC)=.FALSE.
            GFLAG(IRALTRC)=.FALSE.
            MARSTYPE='an'
          ENDIF

!         change the output flags need to reset mapping
          CALL MPCRTBL

        ENDIF

!       NEED TO TEMPORARLY RESET THE IFS FORECAST STEP
!       IF THE IFS GRIB HEADER IS USED SUCH THAT IT POINTS TO
!       THE START OF THE RUN.
        IF (LGRHDIFS) LRSTST0=.TRUE.

!       COMPUTE OUTPUT PARAMETERS AND PRINT OUT NORMS
        IF (NIPRMOUT > 0) THEN
          ALLOCATE(BOUTST0(NPROMA_WAM, NIPRMOUT, NCHNK))
          CALL OUTBS (MIJ, FL1, XLLWS,                             &
     &                WVPRPT, WVENVI, FF_NOW, INTFLDS, NEMO2WAM,   &
     &                BOUTST0)
        ENDIF

        IF ( .NOT. LRESTARTED ) THEN
          IF (IREST == 1 .AND. MARSTYPE /= 'an' .AND. LGRIBOUT) THEN
            CALL OUTSPEC(FL1, FF_NOW)
            WRITE(IU06,*) ''
            WRITE(IU06,*) ' OUTSTEP0: OUTSPEC WAS CALLED !'
            LLFLUSH = .TRUE.
          ENDIF

          IF (NIPRMOUT > 0 ) THEN
            CALL OUTWINT(BOUTST0)
            WRITE(IU06,*) ''
            WRITE(IU06,*) ' OUTSTEP0: OUTWINT WAS CALLED !'
            LLFLUSH = .TRUE.
          ENDIF

        ENDIF

        IF (ALLOCATED(BOUTST0)) DEALLOCATE(BOUTST0)

        IF (LFDB .AND. LLFLUSH) THEN
           CALL GSTATS(1976,0)
           CALL WAM_MULTIO_FLUSH()
           CALL GSTATS(1976,1)
           LLFLUSH=.FALSE.
        ENDIF

        IF (LFRSTFLD) THEN
          MARSTYPE=MARSTYPEBAK
          FFLAG=FFLAGBAK
          GFLAG=GFLAGBAK
          LFRSTFLD=.FALSE.
!         change the output flags need to reset mapping
          CALL MPCRTBL
        ENDIF

        CDTINTT=CDTINTTBAK

        IF (LGRHDIFS) LRSTST0=.FALSE.

      ENDIF


IF (LHOOK) CALL DR_HOOK('OUTSTEP0',1,ZHOOK_HANDLE)

END SUBROUTINE OUTSTEP0
