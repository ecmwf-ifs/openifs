! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE OUTBLOCK (KIJS, KIJL, MIJ,                 &
 &                   FL1, XLLWS,                      & 
 &                   WAVNUM, CINV, CGROUP,            &
 &                   DEPTH, UCUR, VCUR, IODP,         &
 &                   ALTWH, CALTWH, RALTCOR,          &
 &                   USTOKES, VSTOKES, STRNMS,        &
 &                   TAUXD, TAUYD, TAUOCXD,           &
 &                   TAUOCYD, TAUOC, PHIOCD,          &
 &                   PHIEPS, PHIAW,                   &
 &                   AIRD, WDWAVE, CICOVER,           &
 &                   WSWAVE, WSTAR,                   &
 &                   UFRIC, TAUW,                     &
 &                   Z0M, Z0B, CHRNCK,                &
 &                   CITHICK,                         &
 &                   NEMOSST, NEMOCICOVER,            &
 &                   NEMOCITHICK, NEMOUCUR, NEMOVCUR, &
 &                   BOUT)

! ----------------------------------------------------------------------

!**** *OUTBLOCK* - MODEL OUTPUT

!*    PURPOSE.
!     --------

!       PREPARE OUTPUT INTEGRATED PARAMETERS.


!**   INTERFACE.
!     ----------
!      *CALL*OUTBLOCK (KIJS, KIJL, MIJ,
!     &                FL1, XLLWS,
!     &                WVPRPT,
!     &                WVENVI, FF_NOW, INTFLDS, NEMO2WAM,
!     &                BOUT)
!      *KIJS*    - INDEX OF FIRST LOCAL GRIDPOINT
!      *KIJL*    - INDEX OF LAST LOCAL GRIDPOINT
!      *MIJ*     - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!      *FL1*     - INPUT SPECTRUM.
!      *XLLWS*   - WINDSEA MASK FROM INPUT SOURCE TERM
!      *WVPRPT*  - WAVE PROPERTIES
!      *WVENVI*  - WAVE ENVORONEMENT
!      *FF_NOW*  - FORCING FIELDS
!      *INTFLDS* - INTEGRATED/DERIVED PARAMETERS
!      *NEMO2WAM*- FIELDS FRON OCEAN MODEL to WAM
!      *BOUT*    - OUTPUT PARAMETER BUFFER

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
     &                         INTGT_PARAM_FIELDS, OCEAN2WAVE

      USE YOWCOUT  , ONLY : NTRAIN   ,IPFGTBL  ,LSECONDORDER,           &
     &            NIPRMOUT, ITOBOUT
      USE YOWCOUP  , ONLY : LWNEMOCOUSTRN
      USE YOWFRED  , ONLY : FR, TH , DFIM, DELTH, COSTH, SINTH, XKMSS_CUTOFF


      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZMISS    ,DEG      ,EPSUS    ,EPSU10, G, ZPI
      USE YOWSTAT  , ONLY : IREFRA

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "cal_second_order_spec.intfb.h"
#include "cimsstrn.intfb.h"
#include "ctcor.intfb.h"
#include "femean.intfb.h"
#include "intpol.intfb.h"
#include "kurtosis.intfb.h"
#include "meansqs.intfb.h"
#include "mwp1.intfb.h"
#include "mwp2.intfb.h"
#include "outbeta.intfb.h"
#include "outsetwmask.intfb.h"
#include "dominant_period.intfb.h"
#include "se10mean.intfb.h"
#include "sebtmean.intfb.h"
#include "sepwisw.intfb.h"
#include "sthq.intfb.h"
#include "wdirspread.intfb.h"
#include "weflux.intfb.h"
#include "w_maxh.intfb.h"
#include "halphap.intfb.h"
#include "alphap_tail.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: MIJ
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1, XLLWS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NFRE), INTENT(IN) :: WAVNUM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NFRE), INTENT(IN) :: CINV
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NFRE), INTENT(IN) :: CGROUP

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: UCUR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: VCUR
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: IODP

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: AIRD, WDWAVE, CICOVER, WSWAVE, WSTAR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: UFRIC, TAUW, Z0M, Z0B, CHRNCK, CITHICK
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: ALTWH, CALTWH, RALTCOR, USTOKES, VSTOKES, STRNMS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUOC, PHIOCD
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: PHIEPS, PHIAW
      REAL(KIND=JWRO), DIMENSION(KIJS:KIJL), INTENT(IN) :: NEMOSST, NEMOCICOVER, NEMOCITHICK, NEMOUCUR, NEMOVCUR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NIPRMOUT), INTENT(OUT) :: BOUT


      INTEGER(KIND=JWIM), PARAMETER :: NTEWH=6
      INTEGER(KIND=JWIM) :: IJ, K, M, ITG, IR, ITR, IH
      INTEGER(KIND=JWIM) :: IRA
      
      REAL(KIND=JWRB) :: SIG
      REAL(KIND=JWRB) :: GOZPI 
      REAL(KIND=JWRB) :: XMODEL_CUTOFF
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(0:NTEWH) :: TEWH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EM, FM, DP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: C3, C4, BF, QP, HMAX, TMAX
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CMAX_F, HMAX_N, CMAX_ST, HMAX_ST, PHIST
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ETA_M, R, XNSLC, SIG_TH, EPS, XNU
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: FLD1, FLD2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ESWELL ,FSWELL ,THSWELL, P1SWELL, P2SWELL, SPRDSWELL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ESEA   ,FSEA   ,THWISEA, P1SEA  , P2SEA  , SPRDSEA
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CHARNOCK, BETAHQ, CDATM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: HALP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NTRAIN) :: EMTRAIN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NTRAIN) :: THTRAIN, PMTRAIN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: COSWDIF

!     *FL2ND*  SPECTRUM with second order effect added if LSECONDORDER is true .
!            and in the absolute frame of reference if currents are used 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: FL2ND

      LOGICAL :: LLPEAKF

      DATA TEWH /10._JWRB,12._JWRB,14._JWRB,17._JWRB,21._JWRB,25._JWRB,30._JWRB/

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OUTBLOCK',0,ZHOOK_HANDLE)

!
!*    1. COMPUTE MEAN PARAMETERS.
!        ------------------------

!     PREPARE THE WAVE SPECTRA THAt SHOULD BE USED FOR OUTPUT

      LLPEAKF = .FALSE.

      IRA=1
      SIG = 1._JWRB
      GOZPI=G/ZPI

      IF (IREFRA == 2 .OR. IREFRA == 3) THEN
        CALL INTPOL (KIJS, KIJL, FL1, FL2ND, WAVNUM, UCUR, VCUR, IRA)
      ELSE
        FL2ND(KIJS:KIJL,:,:) = FL1(KIJS:KIJL,:,:)
      ENDIF
      IF (LSECONDORDER) CALL CAL_SECOND_ORDER_SPEC(KIJS, KIJL, FL2ND, WAVNUM, DEPTH, SIG)


!     COMPUTE MEAN PARAMETERS

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          COSWDIF(IJ,K) = COS(TH(K)-WDWAVE(IJ))
        ENDDO
      ENDDO

      CALL FEMEAN (KIJS, KIJL, FL2ND, EM, FM)

      CALL DOMINANT_PERIOD (KIJS, KIJL, FL1, DP)

      CALL KURTOSIS(KIJS, KIJL, FL1,                     &
     &              DEPTH,                               &
     &              C3, C4, BF, QP, HMAX, TMAX,          &
     &              ETA_M, R, XNSLC, SIG_TH, EPS, XNU)

!     WIND/SWELL PARAMETERS
      CALL SEPWISW (KIJS, KIJL, MIJ, FL1, XLLWS, CINV,             &
     &              UFRIC, WSWAVE, WDWAVE, COSWDIF,  &
     &              ESWELL, FSWELL, THSWELL, P1SWELL, P2SWELL, SPRDSWELL, &
     &              ESEA, FSEA, THWISEA, P1SEA, P2SEA, SPRDSEA,           &
     &              EMTRAIN, THTRAIN, PMTRAIN)


!     LOAD THE OUTPUT BUFFER:
      IR=0

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
!       SIGNIFICANT WAVE HEIGHT CONVERSION
        BOUT(KIJS:KIJL,ITOBOUT(IR))=4._JWRB*SQRT(MAX(EM(KIJS:KIJL),0._JWRB))
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        ITG=ITOBOUT(IR)
        CALL STHQ (KIJS, KIJL, FL2ND, BOUT(KIJS,ITG))
!       CONVERT DIRECTIONS TO DEGREES AND METEOROLOGICAL CONVENTION
        BOUT(KIJS:KIJL,ITG)=MOD(DEG*BOUT(KIJS:KIJL,ITG)+180._JWRB,360._JWRB)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
!       CONVERSION TO PERIOD
        DO IJ=KIJS,KIJL
          IF (FM(IJ) > 0._JWRB) THEN
            BOUT(IJ,ITOBOUT(IR))=1._JWRB/FM(IJ)
          ELSE
            BOUT(IJ,ITOBOUT(IR))=ZMISS
          ENDIF
        ENDDO
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=UFRIC(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
!       CONVERT DIRECTIONS TO DEGREES AND METEOROLOGICAL CONVENTION
        BOUT(KIJS:KIJL,ITOBOUT(IR))=MOD(DEG*WDWAVE(KIJS:KIJL)+180._JWRB,360._JWRB)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
!       CONVERSION TO PERIOD
        DO IJ=KIJS,KIJL
          IF (DP(IJ) > 0.0_JWRB) THEN
            BOUT(IJ,ITOBOUT(IR))=DP(IJ)
          ELSE
            BOUT(IJ,ITOBOUT(IR))=ZMISS
          ENDIF
        ENDDO
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
!!      if the numerical computation of TAU and CD changes, a similar
!!      modification has to be put in buildstress where the friction
!!      velocity is determined from U10 and CD.
!!      Because of the limited numerical resolution when encoding in grib, the maximum value for Cd is set to 0.01
!!!!!!!!!!        BOUT(KIJS:KIJL,ITOBOUT(IR))=MIN(MAX(UFRIC(KIJS:KIJL)**2,EPSUS)/MAX(WSWAVE(KIJS:KIJL)**2,EPSU10**2), 0.01_JWRB)
!!! output the drag coeffient that is consistent with the Charnock parameter that is returned to the atmosphere model:
        CALL OUTBETA (KIJS, KIJL,                      &
     &                WSWAVE, UFRIC, Z0M, Z0B, CHRNCK, &
     &                CHARNOCK, BETAHQ, CD=CDATM)

        BOUT(KIJS:KIJL,ITOBOUT(IR))=MIN(CDATM(KIJS:KIJL), 0.01_JWRB)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=TAUW(KIJS:KIJL)/MAX(UFRIC(KIJS:KIJL)**2,EPSUS)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        CALL MEANSQS (XKMSS_CUTOFF, KIJS, KIJL, FL1, WAVNUM, UFRIC, COSWDIF, BOUT(KIJS,ITOBOUT(IR)))
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=WSWAVE(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
!       WINDSEA SIGNIFICANT WAVE HEIGHT CONVERSION
        BOUT(KIJS:KIJL,ITOBOUT(IR))=4._JWRB*SQRT(MAX(ESEA(KIJS:KIJL),0._JWRB))
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
!       TOTAL SWELL SIGNIFICANT WAVE HEIGHT CONVERSION
        BOUT(KIJS:KIJL,ITOBOUT(IR))=4._JWRB*SQRT(MAX(ESWELL(KIJS:KIJL),0._JWRB))
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
!       CONVERT DIRECTIONS TO DEGREES AND METEOROLOGICAL CONVENTION
        BOUT(KIJS:KIJL,ITOBOUT(IR))=MOD(DEG*THWISEA(KIJS:KIJL)+180._JWRB,360._JWRB)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
!       CONVERT DIRECTIONS TO DEGREES AND METEOROLOGICAL CONVENTION
        BOUT(KIJS:KIJL,ITOBOUT(IR))=MOD(DEG*THSWELL(KIJS:KIJL)+180._JWRB,360._JWRB)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
!       CONVERSION TO PERIOD
        DO IJ=KIJS,KIJL
          IF (FSEA(IJ) > 0._JWRB) THEN
            BOUT(IJ,ITOBOUT(IR))=1._JWRB/FSEA(IJ)
          ELSE
            BOUT(IJ,ITOBOUT(IR))=ZMISS
          ENDIF
        ENDDO
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
!       CONVERSION TO PERIOD
        DO IJ=KIJS,KIJL
          IF (FSWELL(IJ) > 0._JWRB) THEN
            BOUT(IJ,ITOBOUT(IR))=1._JWRB/FSWELL(IJ)
          ELSE
            BOUT(IJ,ITOBOUT(IR))=ZMISS
          ENDIF
        ENDDO
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=ALTWH(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=CALTWH(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=RALTCOR(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        CALL MWP1 (KIJS, KIJL, FL2ND, BOUT(KIJS,ITOBOUT(IR)))
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        CALL MWP2 (KIJS, KIJL, FL2ND, BOUT(KIJS,ITOBOUT(IR)))
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        CALL WDIRSPREAD (KIJS, KIJL, FL2ND, EM, LLPEAKF, BOUT(KIJS,ITOBOUT(IR)))
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=P1SEA(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=P1SWELL(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=P2SEA(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=P2SWELL(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=SPRDSEA(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=SPRDSWELL(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=C4(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=BF(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=QP(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=DEPTH(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=HMAX(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=TMAX(KIJS:KIJL)
      ENDIF

!     SURFACE STOKES DRIFT U and V
      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=USTOKES(KIJS:KIJL)
      ENDIF
      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=VSTOKES(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=UCUR(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=VCUR(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=PHIEPS(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=PHIAW(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=TAUOC(KIJS:KIJL)
      ENDIF

      DO ITR=1,NTRAIN
        IR=IR+1
        IF (IPFGTBL(IR) /= 0) THEN
          BOUT(KIJS:KIJL,ITOBOUT(IR))=4._JWRB*SQRT(MAX(EMTRAIN(KIJS:KIJL,ITR),0._JWRB))
        ENDIF

        IR=IR+1
        IF (IPFGTBL(IR) /= 0) THEN
          BOUT(KIJS:KIJL,ITOBOUT(IR))=MOD(DEG*THTRAIN(KIJS:KIJL,ITR)+180._JWRB,360._JWRB)
        ENDIF

        IR=IR+1
        IF (IPFGTBL(IR) /= 0) THEN
          BOUT(KIJS:KIJL,ITOBOUT(IR))=PMTRAIN(KIJS:KIJL,ITR)
        ENDIF
      ENDDO

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        IF (LWNEMOCOUSTRN) THEN
          BOUT(KIJS:KIJL,ITOBOUT(IR))=STRNMS(KIJS:KIJL)
        ELSE
          CALL CIMSSTRN (KIJS, KIJL, FL1, WAVNUM, DEPTH, CITHICK, BOUT(KIJS,ITOBOUT(IR)))
        ENDIF
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        CALL SE10MEAN (KIJS, KIJL, FL2ND, FLD1)
        BOUT(KIJS:KIJL,ITOBOUT(IR))=4._JWRB*SQRT(MAX(FLD1(KIJS:KIJL),0._JWRB))
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=AIRD(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=WSTAR(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=CICOVER(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=CITHICK(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=C3(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=NEMOSST(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=NEMOCICOVER(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=NEMOCITHICK(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=NEMOUCUR(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=NEMOVCUR(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0 .OR. IPFGTBL(IR+1) /= 0) THEN
           CALL WEFLUX (KIJS, KIJL, FL1, CGROUP,    &
     &                  NFRE, NANG, DFIM, DELTH,    &
     &                  COSTH, SINTH,               &
     &                  FLD1, FLD2)
      ENDIF
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=FLD1(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
!       CONVERT DIRECTIONS TO DEGREES AND METEOROLOGICAL CONVENTION
        BOUT(KIJS:KIJL,ITOBOUT(IR))=MOD(DEG*FLD2(KIJS:KIJL)+180._JWRB,360._JWRB)
      ENDIF

      DO IH=1,NTEWH
        IR=IR+1
        IF (IPFGTBL(IR) /= 0) THEN
          CALL SEBTMEAN (KIJS, KIJL, FL2ND, TEWH(IH-1), TEWH(IH), BOUT(KIJS,ITOBOUT(IR)))
!         SIGNIFICANT WAVE HEIGHT CONVERSION
          BOUT(KIJS:KIJL,ITOBOUT(IR))=4._JWRB*SQRT(MAX(BOUT(KIJS:KIJL,ITOBOUT(IR)),0._JWRB))
        ENDIF
      ENDDO

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=ETA_M(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=R(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=XNSLC(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=TAUXD(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=TAUYD(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=TAUOCXD(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=TAUOCYD(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
!      !!! make the energy flux positive
        BOUT(KIJS:KIJL,ITOBOUT(IR))=MAX(-PHIOCD(KIJS:KIJL),0.0_JWRB)
      ENDIF

!!    alternative ways to determine wave height extremes 
      IF (IPFGTBL(IR  ) /= 0 .OR. IPFGTBL(IR+1) /= 0  .OR. &
&         IPFGTBL(IR+2) /= 0 .OR. IPFGTBL(IR+3) /= 0 ) THEN
        CALL W_MAXH (KIJS, KIJL, FL1, DEPTH, WAVNUM,            &
     &               CMAX_F, HMAX_N, CMAX_ST, HMAX_ST, PHIST)
      ENDIF 

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=CMAX_F(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=HMAX_N(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=CMAX_ST(KIJS:KIJL)
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=HMAX_ST(KIJS:KIJL)
      ENDIF
!!


!     COMPUTE OUTPUT EXTRA FIELDS
!     add necessary code to compute the extra output fields
!!!for testing
      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        CALL CTCOR (KIJS, KIJL, FL1, BOUT(KIJS,ITOBOUT(IR)))
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        XMODEL_CUTOFF=(ZPI*FR(NFRE))**2/G
        CALL MEANSQS (XMODEL_CUTOFF, KIJS, KIJL, FL1, WAVNUM, UFRIC, COSWDIF, BOUT(KIJS,ITOBOUT(IR)))
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=0._JWRB
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=0._JWRB
      ENDIF

      IR=IR+1
      IF (IPFGTBL(IR) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(IR))=0._JWRB
      ENDIF



!     APPLY SEA ICE MASK AND SEA MASK IF NECESSARY
      CALL OUTSETWMASK (KIJS, KIJL, IODP(KIJS:KIJL), CICOVER, BOUT)

IF (LHOOK) CALL DR_HOOK('OUTBLOCK',1,ZHOOK_HANDLE)

END SUBROUTINE OUTBLOCK
