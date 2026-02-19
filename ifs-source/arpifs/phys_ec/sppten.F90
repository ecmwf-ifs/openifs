! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SPPTEN (YDMODEL,YGFL,KIDIA,KFDIA,KLON,KLEV,KSDTN2D,PTSPHY,&
 & PTSL, PQSL, PA, PAP, PAPH, &
 & PDYN_U, PDYN_V, PDYN_T, PDYN_Q, &
 & PUNP_U ,PUNP_V ,PUNP_T ,PUNP_Q, &
 & PPHY_U ,PPHY_V ,PPHY_T ,PPHY_Q, &
 & PMULNOISE, &
 & PTENU  ,PTENV  ,PTENT  ,PTENQ, PTENGFL, KTRAC, PTENC, PPERT)

!**** *SPPTEN * - Add stochastically perturbed parametrization tendencies

!     PURPOSE.
!     --------
!       Add stochastically perturbed parametrization tendencies (SPPT)

!**   INTERFACE.
!     ----------
!        *CALL* *SPPTEN(...)*

!     INPUT ARGUMENTS.

!     KIDIA   : start of horizontal loop
!     KFDIA   : end   of horizontal loop
!     KLON    : horizontal dimension
!     KSDTN2D : number of independent (2D) patterns
!     KLEV    : end of vertical loop and vertical dimension
!     PMULNOISE: multiplicative noise (minus 1) applied to U, V, T, Q
!     PDYN_x  : U,V,T,Q tendencies: dynamics ONLY == NOT to be perturbed
!     PUNP_x  : U,V,T,Q tendencies: physics, NOT to be perturbed (LRADCLR_SDT/iSPPT options)
!     PPHY_x  : U,V,T,Q tendencies: physics, to be perturbed
!     PA      : cloud fraction
!     PAP     : full-level pressure
!     PAPH    : half-level pressure
!     KTRAC   : tracer dimension
!
!     OUTPUT ARGUMENTS.

!     PTENx   : U,V,T,Q tendencies: totals (including perturbations)
!     PPERT%PSPPTGFIX : Initial (unperturbed) and final (perturbed) physics tendencies modified by SPPT
!

!     EXTERNALS.  NONE
!     ---------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!     14-Dec-2009  Martin Leutbecher

!     MODIFICATIONS.
!     --------------
!    Aug2011 F.Bouttier: no supersat removal with Arome physics (it's already done elsewhere)
!    Dec2015 S.Massart: perturbation of GFL tendencies
!    Dec2015 A.Innes: CHEM tracers
!    Jan2016 S.Lang    : modifications for SPPTGFIX
!    Jan-2016, SJ Lock : Enabled independent perturbation patterns (iSPPT)
!    Feb-2016, SJ Lock : LSPPTGFL, modifications to CHEM perturbations
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RETV     ,RTT      ,RLSTT    ,RLVTT
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 &                    R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 &                    RALVDCP  ,RALSDCP  ,RTWAT    ,&
 &                    RTICE    ,RTICECU  ,&
 &                    RTWAT_RTICE_R      ,RTWAT_RTICECU_R  ,&
 &                    RKOOP1   ,RKOOP2
USE YOECLDP  , ONLY : TECLDP
USE YOMCT0   , ONLY : LAROME
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOMPHYDER, ONLY : PERTURB_TYPE
USE TYPE_MODEL,ONLY : MODEL

IMPLICIT NONE

TYPE(MODEL)               ,INTENT (INOUT)          :: YDMODEL
TYPE(TYPE_GFLD)           ,INTENT (INOUT)          :: YGFL
TYPE(PERTURB_TYPE)        ,INTENT (INOUT),OPTIONAL :: PPERT

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KSDTN2D
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PTSL(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PQSL(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PA(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PAP(KLON,KLEV)    ! full level pressure
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PAPH(KLON,KLEV+1) ! half level pressure
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PDYN_U(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PDYN_V(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PDYN_T(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PDYN_Q(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PUNP_U(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PUNP_V(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PUNP_T(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PUNP_Q(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PPHY_U(KLON,KLEV,KSDTN2D)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PPHY_V(KLON,KLEV,KSDTN2D)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PPHY_T(KLON,KLEV,KSDTN2D)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PPHY_Q(KLON,KLEV,KSDTN2D)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PMULNOISE(KLON,KSDTN2D)
REAL(KIND=JPRB)   ,INTENT(OUT),OPTIONAL   :: PTENU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT),OPTIONAL   :: PTENV(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT),OPTIONAL   :: PTENT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT),OPTIONAL   :: PTENQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL    :: PTENGFL(KLON,KLEV,YGFL%NDIM1)
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL    :: KTRAC
REAL(KIND=JPRB)   ,INTENT(INOUT),OPTIONAL :: PTENC(:,:,:)

! Local variables
INTEGER(KIND=JPIM) :: JK, JL, J2D
INTEGER(KIND=JPIM), PARAMETER :: ISPECNOPT = -1  ! unperturbed species
INTEGER(KIND=JPIM), PARAMETER :: ISPECPT = 1  ! marker for perturbed species
INTEGER(KIND=JPIM), DIMENSION(YGFL%NCHEM) :: ICHEM

REAL(KIND=JPRB) :: ZTP1(KLON,KLEV), ZQP1(KLON,KLEV)
REAL(KIND=JPRB) :: ZQS, ZFACT, ZZR
REAL(KIND=JPRB) :: ZPERT, ZDELQP, ZRED1, ZRED2, ZQSC, ZDQSDTC, ZRP
REAL(KIND=JPRB) :: ZDELTP,ZDELQAMX, ZDELQP2, ZDQX, ZXKOOP, ZXRETV
REAL(KIND=JPRB) :: ZMEAN(KLON), ZTHICK(KLON), ZTHICK1, ZREDMIN, ZRTHICK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZRPBOT(KLON)
REAL(KIND=JPRB), DIMENSION(KLON,KLEV) :: ZRED, ZREDX
REAL(KIND=JPRB), DIMENSION(KLON,KLEV,KSDTN2D) :: ZTAPER
REAL(KIND=JPRB) :: ZSIG, ZPRESS, ZPRESSLIMTOP, ZSIGMALIMBOT, ZTAPER_ST0, ZTAPER_BL0
LOGICAL :: LLDONTPERT, LLSUPSAT, LLVINLIM
LOGICAL,  DIMENSION(KLON,KLEV) :: LLPERT
LOGICAL, DIMENSION(KLON) :: LLTHICKENOUGH
LOGICAL :: LLALLTHICKENOUGH
INTEGER(KIND=JPIM)  :: IKS
INTEGER(KIND=JPIM), DIMENSION(1) :: IKS1
REAL(KIND=JPRB), DIMENSION(KLEV) :: ZRX

INTEGER(KIND=JPIM) :: JEXT,ITRC
REAL(KIND=JPRB) :: ZPERT_GFL
REAL(KIND=JPRB), DIMENSION(:), ALLOCATABLE :: ZDZPGFL

#include "abor1.intfb.h"
#include "fcttre.func.h"
#include "fccld.func.h"

!-----------------------------------------------------------------------

!     1.  ADD STOCHASTIC PERTURBATION TO PHYSICS TENDENCIES IF REQUIRED
!         (FOR EPS, if LSPSDT=.T. )
!         ----------------------------------------------------------

!     The perturbation is defined as:  dp/p=1+a*(1-2*r)
!     where:
!     p  : tendency
!     dp : perturbed tendency
!     a  : perturbation amplitude
!     r  : random number, 0<r<1

IF (LHOOK) CALL DR_HOOK('SPPTEN',0,ZHOOK_HANDLE)
ASSOCIATE(YDSPPT_CONFIG=>YDMODEL%YRML_GCONF%YRSPPT_CONFIG, &
  &  NSSOPT=>YDMODEL%YRML_PHY_EC%YRECLDP%NSSOPT,NGEMS=>YGFL%NGEMS,NGHG=>YGFL%NGHG, &
  &  YGHG=>YGFL%YGHG, NAERO=>YGFL%NAERO,YAERO=>YGFL%YAERO,NCHEM=>YGFL%NCHEM, &
  &  YCHEM=>YGFL%YCHEM,LSPPTGFL=>YGFL%LSPPTGFL)

IF (LSPPTGFL) THEN
  !for CHEM perturbations
  ALLOCATE(ZDZPGFL(KTRAC))
ENDIF

ZPRESSLIMTOP=0.5_JPRB*(YDSPPT_CONFIG%XPRESSTOP_ST0 + YDSPPT_CONFIG%XPRESSBOT_ST0)
ZSIGMALIMBOT=0.5_JPRB*(YDSPPT_CONFIG%XSIGMATOP + YDSPPT_CONFIG%XSIGMABOT)
ZRTHICK=1._JPRB/YDSPPT_CONFIG%XMEANRED_THICKNESS_SDT
ZZR=1._JPRB/PTSPHY
IF (YDSPPT_CONFIG%LTAPER_BL0) THEN
  DO JL=KIDIA,KFDIA
    ZRPBOT(JL) = 1._JPRB/PAP(JL,KLEV)
  ENDDO
ENDIF

! Record the indices of the species to be perturbed later
ICHEM(:) = ISPECNOPT ! default - no perturbations
DO JEXT=1,NCHEM
  SELECT CASE (TRIM(YCHEM(JEXT)%CNAME))
    CASE ("O3", "CO", "NO2", "CH2O")  ! species to perturb
      ICHEM(JEXT) = ISPECPT
  END SELECT
ENDDO

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    IF (YDSPPT_CONFIG%LSPPTGFIX) THEN
      ! Accumulate and store the unperturbed physics tendencies
      PPERT%PSPPTGFIX(JL,JK,1,1)=PUNP_Q(JL,JK)
      PPERT%PSPPTGFIX(JL,JK,2,1)=PUNP_T(JL,JK)
      PPERT%PSPPTGFIX(JL,JK,3,1)=PUNP_U(JL,JK)
      PPERT%PSPPTGFIX(JL,JK,4,1)=PUNP_V(JL,JK)
      DO J2D=1,KSDTN2D
        PPERT%PSPPTGFIX(JL,JK,1,1)=PPERT%PSPPTGFIX(JL,JK,1,1)+PPHY_Q(JL,JK,J2D)
        PPERT%PSPPTGFIX(JL,JK,2,1)=PPERT%PSPPTGFIX(JL,JK,2,1)+PPHY_T(JL,JK,J2D)
        PPERT%PSPPTGFIX(JL,JK,3,1)=PPERT%PSPPTGFIX(JL,JK,3,1)+PPHY_U(JL,JK,J2D)
        PPERT%PSPPTGFIX(JL,JK,4,1)=PPERT%PSPPTGFIX(JL,JK,4,1)+PPHY_V(JL,JK,J2D)
      ENDDO
    ENDIF

    !Extract tendencies for CHEM perturbations
    IF (LSPPTGFL) THEN
      ITRC=0
      DO JEXT=1,NGHG
        ITRC=ITRC+1
        ZDZPGFL(ITRC)=PTENC(JL,JK, ITRC)-PTENGFL(JL,JK,YGHG(JEXT)%MP1)
      ENDDO
      DO JEXT=1,NAERO
        ITRC=ITRC+1
        ZDZPGFL(ITRC)=PTENC(JL,JK, ITRC)-PTENGFL(JL,JK,YAERO(JEXT)%MP1)
      ENDDO
      DO JEXT=1,NCHEM
        ITRC=ITRC+1
        ZDZPGFL(ITRC)=PTENC(JL,JK, ITRC)-PTENGFL(JL,JK,YCHEM(JEXT)%MP1)
      ENDDO
    ENDIF

    ! calculate tapers
    ZTAPER_ST0=1._JPRB
    ZTAPER_BL0=1._JPRB
    IF (YDSPPT_CONFIG%LTAPER_BL0) THEN
      !
      !   taper perturbations to zero in the boundary layer
      !
      ZSIG=PAP(JL,JK)*ZRPBOT(JL) ! pseudo vertical coordinate for tapering
      IF (ZSIG > YDSPPT_CONFIG%XSIGMABOT) THEN
        !
        !   below transition layer
        !
        ZTAPER_BL0=0._JPRB
      ELSE
        IF (ZSIG > YDSPPT_CONFIG%XSIGMATOP) THEN
          !
          !   in transition layer
          !
          ZTAPER_BL0= ((YDSPPT_CONFIG%XTAPER3*ZSIG + YDSPPT_CONFIG%XTAPER2)*ZSIG &
                       & + YDSPPT_CONFIG%XTAPER1)*ZSIG + YDSPPT_CONFIG%XTAPER0
        ENDIF
      ENDIF
    ENDIF
    IF (YDSPPT_CONFIG%LTAPER_ST0) THEN
      ZPRESS=PAP(JL,JK)
      !
      !   taper perturbations to zero above tropopause
      !
      IF (ZPRESS < YDSPPT_CONFIG%XPRESSTOP_ST0) THEN
        !
        !   above transition layer
        !
        ZTAPER_ST0=0._JPRB
      ELSE
        IF (ZPRESS < YDSPPT_CONFIG%XPRESSBOT_ST0) THEN
          !
          !   in transition layer
          !
          ZTAPER_ST0 = ((YDSPPT_CONFIG%XTAPER3_ST0*ZPRESS + YDSPPT_CONFIG%XTAPER2_ST0)*ZPRESS &
                        & + YDSPPT_CONFIG%XTAPER1_ST0)*ZPRESS + YDSPPT_CONFIG%XTAPER0_ST0
        ENDIF
      ENDIF
    ENDIF

    ! Accumulate the total perturbed tendencies:
    !  = unperturbed tendencies (dynamics and optional physics)
    PTENU(JL,JK)=PDYN_U(JL,JK) + PUNP_U(JL,JK)
    PTENV(JL,JK)=PDYN_V(JL,JK) + PUNP_V(JL,JK)
    PTENT(JL,JK)=PDYN_T(JL,JK) + PUNP_T(JL,JK)
    PTENQ(JL,JK)=PDYN_Q(JL,JK) + PUNP_Q(JL,JK)
    ! Also for the supersaturation check
    ZTP1 (JL,JK)=PTSL(JL,JK) + PUNP_T(JL,JK)
    ZQP1 (JL,JK)=PQSL(JL,JK) + PUNP_Q(JL,JK)
    ! plus the perturbed physics tendencies
    DO J2D=1,KSDTN2D
      ZTAPER(JL,JK,J2D) =  1.0_JPRB
      IF (YDSPPT_CONFIG%LTAPER_BL0 .AND. YDSPPT_CONFIG%LTAPER_BLI(J2D)) ZTAPER(JL,JK,J2D) = ZTAPER(JL,JK,J2D)*ZTAPER_BL0
      IF (YDSPPT_CONFIG%LTAPER_ST0 .AND. YDSPPT_CONFIG%LTAPER_STI(J2D)) ZTAPER(JL,JK,J2D) = ZTAPER(JL,JK,J2D)*ZTAPER_ST0
      ZPERT = ZTAPER(JL,JK,J2D)*PMULNOISE(JL,J2D)
      ZFACT = 1._JPRB+ZPERT
      PTENU(JL,JK)=PTENU(JL,JK) + ZFACT*PPHY_U(JL,JK,J2D)
      PTENV(JL,JK)=PTENV(JL,JK) + ZFACT*PPHY_V(JL,JK,J2D)
      PTENT(JL,JK)=PTENT(JL,JK) + ZFACT*PPHY_T(JL,JK,J2D)
      PTENQ(JL,JK)=PTENQ(JL,JK) + ZFACT*PPHY_Q(JL,JK,J2D)
      ! for the supersaturation check
      ZTP1 (JL,JK)=ZTP1 (JL,JK) + ZPERT*PPHY_T(JL,JK,J2D)*PTSPHY
      ZQP1 (JL,JK)=ZQP1 (JL,JK) + ZPERT*PPHY_Q(JL,JK,J2D)*PTSPHY
    ENDDO

    !Perturb CHEM tendencies
    IF (LSPPTGFL) THEN
      J2D=1                        !temporary fix: always use pattern number 1
      ZPERT_GFL=PMULNOISE(JL,J2D)
      ZFACT = 1._JPRB+ZPERT_GFL
      ITRC=0
      DO JEXT=1,NGHG
        ITRC=ITRC+1
        PTENC(JL,JK, ITRC) = PTENGFL(JL,JK,YGHG(JEXT)%MP1) &
        &                   +  ZFACT * ZDZPGFL(ITRC)
      ENDDO
      DO JEXT=1,NAERO
        ITRC=ITRC+1
        PTENC(JL,JK, ITRC) = PTENGFL(JL,JK,YAERO(JEXT)%MP1) &
        &                   +  ZFACT * ZDZPGFL(ITRC)
      ENDDO
      DO JEXT=1,NCHEM     !It crashes if all CHEM tendencies are perturbed
        ITRC=ITRC+1
        IF (ICHEM(JEXT) == ISPECPT) THEN    ! O3, CO, NO2 and CH2O
          PTENC(JL,JK, ITRC) = PTENGFL(JL,JK,YCHEM(JEXT)%MP1) &
        &                   +  ZFACT * ZDZPGFL(ITRC)
        ELSE
          PTENC(JL,JK, ITRC) = PTENGFL(JL,JK,YCHEM(JEXT)%MP1) &
        &                            + ZDZPGFL(ITRC)
        ENDIF
      ENDDO
    ENDIF

    ! ------------------------------
    !     Supersaturation check
    ! ------------------------------
    ZQS=FOEEWM(ZTP1(JL,JK))/PAP(JL,JK)
    ZQS=MIN(0.5_JPRB,ZQS)
    ZQS=ZQS/(1.0_JPRB-RETV*ZQS)

    !
    ! supersaturation adjustments
    !
    IF (ZTP1(JL,JK)<RTICE.AND.NSSOPT>0.AND..NOT.(LAROME)) &
    & ZQS=ZQS*(PA(JL,JK)+FOKOOP(ZTP1(JL,JK))*(1.0_JPRB-PA(JL,JK)))

    LLDONTPERT=.FALSE.
    IF (YDSPPT_CONFIG%LQPERTLIMIT2) THEN
      LLDONTPERT=LLDONTPERT.OR.(ZQP1(JL,JK)< 0._JPRB)
    ENDIF
    IF (YDSPPT_CONFIG%NQSAT_SDT==0) THEN
      LLDONTPERT=LLDONTPERT.OR.(ZQP1(JL,JK)> ZQS)
    ENDIF
    LLPERT(JL,JK)=.NOT.LLDONTPERT
    IF (LLDONTPERT) THEN  ! accumulate total unperturbed tendencies
      PTENQ(JL,JK)=PDYN_Q(JL,JK) + PUNP_Q(JL,JK)
      PTENT(JL,JK)=PDYN_T(JL,JK) + PUNP_T(JL,JK)
      DO J2D=1,KSDTN2D
        PTENQ(JL,JK)=PTENQ(JL,JK) + PPHY_Q(JL,JK,J2D)
        PTENT(JL,JK)=PTENT(JL,JK) + PPHY_T(JL,JK,J2D)
      ENDDO
    ENDIF
    !
    !   determine moderation factor using linearised
    !        saturation humidity function
    !
    LLSUPSAT=(ZQP1(JL,JK) > ZQS)
    LLVINLIM=.TRUE.
    IF (YDSPPT_CONFIG%LTAPER_BL0) LLVINLIM=LLVINLIM.AND.(ZSIG<ZSIGMALIMBOT)
    IF (YDSPPT_CONFIG%LTAPER_ST0) LLVINLIM=LLVINLIM.AND.(ZPRESS>ZPRESSLIMTOP)
    ZRED(JL,JK) = 1._JPRB
    IF ((YDSPPT_CONFIG%NQSAT_SDT==3).AND.LLPERT(JL,JK).AND.LLSUPSAT.AND.LLVINLIM) THEN
      ZRP=1._JPRB/PAP(JL,JK)
      ZQSC=FOEEWM(PTSL(JL,JK))*ZRP ! satur. humidity at final unperturbed point
      ZQSC=MIN(0.5_JPRB,ZQSC)
      ZXRETV=1.0_JPRB/(1.0_JPRB-RETV*ZQSC)
      ZQSC=ZQSC*ZXRETV
      !
      ! d(sat. hum.)/dT at final unpert. point
      !
      ZDQSDTC=FOE_DEWM_DT(PTSL(JL,JK))*ZRP*ZXRETV
      !
      ! adjust for homogeneous nucleation at cold temperatures
      !
      IF (PTSL(JL,JK)<RTICE.AND.NSSOPT>0) THEN
        ZXKOOP =PA(JL,JK)+FOKOOP(PTSL(JL,JK))*(1.0_JPRB-PA(JL,JK))
        ZQSC   =ZQSC   *ZXKOOP
        ZDQSDTC=ZDQSDTC*ZXKOOP
      ENDIF

      ZRED1=1._JPRB
      ZRED2=1._JPRB
      IF (ZQP1(JL,JK)> ZQS ) THEN
        ZDELQP=ZQP1(JL,JK)-PQSL(JL,JK) ! pos. Q pert.
        ZDELTP=ZTP1(JL,JK)-PTSL(JL,JK) ! pos. T pert.
        ZDELQP2=ZDQSDTC*ZDELTP  ! change in sat. hum. due to T pert
        ZDELQAMX=MAX(ABS(ZDELQP),ABS(ZDELQP2),1.E-9_JPRB)
        ZDQX=ZDELQP - ZDELQP2
        IF (ABS(ZDQX)>=0.001_JPRB*ZDELQAMX) THEN
          ZRED1 = (ZQSC - PQSL(JL,JK))/ZDQX
        ELSE
          ZRED1 =0._JPRB
        ENDIF
      ENDIF
      ZRED(JL,JK) = MIN( 1._JPRB, ZRED1, ZRED2)
      ZRED(JL,JK) = MAX( 0._JPRB, ZRED(JL,JK))
    ENDIF
  ENDDO
ENDDO
!
!  vertical consistency
!
IF (YDSPPT_CONFIG%NQSAT_SDT==3) THEN
  ZTHICK(:)=1.E-12_JPRB
  ZMEAN(:) =1.E-12_JPRB
  LLTHICKENOUGH(:)=.FALSE.
  ZREDX(:,:) =  ZRED(:,:)
  JKLOOP: DO JK=1,KLEV
    JLLOOP: DO JL=KIDIA,KFDIA
      !
      !   Is layer thick enough?
      !
      IF (LLTHICKENOUGH(JL)) THEN
        CYCLE JLLOOP
      ELSE
        ZRX(:) = ZREDX(JL,:)
        IKS1 = MINLOC(ZRX)
        IKS  = IKS1(1)
        ZREDX(JL, IKS) = 99._JPRB
        ZTHICK1=PAPH(JL,IKS+1)-PAPH(JL,IKS)
        ZREDMIN  =ZRED(JL,IKS)
        ZMEAN(JL)  = ZMEAN(JL)  + ZTHICK1*ZREDMIN
        ZTHICK(JL) = ZTHICK(JL) + ZTHICK1
        IF ((ZREDMIN>0.999999999_JPRB).OR.&
             &(ZTHICK(JL)>=YDSPPT_CONFIG%XMEANRED_THICKNESS_SDT)) THEN
          LLTHICKENOUGH(JL)=.TRUE.
          ZTHICK1=YDSPPT_CONFIG%XMEANRED_THICKNESS_SDT-ZTHICK(JL)
          ZMEAN(JL)=ZMEAN(JL) + ZTHICK1*ZREDMIN
          ZTHICK(JL)=YDSPPT_CONFIG%XMEANRED_THICKNESS_SDT
        ENDIF
      ENDIF
    ENDDO JLLOOP
    LLALLTHICKENOUGH=ALL(LLTHICKENOUGH(KIDIA:KFDIA))
    IF (LLALLTHICKENOUGH) THEN
      EXIT JKLOOP
    ENDIF
  ENDDO JKLOOP
  IF (.NOT.LLALLTHICKENOUGH) THEN
    CALL ABOR1('XMEANRED_THICKNESS_SDT is too large; LLALLTHICKENOUGH=FALSE despite JK=KLEV')
  ENDIF
  DO JL=KIDIA,KFDIA
    ZMEAN(JL) = ZMEAN(JL)*ZRTHICK
    ZRED(JL,:)= ZMEAN(JL)
  ENDDO
  ! Accumulate the total perturbed tendencies,
  !   applying reduction factor for T/Q perturbations
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      IF ((ZRED(JL,JK) /= 1._JPRB).AND.LLPERT(JL,JK)) THEN
        ! = unperturbed tendencies (dynamics and optional physics)
        PTENT(JL,JK)=PDYN_T(JL,JK) + PUNP_T(JL,JK)
        PTENQ(JL,JK)=PDYN_Q(JL,JK) + PUNP_Q(JL,JK)
        ! + perturbed physics tendencies
        DO J2D=1,KSDTN2D
          ZPERT= ZRED(JL,JK)*ZTAPER(JL,JK,J2D)*PMULNOISE(JL,J2D)
          ZFACT=1._JPRB + ZPERT
          PTENT(JL,JK)=PTENT(JL,JK) + ZFACT*PPHY_T(JL,JK,J2D)
          PTENQ(JL,JK)=PTENQ(JL,JK) + ZFACT*PPHY_Q(JL,JK,J2D)
        ENDDO
      ENDIF
    ENDDO
  ENDDO
ENDIF

IF (YDSPPT_CONFIG%LSPPTGFIX) THEN
  ! Store total perturbed physics tendencies
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PPERT%PSPPTGFIX(JL,JK,1,2)=PTENQ(JL,JK)-PDYN_Q(JL,JK)
      PPERT%PSPPTGFIX(JL,JK,2,2)=PTENT(JL,JK)-PDYN_T(JL,JK)
      PPERT%PSPPTGFIX(JL,JK,3,2)=PTENU(JL,JK)-PDYN_U(JL,JK)
      PPERT%PSPPTGFIX(JL,JK,4,2)=PTENV(JL,JK)-PDYN_V(JL,JK)
    ENDDO
  ENDDO
ENDIF

IF (LSPPTGFL) THEN
  DEALLOCATE(ZDZPGFL)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPPTEN',1,ZHOOK_HANDLE)
END SUBROUTINE SPPTEN
