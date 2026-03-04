! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#ifdef RS6K
@PROCESS NOSTRICT
#endif
SUBROUTINE VDFEXCU(YDEGWD,YDEPHY,YDSPP_CONFIG,KIDIA  , KFDIA  , KLON   , KLEV   , &
                  &KINV   , KCBASE , KCTOP  , PLAM,  PZ0MM  , &
                  &PHRLW  , PHRSW  , PUM1   , PVM1   , PTM1   , PQM1   , &
                  &PAPHM1 , PGEOM1 , PGEOH  , PGELAT , PCPTGZ , &
                  &PKMFL  , PKHFL  , PKQFL  , PKHVFL , PGP2DSPP,PCFM   , PCFH   , &
                  &PZINV  , PKH    , PKM    , PRI    , PZCLDBASE , KPBLTYPE )
!     ------------------------------------------------------------------

!**   *VDFEXCU* - DETERMINES THE EXCHANGE COEFFICIENTS BETWEEN THE
!                 UPPER MODEL LEVELS WITH STABILITY AS A FUNCTION OF
!                 OBUKHOV-L

!     PURPOSE
!     -------

!     DETERMINE EXCHANGE COEFFICIENTS BETWEEN THE UPPER MODEL LEVELS

!     INTERFACE
!     ---------

!     *VDFEXCU* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):

!     *PUM1*         X-VELOCITY COMPONENT AT T-1
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1
!     *PTM1*         TEMPERATURE AT T-1
!     *PQM1*         SPECIFIC HUMUDITY AT T-1
!     *PAPHM1*       PRESSURE AT HALF LEVELS AT T-1
!     *PGEOM1*       GEOPOTENTIAL AT T-1
!     *PCPTGZ*       DRY STATIC ENERGY
!     *PKMFL*        KINEMATIC MOMENTUM FLUX
!     *PKHFL*        KINEMATIC HEAT FLUX
!     *PKQFL*        KINEMATIC MOISTURE FLUX
!     *PKHVFL*       KINEMATIC BUOYANCY FLUX
!     *PGP2DSPP*     Standard stochastic variable (mean=0, SD=1)
!     *PZINV*        INVERSION HEIGHT                  [M]

!     OUTPUT PARAMETERS (REAL):

!     *PCFM*         PROP. TO EXCH. COEFF. FOR MOMENTUM (C-STAR IN DOC.)
!     *PCFH*         PROP. TO EXCH. COEFF. FOR HEAT     (C-STAR IN DOC.)
!                    (ONLY PCFM(*,1:KLEV-1) AND
!                          PCFH(*,1:KLEV-1) ARE COMPUTED)
!     *PKH*          TURB. DIFF. COEFF. FOR HEAT ABOVE SURF. LAY.  (M2/S)
!     *PKM*          TURB. DIFF. COEFF. FOR MOMENTUM               (M2/S)
!     *PRI*          RICHARDSON NUMBER

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     AUTHOR.
!     -------
!      A.C.M. BELJAARS  26/03/90.  Original

!     MODIFICATIONS.
!     --------------
!      J.HAGUE          13/01/2003 MASS Vector Functions
!      M. Ko"hler        3/12/2004 Moist Advection-Diffusion incl.
!                                  K,cloud and cloud top entrainment
!      P. Lopez         02/06/2005 Removed option for linearized
!                                  physics (now called separately)
!      M. Ko"hler        1/11/2007 reduced K diffusion above surface or mixed layer
!      M. Ko"hler        2/15/2008 added shear for K and Ri calculation
!      N. Semane+P.Becht 8/06/2012 scaling for small planet
!      A.Beljaars+I.Sandu 1/11/2012 smooth reduction in diffusion above tropopause 
!      +T.Stockdale+P.Bechtold (optional,by default : not active)
!      N.Semane+P.Bechtold 04-10-2012 Add RPLRG/RPLDARE factors for small planet
!      I. Sandu+A.Beljaars 15/3/2013 changed treatment of diffusion in stable conditions,i.e.
!                              LTG functions allover, asymptotic mixing length proportional to PBL
!                              height within stable boundary layers and equal to 30m in free-shear layers,
!                              removed non-resolved shear term
!      F. Vana  05-Mar-2015  Support for single precision
!      F. Vana  17-Dec-2015  More support for single precision
!      M. Leutbecher & S.-J. Lock (Jan 2016) Introduced SPP scheme (LSPP)
!      M. Leutbecher & S. Lang (Oct 2020) SPP abstraction and revision
!     -----------------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST    , ONLY : RG, RD, RCPD, RETV, RLVTT
USE PARPHY    , ONLY : RKAP,REPDU2
USE YOEVDFS   , ONLY : RCHBA, RCHBB, RCHBD, RCHB23A, RCHBBCD, RCHBCD, &
 &                     RCHETA, RCHETB, RCDHALF, RCDHPI2
USE YOEGWD    , ONLY : TEGWD
USE YOEPHY    , ONLY : TEPHY
USE YOMDYNCORE, ONLY : RPLRG, RPLDARE

USE SPP_MOD        , ONLY : TSPP_CONFIG
USE SPP_GEN_MOD    , ONLY : SPP_PERT

IMPLICIT NONE

TYPE(TEGWD)       ,INTENT(INOUT) :: YDEGWD
TYPE(TEPHY)       ,INTENT(INOUT) :: YDEPHY
TYPE(TSPP_CONFIG) ,INTENT(IN)    :: YDSPP_CONFIG
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KINV(KLON)
INTEGER(KIND=JPIM),INTENT(IN)    :: KCBASE(KLON)
INTEGER(KIND=JPIM),INTENT(IN)    :: KCTOP(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAM
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0MM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRLW(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRSW(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAT(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTGZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP2DSPP(KLON,YDSPP_CONFIG%SM%NRFTOTAL)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKMFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKHFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKQFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKHVFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCFM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCFH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZINV(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZCLDBASE(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPBLTYPE(KLON) 

REAL(KIND=JPRD) ::    ZRI(KLON),ZMGEOM(KLON),ZUST(KLON),&
                    & ZDTV(KLON),ZL(KLON),ZPHIM(KLON),&
                    & ZPHIH(KLON)
REAL(KIND=JPRB) ::    ZP(KLON), ZDU2(KLON), ZPRKAP(KLON)

INTEGER(KIND=JPIM) :: JK, JL
INTEGER(KIND=JPIM) :: IINV(KLON)

! A bunch of SPP variables
LOGICAL            :: LLPERT_VDEXC ! SPP perturbation on?
LOGICAL            :: LLPERT_RKAP, LLPERT_RKAP1, LLPERT_RKAP2, LLPERT_RKAP3 ! SPP perturbation on?
INTEGER(KIND=JPIM) :: IPN        ! SPP perturbation pointer
INTEGER(KIND=JPIM) :: IP1, IPRKAP, IPRKAP1, IPRKAP2, IPRKAP3 ! SPP random field pointer
TYPE(SPP_PERT)     :: PN1, PNRKAP, PNRKAP1, PNRKAP2, PNRKAP3 ! SPP pertn. configs

REAL(KIND=JPRB) ::    ZENTRSFC, ZENTRRAD, ZKLEN2, &
                    & ZCB, ZCD, ZCFNC1, ZCONS13, ZRG, &
                    & ZDH, ZDL, ZDRORO, ZEPS, &
                    & ZPHIKH, ZPHIKM, ZSCF, &
                    & ZZ, ZWTVENTR, ZKH, &
                    & ZML, ZBASE, ZVSC, ZKCLD, ZDRADFLX(KLON), &
                    & ZREPUST,ZRA, ZLAT, ZETA, ZRCPD
REAL(KIND=JPRB) ::    ZZH, ZIFLTGM, ZIFLTGH, ZIFMOM(KLON), ZIFMOH(KLON), ZDUDZ(KLON)
LOGICAL ::            LLDONE(KLON)
REAL(KIND=JPRB) ::    ZPBLHEIGHT(KLON),ZRIB(KLON),ZSVBOT(KLON),ZRILEV,&
                    & ZRICRI,ZKLENT(KLON,KLEV),ZSV, ZZNORM, ZLMIN

REAL(KIND=JPRB) ::    ZEPSILON

REAL(KIND=JPHOOK) ::    ZHOOK_HANDLE

#include "surf_inq.h"

#include "fcvdfs.func.h"

!     ------------------------------------------------------------------

!*         1.     INITIALIZE CONSTANTS
!                 --------------------

IF (LHOOK) CALL DR_HOOK('VDFEXCU',0,ZHOOK_HANDLE)

ZEPSILON=100._JPRB*EPSILON(ZEPSILON)

ZENTRSFC  = 0.2_JPRB       ! factor for surface based top entrainment 
ZENTRRAD  = 0.2_JPRB       ! factor for radiative based top entrainment 
ZCD       = 1.0_JPRB
ZCB       = 5.0_JPRB
ZEPS      = 1.E-10_JPRB
ZRICRI   = 0.25_JPRB

CALL SURF_INQ(YDEPHY%YSURF,PREPUST=ZREPUST)

! optimization
ZRG       = 1.0_JPRB/RG
ZCONS13   = 1.0_JPRB/3._JPRB
ZRCPD     = 1.0_JPRB/RCPD

ZLMIN=30.0_JPRB

!  Prepare SPP -------------------------------------------------------------

IF (YDSPP_CONFIG%LSPP) THEN
  IPN = YDSPP_CONFIG%PPTR%VDEXC
  LLPERT_VDEXC= IPN > 0
  IF (LLPERT_VDEXC) THEN
    PN1  = YDSPP_CONFIG%SM%PN(IPN)
    IP1 = PN1%MP
  ENDIF

  IPN = YDSPP_CONFIG%PPTR%RKAP
  LLPERT_RKAP= IPN > 0
  IF (LLPERT_RKAP) THEN
    PNRKAP  = YDSPP_CONFIG%SM%PN(IPN)
    IPRKAP = PNRKAP%MP
  ENDIF

  IPN = YDSPP_CONFIG%PPTR%RKAP1
  LLPERT_RKAP1= IPN > 0
  IF (LLPERT_RKAP1) THEN
    PNRKAP1  = YDSPP_CONFIG%SM%PN(IPN)
    IPRKAP1 = PNRKAP1%MP
  ENDIF

  IPN = YDSPP_CONFIG%PPTR%RKAP2
  LLPERT_RKAP2= IPN > 0
  IF (LLPERT_RKAP2) THEN
    PNRKAP2  = YDSPP_CONFIG%SM%PN(IPN)
    IPRKAP2 = PNRKAP2%MP
  ENDIF

  IPN = YDSPP_CONFIG%PPTR%RKAP3
  LLPERT_RKAP3= IPN > 0
  IF (LLPERT_RKAP3) THEN
    PNRKAP3  = YDSPP_CONFIG%SM%PN(IPN)
    IPRKAP3 = PNRKAP3%MP
  ENDIF

ELSE
  LLPERT_RKAP  =.FALSE.
  LLPERT_RKAP1 =.FALSE.
  LLPERT_RKAP2 =.FALSE.
  LLPERT_RKAP3 =.FALSE.
  LLPERT_VDEXC =.FALSE.
ENDIF

IF (LLPERT_VDEXC) THEN
  DO JL=KIDIA,KFDIA
    ZP(JL)= EXP(PN1%MU(1)+PN1%XMAG(1)*PGP2DSPP(JL, IP1))
  ENDDO
ENDIF

DO JL=KIDIA,KFDIA
  IF (KPBLTYPE(JL) == 0 .AND. LLPERT_RKAP) THEN
    ZPRKAP(JL) = EXP(PNRKAP%MU(1)+PNRKAP%XMAG(1)*PGP2DSPP(JL, IPRKAP))
  ELSEIF (KPBLTYPE(JL) == 1 .AND. LLPERT_RKAP1) THEN
    ZPRKAP(JL) = EXP(PNRKAP1%MU(1)+PNRKAP1%XMAG(1)*PGP2DSPP(JL, IPRKAP1))
  ELSEIF (KPBLTYPE(JL) == 2 .AND. LLPERT_RKAP2) THEN
    ZPRKAP(JL) = EXP(PNRKAP2%MU(1)+PNRKAP2%XMAG(1)*PGP2DSPP(JL, IPRKAP2))
  ELSEIF (KPBLTYPE(JL) == 3 .AND. LLPERT_RKAP3) THEN
    ZPRKAP(JL) = EXP(PNRKAP3%MU(1)+PNRKAP3%XMAG(1)*PGP2DSPP(JL, IPRKAP3))
  ELSE ! do not perturb
    ZPRKAP(JL) = 1.0_JPRB
  ENDIF
ENDDO

!     ------------------------------------------------------------------

!*         2.     PREPARE SCALING COEFFICIENTS FOR MIXED LAYER
!                 --------------------------------------------

DO JL=KIDIA,KFDIA
  ZL(JL)  = 1._JPRD
  ZUST  (JL)=MAX(SQRT(PKMFL(JL)),ZREPUST)
!shifts the inversion level by one as pgeoh defined from 0 to klev, while in vdfh is from 1 to klev+1
  IINV(JL)=KINV(JL)-1
  IF (PKHVFL(JL)  <  0.0_JPRB) THEN
    ZL(JL)  = ZUST (JL)**3*PTM1(JL,KLEV)/(RKAP * ZPRKAP(JL) *RG*(PKHVFL(JL)-ZEPS))
    ZL(JL)  = SIGN(MAX(ZEPSILON,ABS(ZL(JL))),ZL(JL))
  ENDIF
ENDDO

!          Calculate PBL cloud top radiative flux jump [Km/s] (cooling)
!          for top-driven K and entrainment velocity formulations.

DO JL=KIDIA,KFDIA
  ZDRADFLX(JL) = 0.0_JPRB
  IF ( KPBLTYPE(JL) == 2 ) THEN
     JK =IINV(JL) 
     ZDRADFLX(JL) = -(MIN(PHRLW(JL,JK)+PHRSW(JL,JK),PHRLW(JL,JK-1)+PHRSW(JL,JK-1)))*(PGEOH(JL,JK-1)-PGEOH(JL,JK))*ZRG
     ZDRADFLX(JL) = MAX( ZDRADFLX(JL)/RPLDARE, 0.0_JPRB )    !safety against rad. heating cases
  ENDIF
ENDDO

!     ------------------------------------------------------------------

!*         3.     VERTICAL LOOP
!                 -------------
!   compute pbl height in stable boundary layers
 DO JL=KIDIA,KFDIA
   LLDONE(JL)=.FALSE.
   ZPBLHEIGHT(JL)=0.0_JPRB
   ZRIB(JL)=0.0_JPRB
   ZSVBOT(JL)=RCPD*PTM1(JL,KLEV)*(1.0_JPRB+RETV*PQM1(JL,KLEV))+PGEOM1(JL,KLEV)+ZEPS
 ENDDO
 DO JK = KLEV-1, 1, -1
    DO JL=KIDIA,KFDIA
      IF (.NOT. LLDONE(JL) .AND. PKHVFL(JL)  >  0.0_JPRB) THEN
        ZSV=RCPD*PTM1(JL,JK)*(1.0_JPRB+RETV*PQM1(JL,JK))+PGEOM1(JL,JK)
!   pbl height diag, which considers the winds close to the surf eq to 0
        ZDU2(JL)=MAX(REPDU2, PUM1(JL,JK)**2+PVM1(JL,JK)**2)  
        ZDRORO=(ZSV-ZSVBOT(JL))&
         & /(ZSV-PGEOM1(JL,JK))  
        ZRILEV=(PGEOM1(JL,JK)-PGEOM1(JL,KLEV))*ZDRORO/ZDU2(JL)
!
          IF (ZRILEV  >  ZRICRI) THEN
           ZPBLHEIGHT(JL)=( (ZRILEV-ZRICRI)*PGEOM1(JL,JK+1)&
           & +(ZRICRI-ZRIB(JL))*PGEOM1(JL,JK) )/&
           & ((ZRILEV-ZRIB(JL))*RG)  
           ZRIB(JL)=ZRILEV
           LLDONE(JL)=.TRUE.
          ELSE
           ZRIB(JL)=ZRILEV
          ENDIF
      ENDIF
    ENDDO
 ENDDO

!initialization
DO JK = KLEV, 1, -1
  DO JL=KIDIA,KFDIA
    ZKLENT(JL,JK) =0.0_JPRB
  ENDDO
ENDDO

!***
DO JK = KLEV-1, 1, -1
!***

  DO JL=KIDIA,KFDIA
    PCFM(JL,JK)=0.0_JPRB
    PCFH(JL,JK)=0.0_JPRB
    PKH(JL,JK) =0.0_JPRB
    PKM(JL,JK) =0.0_JPRB
  ENDDO

!          COMPUTE RI-NUMBER

  DO JL=KIDIA,KFDIA
    ZMGEOM(JL)=PGEOM1(JL,JK)-PGEOM1(JL,JK+1)
    ZDU2(JL)=MAX(REPDU2,(PUM1(JL,JK)-PUM1(JL,JK+1))**2&
                     & +(PVM1(JL,JK)-PVM1(JL,JK+1))**2)
    ZDUDZ(JL)= SQRT(ZDU2(JL)) /ZMGEOM(JL)*RG
    ZRI(JL)=PRI(JL,JK)
  ENDDO

  DO JL = KIDIA, KFDIA

!          DIMENSIONLESS COEFFICIENTS MULTIPLIED BY PRESSURE
!          THICKNESSES FOR MOMENTUM AND HEAT EXCHANGE.

    ZZH     = 0.5_JPRB * ZRG * (PGEOM1(JL,JK)+PGEOM1(JL,JK+1)) + PZ0MM(JL)

!          STATICALLY STABLE

    IF ( ZRI(JL) > 0.0_JPRB ) THEN
!     ASYMPTOTIC MIXING LENGTH FOR STABLE SITUATIONS
      IF (PGEOM1(JL,JK)*ZRG <= ZPBLHEIGHT(JL)) THEN
        IF (LLPERT_VDEXC) THEN
          ZKLENT(JL,JK)=MAX(ZLMIN,ZP(JL)*0.1_JPRB*ZPBLHEIGHT(JL)*RPLRG)
        ELSE
          ZKLENT(JL,JK)=MAX(ZLMIN,0.1_JPRB*ZPBLHEIGHT(JL)*RPLRG)
        ENDIF
        ZKLENT(JL,JK)=MIN(300.0_JPRB,ZKLENT(JL,JK))
      ELSE
        ZKLENT(JL,JK)=ZLMIN
      ENDIF
      ZKLENT(JL,JK)=ZKLENT(JL,JK)/RPLRG

      ZKLEN2  = RKAP * ZPRKAP(JL) * ZZH * ZKLENT(JL,JK) / ( RKAP * ZPRKAP(JL) * ZZH + ZKLENT(JL,JK))

!          COMPUTE STABILITY FUNCTIONS
      ZSCF    = SQRT(1.0_JPRB+ZCD*ZRI(JL))
      ZIFLTGM = 1.0_JPRB / (1.0_JPRB + 2.0_JPRB *ZCB * ZRI(JL)/ZSCF) !F(LTG),M
      ZIFLTGH = 1.0_JPRB / (1.0_JPRB + 2.0_JPRB *ZCB * ZRI(JL)*ZSCF) !F(LTG),H

      PCFM(JL,JK) = ZDUDZ(JL) * ZKLEN2**2 * ZIFLTGM
      PCFH(JL,JK) = ZDUDZ(JL) * ZKLEN2**2 * ZIFLTGH
      ZPHIM(JL) = 0.0_JPRB
      ZPHIH(JL) = 0.0_JPRB
!          STATICALLY UNSTABLE

    ELSE
!     ASYMPTOTIC MIXING LENGTH FOR UNSTABLE SITUATIONS
      ZKLENT(JL,JK)=PLAM
      ZKLEN2  = RKAP * ZPRKAP(JL) * ZZH * ZKLENT(JL,JK) / ( RKAP * ZPRKAP(JL) * ZZH + ZKLENT(JL,JK) )

!          COMPUTE STABILITY FUNCTIONS
      ZETA  = ZRI(JL)
      ZPHIM(JL) = PHIMU(ZETA)
      ZPHIH(JL) = PHIHU(ZETA)
!JJJ Forcing compiler to create reciprocal may not give most accurate computation, 
!JJJ but gives bit reproducibility (in T159) with previous 38r2 version
      ZIFMOM(JL)  = 1.0_JPRB / (ZPHIM(JL)**2)                              !F(MO),M
      ZIFMOH(JL)  = 1.0_JPRB / (ZPHIM(JL)*ZPHIH(JL))                       !F(MO),H
      PCFM(JL,JK) = ZDUDZ(JL) * ZKLEN2**2 * ZIFMOM(JL)
      PCFH(JL,JK) = ZDUDZ(JL) * ZKLEN2**2 * ZIFMOH(JL)
    ENDIF

!          ADD MIXED LAYER PARAMETRIZATION (SURF. AND OUTER LAYER)

!    IF (PGEOH(JL,JK-1)*ZRG <= PZINV(JL) ) THEN  ! up to level below entr. level
    IF (JK > IINV(JL) .AND. IINV(JL) < (KLEV-1) ) THEN  ! up to level below entr. level


      ZZ     = PGEOH(JL,JK)*ZRG
      ZDH    = ZZ/MAX(ZEPSILON,(PGEOH(JL,IINV(JL))*ZRG))
      ZDL    = ZZ/ZL(JL)
      ZPHIKH = (1.0_JPRB-39._JPRB*ZDL)**(-ZCONS13)
      ZPHIKM = (1.0_JPRB-15._JPRB*ZDL)**(-ZCONS13)

!          K,surface

      ZZNORM= ZZ * (1.0_JPRB-ZDH)**2
      IF ( KPBLTYPE(JL) >= 2 .AND. KCBASE(JL)>KCTOP(JL)+1)THEN
        IF(RCPD*PKHFL(JL)<-0.1_JPRB.AND.-RCPD*PKHFL(JL)<-0.25_JPRB*RLVTT*PKQFL(JL)) THEN
           ZZNORM=ZZ * (1.0_JPRB-ZDH)**4
           ZDH   =MIN(1.0_JPRB,ZZ/MAX(ZEPSILON,PGEOH(JL,KCBASE(JL)-1)*ZRG) )
           ZZNORM=ZZNORM + ZZ * (1.0_JPRB-ZDH)**2
        ENDIF
      ENDIF

      PCFH(JL,JK) = RKAP * ZPRKAP(JL)/ ZPHIKH * ZUST(JL) * ZZNORM
      PCFM(JL,JK) = RKAP * ZPRKAP(JL)/ ZPHIKM * ZUST(JL) * ZZNORM

!          add cloud-top driven K (Lock et al. 2000, MWR p3187f, equ. 5)
!          (using simplified radiative velocity scale as in Lock, 1998, equ. 12
!          and ignore buoyancy reversal velocity scale)
!          apply K-cloud throughout full PBL (only for stratocumulus)

      ZML   = PGEOH(JL,IINV(JL))*ZRG               ! mixing depth
      ZBASE = 0.0_JPRB                             ! mixing base

      IF ( KPBLTYPE(JL) == 2 ) THEN  
        ZVSC  = ( RG / PTM1(JL,JK) * ZML * ZDRADFLX(JL) ) ** ZCONS13 
        ZKCLD = 0.85_JPRB * ZPRKAP(JL) * RKAP * ZVSC &
            & * (ZZ-ZBASE) ** 2 / ZML &
            & * ( 1.0_JPRB - (ZZ-ZBASE) / ZML ) ** 0.5_JPRB  
        PCFH(JL,JK) = PCFH(JL,JK) + ZKCLD
        PCFM(JL,JK) = PCFM(JL,JK) + ZKCLD

      ENDIF

    ENDIF

!          ADD ENTRAINMENT
!          entrainment velocity * T,v jump due to lw-radiative cooling & sfc flux
!          (Lock & MacVean, 1999, equ. 11)

!   IF ( PGEOH(JL,JK)*ZRG <= PZINV(JL)  .AND.  PZINV(JL) < PGEOH(JL,JK-1)*ZRG ) THEN
    IF ( JK == IINV(JL) .AND. IINV(JL) < (KLEV-1) ) THEN
      ZWTVENTR = 0.0_JPRB
      ZWTVENTR = - ZENTRSFC * PKHVFL(JL)  !sfc flux top entrainment

      IF ( KPBLTYPE(JL) == 2 ) THEN
        ZWTVENTR = ZWTVENTR + ZENTRRAD * ZDRADFLX(JL) !radiation flux jump
      ENDIF
      ZWTVENTR = MAX(0.0_JPRB,ZWTVENTR)

      ZDTV(JL)=(PCPTGZ(JL,JK-1)-PCPTGZ(JL,JK))*ZRCPD&
     & + RETV*0.5_JPRB * (PQM1(JL,JK-1)-PQM1(JL,JK))&
     & * (PTM1(JL,JK)+PTM1(JL,JK-1))

      IF (ZDTV(JL) > 1.E-3_JPRB) THEN
        ZKH     = ZWTVENTR * (PGEOM1(JL,JK-1)-PGEOM1(JL,JK))*ZRG/ ZDTV(JL) 
      ELSE
        ZKH=0.0_JPRB
      ENDIF
      PCFH(JL,JK) = ZKH
      PCFM(JL,JK) = ZKH

    ENDIF

!          LIMIT K TO 1000000 M2/S FOR SAFETY

    PCFH(JL,JK) = MIN( PCFH(JL,JK), 1000000.0_JPRB )
    PCFM(JL,JK) = MIN( PCFM(JL,JK), 1000000.0_JPRB )


    IF(YDEGWD%LRDIFF_STRATO) THEN
      IF(YDEGWD%NDIFF_STRATO == 1) THEN
! Reduced diffusion in stratosphere just a function of pressure (quadratic from 80-120hPa)
        IF (ZRI(JL) > 0.25_JPRB) THEN
          IF(PAPHM1(JL,JK)>8000.0_JPRB) THEN
            ZRA=(PAPHM1(JL,JK)-8000.0_JPRB)*(PAPHM1(JL,JK)-8000.0_JPRB)/16.0E6_JPRB
            ZRA=MAX(MIN(1.0_JPRB,ZRA),0.002_JPRB)
          ELSE
            ZRA=0.002_JPRB
          ENDIF
          PCFH(JL,JK) = PCFH(JL,JK)*ZRA
          PCFM(JL,JK) = PCFM(JL,JK)*ZRA
        ENDIF
      ELSEIF(YDEGWD%NDIFF_STRATO == 2) THEN
! Reduced diffusion in stratosphere only in tropical lower stratosphere
        IF(ABS(PGELAT(JL))<0.5_JPRB) THEN
        IF (ZRI(JL) > 0.25_JPRB) THEN
          IF(PAPHM1(JL,JK)>8000.0_JPRB) THEN
            ZRA=(PAPHM1(JL,JK)-8000.0_JPRB)*(PAPHM1(JL,JK)-8000.0_JPRB)/16.0E6_JPRB
          ELSEIF(PAPHM1(JL,JK)>1200.0_JPRB) THEN
            ZRA=0.0_JPRB
          ELSEIF(PAPHM1(JL,JK)>800.0_JPRB) THEN
            ZRA=(PAPHM1(JL,JK)-1200.0_JPRB)*(PAPHM1(JL,JK)-1200.0_JPRB)/16.0E4_JPRB
          ELSE
            ZRA=1.0_JPRB
          ENDIF
          ZRA=MAX(MIN(1.0_JPRB,ZRA),(1.0_JPRB+ZRI(JL))**(-4))
          PCFH(JL,JK) = PCFH(JL,JK)*ZRA
          PCFM(JL,JK) = PCFM(JL,JK)*ZRA
        ENDIF
        ENDIF
      ELSEIF(YDEGWD%NDIFF_STRATO == 5) THEN
! Reduced diffusion in tropical lower stratosphere, smooth function of latitude
        IF (ZRI(JL) > 0.25_JPRB) THEN
          IF(PAPHM1(JL,JK)>12000.0_JPRB) THEN
            ZRA=1.0_JPRB
          ELSEIF(PAPHM1(JL,JK)>8000.0_JPRB) THEN
            ZRA=(PAPHM1(JL,JK)-8000.0_JPRB)*(PAPHM1(JL,JK)-8000.0_JPRB)/16.0E6_JPRB
          ELSEIF(PAPHM1(JL,JK)>1200.0_JPRB) THEN
            ZRA=0.0_JPRB
          ELSEIF(PAPHM1(JL,JK)>800.0_JPRB) THEN
            ZRA=(PAPHM1(JL,JK)-1200.0_JPRB)*(PAPHM1(JL,JK)-1200.0_JPRB)/16.0E4_JPRB
          ELSE
            ZRA=1.0_JPRB
          ENDIF
! Restrict reduction to tropics (in extratropics, shear associated with wave activity and noise)
          ZLAT=ABS(PGELAT(JL))*57.296_JPRB
          IF(ZLAT < 20.0_JPRB) THEN
            ZRA=ZRA
          ELSEIF(ZLAT < 30.0_JPRB) THEN
            ZRA=1.0_JPRB-(1.0_JPRB-0.01_JPRB*(ZLAT-20.0_JPRB)*(ZLAT-20.0_JPRB))*(1.0_JPRB-ZRA)
          ELSE
            ZRA=1.0
          ENDIF
          ZRA=MAX(ZRA,(1.0_JPRB+ZRI(JL))**(-4))
          PCFH(JL,JK) = PCFH(JL,JK)*ZRA
          PCFM(JL,JK) = PCFM(JL,JK)*ZRA
        ENDIF
      ENDIF
    ENDIF

!          DIFFUSION COEFFICIENT FOR HEAT FOR POSTPROCESSING ONLY IN (M2/S)

    PKH(JL,JK) = PCFH(JL,JK)*RPLDARE
    PKM(JL,JK) = PCFM(JL,JK)*RPLDARE

!          SCALE DIFFUSION COEFFICIENTS FOR USE IN VDFDIFH/M

    ZCFNC1 = RG * PAPHM1(JL,JK)&
         & /( 0.5_JPRB*RD * ZMGEOM(JL)&
         &    *( PTM1(JL,JK  )*(1.0_JPRB+RETV*PQM1(JL,JK  ))&
         &    +  PTM1(JL,JK+1)*(1.0_JPRB+RETV*PQM1(JL,JK+1)))) 
    PCFH(JL,JK) = PCFH(JL,JK) * ZCFNC1 * RPLDARE
    PCFM(JL,JK) = PCFM(JL,JK) * ZCFNC1 * RPLDARE

  ENDDO

!***
ENDDO
!***

IF (LHOOK) CALL DR_HOOK('VDFEXCU',1,ZHOOK_HANDLE)
END SUBROUTINE VDFEXCU
