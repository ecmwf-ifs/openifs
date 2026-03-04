! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFEXCUS (YDEPHLI, YDEPHY, &
 & KIDIA, KFDIA, KLON, KLEV, &
 & KINV , KCBASE, KCTOP, KPBLTYPE,&
 & PTMST  , PVDIFTS, &
 & PZ0MM  , PUM1   , PVM1   , PTM1  , PQM1 , &
 & PAPHM1 , PGEOM1 , PGEOH , PCPTGZ, &
 & PKMFL  , PKHFL  , PKQFL , PZINV , PCFM , PCFH)
!     ------------------------------------------------------------------

!**   *VDFEXCUS* - DETERMINES THE EXCHANGE COEFFICIENTS BETWEEN THE
!                   UPPER MODEL LEVELS WITH STABILITY AS A FUNCTION OF
!                   RICHARDSON NUMBER
!                   (Nonlinear version for trajectory in adjoint)

!     J.F. MAHFOUF    02/10/95     E.C.M.W.F. 

!     Modified by
!     -----------

!     P. LOPEZ        25/02/05     nonlinear version for
!                                  trajectory in adjoint
!     M. Ko"hler      25/09/07     add mixed layer
!                                  reduce K towards MO higher up 
!     M. Janiskova    19/12/07     modified wind shear for Ri computation
!     M. Janiskova    14/01/11     new transition from LTG to MO in stable
!                                  BL close to surface (as in vdfexcu)
!     P. Lopez        31/07/13     Revision to match full non-linear version changes
!     P. Lopez        15/01/21     use new inversion level from VDFHGHTNS & changed mixed layer.

!     PURPOSE
!     -------

!     DETERMINE EXCHANGE COEFFICIENTS BETWEEN THE UPPER MODEL LEVELS

!     INTERFACE
!     ---------

!     *VDFEXCUS* IS CALLED BY *VDFMAINS*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     *KINV*         Inversion level (from 1 to KLEV+1)
!     *KCBASE*       Cloud base level
!     *KCTOP*        Cloud top level
!     *KPBLTYPE*     0: stable PBL
!                    1: dry convective PBL (no cloud below parcel top)
!                    2: cloudy PBL ("stratocumulus")
!                    3: dry PBL below convection ("cumulus")

!     INPUT PARAMETERS (REAL):

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)
!     *PUM1*         X-VELOCITY COMPONENT AT T-1  (Trajectory)
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1  (Trajectory)
!     *PTM1*         TEMPERATURE AT T-1           (Trajectory)
!     *PQM1*         SPECIFIC HUMUDITY AT T-1     (Trajectory)
!     *PAPHM1*       PRESSURE AT T-1              (Trajectory)
!     *PGEOM1*       GEOPOTENTIAL AT T-1          (Trajectory)
!     *PGEOH*        GEOPOTENTIAL AT HALF LEVELS  (Trajectory)
!     *PCPTGZ*       DRY STATIC ENERGY            (Trajectory)
!     *PKMFL*        KINEMATIC MOMENTUM FLUX   
!     *PKHFL*        KINEMATIC HEAT FLUX       
!     *PKQFL*        KINEMATIC MOISTURE FLUX  
!     *PZINV*        INVERSION HEIGHT         

!     OUTPUT PARAMETERS (REAL):

!     *PCFM*         PROP. TO EXCH. COEFF. FOR MOMENTUM     (Trajectory)
!     *PCFH*         PROP. TO EXCH. COEFF. FOR HEAT         (Trajectory)
!                    (ONLY PCFM(*,1:KLEV-1) AND
!                          PCFH(*,1:KLEV-1) ARE COMPUTED)

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     THE EXCHANGE COEFFICIENTS ARE EXPRESSED ACCORDING TO 
!     LOUIS et al. (1982) PROPOSAL NEAR THE SURFACE.  HIGHER
!     UP MONIN-OBUKOV VALUES ARE USED WITH A LENGTH SCALE OF 
!     300 M (150M).
!     POSSIBILITY TO USE LOUIS et al. (1982) PROPOSAL IN STABLE
!     SITUATIONS STILL KEPT.

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG, RD, RCPD, RETV, RLVTT
USE PARPHY   , ONLY : REPDU2,RKAP
USE YOEVDFS  , ONLY : RCDHALF
USE YOEPHLI  , ONLY : TEPHLI
USE YOEPHY   , ONLY : TEPHY

IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

TYPE(TEPHLI)      ,INTENT(IN)    :: YDEPHLI
TYPE(TEPHY)       ,INTENT(IN)    :: YDEPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KINV(KLON)
INTEGER(KIND=JPIM),INTENT(IN)    :: KCBASE(KLON)
INTEGER(KIND=JPIM),INTENT(IN)    :: KCTOP(KLON)
INTEGER(KIND=JPIM),INTENT(IN)    :: KPBLTYPE(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVDIFTS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0MM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTGZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKMFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKHFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKQFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZINV(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCFM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCFH(KLON,KLEV) 

!*         0.2    LOCAL VARIABLES

INTEGER(KIND=JPIM) :: JK, JL
INTEGER(KIND=JPIM) :: IINV(KLON)

REAL(KIND=JPRB) :: ZRI(KLON),ZDU1(KLON),ZDU2(KLON),ZMGEOM(KLON),ZUST(KLON), &
                 & ZKHVFL(KLON),ZL(KLON),ZLI(KLON)

REAL(KIND=JPRB) :: ZENTRSFC, ZWTVENTR, &
                 & ZCB, ZCD, ZCFNC1, ZDRORO, ZSCF, ZEPS, &
                 & ZZ, ZDH, ZDL, ZETA, ZPHIKH, ZPHIKM, ZPHIKHI, ZPHIKMI,&
                 & ZDHZZ, ZKH, ZREPUST, ZCONS13, ZRG, &
                 & ZIFMOM, ZIFMOH, ZDUDZ, ZDUDZ1, &
                 & ZZH, ZIFLTGM, ZIFLTGH, &
                 & ZDTV, ZDTVI, ZZNORM, ZDENOM, ZFACT2

REAL(KIND=JPRB) :: Z1S, ZCC, ZCONS10
REAL(KIND=JPRB) :: ZKLEN, ZKLEN2
REAL(KIND=JPRB) :: ZKLENT, ZKLENT2
REAL(KIND=JPRB) :: ZSQRT1, ZSQRT2

REAL(KIND=JPRB) :: ZPBLHEIGHT(KLON),ZRIB(KLON),ZSVBOT(KLON),ZRILEV,&
                 & ZRICRI,ZSV
LOGICAL :: LLDONE(KLON)

REAL(KIND=JPHOOK) ::    ZHOOK_HANDLE

#include "surf_inq.h"

!     ------------------------------------------------------------------

!*         1.     INITIALIZE CONSTANTS
!                 --------------------

IF (LHOOK) CALL DR_HOOK('VDFEXCUS',0,ZHOOK_HANDLE)
ASSOCIATE(RLPCC=>YDEPHLI%RLPCC, &
 & YSURF=>YDEPHY%YSURF)
 
! Various constants

ZENTRSFC  = 0.2_JPRB       ! factor for surface based top entrainment 
ZCD       = 1.0_JPRB
ZCB       = 5.0_JPRB
ZEPS      = 1.E-10_JPRB
ZRICRI    = 0.25_JPRB
ZCC       = RLPCC
ZCONS10   = (PVDIFTS*PTMST*RG**2)/(0.5_JPRB*RD)

CALL SURF_INQ(YSURF,PREPUST=ZREPUST)

! optimization
ZRG       = 1.0_JPRB/RG
ZCONS13   = 1.0_JPRB/3._JPRB


!     ------------------------------------------------------------------

!*         2.     PREPARE SCALING COEFFICIENTS FOR MIXED LAYER
!                 --------------------------------------------

DO JL=KIDIA,KFDIA
  ZL(JL) = 1.0_JPRB
  ZUST (JL)=MAX(SQRT(PKMFL(JL)),ZREPUST)
  ZKHVFL(JL)=PKHFL(JL)+RETV*PTM1(JL,KLEV)*PKQFL(JL)
  ! Shifts the inversion level by one as pgeoh defined from 0 to klev, while in vdfh is from 1 to klev+1.
  IINV(JL) = KINV(JL) - 1
  IF (ZKHVFL(JL) < 0.0_JPRB) THEN
    ZL(JL) = ZUST(JL)**3*PTM1(JL,KLEV)/(RKAP*RG*(ZKHVFL(JL)-ZEPS))
  ENDIF
  ZLI(JL) = 1.0_JPRB/ZL(JL)
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
    IF (.NOT. LLDONE(JL) .AND. ZKHVFL(JL)  >  0.0_JPRB) THEN
      ZSV=RCPD*PTM1(JL,JK)*(1.0_JPRB+RETV*PQM1(JL,JK))+PGEOM1(JL,JK)
!   pbl height diag, which considers the winds close to the surf eq to 0
      ZDU2(JL)=MAX(REPDU2, PUM1(JL,JK)**2+PVM1(JL,JK)**2)  
      ZDRORO=(ZSV-ZSVBOT(JL))&
       & /(ZSV-PGEOM1(JL,JK))  
      ZRILEV=(PGEOM1(JL,JK)-PGEOM1(JL,KLEV))*ZDRORO/ZDU2(JL)

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

!***
DO JK = KLEV-1, 1, -1
!***

  DO JL = KIDIA, KFDIA
    PCFM(JL,JK) = 0.0_JPRB
    PCFH(JL,JK) = 0.0_JPRB
  ENDDO

!           COMPUTE RI-NUMBER

  DO JL = KIDIA, KFDIA
    ZMGEOM(JL) = PGEOM1(JL,JK)-PGEOM1(JL,JK+1)
    Z1S = (PUM1(JL,JK)-PUM1(JL,JK+1))**2+(PVM1(JL,JK)-PVM1(JL,JK+1))**2
    ZDU1(JL) = MAX(REPDU2,Z1S)
    ZSQRT1 = SQRT(ZDU1(JL))
    ZDUDZ1 = ZSQRT1/ZMGEOM(JL)*RG              !shear
    ZDU2(JL) = (ZDUDZ1*ZMGEOM(JL)*ZRG) **2

    ZDRORO = 2.0_JPRB*(PCPTGZ(JL,JK)-PCPTGZ(JL,JK+1))/&
     & (PCPTGZ(JL,JK)+PCPTGZ(JL,JK+1)-PGEOM1(JL,JK)-&
     & PGEOM1(JL,JK+1))+&
     & RETV*(PQM1(JL,JK)-PQM1(JL,JK+1))
    ZRI(JL) = ZMGEOM(JL)*ZDRORO/ZDU2(JL)
  ENDDO

  DO JL = KIDIA, KFDIA

    ZSQRT2 = SQRT(ZDU2(JL))
    ZDUDZ = ZSQRT2/ZMGEOM(JL)*RG    !shear
    ZZH = 0.5_JPRB * ZRG * (PGEOM1(JL,JK)+PGEOM1(JL,JK+1)) + PZ0MM(JL)

!          COMPUTE EXCHANGE COEFFICIENTS FOR MOMENTUM AND HEAT EXCHANGE 
!          IN STABLE OR UNSTABLE REGIME (BASED ON SIGN OF RICHARDSON NUMBER).

    IF (ZRI(JL) > 0.0_JPRB) THEN  ! statically stable

!     ASYMPTOTIC MIXING LENGTH FOR STABLE SITUATIONS
      IF (PGEOM1(JL,JK)*ZRG <= ZPBLHEIGHT(JL)) THEN
        ZKLEN=120.0_JPRB
        ZKLENT=60.0_JPRB
      ELSE
        ZKLEN =60.0_JPRB
        ZKLENT=60.0_JPRB
      ENDIF
      ZKLEN2  = RKAP * ZZH * ZKLEN  / ( RKAP * ZZH + ZKLEN )
      ZKLENT2 = RKAP * ZZH * ZKLENT / ( RKAP * ZZH + ZKLENT )

!          COMPUTE STABILITY FUNCTIONS
      ZSCF    = SQRT(1.0_JPRB+ZCD*ZRI(JL))
!NL      ZIFLTGM = 1.0_JPRB / (1.0_JPRB + 2.0_JPRB *ZCB * ZRI(JL)/ZSCF) !F(LTG),M
      ZIFLTGM = ZSCF/(ZSCF+2.0_JPRB*ZCB*ZRI(JL))                     !F(LTG),M
      ZIFLTGH = 1.0_JPRB / (1.0_JPRB + 2.0_JPRB *ZCB * ZRI(JL)*ZSCF) !F(LTG),H

      PCFM(JL,JK) = ZDUDZ * ZKLEN2**2  * ZIFLTGM
      PCFH(JL,JK) = ZDUDZ * ZKLENT2**2 * ZIFLTGH

    ELSE  ! statically unstable

      ZKLENT=150.0_JPRB
      ZKLEN2 = RKAP * ZZH * ZKLENT / ( RKAP * ZZH + ZKLENT )

!          COMPUTE STABILITY FUNCTIONS
      ZETA = ZRI(JL)
      ZIFMOM = SQRT(1.0_JPRB-RCDHALF*ZETA)              !F(MO),M
      ZIFMOH = (1.0_JPRB-RCDHALF*ZETA)**0.75_JPRB       !F(MO),H
      PCFM(JL,JK) = ZDUDZ * ZKLEN2**2 * ZIFMOM
      PCFH(JL,JK) = ZDUDZ * ZKLEN2**2 * ZIFMOH

    ENDIF

!           ADD MIXED LAYER PARAMETRIZATION (SURF. AND OUTER LAYER)

    IF (JK > IINV(JL) .AND. IINV(JL) < (KLEV-1)) THEN  ! up to level below entr.lev.

      ZZ     = PGEOH(JL,JK)*ZRG
      ZDH  = ZZ /(PGEOH(JL,IINV(JL))*ZRG)
      ZDL    = ZZ*ZLI(JL)
      ZPHIKH = (1.0_JPRB-39._JPRB*ZDL)**(-ZCONS13)
      ZPHIKM = (1.0_JPRB-15._JPRB*ZDL)**(-ZCONS13)

!           K,surface

      ZZNORM = ZZ * (1.0_JPRB - ZDH)**2
      IF (KPBLTYPE(JL) >= 2 .AND. KCBASE(JL) > KCTOP(JL)+1) THEN
        IF (RCPD*PKHFL(JL) < -0.1_JPRB .AND. -RCPD*PKHFL(JL) < -0.25_JPRB*RLVTT*PKQFL(JL)) THEN
           ZZNORM = ZZ * (1.0_JPRB - ZDH)**4
           ZDH = MIN(1.0_JPRB, ZZ/(PGEOH(JL,KCBASE(JL)-1)*ZRG))
           ZZNORM = ZZNORM + ZZ * (1.0_JPRB - ZDH)**2 
        ENDIF
      ENDIF
      ZDHZZ = ZUST(JL) * ZZNORM

      ZPHIKHI = RKAP / ZPHIKH
      PCFH(JL,JK) = ZPHIKHI * ZDHZZ
      ZPHIKMI = RKAP / ZPHIKM
      PCFM(JL,JK) = ZPHIKMI * ZDHZZ

    ENDIF

!           ADD ENTRAINMENT
!           entrainment velocity * T,v jump due to sfc flux

    IF (JK == IINV(JL) .AND. IINV(JL) < (KLEV-1)) THEN

      ZWTVENTR = - ZENTRSFC * ZKHVFL(JL)              !sfc flux
      ZWTVENTR = MAX(0.0_JPRB,ZWTVENTR)
      
      ZDTV = (PCPTGZ(JL,JK-1) - PCPTGZ(JL,JK)) / RCPD &
         & +  RETV * 0.5_JPRB * (PQM1(JL,JK-1) - PQM1(JL,JK)) &
         & *  (PTM1(JL,JK) + PTM1(JL,JK-1))

      IF (ZDTV > 1.E-3_JPRB) THEN
        ZDTVI = 1.0_JPRB / ZDTV
        ZKH = ZWTVENTR * ZRG * (PGEOM1(JL,JK-1) - PGEOM1(JL,JK)) * ZDTVI
      ELSE
        ZKH = 0.0_JPRB
      ENDIF

      PCFH(JL,JK) = ZKH
      PCFM(JL,JK) = ZKH
    ENDIF

!           SCALE DIFFUSION COEFFICIENTS FOR USE IN VDFDIFH/M

    ZDENOM = ZMGEOM(JL) &
          & * (PTM1(JL,JK)   * (1.0_JPRB + RETV * PQM1(JL,JK)) &
          & +  PTM1(JL,JK+1) * (1.0_JPRB + RETV * PQM1(JL,JK+1)))
    ZFACT2 = 1.0_JPRB / ZDENOM
    ZCFNC1 = ZCONS10 * PAPHM1(JL,JK) * ZFACT2
    PCFH(JL,JK) = PCFH(JL,JK) * ZCFNC1
    PCFM(JL,JK) = PCFM(JL,JK) * ZCFNC1

  ENDDO

!***
ENDDO
!***

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('VDFEXCUS',1,ZHOOK_HANDLE)
END SUBROUTINE VDFEXCUS
