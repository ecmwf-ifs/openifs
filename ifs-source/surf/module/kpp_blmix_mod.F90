MODULE KPP_BLMIX_MOD
CONTAINS
SUBROUTINE KPP_BLMIX &
  & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KLEVOP1  ,&
  &   LDKPPCAL ,PZO      ,PHO      ,PHO_INV  ,PUSTAR   ,&
  &   PBFSFC   ,PHBL     ,PSTABLE  ,PCASEA   ,PDIFM    ,&
  &   PDIFS    ,PDIFT    ,KBL      ,PGHAT    ,PSIGMA   ,&
  &   PWM      ,PWS      ,YDOCEAN_ML)

! (C) Copyright 1994- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! Purpose :
! -------
!   This routine sets computes mixing coefficients within boundary 
!   layer depend on surface forcing and the magnitude and gradient 
!   of interior mixing below the boundary layer ("matching").

! Interface :
! ---------
!   Call *KPP_BLMIX* from *KPP_KPPMIX*

! Method :
! ------

! Externals :
! ---------

! Reference :
! ---------
! Large, W. G., J. C. McWilliams, S. C. Doney (1994), Rev. Geophys.

! Modifications :
! -------------
!     06-Jun-1994  Bill Large
!            2002  Steve Woolnough, Reading Univ.
!     07-Oct-2008  Yuhei Takaya,    E.C.M.W.F.    Implemented to IFS.
! End Modifications :
!---------------------------------------------------------------------

USE PARKIND1,     ONLY : JPIM, JPRB
USE YOMHOOK,      ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_OCEAN_ML, ONLY : TOCEAN_ML

USE KPP_WSCALE_MOD

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVO
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVOP1
INTEGER(KIND=JPIM),INTENT(IN) :: KBL(KLON)         ! idx of 1st lev. below PHBL
REAL(KIND=JPRB),INTENT(IN)    :: PZO(KLON,KLEVOP1) ! vertical grid (<=0) d  (m)
REAL(KIND=JPRB),INTENT(IN)    :: PHO(KLON,KLEVOP1) ! layer thicknesses      (m)
REAL(KIND=JPRB),INTENT(IN)    :: PHO_INV(KLON,KLEVOP1) ! 1/PHO
REAL(KIND=JPRB),INTENT(IN)    :: PUSTAR(KLON)      ! sfc friction vel. (m/s)
REAL(KIND=JPRB),INTENT(IN)    :: PBFSFC(KLON)      ! sfc buoy forcing (m^2/s^3)
REAL(KIND=JPRB),INTENT(IN)    :: PHBL(KLON)        ! boundary layer depth   (m)
REAL(KIND=JPRB),INTENT(IN)    :: PSTABLE(KIDIA:KFDIA) ! = 1 in stable forcing
REAL(KIND=JPRB),INTENT(IN)    :: PCASEA(KLON)      ! = 1 in case A
REAL(KIND=JPRB),INTENT(OUT)   :: PGHAT(KLON,KLEVO) ! nonlocal scalar transport
REAL(KIND=JPRB),INTENT(OUT)   :: PSIGMA(KLON)      ! normalized depth (d / hbl)
REAL(KIND=JPRB),INTENT(OUT)   :: PWS(KIDIA:KFDIA)  ! turb. scalar scales (m/s)
REAL(KIND=JPRB),INTENT(OUT)   :: PWM(KIDIA:KFDIA)  ! turb. vel. scales (m/s)
REAL(KIND=JPRB),INTENT(INOUT) :: PDIFM(KLON,0:KLEVO) ! visc. coef. (m^2/s)
REAL(KIND=JPRB),INTENT(INOUT) :: PDIFS(KLON,0:KLEVO) ! scalar diff.(m^2/s)
REAL(KIND=JPRB),INTENT(INOUT) :: PDIFT(KLON,0:KLEVO) ! temp. diff. (m^2/s)

LOGICAL,INTENT(IN) :: LDKPPCAL(KLON) 

TYPE(TOCEAN_ML),INTENT(INOUT) :: YDOCEAN_ML

INTEGER(KIND=JPIM) :: JZ
INTEGER(KIND=JPIM) :: JL
INTEGER(KIND=JPIM) :: IN
INTEGER(KIND=JPIM) :: INP1
INTEGER(KIND=JPIM) :: ILEVMAX

REAL(KIND=JPRB) :: ZDKM1(KIDIA:KFDIA,YDOCEAN_ML%NDIFFO)   ! boundary layer difs at kbl-1
REAL(KIND=JPRB) :: ZBLMC(KIDIA:KFDIA,KLEVO,YDOCEAN_ML%NDIFFO) ! BL mixing coef. (m^2/s)
REAL(KIND=JPRB) :: ZGAT1(KIDIA:KFDIA,YDOCEAN_ML%NDIFFO)
REAL(KIND=JPRB) :: ZDAT1(KIDIA:KFDIA,YDOCEAN_ML%NDIFFO)
REAL(KIND=JPRB) :: ZCG
REAL(KIND=JPRB) :: ZR(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZDVDZDN(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZDVDZUP(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZDIFMH(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZDIFTH(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZDIFSH(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZDIFMP(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZDIFTP(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZDIFSP(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZDELHAT(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZF1(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZSIG(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZGM(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZGT(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZGS(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZA1(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZA2(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZA3(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZDEL(KIDIA:KFDIA) 
REAL(KIND=JPRB) :: ZDEL2(KIDIA:KFDIA) 
REAL(KIND=JPRB) :: ZDKMP5(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZDSTAR(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZHBL_INV(KIDIA:KFDIA)
REAL(KIND=JPRB) :: Z1CASEA(KIDIA:KFDIA)
REAL(KIND=JPRB) :: Z1DEL(KIDIA:KFDIA)
REAL(KIND=JPRB) :: Z1DEL2(KIDIA:KFDIA)
REAL(KIND=JPRB) :: Z1STABLE(KIDIA:KFDIA) 
REAL(KIND=JPRB) :: ZWM_INV(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZWS_INV(KIDIA:KFDIA)
REAL(KIND=JPRB) :: Z1R(KIDIA:KFDIA)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('KPP_BLMIX_MOD:KPP_BLMIX',0,ZHOOK_HANDLE)
ASSOCIATE(NDIFFO=>YDOCEAN_ML%NDIFFO, RC1_KPP=>YDOCEAN_ML%RC1_KPP, &
 & RCSTAR_KPP=>YDOCEAN_ML%RCSTAR_KPP, RCS_KPP=>YDOCEAN_ML%RCS_KPP, &
 & REPSILON_KPP=>YDOCEAN_ML%REPSILON_KPP, REPS_KPP=>YDOCEAN_ML%REPS_KPP, &
 & RVONK=>YDOCEAN_ML%RVONK)

ZCG = RCSTAR_KPP * RVONK &
    & * (RCS_KPP * RVONK * REPSILON_KPP)**(1.0_JPRB/3.0_JPRB)
 

! 1. Compute nondimensional shape function
!    -------------------------------------

! 1.1 Compute velocity scales at PHBL
!     -------------------------------

DO JL = KIDIA, KFDIA
  IF( LDKPPCAL(JL) ) THEN
    Z1STABLE(JL) = 1.0_JPRB - PSTABLE(JL)
    PSIGMA(JL) = PSTABLE(JL) + Z1STABLE(JL) * REPSILON_KPP
  ENDIF
ENDDO

CALL KPP_WSCALE &
 & ( KIDIA    ,KFDIA    ,KLON     ,LDKPPCAL ,PSIGMA   ,&
 &   PHBL     ,PUSTAR   ,PBFSFC   ,PWM      ,PWS      ,&
 &   YDOCEAN_ML)

! 1.2 find the interior viscosities and derivatives at PHBL 
!     -----------------------------------------------------

ILEVMAX=0

DO JL = KIDIA, KFDIA
  IF( LDKPPCAL(JL) ) THEN

    IN = INT( PCASEA(JL) + REPS_KPP )
    IN = IN* ( KBL(JL) - 1 ) + ( 1 - IN ) * KBL(JL)
    INP1 = MIN( IN+1, KLEVO )  !used for diff.

    ZWM_INV(JL) = 1.0_JPRB /  ( PWM(JL) + REPS_KPP )
    ZWS_INV(JL) = 1.0_JPRB /  ( PWS(JL) + REPS_KPP )

    ZDELHAT(JL) = 0.5_JPRB * PHO(JL,IN) - PZO(JL,IN) - PHBL(JL)
    ZR(JL)      = 0.5_JPRB * ( 1.0_JPRB - ZDELHAT(JL) * PHO_INV(JL,IN) )
    Z1R(JL) =  0.5_JPRB * ( 1.0_JPRB - ZR(JL) )
    ZDVDZUP(JL) = ( PDIFM(JL,IN-1) - PDIFM(JL,IN) ) * PHO_INV(JL,IN) 
    ZDVDZDN(JL) = ( PDIFM(JL,IN) - PDIFM(JL,INP1) )  * PHO_INV(JL,IN+1)
    ZDIFMP(JL)  = ( Z1R(JL)*( ZDVDZUP(JL) + ABS(ZDVDZUP(JL)) ) &
                & + ZR(JL) * ( ZDVDZDN(JL) + ABS(ZDVDZDN(JL) ) ) )
    ZDVDZUP(JL) = ( PDIFS(JL,IN-1) - PDIFS(JL,IN) ) * PHO_INV(JL,IN)
    ZDVDZDN(JL) = ( PDIFS(JL,IN)   - PDIFS(JL,INP1) )  * PHO_INV(JL,IN+1)
    ZDIFSP(JL)  =  ( Z1R(JL)*( ZDVDZUP(JL) + ABS(ZDVDZUP(JL)) ) &
                & + ZR(JL) * ( ZDVDZDN(JL) + ABS(ZDVDZDN(JL)) ) )
    ZDVDZUP(JL) = ( PDIFT(JL,IN-1) - PDIFT(JL,IN) ) * PHO_INV(JL,IN)
    ZDVDZDN(JL) = ( PDIFT(JL,IN)   - PDIFT(JL,INP1) ) * PHO_INV(JL,IN+1)
    ZDIFTP(JL)  = ( Z1R(JL)*( ZDVDZUP(JL) + ABS(ZDVDZUP(JL)) ) &
                & + ZR(JL) * ( ZDVDZDN(JL) + ABS(ZDVDZDN(JL)) ) )
    ZDIFMH(JL)  = PDIFM(JL,IN) + ZDIFMP(JL) * ZDELHAT(JL)
    ZDIFSH(JL)  = PDIFS(JL,IN) + ZDIFSP(JL) * ZDELHAT(JL)
    ZDIFTH(JL)  = PDIFT(JL,IN) + ZDIFTP(JL) * ZDELHAT(JL)

    ZHBL_INV(JL) = 1.0_JPRB / PHBL(JL) 
    ZF1(JL) = PSTABLE(JL) * RC1_KPP * PBFSFC(JL) &
            & / ( PUSTAR(JL)**4 + REPS_KPP ) 
    ZGAT1(JL,1) = ZDIFMH(JL) * ZHBL_INV(JL) * ZWM_INV(JL)
    ZDAT1(JL,1) = -ZDIFMP(JL) *  ZWM_INV(JL) &
                & + ZF1(JL) * ZDIFMH(JL)
    ZDAT1(JL,1) = MIN( ZDAT1(JL,1), 0.0_JPRB ) 
    ZGAT1(JL,2) = ZDIFSH(JL)  * ZHBL_INV(JL) * ZWS_INV(JL) 
    ZDAT1(JL,2) = -ZDIFSP(JL) * ZWS_INV(JL) &
                & + ZF1(JL) * ZDIFSH(JL) 
    ZDAT1(JL,2) = MIN( ZDAT1(JL,2), 0.0_JPRB ) 
    ZGAT1(JL,3) = ZDIFTH(JL) * ZHBL_INV(JL) * ZWS_INV(JL)
    ZDAT1(JL,3) = -ZDIFTP(JL) * ZWS_INV(JL) &
                & + ZF1(JL) * ZDIFTH(JL) 
    ZDAT1(JL,3) = MIN( ZDAT1(JL,3), 0.0_JPRB ) 

    ILEVMAX = MAX(ILEVMAX,KBL(JL)-1) !for the next loop 

  ENDIF
ENDDO

DO JZ = 1, ILEVMAX 

! 1.3 Compute turbulent velocity scales on the interfaces
!     ---------------------------------------------------

  DO JL = KIDIA, KFDIA
    IF( LDKPPCAL(JL) ) THEN
      ZSIG(JL) = ( -PZO(JL,JZ) + 0.5_JPRB * PHO(JL,JZ) ) * ZHBL_INV(JL)
      PSIGMA(JL) = PSTABLE(JL) * ZSIG(JL) &
                 & + Z1STABLE(JL) * MIN( ZSIG(JL), REPSILON_KPP )
    ENDIF
  ENDDO

  CALL KPP_WSCALE &
   & ( KIDIA    ,KFDIA    ,KLON     ,LDKPPCAL ,PSIGMA   ,&
   &   PHBL     ,PUSTAR   ,PBFSFC   ,PWM      ,PWS      ,&
   &   YDOCEAN_ML)

! 1.4 Compute the dimensionless shape functions at the interfaces
!     -----------------------------------------------------------

  DO JL = KIDIA, KFDIA
    IF( LDKPPCAL(JL) ) THEN

      ZSIG(JL) = ( -PZO(JL,JZ) + 0.5_JPRB * PHO(JL,JZ) ) * ZHBL_INV(JL)
      ZA1(JL) = ZSIG(JL) - 2.0_JPRB
      ZA2(JL) = 3.0_JPRB - 2.0_JPRB * ZSIG(JL)
      ZA3(JL) = ZSIG(JL) - 1.0_JPRB
      ZGM(JL) = ZA1(JL) + ZA2(JL) * ZGAT1(JL,1) + ZA3(JL) * ZDAT1(JL,1) 
      ZGS(JL) = ZA1(JL) + ZA2(JL) * ZGAT1(JL,2) + ZA3(JL) * ZDAT1(JL,2)
      ZGT(JL) = ZA1(JL) + ZA2(JL) * ZGAT1(JL,3) + ZA3(JL) * ZDAT1(JL,3)

!     compute boundary layer diffusivities at the interfaces
      ZBLMC(JL,JZ,1) = PHBL(JL) * PWM(JL) * ZSIG(JL) &
                     & * ( 1.0_JPRB + ZSIG(JL) * ZGM(JL) )
      ZBLMC(JL,JZ,2) = PHBL(JL) * PWS(JL) * ZSIG(JL) &
                     & * ( 1.0_JPRB + ZSIG(JL) * ZGS(JL) )
      ZBLMC(JL,JZ,3) = PHBL(JL) * PWS(JL) * ZSIG(JL) &
                     & * ( 1.0_JPRB + ZSIG(JL) * ZGT(JL) )

!     nonlocal transport term = ghats * <ws>o
      PGHAT(JL,JZ) = Z1STABLE(JL) * ZCG &
                   & / ( PWS(JL) * PHBL(JL) + REPS_KPP )

    ENDIF
  ENDDO
ENDDO !JZ

! 1.5 Find diffusivities at KBL-1 grid level 
!     --------------------------------------

DO JL = KIDIA, KFDIA
  IF( LDKPPCAL(JL) ) THEN
    ZSIG(JL)   = -PZO(JL,KBL(JL)-1) * ZHBL_INV(JL)
    PSIGMA(JL) = PSTABLE(JL) * ZSIG(JL) &
               & + Z1STABLE(JL) * MIN(ZSIG(JL), REPSILON_KPP)
  ENDIF
ENDDO

CALL KPP_WSCALE &
 & ( KIDIA    ,KFDIA    ,KLON     ,LDKPPCAL ,PSIGMA   ,&
 &   PHBL     ,PUSTAR   ,PBFSFC   ,PWM      ,PWS      ,&
 &   YDOCEAN_ML)

DO JL = KIDIA, KFDIA
  IF( LDKPPCAL(JL) ) THEN
    ZSIG(JL) = -PZO(JL,KBL(JL)-1) * ZHBL_INV(JL)
    ZA1(JL) = ZSIG(JL) - 2.0_JPRB
    ZA2(JL) = 3.0_JPRB - 2.0_JPRB * ZSIG(JL)
    ZA3(JL) = ZSIG(JL) - 1.0_JPRB
    ZGM(JL) = ZA1(JL) + ZA2(JL) * ZGAT1(JL,1) + ZA3(JL) * ZDAT1(JL,1)
    ZGS(JL) = ZA1(JL) + ZA2(JL) * ZGAT1(JL,2) + ZA3(JL) * ZDAT1(JL,2)
    ZGT(JL) = ZA1(JL) + ZA2(JL) * ZGAT1(JL,3) + ZA3(JL) * ZDAT1(JL,3)
    ZDKM1(JL,1) = PHBL(JL) * PWM(JL) * ZSIG(JL) &
                & * ( 1.0_JPRB + ZSIG(JL) * ZGM(JL) )
    ZDKM1(JL,2) = PHBL(JL) * PWS(JL) * ZSIG(JL) &
                & * ( 1.0_JPRB + ZSIG(JL) * ZGS(JL) )
    ZDKM1(JL,3) = PHBL(JL) * PWS(JL) * ZSIG(JL) &
                & * ( 1.0_JPRB + ZSIG(JL) * ZGT(JL) )
  ENDIF
ENDDO

! 2. Interior-Boundary Layer Matching
!    --------------------------------

! subroutine ENHANCE is inserted here

DO JL = KIDIA, KFDIA
  IF( LDKPPCAL(JL) ) THEN
!    DO JZ = 1, KLEVO-1
!      IF( JZ  ==  (KBL(JL) - 1) ) THEN
        JZ = KBL(JL) - 1_JPIM
        ZDEL(JL) = ( PHBL(JL) + PZO(JL,JZ) ) / ( PZO(JL,JZ) - PZO(JL,JZ+1) )
        ZDEL2(JL) = ZDEL(JL) * ZDEL(JL)
        Z1DEL(JL) = 1.0_JPRB - ZDEL(JL) 
        Z1DEL2(JL) = Z1DEL(JL) * Z1DEL(JL)
        Z1CASEA(JL)= 1.0_JPRB - PCASEA(JL)
        ZDKMP5(JL) = PCASEA(JL) * PDIFM(JL,JZ) &
                   & + Z1CASEA(JL) * ZBLMC(JL,JZ,1)
        ZDSTAR(JL) = Z1DEL2(JL) * ZDKM1(JL,1) &
                   & + ZDEL2(JL) * ZDKMP5(JL)
        ZBLMC(JL,JZ,1) = Z1DEL(JL) * PDIFM(JL,JZ) &
                       & + ZDEL(JL) * ZDSTAR(JL)
        ZDKMP5(JL) = PCASEA(JL) * PDIFS(JL,JZ) &
                   & + Z1CASEA(JL) * ZBLMC(JL,JZ,2)
        ZDSTAR(JL) = Z1DEL2(JL) * ZDKM1(JL,2) &
                   & + ZDEL2(JL) * ZDKMP5(JL)
        ZBLMC(JL,JZ,2) = Z1DEL(JL) * PDIFS(JL,JZ) &
                       & + ZDEL(JL) * ZDSTAR(JL)
        ZDKMP5(JL) = PCASEA(JL) * PDIFT(JL,JZ) &
                   & + Z1CASEA(JL) * ZBLMC(JL,JZ,3)
        ZDSTAR(JL) = Z1DEL2(JL) * ZDKM1(JL,3) &
                   & + ZDEL2(JL) * ZDKMP5(JL)
        ZBLMC(JL,JZ,3) = Z1DEL(JL) * PDIFT(JL,JZ) &
                       & + ZDEL(JL) * ZDSTAR(JL)
        PGHAT(JL,JZ) = Z1CASEA(JL) * PGHAT(JL,JZ)
!      ENDIF
!    ENDDO

! 3. Set boundery layer diffusivity
!    ------------------------------
    DO JZ = 1, KLEVO
      IF( JZ < KBL(JL) ) THEN
!        PDIFM(JL,JZ) = ZBLMC(JL,JZ,1)
!        PDIFS(JL,JZ) = ZBLMC(JL,JZ,2)
!        PDIFT(JL,JZ) = ZBLMC(JL,JZ,3)
         PDIFM(JL,JZ) = MAX(ZBLMC(JL,JZ,1),PDIFM(JL,JZ))  !Y.T. 05/04/2009
         PDIFS(JL,JZ) = MAX(ZBLMC(JL,JZ,2),PDIFS(JL,JZ))
         PDIFT(JL,JZ) = MAX(ZBLMC(JL,JZ,3),PDIFT(JL,JZ))
      ELSE
        PGHAT(JL,JZ) = 0.0_JPRB
      ENDIF
    ENDDO

  ENDIF
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('KPP_BLMIX_MOD:KPP_BLMIX',1,ZHOOK_HANDLE)
END SUBROUTINE KPP_BLMIX
END MODULE KPP_BLMIX_MOD
