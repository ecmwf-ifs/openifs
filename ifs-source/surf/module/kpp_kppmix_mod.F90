MODULE KPP_KPPMIX_MOD
CONTAINS
SUBROUTINE KPP_KPPMIX &
  & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KLEVOP1  ,&
  &   LDINIT_KPP, LDKPPCAL ,LDSICE ,PTAUX    ,PTAUY    ,&
  &   PSWF     ,PLWF     ,PLHF     ,PSHF     ,PRAIN    ,&
  &   PSNOW    ,PZO      ,PHO      ,PHO_INV  ,PDO      ,&
  &   PSFLUX   ,PUO      ,PXO      ,PHBL     ,KBL      ,&
  &   PF       ,PDIFM    ,PDIFS    ,PDIFT    ,PGHAT    ,&
  &   PRHO     ,PCP      ,PWU      ,PWX      ,PWXNT    ,&
  &   PSWDK_SAVE, PTALPHA ,PSBETA  ,POCDEPTH ,PRHOH2O  ,&
  &   PRHO_INV ,YDCST    ,YDOCEAN_ML )

! (C) Copyright 1994- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! Purpose :
! -------
!   This routine computes coefficients for vertical mixing.

! Interface :
! ---------

! Method :
! ------
!   K-Profile parameterization

! Externals :
! ---------

! Reference :
! ---------
!   Large, W. G., J. C. McWilliams, S. C. Doney (1994), Rev. Geophys.

! Modifications :
! -------------
!     06-Jun-1994  Bill Large
!            2002  Steve Woolnough, Reading Univ.
!     07-Oct-2008  Yuhei Takaya,    E.C.M.W.F.    Implemented to IFS.
! End Modifications :
!---------------------------------------------------------------------


USE PARKIND1,     ONLY : JPIM, JPRB
USE YOMHOOK,      ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_CST,      ONLY : TCST
USE YOS_OCEAN_ML, ONLY : TOCEAN_ML, NVELO, NSCLRO

USE KPP_ABK80_MOD
USE KPP_CPSW_MOD
USE KPP_SWFRAC_MOD
USE KPP_INTERIOR_MIX_MOD
USE KPP_BLDEPTH_MOD
USE KPP_BLMIX_MOD

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVO
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVOP1
INTEGER(KIND=JPIM),INTENT(INOUT) :: KBL(KLON)
REAL(KIND=JPRB),INTENT(INOUT) :: PTAUX(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PTAUY(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PSWF(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PLWF(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PLHF(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PSHF(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PRAIN(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PSNOW(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PZO(KLON,KLEVOP1)
REAL(KIND=JPRB),INTENT(IN)    :: PDO(KLON,0:KLEVO)
REAL(KIND=JPRB),INTENT(IN)    :: PHO(KLON,KLEVOP1)
REAL(KIND=JPRB),INTENT(IN)    :: PHO_INV(KLON,KLEVOP1)
REAL(KIND=JPRB),INTENT(IN)    :: PF(KIDIA:KFDIA)
REAL(KIND=JPRB),INTENT(IN)    :: PUO(KIDIA:KFDIA,KLEVOP1,NVELO)
REAL(KIND=JPRB),INTENT(IN)    :: PXO(KIDIA:KFDIA,KLEVOP1,NSCLRO)
REAL(KIND=JPRB),INTENT(IN)    :: POCDEPTH(KLON) ! ocean depth
REAL(KIND=JPRB),INTENT(INOUT) :: PHBL(KLON) 
REAL(KIND=JPRB),INTENT(INOUT) :: PWX(KLON,0:KLEVO,NSCLRO+1)
REAL(KIND=JPRB),INTENT(INOUT) :: PWXNT(KLON,0:KLEVO,NSCLRO)
REAL(KIND=JPRB),INTENT(INOUT) :: PWU(KLON,0:KLEVO,NVELO)
REAL(KIND=JPRB),INTENT(INOUT) :: PSWDK_SAVE(KLON,0:KLEVO)
REAL(KIND=JPRB),INTENT(OUT)   :: PSFLUX(KLON,6)
REAL(KIND=JPRB),INTENT(OUT)   :: PRHO(KLON,0:KLEVOP1)
REAL(KIND=JPRB),INTENT(OUT)   :: PRHO_INV(KIDIA:KFDIA)
REAL(KIND=JPRB),INTENT(OUT)   :: PCP(KLON,0:KLEVOP1)
REAL(KIND=JPRB),INTENT(OUT)   :: PDIFM(KLON,0:KLEVO)
REAL(KIND=JPRB),INTENT(OUT)   :: PDIFS(KLON,0:KLEVO)
REAL(KIND=JPRB),INTENT(OUT)   :: PDIFT(KLON,0:KLEVO)
REAL(KIND=JPRB),INTENT(OUT)   :: PGHAT(KLON,KLEVO)
REAL(KIND=JPRB),INTENT(OUT)   :: PRHOH2O(KLON)
REAL(KIND=JPRB),INTENT(OUT)   :: PTALPHA(KIDIA:KFDIA,0:KLEVOP1)
REAL(KIND=JPRB),INTENT(OUT)   :: PSBETA(KIDIA:KFDIA,0:KLEVOP1)

LOGICAL,INTENT(IN) :: LDKPPCAL(KLON)
LOGICAL,INTENT(IN) :: LDSICE(KLON)
LOGICAL,INTENT(IN) :: LDINIT_KPP

TYPE(TCST),     INTENT(IN)    :: YDCST
TYPE(TOCEAN_ML),INTENT(INOUT) :: YDOCEAN_ML

INTEGER(KIND=JPIM) :: JL
INTEGER(KIND=JPIM) :: JZ
INTEGER(KIND=JPIM) :: JS
INTEGER(KIND=JPIM) :: JK
INTEGER(KIND=JPIM) :: JKL

REAL(KIND=JPRB) :: ZRHOB(KLON)
REAL(KIND=JPRB) :: ZSHSQ(KIDIA:KFDIA,KLEVO) ! (local velocity shear)^2 (m/s)^2
REAL(KIND=JPRB) :: ZDVSQ(KIDIA:KFDIA,KLEVO) ! (velocity shear re sfc)^2 (m/s)^2
REAL(KIND=JPRB) :: ZDBLOC(KIDIA:KFDIA,KLEVO)! local delta buoyancy     (m/s^2)
REAL(KIND=JPRB) :: ZRITOP(KIDIA:KFDIA,KLEVO)! numerator of bulk Ri     (m/s)^2
REAL(KIND=JPRB) :: ZB0(KIDIA:KFDIA)   
REAL(KIND=JPRB) :: ZB0SOL(KIDIA:KFDIA) 
REAL(KIND=JPRB) :: ZBFSFC(KLON)  ! bo+radiation absorbed to d=hbf*hbl
REAL(KIND=JPRB) :: ZSTABLE(KIDIA:KFDIA) ! =1 in stable forcing
                                        ! =0 unstable
REAL(KIND=JPRB) :: ZCASEA(KLON)         ! =1 in stable forcing
REAL(KIND=JPRB) :: ZSIGMA(KLON)         ! sigma for profile function
REAL(KIND=JPRB) :: ZWM(KIDIA:KFDIA)     ! momentum velocity scale
REAL(KIND=JPRB) :: ZWS(KIDIA:KFDIA)     ! scalar velocity scale
REAL(KIND=JPRB) :: ZSWDK(KLON)
REAL(KIND=JPRB) :: ZBUOY(KIDIA:KFDIA,KLEVOP1)
REAL(KIND=JPRB) :: ZALPHADT(KIDIA:KFDIA,KLEVO)   ! alpha * dt across interfaces
REAL(KIND=JPRB) :: ZBETADS(KIDIA:KFDIA,KLEVO)    ! beta * ds across interfaces
REAL(KIND=JPRB) :: ZDUO(KIDIA:KFDIA,KLEVO,NVELO) ! PUO(JZ)-PUO(JZ+1)
REAL(KIND=JPRB) :: ZDXO(KIDIA:KFDIA,KLEVO,NSCLRO)! PXO(JZ)-PXO(JZ+1)
REAL(KIND=JPRB) :: ZDZO(KIDIA:KFDIA,KLEVO)       ! PZO(JZ)-PZO(JZ+1)
REAL(KIND=JPRB) :: ZALPHA(KLON)            
REAL(KIND=JPRB) :: ZBETA(KLON)
REAL(KIND=JPRB) :: ZSIG(KLON,KLEVOP1)
REAL(KIND=JPRB) :: ZSIG0(KLON,KLEVOP1)
REAL(KIND=JPRB) :: ZTAU(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZUSTAR(KLON)
REAL(KIND=JPRB) :: ZREF(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZWZ(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZBREF(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZUREF(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZVREF(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZDEL(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZDLIMIT
REAL(KIND=JPRB) :: ZVLIMIT
REAL(KIND=JPRB) :: ZTMP1(KIDIA:KFDIA,KLEVOP1)
REAL(KIND=JPRB) :: ZTMP2(KIDIA:KFDIA,KLEVOP1)
REAL(KIND=JPRB) :: ZTMP3(KLON)
REAL(KIND=JPRB) :: Z1
REAL(KIND=JPRB) :: Z0

REAL(KIND=JPRB) :: Z1000
REAL(KIND=JPRB) :: ZRG1000
REAL(KIND=JPRB) :: ZRLMLT_INV
REAL(KIND=JPRB) :: ZRLVTT_INV
REAL(KIND=JPRB) :: ZRRAMSICE_INV
REAL(KIND=JPRB) :: ZCP_INV(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZREF_INV(KIDIA:KFDIA)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

LOGICAL :: LLALPHA
LOGICAL :: LLBETA

!---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('KPP_KPPMIX_MOD:KPP_KPPMIX',0,ZHOOK_HANDLE)
ASSOCIATE(RG=>YDCST%RG, RLMLT=>YDCST%RLMLT, RLVTT=>YDCST%RLVTT, &
 & LKPP_KPP=>YDOCEAN_ML%LKPP_KPP, LNBFLX_KPP=>YDOCEAN_ML%LNBFLX_KPP, &
 & RDIFM_BOT=>YDOCEAN_ML%RDIFM_BOT, RDIFS_BOT=>YDOCEAN_ML%RDIFS_BOT, &
 & REPSILON_KPP=>YDOCEAN_ML%REPSILON_KPP, REPS_KPP=>YDOCEAN_ML%REPS_KPP, &
 & RRAMSICE=>YDOCEAN_ML%RRAMSICE, RSICE=>YDOCEAN_ML%RSICE, &
 & RTICE=>YDOCEAN_ML%RTICE)

Z1=1.0_JPRB
Z0=0.0_JPRB
Z1000=1000.0_JPRB
ZRG1000= RG / 1000.0_JPRB
ZRLMLT_INV = Z1 / RLMLT
ZRLVTT_INV = Z1 / RLVTT
ZRRAMSICE_INV = Z1 / RRAMSICE


! Reset the mixing coefficients
DO JL = KIDIA, KFDIA
  IF( LDKPPCAL(JL) ) THEN
    DO JZ=1,KLEVO
      PGHAT(JL,JZ) = Z0
    ENDDO

! Reset fluxes
    PWX(JL,:,:)   = Z0
    PWU(JL,:,:)   = Z0
    DO JS=2,NSCLRO
      PWXNT(JL,:,JS) = Z0
    ENDDO

    DO JZ = 1, KLEVOP1
      ZTMP2(JL,JZ) = -PZO(JL,JZ) !pressure  for next KPP_ABK80
    ENDDO
    ZTMP1(JL,1) = Z0             !fresh water salinity  for next KPP_ABK80
    ZTMP1(JL,2) = RSICE          !salinity of brine for next KPP_ABK80

  ENDIF
ENDDO

! 1. Calculate density of fresh water and brine in surface layer
!    -----------------------------------------------------------

LLALPHA = .TRUE.
LLBETA  = .TRUE.

! density of fresh water at surface (PRHOH2O)
CALL KPP_ABK80 &
 & ( KIDIA    ,KFDIA    ,KLON     ,1_JPIM   ,LDKPPCAL ,&
 &   ZTMP1(KIDIA,1) ,PXO(KIDIA,1,1) ,ZTMP2(KIDIA,1) ,PTALPHA ,PSBETA ,& !PTALPHA,PSBETA dmy
 &   ZSIG0    ,ZSIG     ,PRHOH2O  ,LLALPHA  ,LLBETA   )

! density of brine at surface (ZRHOB)
CALL KPP_ABK80 &
 & ( KIDIA    ,KFDIA    ,KLON     ,1_JPIM   ,LDKPPCAL ,&
 &   ZTMP1(KIDIA,2) ,PXO(KIDIA,1,1) ,ZTMP2(KIDIA,1) ,PTALPHA ,PSBETA ,&!PTALPHA,PSBETA dmy
 &   ZSIG0    ,ZSIG     ,ZRHOB    ,LLALPHA  ,LLBETA   )

! 2. Calculate buoyancy profile (m/s**2) ,buoyancy gradients alpha and beta
!    ----------------------------------------------------------------------

CALL KPP_ABK80 & 
 & ( KIDIA    ,KFDIA    ,KLON     ,KLEVOP1  ,LDKPPCAL         ,&
 &   PXO(KIDIA,1,2), PXO(KIDIA,1,1), ZTMP2 ,PTALPHA(KIDIA,1) ,PSBETA(KIDIA,1) ,&
 &   ZSIG0    ,ZSIG     ,PRHO(1,1),LLALPHA          ,LLBETA   )

CALL KPP_CPSW &
 & ( KIDIA    ,KFDIA    ,KLON     ,KLEVOP1  ,LDKPPCAL,&
 &   PXO(KIDIA,1,2), PXO(KIDIA,1,1),ZTMP2 ,PCP(1,1))

DO JL = KIDIA, KFDIA
  IF( LDKPPCAL(JL) ) THEN

    DO JZ = 1, KLEVOP1
!      ZBUOY(JL,JZ)   = - RG * ZSIG0(JL) / 1000.0_JPRB
      ZBUOY(JL,JZ)   = - ZRG1000 * ZSIG0(JL,JZ)
    ENDDO

    PRHO(JL,0) = PRHO(JL,1)
    PCP(JL,0) = PCP(JL,1)
    PTALPHA(JL,0) = PTALPHA(JL,1)
    PSBETA(JL,0) = PSBETA(JL,1)  
    PRHO_INV(JL)= Z1 / PRHO(JL,0)
    ZCP_INV(JL)= Z1 / PCP(JL,0) 

! 3. Set fluxes
!    ----------

    IF ( ABS(PTAUX(JL)) < REPS_KPP .AND. ABS(PTAUY(JL)) < REPS_KPP ) THEN
      PTAUX(JL)=REPS_KPP
    ENDIF

!   positive downward 
    PSFLUX(JL,1) = PTAUX(JL)
    PSFLUX(JL,2) = PTAUY(JL)
    PSFLUX(JL,3) = PSWF(JL)
    PSFLUX(JL,4) = PLWF(JL) + PLHF(JL) + PSHF(JL) - PSNOW(JL) * RLMLT
    PSFLUX(JL,5) = Z0     ! melting of sea-ice = 0.0
!    PSFLUX(JL,6) = ( PRAIN(JL) + PSNOW(JL) + PLHF(JL)/RLVTT ) ! P-E
    PSFLUX(JL,6) = ( PRAIN(JL) + PSNOW(JL) + PLHF(JL)*ZRLVTT_INV ) ! P-E

    IF( LDSICE(JL) ) THEN
      ! Melting of sea-ice is proportional to the difference 
      ! btw brine and mdl temp. 
      PSFLUX(JL,3) = 0.0_JPRB 
!      PSFLUX(JL,4) = ( RTICE - PXO(JL,1,1) ) * PRHO(JL,0) &
!                   & * PCP(JL,0) * PDO(JL,0) / RRAMSICE
      PSFLUX(JL,4) = ( RTICE - PXO(JL,1,1) ) * PRHO(JL,0) &
                   & * PCP(JL,0) * PDO(JL,1) * ZRRAMSICE_INV
!      PSFLUX(JL,5) = - MIN( PSFLUX(JL,4), 0.0_JPRB ) / RLMLT
      PSFLUX(JL,5) = - MIN( PSFLUX(JL,4), 0.0_JPRB ) * ZRLMLT_INV
    ENDIF

    ZTMP3(JL) = -Z1 ! for next KPP_SWFRAC

  ENDIF
ENDDO

DO JZ=0,KLEVO
  IF( LDINIT_KPP ) THEN
    CALL KPP_SWFRAC &
     & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,JZ       ,&
     &   LDINIT_KPP, LDKPPCAL ,ZTMP3, PDO(1,JZ), PSWDK_SAVE ,&
     &   ZSWDK, YDOCEAN_ML )
  ELSE
    ZSWDK(KIDIA:KFDIA)=PSWDK_SAVE(KIDIA:KFDIA,JZ)
  ENDIF

  DO JL = KIDIA, KFDIA
    IF( LDKPPCAL(JL) ) THEN
!      PWXNT(JL,JZ,1) = - PSFLUX(JL,3) * ZSWDK(JL) &
!                     &   / ( PRHO(JL,0) * PCP(JL,0) )
      PWXNT(JL,JZ,1) = - PSFLUX(JL,3) * ZSWDK(JL) &
                     &   *  PRHO_INV(JL) * ZCP_INV(JL)
    ENDIF
  ENDDO

ENDDO

DO JL = KIDIA, KFDIA
  IF( LDKPPCAL(JL) ) THEN

!   calculate kinematic surface momentum fluxes
!    PWU(JL,0,1) = - PSFLUX(JL,1) / PRHO(JL,0)
    PWU(JL,0,1) = - PSFLUX(JL,1) * PRHO_INV(JL)
!    PWU(JL,0,2) = - PSFLUX(JL,2) / PRHO(JL,0)
    PWU(JL,0,2) = - PSFLUX(JL,2) * PRHO_INV(JL)
    ZTAU(JL)    = SQRT( PSFLUX(JL,1)*PSFLUX(JL,1) &
                      & + PSFLUX(JL,2)*PSFLUX(JL,2) )
!    ZUSTAR(JL)   = SQRT( ZTAU(JL) / PRHO(JL,0) )
    ZUSTAR(JL)   = SQRT( ZTAU(JL) * PRHO_INV(JL) )

!   total turbulent kinematic temperature flux (c m/s)
!    PWX(JL,0,1)  = - PSFLUX(JL,4) / PRHO(JL,0) / PCP(JL,0)
    PWX(JL,0,1)  = - PSFLUX(JL,4) * PRHO_INV(JL) * ZCP_INV(JL)

!   total turbulent kinematic salinity flux (o/oo m/s)
!    PWX(JL,0,2) = PXO(JL,1,2) * PSFLUX(JL,6) / PRHOH2O(JL) &
!                & + ( PXO(JL,1,2) - RSICE ) * PSFLUX(JL,5) / ZRHOB(JL)
    PWX(JL,0,2) = ( PXO(JL,1,2)*PSFLUX(JL,6)*ZRHOB(JL) &
                &  + (PXO(JL,1,2)-RSICE)*PSFLUX(JL,5)*PRHOH2O(JL) ) &
                & / ( PRHOH2O(JL)*ZRHOB(JL) )

!   calculate total kinematic surface buoyancy flux (m**2/s**3)
    ZB0(JL) = - RG * ( PTALPHA(JL,0) * PWX(JL,0,1) &
            &         - PSBETA(JL,0) * PWX(JL,0,2) )
    PWX(JL,0,NSCLRO+1) =  - ZB0(JL) 
!    ZB0SOL(JL) = RG * PTALPHA(JL,0) * PSFLUX(JL,3) &
!               & / ( PRHO(JL,0) * PCP(JL,0) )
    ZB0SOL(JL) = RG * PTALPHA(JL,0) * PSFLUX(JL,3) &
               & * PRHO_INV(JL) * ZCP_INV(JL)

!   calculate temperature and salt contributions of buoyancy gradients
!   on interfaces for double diffusion


    DO JZ = 1, KLEVO
      ZDXO(JL,JZ,1) =  PXO(JL,JZ,1) - PXO(JL,JZ+1,1)
      ZDXO(JL,JZ,2) =  PXO(JL,JZ,2) - PXO(JL,JZ+1,2)
      ZDUO(JL,JZ,1) =  PUO(JL,JZ,1) - PUO(JL,JZ+1,1)
      ZDUO(JL,JZ,2) =  PUO(JL,JZ,2) - PUO(JL,JZ+1,2)
      ZDZO(JL,JZ) =  PZO(JL,JZ) - PZO(JL,JZ+1)
!      ZALPHADT(JL,JZ) = 0.5_JPRB * ( PTALPHA(JL,JZ) + PTALPHA(JL,JZ+1) ) &
!                      & * ( PXO(JL,JZ,1) - PXO(JL,JZ+1,1) )
      ZALPHADT(JL,JZ) = 0.5_JPRB * ( PTALPHA(JL,JZ) + PTALPHA(JL,JZ+1) ) &
                      & * ZDXO(JL,JZ,1)
!      ZBETADS(JL,JZ) = 0.5_JPRB * ( PSBETA(JL,JZ) + PSBETA(JL,JZ+1) )  &
!                     & * ( PXO(JL,JZ,2) - PXO(JL,JZ+1,2) )
      ZBETADS(JL,JZ) = 0.5_JPRB * ( PSBETA(JL,JZ) + PSBETA(JL,JZ+1) )  &
                     & * ZDXO(JL,JZ,2)
    ENDDO

! 4. Compute buoyancy and shear profiles
!    -----------------------------------

    DO JZ = 1, KLEVO

      ZREF(JL) =  REPSILON_KPP * PZO(JL,JZ)
      ZREF_INV(JL) = 1.0_JPRB / ZREF(JL)

!     compute reference buoyancy and velocity
      ZWZ(JL)    = MAX( PZO(JL,1), ZREF(JL) )
!      ZUREF(JL)  = PUO(JL,1,1)  * ZWZ(JL) / ZREF(JL)
      ZUREF(JL)  = PUO(JL,1,1)  * ZWZ(JL) * ZREF_INV(JL)
!      ZVREF(JL)  = PUO(JL,1,2)  * ZWZ(JL) / ZREF(JL)
      ZVREF(JL)  = PUO(JL,1,2)  * ZWZ(JL) * ZREF_INV(JL)
!      ZBREF(JL)  = ZBUOY(JL,1) * ZWZ(JL) / ZREF(JL)
      ZBREF(JL)  = ZBUOY(JL,1) * ZWZ(JL) * ZREF_INV(JL)

      JKLOOP: DO JK = 1, KLEVO
        IF( ZREF(JL) >= PZO(JL,JK) ) EXIT JKLOOP
!        ZWZ(JL) = MIN( PZO(JL,JK)-PZO(JL,JK+1), PZO(JL,JK)-ZREF(JL) )
        ZWZ(JL) = MIN( ZDZO(JL,JK), PZO(JL,JK)-ZREF(JL) )
!        ZDEL(JL) = 0.5_JPRB * ZWZ(JL) / (PZO(JL,JK) - PZO(JL,JK+1))
        ZDEL(JL) = 0.5_JPRB * ZWZ(JL) / ZDZO(JL,JK)
!        ZUREF(JL) = ZUREF(JL) - ZWZ(JL) * ( PUO(JL,JK,1) + ZDEL(JL)  &
!                  & * ( PUO(JL,JK+1,1) - PUO(JL,JK,1))      ) / ZREF(JL)
        ZUREF(JL) = ZUREF(JL) - ZWZ(JL) * ( PUO(JL,JK,1) - ZDEL(JL)  &
                  & * ZDUO(JL,JK,1) ) * ZREF_INV(JL)
!        ZVREF(JL) = ZVREF(JL) - ZWZ(JL) * ( PUO(JL,JK,2) + ZDEL(JL)  &
!                  & * ( PUO(JL,JK+1,2) - PUO(JL,JK,2))      ) / ZREF(JL)
        ZVREF(JL) = ZVREF(JL) - ZWZ(JL) * ( PUO(JL,JK,2) - ZDEL(JL)  &
                  & * ZDUO(JL,JK,2) ) * ZREF_INV(JL)
!        ZBREF(JL) = ZBREF(JL) - ZWZ(JL) * ( ZBUOY(JL,JK) + ZDEL(JL) &
!                  & * ( ZBUOY(JL,JK+1) - ZBUOY(JL,JK))    ) / ZREF(JL)
        ZBREF(JL) = ZBREF(JL) - ZWZ(JL) * ( ZBUOY(JL,JK) + ZDEL(JL) &
                  & * ( ZBUOY(JL,JK+1) - ZBUOY(JL,JK))    ) * ZREF_INV(JL)
      ENDDO JKLOOP

      ZRITOP(JL,JZ) = ( ZREF(JL) - PZO(JL,JZ) )              &
                    & * ( ZBREF(JL) - ZBUOY(JL,JZ) )
      ZDBLOC(JL,JZ) = ZBUOY(JL,JZ) - ZBUOY(JL,JZ+1)
      ZDVSQ(JL,JZ)  = ( ZUREF(JL) - PUO(JL,JZ,1) )**2        &
                    & + ( ZVREF(JL) - PUO(JL,JZ,2) )**2
!      ZSHSQ(JL,JZ)  = ( PUO(JL,JZ,1) - PUO(JL,JZ+1,1) )**2   &
!                    & + ( PUO(JL,JZ,2) - PUO(JL,JZ+1,2) )**2
      ZSHSQ(JL,JZ)  = ZDUO(JL,JZ,1)*ZDUO(JL,JZ,1) &
                    & + ZDUO(JL,JZ,2)*ZDUO(JL,JZ,2)

    ENDDO

  ENDIF
ENDDO

! 5. Compute interior diffusivities 
!    ------------------------------

CALL KPP_INTERIOR_MIX &
 & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KLEVOP1  ,&
 &   LDKPPCAL ,ZSHSQ    ,ZDBLOC   ,PZO      ,ZDZO     ,&
 &   ZALPHADT ,ZBETADS  ,YDOCEAN_ML,&
 &   PDIFM    ,PDIFS    ,PDIFT    )


! 6. Compute non local KPP mixing
!    ----------------------------

IF(LKPP_KPP) THEN

! 6.1 Diagnose the new boundary layer depth
!     -------------------------------------

  CALL KPP_BLDEPTH &
   & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KLEVOP1  ,&
   &   LDINIT_KPP, LDKPPCAL ,PZO    ,ZDZO     ,PHO      ,&
   &   PSWDK_SAVE,ZDVSQ   ,ZDBLOC   ,ZRITOP   ,ZUSTAR   ,&
   &   ZB0      ,ZB0SOL   ,PF       ,PHBL     ,ZBFSFC   ,&
   &   ZSTABLE  ,ZCASEA   ,KBL      ,ZSIGMA   ,ZWM      ,&
   &   ZWS      ,POCDEPTH ,YDOCEAN_ML)

! 6.2 Compute boundary layer mixing
!     ----------------------------- 

  CALL KPP_BLMIX &
   & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KLEVOP1  ,&
   &   LDKPPCAL ,PZO      ,PHO      ,PHO_INV  ,ZUSTAR   ,&
   &   ZBFSFC   ,PHBL     ,ZSTABLE  ,ZCASEA   ,PDIFM    ,&
   &   PDIFS    ,PDIFT    ,KBL      ,PGHAT    ,ZSIGMA   ,&
   &   ZWM      ,ZWS      ,YDOCEAN_ML)

ENDIF ! LKPP_KPP


! 7. Set bottom layer diffusivity
!    ---------------------------- 

! Reset diffusivities for no bottom flux option

IF(LNBFLX_KPP) THEN 

  ZDLIMIT = Z0
  ZVLIMIT = Z0

  DO JL = KIDIA, KFDIA
    IF( LDKPPCAL(JL) ) THEN
      DO JS=1,NSCLRO
        PWXNT(JL,KLEVO,JS) = Z0
      ENDDO
    ENDIF
  ENDDO

ELSE ! Limit the bottom diffusity and viscosity

  ZDLIMIT = RDIFS_BOT
  ZVLIMIT = RDIFM_BOT

ENDIF

DO JL = KIDIA, KFDIA
  IF( LDKPPCAL(JL) ) THEN
    PDIFM(JL,KLEVO) = ZVLIMIT
    PDIFS(JL,KLEVO) = ZDLIMIT
    PDIFT(JL,KLEVO) = ZDLIMIT
    PGHAT(JL,KLEVO) = Z0
  ENDIF
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('KPP_KPPMIX_MOD:KPP_KPPMIX',1,ZHOOK_HANDLE)

END SUBROUTINE KPP_KPPMIX
END MODULE KPP_KPPMIX_MOD
