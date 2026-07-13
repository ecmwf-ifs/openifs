MODULE OCEAN_ML_DRIVER_MOD

CONTAINS

SUBROUTINE OCEAN_ML_DRIVER &
 & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KSTART   ,&
 &   KSTEP    ,LDOCN_KPP,LDSICE   ,PSSRFLTI ,PSLRFL   ,&
 &   PAHFSTI  ,PEVAPTI  ,PUSTRTI  ,PVSTRTI  ,PUSTRC   ,&
 &   PVSTRC   ,PRSFC    ,PRSFL    ,PSSFC    ,PSSFL    ,&
 &   PGEMU    ,PZO      ,PHO      ,PHO_INV  ,PDO      ,&
 &   POCDEPTH ,PUO0     ,PVO0     ,PUOC     ,PVOC     ,&
 &   PTO0     ,PSO0     ,PUOE1    ,PVOE1    ,PTOE1    ,&
 &   PSOE1    ,PADVT    ,PADVS    ,PTSTP    ,PDIFM    ,&
 &   PDIFT    ,PDIFS    ,PTRI0    ,PTRI1    ,PSWDK_SAVE,&
 &   YDCST    ,YDOCEAN_ML)


! (C) Copyright 1994- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! -------
!   This routine is the main driver for the ocean mixed layer model (KPP).

! Interface :
! ---------

! Method :
! ------
!   K-profile parameterization

! Externals :
! ---------

! Reference :
! ---------
!   Large, W. G., J. C. McWilliams, S. C. Doney (1994), Rev. Geophys.

! Modifications :
! -------------
!   06-Jun-1994  Bill Large
!          2002  Steve Woolnough, Reading Univ.
!   07-Oct-2008  Yuhei Takaya,    E.C.M.W.F.    Implemented to IFS.
! End Modifications :
!---------------------------------------------------------------------

USE PARKIND1,     ONLY : JPIM, JPRB
USE YOMHOOK,      ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_CST,      ONLY : TCST
USE YOS_OCEAN_ML, ONLY : TOCEAN_ML, NVELO, NSCLRO
USE KPP_SWFRAC_MOD
USE KPP_KPPMIX_MOD
USE KPP_OCNINT_MOD

IMPLICIT NONE

INTEGER(KIND = JPIM), INTENT(IN) :: KIDIA
INTEGER(KIND = JPIM), INTENT(IN) :: KFDIA
INTEGER(KIND = JPIM), INTENT(IN) :: KLON
INTEGER(KIND = JPIM), INTENT(IN) :: KLEVO
INTEGER(KIND = JPIM), INTENT(IN) :: KSTART
INTEGER(KIND = JPIM), INTENT(IN) :: KSTEP
INTEGER(KIND = JPIM) :: IBLS(KLON,0:1)     

REAL(KIND = JPRB), INTENT(INOUT)  :: PZO(KLON,KLEVO+1)
REAL(KIND = JPRB), INTENT(INOUT)  :: PDO(KLON,0:KLEVO)
REAL(KIND = JPRB), INTENT(INOUT)  :: PHO(KLON,KLEVO+1)
REAL(KIND = JPRB), INTENT(INOUT)  :: PHO_INV(KLON,KLEVO+1)

REAL(KIND = JPRB), INTENT(IN)    :: PSSRFLTI(KLON) ! Net SW flx at sfc [W m^-2]
REAL(KIND = JPRB), INTENT(IN)    :: PSLRFL(KLON)   ! Net LW flx at sfc [W m^-2]
REAL(KIND = JPRB), INTENT(IN)    :: PAHFSTI(KLON)  ! SH flux downwards [W m^-2]
REAL(KIND = JPRB), INTENT(IN)    :: PEVAPTI(KLON)  ! surf. moisture flux
REAL(KIND = JPRB), INTENT(IN)    :: PUSTRTI(KLON)  ! X-momentum flux [N m^-2]
REAL(KIND = JPRB), INTENT(IN)    :: PVSTRTI(KLON)  ! Y-momentum flux [N m^-2]
REAL(KIND = JPRB), INTENT(IN)    :: PUSTRC(KLON)   ! taux climatology 
REAL(KIND = JPRB), INTENT(IN)    :: PVSTRC(KLON)   ! tauy climatology
REAL(KIND = JPRB), INTENT(IN)    :: POCDEPTH(KLON) ! Depth of ML model [m]
REAL(KIND = JPRB), INTENT(IN)    :: PRSFC(KLON)    ! convective rain
REAL(KIND = JPRB), INTENT(IN)    :: PRSFL(KLON)    ! large scale rain
REAL(KIND = JPRB), INTENT(IN)    :: PSSFC(KLON)    ! convective snow
REAL(KIND = JPRB), INTENT(IN)    :: PSSFL(KLON)    ! large scale snow
REAL(KIND = JPRB), INTENT(IN)    :: PGEMU(KLON)    ! The SIN of latitude
REAL(KIND = JPRB), INTENT(IN)    :: PTSTP          ! The model time step [s]
REAL(KIND = JPRB), INTENT(INOUT) :: PADVT(KLON,KLEVO+1)  ! advection 
REAL(KIND = JPRB), INTENT(INOUT) :: PADVS(KLON,KLEVO+1)  ! advection 
REAL(KIND = JPRB), INTENT(OUT)   :: PDIFM(KLON,0:KLEVO)  ! visc. coef.
REAL(KIND = JPRB), INTENT(OUT)   :: PDIFS(KLON,0:KLEVO)  ! diff. of sal.
REAL(KIND = JPRB), INTENT(OUT)   :: PDIFT(KLON,0:KLEVO)  ! diff. of temp.
REAL(KIND = JPRB), INTENT(INOUT) :: PUO0(KLON,KLEVO+1)   ! velocity
REAL(KIND = JPRB), INTENT(INOUT) :: PVO0(KLON,KLEVO+1)   ! velocity
REAL(KIND = JPRB), INTENT(INOUT) :: PUOC(KLON,KLEVO+1)   ! velocity climatology
REAL(KIND = JPRB), INTENT(INOUT) :: PVOC(KLON,KLEVO+1)   ! velocity climatology
REAL(KIND = JPRB), INTENT(INOUT) :: PTO0(KLON,KLEVO+1)   ! scalars
REAL(KIND = JPRB), INTENT(INOUT) :: PSO0(KLON,KLEVO+1)   ! scalars
REAL(KIND = JPRB), INTENT(INOUT) :: PUOE1(KLON,KLEVO+1)  ! tendency of velocity
REAL(KIND = JPRB), INTENT(INOUT) :: PVOE1(KLON,KLEVO+1)  ! tendency of velocity
REAL(KIND = JPRB), INTENT(INOUT) :: PTOE1(KLON,KLEVO+1)  ! tendency of scalars
REAL(KIND = JPRB), INTENT(INOUT) :: PSOE1(KLON,KLEVO+1)  ! tendency of scalars
REAL(KIND = JPRB), INTENT(INOUT) :: PTRI0(KLON,0:KLEVO)  !array for diff. eq.
REAL(KIND = JPRB), INTENT(INOUT) :: PTRI1(KLON,0:KLEVO)  !array for diff. eq.
REAL(KIND = JPRB), INTENT(INOUT) :: PSWDK_SAVE(KLON,0:KLEVO) !coef of radiation

REAL(KIND = JPRB) :: ZUO0(KIDIA:KFDIA,KLEVO+1,NVELO) ! velocity
REAL(KIND = JPRB) :: ZXO0(KIDIA:KFDIA,KLEVO+1,NSCLRO)! scalars

REAL(KIND = JPRB) :: ZUOS(KIDIA:KFDIA,KLEVO+1,NVELO,0:1)  
REAL(KIND = JPRB) :: ZXOS(KIDIA:KFDIA,KLEVO+1,NSCLRO,0:1) 
REAL(KIND = JPRB) :: ZADV(KIDIA:KFDIA,KLEVO,NSCLRO)    ! advection 
REAL(KIND = JPRB) :: ZHBLS(KLON,0:1)               
REAL(KIND = JPRB) :: ZGHAT(KLON,KLEVO)
REAL(KIND = JPRB) :: ZRHO(KLON,0:KLEVO+1)
REAL(KIND = JPRB) :: ZCP(KLON,0:KLEVO+1)
REAL(KIND = JPRB) :: ZTALPHA(KIDIA:KFDIA,0:KLEVO+1)
REAL(KIND = JPRB) :: ZSBETA(KIDIA:KFDIA,0:KLEVO+1)
REAL(KIND = JPRB) :: ZWX(KLON,0:KLEVO,NSCLRO+1)
REAL(KIND = JPRB) :: ZWXNT(KLON,0:KLEVO,NSCLRO)
REAL(KIND = JPRB) :: ZWU(KLON,0:KLEVO,NVELO)

LOGICAL, INTENT(IN) ::  LDOCN_KPP(KLON)  ! =.TRUE.  grid points coupled to 
                                         !          the ocean mixed layer
                                         ! =.FALSE. not computed
LOGICAL, INTENT(IN) ::  LDSICE(KLON)     ! =.TRUE.  sea-ice grid 

TYPE(TCST),      INTENT(IN)    :: YDCST
TYPE(TOCEAN_ML), INTENT(INOUT) :: YDOCEAN_ML

REAL(KIND = JPRB) :: ZF(KIDIA:KFDIA)     ! Coriolis parameter
REAL(KIND = JPRB) :: ZTAUX(KLON)
REAL(KIND = JPRB) :: ZTAUY(KLON)
REAL(KIND = JPRB) :: ZSWF(KLON)
REAL(KIND = JPRB) :: ZLWF(KLON)
REAL(KIND = JPRB) :: ZLHF(KLON)
REAL(KIND = JPRB) :: ZSHF(KLON)
REAL(KIND = JPRB) :: ZRAIN(KLON)         ! Precipitation [kg m^-2]
REAL(KIND = JPRB) :: ZSNOW(KLON)         ! Snow [kg m^-2]
REAL(KIND = JPRB) :: ZSFLUX(KLON,6)      ! surface flux
REAL(KIND = JPRB) :: ZDTO                ! Time step of KPP
REAL(KIND = JPRB) :: ZDIV
REAL(KIND = JPRB) :: ZFRAC    
REAL(KIND = JPRB) :: ZHBLEPS(KIDIA:KFDIA)              
REAL(KIND = JPRB) :: ZRHOH2O(KLON)
REAL(KIND = JPRB) :: ZRHO_INV(KIDIA:KFDIA)              
REAL(KIND = JPRB) :: ZSINMIN
REAL(KIND = JPRB) :: ZDELTAZ(KIDIA:KFDIA)

REAL(KIND = JPRB) :: ZDZB(KLEVO)
REAL(KIND = JPRB) :: ZSUMH(KIDIA:KFDIA)
REAL(KIND = JPRB) :: ZDFAC(KIDIA:KFDIA)
REAL(KIND = JPRB) :: ZSK(KIDIA:KFDIA)
REAL(KIND = JPRB) :: ZHSUM(KIDIA:KFDIA)
REAL(KIND = JPRB) :: ZSWDK(KLON)
REAL(KIND = JPRB) :: ZTMP1(KLON)
REAL(KIND = JPRB) :: ZSUM(KIDIA:KFDIA)
REAL(KIND = JPRB) :: ZA
REAL(KIND = JPRB) :: ZB(KIDIA:KFDIA)
REAL(KIND = JPRB) :: ZX1
REAL(KIND = JPRB) :: ZX2

REAL(KIND = JPRB) :: ZROMEGA2

INTEGER(KIND = JPIM) :: ILEVOP1 
INTEGER(KIND = JPIM) :: IKBL(KLON)       ! KBL
INTEGER :: INEW   ! = 1
INTEGER :: IOLD   ! = 0
INTEGER :: ICONV(KLON)                   ! counter for convergence test

INTEGER :: JT
INTEGER :: JITER
INTEGER :: JZ
INTEGER :: JL
INTEGER :: JV
INTEGER :: JS

LOGICAL :: LLKPPCAL(KLON)                 
LOGICAL :: LLCONV
LOGICAL :: LLINIT_KPP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!--------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OCEAN_ML_DRIVER_MOD:OCEAN_ML_DRIVER',0,ZHOOK_HANDLE)
ASSOCIATE(RG=>YDCST%RG, RLVTT=>YDCST%RLVTT, ROMEGA=>YDCST%ROMEGA, &
 & RPI=>YDCST%RPI, &
 & LCDIAG_KPP=>YDOCEAN_ML%LCDIAG_KPP, LDD_KPP=>YDOCEAN_ML%LDD_KPP, &
 & LGEO_KPP=>YDOCEAN_ML%LGEO_KPP, LVL_NEW_KPP=>YDOCEAN_ML%LVL_NEW_KPP, &
 & NITERMAX_KPP=>YDOCEAN_ML%NITERMAX_KPP, NOCNSTEP_KPP=>YDOCEAN_ML%NOCNSTEP_KPP, &
 & RDSCALE_KPP=>YDOCEAN_ML%RDSCALE_KPP, RLAMBDA_KPP=>YDOCEAN_ML%RLAMBDA_KPP, &
 & RTOLFAC_KPP=>YDOCEAN_ML%RTOLFAC_KPP, RVC_KPP=>YDOCEAN_ML%RVC_KPP, &
 & RVDEL_KPP=>YDOCEAN_ML%RVDEL_KPP, RVD_KPP=>YDOCEAN_ML%RVD_KPP, &
 & RVE_KPP=>YDOCEAN_ML%RVE_KPP)

! 1. Initial Setting
!    ---------------

IOLD = 0  
INEW = 1  
ILEVOP1 = KLEVO + 1

!IF ( NSTEP == NSTART ) THEN
IF ( KSTEP == KSTART ) THEN
  LLINIT_KPP = .TRUE.
ELSE
  LLINIT_KPP = .FALSE.
ENDIF

ZDTO = PTSTP / REAL ( NOCNSTEP_KPP ) 

ZFRAC = 1.0_JPRB - RLAMBDA_KPP

IF(LGEO_KPP) THEN
  DO JZ = 1, KLEVO+1
    DO JL = KIDIA, KFDIA
      ZUO0(JL,JZ,1) = PUO0(JL,JZ) + PUOC(JL,JZ)
      ZUO0(JL,JZ,2) = PVO0(JL,JZ) + PVOC(JL,JZ)
      ZXO0(JL,JZ,1) = PTO0(JL,JZ)
      ZXO0(JL,JZ,2) = PSO0(JL,JZ)
    ENDDO
  ENDDO
ELSE
  DO JZ = 1, KLEVO+1
    DO JL = KIDIA, KFDIA
      ZUO0(JL,JZ,1) = PUO0(JL,JZ)
      ZUO0(JL,JZ,2) = PVO0(JL,JZ)
      ZXO0(JL,JZ,1) = PTO0(JL,JZ)
      ZXO0(JL,JZ,2) = PSO0(JL,JZ)
    ENDDO
  ENDDO
ENDIF

DO JZ = 1, KLEVO
  DO JL = KIDIA, KFDIA
    ZADV(JL,JZ,1) = PADVT(JL,JZ)
    ZADV(JL,JZ,2) = PADVS(JL,JZ)
  ENDDO
ENDDO

! 1.1. Accummulate Flux
!      ----------------

ZSINMIN = SIN( 0.5_JPRB * RPI / 180.0_JPRB )
ZROMEGA2 = ROMEGA * 2.0_JPRB

DO JL = KIDIA, KFDIA
  IF( LDOCN_KPP(JL) ) THEN
    ZTAUX(JL) = PUSTRTI(JL)
    ZTAUY(JL) = PVSTRTI(JL) 
    ZSWF(JL)  = PSSRFLTI(JL)
    ZLWF(JL)  = PSLRFL(JL)        !ZLWF: downward positive 
    ZLHF(JL)  = PEVAPTI(JL)*RLVTT !ZLHF: latent heat, downward positive
!    ZLHF(JL)  = PEVAPTI(JL)  !ZLHF: latent heat, downward positive
    ZSHF(JL)  = PAHFSTI(JL) ! ZSHF: downward positive
    ZRAIN(JL) = PRSFC(JL) + PRSFL(JL)
    ZSNOW(JL) = PSSFC(JL) + PSSFL(JL)

    ZF(JL)= ZROMEGA2 * PGEMU(JL)

  ENDIF
ENDDO

IF( LLINIT_KPP ) THEN

! 1.2. Set coordinates
!      ---------------

PHO(KIDIA:KFDIA,:)=0.0_JPRB
PZO(KIDIA:KFDIA,:)=0.0_JPRB
PDO(KIDIA:KFDIA,:)=0.0_JPRB

IF(LVL_NEW_KPP) THEN

  ZA = RVDEL_KPP * REAL(KLEVO)
  DO JL = KIDIA, KFDIA
    ZB(JL) = POCDEPTH(JL) - ZA - RVD_KPP
  ENDDO
  ZSUM(KIDIA:KFDIA) = 0.0_JPRB

  DO JZ=1,KLEVO
    ZX1 = REAL(JZ)/REAL(KLEVO)
    ZX2 = REAL(JZ-1)/REAL(KLEVO-1)
    DO JL = KIDIA, KFDIA
      PDO(JL,JZ) = ZA*ZX1 + ZB(JL)*(ZX2**RVC_KPP) &
                 & + RVD_KPP*(ZX2**RVE_KPP)
      PHO(JL,JZ) = PDO(JL,JZ) - ZSUM(JL)
      PZO(JL,JZ) = -( ZSUM(JL) + 0.5_JPRB*PHO(JL,JZ) )
      PHO_INV(JL,JZ) = 1.0_JPRB / PHO(JL,JZ)
      ZSUM(JL)   = ZSUM(JL) + PHO(JL,JZ)
    ENDDO
  ENDDO

ELSE ! Stretched grid of Large et al. (1994)

  IF ( RDSCALE_KPP > 0.0_JPRB ) THEN
    DO JL = KIDIA, KFDIA
      ZSUMH(JL) = 0.0_JPRB
      ZDFAC(JL) = 1.0_JPRB - EXP( -RDSCALE_KPP )
    ENDDO

    DO JZ = 1, KLEVO
      DO JL = KIDIA, KFDIA
        ZSK(JL) = - ( REAL( JZ ) - 0.5_JPRB ) / REAL( KLEVO )
        PHO(JL,JZ) = POCDEPTH(JL) * ZDFAC(JL) / REAL( KLEVO ) &
                   & / RDSCALE_KPP / ( 1.0_JPRB + ZSK(JL) * ZDFAC(JL) )
        ZSUMH(JL) = ZSUMH(JL) + PHO(JL,JZ)
      ENDDO
    ENDDO
  ENDIF

  ! To ensure the PZO(KLEVO) = - POCDEPTH

  ZHSUM(KIDIA:KFDIA) = 0.0_JPRB
  PDO(KIDIA:KFDIA,0) = 0.0_JPRB

  DO JZ = 1, KLEVO
    DO JL = KIDIA, KFDIA
      IF( RDSCALE_KPP > 0.0_JPRB ) THEN
        PHO(JL,JZ) = PHO(JL,JZ) * POCDEPTH(JL) / ZSUMH(JL)
      ELSE
        PHO(JL,JZ) = POCDEPTH(JL) / REAL(KLEVO)
      ENDIF
      PHO_INV(JL,JZ) = 1.0_JPRB / PHO(JL,JZ)
      PZO(JL,JZ) =  - ( ZHSUM(JL) + 0.5_JPRB * PHO(JL,JZ) )
      ZHSUM(JL) = ZHSUM(JL) + PHO(JL,JZ)
      PDO(JL,JZ) = ZHSUM(JL)
!     WRITE(NULOUT,'(A4,I3,A4,I3,3(3X,A4,F9.3))') &
!&         ' JL=',JL,' JZ=',JZ,' PHO=',PHO(JL,JZ),&
!&         ' PDO=',PDO(JL,JZ),' PZO=',PZO(JL,JZ)
    ENDDO
  ENDDO

ENDIF

DO JL = KIDIA, KFDIA
  PHO(JL,KLEVO+1) = 1.E-10_JPRB ! small number
  PHO_INV(JL,KLEVO+1) = 1.0_JPRB / PHO(JL,KLEVO+1)
  PZO(JL,KLEVO+1) = - POCDEPTH(JL)
ENDDO

! 1.3. Radiation table
!      ---------------

ZTMP1(KIDIA:KFDIA) = -1.0_JPRB
DO JZ = 0, KLEVO
  CALL KPP_SWFRAC &
   & ( KIDIA    ,KFDIA     ,KLON       ,KLEVO    ,JZ       ,&
   &   LLINIT_KPP,LDOCN_KPP,ZTMP1      ,PDO(1,JZ),PSWDK_SAVE,&
   &   ZSWDK    ,YDOCEAN_ML )  ! PDO >= 0
ENDDO

! 1.4. Set trigonal matrix for implicit scheme
!      ---------------------------------------

  DO JL = KIDIA, KFDIA
    IF( LDOCN_KPP(JL) ) THEN
      DO JZ = 1, KLEVO
        ZDZB(JZ) = PZO(JL,JZ) - PZO(JL,JZ+1)
      ENDDO
!      PTRI1(JL,0) = ZDTO / PHO(JL,1)
!      PTRI1(JL,1) = ZDTO / PHO(JL,1) / ZDZB(1) 
      PTRI1(JL,0) = ZDTO * PHO_INV(JL,1)
      PTRI1(JL,1) = ZDTO * PHO_INV(JL,1) / ZDZB(1) 
      DO JZ = 2, KLEVO
!        PTRI1(JL,JZ) = ZDTO / PHO(JL,JZ) / ZDZB(JZ)
        PTRI1(JL,JZ) = ZDTO * PHO_INV(JL,JZ) / ZDZB(JZ)
                      ! dt / h(k) / {dzb(k-1)=dzabove}
!        PTRI0(JL,JZ) = ZDTO / PHO(JL,JZ) / ZDZB(JZ-1) 
        PTRI0(JL,JZ) = ZDTO * PHO_INV(JL,JZ) / ZDZB(JZ-1) 
                      ! dt / h(k) / {dzb(k  )=dzbelow}
      ENDDO
    ENDIF
  ENDDO

  LLINIT_KPP=.FALSE.

ENDIF ! LLINIT_KPP

! 2. Time Integration
!    ----------------

DO JT = 1, NOCNSTEP_KPP 
! NOCNSTEP_KPP should be 1 for IFS, because the time integration 
! is done in POSTPHY.  If ZUO,ZXO are updated (integrated) in JT 
! loop, we can activate NOCNSTEP_KPP.

! 2.1. Estimate new profiles
!      ---------------------

! Estimate new profiles by extraporation for Euler-Backward scheme

    DO JZ = 1, ILEVOP1 
      DO JL = KIDIA, KFDIA
        IF( LDOCN_KPP(JL) ) THEN
          ZUOS(JL,JZ,1,INEW) = ZUO0(JL,JZ,1) 
          ZUOS(JL,JZ,1,IOLD) = ZUOS(JL,JZ,1,INEW)
          ZUOS(JL,JZ,2,INEW) = ZUO0(JL,JZ,2) 
          ZUOS(JL,JZ,2,IOLD) = ZUOS(JL,JZ,2,INEW)
          ZXOS(JL,JZ,1,INEW) = ZXO0(JL,JZ,1) !+ PTOE1(JL,JZ)*ZDTO
          ZXOS(JL,JZ,1,IOLD) = ZXOS(JL,JZ,1,INEW) 
          ZXOS(JL,JZ,2,INEW) = ZXO0(JL,JZ,2) !+ PSOE1(JL,JZ)*ZDTO
          ZXOS(JL,JZ,2,IOLD) = ZXOS(JL,JZ,2,INEW) 
        ENDIF
      ENDDO
    ENDDO

    DO JL = KIDIA, KFDIA
      IF( LDOCN_KPP(JL) ) THEN !Bug fix
        ZHBLS(JL,INEW) = 100000.0_JPRB  !initial setting
        IBLS(JL,INEW) = 1               !not used 
      ENDIF
    ENDDO

! 2.2 Iteration of Euler-Backward Scheme
!     ----------------------------------

  DO JL = KIDIA, KFDIA
    ICONV(JL) = 0
    IF( LDOCN_KPP(JL) ) THEN
      LLKPPCAL(JL) = .TRUE.  ! Compute KPP scheme, solution is not convergent yet
    ELSE
      LLKPPCAL(JL) = .FALSE. 
    ENDIF
  ENDDO

  DO JITER = 1, NITERMAX_KPP

    DO JZ = 1, ILEVOP1
      DO JL = KIDIA, KFDIA
        IF( LLKPPCAL(JL) ) THEN
          ZUOS(JL,JZ,1,INEW) = RLAMBDA_KPP * ZUOS(JL,JZ,1,IOLD) &
                              & + ZFRAC * ZUOS(JL,JZ,1,INEW) 
          ZUOS(JL,JZ,1,IOLD) = ZUOS(JL,JZ,1,INEW)
          ZUOS(JL,JZ,2,INEW) = RLAMBDA_KPP * ZUOS(JL,JZ,2,IOLD) &
                              & + ZFRAC * ZUOS(JL,JZ,2,INEW) 
          ZUOS(JL,JZ,2,IOLD) = ZUOS(JL,JZ,2,INEW)
          ZXOS(JL,JZ,1,INEW) = RLAMBDA_KPP * ZXOS(JL,JZ,1,IOLD) &
                                & + ZFRAC * ZXOS(JL,JZ,1,INEW)
          ZXOS(JL,JZ,1,IOLD) = ZXOS(JL,JZ,1,INEW)
          ZXOS(JL,JZ,2,INEW) = RLAMBDA_KPP * ZXOS(JL,JZ,2,IOLD) &
                                & + ZFRAC * ZXOS(JL,JZ,2,INEW)
          ZXOS(JL,JZ,2,IOLD) = ZXOS(JL,JZ,2,INEW)
        ENDIF
      ENDDO
    ENDDO
    DO JL = KIDIA, KFDIA
      IF( LLKPPCAL(JL) ) THEN
        ZHBLS(JL,IOLD) = ZHBLS(JL,INEW)
        IBLS(JL,IOLD)  = IBLS(JL,INEW)
      ENDIF
    ENDDO

! 2.2. Compute vertical diffusion coefficients   
!      ---------------------------------------

    CALL KPP_KPPMIX &
     & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,ILEVOP1  ,&
     &   LLINIT_KPP, LLKPPCAL ,LDSICE   ,ZTAUX  ,ZTAUY    ,&
     &   ZSWF     ,ZLWF     ,ZLHF     ,ZSHF     ,ZRAIN    ,&
     &   ZSNOW    ,PZO      ,PHO      ,PHO_INV  ,PDO      ,&
     &   ZSFLUX   ,ZUOS(KIDIA,1,1,INEW),ZXOS(KIDIA,1,1,INEW),ZHBLS(1,INEW),IBLS(1,INEW),&
     &   ZF       ,PDIFM    ,PDIFS    ,PDIFT    ,ZGHAT    ,&
     &   ZRHO     ,ZCP      ,ZWU      ,ZWX      ,ZWXNT    ,&
     &   PSWDK_SAVE ,ZTALPHA,ZSBETA   ,POCDEPTH ,ZRHOH2O  ,&
     &   ZRHO_INV ,YDCST    ,YDOCEAN_ML )

    IF(LGEO_KPP) THEN 
      DO JL = KIDIA, KFDIA
        IF( LLKPPCAL(JL) ) THEN
          DO JZ = 1, ILEVOP1
            ZUOS(JL,JZ,1,INEW) = ZUOS(JL,JZ,1,INEW) - PUOC(JL,JZ)
            ZUOS(JL,JZ,2,INEW) = ZUOS(JL,JZ,2,INEW) - PVOC(JL,JZ)
            ZUO0(JL,JZ,1) = ZUO0(JL,JZ,1) - PUOC(JL,JZ)
            ZUO0(JL,JZ,2) = ZUO0(JL,JZ,2) - PVOC(JL,JZ)
          ENDDO
          ZWU(JL,0,1) = ZWU(JL,0,1) + PUSTRC(JL)*ZRHO_INV(JL) !caution for sign, ZWU = -TAUX
          ZWU(JL,0,2) = ZWU(JL,0,2) + PVSTRC(JL)*ZRHO_INV(JL) !                  ZWU = -TAUY
        ENDIF
      ENDDO
    ENDIF

    CALL KPP_OCNINT &
     & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,ILEVOP1  ,&
     &   LLKPPCAL ,PZO      ,PHO      ,PHO_INV  ,PDO      ,&
     &   PTRI0    ,PTRI1    ,IBLS(1,INEW), ZUO0 ,ZXO0     ,&
     &   ZF       ,ZCP      ,ZRHO     ,PDIFM    ,PDIFT    ,&
     &   PDIFS    ,ZGHAT    ,ZADV     ,ZUOS(KIDIA,1,1,INEW),ZXOS(KIDIA,1,1,INEW),&
     &   ZWX      ,ZWXNT    ,ZWU      ,ZDTO     ,YDOCEAN_ML )


    IF(LGEO_KPP) THEN
      DO JL = KIDIA, KFDIA
        IF( LLKPPCAL(JL) ) THEN
          DO JZ = 1, ILEVOP1
            ZUOS(JL,JZ,1,INEW) = ZUOS(JL,JZ,1,INEW) + PUOC(JL,JZ)
            ZUOS(JL,JZ,2,INEW) = ZUOS(JL,JZ,2,INEW) + PVOC(JL,JZ)
            ZUO0(JL,JZ,1) = ZUO0(JL,JZ,1) + PUOC(JL,JZ)
            ZUO0(JL,JZ,2) = ZUO0(JL,JZ,2) + PVOC(JL,JZ)
          ENDDO
          ZWU(JL,0,1) = ZWU(JL,0,1) - PUSTRC(JL)*ZRHO_INV(JL)
          ZWU(JL,0,2) = ZWU(JL,0,2) - PVSTRC(JL)*ZRHO_INV(JL)
        ENDIF
      ENDDO
    ENDIF


!  2.3. Check convergence
!       -----------------

    LLCONV = .TRUE.

    DO JL = KIDIA, KFDIA
      IF( LLKPPCAL(JL) ) THEN

        ZHBLEPS(JL) = RTOLFAC_KPP * PHO(JL,IBLS(JL,INEW))
        IF( IBLS(JL,INEW) == ILEVOP1 ) THEN
          ZHBLEPS(JL) = RTOLFAC_KPP * PHO(JL,KLEVO)
        ENDIF
        IF( ABS( ZHBLS(JL,INEW)-ZHBLS(JL,IOLD) ) > ZHBLEPS(JL) ) THEN
          ICONV(JL) = 0
        ELSE
          ICONV(JL) = ICONV(JL) + 1
        ENDIF

! The criteria for convergence is based on the boundary layer depth.
! Sometimes, this doesn't work after convective vertical mixing, but 
! the solution seems to be converged, so the original source codes
! was not changed.  Y.T.

        IF( ICONV(JL) > 0 ) THEN
          LLKPPCAL(JL) = .FALSE.
        ELSEIF( ( JITER > NITERMAX_KPP - 2 )  .AND.  &
              & ( ZHBLS(JL,INEW) < ZHBLS(JL,IOLD) )  ) THEN
!       select smaller HBL for safety.
          LLKPPCAL(JL) = .FALSE.
        ENDIF

        IF( LLKPPCAL(JL) ) LLCONV = .FALSE.

      ENDIF
    ENDDO

    IF( LLCONV ) EXIT ! Solutions at all grid points are converged.

  ENDDO ! JITER

! 3. Tendency
!    --------

  ZDIV = 1.0_JPRB / ZDTO

  DO JZ = 1, ILEVOP1
    DO JL = KIDIA, KFDIA
      IF( LDOCN_KPP(JL) ) THEN
        PUOE1(JL,JZ) = ( ZUOS(JL,JZ,1,INEW) - ZUO0(JL,JZ,1) ) * ZDIV
        PVOE1(JL,JZ) = ( ZUOS(JL,JZ,2,INEW) - ZUO0(JL,JZ,2) ) * ZDIV
        PTOE1(JL,JZ) = ( ZXOS(JL,JZ,1,INEW) - ZXO0(JL,JZ,1) ) * ZDIV
        PSOE1(JL,JZ) = ( ZXOS(JL,JZ,2,INEW) - ZXO0(JL,JZ,2) ) * ZDIV
      ELSE
        PUOE1(JL,JZ) = 0.0_JPRB
        PVOE1(JL,JZ) = 0.0_JPRB
        PTOE1(JL,JZ) = 0.0_JPRB
        PSOE1(JL,JZ) = 0.0_JPRB
      ENDIF
    ENDDO
  ENDDO

  DO JZ = 1, KLEVO
    DO JL = KIDIA, KFDIA
      IF( LDOCN_KPP(JL) ) THEN
        PADVT(JL,JZ) = ZADV(JL,JZ,1)  
        PADVS(JL,JZ) = ZADV(JL,JZ,2) 
      ENDIF
    ENDDO
  ENDDO

! 4. Compute diagnostic fluxes
!    -------------------------

! ZWX(0,*),ZWU(0,*) is computed in KPP_KPPMIX. 
! Caution, fluxes are computed with t0 step for just only diagnostics.

  IF( LCDIAG_KPP ) THEN
    DO JZ = 1, KLEVO
      DO JL = KIDIA, KFDIA
        IF( LDOCN_KPP(JL) ) THEN
          ZDELTAZ(JL) = 0.5_JPRB * ( PHO(JL,JZ) + PHO(JL,JZ+1) )
          DO JS=1,NSCLRO 
            ZWX(JL,JZ,JS) = -PDIFS(JL,JZ) &
                        & * ( ( ZXO0(JL,JZ,JS) - ZXO0(JL,JZ+1,JS) ) / ZDELTAZ(JL) &
                        &      - ZGHAT(JL,JZ) * ZWX(JL,0,JS) )
          ENDDO
          IF(LDD_KPP) THEN
            ZWX(JL,JZ,1) = -PDIFT(JL,JZ) &
                      & * ( ( ZXO0(JL,JZ,1) - ZXO0(JL,JZ+1,1) ) / ZDELTAZ(JL) &
                      &     - ZGHAT(JL,JZ) * ZWX(JL,0,1) )
          ENDIF 
          DO JV=1,NVELO
            ZWU(JL,JZ,JV)= -PDIFM(JL,JZ) &
                         & * ( ZUO0(JL,JZ,JV) - ZUO0(JL,JZ+1,JV) ) / ZDELTAZ(JL)
          ENDDO
          ZWX(JL,JZ,NSCLRO+1) = RG * ( ZTALPHA(JL,JZ) * ZWX(JL,JZ,1) &
                            &        - ZSBETA(JL,JZ) * ZWX(JL,JZ,2) )
        ENDIF
      ENDDO
    ENDDO
  ENDIF

!  WRITE(NULOUT,*) 'OCEAN MIXED LAYER MODEL, NO OF ITERATION =',JITER

ENDDO ! JT loop

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('OCEAN_ML_DRIVER_MOD:OCEAN_ML_DRIVER',1,ZHOOK_HANDLE)

END SUBROUTINE OCEAN_ML_DRIVER
END MODULE OCEAN_ML_DRIVER_MOD
