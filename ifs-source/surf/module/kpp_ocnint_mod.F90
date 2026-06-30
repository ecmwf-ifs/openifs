MODULE KPP_OCNINT_MOD
CONTAINS
SUBROUTINE KPP_OCNINT &
  & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KLEVOP1  ,&
  &   LDKPPCAL ,PZO      ,PHO      ,PHO_INV  ,PDO      ,&
  &   PTRI0    ,PTRI1    ,KBL      ,PUO      ,PXO      ,&
  &   PF       ,PCP      ,PRHO     ,PDIFM    ,PDIFT    ,&
  &   PDIFS    ,PGHAT    ,PADV     ,PUN      ,PXN      ,&
  &   PWX      ,PWXNT    ,PWU      ,PDTO     ,YDOCEAN_ML )

! (C) Copyright 1994- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! Purpose :
! -------
!   This routine calculates time integration of the ocean model.

! Interface :
! ---------

! Method :
! ------
!   Integrate the ocean model by Euler-Backward (implicit) scheme

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
USE YOS_OCEAN_ML, ONLY : TOCEAN_ML, NVELO, NSCLRO

USE KPP_TRIDCOF_MOD
USE KPP_TRIDMAT_MOD
USE KPP_TRIDRHS_MOD

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVO
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVOP1
INTEGER(KIND=JPIM),INTENT(IN) :: KBL(KLON)              ! index of gridpoint 
                                                        ! just below h
REAL(KIND=JPRB),INTENT(IN) :: PUO(KIDIA:KFDIA,KLEVOP1,NVELO)   ! old velocity profiles
REAL(KIND=JPRB),INTENT(IN) :: PXO(KIDIA:KFDIA,KLEVOP1,NSCLRO)  ! old scalar profiles
REAL(KIND=JPRB),INTENT(IN) :: PZO(KLON,KLEVOP1)         ! vertical layer grids
REAL(KIND=JPRB),INTENT(IN) :: PHO(KLON,KLEVOP1)         ! layer thickness
REAL(KIND=JPRB),INTENT(IN) :: PHO_INV(KLON,KLEVOP1)     ! 1/PHO
REAL(KIND=JPRB),INTENT(IN) :: PDO(KLON,0:KLEVO)         ! layer interface depth
REAL(KIND=JPRB),INTENT(IN) :: PCP(KLON,0:KLEVOP1)       ! heat capacitiy
REAL(KIND=JPRB),INTENT(IN) :: PRHO(KLON,0:KLEVOP1)      ! density
REAL(KIND=JPRB),INTENT(IN) :: PF(KIDIA:KFDIA)           ! Coriolis parameter
REAL(KIND=JPRB),INTENT(IN) :: PDIFM(KLON,0:KLEVO)       ! viscosity
REAL(KIND=JPRB),INTENT(IN) :: PDIFT(KLON,0:KLEVO)       ! diffusivity of heat
REAL(KIND=JPRB),INTENT(IN) :: PDIFS(KLON,0:KLEVO)       ! diffusivity of saliity
REAL(KIND=JPRB),INTENT(IN) :: PTRI0(KLON,0:KLEVO)       ! array for diff. eq.
REAL(KIND=JPRB),INTENT(IN) :: PTRI1(KLON,0:KLEVO)       ! array for diff. eq.
REAL(KIND=JPRB),INTENT(IN) :: PADV(KIDIA:KFDIA,KLEVO,NSCLRO)   ! advection correction
REAL(KIND=JPRB),INTENT(IN) :: PDTO                      ! time step for KPP
REAL(KIND=JPRB),INTENT(IN) :: PWX(KLON,0:KLEVO,NSCLRO+1)! scalar fluxes
REAL(KIND=JPRB),INTENT(IN) :: PWXNT(KLON,0:KLEVO,NSCLRO)! non-turbulent fluxes
REAL(KIND=JPRB),INTENT(IN) :: PWU(KLON,0:KLEVO,NVELO)   ! momentum fluxes
REAL(KIND=JPRB),INTENT(IN) :: PGHAT(KLON,KLEVO)         ! ghat
REAL(KIND=JPRB),INTENT(OUT) :: PUN(KIDIA:KFDIA,KLEVOP1,NVELO)  ! new velocity profiles
REAL(KIND=JPRB),INTENT(OUT) :: PXN(KIDIA:KFDIA,KLEVOP1,NSCLRO) ! new scalar profiles

LOGICAL,INTENT(IN) :: LDKPPCAL(KLON)            

TYPE(TOCEAN_ML), INTENT(IN) :: YDOCEAN_ML

INTEGER(KIND=JPIM) :: JZ                   ! loop control variable
INTEGER(KIND=JPIM) :: JL                   ! loop control variable
INTEGER(KIND=JPIM) :: JS                   ! loop control variable

REAL(KIND=JPRB) :: ZCU (KIDIA:KFDIA,KLEVO) ! upper coeff for (k-1) 
                                           ! on k line of trid. mtx
REAL(KIND=JPRB) :: ZCC (KIDIA:KFDIA,KLEVO) ! central ...     (k  ) ..
REAL(KIND=JPRB) :: ZCL (KIDIA:KFDIA,KLEVO) ! lower .....     (k-1) ..
REAL(KIND=JPRB) :: ZRHS(KIDIA:KFDIA,KLEVO) ! right-hand-side terms
REAL(KIND=JPRB) :: ZNTFLX(KIDIA:KFDIA,0:KLEVO,NSCLRO) ! non-turbulent flux
REAL(KIND=JPRB) :: ZGHATFLUX(KIDIA:KFDIA)  ! surface flux incl. solar
REAL(KIND=JPRB) :: ZSTURFLUX(KIDIA:KFDIA)  ! surface flux
REAL(KIND=JPRB) :: ZFHF(KIDIA:KFDIA)       ! (Corioris Parameter)/2  
REAL(KIND=JPRB) :: ZDTOFHF(KIDIA:KFDIA)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('KPP_OCNINT_MOD:KPP_OCNINT',0,ZHOOK_HANDLE)

! 1. Calculate u and v solution of tridiagonal matrix
!    ------------------------------------------------

DO JL = KIDIA, KFDIA
  IF ( LDKPPCAL(JL) ) THEN
    ZFHF(JL) = 0.5_JPRB * PF(JL)
    ZDTOFHF(JL) = PDTO * ZFHF(JL)
  ENDIF
ENDDO

! 1.1 Set coefficients of tridiagonal matrix
!     --------------------------------------

CALL KPP_TRIDCOF &
 & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,LDKPPCAL ,&
 &   PDIFM    ,PTRI0    ,PTRI1    ,ZCU      ,ZCC      ,&
 &   ZCL      )


! 1.2 Set U right hand side and solve matrix
!     --------------------------------------

DO JL = KIDIA, KFDIA
  IF ( LDKPPCAL(JL) ) THEN
!   top layer
!    ZRHS(JL,1) = PUO(JL,1,1) &
!               & + PDTO * ( ZFHF(JL) * ( PUO(JL,1,2) + PUN(JL,1,2) ) &
!               &            - PWU(JL,0,1) / PHO(JL,1) ) 
    ZRHS(JL,1) = PUO(JL,1,1) &
               & + PDTO * ( ZFHF(JL) * ( PUO(JL,1,2) + PUN(JL,1,2) ) &
               &            - PWU(JL,0,1) * PHO_INV(JL,1) ) 
!   bottom
!    ZRHS(JL,KLEVO) = PUO(JL,KLEVO,1) &
!                & + PDTO * ZFHF(JL) * ( PUO(JL,KLEVO,2) + PUN(JL,KLEVO,2) ) &
!                & + PTRI1(JL,KLEVO) * PDIFM(JL,KLEVO) * PUO(JL,KLEVO+1,1)
    ZRHS(JL,KLEVO) = PUO(JL,KLEVO,1) &
                & + ZDTOFHF(JL) * ( PUO(JL,KLEVO,2) + PUN(JL,KLEVO,2) ) &
                & + PTRI1(JL,KLEVO) * PDIFM(JL,KLEVO) * PUO(JL,KLEVO+1,1)
    DO JZ = 2, KLEVO-1
      ZRHS(JL,JZ) = PUO(JL,JZ,1) &
                  & + ZDTOFHF(JL) * ( PUO(JL,JZ,2) + PUN(JL,JZ,2) ) 
    ENDDO
  ENDIF
ENDDO

CALL KPP_TRIDMAT &
 & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KLEVOP1  ,&
 &   LDKPPCAL ,ZCU      ,ZCC      ,ZCL      ,ZRHS     ,&
 &   PUO(KIDIA,1,1)     ,YDOCEAN_ML,&
 &   PUN(KIDIA,1,1) )


! 1.3 Set V right hand side and solve matrix
!     --------------------------------------

DO JL = KIDIA, KFDIA
  IF ( LDKPPCAL(JL) ) THEN
!    ZRHS(JL,1) = PUO(JL,1,2) &
!               & - PDTO * ( ZFHF(JL)*(PUO(JL,1,1)+PUN(JL,1,1)) &
!               &            + PWU(JL,0,2) / PHO(JL,1) )
    ZRHS(JL,1) = PUO(JL,1,2) &
               & - PDTO * ( ZFHF(JL)*(PUO(JL,1,1)+PUN(JL,1,1)) &
               &            + PWU(JL,0,2) * PHO_INV(JL,1) )
!    ZRHS(JL,KLEVO) = PUO(JL,KLEVO,2) &
!                   & - PDTO*ZFHF(JL)*(PUO(JL,KLEVO,1)+PUN(JL,KLEVO,1)) &
!                   & + PTRI1(JL,KLEVO)*PDIFM(JL,KLEVO)*PUO(JL,KLEVO+1,2)
    ZRHS(JL,KLEVO) = PUO(JL,KLEVO,2) &
                   & - ZDTOFHF(JL)*(PUO(JL,KLEVO,1)+PUN(JL,KLEVO,1)) &
                   & + PTRI1(JL,KLEVO)*PDIFM(JL,KLEVO)*PUO(JL,KLEVO+1,2)
    DO JZ=2,KLEVO-1
      ZRHS(JL,JZ) = PUO(JL,JZ,2) &
                  & - ZDTOFHF(JL) * ( PUO(JL,JZ,1) + PUN(JL,JZ,1) ) 
    ENDDO
  ENDIF
ENDDO

CALL KPP_TRIDMAT &
 & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KLEVOP1  ,&
 &   LDKPPCAL ,ZCU      ,ZCC      ,ZCL      ,ZRHS     ,&
 &   PUO(KIDIA,1,2) ,    YDOCEAN_ML,&
 &   PUN(KIDIA,1,2) )


! 2. Calculate T and S solution of tridiagonal matrix
!    ------------------------------------------------

! 2.1 Set T right hand side and solve matrix
!     --------------------------------------

! scalar solutions of tridiagonal matrix
!     temperature (different from other scalars because of ghat-term
!                  and double diffusion)

DO JL = KIDIA, KFDIA
  IF ( LDKPPCAL(JL) ) THEN
    ZGHATFLUX(JL) = PWX(JL,0,1)
    ZSTURFLUX(JL) = PWX(JL,0,1)
    ZNTFLX(JL,0,1) = PWXNT(JL,0,1)
    DO JZ = 1, KLEVO
      ZNTFLX(JL,JZ,1) = PWXNT(JL,JZ,1)
    ENDDO
  ENDIF
ENDDO

CALL KPP_TRIDCOF &
 & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,LDKPPCAL ,&
 &   PDIFT    ,PTRI0    ,PTRI1    ,ZCU      ,ZCC      ,&
 &   ZCL      )

CALL KPP_TRIDRHS &
 & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KLEVOP1  ,&
 &   LDKPPCAL ,PHO      ,PHO_INV  ,PTRI0    ,PTRI1    ,&
 &   PXO(KIDIA,1,1), ZNTFLX(KIDIA,0,1), PDIFT ,PGHAT,ZSTURFLUX,&
 &   ZGHATFLUX,ZRHS     ,PDTO     ,PADV     ,PRHO     ,&
 &   PCP      ,1_JPIM   ,YDOCEAN_ML)

CALL KPP_TRIDMAT &
 & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KLEVOP1  ,&
 &   LDKPPCAL ,ZCU      ,ZCC      ,ZCL      ,ZRHS     ,&
 &   PXO(KIDIA,1,1)     ,YDOCEAN_ML,&
 &   PXN(KIDIA,1,1) )

! 2.2 Set S right hand side and solve matrix
!     --------------------------------------

CALL KPP_TRIDCOF &
 & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,LDKPPCAL ,&
 &   PDIFS    ,PTRI0    ,PTRI1    ,ZCU      ,ZCC      ,&
 &   ZCL      )

DO JS = 2, NSCLRO

   DO JL = KIDIA, KFDIA
     IF ( LDKPPCAL(JL) ) THEN
       DO JZ = 0, KLEVO !KLEVOP1
         ZNTFLX(JL,JZ,JS) = PWXNT(JL,JZ,JS)
       ENDDO
       ZGHATFLUX(JL) = PWX(JL,0,JS) 
       ZSTURFLUX(JL) = PWX(JL,0,JS)
     ENDIF
   ENDDO

   CALL KPP_TRIDRHS &
    & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KLEVOP1  ,&
    &   LDKPPCAL ,PHO      ,PHO_INV  ,PTRI0    ,PTRI1    ,&
    &   PXO(KIDIA,1,JS),ZNTFLX(KIDIA,0,JS),PDIFS ,PGHAT,ZSTURFLUX,&
    &   ZGHATFLUX,ZRHS     ,PDTO     ,PADV     ,PRHO     ,&
    &   PCP      ,JS       ,YDOCEAN_ML)

   CALL KPP_TRIDMAT &
    & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KLEVOP1  ,&
    &   LDKPPCAL ,ZCU      ,ZCC      ,ZCL      ,ZRHS     ,&
    &   PXO(KIDIA,1,JS)    ,YDOCEAN_ML,&
    &   PXN(KIDIA,1,JS) )

ENDDO

IF (LHOOK) CALL DR_HOOK('KPP_OCNINT_MOD:KPP_OCNINT',1,ZHOOK_HANDLE)
 
END SUBROUTINE KPP_OCNINT
END MODULE KPP_OCNINT_MOD
