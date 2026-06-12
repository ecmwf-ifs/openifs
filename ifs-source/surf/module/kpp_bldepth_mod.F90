MODULE KPP_BLDEPTH_MOD
CONTAINS 
SUBROUTINE KPP_BLDEPTH &
  & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KLEVOP1  ,&
  &   LDINIT_KPP, LDKPPCAL ,PZO    ,PDZO     ,PHO      ,&
  &   PSWDK_SAVE, PDVSQ  ,PDBLOC   ,PRITOP   ,PUSTAR   ,&
  &   PB0      ,PB0SOL   ,PF       ,PHBL     ,PBFSFC   ,&
  &   PSTABLE  ,PCASEA   ,KBL      ,PSIGMA   ,PWM      ,&
  &   PWS      ,POCDEPTH ,YDOCEAN_ML )

! (C) Copyright 1994- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! Purpose :
! -------
!   This routine calculates boundary layer depth.

! Interface :
! ---------
!   Call *KPP_BLDEPTH* from *KPP_KPPMIX*

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

USE KPP_SWFRAC_MOD
USE KPP_WSCALE_MOD

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)  :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)  :: KLON
INTEGER(KIND=JPIM),INTENT(IN)  :: KLEVO
INTEGER(KIND=JPIM),INTENT(IN)  :: KLEVOP1
INTEGER(KIND=JPIM),INTENT(OUT) :: KBL(KLON)      ! Layer Idx of the first 
                                                 ! level below PHBL 
REAL(KIND=JPRB),INTENT(IN) :: PZO(KLON,KLEVOP1)  ! vertical grid (<= 0)     (m)
REAL(KIND=JPRB),INTENT(IN) :: PDZO(KIDIA:KFDIA,KLEVO) ! PZO(JZ)-PZO(JZ+1)
REAL(KIND=JPRB),INTENT(IN) :: PHO(KLON,KLEVOP1)  ! layer thicknesses        (m)
REAL(KIND=JPRB),INTENT(INOUT) :: PSWDK_SAVE(KLON,0:KLEVO) ! coef. for radiation
REAL(KIND=JPRB),INTENT(IN) :: PDVSQ(KIDIA:KFDIA,KLEVO)
                                                 ! (vel shear re sfc)^2 (m/s)^2
REAL(KIND=JPRB),INTENT(IN) :: PDBLOC(KIDIA:KFDIA,KLEVO) 
                                                 ! local delta buoyancy (m/s^2)
REAL(KIND=JPRB),INTENT(IN) :: PRITOP(KIDIA:KFDIA,KLEVO) 
                                                 ! numerator of bulk Ri (m/s)^2
REAL(KIND=JPRB),INTENT(IN) :: PUSTAR(KLON)       ! sfc friction velocity  (m/s)
REAL(KIND=JPRB),INTENT(IN) :: PB0(KIDIA:KFDIA)   ! sfc turb. buoy. forcing 
                                                 !                    (m^2/s^3)
REAL(KIND=JPRB),INTENT(IN) :: PB0SOL(KIDIA:KFDIA)! rad. buoy. forcing (m^2/s^3)
REAL(KIND=JPRB),INTENT(IN) :: PF(KIDIA:KFDIA)    ! Coriolis parameter     (1/s)
REAL(KIND=JPRB),INTENT(IN) :: POCDEPTH(KLON)     ! Ocean depth (>=0)        (m)
REAL(KIND=JPRB),INTENT(OUT) :: PHBL(KLON)        ! boundary layer depth     (m)
REAL(KIND=JPRB),INTENT(OUT) :: PBFSFC(KLON)      ! bo+radiation absorbed 
                                                 ! to d=hbf*hbl       (m^2/s^3)
REAL(KIND=JPRB),INTENT(OUT) :: PSTABLE(KIDIA:KFDIA) ! =1 in stable forcing
                                                 ! =0 unstable
REAL(KIND=JPRB),INTENT(OUT) :: PCASEA(KLON)      ! =1 case A: 
                                                 !    PHBL is near to PZO(KBL)
                                                 ! =0 case B: 
                                                 !    PHBL is near to PZO(KBL-1)
REAL(KIND=JPRB),INTENT(OUT) :: PSIGMA(KLON)      ! normalized depth (PD/PHBL)
REAL(KIND=JPRB),INTENT(OUT) :: PWM(KIDIA:KFDIA)  ! turbulent vel. scales  (m/s)
REAL(KIND=JPRB),INTENT(OUT) :: PWS(KIDIA:KFDIA)  ! turbulent scalar scales(m/s)

LOGICAL,INTENT(IN) :: LDKPPCAL(KLON)
LOGICAL,INTENT(IN) :: LDINIT_KPP

TYPE(TOCEAN_ML),INTENT(INOUT) :: YDOCEAN_ML

INTEGER(KIND=JPIM) :: ITMP
INTEGER(KIND=JPIM) :: IA
INTEGER(KIND=JPIM) :: IU
INTEGER(KIND=JPIM) :: JL
INTEGER(KIND=JPIM) :: JZ

REAL(KIND=JPRB) :: ZRIB(KIDIA:KFDIA,2)   ! Bulk Richardson number
REAL(KIND=JPRB) :: ZDMO(KIDIA:KFDIA,2)   ! Monin-Obukhov depth          (m)
REAL(KIND=JPRB) :: ZDEK(KIDIA:KFDIA)     ! Ekman depth                  (m)
REAL(KIND=JPRB) :: ZHMONOB(KIDIA:KFDIA)  ! Monin-Obukhov depth          (m)
REAL(KIND=JPRB) :: ZHRI(KIDIA:KFDIA)     ! Depth of critical Ri number  (m)
REAL(KIND=JPRB) :: ZHMIN(KIDIA:KFDIA)                  
REAL(KIND=JPRB) :: ZHMIN2(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZVTSQ(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZBVSQ(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZDEKMAN(KIDIA:KFDIA)
REAL(KIND=JPRB) :: Z1        ! 1.0
REAL(KIND=JPRB) :: ZTMP1(KLON)
REAL(KIND=JPRB) :: ZRVONK_INV
REAL(KIND=JPRB) :: Z1STABLE(KIDIA:KFDIA)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('KPP_BLDEPTH_MOD:KPP_BLDEPTH',0,ZHOOK_HANDLE)
ASSOCIATE(RCEKMAN_KPP=>YDOCEAN_ML%RCEKMAN_KPP, &
 & RCMONOB_KPP=>YDOCEAN_ML%RCMONOB_KPP, RCS_KPP=>YDOCEAN_ML%RCS_KPP, &
 & RCV_KPP=>YDOCEAN_ML%RCV_KPP, REPSILON_KPP=>YDOCEAN_ML%REPSILON_KPP, &
 & REPS_KPP=>YDOCEAN_ML%REPS_KPP, RICR_KPP=>YDOCEAN_ML%RICR_KPP, &
 & RVONK=>YDOCEAN_ML%RVONK, RVTC=>YDOCEAN_ML%RVTC)

Z1 = 1.0_JPRB
ZRVONK_INV = 1.0_JPRB / RVONK

! eq. (23) in LMD94
!ZVTC =  RCV_KPP * SQRT( 0.2_JPRB / RCS_KPP / REPSILON_KPP ) &
!     &  / RVONK**2 / RICR_KPP

! Indices for array rib(i,k), the bulk richardson number.
IA = 1
IU = 2

! Initialize hbl and kbl to bottomed out values

DO JL = KIDIA, KFDIA
  IF ( LDKPPCAL(JL) ) THEN
    ZRIB(JL,IA) = 0.0_JPRB
    ZDMO(JL,IA) = -PZO(JL,KLEVOP1)
    KBL(JL)     = KLEVO
    PHBL(JL)    = -PZO(JL,KLEVO)
    ZDEK(JL)    = RCEKMAN_KPP * PUSTAR(JL) / ( ABS(PF(JL)) + REPS_KPP )
    ZTMP1(JL)   = Z1 !used next KPP_SWFAC
  ENDIF
ENDDO

DO JZ = 2, KLEVO

  IF( LDINIT_KPP ) THEN
    CALL KPP_SWFRAC &
     & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO      ,JZ       ,&
     &   LDINIT_KPP, LDKPPCAL ,ZTMP1  ,PZO(1,JZ) ,PSWDK_SAVE,&
     &   PBFSFC   ,YDOCEAN_ML) 
  ELSE
    PBFSFC(KIDIA:KFDIA)=PSWDK_SAVE(KIDIA:KFDIA,JZ)
  ENDIF

  DO JL = KIDIA, KFDIA
    IF ( LDKPPCAL(JL) ) THEN

      IF( KBL(JL) >= KLEVO ) THEN
!       use PCASEA as temporary array for next call to WSCALE
        PCASEA(JL)  = -PZO(JL,JZ)
 
!       compute bfsfc= bo + radiative contribution down to hbf * hbl
        PBFSFC(JL)  = PB0(JL) + PB0SOL(JL) * ( Z1 - PBFSFC(JL) )
        PSTABLE(JL) = 0.5_JPRB + SIGN( 0.5_JPRB, PBFSFC(JL) + REPS_KPP )
        Z1STABLE(JL) =  Z1-PSTABLE(JL)
        PSIGMA(JL)  = PSTABLE(JL) +  Z1STABLE(JL) * REPSILON_KPP

      ENDIF

    ENDIF
  ENDDO !JL loop

! compute velocity scales at sigma, for hbl= casea = -zgrid(kl)
  CALL KPP_WSCALE &
   & ( KIDIA    ,KFDIA    ,KLON     ,LDKPPCAL ,PSIGMA   ,&
   &   PCASEA   ,PUSTAR   ,PBFSFC   ,PWM      ,PWS      ,&
   &   YDOCEAN_ML )

  DO JL = KIDIA, KFDIA
    IF ( LDKPPCAL(JL) ) THEN

      IF( KBL(JL) >= KLEVO ) THEN
 
!       compute the turbulent shear contribution to rib
!        ZBVSQ(JL) = 0.5_JPRB &
!                  & * ( PDBLOC(JL,JZ-1) / ( PZO(JL,JZ-1) - PZO(JL,JZ)   ) &
!                  &     + PDBLOC(JL,JZ) / ( PZO(JL,JZ)   - PZO(JL,JZ+1) ) )
!        ZBVSQ(JL) = 0.5_JPRB &
!                  & * ( PDBLOC(JL,JZ-1) / PDZO(JL,JZ-1) &
!                  &     + PDBLOC(JL,JZ) / PDZO(JL,JZ)   )
        ZBVSQ(JL) = 0.5_JPRB &
                  & * ( PDBLOC(JL,JZ-1)*PDZO(JL,JZ) + PDBLOC(JL,JZ)*PDZO(JL,JZ-1) ) &
                  &   / ( PDZO(JL,JZ)*PDZO(JL,JZ-1) )   
        ZVTSQ(JL) = - PZO(JL,JZ) * PWS(JL) * SQRT(ABS(ZBVSQ(JL))) * RVTC

!       compute bulk richardson number at new level, dunder
        ZRIB(JL,IU) = PRITOP(JL,JZ) / ( PDVSQ(JL,JZ) + ZVTSQ(JL) + REPS_KPP )
        ZRIB(JL,IU) = MAX( ZRIB(JL,IU), ZRIB(JL,IA) + REPS_KPP )

!       linear interpolate to find hbl where rib = ricr
!        ZHRI(JL) = -PZO(JL,JZ-1) + ( PZO(JL,JZ-1) - PZO(JL,JZ) )     &
!                 &                * ( RICR_KPP - ZRIB(JL,IA) )    &
!                 &                / ( ZRIB(JL,IU) - ZRIB(JL,IA) )
        ZHRI(JL) = -PZO(JL,JZ-1) +  PDZO(JL,JZ-1) &
                 &                * ( RICR_KPP - ZRIB(JL,IA) )    &
                 &                / ( ZRIB(JL,IU) - ZRIB(JL,IA) )

!       compute the monin obukov length scale 
!        ZDMO(JL,IU) = RCMONOB_KPP * PUSTAR(JL) * PUSTAR(JL) * PUSTAR(JL) &
!                    & / RVONK / ( ABS(PBFSFC(JL)) + REPS_KPP )
        ZDMO(JL,IU) = RCMONOB_KPP * PUSTAR(JL) * PUSTAR(JL) * PUSTAR(JL) &
                    & * ZRVONK_INV / ( ABS(PBFSFC(JL)) + REPS_KPP )
!        ZDMO(JL,IU) = PSTABLE(JL) * ZDMO(JL,IU) &
!                    & - (Z1-PSTABLE(JL))*PZO(JL,KLEVOP1) !Changed from LMD
        ZDMO(JL,IU) = PSTABLE(JL) * ZDMO(JL,IU) &
                    & - Z1STABLE(JL) * PZO(JL,KLEVOP1) !Changed from LMD
        IF( -PZO(JL,JZ) >= ZDMO(JL,IU) ) THEN
!          ZHMONOB(JL) = ( ZDMO(JL,IU) - ZDMO(JL,IA) ) &
!                      & / ( PZO(JL,JZ-1) - PZO(JL,JZ) )
          ZHMONOB(JL) = ( ZDMO(JL,IU) - ZDMO(JL,IA) ) / PDZO(JL,JZ-1)
          ZHMONOB(JL) = ( ZDMO(JL,IU) + ZHMONOB(JL)*PZO(JL,JZ) ) &
                      & / ( Z1 - ZHMONOB(JL) )
        ELSE
          ZHMONOB(JL) = -PZO(JL,KLEVOP1)
        ENDIF

!       Added by Y. Takaya 
!        ZHMONOB(JL) = PSTABLE(JL) * ZHMONOB(JL) &
!                    & - ( Z1 - PSTABLE(JL) ) * PZO(JL,KLEVOP1)
        ZHMONOB(JL) = PSTABLE(JL) * ZHMONOB(JL) &
                    & - Z1STABLE(JL) * PZO(JL,KLEVOP1)
 
!       compute the ekman depth 
!        ZDEKMAN(JL) = PSTABLE(JL) * ZDEK(JL) &
!                    & - ( Z1 - PSTABLE(JL) ) * PZO(JL,KLEVOP1)
        ZDEKMAN(JL) = PSTABLE(JL) * ZDEK(JL) &
                    & - Z1STABLE(JL) * PZO(JL,KLEVOP1)
      
!       compute boundary layer depth
        ZHMIN(JL) = MIN( ZHRI(JL), ZHMONOB(JL), ZDEKMAN(JL), POCDEPTH(JL) )

        IF( ZHMIN(JL)  <  -PZO(JL,JZ) ) THEN

!         code below added by sjw 09/07/04 to solve problems where hek 
!         less than zgrid(kl-1) giving negative diffusions
!         if this occurs too often then we need to rethink this fix

!          IF ( .NOT. LDINIT_KPP ) THEN
            IF ( ZHMIN(JL)  <  -PZO(JL,JZ-1) ) THEN
              ZHMIN2(JL) = MIN( ZHRI(JL) , POCDEPTH(JL) )
              IF ( ZHMIN2(JL)  <  -PZO(JL,JZ) ) THEN
               ZHMIN(JL) = ZHMIN2(JL)
              ENDIF
            ENDIF
!          ENDIF

          PHBL(JL) = ZHMIN(JL)
          KBL(JL) =  JZ

        ENDIF

      ENDIF 

    ENDIF
  ENDDO !JL loop

  ITMP  = IA
  IA    = IU
  IU    = ITMP

ENDDO !JZ loop

DO JL = KIDIA, KFDIA
  ZTMP1(JL)=-Z1
ENDDO

CALL KPP_SWFRAC &
 & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,-1_JPIM  ,&
 &   LDINIT_KPP,LDKPPCAL ,ZTMP1    ,PHBL     ,PSWDK_SAVE ,&
 &   PBFSFC    ,YDOCEAN_ML )

DO JL = KIDIA, KFDIA
  IF ( LDKPPCAL(JL) ) THEN
    PBFSFC(JL)  = PB0(JL) + PB0SOL(JL) * ( Z1 - PBFSFC(JL) )
    PSTABLE(JL) = 0.5_JPRB + SIGN( 0.5_JPRB, PBFSFC(JL) ) !1:stable,0:unstable 
    PBFSFC(JL)  = PBFSFC(JL) + PSTABLE(JL) * REPS_KPP  ! ensures bfsfc never=0 
    PCASEA(JL)  = 0.5_JPRB  &
                & - SIGN( 0.5_JPRB , &
                &         PZO(JL,KBL(JL))+0.5_JPRB*PHO(JL,KBL(JL))+PHBL(JL) )
  ENDIF
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('KPP_BLDEPTH_MOD:KPP_BLDEPTH',1,ZHOOK_HANDLE)

END SUBROUTINE KPP_BLDEPTH
END MODULE KPP_BLDEPTH_MOD
