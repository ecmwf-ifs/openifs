MODULE KPP_INTERIOR_MIX_MOD
CONTAINS
SUBROUTINE KPP_INTERIOR_MIX &
  & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KLEVOP1  , &
  &   LDKPPCAL ,PSHSQ    ,PDBLOC   ,PZO      ,PDZO     , &
  &   PALPHADT ,PBETADS  ,YDOCEAN_ML, &
  &   PDIFM    ,PDIFS    ,PDIFT    )

! (C) Copyright 1994- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! Purpose :
! -------
!   This routine calculates interior viscosity/diffusivity 
!   coefficients dependent on local Richardson Number and
!   due to background internal wave mixing and double 
!   diffusion.

! Interface :
! ---------
!   Call *KPP_INTERIOR_MIX* from *KPP_KPPMIX*

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

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVO
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVOP1

REAL(KIND=JPRB),INTENT(IN)  :: PSHSQ(KIDIA:KFDIA,KLEVO) 
                                                  !(local vel. shear)^2 (m/s)^2
REAL(KIND=JPRB),INTENT(IN)  :: PDBLOC(KIDIA:KFDIA,KLEVO) 
                                                  ! local delta buoy    (m/s^2)
REAL(KIND=JPRB),INTENT(IN)  :: PZO(KLON,KLEVOP1)  ! vertical grid (<= 0) (m)
REAL(KIND=JPRB),INTENT(IN)  :: PDZO(KIDIA:KFDIA,KLEVO)      ! PZO(JZ)-PZO(JZ+1)
REAL(KIND=JPRB),INTENT(IN)  :: PALPHADT(KIDIA:KFDIA,KLEVO)  ! alpha*dt  
REAL(KIND=JPRB),INTENT(IN)  :: PBETADS(KIDIA:KFDIA,KLEVO)   ! beta*ds
REAL(KIND=JPRB),INTENT(OUT) :: PDIFM(KLON,0:KLEVO) ! visc. coef. (m^2/s)
REAL(KIND=JPRB),INTENT(OUT) :: PDIFS(KLON,0:KLEVO) ! scalar diff.(m^2/s)
REAL(KIND=JPRB),INTENT(OUT) :: PDIFT(KLON,0:KLEVO) ! temp. diff. (m^2/s)

LOGICAL,INTENT(IN) :: LDKPPCAL(KLON)

TYPE(TOCEAN_ML),INTENT(IN) :: YDOCEAN_ML

INTEGER :: JL                           ! loop control
INTEGER :: JZ                           ! loop control
INTEGER :: JM                           ! loop control

REAL(KIND=JPRB) :: ZRIG(KIDIA:KFDIA)    ! local richardson number
REAL(KIND=JPRB) :: ZFRI(KIDIA:KFDIA)    ! function of Ri_g
REAL(KIND=JPRB) :: ZRATIO(KIDIA:KFDIA)  ! Ri_g / Ri_0
REAL(KIND=JPRB) :: ZTMP(KIDIA:KFDIA)    ! temporary variable for smoothing 
REAL(KIND=JPRB) :: ZWT(KIDIA:KFDIA)     ! weight for smoothing
REAL(KIND=JPRB) :: ZRRHO(KIDIA:KFDIA)   ! double diffusion parameter
REAL(KIND=JPRB) :: ZDIFFDD(KIDIA:KFDIA) ! double diffusion diffusivity scale
REAL(KIND=JPRB) :: ZPRANDTL(KIDIA:KFDIA)! prandtl number
REAL(KIND=JPRB) :: Z0                   ! 0.0
REAL(KIND=JPRB) :: Z1                   ! 1.0
REAL(KIND=JPRB) :: Z2                   ! 2.0
REAL(KIND=JPRB) :: ZRRIINFTY_INV
REAL(KIND=JPRB) :: ZRRRHO01
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('KPP_INTERIOR_MIX_MOD:KPP_INTERIOR_MIX',0,ZHOOK_HANDLE)
ASSOCIATE(LDD_KPP=>YDOCEAN_ML%LDD_KPP, LRI_KPP=>YDOCEAN_ML%LRI_KPP, &
 & LSF_NEW_KPP=>YDOCEAN_ML%LSF_NEW_KPP, NMRI_KPP=>YDOCEAN_ML%NMRI_KPP, &
 & RDIFM_IW=>YDOCEAN_ML%RDIFM_IW, RDIFM_MAX=>YDOCEAN_ML%RDIFM_MAX, &
 & RDIFS_IW=>YDOCEAN_ML%RDIFS_IW, RDIFS_MAX=>YDOCEAN_ML%RDIFS_MAX, &
 & RDMOL=>YDOCEAN_ML%RDMOL, RDSFMAX=>YDOCEAN_ML%RDSFMAX, &
 & REPS_KPP=>YDOCEAN_ML%REPS_KPP, RRIINFTY=>YDOCEAN_ML%RRIINFTY, &
 & RRRHO0=>YDOCEAN_ML%RRRHO0, RRRHO0_D06=>YDOCEAN_ML%RRRHO0_D06)

Z0 = 0.0_JPRB
Z1 = 1.0_JPRB
Z2 = 2.0_JPRB

! 1. Compute interior gradient Ri at all interfaces, except surface
!    --------------------------------------------------------------

IF(LRI_KPP) THEN

  ! eq. (27) in LMD94
  DO JL = KIDIA, KFDIA
    IF( LDKPPCAL(JL) ) THEN

      DO JZ= 1, KLEVO
        PDIFM(JL,JZ) = PDBLOC(JL,JZ) * PDZO(JL,JZ) &
                     & / ( PSHSQ(JL,JZ) + REPS_KPP )
      ENDDO

! 2. Vertically smooth MRI times
!    ---------------------------

! PDIFS is used for work valiable here, top(0) value is used as a dummy

      DO JM = 1, NMRI_KPP

        PDIFS(JL,0)       = Z0
        PDIFM(JL,0)       = Z0

        DO JZ=1,KLEVO
          IF( ( PDIFM(JL,JZ) < Z0 ) .OR. ( PDIFM(JL,JZ) > RRIINFTY ) ) THEN
            PDIFS(JL,JZ) = Z0
          ELSE
            PDIFS(JL,JZ) = Z1
          ENDIF
        ENDDO

        DO JZ=1,KLEVO-1
          ZTMP(JL)     = PDIFM(JL,JZ)
          PDIFM(JL,JZ) = PDIFS(JL,JZ-1)*PDIFM(JL,0) + Z2 * PDIFM(JL,JZ) &
                       & + PDIFS(JL,JZ+1)*PDIFM(JL,JZ+1) 
          ZWT(JL)      = PDIFS(JL,JZ-1) + Z2 + PDIFS(JL,JZ+1)
          PDIFM(JL,JZ) = PDIFM(JL,JZ) / ZWT(JL)
          PDIFM(JL,0)  = ZTMP(JL) !value at JZ before smoothing
        ENDDO

        !--- KLEVO
        DO JZ=1,KLEVO-1
          ZTMP(JL)     = PDIFM(JL,JZ)
          PDIFM(JL,JZ) = PDIFS(JL,JZ-1)*PDIFM(JL,0) + Z2 * PDIFM(JL,JZ) 
          ZWT(JL)      = PDIFS(JL,JZ-1) + Z2 
          PDIFM(JL,JZ) = PDIFM(JL,JZ) / ZWT(JL)
        ENDDO

      ENDDO ! JM loop

! 3. Effective overall interior diffusivity
!    --------------------------------------

! eq. (28) in LMD94

      ZRRIINFTY_INV = 1.0_JPRB / RRIINFTY
      DO JZ= 1, KLEVO
        ZRIG(JL)   = MAX( PDIFM(JL,JZ) , Z0 )
        ZRATIO(JL) = MIN( ZRIG(JL) * ZRRIINFTY_INV , Z1 )
        ZFRI(JL)   = Z1 - ZRATIO(JL) * ZRATIO(JL)
        ZFRI(JL)   = ZFRI(JL) * ZFRI(JL) * ZFRI(JL)

!       eq. (25) in LMD94
        PDIFM(JL,JZ) = RDIFM_IW + ZFRI(JL) * RDIFM_MAX 
        PDIFS(JL,JZ) = RDIFS_IW + ZFRI(JL) * RDIFS_MAX 
        PDIFT(JL,JZ) = PDIFS(JL,JZ)
      ENDDO

    ENDIF
  ENDDO

ELSE
! Reset the mixing coefficients
  DO JL = KIDIA, KFDIA
    IF( LDKPPCAL(JL) ) THEN
      DO JZ=0,KLEVOP1
        PDIFM(JL,JZ) = Z0
        PDIFT(JL,JZ) = Z0
        PDIFS(JL,JZ) = Z0
      ENDDO
    ENDIF
  ENDDO
ENDIF

! 4. Sets coefficients for double-diffusive mixing
!    ---------------------------------------------

IF(LDD_KPP) THEN  

  IF(LSF_NEW_KPP) THEN
    ZRRRHO01 = 1.0_JPRB / ( RRRHO0_D06 - Z1 )
  ELSE
    ZRRRHO01 = 1.0_JPRB / ( RRRHO0 - Z1 )
  ENDIF

  DO JL = KIDIA, KFDIA
    IF( LDKPPCAL(JL) ) THEN

      DO JZ = 1, KLEVO

        IF( ( PALPHADT(JL,JZ) > PBETADS(JL,JZ) ) &
          & .AND. ( PBETADS(JL,JZ) > Z0 ) ) THEN

          ! salt fingering case

          IF(LSF_NEW_KPP) THEN
            !Danabasoglu et al. 2006
            ZRRHO(JL) = MIN( PALPHADT(JL,JZ) / PBETADS(JL,JZ) , RRRHO0_D06 )
            ZDIFFDD(JL) = Z1 - ( (ZRRHO(JL)-Z1) * ZRRRHO01 )**3
            PDIFT(JL,JZ) = PDIFT(JL,JZ) + ZDIFFDD(JL) * 0.7_JPRB
            PDIFS(JL,JZ) = PDIFS(JL,JZ) + ZDIFFDD(JL)
          ELSE
            ! eq.(31) in LMD94
            ZRRHO(JL) = MIN( PALPHADT(JL,JZ) / PBETADS(JL,JZ) , RRRHO0 )
!            ZDIFFDD(JL) = Z1 - ( (ZRRHO(JL)-Z1) / ( RRRHO0-Z1 ) )**2
            ZDIFFDD(JL) = Z1 - ( (ZRRHO(JL)-Z1) * ZRRRHO01 )**2
            ZDIFFDD(JL) = RDSFMAX * ZDIFFDD(JL) * ZDIFFDD(JL) * ZDIFFDD(JL)
            PDIFT(JL,JZ) = PDIFT(JL,JZ) &
                         & + ZDIFFDD(JL) * 0.8_JPRB / ZRRHO(JL) 
                                                      !Changed from LMD94
            PDIFS(JL,JZ) = PDIFS(JL,JZ) + ZDIFFDD(JL)
          ENDIF

        ELSE IF ( ( PALPHADT(JL,JZ) < Z0 ) .AND. ( PBETADS(JL,JZ) < Z0 ) &
                & .AND. ( PALPHADT(JL,JZ) < PBETADS(JL,JZ) ) ) THEN

          ! diffusive convection
          ZRRHO(JL) = PALPHADT(JL,JZ) / PBETADS(JL,JZ)
          ! eq.(32) in LMD94
          ZDIFFDD(JL)  = RDMOL * 0.909_JPRB &
                       & *EXP(4.6_JPRB*EXP(-0.54_JPRB*(Z1/ZRRHO(JL)-Z1)))
          IF ( ZRRHO(JL) < 0.5_JPRB ) THEN
            ZPRANDTL(JL) = 0.15_JPRB * ZRRHO(JL)
          ELSE
            ZPRANDTL(JL) = ( 1.85_JPRB - 0.85_JPRB/ZRRHO(JL) ) * ZRRHO(JL)
          ENDIF

          PDIFT(JL,JZ) = PDIFT(JL,JZ) + ZDIFFDD(JL)
          PDIFS(JL,JZ) = PDIFS(JL,JZ) + ZPRANDTL(JL) * ZDIFFDD(JL)

        ENDIF

      ENDDO !JZ

    ENDIF
  ENDDO !JL

ENDIF

! 5. Reset surface values, fill the bottom KLEVOP1 coefficients for KPP_BLMIX
!    ------------------------------------------------------------------------

DO JL = KIDIA, KFDIA
  IF( LDKPPCAL(JL) ) THEN
    PDIFM(JL,0) = Z0
    PDIFT(JL,0) = Z0
    PDIFS(JL,0) = Z0
  ENDIF
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('KPP_INTERIOR_MIX_MOD:KPP_INTERIOR_MIX',1,ZHOOK_HANDLE)

END SUBROUTINE KPP_INTERIOR_MIX
END MODULE KPP_INTERIOR_MIX_MOD
