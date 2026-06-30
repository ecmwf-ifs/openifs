MODULE KPP_WSCALE_MOD
CONTAINS
SUBROUTINE KPP_WSCALE &
 & ( KIDIA    ,KFDIA    ,KLON     ,LDKPPCAL ,PSIGMA   ,&
 &   PHBL     ,PUSTAR   ,PBFSFC   ,PWM      ,PWS      ,&
 &   YDOCEAN_ML )

! (C) Copyright 1994- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! Purpose :
! -------
!   This routine compute turbulent velocity scales.

! Interface :
! ---------
!   Call *KPP_WSCALE* from *KPP_BLDEPTH*

! Method :
! ------
!   This subroutine uses a 2d-lookup table for wm and ws as functions of 
!   ustar and zetahat (=vonk*sigma*hbl*bfsfc) for computational efficiency.

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

REAL(KIND=JPRB),INTENT(IN)  :: PSIGMA(KLON) ! normalized depth          (d/hbl)
REAL(KIND=JPRB),INTENT(IN)  :: PHBL(KLON)   ! boundary layer depth          (m)
REAL(KIND=JPRB),INTENT(IN)  :: PUSTAR(KLON) ! surface friction velocity   (m/s)
REAL(KIND=JPRB),INTENT(IN)  :: PBFSFC(KLON) ! total sfc buoyancy flux (m^2/s^3)
REAL(KIND=JPRB),INTENT(OUT) :: PWM(KIDIA:KFDIA)! turb. vel. scales at sigma
REAL(KIND=JPRB),INTENT(OUT) :: PWS(KIDIA:KFDIA)! ...        scalar ...

LOGICAL,INTENT(IN) :: LDKPPCAL(KLON)

TYPE(TOCEAN_ML), INTENT(INOUT) :: YDOCEAN_ML

INTEGER(KIND=JPIM) :: JL       ! loop control 
INTEGER(KIND=JPIM) :: JZ       ! ...
INTEGER(KIND=JPIM) :: JU       ! ...

INTEGER(KIND=JPIM) :: IZ       ! index for lookup table
INTEGER(KIND=JPIM) :: IZP1     ! ...
INTEGER(KIND=JPIM) :: IU       ! ...
INTEGER(KIND=JPIM) :: IUP1     ! ...

REAL(KIND=JPRB) :: ZZFRAC      ! temporary variable for lookup table
REAL(KIND=JPRB) :: ZUFRAC      ! temporary variable for lookup table
REAL(KIND=JPRB) :: ZFFRAC      ! temporary variable for lookup table
REAL(KIND=JPRB) :: ZWAM        ! temporary variable for lookup table
REAL(KIND=JPRB) :: ZWBM        ! temporary variable for lookup table
REAL(KIND=JPRB) :: ZWAS        ! temporary variable for lookup table
REAL(KIND=JPRB) :: ZWBS        ! temporary variable for lookup table
REAL(KIND=JPRB) :: ZUCUBE      ! ustar**3
REAL(KIND=JPRB) :: ZDELTAZ     ! delta zehat in table
REAL(KIND=JPRB) :: ZDELTAU     ! delta ustar in table
REAL(KIND=JPRB) :: ZDELTAZ_INV
REAL(KIND=JPRB) :: ZDELTAU_INV
REAL(KIND=JPRB) :: ZEHAT       ! = zeta *  ustar**3
REAL(KIND=JPRB) :: ZETA        ! = stability parameter d/l
REAL(KIND=JPRB) :: ZUSTA       ! ustar for lookup table
REAL(KIND=JPRB) :: ZZDIFF      ! temporary variable for looup table
REAL(KIND=JPRB) :: ZUDIFF      ! temporary variable for looup table
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('KPP_WSCALE_MOD:KPP_WSCALE',0,ZHOOK_HANDLE)
ASSOCIATE(LINIT_WSCALE=>YDOCEAN_ML%LINIT_WSCALE, NUTBLO=>YDOCEAN_ML%NUTBLO, &
 & NZTBLO=>YDOCEAN_ML%NZTBLO, RAM_KPP=>YDOCEAN_ML%RAM_KPP, &
 & RAS_KPP=>YDOCEAN_ML%RAS_KPP, RC1_KPP=>YDOCEAN_ML%RC1_KPP, &
 & RC2_KPP=>YDOCEAN_ML%RC2_KPP, RC3_KPP=>YDOCEAN_ML%RC3_KPP, &
 & RCM_KPP=>YDOCEAN_ML%RCM_KPP, RCS_KPP=>YDOCEAN_ML%RCS_KPP, &
 & REPS_KPP=>YDOCEAN_ML%REPS_KPP, RUMAX_KPP=>YDOCEAN_ML%RUMAX_KPP, &
 & RUMIN_KPP=>YDOCEAN_ML%RUMIN_KPP, RVONK=>YDOCEAN_ML%RVONK, &
 & RWMT=>YDOCEAN_ML%RWMT, RWST=>YDOCEAN_ML%RWST, RZETAM=>YDOCEAN_ML%RZETAM, &
 & RZETAS=>YDOCEAN_ML%RZETAS, RZETMAX_KPP=>YDOCEAN_ML%RZETMAX_KPP, &
 & RZETMIN_KPP=>YDOCEAN_ML%RZETMIN_KPP)

! 1. Construct the wm and ws lookup tables
!    -------------------------------------

ZDELTAZ = ( RZETMAX_KPP - RZETMIN_KPP ) / REAL(NZTBLO+1) 
ZDELTAZ_INV = 1.0_JPRB / ZDELTAZ 
ZDELTAU = ( RUMAX_KPP - RUMIN_KPP) / REAL(NUTBLO+1)
ZDELTAU_INV = 1.0_JPRB / ZDELTAU

IF(LINIT_WSCALE) THEN

  DO JZ=0, NZTBLO+1

    ZEHAT = ZDELTAZ * REAL(JZ) + RZETMIN_KPP

    DO JU=0, NUTBLO+1

      ZUSTA = ZDELTAU * REAL(JU) + RUMIN_KPP
      ZETA = ZEHAT / ( ZUSTA**3 + REPS_KPP )

      IF(ZEHAT >= 0.0_JPRB) THEN
        RWMT(JZ,JU) = RVONK * ZUSTA / ( 1.0_JPRB + RC1_KPP * ZETA )
        RWST(JZ,JU) = RWMT(JZ,JU)

      ELSE
        IF(ZETA > RZETAM) THEN
          RWMT(JZ,JU) = RVONK * ZUSTA &
                      & * ABS( 1.0_JPRB - RC2_KPP*ZETA )**(1.0_JPRB/4.0_JPRB)
        ELSE
          RWMT(JZ,JU) = RVONK * ( RAM_KPP*ZUSTA**3.0_JPRB &
                      &          - RCM_KPP*ZEHAT )**(1.0_JPRB/3.0_JPRB)
        ENDIF

        IF(ZETA > RZETAS) THEN
          RWST(JZ,JU) = RVONK * ZUSTA &
                      & * ABS(1.0_JPRB-RC3_KPP*ZETA)**(1.0_JPRB/2.0_JPRB)
        ELSE
          RWST(JZ,JU) = RVONK * ( RAS_KPP*ZUSTA**3 &
                      &          - RCS_KPP*ZEHAT )**(1.0_JPRB/3.0_JPRB)
        ENDIF

      ENDIF   

    ENDDO

  ENDDO

ENDIF ! LINIT_WSCALE

! 2. Use lookup table for turbulent scale
!    ------------------------------------

! use the table for zehat < zmax  only
! otherwise use stable formulae

DO JL = KIDIA, KFDIA
  IF ( LDKPPCAL(JL) ) THEN

    ZEHAT = RVONK * PSIGMA(JL) * PHBL(JL) * PBFSFC(JL)

    IF ( ZEHAT <= RZETMAX_KPP ) THEN

      ZZDIFF  = ZEHAT - RZETMIN_KPP
      IZ = INT( ZZDIFF * ZDELTAZ_INV )
      IZ = MIN( IZ, NZTBLO )
      IZ = MAX( IZ, 0  )
      IZP1 = IZ + 1

      ZUDIFF = PUSTAR(JL) - RUMIN_KPP
      IU = INT( ZUDIFF * ZDELTAU_INV )
      IU = MIN( IU, NUTBLO )
      IU = MAX( IU, 0 )
      IUP1 = IU + 1

      ZZFRAC = ZZDIFF * ZDELTAZ_INV - REAL(IZ)
      ZUFRAC = ZUDIFF * ZDELTAU_INV - REAL(IU)

      ZFFRAC= 1.0_JPRB - ZZFRAC
      ZWAM = ZFFRAC * RWMT(IZ,IUP1) + ZZFRAC * RWMT(IZP1,IUP1)
      ZWBM = ZFFRAC * RWMT(IZ,IU) + ZZFRAC * RWMT(IZP1,IU)
      ZWAS = ZFFRAC * RWST(IZ,IUP1) + ZZFRAC*RWST(IZP1,IUP1)
      ZWBS = ZFFRAC * RWST(IZ,IU) + ZZFRAC*RWST(IZP1,IU)

      ZFFRAC= 1.0_JPRB - ZUFRAC
      PWM(JL)  = ZFFRAC * ZWBM  + ZUFRAC * ZWAM
      PWS(JL)  = ZFFRAC * ZWBS  + ZUFRAC * ZWAS

    ELSE
 
      ZUCUBE = PUSTAR(JL) * PUSTAR(JL) * PUSTAR(JL)
      PWM(JL) = RVONK * PUSTAR(JL) * ZUCUBE &
              & / ( ZUCUBE + RC1_KPP * ZEHAT )
      PWS(JL) = PWM(JL)

    ENDIF   

  ENDIF
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('KPP_WSCALE_MOD:KPP_WSCALE',1,ZHOOK_HANDLE)

END SUBROUTINE KPP_WSCALE
END MODULE KPP_WSCALE_MOD
