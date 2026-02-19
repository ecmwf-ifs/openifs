! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE DDR_LAMINAR_RES_GC(PTS, PRESS, PXM,PUSTAR, PXDIMO, PWRB)

!!****  *SUBROUTINE COMP_LAMINAR_RES * - subroutine to compute laminar resistances for wesely scheme
!!       
!!
!!    INTERFACE.
!!    ----------
!!    Called from depvel_wsl.F90
!!  
!!    PURPOSE
!!    -------
!!    The purpose of this routine is to compute laminar resistances.
!!   
!!    METHOD
!!    ------
!!
!! INPUTS:
!! -------
!! PTS       : Surface temperature
!! PXM       : Molecular weight [kg/mole]
!! PRESS     : Surface pressure
!! PUSTAR    : U*
!! PXDIMO    : Molecular diffusivities (no longer needed)
!!
!! OUTPUTS:
!! -------
!! PWRB      : Laminar resistance
!!
!! LOCAL:
!! -------
!!
!!    AUTHOR
!!    ------
!!    M. Michou 
!!
!!
!!    MODIFICATIONS
!!    -------------
!!    Original   July 1999
!!    J Flemming 23.4.2018   Adapted for IFS 
!!
!!-------------------------------------------------------------------------------
  

USE PARKIND1 ,ONLY : JPIM, JPRB
USE DRYDEP_PAR   ,ONLY :PPRMAX, PPKARMAN
USE YOMCST, ONLY : R, RPI
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

REAL (KIND=JPRB), INTENT(IN)  :: PXDIMO
REAL (KIND=JPRB), INTENT(IN)  :: PUSTAR ! U*
REAL (KIND=JPRB), INTENT(IN)  :: PTS ! Surface temperature
REAL (KIND=JPRB), INTENT(IN)  :: PRESS ! Surface pressure
REAL (KIND=JPRB), INTENT(IN)  :: PXM ! Molecular weight of gas
REAL (KIND=JPRB), INTENT(OUT) :: PWRB ! Aerodynamic resistance

! DAIR is the thermal diffusivity of air; value of 0.02 cm2 s-1 cited on p. 16,476 of
! Jacob et al. [1992]
REAL(KIND=JPRB), PARAMETER  :: ZDAIR=0.2E-4 ! assumed units in m2 s-1
REAL(KIND=JPRB), PARAMETER  :: ZXMAIR  = 28.8E-3_JPRB ! Moist air molec wt
REAL(KIND=JPRB), PARAMETER  :: ZRADAIR = 1.2E-10_JPRB
REAL(KIND=JPRB), PARAMETER  :: ZRADX   = 1.5E-10_JPRB
REAL(KIND=JPRB), PARAMETER  :: ZAVOGAD = 6.023E+23_JPRB ! Avogadros constant
!! Local variables
REAL(KIND=JPRB)             :: ZAIRDEN, ZZ, ZDIAM, ZFRPATH, ZSPEED, ZDIFF_G

REAL(KIND=JPHOOK) ::  ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DDR_LAMINAR_RES_GC',0,ZHOOK_HANDLE)

! Calculate molecular diffusivity of H20 for calculation
! Air density [molec/m3]
ZAIRDEN = ( PRESS * ZAVOGAD ) / ( R * PTS )
! DIAM is the collision diameter for gas X with air.
ZDIAM   = ZRADX + ZRADAIR
!Calculate the mean free path for gas X in air:
! eq. 8.5 of Seinfeld [1986];
ZZ      = PXM  / ZXMAIR
ZFRPATH = 1._JPRB /( RPI * SQRT( 1._JPRB + ZZ ) * ZAIRDEN * ( ZDIAM**2_JPIM ) )

! Calculate average speed of gas X; eq. 15.47 of Levine [1988]
ZSPEED  = SQRT( 8_JPRB * R * PTS / ( RPI * PXM ) )

! Calculate diffusion coefficient of gas X in air;
! eq. 8.9 of Seinfeld [1986]
ZDIFF_G = ( 3_JPRB * RPI / 32_JPRB ) * ( 1_JPRB + ZZ ) * ZFRPATH * ZSPEED

! PXDIMO = Ratio between D_H2O / D_I MOLECULAR diffusivity
! Therefore D_I MOLECULAR = D_H20 / PXDIMO
PWRB  = (2.0_JPRB / (PPKARMAN * PUSTAR)) * (ZDAIR/ZDIFF_G) ** (2.0_JPRB / 3.0_JPRB)

! Resistances limited to a maximum value
IF (PWRB > PPRMAX) PWRB = PPRMAX

IF (LHOOK) CALL DR_HOOK('DDR_LAMINAR_RES_GC',1,ZHOOK_HANDLE)

END SUBROUTINE DDR_LAMINAR_RES_GC
