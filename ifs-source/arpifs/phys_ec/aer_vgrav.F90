! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      SUBROUTINE AER_VGRAV(KIK,PDP,PSIGMA,PDVISC,PMFPA,PRHOP,PVGRAV_AV)

!**** *AER_VGRAV* -  COMPUTES THE AVERAGE GRAVITATIONAL SETTLING
!                    VELOCITY

!**   DESCRIPTION 
!     ----------
!
! Calculates the average (over kth moment) gravitational
! settling velocity for a log-normally distributed aerosol
! population with geometric mean diameter DP and geometric
! mean standard deviation SIGMA following method in Regional
! Particulate Model as described by Binkowski & Shankar (1995).

!**   INTERFACE.
!     ----------
!          *AER_VGRAV* IS CALLED FROM *AER_DRYDEPVEL*.

! INPUTS:
! -------
! KIK              : Index of moment for calculation
! PDP             : Geometric mean particle diameter for mode (m)
! PSIGMA          : Geometric standard deviation for mode
! PDVISC          : Dynamic viscosity of air (kg m-1 s-1)
! PMFPA           : Mean free path of air (m)
! PRHOP           : Particle density (kg.m-3)


! OUTPUTS:
! --------
! VGRAV_AV : Avg. grav. settling velocity (m s-1)

! LOCAL VARIABLES
! ---------------
! ZKNG        : Knudsen number for geo. mean sized particle
! ZPREF       : Prefactor term to expression
! ZLNSQSG     : ln(SIGMA)*ln(SIGMA)

USE YOMCST,  ONLY: RG
USE PARKIND1, ONLY: JPRB, JPIM
USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

! Subroutine interface
INTEGER(KIND=JPIM), INTENT(IN) :: KIK
REAL(KIND=JPRB), INTENT(IN)    :: PDP
REAL(KIND=JPRB), INTENT(IN)    :: PSIGMA
REAL(KIND=JPRB), INTENT(IN)    :: PDVISC
REAL(KIND=JPRB), INTENT(IN)    :: PMFPA
REAL(KIND=JPRB), INTENT(IN)    :: PRHOP

REAL(KIND=JPRB), INTENT(OUT)    :: PVGRAV_AV

! Local variables
REAL(KIND=JPRB) :: ZKNG
REAL(KIND=JPRB) :: ZPREF
REAL(KIND=JPRB) :: ZLNSQSG

REAL(KIND=JPHOOK)               :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('AER_VGRAV',0,ZHOOK_HANDLE)

ZLNSQSG=LOG(PSIGMA)*LOG(PSIGMA)
!ZKNG=2.0_JPRB*PMFPA/PDP
ZKNG=(1._JPRB +( 2.0_JPRB*PMFPA/PDP)*(1.257_JPRB + 0.4_JPRB*exp(-1.1_JPRB * PDP/ PMFPA)))
ZPREF=PRHOP*PDP*PDP*RG/(18.0_JPRB*PDVISC)
PVGRAV_AV=ZPREF*ZKNG      ! reduced by 50_JPRB %


IF (LHOOK) CALL DR_HOOK('AER_VGRAV',1,ZHOOK_HANDLE)
END SUBROUTINE AER_VGRAV
