! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

     SUBROUTINE AER_DCOFF(KIK,PDP,PSIGMA,PT,PDVISC,PMFPA,PDCOFF_PAR_AV_K)

!**** *AER_DCOFF* -  COMPUTES THE AVERAGE DIFFUSION COEFFICIENT
!                    VELOCITY

!**   DESCRIPTION 
!     ----------
!
! Calculates the average (over kth moment) diffusion coefficient
! for a log-normally distributed aerosol population with
! geometric mean diameter DG, geometric mean standard
! deviation SIGMA following method in Regional Particulate
! Model as described by Binkowski & Shankar (1995).
! Also follows expression in Seinfeld & Pandis pg 474
! (Stokes-Einstein relation with slip-flow correction)

!**   INTERFACE.
!     ----------
!          *AER_DCOFF* IS CALLED FROM *AER_DRYDEPVEL*.

! INPUTS:
! -------
! KIK              : Index of moment for calculation
! PDP             : Geometric mean particle diameter for mode (m)
! PSIGMA          : Geometric standard deviation for mode
! PT              : Temperature of air (K)
! PDVISC          : Dynamic viscosity of air (kg m-1 s-1)
! PMFPA           : Mean free path of air (m)


! OUTPUTS:
! --------
! DCOFF_PAR_AV_K : Avg. ptcl diffusion coefficient (m^2 s-1)

! LOCAL VARIABLES
! ---------------
! ZKNG            : Knudsen number for geo. mean sized particle
! ZLNSQSG         : ln(SIGMA)*ln(SIGMA)
! ZPREF           : Prefactor term to expression
! ZTERM1,ZTERM2    : Terms in average diff coeff expression

!     EXTERNALS.
!     ----------

!     AUTHOR.
!     -------
!        SAMUEL REMY   *CNRS-IPSL*
!        GLOMAP community
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2016-09-15


!--------------------------------------------------------------------
!
! Purpose
! -------
!
!
!--------------------------------------------------------------------
USE YOMCST,   ONLY: RKBOL, RPI
USE YOMHOOK,          ONLY: LHOOK, DR_HOOK, JPHOOK
USE PARKIND1,         ONLY: JPRB, JPIM
IMPLICIT NONE
!
! Subroutine interface
INTEGER(KIND=JPIM), INTENT(IN) :: KIK
REAL(KIND=JPRB), INTENT(IN)    :: PDP
REAL(KIND=JPRB), INTENT(IN)    :: PSIGMA
REAL(KIND=JPRB), INTENT(IN)    :: PT
REAL(KIND=JPRB), INTENT(IN)    :: PDVISC
REAL(KIND=JPRB), INTENT(IN)    :: PMFPA
REAL(KIND=JPRB), INTENT(OUT)   :: PDCOFF_PAR_AV_K

! Local variables
REAL(KIND=JPRB)    :: ZLNSQSG
REAL(KIND=JPRB)    :: ZTERM1
REAL(KIND=JPRB)    :: ZTERM2
REAL(KIND=JPRB)    :: ZTERM3
REAL(KIND=JPRB)    :: ZTERM4
REAL(KIND=JPRB)    :: ZKNG
REAL(KIND=JPRB)    :: ZPREF

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('AER_DCOFF',0,ZHOOK_HANDLE)

ZLNSQSG=LOG(PSIGMA)*LOG(PSIGMA)
ZTERM1=(-2.0_JPRB*FLOAT(KIK)+1.0_JPRB)/2.0_JPRB
ZTERM2=(-4.0_JPRB*FLOAT(KIK)+4.0_JPRB)/2.0_JPRB
ZTERM3=EXP(ZTERM1*ZLNSQSG)
ZTERM4=1.246_JPRB*EXP(ZTERM2*ZLNSQSG)

ZKNG=2.0_JPRB*PMFPA/PDP
ZPREF=RKBOL*PT/(3.0_JPRB*RPI*PDVISC*PDP)
!PDCOFF_PAR_AV_K=ZPREF*(ZTERM3+ZTERM4*ZKNG)
PDCOFF_PAR_AV_K=ZPREF

IF (LHOOK) CALL DR_HOOK('AER_DCOFF',1,ZHOOK_HANDLE)
END SUBROUTINE AER_DCOFF
