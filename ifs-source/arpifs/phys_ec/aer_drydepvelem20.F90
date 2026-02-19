! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      SUBROUTINE AER_DRYDEPVELEM20(KSEASON_WE,KVEG_ZH,PRHOP,PWETD,PSIGMA, &
      &    PZ0M,PCI,PUSTR,PDZ,PT, PRHO, PVDEP)

!**** *AER_DRYDEPVEL* -  ROUTINE FOR PARAMETRIZATION OF DRY DEPOSITION VELOCITY

!**   DESCRIPTION 
!     ----------
! Calculates aerosol dry deposition (and sedimentation).
! Based on the parameterisation of Zhang et al (2001) which
! uses the method in the model of Slinn (1982).
!
! Evaluate deposition velocity in lowest level as:
!
! V_dep = 1/(A_r + S_r)
!
! where V_dep is the deposition velocity
!       V_g   is the gravitational velocity = rho_p*Dp^2*g*CF/(18*DVISC)
!       CF    is the Cunningham slip correction
!       DVISC is the dynamic viscosity
!       A_r   is the aerodynamic resitance
!       S_r   is the surface resistance
!       Dp    is the particle diameter
!       rho_p is the particle density
!       g     is the gravitational acceleration
!
! Evaluate S_r=1/{ 3 * ustar * (EB + EIM + EIN) }
!
! following parameterization by Zhang et al (2001) where
!
! EB,EIM,EIN are collection efficiencies for Brownian diffusion,
! impaction and interception respectively.
!
! EB = Sc^-YR where Sc is the particle Schmidt number = nu/D
!                                where nu = kinematic viscosity of air
!                                      D =  particle diffusion coeff.
!
!  and YR is surface-dependent constant, values as in Table 3 (Zhang01)
!         0.50 over water          (Land use category 13-14)
!         0.56 over forest         (Land use category  1- 5)
!         0.54 over grass and ice  (Land use category  6,12)
!
! EIM = { St/(ALPHA+St) }^2
!
!    where St is the Stokes number = V_g * ustar^2 / DVISC  (z0<1mm)
!                                  = V_g * ustar   / (g*CR) (z0>1mm)
!
!                                    [smooth & rough flow regimes]
!
!      and ALPHA,CR are surface-dependent constant, values as Table 3:
!         ALPHA=100.0, CR=0.0 over water [only divide by CR for veg]
!         ALPHA= 50.0, CR=0.0 over ice   [only divide by CR for veg]
!         ALPHA=  1.0, CR=0.005 over grass
!         ALPHA=  1.2, CR=0.002 over forest
!
! EIN = 0.5*Dp/CR


!**   INTERFACE.
!     ----------
!          *AER_DRYDEPVEL* IS CALLED FROM *AER_DRYDEP*.

! INPUTS:
! -------
! KSEASON_WE : Weseley et al season    (index)
! KVEG_ZH    : Zhand et al land class  (index)
! PRHOP      : PARTICLE DENSITY        (kg m-3)
! PWETD      : PARTICLE WET DIAMETER   (m)
! PSIGMA     : PARTICLE DISTRIBUTION STDEV 
! PZ0M       : ROUGHNESS LENGTH         (m)
! PCI        : SEA-ICE MASK
! PUST       : FRICTION VELOCITY        (m.s-1)
! PDZ        : DELTA Z                  (m)
! PT         : TEMPERATURE              (K)
! PRHO       : AIR DENSITY              (kg m-3)

! OUTPUTS:
! --------
! PVDEP      : DRY DEPOSITION VELOCITY  (m.s-1)

! LOCAL VARIABLES:
! ---------------
! ZPS_AV_3    : 3rd moment avg particle Schmidt Number
! ZVISC      :  viscosity of air (m2 s-1)
! ZKVISC      : Kinematic viscosity of air (m2 s-1)
! ZVGRAV_AV_3 : 3rd moment avg grav. settling vel. (m/s)
! ZVDEP_AV_3  : 3rd moment avg deposition velocity (m/s)
! ZDCOEF_AV_3 : 3rd moment avg particle diffusion coefficient(m2/s)
! ZSN_AV_3    : 3rd moment avg Stokes number
! ZSR_AV_3    : 3rd moment avg surface resistance
! ZEB_AV_3    : 3rd moment avg collection eff. for Brownian diffusion
! ZEIM_AV_3   : 3rd moment avg collection eff. for impaction
! ZEIN        : Collection eff. for interception
! ZAR         : Aerodynamic resistance
! ZCR,ZY,ZALPHA: aerosol deposition coefficients
!     [vary with land category & input via DATA statements]
! ZCR        : Characteristic radius of collectors (m)
! ZY         : Parameter for calculating Brownian diffusion
! ZALPHA     : Parameter for calculating EIM
!

!     EXTERNALS.
!     ----------
!    aer_vgrav.F90
!    aer_dcoff.F90

!     AUTHOR.
!     -------
!        SAMUEL REMY   *CNRS-IPSL*
!        GLOMAP community
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2016-09-15
!        SR         2017-01-17   externalize land class attribution

! 
!     REFERENCES.
!     ----------
! Slinn, Atmos. En., 1982, 16, 1785-1794
! Zhang et al, Atmos. En., 2001, 35, 549-560
!
!----------------------------------------------------------------------
USE YOMCST ,             ONLY : RG, RPI, RKBOL
USE PARKIND1,            ONLY: JPRB, JPIM
USE YOMHOOK,             ONLY: LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE


! .. Subroutine interface
INTEGER(KIND=JPIM), INTENT(IN) :: KSEASON_WE
INTEGER(KIND=JPIM), INTENT(IN) :: KVEG_ZH
REAL(KIND=JPRB), INTENT(IN) :: PRHOP
REAL(KIND=JPRB), INTENT(IN) :: PWETD
REAL(KIND=JPRB), INTENT(IN) :: PSIGMA
REAL(KIND=JPRB), INTENT(IN) :: PZ0M
REAL(KIND=JPRB), INTENT(IN) :: PCI
REAL(KIND=JPRB), INTENT(IN) :: PUSTR
REAL(KIND=JPRB), INTENT(IN) :: PDZ
REAL(KIND=JPRB), INTENT(IN) :: PT
REAL(KIND=JPRB), INTENT(IN) :: PRHO
REAL(KIND=JPRB), INTENT(OUT) :: PVDEP
!
!    Local Variables
REAL(KIND=JPRB)    :: ZPS_AV_3
REAL(KIND=JPRB)    :: ZKVISC, ZVISC,ZVBA,ZMFPA,ZKARMN
REAL(KIND=JPRB)    :: ZVGRAV_AV_3
REAL(KIND=JPRB)    :: ZDCOEF_AV_3
REAL(KIND=JPRB)    :: ZEB_AV_3
REAL(KIND=JPRB)    :: ZEIM_AV_3
REAL(KIND=JPRB)    :: ZEIN
REAL(KIND=JPRB)    :: ZSN_AV_3
REAL(KIND=JPRB)    :: ZAR
REAL(KIND=JPRB)    :: ZSR_AV_3
REAL(KIND=JPRB) :: ZCR(15,5)
REAL(KIND=JPRB), PARAMETER :: ZALPHA(15) = (/1.0_JPRB,0.6_JPRB,1.1_JPRB,0.8_JPRB,0.8_JPRB,         &
 &    1.2_JPRB,1.2_JPRB,50.0_JPRB,50.0_JPRB,1.3_JPRB,2.0_JPRB,50._JPRB,100._JPRB,100._JPRB,1.5_JPRB/)
REAL(KIND=JPRB), PARAMETER :: ZIMA=4.78E-26_JPRB

REAL(KIND=JPHOOK)             :: ZHOOK_HANDLE


#include "aer_vgrav.intfb.h"
#include "aer_dcoff.intfb.h"

IF (LHOOK) CALL DR_HOOK('AER_DRYDEPVELEM20',0,ZHOOK_HANDLE)


ZCR(:,1)=(/8E-4_JPRB,2e-3_JPRB,8E-4_JPRB,2E-3_JPRB,2E-3_JPRB,8E-4_JPRB,8E-4_JPRB,8E-4_JPRB,8E-4_JPRB,4E-3_JPRB,4E-3_JPRB,&
         & 4E-3_JPRB,4E-3_JPRB,4E-3_JPRB,4E-3_JPRB/)
ZCR(:,2)=(/8E-4_JPRB,2e-3_JPRB,8E-4_JPRB,2E-3_JPRB,2E-3_JPRB,8E-4_JPRB,8E-4_JPRB,8E-4_JPRB,8E-4_JPRB,4E-3_JPRB,4E-3_JPRB,&
         & 4E-3_JPRB,4E-3_JPRB,4E-3_JPRB,4E-3_JPRB/)
ZCR(:,3)=(/8E-4_JPRB,2e-3_JPRB,2E-3_JPRB,4E-3_JPRB,2E-3_JPRB,2E-3_JPRB,8E-4_JPRB,8E-4_JPRB,8E-4_JPRB,4E-3_JPRB,4E-3_JPRB,&
         & 4E-3_JPRB,4E-3_JPRB,4E-3_JPRB,4E-3_JPRB/)
ZCR(:,4)=(/8E-4_JPRB,2e-3_JPRB,2E-3_JPRB,4E-3_JPRB,2E-3_JPRB,2E-3_JPRB,8E-4_JPRB,8E-4_JPRB,8E-4_JPRB,4E-3_JPRB,4E-3_JPRB,&
         & 4E-3_JPRB,4E-3_JPRB,4E-3_JPRB,4E-3_JPRB/)
ZCR(:,5)=(/8E-4_JPRB,2e-3_JPRB,8E-4_JPRB,2E-3_JPRB,2E-3_JPRB,8E-4_JPRB,8E-4_JPRB,8E-4_JPRB,8E-4_JPRB,4E-3_JPRB,4E-3_JPRB,&
         & 4E-3_JPRB,4E-3_JPRB,4E-3_JPRB,4E-3_JPRB/)


!ZCR(:,1)=(/2E-3_JPRB,5e-3_JPRB,2E-4_JPRB,5E-3_JPRB,5E-3_JPRB,2E-3_JPRB,2E-3_JPRB,8E-4_JPRB,8E-4_JPRB,1E-2_JPRB,1E-2_JPRB,&
!         & 4E-3_JPRB,4E-3_JPRB,4E-3_JPRB,1E-2_JPRB/)
!ZCR(:,2)=(/2E-3_JPRB,5e-3_JPRB,2E-4_JPRB,5E-3_JPRB,5E-3_JPRB,2E-3_JPRB,2E-3_JPRB,8E-4_JPRB,8E-4_JPRB,1E-2_JPRB,1E-2_JPRB,&
!         & 4E-3_JPRB,4E-3_JPRB,4E-3_JPRB,1E-2_JPRB/)
!ZCR(:,3)=(/2E-3_JPRB,5e-3_JPRB,5E-4_JPRB,1E-2_JPRB,5E-3_JPRB,5E-3_JPRB,5E-3_JPRB,8E-4_JPRB,8E-4_JPRB,1E-2_JPRB,1E-2_JPRB,&
!         & 4E-3_JPRB,4E-3_JPRB,4E-3_JPRB,1E-2_JPRB/)
!ZCR(:,4)=(/2E-3_JPRB,5e-3_JPRB,5E-4_JPRB,1E-2_JPRB,5E-3_JPRB,5E-3_JPRB,5E-3_JPRB,8E-4_JPRB,8E-4_JPRB,1E-2_JPRB,1E-2_JPRB,&
!         & 4E-3_JPRB,4E-3_JPRB,4E-3_JPRB,1E-2_JPRB/)
!ZCR(:,5)=(/2E-3_JPRB,5e-3_JPRB,2E-4_JPRB,5E-3_JPRB,5E-3_JPRB,2E-3_JPRB,2E-3_JPRB,8E-4_JPRB,8E-4_JPRB,1E-2_JPRB,1E-2_JPRB,&
!         & 4E-3_JPRB,4E-3_JPRB,4E-3_JPRB,1E-2_JPRB/)


PVDEP=0.0_JPRB
! dynamic viscosity of air (kg m^-1 s^-1)
ZVISC = 1.83E-5_JPRB*(416.16_JPRB/(PT+120.0_JPRB))*(SQRT(PT/296.16_JPRB)**3.0_JPRB)
! mean free speed of air molecules
ZVBA  = SQRT(8.0_JPRB*RKBOL*PT)/(RPI*ZIMA) 
 ! mean free path of air molecules
ZMFPA   = 2.0_JPRB*ZVISC/(PRHO*ZVBA)

ZKARMN=0.4_JPRB
! .. Calculate aerodynamic resistance
ZAR=LOG(PDZ/PZ0M)/(ZKARMN*PUSTR)

!       Calculate 3rd moment avg. grav. settling velocities
CALL AER_VGRAV(3,PWETD,PSIGMA,ZVISC,ZMFPA,PRHOP,ZVGRAV_AV_3)

!       Calculate 3rd moment avg particle diffusion coeffs
CALL AER_DCOFF(3,PWETD,PSIGMA,PT,ZVISC,ZMFPA,ZDCOEF_AV_3)

!      Calculate kinematic viscosity of air
ZKVISC=ZVISC/PRHO

!      Calculate 0th and 3rd moment avg. particle Schmidt number
ZPS_AV_3=ZKVISC/ZDCOEF_AV_3
!      Calculate particle collection efficiencies
!       -- For Brownian Diffusion
ZEB_AV_3=0.2_JPRB*ZPS_AV_3**(-0.66_JPRB)

!       -- For Impaction
IF (KVEG_ZH==8 .OR. KVEG_ZH==9 .OR. KVEG_ZH==12 .OR. KVEG_ZH==13 .OR. KVEG_ZH==14 ) THEN
!        Calculate stokes number for smooth surfaces
  ZSN_AV_3=ZVGRAV_AV_3*PUSTR*PUSTR/ZVISC
ELSE
!        Calculate stokes number for vegetated surfaces
  ZSN_AV_3=ZVGRAV_AV_3*PUSTR/(RG*ZCR(KVEG_ZH,KSEASON_WE))
ENDIF

ZEIM_AV_3=0.4_JPRB*(ZSN_AV_3/(ZALPHA(KVEG_ZH)+ZSN_AV_3))**1.7_JPRB

!       -- For Interception
IF (KVEG_ZH==8 .OR. KVEG_ZH==9 .OR. KVEG_ZH==12 .OR. KVEG_ZH==13 .OR. KVEG_ZH==14 ) THEN
   ZEIN=0.0_JPRB
ELSE
   ZEIN=2.5_JPRB*(PWETD/ZCR(KVEG_ZH,KSEASON_WE))**0.8_JPRB
ENDIF

!       Calculate surface resistance
ZSR_AV_3=1.0_JPRB/(3.0_JPRB*PUSTR*(ZEB_AV_3+ZEIM_AV_3+ZEIN))
!       Calculate deposition velocity
PVDEP=1.0_JPRB/(ZAR+ZSR_AV_3)

IF (LHOOK) CALL DR_HOOK('AER_DRYDEPVELEM20',1,ZHOOK_HANDLE)
END SUBROUTINE AER_DRYDEPVELEM20
