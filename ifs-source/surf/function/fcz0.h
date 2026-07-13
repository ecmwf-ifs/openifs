! (C) Copyright 2013- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!     ------------------------------------------------------------------

!     *FCZ0** CONTAINS STATEMENT FUNCTIONS DESCRIBING ROUGNESS LENGTHS

!     JEAN BIDLOT    E.C.M.W.F.      24/10/2013.
!     JEAN BIDLOT    E.C.M.W.F.      10/04/2015. sea-state dependent latent and sensible transfer coefficients
!                                                after Janssen( Tech memo 239, 1997) 

!     ------------------------------------------------------------------

!     FCZ0 DESCRIBES THE FUNCTIONAl DEPENDENCE OF THE ROUGHNESS LENGTHS 
!     ON THE ENVIROMENT 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: PRG     ! GRAVITY

! AERODYNAMIC ROUGHNESS LENGTH OVER SEAS AND OCEANS:
! --------------------------------------------------
REAL(KIND=JPRB) :: PCHARNOCK ! CHARNOCK PARAMETER
REAL(KIND=JPRB) :: PUST2     ! FRICTION VELOCITY SQUARED  

REAL(KIND=JPRB) :: PZ0SEA    ! ROUGHNESS LENGTH FOR MOMEMTUM OVER SEA DUE TO WAVES
PZ0SEA(PRG, PCHARNOCK, PUST2) = (PCHARNOCK/PRG)*PUST2


! ROUGHNESS LENGTH FOR HEAT AND HUMIDITY OVER SEAS AND OCEANS:
! ------------------------------------------------------------
REAL(KIND=JPRB), PARAMETER :: Z0HQMIN=0.0000001_JPRB
REAL(KIND=JPRB) :: PZ0       ! ROUGHNESS LENGTH FOR MOMEMTUM OVER SEA 
REAL(KIND=JPRB) :: PZ0HQ     ! ROUGHNESS LENGTH FOR MOMEMTUM OVER SEA WITHOUT THE OCEAN WAVES CONTRIBUTION
REAL(KIND=JPRB) :: PZN       ! ROUGHNESS LENGTH FOR HEAT OR HUMIDITY OVER SEA 
REAL(KIND=JPRB) :: PZP
REAL(KIND=JPRB) :: PZM

REAL(KIND=JPRB) :: PZPLUS    ! LARGEST ROOT OF Z**2 + (PZN+2*PZ0HQ)*Z + PZN*PZ0HQ = 0 
PZPLUS(PZ0HQ, PZN) = -(PZ0HQ+0.5_JPRB*PZN)+SQRT(PZ0HQ**2+0.25_JPRB*PZN**2) 

REAL(KIND=JPRB) :: PZMINS    ! SMALLEST ROOT OF Z**2 + (PZN+2*PZ0HQ)*Z + PZN*PZ0HQ = 0 
PZMINS(PZ0HQ, PZN) = -(PZ0HQ+0.5_JPRB*PZN)-SQRT(PZ0HQ**2+0.25_JPRB*PZN**2) 

REAL(KIND=JPRB) :: PZZ       ! USEFUL FUNCTION
PZZ(PZ0HQ, PZP, PZM) = ABS(PZP)**((PZ0HQ+PZP)/(PZP-PZM))

REAL(KIND=JPRB) :: PZNSEA    ! ROUGHNESS LENGTH FOR HEAT OR HUMIDITY OVER SEA AND OCEANS 
PZNSEA(PZ0HQ, PZN, PZ0, PUST2) = MAX( PZZ(PZ0HQ, PZPLUS(PZ0HQ,PZN), PZMINS(PZ0HQ,PZN)) * &
                               &      PZZ(PZ0HQ, PZMINS(PZ0HQ,PZN), PZPLUS(PZ0HQ,PZN)),  &
                               &      Z0HQMIN)

! AERODYNAMIC ROUGHNESS LENGTH OVER SEA ICE:
! ------------------------------------------
REAL(KIND=JPRB), PARAMETER :: AZ0ICE=6.05E-3_JPRB
REAL(KIND=JPRB), PARAMETER :: BZ0ICE=17._JPRB
REAL(KIND=JPRB), PARAMETER :: CZ0ICE=0.93E-3_JPRB

REAL(KIND=JPRB) :: PCICCTR ! SEA ICE CONCENTRATION
REAL(KIND=JPRB) :: PZ0MIN  ! MINIMUM VALUE FOR PZ0ICE 
REAL(KIND=JPRB) :: PZ0ICE  ! ROUGHNESS LENGTH FOR MOMEMTUM OVER SEA ICE 

PZ0ICE(PZ0MIN, PCICCTR)=MAX(PZ0MIN,(1._JPRB-PCICCTR)*CZ0ICE+AZ0ICE*EXP(-BZ0ICE*(PCICCTR-0.5_JPRB)**2))

!     ------------------------------------------------------------------
