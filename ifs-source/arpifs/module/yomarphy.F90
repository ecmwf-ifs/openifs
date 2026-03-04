! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE YOMARPHY

USE PARKIND1, ONLY : JPIM

IMPLICIT NONE

SAVE

TYPE :: TARPHY
!*
!    -----------------------------------------------------------------

!    VARIABLES DE CONTROLE DE LA PHYSIQUE AROME :
!    CONTROL VARIABLES FOR THE AROME PHYSICS:

LOGICAL :: LMPA        ! Global switch for AROME physics coming from MesoNH (turbulence, microphysic, shallow-convection)
LOGICAL :: LMSE        ! Global switch for surfex use

LOGICAL :: LMICRO      ! Global switch for AROME microphysics 
LOGICAL :: LTURB       ! Global switch for AROME turbulence
LOGICAL :: LKFBCONV    ! Global switch for Kain-Fritch Convection
LOGICAL :: LKFBD       ! Control key to KFB deep convection
LOGICAL :: LKFBS       ! Control key to KFB shallow convection
LOGICAL :: LMFSHAL     ! Control key for Shallow Convection Mass Flux scheme (Pergaud et al, 2008)
LOGICAL :: LUSECHEM    ! Contol key for calling the Gas Chemistry scheme
LOGICAL :: LORILAM     ! Contol key for calling the Aerosol Chemistry scheme 
LOGICAL :: LRDUST      ! Contol key for calling the desertic aerosols scheme
LOGICAL :: LRCO2       ! Contol key for calling the CO2 gas concentration
LOGICAL :: LINITCHEM   ! Contol key to initialise  Gas Chemistry concentration
LOGICAL :: LINITORILAM ! Contol key to initialise  Aerosol Chemistry concentration
LOGICAL :: LINITDUST   ! Contol key to initialise  Desertic aerosol concentration
LOGICAL :: LRDEPOS     ! Contol key for calling aerosol scavenging
LOGICAL :: LBUFLUX     ! If TRUE fluxes are calculated in AROEND_BUDGET, if FALSE, tendencies remain
CHARACTER(LEN=1) :: CCOUPLING ! Type of SURFEX coupling. E - explicit, I - implicit
LOGICAL :: LMDUST      ! Contol key for calling the desertic aerosols IN ALADIN
LOGICAL :: LSURFEX_CRITICAL ! Run SURFEX calculations in a critical section
LOGICAL :: LSURFEX_KFROM ! Use TRUE/FALSE in GOTO_SURFEX from within OpenMP (debug purposes)
INTEGER(KIND=JPIM) :: NSURFEX_ITER  !Run SURFEX calculations several times (test harness)

END TYPE TARPHY


!    -------------------------------------------------------------------
END MODULE YOMARPHY

