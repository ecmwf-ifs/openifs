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

MODULE YOMCFU

USE PARKIND1,    ONLY : JPIM, JPRB
USE TYPE_FLUXES, ONLY : FLUXES_DESCRIPTOR
USE PTRGFU, ONLY : TCFUPTR

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------

!*    Contains variables to control activation of cumulated fluxes.

INTEGER(KIND=JPIM), PARAMETER :: JPFUST=40    ! maximum number of timesteps where CFU can be activated
INTEGER(KIND=JPIM), PARAMETER :: JPMXCFU=99   ! maximum number of CFU fields

TYPE :: TCFU_KEYS

LOGICAL :: LCUMFU  = .FALSE.                  ! controls switch on/off all CFU
LOGICAL :: LSTRD   = .FALSE.                  ! activates gravity wave drag momentum CFU if .T.
LOGICAL :: LSTRC   = .FALSE.                  ! activates contribution of convection to U, V, q and (cp T) CFU if .T.
LOGICAL :: LSTRT   = .FALSE.                  ! activates contribution of turbulence to U, V, q and (cp T) CFU if .T.
LOGICAL :: LFPLC   = .FALSE.                  ! activates convective precipitation CFU if .T.
LOGICAL :: LFPLCG  = .FALSE.                  ! activates convective graupels CFU if .T.
LOGICAL :: LFPLCH  = .FALSE.                  ! activates convective hail CFU if .T.
LOGICAL :: LFPLS   = .FALSE.                  ! activates stratiform precipitation CFU if .T.
LOGICAL :: LFPLSG  = .FALSE.                  ! activates stratiform graupels CFU if .T.
LOGICAL :: LFPLSH  = .FALSE.                  ! activates stratiform hail CFU if .T.
LOGICAL :: LFR     = .FALSE.                  ! activates radiation CFU if .T.
LOGICAL :: LAMIP   = .FALSE.                  ! activates AMIP output if .T.
LOGICAL :: LRAYS   = .FALSE.                  ! activates more radiative CFU if .T.
LOGICAL :: LRAYD   = .FALSE.                  ! activates downwards surface radiative CFU if .T.
LOGICAL :: LNEBTT  = .FALSE.                  ! activates total cloudiness CFU if .T.
LOGICAL :: LFSF    = .FALSE.                  ! activates surface CFU if .T.
LOGICAL :: LFSOIL  = .FALSE.                  ! activates soil CFU if .T.
LOGICAL :: LNEBPAR = .FALSE.                  ! activates partial cloudiness CFU if .T.
LOGICAL :: LTSTRD  = .FALSE.                  ! activates gravity wave drag momentum CFU at all levels if .T.
LOGICAL :: LTSTRC  = .FALSE.                  ! activates contribution of convection to U, V, q and (cp T) CFU at all levels if .T.
LOGICAL :: LTSTRT  = .FALSE.                  ! activates contribution of turbulence to U, V, q and (cp T) CFU at all levels if .T.
LOGICAL :: LTFPLC  = .FALSE.                  ! activates convective precipitation CFU at all levels if .T.
LOGICAL :: LTFPLS  = .FALSE.                  ! activates stratiform precipitation CFU at all levels if .T.
LOGICAL :: LTFR    = .FALSE.                  ! activates radiation CFU at all levels if .T.
LOGICAL :: LTNEB   = .FALSE.                  ! activates cloudiness CFU at all levels if .T.
LOGICAL :: LFDUTP  = .FALSE.                  ! activates filtered duration of total precipitations CFU if .T.
LOGICAL :: LMOON   = .FALSE.                  ! activates moon radiation CFU if .T.
LOGICAL :: LFRRC   = .FALSE.                  ! activates clear sky radiation calculation if .T.
LOGICAL :: LFLASH  = .FALSE.                  ! activates diagnostics of lightning

END TYPE TCFU_KEYS

TYPE, EXTENDS(TCFU_KEYS) :: TCFU

INTEGER(KIND=JPIM) :: NCFUTS(0:JPFUST)        ! array containing flux accumulation write-up steps
INTEGER(KIND=JPIM) :: NFRRC                   ! frequency for clear sky radiation calculation
INTEGER(KIND=JPIM) :: NFRCFU                  ! frequency of write up of flux diagnostics
INTEGER(KIND=JPIM) :: NFDCFU                  ! total number of fields in buffer
INTEGER(KIND=JPIM) :: NTYPCFU                 ! number of fluxes types in buffer
INTEGER(KIND=JPIM) :: NMTFLASH                ! method used to compute lightening density
REAL(KIND=JPRB) :: CALFLASH1,CALFLASH2        ! calibration factor for lightening density


LOGICAL :: LREACFU = .FALSE.                  ! read first input on historic file if .T.

TYPE(FLUXES_DESCRIPTOR) :: TYPE_CFU(JPMXCFU)  ! contains the fluxes descriptor for the CFU

TYPE(TCFUPTR) :: YCFUPT

REAL (KIND=JPRB), ALLOCATABLE :: GFUBUF (:,:,:)   ! Buffer for cumulative diagnostics

END TYPE TCFU

!     ------------------------------------------------------------
END MODULE YOMCFU
