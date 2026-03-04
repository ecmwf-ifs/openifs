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

TYPE :: TCFU
INTEGER(KIND=JPIM) :: NCFUTS(0:JPFUST)        ! array containing flux accumulation write-up steps
INTEGER(KIND=JPIM) :: NFRRC                   ! frequency for clear sky radiation calculation
INTEGER(KIND=JPIM) :: NFRCFU                  ! frequency of write up of flux diagnostics
INTEGER(KIND=JPIM) :: NFDCFU                  ! total number of fields in buffer
INTEGER(KIND=JPIM) :: NTYPCFU                 ! number of fluxes types in buffer
INTEGER(KIND=JPIM) :: NMTFLASH                ! method used to compute lightening density
REAL(KIND=JPRB)    :: CALFLASH1,CALFLASH2     ! calibration factor for lightening density

TYPE(FLUXES_DESCRIPTOR) :: TYPE_CFU(JPMXCFU)  ! contains the fluxes descriptor for the CFU

LOGICAL :: LCUMFU                             ! controls switch on/off all CFU
LOGICAL :: LREACFU                            ! read first input on historic file if .T.
LOGICAL :: LSTRD                              ! activates gravity wave drag momentum CFU if .T.
LOGICAL :: LSTRC                              ! activates contribution of convection to U, V, q and (cp T) CFU if .T.
LOGICAL :: LSTRT                              ! activates contribution of turbulence to U, V, q and (cp T) CFU if .T.
LOGICAL :: LFPLC                              ! activates convective precipitation CFU if .T.
LOGICAL :: LFPLCG                             ! activates convective graupels CFU if .T.
LOGICAL :: LFPLCH                             ! activates convective hail CFU if .T.
LOGICAL :: LFPLS                              ! activates stratiform precipitation CFU if .T.
LOGICAL :: LFPLSG                             ! activates stratiform graupels CFU if .T.
LOGICAL :: LFPLSH                             ! activates stratiform hail CFU if .T.
LOGICAL :: LFR                                ! activates radiation CFU if .T.
LOGICAL :: LAMIP                              ! activates AMIP output if .T.
LOGICAL :: LRAYS                              ! activates more radiative CFU if .T.
LOGICAL :: LRAYD                              ! activates downwards surface radiative CFU if .T.
LOGICAL :: LNEBTT                             ! activates total cloudiness CFU if .T.
LOGICAL :: LFSF                               ! activates surface CFU if .T.
LOGICAL :: LFSOIL                             ! activates soil CFU if .T.
LOGICAL :: LNEBPAR                            ! activates partial cloudiness CFU if .T.
LOGICAL :: LTSTRD                             ! activates gravity wave drag momentum CFU at all levels if .T.
LOGICAL :: LTSTRC                             ! activates contribution of convection to U, V, q and (cp T) CFU at all levels if .T.
LOGICAL :: LTSTRT                             ! activates contribution of turbulence to U, V, q and (cp T) CFU at all levels if .T.
LOGICAL :: LTFPLC                             ! activates convective precipitation CFU at all levels if .T.
LOGICAL :: LTFPLS                             ! activates stratiform precipitation CFU at all levels if .T.
LOGICAL :: LTFR                               ! activates radiation CFU at all levels if .T.
LOGICAL :: LTNEB                              ! activates cloudiness CFU at all levels if .T.
LOGICAL :: LFDUTP                             ! activates filtered duration of total precipitations CFU if .T.
LOGICAL :: LMOON                              ! activates moon radiation CFU if .T.
LOGICAL :: LFRRC                              ! activates clear sky radiation calculation if .T.
LOGICAL :: LFLASH                             ! activates diagnostics of lightning

TYPE(TCFUPTR) :: YCFUPT

REAL (KIND=JPRB), ALLOCATABLE :: GFUBUF (:,:,:)   ! Buffer for cumulative diagnostics

END TYPE TCFU

!!TYPE(TCFU), POINTER :: YRCFU => NULL()

!     ------------------------------------------------------------
END MODULE YOMCFU
