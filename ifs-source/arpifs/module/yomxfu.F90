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

MODULE YOMXFU

USE PARKIND1  ,ONLY : JPIM, JPRB

USE TYPE_FLUXES, ONLY : FLUXES_DESCRIPTOR
USE PTRXFU, ONLY : TXFUPTR

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------

!*    Contains variables to control activation of instantaneous fluxes.

INTEGER(KIND=JPIM), PARAMETER :: JPFUXT=240  ! maximum number of timesteps where XFU can be activated
INTEGER(KIND=JPIM), PARAMETER :: JPMXXFU=201 ! maximum number of XFU fields

TYPE :: TXFU
TYPE(FLUXES_DESCRIPTOR) :: TYPE_XFU(JPMXXFU) ! contains the fluxes descriptor for the XFU

REAL(KIND=JPRB),ALLOCATABLE:: RMWINDCALC(:)  ! needed for mean wind calculation
REAL(KIND=JPRB),ALLOCATABLE:: RMNWINDCALC(:) ! needed for mean neutral wind calculation
INTEGER(KIND=JPIM) :: MEANPERIOD             ! period (in seconds) for the mean calculation
INTEGER(KIND=JPIM) :: NMEANSTEPS             ! number of timesteps involved in mean calculation
INTEGER(KIND=JPIM) :: NXGSTPERIOD            ! period for maximum gusts
INTEGER(KIND=JPIM) :: NXGSTPERIOD2           ! period for second maximum gusts
INTEGER(KIND=JPIM) :: NVISIPERIOD            ! period for visibilities
INTEGER(KIND=JPIM) :: NVISIPERIOD2           ! period for second visibilities
INTEGER(KIND=JPIM) :: NXGSTTS                ! number of timesteps involved in max gust calculation
INTEGER(KIND=JPIM) :: NTYPXFU                ! number of fluxes types in buffer
INTEGER(KIND=JPIM) :: NXFUTS(0:JPFUXT)       ! array containing flux accumulation write-up steps
INTEGER(KIND=JPIM) :: NFRXFU                 ! frequency of write up of flux diagnostics
INTEGER(KIND=JPIM) :: NRAZTS(0:JPFUXT)       ! array containing instantaneous flux reset steps
INTEGER(KIND=JPIM) :: NFRRAZ                 ! frequency of reset of flux diagnostics
INTEGER(KIND=JPIM) :: N1RAZ                  ! over-riding switch for instantaneous flux reset (0 = false)
INTEGER(KIND=JPIM) :: NFDXFU                 ! total number of fields in buffer

LOGICAL :: LXFU                              ! controls switch on/off all XFU
LOGICAL :: LRESET                            ! reset extreme temperatures to zero
LOGICAL :: LRESET_GST                        ! reset Gust calculation
LOGICAL :: LRESET_GST2                       ! reset Gust2 calculation
LOGICAL :: LRESET_PRECIP                     ! reset Precips type calcultation
LOGICAL :: LRESET_PRECIP2                    ! reset Precips type calcultation
LOGICAL :: LRESET_VISI                       ! reset visibilities calculations
LOGICAL :: LRESET_VISI2                      ! reset visibilities calculations

LOGICAL :: LREAXFU                           ! read first input on historic file if .T.
LOGICAL :: LXTRD                             ! activates gravity wave drag momentum XFU if .T.
LOGICAL :: LXTRC                             ! activates contribution of convection to U, V, q and (cp T) XFU if .T.
LOGICAL :: LXTRT                             ! activates contribution of turbulence to U, V, q and (cp T) XFU if .T.
LOGICAL :: LXPLC                             ! activates convective precipitation XFU if .T.
LOGICAL :: LXPLCG                            ! activates convective graupels CFU if .T.
LOGICAL :: LXPLCH                            ! activates convective hail CFU if .T.
LOGICAL :: LXPLS                             ! activates stratiform precipitation XFU if .T.
LOGICAL :: LXPLSG                            ! activates stratiform graupels CFU if .T.
LOGICAL :: LXPLSH                            ! activates stratiform hail CFU if .T.
LOGICAL :: LXR                               ! activates radiation XFU if .T.
LOGICAL :: LXNEBTT                           ! activates total cloudiness XFU if .T.
LOGICAL :: LXNEBPA                           ! activates partial cloudiness XFU if .T.
LOGICAL :: LXCLS                             ! activates U, V, T, q and relative humidity at 2 or 10 m (time t-dt) if .T.
LOGICAL :: LXMWINDCLS                        ! activates mean of U and V at 10 m if .T., also NU/NV if LXNUVCLS
LOGICAL :: LXNUVCLS                          ! activates neutral U and V at 10 m (time t-dt) if .T.
LOGICAL :: LXTTCLS                           ! activates extreme temperatures at 2 m if .T.
LOGICAL :: LXHHCLS                           ! activates extreme relative moistures at 2 m if .T
LOGICAL :: LXTPWCLS                          ! activates T'w at 2 m if .T
LOGICAL :: LXSOIL                            ! activates soil XFU if .T.
LOGICAL :: LTXTRD                            ! activates gravity wave drag momentum XFU at all levels if .T.
LOGICAL :: LTXTRC                            ! activates contribution of convection to U, V, q and (cp T) XFU if .T.
LOGICAL :: LTXTRT                            ! activates contribution of turbulence to U, V, q and (cp T) XFU if .T.
LOGICAL :: LTXR                              ! activates radiation XFU at all levels if .T.
LOGICAL :: LTXNEB                            ! activates cloudiness XFU at all levels if .T.
LOGICAL :: LTXQICE                           ! total ice water content at all levels
LOGICAL :: LTXQLI                            ! total liquid water content at all levels
LOGICAL :: LXICV                             ! activates indices of convection
                                             ! (CAPE and moisture convergence) XFU at all levels if .T.
LOGICAL :: LXCTOP                            ! activates pressure of top deep convection
LOGICAL :: LXCLP                             ! activates height (in meters) of PBL XFU at all levels if .T.
LOGICAL :: LXVEIN                            ! activates ventilation index
LOGICAL :: LXTGST                            ! activates gusts as U and V components XFU at all levels if .T.
LOGICAL :: LXXGST                            ! activates extreme gusts as U and V components XFU at all levels if .T.
LOGICAL :: LXXGST2                           ! activates extreme gusts2 as U and V components XFU at all levels if .T.
LOGICAL :: LXQCLS                            ! activates specific moisture at 2 meters
LOGICAL :: LXTHW                             ! activates "theta'_w" surface flux
LOGICAL :: LXXDIAGH                          ! activates extreme value of hail diagnostic
LOGICAL :: LXMRT                             ! activates mean radiant temperature
LOGICAL :: LXVISI                            ! activates visibilities diagnostic
LOGICAL :: LXVISI2                           ! activates visibilities diagnostic


TYPE(TXFUPTR) :: YXFUPT

REAL (KIND=JPRB), ALLOCATABLE :: XFUBUF(:,:,:)  ! Buffer for instantaneous diagnostics

END TYPE TXFU

!!TYPE(TXFU), POINTER :: YRXFU => NULL()

!     ------------------------------------------------------------
END MODULE YOMXFU
