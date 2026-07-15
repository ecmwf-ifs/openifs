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

INTEGER(KIND=JPIM), PARAMETER :: JPFUXT=250  ! maximum number of timesteps where XFU can be activated
INTEGER(KIND=JPIM), PARAMETER :: JPMXXFU=201 ! maximum number of XFU fields

TYPE :: TXFU_KEYS

LOGICAL :: LXFU=.FALSE.                      ! controls switch on/off all XFU
LOGICAL :: LXTRD=.FALSE.                     ! activates gravity wave drag momentum XFU if .T.
LOGICAL :: LXTRC=.FALSE.                     ! activates contribution of convection to U, V, q and (cp T) XFU if .T.
LOGICAL :: LXTRT=.FALSE.                     ! activates contribution of turbulence to U, V, q and (cp T) XFU if .T.
LOGICAL :: LXPLC=.FALSE.                     ! activates convective precipitation XFU if .T.
LOGICAL :: LXPLCG=.FALSE.                    ! activates convective graupels CFU if .T.
LOGICAL :: LXPLCH=.FALSE.                    ! activates convective hail CFU if .T.
LOGICAL :: LXPLS=.FALSE.                     ! activates stratiform precipitation XFU if .T.
LOGICAL :: LXPLSG=.FALSE.                    ! activates stratiform graupels CFU if .T.
LOGICAL :: LXPLSH=.FALSE.                    ! activates stratiform hail CFU if .T.
LOGICAL :: LXR=.FALSE.                       ! activates radiation XFU if .T.
LOGICAL :: LXNEBTT=.FALSE.                   ! activates total cloudiness XFU if .T.
LOGICAL :: LXNEBPA=.FALSE.                   ! activates partial cloudiness XFU if .T.
LOGICAL :: LXCLS=.FALSE.                     ! activates U, V, T, q and relative humidity at 2 or 10 m (time t-dt) if .T.
LOGICAL :: LXMWINDCLS=.FALSE.                ! activates mean of U and V at 10 m if .T., also NU/NV if LXNUVCLS
LOGICAL :: LXNUVCLS=.FALSE.                  ! activates neutral U and V at 10 m (time t-dt) if .T.
LOGICAL :: LXTTCLS=.FALSE.                   ! activates extreme temperatures at 2 m if .T.
LOGICAL :: LXHHCLS=.FALSE.                   ! activates extreme relative moistures at 2 m if .T
LOGICAL :: LXTPWCLS=.FALSE.                  ! activates T'w at 2 m if .T
LOGICAL :: LXSOIL=.FALSE.                    ! activates soil XFU if .T.
LOGICAL :: LTXTRD=.FALSE.                    ! activates gravity wave drag momentum XFU at all levels if .T.
LOGICAL :: LTXTRC=.FALSE.                    ! activates contribution of convection to U, V, q and (cp T) XFU if .T.
LOGICAL :: LTXTRT=.FALSE.                    ! activates contribution of turbulence to U, V, q and (cp T) XFU if .T.
LOGICAL :: LTXR=.FALSE.                      ! activates radiation XFU at all levels if .T.
LOGICAL :: LTXNEB=.FALSE.                    ! activates cloudiness XFU at all levels if .T.
LOGICAL :: LTXQICE=.FALSE.                   ! total ice water content at all levels
LOGICAL :: LTXQLI=.FALSE.                    ! total liquid water content at all levels
LOGICAL :: LXICV=.FALSE.                     ! activates indices of convection
                                             ! (CAPE and moisture convergence) XFU at all levels if .T.
LOGICAL :: LXCTOP=.FALSE.                    ! activates pressure of top deep convection
LOGICAL :: LXCLP=.FALSE.                     ! activates height (in meters) of PBL XFU at all levels if .T.
LOGICAL :: LXVEIN=.FALSE.                    ! activates ventilation index
LOGICAL :: LXTGST=.FALSE.                    ! activates gusts as U and V components XFU at all levels if .T.
LOGICAL :: LXXGST=.FALSE.                    ! activates extreme gusts as U and V components XFU at all levels if .T.
LOGICAL :: LXXGST2=.FALSE.                   ! activates extreme gusts2 as U and V components XFU at all levels if .T.
LOGICAL :: LXQCLS=.FALSE.                    ! activates specific moisture at 2 meters
LOGICAL :: LXTHW=.FALSE.                     ! activates "theta'_w" surface flux
LOGICAL :: LXXDIAGH=.FALSE.                  ! activates extreme value of hail diagnostic
LOGICAL :: LXMRT=.FALSE.                     ! activates mean radiant temperature
LOGICAL :: LXVISI=.FALSE.                    ! activates visibilities diagnostic
LOGICAL :: LXVISI2=.FALSE.                   ! activates visibilities diagnostic
LOGICAL :: LXPRECIPS1=.FALSE.                ! activates precipitations types nr 1 diagnostic
LOGICAL :: LXPRECIPS2=.FALSE.                ! activates precipitations types nr 2 diagnostic

END TYPE TXFU_KEYS

TYPE, EXTENDS(TXFU_KEYS) :: TXFU

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

LOGICAL :: LRESET                            ! reset extreme temperatures to zero
LOGICAL :: LRESET_GST                        ! reset Gust calculation
LOGICAL :: LRESET_GST2                       ! reset Gust2 calculation
LOGICAL :: LRESET_PRECIP                     ! reset Precips type calcultation
LOGICAL :: LRESET_PRECIP2                    ! reset Precips type calcultation
LOGICAL :: LRESET_VISI                       ! reset visibilities calculations
LOGICAL :: LRESET_VISI2                      ! reset visibilities calculations

LOGICAL :: LREAXFU                           ! read first input on historic file if .T.

TYPE(TXFUPTR) :: YXFUPT

REAL (KIND=JPRB), ALLOCATABLE :: XFUBUF(:,:,:)  ! Buffer for instantaneous diagnostics

END TYPE TXFU

!     ------------------------------------------------------------
END MODULE YOMXFU
