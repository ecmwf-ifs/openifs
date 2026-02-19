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

MODULE YOMCPL

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!     -------------------------------------------------------
!**** *YOMCPL * - MODULE INCLUDING THE COUPLED FIELDS
!                  (OPA or SLAB)
!     AUTHOR.
!     -------
!        /CNRM/GMGEC/EAC/JPh Piedelievre 2000-11-30

!     MODIFICATIONS.
!     --------------
!        03-04-08 : /CNRM/GMGEC/EAC/JPh Piedelievre
!        27-Apr-2009, E.Maisonnave: Add 10m wind module and sublimation
!     -------------------------------------------------------

!*    Flux de couplage

! FCSTSU   STRESS SURFACE DE U                ; OPA, SLAB
! FCSTSV   STRESS SURFACE DE V                ; OPA, SLAB
! FCCHAS   FLUX CHALEUR SURFACE               ; OPA, SLAB
! FCRSOS   FLUX RAYT SOLAIRE SURFACE          ; OPA, SLAB
! FCHUMS   FLUX HUMIDITE SURFACE              ; OPA
! FCRUIS   RUISSELLEMENT TOTAL                ; OPA
! FCCHLL   FLUX CHALEUR LATENTE EAU LIQUIDE   ;      SLAB
! FCCHLN   FLUX CHALEUR LATENTE NEIGE         ;      SLAB
! FCCHSS   FLUX CHALEUR SENSIBLE              ; OPA, SLAB
! FCHUML   FLUX PRECIPITATIONS EAU LIQUIDE    ; OPA
! FCHUMN   FLUX PRECIPITATIONS NEIGE          ; OPA
! FCTSUR2  TEMPERATURE DE SURFACE             ; OPA
! FCTSTS   TSURF2*TSUR2                       ; OPA
! FCTSFL   TSUR2*FCHA                         ; OPA
! FCWIMO   MODULE DU VENT A 10M               ; OPA
! FCSUBL   SUBLIMATION NEIGE                  ; OPA


REAL(KIND=JPRB),ALLOCATABLE:: FCSTSU(:)
REAL(KIND=JPRB),ALLOCATABLE:: FCSTSV(:)
REAL(KIND=JPRB),ALLOCATABLE:: FCCHAS(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: FCRSOS(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: FCHUMS(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: FCRUIS(:)
REAL(KIND=JPRB),ALLOCATABLE:: FCCHLL(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: FCCHLN(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: FCCHSS(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: FCHUML(:)
REAL(KIND=JPRB),ALLOCATABLE:: FCHUMN(:)
REAL(KIND=JPRB),ALLOCATABLE:: FCTSUR2(:)
REAL(KIND=JPRB),ALLOCATABLE:: FCTSTS(:)
REAL(KIND=JPRB),ALLOCATABLE:: FCTSFL(:)
REAL(KIND=JPRB),ALLOCATABLE:: FCWIMO(:)
REAL(KIND=JPRB),ALLOCATABLE:: FCSUBL(:)

!     -----------------------------------------------------------------
END MODULE YOMCPL
