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

MODULE YEMGEO

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

TYPE :: TEGEO

!     ------------------------------------------------------------------

!*    Defining specific geometry variables for LAM model.
!     These variables are set-up in SUEGEM_NAML and should not be modified elsewhere.

!     IMPORTANT REMARK ABOUT UNITS OF ANGLES (ELON.., ELAT.. RLON_ACAD, RLAT_ACAD):
!      - Namelist values are in degrees.
!      - Printings done in SUEGEM_NAML are in degrees.
!      - Final values are computed in SUEGEM_NAML.
!      - When leaving SUEGEM_NAML, values SHOULD be in radians.
!      - Other routines using these values assume that they are in radians.

!   LMAP :
!        .T. : COORDINATES OF THE DOMAIN IN NAMELIST
!        .F. : WAVELENGHTS OF THE DOMAIN IN NAMELIST
!   ERPK                     : K projection parameter and definition
!         ERPK=10   -- projection type self determined by minimizing
!                      the variation of the map factor
!         ERPK=1    -- polar stereographic projection
!   0 < ERPK < 1    -- lambert conformal projection with the cone
!                      parameter ERPK
!         ERPK=0    -- mercator conformal projection
!         ERPK<0    -- no projection
!   ELON1                    : LA1G geographic longitude of the south-west
!                              corner of useful domain
!   ELAT1                    : FI1G geographic latitude of the south-west
!                              corner of useful domain
!   ELON2                    : LA1G geographic longitude of the north-east
!                              corner of useful domain
!   ELAT2                    : FI1G geographic latitude of the north-east
!                              corner of useful domain
!   ELONC                    : geographic longitude of the centre of domain
!   ELATC                    : geographic latitude of the centre of domain
!   ELON0                    : LA0 geographic longitude of reference for
!                              the projection
!                              (it is not required when NGIV0=1 or 3)
!   ELAT0                    : FI0 geographic latitude of reference for
!                              the projection
!                              (it is not required when NGIV0=2 or 3)
!                                     unit on y axis
!   EDELX                    : grid size in m along x
!   EDELY                    : grid size in m along y
!   ELX                      : wavelenght of the domain in x direction
!   ELY                      : wavelenght of the domain in y direction
!   EXWN                     : wavenumber in x direction
!   EYWN                     : wavenumber in y direction

!   RCORI_ACAD               : Coriolis parameter for academic experiments
!                              (LMAP=.F.)
!   RLAT_ACAD                : Latitude for academic experiments
!                              (LMAP=.F.)
! if both RCORI_ACAD and RLAT_ACAD are specified in namelist,
! RLAT_ACAD is used for radiation only and RCORI_ACAD is set to 0.

!   RLON_ACAD                : Longitude for academic experiments
!                              (LMAP=.F.), for radiation only

!   LMRT                     : key for using Mercator Rotated Tilted 
!                              with new eggx cadre format
!                  .T.   Mercator Rotated Tilted
!                  .F.   other projection

!   LREDEL_IN_METRES         : T: EDELX, EDELY, ELX, ELY are in metres.
!                              F: EDELX, EDELY, ELX, ELY are in radians.

REAL(KIND=JPRB) :: ERPK
REAL(KIND=JPRB) :: ELON1
REAL(KIND=JPRB) :: ELAT1
REAL(KIND=JPRB) :: ELON2
REAL(KIND=JPRB) :: ELAT2
REAL(KIND=JPRB) :: ELON0
REAL(KIND=JPRB) :: ELAT0
REAL(KIND=JPRB) :: ELONC
REAL(KIND=JPRB) :: ELATC
REAL(KIND=JPRB) :: EDELX
REAL(KIND=JPRB) :: EDELY
REAL(KIND=JPRB) :: ELX
REAL(KIND=JPRB) :: ELY
REAL(KIND=JPRB) :: EXWN
REAL(KIND=JPRB) :: EYWN
REAL(KIND=JPRB) :: RCORI_ACAD
REAL(KIND=JPRB) :: RLAT_ACAD
REAL(KIND=JPRB) :: RLON_ACAD
LOGICAL :: LMAP
LOGICAL :: LMRT
LOGICAL :: LREDEL_IN_METRES

END TYPE TEGEO

!     ------------------------------------------------------------------
END MODULE YEMGEO
