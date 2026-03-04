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

MODULE DRYDEP_PAR_GC

!!****  * MODULE DRYDEP_PAR_GC
!!
!!     declare parameterization constants for Geos-Chem version of dry
!!     deposition
!!
!!    AUTHOR
!!    ------
!!    D. Finch (Univ Edinburgh)
!!
!!
!!    MODIFICATIONS
!!    -------------
!!    Original   Oct 2020a
!!
!!-------------------------------------------------------------------------------

USE  PARKIND1 ,ONLY : JPIM, JPRB

IMPLICIT NONE
SAVE

  INTEGER (KIND=JPIM), PARAMETER :: JPVEGNB_GC = 11_JPIM  ! number of Wesely vegetation types

  REAL (KIND=JPRB), PARAMETER :: ZLUSW = 50.0_JPRB
  REAL (KIND=JPRB), PARAMETER :: ZGSS_W = 50.0_JPRB

  REAL (KIND=JPRB), PARAMETER :: ZEPS = 1.E-5_JPRB
  REAL (KIND=JPRB), PARAMETER :: ZEPS1 = 1.E-4_JPRB
  REAL (KIND=JPRB), PARAMETER :: ZEPS2 = 10.0_JPRB

! DRY DEPOSITION TYPES:
!      GEOS-CHEM UPDATED
!  1)      SNOW/ICE
!  2)      DECIDUOUS FOREST
!  3)      CONIFEROUS FOREST
!  4)      AGRICULUTRAL LAND
!  5)      SHRUB/GRASSLAND
!  6)      AMAZON FOREST
!  7)      TUNDRA
!  8)      DESERT
!  9)      WETLAND
! 10)      URBAN
! 11)      WATER

! GEOS-Chem dry dep inputs 
!Surface resistance components for the "deposition land
!* types" are from Wesely [1989] except for tropical forests [Jacob and Wofsy,
!* 1990] and for tundra [Jacob et al., 1992].  All surface resistance
!* components are normalized to a leaf area index of unity.
REAL (KIND=JPRB), DIMENSION(JPVEGNB_GC) :: IRLU = &
& (/9999, 9000, 9000, 9000, 9000, 1000, 4000, 9999, 9000, 9999, 9999/)


! GEOS-Chem dry dep inputs : in-canopy resistance
REAL (KIND=JPRB), DIMENSION(JPVEGNB_GC) :: IRAC = &
& (/0, 2000, 2000, 200, 100, 2000, 0, 0, 300, 100, 0/)

! GEOS-Chem dry dep inputs : Ground, sulfur
REAL (KIND=JPRB), DIMENSION(JPVEGNB_GC) :: IRGSS = &
& (/100, 500, 500, 150, 350, 200, 340, 1000, 0, 400, 0/)

! GEOS-Chem dry dep inputs : ground, O3
REAL (KIND=JPRB), DIMENSION(JPVEGNB_GC) :: IRGSO = &
& (/9999, 200, 200, 150, 200, 200, 340, 400, 1000, 300, 2000/)

! GEOS-Chem dry dep inputs : lower canopy, sulfur
REAL (KIND=JPRB), DIMENSION(JPVEGNB_GC) :: IRCLS = &
& (/9999, 2000, 2000, 2000, 2000, 9999, 9999, 9999, 2500, 9999, 9999/)

! GEOS-Chem dry dep inputs : lower canopy for O3
REAL (KIND=JPRB), DIMENSION(JPVEGNB_GC) :: IRCLO = &
& (/1000, 1000, 1000, 1000, 1000, 9999, 9999, 9999, 1000, 9999, 9999/)

! GEOS-Cheme dry dep inputs: Stomatal
REAL (KIND=JPRB), DIMENSION(JPVEGNB_GC) :: IRI = &
& (/9999, 200, 400, 200, 200, 200, 200, 9999, 200, 9999, 9999/)

! Baldocchi polynomial coefficients for dry deposition
INTEGER (KIND=JPIM),PARAMETER  :: NPOLY = 20
REAL (KIND=JPRB), DIMENSION(NPOLY) :: COEFF = &
& (/-0.358, 3.02, 3.85, -0.0978, -3.66, 12., 0.252, -7.8, 0.226, 0.274, &
&    1.14, -2.19, 0.261, -4.62, 0.685, -0.254, 4.37, -0.266, -0.159,   &
&   -0.20/)

REAL(KIND=JPRB) , PARAMETER               :: PPRMAX = 1E+5 ! resistance maximum value
! used for surface, aerodynamic and laminar
! resistances

REAL (KIND=JPRB) , PARAMETER  :: PPDIMOH2O=0.15/0.595*1E-4  ! molecular diffusivity of H2O (m2/s)

END MODULE DRYDEP_PAR_GC
