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

!> geometry derived type for the IFS model
!> separated from geometry_mod to break dependency issues with sugeometry

MODULE TYPE_GEOMETRY

USE YOMVERT    , ONLY : TVAB, TVETA, TVFE, TVERTICAL_GEOM
USE YOMSTA     , ONLY : TSTA
USE YOMLAP     , ONLY : TLAP
USE YOMLEG     , ONLY : TCSGLEG
USE YOMDIM     , ONLY : TDIM
USE YOMDIMV    , ONLY : TDIMV
USE YOMMP      , ONLY : TMP
USE YOMGEM     , ONLY : TGEM
USE YOMCSGEOM  , ONLY : TCSGEOM
USE YOMGSGEOM  , ONLY : TGSGEOM, TGSGEOM_BLOCKED
USE YOMOROG    , ONLY : TOROG,   TOROG_BLOCKED
USE TYPE_SPGEOM, ONLY : TSPGEOM
#ifdef WITH_ATLAS
USE YOMATLAS   , ONLY : TATLAS
#endif
! LAM specific objects
USE YEMDIM     , ONLY : TEDIM
USE YEMGEO     , ONLY : TEGEO
USE YEMMP      , ONLY : TEMMP
USE YEMLAP     , ONLY : TLEP
USE YEMGSL     , ONLY : TEGSL
!
USE YEMLBC_GEO , ONLY : TELBC_GEO
!
USE YOMCVER , ONLY : TCVER


IMPLICIT NONE 

PRIVATE

TYPE, PUBLIC :: GEOMETRY
  LOGICAL,POINTER  :: LNONHYD_GEOM => NULL()
  TYPE(TVERTICAL_GEOM) :: YRVERT_GEOM

  TYPE(TVAB),POINTER     :: YRVAB => NULL()    ! Geometry restricted vertical hybrid coordinates 0->NFLEVG
  TYPE(TVETA),POINTER    :: YRVETA => NULL()
  TYPE(TVFE),POINTER     :: YRVFE => NULL()
  TYPE(TCVER),POINTER    :: YRCVER => NULL()
  TYPE(TSTA)    :: YRSTA
  TYPE(TLAP)    :: YRLAP
  TYPE(TCSGLEG) :: YRCSGLEG

  TYPE(TDIM)    :: YRDIM
  TYPE(TDIMV)   :: YRDIMV
  TYPE(TMP)     :: YRMP
  TYPE(TGEM)    :: YRGEM


  TYPE(TCSGEOM), POINTER :: YRCSGEOM(:) => NULL()
  TYPE(TCSGEOM)          :: YRCSGEOM_NB

  TYPE(TGSGEOM), POINTER :: YRGSGEOM(:) => NULL() 
  TYPE(TGSGEOM)          :: YRGSGEOM_NB
  TYPE(TGSGEOM_BLOCKED)  :: YRGSGEOM_B

  LOGICAL                    :: AD_GEOM = .FALSE.
  TYPE(TCSGEOM), POINTER     :: YRCSGEOMAD(:) => NULL()
  TYPE(TCSGEOM)              :: YRCSGEOMAD_NB
  TYPE(TGSGEOM), POINTER     :: YRGSGEOMAD(:) => NULL()
  TYPE(TGSGEOM)              :: YRGSGEOMAD_NB


  TYPE(TOROG), ALLOCATABLE :: YROROG(:)
  TYPE(TOROG_BLOCKED)      :: YROROG_B

  TYPE(TSPGEOM) :: YSPGEOM

#ifdef WITH_ATLAS
  TYPE(TATLAS),  POINTER :: YRATLAS => NULL()
#endif

  TYPE(TVAB)    :: YVABIO ! I/O full vertical hybrid coordinates 0->NIOLEVG (NIOLEVG >=NFLEVG)

  ! LAM specific objects
  TYPE(TEDIM) :: YREDIM
  TYPE(TEGEO) :: YREGEO
  TYPE(TEMMP) :: YREMP
  TYPE(TLEP)  :: YRELAP
  TYPE(TEGSL) :: YREGSL
!
  TYPE(TELBC_GEO), POINTER :: YRELBC_GEO => NULL()
!

CONTAINS

FINAL :: GEOMETRY_FINAL

END TYPE GEOMETRY

CONTAINS

SUBROUTINE GEOMETRY_FINAL(THIS)
  TYPE(GEOMETRY) :: THIS
  ! If we don't add this, we may get a internal compiler error
END SUBROUTINE GEOMETRY_FINAL

END MODULE TYPE_GEOMETRY
