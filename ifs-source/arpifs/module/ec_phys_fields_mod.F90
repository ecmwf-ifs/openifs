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

MODULE EC_PHYS_FIELDS_MOD 

USE YOE_TILE_PROP, ONLY : TETILEPROP
USE YOE_PHYS_MWAVE, ONLY : TEPHYSMWAVE
IMPLICIT NONE
PRIVATE

! Fields related to ECMWF physics package

TYPE,PUBLIC  :: TEC_PHYS_FIELDS
TYPE(TETILEPROP) :: YRTILEPROP
TYPE(TEPHYSMWAVE) :: YRPHYSMWAVE
CONTAINS
PROCEDURE :: CREATE
PROCEDURE :: ZERO
PROCEDURE :: COPY
END TYPE TEC_PHYS_FIELDS
CONTAINS
!============================================================================
SUBROUTINE CREATE(SELF,YDGEOMETRY,YDDPHY)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDPHY      , ONLY : TDPHY
CLASS(TEC_PHYS_FIELDS) :: SELF
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDPHY)       ,INTENT(INOUT):: YDDPHY

CALL SELF%YRTILEPROP%CREATE(YDGEOMETRY,YDDPHY)
CALL SELF%YRPHYSMWAVE%CREATE(YDGEOMETRY)

END SUBROUTINE CREATE
!============================================================================
SUBROUTINE ZERO(SELF)
CLASS(TEC_PHYS_FIELDS) :: SELF

CALL SELF%YRTILEPROP%ZERO()
CALL SELF%YRPHYSMWAVE%ZERO()

END SUBROUTINE ZERO
!============================================================================

SUBROUTINE COPY(SELF,RHS)
CLASS(TEC_PHYS_FIELDS) :: SELF
CLASS(TEC_PHYS_FIELDS) :: RHS

CALL SELF%YRTILEPROP%COPY(RHS%YRTILEPROP)
CALL SELF%YRPHYSMWAVE%COPY(RHS%YRPHYSMWAVE)

END SUBROUTINE COPY
!============================================================================

END MODULE EC_PHYS_FIELDS_MOD

