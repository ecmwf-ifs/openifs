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

MODULE YOE_PHYS_MWAVE

USE PARKIND1  ,ONLY : JPRB, JPIM

IMPLICIT NONE


PRIVATE 

!       ----------------------------------------------------------------
!PHYSICS VARIABLES FOR MWAVE OPERATOR
!     D.Salmond      E.C.M.W.F.    23 Jan 2018
!     Moved out these arrays from GFL     
!       ----------------------------------------------------------------

TYPE,PUBLIC :: TEPHYSMWAVE
  LOGICAL :: LPHYS_MWAVE_FILLED_IN=.FALSE.
  REAL(KIND=JPRB),ALLOCATABLE :: PHYS_MWAVE(:,:,:,:)

  CONTAINS
  PROCEDURE :: CREATE
  PROCEDURE :: DESTROY
  PROCEDURE :: ZERO
  PROCEDURE :: COPY
END TYPE TEPHYSMWAVE

INTEGER(KIND=JPIM),PUBLIC :: N_PHYS_MWAVE

CONTAINS
!============================================================================
SUBROUTINE CREATE(SELF,YDGEOMETRY)
USE GEOMETRY_MOD , ONLY : GEOMETRY
IMPLICIT NONE
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
CLASS(TEPHYSMWAVE) :: SELF

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, NFLEVG=>YDDIMV%NFLEVG, &
 & NGPBLKS=>YDDIM%NGPBLKS)
N_PHYS_MWAVE=9 !THIS IS POOR CODING - NEEDS TO BE CLEANED UP.
               !THIS NUMBER MUST BE THE SAME AS IN ifs/module/par_gfl.F90
ALLOCATE(SELF%PHYS_MWAVE(NPROMA,NFLEVG,N_PHYS_MWAVE,NGPBLKS))

SELF%PHYS_MWAVE(:,:,:,:) = 0.0_JPRB

END ASSOCIATE
END ASSOCIATE

END SUBROUTINE CREATE
!============================================================================
SUBROUTINE DESTROY(SELF)
IMPLICIT NONE
CLASS(TEPHYSMWAVE) :: SELF

DEALLOCATE(SELF%PHYS_MWAVE)

END SUBROUTINE DESTROY
!============================================================================

SUBROUTINE ZERO(SELF)
CLASS(TEPHYSMWAVE) :: SELF
IF(ALLOCATED(SELF%PHYS_MWAVE)) SELF%PHYS_MWAVE=0.0_JPRB
END SUBROUTINE ZERO
!============================================================================
SUBROUTINE COPY(SELF,RHS)
CLASS(TEPHYSMWAVE) :: SELF
CLASS(TEPHYSMWAVE) :: RHS
!---------------------------------------------------------------------------
IF(ALLOCATED(SELF%PHYS_MWAVE).AND.ALLOCATED(RHS%PHYS_MWAVE)) THEN
  IF (ALL(SHAPE(SELF%PHYS_MWAVE) == SHAPE(RHS%PHYS_MWAVE))) THEN
    SELF%PHYS_MWAVE=RHS%PHYS_MWAVE
  ELSE
    CALL ABOR1("YOE_PHYS_MWAVE:COPY - DIFFERNT SHAPES")
  ENDIF
ELSE
  CALL ABOR1("YOE_PHYS_MWAVE:COPY - UNALLOCATED ARRAY")
ENDIF
!---------------------------------------------------------------------------
END SUBROUTINE COPY

END MODULE YOE_PHYS_MWAVE
