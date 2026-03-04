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

MODULE YOMSLINT

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMVSPLIP  , ONLY : TVSPLIP
USE YOMVSLETA  , ONLY : TVSLETA
USE YOMHSLMER  , ONLY : THSLMER

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!     USED FOR BIT REPRODUCIBILITY OF SEMI-LAGRANGIAN ADJOINT (SLAD)

!     NADMAP     : USED FOR REMAPPING POINTS IN SL BUFFER SO THAT THEY
!                : ARE UNIQUE
!     NADCORE    : INNER HALO FOR ADJOINT INTERPOLATIONS
!                : ARE UNIQUE
!     LADCORE    : T - IF ONE OF THE NGPTOTAD CORE POINTS, F OTHERWISE
!     NGPTOTAD   : NUMBER OF GRID POINTS IN INNER HALO
!     RSASIGN    : USED FOR CORRECTING SYMMETRIC / ANTISYMMETRIC PROPERTIES
!                : (MIRRORED LATITUDES NEAR POLES)

TYPE :: TSLINT
  TYPE(TVSPLIP) :: YRVSPLIP
  TYPE(TVSLETA) :: YRVSLETA
  TYPE(THSLMER) :: YRHSLMER

!     ------------------------------------------------------------------
CONTAINS
  
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 

END TYPE TSLINT
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------


!!TYPE(TSLINT), POINTER :: YRSLINT => NULL()

CONTAINS 
  
SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TSLINT), INTENT(IN) :: SELF
  INTEGER      , INTENT(IN) :: KDEPTH
  INTEGER      , INTENT(IN) :: KOUTNO

  INTEGER :: IDEPTHLOC

  IDEPTHLOC = KDEPTH + 2
  
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_dyn%yrslint : '
  
END SUBROUTINE PRINT_CONFIGURATION
!     ------------------------------------------------------------------
END MODULE YOMSLINT
