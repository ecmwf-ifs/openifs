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

MODULE YOMSENS

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     --------------------------------------------------

!*    Control for SENSitivity job

!     LGRVOL = .T.  : ACTIVATE GRADIENT NORMALIZATION BY SIDELP TO WRITE
!                        FILES OF 3D DENSITY OF GRADIENT
!     NJROPT        : TYPE OF COST FUNCTION:
!                             1 - QUADRATIC DISTANCE TO A REFERENCE (DEFAULT)
!                             2 - LINEAR INTEGRAL OF PARAMETERS
!     LBSENS        : USE B-matrix in the initial norm

LOGICAL :: LGRVOL
INTEGER(KIND=JPIM) :: NJROPT
LOGICAL :: LBSENS

!     --------------------------------------------------
END MODULE YOMSENS
