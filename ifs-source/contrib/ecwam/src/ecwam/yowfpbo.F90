! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWFPBO

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *FPBOUN* USED FOR THE FINE GRID
!                       ORGANIZATION THE BOUNDARY POINTS.

      INTEGER(KIND=JWIM)              :: IBOUNF
      INTEGER(KIND=JWIM)              :: NBOUNF
      INTEGER(KIND=JWIM), ALLOCATABLE :: IJARF(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: IGARF(:) !!!! obsolete
      INTEGER(KIND=JWIM), ALLOCATABLE :: IBFL(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: IBFR(:)

      REAL(KIND=JWRB), ALLOCATABLE    :: BFW(:)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IBOUNF*    INTEGER   FLAG FOR THE FINE GRID INPUT OF BOUNDARY VALUES.
!                            = 1; THE RUN INCLUDES INPUT BOUNDARY POINTS.
!                            ELSE; NO BOUNDARY POINTS.
!      *NBOUNF*    INTEGER   NUMBER OF FINE GRID BOUNDARY POINTS.
!      *IGARF*     INTEGER   BLOCK INDEX FOR A FINE GRID BOUNDARY POINT.
!      *IJARF*     INTEGER   POINT INDEX FOR A FINE GRID BOUNDARY POINT.
!      *IBFL*      INTEGER   INDEX OF LEFT COARSE GRID OUTPUT POINT.
!      *IBFR*      INTEGER   INDEX OF RIGHT COARSE GRID OUTPUT POINT.
!      *BFW*       REAL      SPACE INTERPOLATION WEIGHT.

! ----------------------------------------------------------------------
      END MODULE YOWFPBO
