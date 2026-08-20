! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE INIWCST(PRPLRADI) 

! ----------------------------------------------------------------------

!**** *INIWCST* -

!*    PURPOSE.
!     --------
!     INITIALISES CONSTANTS

!      *RPLRADI*   REAL      EARTH RADIUS REDUCTION FACTOR FOR SMALL PLANET

!**   INTERFACE.
!     ----------
!     *CALL INIWCST*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : PI       ,CIRC     ,ZPI      ,ZCONST   ,    &
     &            RAD      ,DEG      ,R        ,ZPISQRT  ,ZPI4GM1  ,    &
     &            ZPI4GM2  ,G        ,THREEZPI

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: PRPLRADI

      PI = 4.0_JWRB*ATAN(1.0_JWRB)
      ZPI= 2.0_JWRB*PI
      THREEZPI = 3.0_JWRB*ZPI
      ZPISQRT = SQRT(ZPI)
      ZPI4GM1 = ZPI**4/G
      ZPI4GM2 = ZPI**4/G**2

      ZCONST=1.0_JWRB/(8.0_JWRB*PI*SQRT(2.0_JWRB))
      RAD = PI/180.0_JWRB
      DEG = 180./PI
      R = CIRC/ZPI*PRPLRADI

      END SUBROUTINE INIWCST
