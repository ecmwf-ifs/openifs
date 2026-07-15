! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE YOWSPHERE

! ------------------------------------------------------------------------------
IMPLICIT NONE

SAVE
PRIVATE

PUBLIC :: SPHERICAL_COORDINATE_DISTANCE
PUBLIC :: SPHERICAL_COORDINATE_AREA

INTERFACE SPHERICAL_COORDINATE_DISTANCE
  MODULE PROCEDURE SPHERICAL_COORDINATE_DISTANCE
END INTERFACE

INTERFACE SPHERICAL_COORDINATE_AREA
  MODULE PROCEDURE SPHERICAL_COORDINATE_AREA
END INTERFACE

CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
SUBROUTINE SPHERICAL_COORDINATE_DISTANCE(LON1, LON2, LAT1, LAT2, DIST)
!     Purpose.: computes the distance on a sphere of radius=1 between (LON1,LAT1) and (LON2,LAT2)
!     --------

!        Explicit arguments :  
!        --------------------   

!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!     http://en.wikipedia.org/wiki/Great-circle_distance

!     Author.
!     -------

!     Modifications.
!     --------------
!     --------------------------------------------------------------
USE YOWPCONS , ONLY : RAD
USE PARKIND_WAVE, ONLY : JWRU

IMPLICIT NONE

REAL(KIND=JWRU), INTENT(IN) :: LON1, LON2, LAT1, LAT2
REAL(KIND=JWRU), INTENT(OUT) :: DIST
REAL(KIND=JWRU) :: SLAT, SLON, C1, C2

SLAT = SIN(0.5_JWRU*(LAT1-LAT2)*RAD)**2
SLON = SIN(0.5_JWRU*(LON1-LON2)*RAD)**2
C1=COS(LAT1*RAD)
C2=COS(LAT2*RAD)
DIST = SQRT(MAX(SLAT+C1*C2*SLON,0.0_JWRU))

IF (DIST .ge. 1.0_JWRU) THEN
  DIST=0.0_JWRU
ELSE
  DIST = 2.0_JWRU*ASIN(DIST) 
END IF
END SUBROUTINE SPHERICAL_COORDINATE_DISTANCE
!**********************************************************************
!*                                                                    *
!**********************************************************************
SUBROUTINE SPHERICAL_COORDINATE_AREA(LON1, LON2, LON3, LAT1, LAT2, LAT3, AREA)
USE PARKIND_WAVE, ONLY : JWRU
IMPLICIT NONE
REAL(KIND=JWRU), INTENT(IN) :: LON1, LON2, LON3, LAT1, LAT2, LAT3
REAL(KIND=JWRU), INTENT(OUT) :: AREA
REAL(KIND=JWRU) :: DistA, DistB, DistC, DistS
REAL(KIND=JWRU) :: eTan1, eTan2, eTan3, eTan4
REAL(KIND=JWRU) :: eProd, sqrtProd
CALL SPHERICAL_COORDINATE_DISTANCE(LON1, LON2, LAT1, LAT2, DistA)
CALL SPHERICAL_COORDINATE_DISTANCE(LON1, LON3, LAT1, LAT3, DistB)
CALL SPHERICAL_COORDINATE_DISTANCE(LON2, LON3, LAT2, LAT3, DistC)
DistS=0.5_JWRU*(DistA + DistB + DistC)
eTan1=tan(0.5_JWRU*DistS)
eTan2=tan(0.5_JWRU*(DistS - DistA))
eTan3=tan(0.5_JWRU*(DistS - DistB))
eTan4=tan(0.5_JWRU*(DistS - DistC))
eProd=eTan1*eTan2*eTan3*eTan4
sqrtProd=SQRT(eProd)
AREA=4.0_JWRU*ATAN(sqrtProd)
END SUBROUTINE SPHERICAL_COORDINATE_AREA
!**********************************************************************
!*                                                                    *
!**********************************************************************
END MODULE YOWSPHERE 
