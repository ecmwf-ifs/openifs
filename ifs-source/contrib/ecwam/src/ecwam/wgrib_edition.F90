! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTEGER FUNCTION WGRIB_EDITION (IPARAMID)

! ----------------------------------------------------------------------

!**** *WGRIB_EDITION* - FUNCTION TO DETERMINE WHETHER A WAVE PARAMETER WITH IPARAMID
!                       CAN ONLY BE ENCODED IN GRIB 2 (i.e. WGRIB_EDITION=2).
!                       ALL OTHER PARAMETERS WGRIB_EDITION WILL RETURN 0
!                       AND THE DEFAULT GRIB EDITION WILL NEED TO BE USED.

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

! ----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JWIM), INTENT(IN) :: IPARAMID   ! grib parameter id (6 digit) 

! ----------------------------------------------------------------------

IF ( IPARAMID >= 140131 .AND. IPARAMID <= 140134 ) THEN
  WGRIB_EDITION = 2
ELSE
  WGRIB_EDITION = 0 
ENDIF

END FUNCTION WGRIB_EDITION
