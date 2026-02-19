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

MODULE BASCOE_LBC_MODULE


USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

  ! Boundary conditions at surface for stratospheric species
  ! Which currently don't feature emissions (!)
  ! array of monthly varying, latitude band dependent, values
  ! allocated and loaded loaded from external file
  ! in routine BASCOE_LBC_INI
  INTEGER(KIND=JPIM)    :: NMONTH_LBC, NLATBOUND_LBC
  INTEGER(KIND=JPIM) , DIMENSION(:), ALLOCATABLE    :: MONTH_LBC
  REAL(KIND=JPRB) , DIMENSION(:), ALLOCATABLE       :: XLATBOUND_LBC
  REAL(KIND=JPRB) , DIMENSION(:,:,:), ALLOCATABLE   :: VALUES_LBC


END MODULE BASCOE_LBC_MODULE
