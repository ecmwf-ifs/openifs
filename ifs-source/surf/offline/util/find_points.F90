! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM FIND_POINTS
USE NETCDF
USE GET_POINTS_MOD
IMPLICIT NONE

INTEGER :: NP_TOT,NPOINTS
INTEGER     ,DIMENSION(:),ALLOCATABLE :: POINTS 
CHARACTER(LEN=500) :: INFO_FILE,NAM_FILE
INTEGER :: NP_LAND_OUT,NP_LAKE_OUT,NP_LAKE_FRAC_OUT,NP_OCEAN_OUT

NAMELIST /NAMFP/INFO_FILE

NAM_FILE='input.nam'

OPEN(99,FILE=TRIM(NAM_FILE),STATUS='old')
READ(99,NAMFP)
CLOSE(99)


CALL GET_POINTS(TRIM(INFO_FILE),TRIM(NAM_FILE),NP_TOT,NPOINTS,POINTS,&
              & NP_LAND_OUT,NP_LAKE_OUT,NP_LAKE_FRAC_OUT,NP_OCEAN_OUT)
	    
! PRINT*,'NP_TOT',NP_TOT
PRINT*,NPOINTS
! PRINT*,'NP_LAND',NP_LAND_OUT
! PRINT*,'NP_LAKE',NP_LAKE_OUT
! PRINT*,'NP_LAKE_FRAC',NP_LAKE_FRAC_OUT
! PRINT*,'NP_OCEAN',NP_OCEAN_OUT

END PROGRAM FIND_POINTS






