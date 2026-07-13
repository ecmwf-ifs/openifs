MODULE YOS_DIM
 
USE PARKIND1  ,ONLY : JPIM     ,JPRB

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

IMPLICIT NONE

SAVE

INTEGER(KIND=JPIM), PARAMETER :: JPSOILTY=11  ! number of soil types (CH)
INTEGER(KIND=JPIM), PARAMETER :: JPTEXT=3     ! number of macro-types (TESSEL->2only)

TYPE :: TDIM
INTEGER(KIND=JPIM) :: NCSS    ! number of levels in SOILB
INTEGER(KIND=JPIM) :: NTILES  ! number of surface tiles
INTEGER(KIND=JPIM) :: NMONTH  ! number of months im a year
INTEGER(KIND=JPIM) :: NCSNEC  ! number of snow levels 
END TYPE TDIM

END MODULE YOS_DIM
