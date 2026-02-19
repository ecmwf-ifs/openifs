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

MODULE SUPPDATE_MOD

USE PARKIND1, ONLY : JPIM

IMPLICIT NONE

INTEGER (KIND=JPIM), PARAMETER :: IDATE0 (22) = &
& (/ 0_JPIM,  0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM,  0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, &
&    0_JPIM,  0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM,  0_JPIM, 0_JPIM, 0_JPIM, 0_JPIM /)

INTEGER (KIND=JPIM), TARGET, DIMENSION (22) :: IDATEF_CFNISH = IDATE0

SAVE

END MODULE SUPPDATE_MOD
