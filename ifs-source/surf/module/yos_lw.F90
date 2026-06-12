MODULE YOS_LW
 
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

!       ----------------------------------------------------------------
!*    ** *YOELW* - COEFFICIENTS OF THE LONGWAVE RADIATION TRANSFER
!       ----------------------------------------------------------------

TYPE :: TLW
INTEGER(KIND=JPIM) :: NSIL    ! NUMBER OF SPECTRAL INTERVALS
REAL(KIND=JPRB) :: TSTAND     ! REFERENCE TEMPERATURE FOR TEMPERATURE DEPENDENCE
REAL(KIND=JPRB) :: XP(6,6)    ! POLYNOMIAL COEFFICIENTS OF PLANCK FUNCTION
END TYPE TLW

END MODULE YOS_LW
