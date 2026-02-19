! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE EC_FLUSH(UNIT)
  
!USE, INTRINSIC :: iso_fortran_env, ONLY: FLUSH
USE EC_PARKIND, ONLY : JPIM

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: UNIT

INTEGER(KIND=JPIM) :: IERR

IF (UNIT >= 0) THEN
   FLUSH(UNIT,IOSTAT=IERR,ERR=99) ! F2003 or later
ENDIF

99 CONTINUE

END SUBROUTINE EC_FLUSH
