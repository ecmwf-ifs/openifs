! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

subroutine get_openmp(kopenmp)
use EC_PARKIND, only : JPIM
implicit none
INTEGER(KIND=JPIM), INTENT(out) :: kopenmp
#ifdef _OPENMP
kopenmp = _OPENMP
#else
kopenmp = 0
#endif
end subroutine get_openmp
