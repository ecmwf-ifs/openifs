! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE EC_LUN

USE EC_PARKIND  ,ONLY : JPIM

IMPLICIT NONE

SAVE

PRIVATE :: JPIM

PUBLIC
!     ------------------------------------------------------------------

!*    Logical units used by code

!     NULOUT :   output unit
!     NULERR :   unit number for comparison with reference run

INTEGER(KIND=JPIM) :: NULOUT = 6
INTEGER(KIND=JPIM) :: NULERR = 0

!     ------------------------------------------------------------------
END MODULE EC_LUN
