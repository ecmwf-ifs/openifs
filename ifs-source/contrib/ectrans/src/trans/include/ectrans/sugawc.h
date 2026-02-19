! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE
SUBROUTINE SUGAWC(KDGLG,PMU,PW)

!**** *SUGAWC* - Compute Gaussian latitudes and weights 

!     Purpose.
!     --------
!     Compute Gaussian latitudes and weights.

!**   Interface.
!     ----------
!     CALL SUGAWC(...)

!     Explicit arguments :
!     -------------------- 
!      INPUT:
!       KDGLG    - number of latitudes.

!      OUTPUT:
!       PMU      - sine of Gaussian latitudes.
!       PW       - Gaussian weights.

!     Method.
!     -------

!     Externals.  SUGAW
!     ----------  

!     Author.
!     -------
!        K. Yessad, from SUGAWA and SULEG (trans)
!        Original : May 2012

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRD

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM) ,INTENT(IN)  :: KDGLG
REAL(KIND=JPRD)    ,INTENT(OUT) :: PMU(:)
REAL(KIND=JPRD)    ,INTENT(OUT) :: PW(:)

END SUBROUTINE SUGAWC

END INTERFACE
