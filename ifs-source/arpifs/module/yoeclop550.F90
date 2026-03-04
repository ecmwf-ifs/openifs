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

MODULE YOECLOP550

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!     * YOECLOP550* PARAMETERS FOR CLOUD OPTICAL PROPERTIES
!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: RSA55, RSB55,  RSC55, RSD55

REAL(KIND=JPRB) :: RFA55(2)

!*    SLINGO (1989) WATER CLOUD OPTICAL PROPERTIES

! RSA :  REAL   : C1 IN OPTICAL THICKNESS FORMULA
! RSB :  REAL   : C2 IN OPTICAL THICKNESS FORMULA

!*    ICE CLOUD OPTICAL PROPERTIES DERIVED FROM FU (1996)

! RFA :  REAL   : COEFFICIENTS IN EXPRESSION FOR SW EXTINCTION COEFF.

!     -----------------------------------------------------------------
END MODULE YOECLOP550
