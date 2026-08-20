! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWTEST

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *TESTO* - PRINTER OUTPUT UNIT AND TEST FLAGS.

      INTEGER(KIND=JWIM) :: IU06 = 6
      INTEGER(KIND=JWIM) :: ITEST 
      INTEGER(KIND=JWIM) :: ITESTB 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IU06*      INTEGER   UNIT FOR PRINTER OUTPUT.
!      *ITEST*     INTEGER   TEST OUTPUT LEVEL:
!                             .LE. 0  NO OUTPUT
!                             .GE. I  OUTPUT TO SUB. LEVEL I
!      *ITESTB*    INTEGER   MAX BLOCK NUMBER FOR OUTPUT IN BLOCK LOOPS

! ----------------------------------------------------------------------
      END MODULE YOWTEST
