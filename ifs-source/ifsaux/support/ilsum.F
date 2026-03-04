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

      FUNCTION ILSUM(N,LX,INCX)
      USE PARKIND1, ONLY : JPIM
      IMPLICIT NONE
      INTEGER(KIND=JPIM) :: ILSUM
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: INCX
      INTEGER(KIND=JPIM) :: JL
      LOGICAL :: LX(*)
!
!  FINDS THE NUMBER OF TRUE ELEMENTS OF A LOGICAL VECTOR
!
      ILSUM=0
      IF ( N.LE.0) RETURN
!
      DO 1 JL=1,N
        IF ( LX(1+(JL-1)*INCX) ) ILSUM=ILSUM+1
 1    CONTINUE
!
      ENDFUNCTION ILSUM
