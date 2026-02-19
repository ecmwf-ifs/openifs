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

      FUNCTION RSUM(N,X,INCX)
      USE PARKIND1, ONLY : JPIM, JPRB
      IMPLICIT NONE
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: INCX
      REAL(KIND=JPRB) :: X(N*INCX)
      REAL(KIND=JPRB) :: RSUM
      INTEGER(KIND=JPIM) :: I, J
      REAL(KIND=JPRB) :: ZSUM
!
      ZSUM=0.0_JPRB
      I=1
      DO J=1,N
        ZSUM=ZSUM+X(I)
        I=I+INCX
      ENDDO
      RSUM=ZSUM
      RETURN
      ENDFUNCTION RSUM
