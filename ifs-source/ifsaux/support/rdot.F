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

      FUNCTION RDOT(N,X,INCX,Y,INCY)
      USE PARKIND1, ONLY : JPIM, JPRB
      IMPLICIT NONE
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: INCX
      INTEGER(KIND=JPIM) :: INCY
      REAL(KIND=JPRB) :: X(N*INCX)
      REAL(KIND=JPRB) :: Y(N*INCY)
      REAL(KIND=JPRB) :: RDOT
      INTEGER :: IX,IY,J
      REAL(KIND=JPRB) :: ZDOT
      ZDOT=0.0_JPRB
      IX=1
      IY=1
      DO J=1,N
        ZDOT=ZDOT+X(IX)*Y(IY)
        IX=IX+INCX
        IY=IY+INCY
      ENDDO
      RDOT=ZDOT
      RETURN
      ENDFUNCTION RDOT
