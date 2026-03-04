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

      FUNCTION ISRCHFGE(N,ARRAY,INC,TARGET)
      USE PARKIND1, ONLY : JPIM, JPRB
      IMPLICIT NONE
      INTEGER(KIND=JPIM) :: ISRCHFGE
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: INC
      REAL(KIND=JPRB) :: ARRAY(1+INC*(N-1))
      REAL(KIND=JPRB) :: TARGET
      INTEGER(KIND=JPIM) :: I
!
      IF (N.LE.0) THEN
        ISRCHFGE=0
        RETURN
      ELSE
        DO 1 I=1,N
          IF ( ARRAY(1+INC*(I-1)).GE.TARGET ) THEN
            ISRCHFGE=I
            RETURN
          ENDIF
 1      CONTINUE
        ISRCHFGE=N+1
      ENDIF
      ENDFUNCTION ISRCHFGE
