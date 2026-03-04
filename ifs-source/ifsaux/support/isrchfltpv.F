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

      FUNCTION ISRCHFLTPV(N,ARRAY,INC,TARGET,NBITER)
      USE PARKIND1, ONLY : JPIM, JPRB
      IMPLICIT NONE
      INTEGER(KIND=JPIM) :: ISRCHFLTPV
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: INC
      REAL(KIND=JPRB) :: ARRAY(1+INC*(N-1))
      REAL(KIND=JPRB) :: TARGET
      INTEGER(KIND=JPIM) :: NBITER
      INTEGER   :: I
!
      IF (N.LE.0) THEN
        ISRCHFLTPV=0
        RETURN
      ELSE
        IF (NBITER.EQ.3) THEN
        DO 3 I=1,N-NBITER
          IF ( ARRAY(1+INC*(I-1)).LT.TARGET ) THEN
             IF ( ARRAY(1+INC*(I)).LT.TARGET ) THEN
                IF ( ARRAY(1+INC*(I+1)).LT.TARGET ) THEN
                   IF ( ARRAY(1+INC*(I+2)).LT.TARGET ) THEN
                      ISRCHFLTPV=I
                      RETURN
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
 3      CONTINUE
        DO I=N-NBITER+1,N
          IF ( ARRAY(1+INC*(I-1)).LT.TARGET ) THEN
            ISRCHFLTPV=I
            RETURN
          ENDIF
        ENDDO
        ELSEIF (NBITER.EQ.2) THEN
        DO 2 I=1,N-NBITER
          IF ( ARRAY(1+INC*(I-1)).LT.TARGET ) THEN
             IF ( ARRAY(1+INC*(I)).LT.TARGET ) THEN
                IF ( ARRAY(1+INC*(I+1)).LT.TARGET ) THEN
                      ISRCHFLTPV=I
                      RETURN
                ENDIF
             ENDIF
          ENDIF
 2      CONTINUE
        DO I=N-NBITER+1,N
          IF ( ARRAY(1+INC*(I-1)).LT.TARGET ) THEN
            ISRCHFLTPV=I
            RETURN
          ENDIF
        ENDDO
        ELSEIF (NBITER.EQ.1) THEN
        DO 1 I=1,N-NBITER
          IF ( ARRAY(1+INC*(I-1)).LT.TARGET ) THEN
             IF ( ARRAY(1+INC*(I)).LT.TARGET ) THEN
                      ISRCHFLTPV=I
                      RETURN
             ENDIF
          ENDIF
 1      CONTINUE
        DO I=N-NBITER+1,N
          IF ( ARRAY(1+INC*(I-1)).LT.TARGET ) THEN
            ISRCHFLTPV=I
            RETURN
          ENDIF
        ENDDO
        ELSEIF (NBITER.EQ.0) THEN
        DO 4 I=1,N-NBITER
          IF ( ARRAY(1+INC*(I-1)).LT.TARGET ) THEN
                      ISRCHFLTPV=I
                      RETURN
          ENDIF
 4      CONTINUE
        ELSE
        CALL ABOR1('ISRCHFLTPV: ERROR WITH NBITER VALUE')
        ENDIF
        ISRCHFLTPV=N+1
      ENDIF
      ENDFUNCTION ISRCHFLTPV
