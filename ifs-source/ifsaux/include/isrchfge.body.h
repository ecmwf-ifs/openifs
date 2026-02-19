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

!define ISRCHFGE_RESULT
!define ISRCHFGE_N
!define ISRCHFGE_INC
!define ISRCHFGE_ARRAY(I)
!define ISRCHFGE_TARGET
      ISRCHFGE_RESULT=-1
      IF (ISRCHFGE_N.LE.0) THEN
        ISRCHFGE_RESULT=0
      ELSE
        DO ISRCHFGE_I=1,ISRCHFGE_N
          IF ( ISRCHFGE_ARRAY(1+ISRCHFGE_INC*(ISRCHFGE_I-1)).GE.ISRCHFGE_TARGET ) THEN
            ISRCHFGE_RESULT=ISRCHFGE_I
            EXIT
          ENDIF
        ENDDO
        IF (ISRCHFGE_RESULT.LT.0) THEN
        ISRCHFGE_RESULT=ISRCHFGE_N+1
        ENDIF
      ENDIF
