! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE NON_CONTI
        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE PARKIND1
        IMPLICIT NONE
CONTAINS
        SUBROUTINE DO_STUFF_WITH_NONCONTIGUOUS_DATA(D)
                REAL(KIND=JPRB) :: D(:,:)
                CLASS(FIELD_2RB), POINTER :: W => NULL()
                REAL(KIND=JPRB), POINTER :: PTR(:,:)
                CALL FIELD_NEW(W, DATA=D)
                CALL W%GET_HOST_DATA_RDWR(PTR)
                PTR=42
                CALL FIELD_DELETE(W)
        END SUBROUTINE DO_STUFF_WITH_NONCONTIGUOUS_DATA
END MODULE NON_CONTI
PROGRAM INIT_WRAPPER_NON_CONTIGUOUS
        ! TEST IF WRAPPER WORKS WITH NON CONTIGUOUS ARRAY
        USE FIELD_MODULE
        USE PARKIND1
        USE NON_CONTI

        IMPLICIT NONE
        REAL(KIND=JPRB), ALLOCATABLE :: D(:,:)
        INTEGER :: I, J

        ALLOCATE(D(10,10))
        D=7
        CALL DO_STUFF_WITH_NONCONTIGUOUS_DATA(D(5:10,1:2))
        DO I=5,10
        DO J=1,2
        IF (D(I,J) /= 42) THEN
                ERROR STOP
        END IF
        END DO
        END DO

END PROGRAM INIT_WRAPPER_NON_CONTIGUOUS
