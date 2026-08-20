! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM WRAPPER_MODIFY_GPU
        ! TEST IF DATA ARE CORRECTLY TRANSFERED ON GPU
        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE PARKIND1
        IMPLICIT NONE
        CLASS(FIELD_2RB), POINTER :: W => NULL()
        REAL(KIND=JPRB), ALLOCATABLE :: D(:,:)
        REAL(KIND=JPRB), POINTER :: D_GPU(:,:), D_CPU(:,:)
        LOGICAL :: RES
        INTEGER :: I,J

        ALLOCATE(D(10,10))
        D=7

        CALL FIELD_NEW(W, DATA=D)
        CALL W%GET_DEVICE_DATA_RDWR(D_GPU)

        !$ACC SERIAL PRESENT(D_GPU) COPY(RES)
        RES=.TRUE.
        DO I=1,10
        DO J=1,10
        D_GPU(J,I)=I*J
        END DO
        END DO
        !$ACC END SERIAL


        ! DATA SHOULD ALL BE 7 BECAUSE THEY HAVEN'T MOVED BACK FROM GPU
        ! YET. BUT IF OPENACC WAS DESACTIVATED THEN THE VALUES HAVE
        ! ALREADY CHANGED SINCE IT ALWAYS HAPPEN ON THE CPU MAIN MEMORY
        IF ( .NOT. ALL(D == 7))THEN
                ERROR STOP
        END IF

        CALL W%GET_HOST_DATA_RDONLY(D_CPU)

        DO I=1,10
        DO J=1,10
        IF (D_CPU(J,I)/=I*J) THEN
                ERROR STOP
        END IF
        END DO
        END DO
        CALL FIELD_DELETE(W)
END PROGRAM WRAPPER_MODIFY_GPU
