! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM INIT_OWNER_GPU
        ! TEST IF DATA ARE CORRECTLY TRANSFERED ON GPU

        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE PARKIND1
        IMPLICIT NONE
        CLASS(FIELD_2RB), POINTER :: O => NULL()
        REAL(KIND=JPRB), POINTER :: PTR_CPU(:,:)
        REAL(KIND=JPRB), POINTER :: PTR_GPU(:,:)
        LOGICAL :: OKAY
        INTEGER :: I,J

        CALL FIELD_NEW(O, LBOUNDS=[10,1], UBOUNDS=[21,11], PERSISTENT=.TRUE.)
        CALL O%GET_HOST_DATA_RDWR(PTR_CPU)
        PTR_CPU=42

        CALL O%GET_DEVICE_DATA_RDONLY(PTR_GPU)
        OKAY=.TRUE.
        !$ACC SERIAL PRESENT(PTR_GPU) COPY(OKAY)
        DO I=10,21
        DO J=1,11
        IF(PTR_GPU(I,J) /= 42) THEN
                OKAY = .FALSE.
        END IF
        END DO
        END DO
        !$ACC END SERIAL

        IF (.NOT. OKAY) THEN
                ERROR STOP
        END IF
        CALL FIELD_DELETE(O)
END PROGRAM INIT_OWNER_GPU
