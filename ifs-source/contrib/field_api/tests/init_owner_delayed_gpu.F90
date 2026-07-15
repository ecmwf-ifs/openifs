! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM INIT_OWNER_DELAYED_GPU
        !TEST IF DELAYED ALLOCATION IS WORKING ON GPU

        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE PARKIND1
        IMPLICIT NONE
        CLASS(FIELD_2RB), POINTER :: O => NULL()
        REAL(KIND=JPRB), POINTER :: PTR(:,:)
        LOGICAL :: OKAY
        INTEGER :: I,J

        CALL FIELD_NEW(O, LBOUNDS=[10,1], UBOUNDS=[21,11], DELAYED=.TRUE.)
        CALL O%GET_DEVICE_DATA_RDWR(PTR)
        !$ACC KERNELS PRESENT(PTR)
        PTR=42
        !$ACC END KERNELS

        OKAY=.TRUE.
        !$ACC SERIAL PRESENT(PTR) COPY(OKAY)
        IF(.NOT. ALL(PTR == 42))THEN
                OKAY=.FALSE.
        ENDIF
        !$ACC END SERIAL

        IF(OKAY .EQV. .FALSE.)THEN
                ERROR STOP
        ENDIF
        CALL FIELD_DELETE(O)
END PROGRAM INIT_OWNER_DELAYED_GPU
