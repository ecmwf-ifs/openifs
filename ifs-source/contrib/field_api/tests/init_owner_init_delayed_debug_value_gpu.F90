! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM INIT_OWNER_INIT_DELAYED_DEBUG_VALUE_GPU
        !TEST THAT ALL DATA ARE INITIALIZED TO THE DEBUG VALUE CHOOSED BY THE
        !USER WHEN COMBINED WITH DELAYED ALLOCATION ON GPU
        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE FIELD_INIT_DEBUG_VALUE_MODULE
        USE PARKIND1
        IMPLICIT NONE
        CLASS(FIELD_2IM), POINTER :: O => NULL()
        INTEGER(KIND=JPIM), POINTER :: PTR(:,:)
        LOGICAL :: OKAY
        integer :: i,j

        USE_INIT_DEBUG_VALUE = .TRUE.
        INIT_DEBUG_VALUE_JPIM=-123456789

        CALL FIELD_NEW(O, LBOUNDS=[1,1], UBOUNDS=[10,10], PERSISTENT=.TRUE., DELAYED=.TRUE.)
        CALL O%GET_DEVICE_DATA_RDONLY(PTR)

        OKAY=.TRUE.
        !$ACC PARALLEL LOOP PRESENT(PTR) COPY(OKAY)
        DO I=1,10
        DO J=1,10
        IF(PTR(I,J) /= -123456789)THEN
                OKAY=.FALSE.
        ENDIF
        END DO
        END DO
        !$ACC END PARALLEL

        IF(OKAY .EQV. .FALSE.)THEN
                ERROR STOP
        ENDIF

        CALL FIELD_DELETE(O)
END PROGRAM INIT_OWNER_INIT_DELAYED_DEBUG_VALUE_GPU
