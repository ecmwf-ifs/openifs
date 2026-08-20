! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM CPU_TO_GPU
        !TEST INITIALISING DATA ON CPU THEN MOVING THEM ON GPU

        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE PARKIND1
        IMPLICIT NONE
        CLASS(FIELD_2IM), POINTER :: O => NULL()
        INTEGER(KIND=JPIM), POINTER :: PTR_HOST(:,:) => NULL()
        INTEGER(KIND=JPIM), POINTER :: PTR_DEV(:,:) => NULL()
        INTEGER ::I,J
        LOGICAL :: OKAY

        CALL FIELD_NEW(O, LBOUNDS=[1,1], UBOUNDS=[11,11], PERSISTENT=.TRUE.)
        CALL O%GET_HOST_DATA_RDWR(PTR_HOST)
        PTR_HOST=7

        IF (.NOT. ALL(PTR_HOST == 7)) THEN
                ERROR STOP
        END IF 

        CALL O%GET_DEVICE_DATA_RDONLY(PTR_DEV)

        OKAY=.TRUE.
        !$ACC PARALLEL PRESENT(PTR_DEV) COPY(OKAY)
        DO I=1,11
        DO J=1,11
          IF(PTR_DEV(I,J)/=7)THEN
                  OKAY=.FALSE.
          ENDIF
        ENDDO
        ENDDO
        !$ACC END PARALLEL
        IF(OKAY .EQV. .FALSE.)THEN
                ERROR STOP
        ENDIF
        CALL FIELD_DELETE(O)
END PROGRAM CPU_TO_GPU
