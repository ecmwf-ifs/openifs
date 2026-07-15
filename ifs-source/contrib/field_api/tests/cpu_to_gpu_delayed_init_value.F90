! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM CPU_TO_GPU_DELAYED_INIT_VALUE
        !TEST THAT DATA OF A FIELD OWNER IS NOT COPIED TO CPU WHEN IT HAS
        !BEEN INITIALISED WITH AN INIT VALUE AFTER A DELAYED ALLOCATION, BUT
        !THAT IT IS TRANSFERED IF WE WANT TO READ IT FROM GPU

        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE PARKIND1
        IMPLICIT NONE
        CLASS(FIELD_2IM), POINTER :: O => NULL()
        INTEGER(KIND=JPIM), POINTER :: PTR_HOST(:,:) => NULL()
        INTEGER(KIND=JPIM), POINTER :: PTR_DEV(:,:) => NULL()
        INTEGER ::I,J
        LOGICAL :: OKAY

        CALL FIELD_NEW(O, LBOUNDS=[1,1],UBOUNDS=[5,5], PERSISTENT=.TRUE., DELAYED=.TRUE., INIT_VALUE=3)

        CALL O%GET_HOST_DATA_RDONLY(PTR_HOST)
        WRITE(*,*)"O%STATS%TRANSFER_CPU_TO_GPU",O%STATS%TRANSFER_CPU_TO_GPU
        WRITE(*,*)"O%STATS%TRANSFER_GPU_TO_CPU",O%STATS%TRANSFER_GPU_TO_CPU
        IF(O%STATS%TRANSFER_GPU_TO_CPU /= 0)THEN
                ERROR STOP
        ENDIF
        DO I=1,5
        DO J=1,5
        IF (PTR_HOST(I,J) /= 3)THEN
                ERROR STOP
        ENDIF
        ENDDO
        ENDDO

        CALL O%GET_DEVICE_DATA_RDWR(PTR_DEV)
        WRITE(*,*)"O%STATS%TRANSFER_CPU_TO_GPU",O%STATS%TRANSFER_CPU_TO_GPU
        WRITE(*,*)"O%STATS%TRANSFER_GPU_TO_CPU",O%STATS%TRANSFER_GPU_TO_CPU
        IF(O%STATS%TRANSFER_CPU_TO_GPU /= 1)THEN
                ERROR STOP
        ENDIF
        OKAY=.TRUE.
        !$ACC PARALLEL PRESENT(PTR_DEV) COPY(OKAY)
        DO I=1,5
        DO J=1,5
        IF (PTR_DEV(I,J) /= 3)THEN
                OKAY=.FALSE.
        ENDIF
        ENDDO
        ENDDO
        !$ACC END PARALLEL
        IF(OKAY .EQV. .FALSE.)THEN
                ERROR STOP
        ENDIF

        CALL FIELD_DELETE(O)
END PROGRAM CPU_TO_GPU_DELAYED_INIT_VALUE
