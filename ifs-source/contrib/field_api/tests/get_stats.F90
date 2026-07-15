! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM GET_STATS
        ! TEST IF STATS ARE CORRECTLY COMPUTED
        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE PARKIND1
        IMPLICIT NONE
        CLASS(FIELD_2RB), POINTER :: W => NULL()
        REAL(KIND=JPRB), ALLOCATABLE :: D(:,:)
        REAL(KIND=JPRB), POINTER :: D_GPU(:,:)
        REAL(KIND=JPRB), POINTER :: D_CPU(:,:)

        ALLOCATE(D(10,10))

        CALL FIELD_NEW(W, DATA=D)
        CALL W%GET_DEVICE_DATA_RDONLY(D_GPU)
        CALL W%GET_DEVICE_DATA_RDONLY(D_GPU)
        ! SHOULD BE 1 BECAUSE THE DATA HAVE NOT CHANGED BETWEEN THE TWO MOVE
        IF (W%STATS%TRANSFER_CPU_TO_GPU /= 1) THEN
                ERROR STOP
        END IF
        IF (W%STATS%TOTAL_TIME_TRANSFER_CPU_TO_GPU == 0) THEN
                ERROR STOP
        END IF

        CALL W%GET_HOST_DATA_RDONLY(D_CPU)
        ! SHOULD BE ZERO, SINCE WE HAVE NOT MODIFED DATA ON THE GPU, THERE IS NO
        ! NEED TO TRANSFER DATA BACK ON CPU, THEY SHOULD BE STILL THE SAME
        IF (W%STATS%TRANSFER_GPU_TO_CPU /= 0) THEN
                ERROR STOP
        END IF

        ! PRETEND WE WILL MODIFIED THE DATA ON THE GPU
        CALL W%GET_DEVICE_DATA_RDWR(D_GPU)
        CALL W%GET_HOST_DATA_RDONLY(D_CPU)
        ! SHOULD BE STILL ONE, BECAUSE WE ALREADY MOVED DATA ON THE GPU BEFORE
        ! AND THEY HAVE NOT BEEN MODIFIED SINCE
        IF (W%STATS%TRANSFER_CPU_TO_GPU /= 1) THEN
                ERROR STOP
        END IF
        IF(W%STATS%TRANSFER_GPU_TO_CPU /= 1 ) THEN
                ERROR STOP
        END IF
        CALL FIELD_DELETE(W)
END PROGRAM GET_STATS
