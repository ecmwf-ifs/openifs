! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM INIT_OWNER_OPENMP
        ! TEST IF OWNER IS REALLY ALLOCATING THE DATA
        ! WHEN PERSITENT IS SET TO FALSE OR IS NOT GIVEN IN ARGUMENT,
        ! THEN THE LAST DIM OF THE FIELD IS THE NUMBER OF OPENMP THREADS

        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE OMP_LIB
        USE PARKIND1
        IMPLICIT NONE
        CLASS(FIELD_2RB), POINTER :: O => NULL()
        REAL(KIND=JPRB), POINTER :: PTR(:,:)

        CALL OMP_SET_NUM_THREADS(4)
        CALL FIELD_NEW(O, LBOUNDS=[10,1], UBOUNDS=[21,11])
        CALL O%GET_HOST_DATA_RDWR(PTR)
        PTR=42

        IF (SIZE(O%PTR,1) /= 12) THEN
                ERROR STOP
        END IF
        !4 BECAUSE WE SET THE NUMBER OF THREADS TO 4
        IF (SIZE(O%PTR,2) /= 4) THEN
                ERROR STOP
        END IF

        IF (.NOT. ALL(O%PTR == 42)) THEN
                ERROR STOP
        END IF 
        CALL FIELD_DELETE(O)
END PROGRAM INIT_OWNER_OPENMP
