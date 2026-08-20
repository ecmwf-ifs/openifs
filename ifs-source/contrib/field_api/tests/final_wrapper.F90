! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM FINAL_WRAPPER
        ! TEST IF DESTROYING THE WRAPPER DOESN'T DESTROY THE WRAPPED DATA

        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE PARKIND1
        IMPLICIT NONE
        CLASS(FIELD_2RB), POINTER :: W => NULL()
        REAL(KIND=JPRB), ALLOCATABLE :: D(:,:)

        ALLOCATE(D(10,10))
        D=7

        CALL FIELD_NEW(W, DATA=D)
        CALL W%FINAL()

        IF (ASSOCIATED(W%PTR)) THEN
                ERROR STOP
        END IF
        IF (ASSOCIATED(W%DEVPTR)) THEN
                ERROR STOP
        END IF

        IF (.NOT. ALLOCATED(D)) THEN
                ERROR STOP
        END IF

        IF (.NOT. ALL(D == 7)) THEN
                ERROR STOP
        END IF 
END PROGRAM FINAL_WRAPPER
