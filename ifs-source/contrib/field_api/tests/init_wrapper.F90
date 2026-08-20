! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM INIT_WRAPPER
        ! TEST IF WRAPPER IS REALLY WRAPPING DATA
        ! IF IT IS WRAPPING THEN MODIFYING THE WRAPED DATA SHOULD BE SEEN IN
        ! WRAPPER TOO

        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE PARKIND1
        IMPLICIT NONE
        CLASS(FIELD_2RB), POINTER :: W => NULL()
        REAL(KIND=JPRB), ALLOCATABLE :: D(:,:)

        ALLOCATE(D(10,10))
        D=7

        CALL FIELD_NEW(W, DATA=D)
        D=42

        IF (.NOT. ALL(W%PTR == 42)) THEN
                ERROR STOP
        END IF 
        CALL FIELD_DELETE(W)
END PROGRAM INIT_WRAPPER
