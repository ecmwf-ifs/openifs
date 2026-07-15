! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM INIT_WRAPPER_MIXED_PRECISION
        ! TEST IF WRAPPER IS REALLY WRAPPING DATA
        ! IF IT IS WRAPPING THEN MODIFYING THE WRAPED DATA SHOULD BE SEEN IN
        ! WRAPPER TOO

        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE PARKIND1
        IMPLICIT NONE
        CLASS(FIELD_2RB), POINTER :: W => NULL()
        CLASS(FIELD_2RD), POINTER :: W_DBLE => NULL()
        REAL(KIND=JPRB), ALLOCATABLE :: D(:,:)
        REAL(KIND=JPRD), ALLOCATABLE :: D_DBLE(:,:)

        ALLOCATE(D(10,10))
        ALLOCATE(D_DBLE(10,10))
        D=7._JPRB
        D_DBLE=7._JPRD

        CALL FIELD_NEW(W, DATA=D)
        CALL FIELD_NEW(W_DBLE, DATA=D_DBLE)
        D=42._JPRB
        D_DBLE=42._JPRD

        IF (.NOT. ALL(W%PTR == 42._JPRB)) THEN
                ERROR STOP
        END IF 
        IF (.NOT. ALL(W_DBLE%PTR == 42._JPRD)) THEN
                ERROR STOP
        END IF 
        CALL FIELD_DELETE(W)
        CALL FIELD_DELETE(W_DBLE)

END PROGRAM INIT_WRAPPER_MIXED_PRECISION
