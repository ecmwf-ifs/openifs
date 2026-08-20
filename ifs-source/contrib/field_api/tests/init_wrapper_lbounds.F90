! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM INIT_WRAPPER_LBOUNDS
        ! TEST IF WRAPPER IS USING THE ORIGINAL LBOUNDS
        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE PARKIND1
        IMPLICIT NONE
        CLASS(FIELD_2RB), POINTER :: W => NULL()
        REAL(KIND=JPRB), ALLOCATABLE :: D(:,:,:)
        INTEGER :: LBOUNDS(2)

        ALLOCATE(D(0:10, 5:10, 100:200))
        CALL FIELD_NEW(W, DATA=D(:,7,:), LBOUNDS=[0, 100])
        LBOUNDS = LBOUND(W%PTR)

        IF (LBOUNDS(1) /= 0 .OR. LBOUNDS(2) /= 100) THEN
                ERROR STOP
        END IF 
        CALL FIELD_DELETE(W)
END PROGRAM INIT_WRAPPER_LBOUNDS
