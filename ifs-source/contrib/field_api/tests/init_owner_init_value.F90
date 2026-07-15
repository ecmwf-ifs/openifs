! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM INIT_OWNER_INIT_VALUE
        !TEST THAT ALL DATA ARE INITIALIZE TO THE CHOOSEN VALUE
        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE PARKIND1
        USE OMP_LIB, ONLY: OMP_GET_MAX_THREADS
        IMPLICIT NONE
        CLASS(FIELD_2IM), POINTER :: O => NULL()
        INTEGER(KIND=JPIM), POINTER :: PTR(:,:)

        CALL FIELD_NEW(O, LBOUNDS=[1,1], UBOUNDS=[10,10], INIT_VALUE=7_JPIM)
        CALL O%GET_HOST_DATA_RDWR(PTR)

        WRITE(*,*)O%PTR(1,1)
        IF (.NOT. ALL(O%PTR == 7)) THEN
                ERROR STOP
        END IF

        CALL FIELD_DELETE(O)
END PROGRAM INIT_OWNER_INIT_VALUE
