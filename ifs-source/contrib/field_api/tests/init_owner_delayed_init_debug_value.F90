! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM INIT_OWNER_DELAYED_INIT_DEBUG_VALUE
        !TEST THAT ALL DATA ARE INITIALIZE TO THE CHOOSEN VALUE
        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE FIELD_INIT_DEBUG_VALUE_MODULE
        USE PARKIND1
        USE OMP_LIB, ONLY: OMP_GET_MAX_THREADS
        IMPLICIT NONE
        CLASS(FIELD_2IM), POINTER :: O => NULL()
        INTEGER(KIND=JPIM), POINTER :: PTR(:,:)

        USE_INIT_DEBUG_VALUE = .TRUE.
        INIT_DEBUG_VALUE_JPIM=128128128

        CALL FIELD_NEW(O, LBOUNDS=[1,1], UBOUNDS=[10,10], DELAYED=.TRUE.)
        CALL O%GET_HOST_DATA_RDWR(PTR)
        IF (.NOT. ALL(PTR == 128128128)) THEN
                ERROR STOP
        END IF

        CALL FIELD_DELETE(O)
END PROGRAM INIT_OWNER_DELAYED_INIT_DEBUG_VALUE
