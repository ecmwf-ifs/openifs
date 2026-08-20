! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM FINAL_OWNER
        ! TEST IF OWNER IS DESTROYED WHEN CALLING FINAL

        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE PARKIND1
        IMPLICIT NONE
        CLASS(FIELD_2RB), POINTER :: O => NULL()
        REAL(KIND=JPRB), POINTER :: PTR(:,:)

        CALL FIELD_NEW(O, [10,1],[21,11])
        CALL O%GET_HOST_DATA_RDWR(PTR)
        PTR=42
        CALL O%FINAL()

        IF (ASSOCIATED(O%PTR)) THEN
                ERROR STOP
        END IF
        IF (ASSOCIATED(O%DEVPTR)) THEN
                ERROR STOP
        END IF
        CALL FIELD_DELETE(O)
END PROGRAM FINAL_OWNER
