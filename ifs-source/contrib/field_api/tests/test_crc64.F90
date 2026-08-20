! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM TEST_CRC64
        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE FIELD_UTIL_MODULE
        USE OMP_LIB
        USE PARKIND1

        CLASS(FIELD_2RB), POINTER :: W => NULL()
        REAL(KIND=JPRB), POINTER :: D(:,:)
        INTEGER :: NPROMA = 24
        INTEGER :: NBLOCKS= 100
        INTEGER :: IBLK,JLON
        REAL(KIND=JPRB), POINTER :: VIEW(:) => NULL()

        INTEGER*8 :: ICRC64_PI, ICRC64_E, ICRC
        DATA ICRC64_PI / Z'33E2D12C420E6E86' /
        DATA ICRC64_E  / Z'817AF6E200A0FC70' /


        ALLOCATE(D(NPROMA, NBLOCKS))

        CALL FIELD_NEW(W, DATA=D)

        D => GET_HOST_DATA_RDWR (W)

        D = 3.14

        ICRC = CRC64 (W)
!       WRITE (*, '(Z16.16)') ICRC
        IF (ICRC /= ICRC64_PI) ERROR STOP

        D => GET_DEVICE_DATA_RDWR (W)

!$acc kernels present (D)

        D = 2.78

!$acc end kernels

        ICRC = CRC64 (W)
!       WRITE (*, '(Z16.16)') ICRC
        IF (ICRC /= ICRC64_E) ERROR STOP
        

        CALL FIELD_DELETE(W)
END PROGRAM TEST_CRC64
