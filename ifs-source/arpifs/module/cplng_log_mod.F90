! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE CPLNG_LOG_MOD

    IMPLICIT NONE

    PRIVATE

    PUBLIC CPLNG_INFO
    PUBLIC CPLNG_WARNING
    PUBLIC CPLNG_ERROR

CONTAINS

SUBROUTINE CPLNG_INFO(CDMSG)

    USE YOMLUN, ONLY: NULOUT

    ! ARGUMENTS
    CHARACTER(LEN=*),INTENT(IN) :: CDMSG

    WRITE (NULOUT,'(/,X,A,/)') 'CPLNG: Info: ' // CDMSG

END SUBROUTINE CPLNG_INFO

SUBROUTINE CPLNG_WARNING(CDMSG)

    USE YOMLUN, ONLY: NULOUT

    ! ARGUMENTS
    CHARACTER(LEN=*),INTENT(IN) :: CDMSG

    WRITE (NULOUT,'(/,X,A,/)') 'CPLNG: Warning: ' // CDMSG

END SUBROUTINE CPLNG_WARNING

SUBROUTINE CPLNG_ERROR(CDMSG)

    USE YOMLUN, ONLY: NULOUT

    ! ARGUMENTS
    CHARACTER(LEN=*),INTENT(IN) :: CDMSG

    WRITE (NULOUT,'(/,X,A)') 'CPLNG: Error: ' // CDMSG
    WRITE (NULOUT,'(X,A)') 'CPLNG: Error: Calling ABOR1 to ABORT!'
    CALL ABOR1('CPLNG: ' // CDMSG)

END SUBROUTINE CPLNG_ERROR

END MODULE CPLNG_LOG_MOD
