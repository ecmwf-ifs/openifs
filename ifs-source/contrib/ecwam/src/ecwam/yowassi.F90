! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE YOWASSI
    !! This module provides interfaces for callback functions that need to
    !! be registered for Data assimilation routines
    !! The idea is that data assimilation can be added via an external library
    !! which provides implementations for the interfaces.
    !! A setup routine needs to point the HANDLERs to the correct implementation

IMPLICIT NONE
PRIVATE

PUBLIC :: WAMASSI, WAM_ODB_OPEN, WAM_ODB_CLOSE
PUBLIC :: WAMASSI_HANDLER, WAM_ODB_OPEN_HANDLER, WAM_ODB_CLOSE_HANDLER

PROCEDURE(WAMASSI),       POINTER :: WAMASSI_HANDLER       => NULL()
PROCEDURE(WAM_ODB_OPEN),  POINTER :: WAM_ODB_OPEN_HANDLER  => NULL()
PROCEDURE(WAM_ODB_CLOSE), POINTER :: WAM_ODB_CLOSE_HANDLER => NULL()

CONTAINS

SUBROUTINE WAMASSI(LDSTOP, LDWRRE, BLK2GLO, WVENVI,   &
    &              WVPRPT, FF_NOW, INTFLDS,  &
    &              WAM2NEMO, NEMO2WAM, FL1)
    USE PARKIND_WAVE, ONLY : JWRB
    USE YOWDRVTYPE  , ONLY : WVGRIDGLO, ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
    &                        INTGT_PARAM_FIELDS, WAVE2OCEAN, OCEAN2WAVE
    USE YOWGRID     , ONLY : NPROMA_WAM, NCHNK
    USE YOWMAP      , ONLY : NIBLO
    USE YOWPARAM    , ONLY : NANG ,NFRE
    USE YOWABORT    , ONLY : WAM_ABORT

! ----------------------------------------------------------------------
    IMPLICIT NONE

    LOGICAL,                  INTENT(IN)    :: LDSTOP
    LOGICAL,                  INTENT(IN)    :: LDWRRE
    TYPE(WVGRIDGLO),          INTENT(IN)    :: BLK2GLO
    TYPE(ENVIRONMENT),        INTENT(INOUT) :: WVENVI
    TYPE(FREQUENCY),          INTENT(IN)    :: WVPRPT
    TYPE(FORCING_FIELDS),     INTENT(INOUT) :: FF_NOW
    TYPE(INTGT_PARAM_FIELDS), INTENT(INOUT) :: INTFLDS
    TYPE(WAVE2OCEAN),         INTENT(INOUT) :: WAM2NEMO
    TYPE(OCEAN2WAVE),         INTENT(IN)    :: NEMO2WAM
    REAL(KIND=JWRB),          INTENT(INOUT) :: FL1(NPROMA_WAM, NANG, NFRE, NCHNK)

    IF( .NOT. ASSOCIATED(WAMASSI_HANDLER) ) THEN
        CALL WAM_ABORT("WAMASSI_HANDLER is not associated. Make sure WAM_SETUP_ASSI was called.", &
            & __FILENAME__, __LINE__)
    ELSE
        CALL WAMASSI_HANDLER(LDSTOP, LDWRRE, BLK2GLO, WVENVI,   &
            &                WVPRPT, FF_NOW, INTFLDS,  &
            &                WAM2NEMO, NEMO2WAM, FL1)
    ENDIF
END SUBROUTINE

SUBROUTINE WAM_ODB_OPEN
    USE YOWABORT    , ONLY : WAM_ABORT
    IF( .NOT. ASSOCIATED(WAM_ODB_OPEN_HANDLER) ) THEN
        CALL WAM_ABORT("WAM_ODB_OPEN_HANDLER is not associated. Make sure WAM_SETUP_ASSI was called.", &
            & __FILENAME__, __LINE__)
    ELSE
        CALL WAM_ODB_OPEN_HANDLER()
    ENDIF
END SUBROUTINE

SUBROUTINE WAM_ODB_CLOSE
    USE YOWABORT    , ONLY : WAM_ABORT
    IF( .NOT. ASSOCIATED(WAM_ODB_CLOSE_HANDLER) ) THEN
        CALL WAM_ABORT("WAM_ODB_CLOSE_HANDLER is not associated. Make sure WAM_SETUP_ASSI was called.", &
            & __FILENAME__, __LINE__)
    ELSE
        CALL WAM_ODB_CLOSE_HANDLER()
    ENDIF
END SUBROUTINE

end module
