! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE GET_PRESET_WGRIB_TEMPLATE(CT, IGRIB_HANDLE, NGRIBV)

!----------------------------------------------------------------------

!**** *GET_PRESET_WGRIB_TEMPLATE* GETS THE WAVE MODEL GRIB TEMPLATE
!     WHICH HAS BEEN CREATED BY THE WAVE MODEL.

!     J. HAWKES    ECMWF NVOEMBER 20017

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!     SUBROUTINE GET_PRESET_WGRIB_TEMPLATE(CT, IGRIB_HANDLE)
!                INPUT:
!                CT           : "I" for INTEGRATED PARAMETERS AND
!                               "2" for INTGRATED PARAMETERS in GRIB 2                            
!                               "S" for SPECTRA
!                OPTIONAL INPUT:
!                NGRIBV         GRIB VERSION TO BE USED (if absent then = NGRIB_VERSION)
!                OUTPUT:
!                IGRIB_HANDLE : GRIB HANDLE THAT WILL BE CREATED. 

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWGRIB_HANDLES, ONLY : NGRIB_HANDLE_WAM_S, NGRIB_HANDLE_WAM_I, NGRIB_HANDLE_WAM_I2
      USE YOWABORT, ONLY : WAM_ABORT

      IMPLICIT NONE

      CHARACTER(LEN=1), INTENT(IN) :: CT 
      INTEGER, INTENT(IN), OPTIONAL :: NGRIBV
      INTEGER(KIND=JWIM), INTENT(OUT) :: IGRIB_HANDLE

      INTEGER(KIND=JWIM) :: IGRIB_VERSION
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-------------------------------------------------------------------

    IF (LHOOK) CALL DR_HOOK('GET_PRESET_WGRIB_TEMPLATE',0,ZHOOK_HANDLE)

    IF ( PRESENT(NGRIBV) ) THEN
      IGRIB_VERSION = NGRIBV  
    ELSE
      IGRIB_VERSION = -1
    ENDIF


    IF(CT == "S") THEN
        IGRIB_HANDLE = NGRIB_HANDLE_WAM_S
    ELSE IF (CT == "I") THEN
      IF ( IGRIB_VERSION == 2 ) THEN
        IGRIB_HANDLE = NGRIB_HANDLE_WAM_I2
      ELSE
        IGRIB_HANDLE = NGRIB_HANDLE_WAM_I
      ENDIF
    ELSE
      CALL WAM_ABORT(' GET_PRESET_WGRIB_TEMPLATE: Value of CT not recognized.', &
        & __FILENAME__, __LINE__)
    END IF

    IF (LHOOK) CALL DR_HOOK('GET_PRESET_WGRIB_TEMPLATE',1,ZHOOK_HANDLE)

END SUBROUTINE GET_PRESET_WGRIB_TEMPLATE
