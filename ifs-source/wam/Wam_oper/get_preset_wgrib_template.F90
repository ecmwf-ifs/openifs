! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#define __FILENAME__ "get_preset_wgrib_template.F90"

SUBROUTINE GET_PRESET_WGRIB_TEMPLATE(CT,IGRIB_HANDLE)

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
!                               "S" for SPECTRA
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
      USE YOWGRIB_HANDLES, ONLY : NGRIB_HANDLE_WAM_S, NGRIB_HANDLE_WAM_I
      USE YOWABORT, ONLY : WAM_ABORT

      IMPLICIT NONE

      CHARACTER(LEN=1), INTENT(IN) :: CT 
      INTEGER(KIND=JWIM), INTENT(OUT) :: IGRIB_HANDLE

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-------------------------------------------------------------------

    IF (LHOOK) CALL DR_HOOK('GET_PRESET_WGRIB_TEMPLATE',0,ZHOOK_HANDLE)

    IF(CT.EQ."S") THEN
        IGRIB_HANDLE = NGRIB_HANDLE_WAM_S
    ELSE IF (CT.EQ."I") THEN
        IGRIB_HANDLE = NGRIB_HANDLE_WAM_I
    ELSE
      CALL WAM_ABORT(' GET_PRESET_WGRIB_TEMPLATE: Value of CT not recognized.', &
        & __FILENAME__, __LINE__)
    END IF

    IF (LHOOK) CALL DR_HOOK('GET_PRESET_WGRIB_TEMPLATE',1,ZHOOK_HANDLE)

END SUBROUTINE GET_PRESET_WGRIB_TEMPLATE 
