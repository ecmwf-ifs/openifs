! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE CHKOOPS(LDUPDATEOOPS)

! ----------------------------------------------------------------------

!**** *CHKOOPS* -

!*    PURPOSE.
!     --------
!      CHECK IF MODEL IS IN OOPS MODE AND ADAPT WHEN WAVE DATA ASSIMILATION IS TAKING PLACE

!**   INTERFACE.
!     ----------
!     *CALL CHKOOPS*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWSTAT  , ONLY : IASSI, IASSI_ORIG 
      USE YOWTEST  , ONLY : IU06
      USE YOWCOUT  , ONLY : LWAMANOUT, LWAMANOUT_ORIG
      USE YOWCOUP  , ONLY : IFSNUPTRA, IFSMUPTRA, IFSCONTEXT

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      LOGICAL, OPTIONAL, INTENT(IN) :: LDUPDATEOOPS

      INTEGER(KIND=JWIM), SAVE :: NTRAJ
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      LOGICAL :: LUPDATEOOPS
      LOGICAL, SAVE :: LFRST_OOPS

      DATA NTRAJ /-1/
      DATA LFRST_OOPS /.TRUE./

! ---------------------------------------------------------------------

!*    1.  THE FIRST CALL TO WAVEMDL PERFORMS INITIALIZATION.

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CHKOOPS',0,ZHOOK_HANDLE)

      LUPDATEOOPS = .FALSE.
      IF(PRESENT(LDUPDATEOOPS)) LUPDATEOOPS = LDUPDATEOOPS

      IF (LUPDATEOOPS .AND. TRIM(IFSCONTEXT) == 'OOPS') THEN
        ! OOPS-IFS may do wave assimilation only in the outer loop MUPTRA - 1
        IF (LFRST_OOPS) THEN
          LFRST_OOPS = .FALSE.
          IASSI_ORIG = IASSI
          LWAMANOUT_ORIG = LWAMANOUT
        ENDIF

        IF( NTRAJ /= IFSNUPTRA ) THEN
          IF (IFSNUPTRA /= IFSMUPTRA - 1) THEN
            IASSI = 0
            LWAMANOUT = .FALSE.
          ELSE
            IASSI = IASSI_ORIG
            LWAMANOUT = LWAMANOUT_ORIG
          ENDIF

          WRITE(IU06,*) ''
          WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++ '
          WRITE(IU06,*) ' SUB. CHKOOPS CALLED FROM ', IFSCONTEXT, &
 &                      ' FOR NUPTRA: ', IFSNUPTRA, ' AND MUPTRA: ', IFSMUPTRA
          WRITE(IU06,*) ' --> IASSI reset to', IASSI 
          WRITE(IU06,*) ' --> LWAMANOUT reset to', LWAMANOUT
          WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++ '

           NTRAJ = IFSNUPTRA
        ENDIF

      ENDIF

      IF (LHOOK) CALL DR_HOOK('CHKOOPS',1,ZHOOK_HANDLE)

      END SUBROUTINE CHKOOPS
