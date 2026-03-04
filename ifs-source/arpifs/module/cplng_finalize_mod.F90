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

MODULE CPLNG_FINALIZE_MOD

   IMPLICIT NONE

   PRIVATE

   PUBLIC CPLNG_FINALIZE

CONTAINS

   SUBROUTINE CPLNG_FINALIZE(LDACTIVE)

      USE PARKIND1, ONLY: JPIM
#ifdef WITH_OASIS
      USE MOD_OASIS
#endif
      USE CPLNG_LOG_MOD
      USE CPLNG_DATA_MOD
      USE YOMMCC, ONLY : TMCC

      ! Argument
      LOGICAL, INTENT(IN) :: LDACTIVE

      INTEGER(KIND=JPIM) :: IERROR
      CHARACTER(LEN=3)   :: CLERRSTR

      ! Early return if CPLNG not active
      IF (.NOT.LDACTIVE) THEN
         CALL CPLNG_WARNING('CPLNG_FINALIZE called but CPLNG_ACTIVE is false.')
         RETURN
      ENDIF

#ifdef WITH_OASIS

      ! Finalize by calling OASIS terminate routine
      CALL OASIS_TERMINATE(IERROR)

      ! Check error code
      IF (IERROR/=OASIS_OK) THEN
         WRITE (CLERRSTR,'(I3)') IERROR
         CALL ABOR1("CPLNG_FINALIZE: OASIS_TERMINATE returns error code: "//CLERRSTR)
      ENDIF

#else

      CALL ABOR1('CPLNG_FINALIZE: SGLEXE NOT DONE YET')

#endif

   END SUBROUTINE CPLNG_FINALIZE

END MODULE CPLNG_FINALIZE_MOD
