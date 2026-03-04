! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE DR_HOOK_END()
  ! Make sure DrHook output is produced before MPI_Finalize (in case it fails)
  IMPLICIT NONE
  EXTERNAL :: c_drhook_prof
  LOGICAL,SAVE :: LL_FIRST_TIME = .TRUE.
  IF( .NOT. LL_FIRST_TIME ) THEN
    LL_FIRST_TIME = .FALSE.
    CALL c_drhook_prof()
  ENDIF
END SUBROUTINE