! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_END_MOD

!**** MPL_END - Terminates the message passing environment

!     Purpose.
!     --------
!     Cleans up all of the MPI state. 
!     Subsequently, no MPI routine can be called

!**   Interface.
!     ----------
!        CALL MPL_END

!        Input required arguments :
!        -------------------------
!           none

!        Input optional arguments :
!        -------------------------
!           LDMEMINFO -  print memory info (default True)

!        Output required arguments :
!        -------------------------
!           none

!        Output optional arguments :
!        -------------------------
!           KERROR    -  return error code.     If not supplied, 
!                        MPL_END aborts when an error is detected.
!     Author.
!     -------
!        D.Dent, M.Hamrud     ECMWF

!     Modifications.
!     --------------
!        Original: 2000-09-01
!        P. Towers     3-Jul-2014 Add call to ec_cray_meminfo

!     ------------------------------------------------------------------

USE EC_PARKIND  ,ONLY : JPIM

USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD

IMPLICIT NONE

PUBLIC MPL_END
PRIVATE

INTEGER :: IERROR

CONTAINS 

SUBROUTINE MPL_END(KERROR,LDMEMINFO)


#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_BUFFER_DETACH => MPI_BUFFER_DETACH8, MPI_FINALIZE => MPI_FINALIZE8
#endif


INTEGER(KIND=JPIM),INTENT(OUT),OPTIONAL :: KERROR
LOGICAL           ,INTENT(IN), OPTIONAL :: LDMEMINFO
INTEGER(KIND=JPIM)                      :: IBUFFMPI,IERROR
LOGICAL                      :: LLMEMINFO
LOGICAL,parameter            :: LLABORT=.TRUE.

#include "ec_mpi_finalize.intfb.h"

IF(MPL_NUMPROC < 1) THEN
  IF(MPL_NUMPROC == -1) THEN
    IF (.NOT.LINITMPI_VIA_MPL) THEN
      ! Neither MPL_INIT_MOD nor MPL_ARG_MOD -modules were called before this
      CALL MPL_MESSAGE(CDMESSAGE=' MPL_END CALLED BEFORE MPL_INIT ')
    ENDIF
!!-- we do not want the following message to appear, since its non-fatal
!!  ELSEIF(MPL_NUMPROC == -2) THEN
!!    CALL MPL_MESSAGE(CDMESSAGE=' MPL_END CALLED MULTIPLE TIMES ')
  ENDIF
  IF(PRESENT(KERROR)) THEN
    IERROR=0
    KERROR=IERROR
  ENDIF
  RETURN
ENDIF

IF (ALLOCATED(MPL_ATTACHED_BUFFER)) THEN
  IF( MPI_IS_FINALIZED() ) THEN
    CALL MPL_MESSAGE(CDMESSAGE='MPL_END -- Cannot call MPI_Buffer_detach() as MPI is already finalized',LDABORT=.FALSE.)
  ELSE
    IBUFFMPI=SIZE(MPL_ATTACHED_BUFFER) * JP_ATTACHED_BUFFER_BYTES ! in bytes
    CALL MPI_BUFFER_DETACH(MPL_ATTACHED_BUFFER,IBUFFMPI,IERROR)
    IF(PRESENT(KERROR)) THEN
      KERROR=IERROR
    ELSE
      IF( IERROR /= 0 )THEN
        CALL MPL_MESSAGE(IERROR,'MPL_END ',LDABORT=LLABORT)
      ENDIF
    ENDIF
  ENDIF
  DEALLOCATE(MPL_ATTACHED_BUFFER)
ENDIF

LLMEMINFO=.TRUE.
IF(PRESENT(LDMEMINFO)) LLMEMINFO=LDMEMINFO
CALL EC_MPI_FINALIZE(IERROR,LINITMPI_VIA_MPL,LLMEMINFO,"mpl_end")

MPL_NUMPROC = -2
LINITMPI_VIA_MPL = .FALSE.

IF(PRESENT(KERROR)) THEN
  KERROR=IERROR
ENDIF

RETURN
END SUBROUTINE MPL_END

FUNCTION MPI_IS_FINALIZED()
  LOGICAL :: MPI_IS_FINALIZED
  LOGICAL :: LLINIT, LLFIN
  INTEGER(KIND=JPIM) :: IERR
  MPI_IS_FINALIZED = .FALSE.
  CALL MPI_INITIALIZED(LLINIT,IERR)
  IF (LLINIT .AND. IERR == 0) THEN
    CALL MPI_FINALIZED(LLFIN,IERR)
    IF( IERR == 0 ) THEN
      MPI_IS_FINALIZED = LLFIN
    ENDIF
  ENDIF
END FUNCTION

END MODULE MPL_END_MOD
