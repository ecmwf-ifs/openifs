! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      PROGRAM chief 

! ----------------------------------------------------------------------

!**** *CHIEF* - SUPERVISES WAVE MODEL EXECUTION.

!     LIANA ZAMBRESKY      GKSS/ECMWF  JUNE 1989
!     H. GUNTHER           ECMWF       JUNE 1990  MODIFIED FOR CYCLE_4.
!     J. BIDLOT            ECMWF       FEBRUARY 1996-97 MESSAGE PASSING
!     J. DOYLE             ECMWF       OCTOBER 1996 ATMOSPHERIC COUPLING
!     J. BIDLOT            ECMWF       APRIL 97 ADD ZDELATM TO CALL TO
!                                               WAVEMDL
!     B. HANSEN            ECMWF       APRIL 97 SIGNAL HANDLING.
!     S. ABDALLA           ECMWF       OCTOBER 2001 MODIFICATION OF THE
!                                                   CALL TO WAVEMDL

!*    PURPOSE.
!     --------

!       THIS PROGRAM SUPERVISES THE EXECUTION OF THE WAM MODEL.

!**   INTERFACE.
!     ----------

!       IN ORDER FOR THE WAM MODEL TO EXECUTE, IT NEEDS
!       FILES FROM ESSENTIALLY FIVE SOURCES.

!       1. THE UNFORMATED FILES CREATED BY THE JOB PREPROC

!       2. USER INPUT FILE

!       3. THE WIND INPUT FILE.

!       4  THE BOUNDARY VALUE INPUT FILES CREATED BY JOB BOUINT.
!          THESE FILES ARE DYNAMICALLY ASSIGNED.

!       5. THE START FILES:
!          THE RESTART FILES HAVE TO BE CREATED BY JOB
!          PRESET, IF A COLD START HAS TO BE DONE.
!          THESE FILES OR FILES FROM A PREVIOUS MODEL RUN
!          ARE AUTOMATICALLY ASSIGNED. (SEE SUB GSFILE).

!       EXPLANATIONS FOR ALL FILES ARE GIVEN IN DETAIL IN SUB INITMDL

!     LIBRARIES.
!     ----------

!         NONE.

!     METHOD.
!     -------

!       THIS VERSION OF THE WAM MODEL HAS BEEN PRODUCED
!       BY MERGING AND CORRECTLY INTERFACING WHAT USED
!       TO BE THE STAND ALONE PROGRAMS:
!               PREWIND AND THE WAM MODEL.
!       PREWIND REFORMATS WINDS INTO THE WAM MODEL BLOCKED
!       STRUCTURE.  STARTING WITH THE INITIAL SEA STATE
!       FILES, THE WAM MODEL CAN THEN INTEGRATE FORWARD
!       IN TIME, DRIVEN BY THE REFORMATTED WINDS.
!       THE SEA STATE AND RESULT FILES ARE SAVED IN REGULAR
!       INTERVALLS. THE SEA STATE FILE SERVE AS THE INITIAL
!       CONDITION FOR A RESTART.

!       EACH CALL OF THE SUB WAVEMDL INTEGRATES FORWARD IN
!       TIME BY ONE WIND INPUT TIMESTEP OR ONE PROPAGATION
!       TIMESTEP, WHAT EVER IS LONGER.
!       IN THE FIRST CALL TO WAVEMDL AN INITIALIZATION IS
!       DONE IN ADDITION.

!     EXTERNALS.
!     ----------

!       *WAVEMDL*   - SUPERVISES THE OVERALL FLOW THROUGH
!                     THE MAIN MODULES: INITMDL, PREWIND
!                     AND WAMODEL.

!     REFERENCE.
!     ----------

!       EACH MODULE IS OF ITSELF THOROUGHLY DOCUMENTED.

! ----------------------------------------------------------------------

      USE YOMHOOK,    ONLY : LHOOK, DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_INIT, MPL_END
      USE YOWNEMOIO,  ONLY : LNEMOIOSERVER, WAMININEMOIO, WAMENDNEMOIO

! ----------------------------------------------------------------------

      IMPLICIT NONE
#ifdef WAM_HAVE_ASSI
#include "wam_setup_assi.intfb.h"
#endif
#include "runwam.intfb.h"

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

!     Initialise NEMOIO, which is enabled with environment variable 'NEMOIOSERVER=1'
      CALL WAMININEMOIO

!     Initialise MPI if not done yet
      CALL MPL_INIT

!     IO-servers skip execution
      IF (.NOT.LNEMOIOSERVER) THEN

!       DR_HOOK must be initialised after MPI
        IF (LHOOK) CALL DR_HOOK('CHIEF',0,ZHOOK_HANDLE)

#ifdef WAM_HAVE_ASSI
!       Set up callback functions for data assimilation
        CALL WAM_SETUP_ASSI
#endif
        CALL RUNWAM

        IF (LHOOK) CALL DR_HOOK('CHIEF',1,ZHOOK_HANDLE)
      ENDIF

!     Close NEMO IO servers
      CALL WAMENDNEMOIO

!     Finalise MPI if not done yet
      CALL MPL_END(LDMEMINFO=.FALSE.)

      END PROGRAM chief
