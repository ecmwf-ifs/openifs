! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE RUNWAM

! ----------------------------------------------------------------------

!**** *RUNWAM* - SUPERVISES WAVE MODEL EXECUTION.

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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU, LWCOU2W, LWCOURNW , LWCOUAST    ,    &
     &                      LWCOUHMF, LWFLUX, LWNEMOCOU            ,    &
     &                      NEMOINIDATE, NEMOINITIME               ,    &
     &                      NEMOITINI,   NEMOITEND                 ,    &
     &                      NEMOTSTEP,   NEMOFRCO                  ,    &
     &                      NEMONSTEP,   NEMOCSTEP, NEMOWSTEP
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWSTAT  , ONLY : CDATEE   ,CDTPRO                       ,    &
     &            IPROPAGS ,LSUBGRID ,IREFRA   ,IDELPRO, TIME_PHYS,     &
     &            TIME_PROPAG
      USE YOWWAMI  , ONLY : CBPLTDT  ,CEPLTDT
      USE YOWALTAS , ONLY : LODBRALT
      USE MPL_MODULE, ONLY : MPL_INIT, MPL_END, MPL_COMM
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWASSI, ONLY : WAM_ODB_OPEN,  WAM_ODB_CLOSE
      USE EC_PARKIND, ONLY : JPRD
      USE YOWTEST, ONLY : IU06
      USE YOWGSTATS, ONLY : WAM_GSTATS_SETUP, WAM_GSTATS_PRINT, &
     &                      WAM_GSTATS_FILE_OPEN, WAM_GSTATS_FILE_CLOSE

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "iniwcst.intfb.h"
#include "mpclose_unit.intfb.h"
#include "wavemdl.intfb.h"
#include "wvalloc.intfb.h"
#include "wvwamdecomp.intfb.h"
#include "wvwaminit.intfb.h"
#include "wvwaminit1.intfb.h"
#include "ec_meminfo.intfb.h"
#include "wam_user_clock.intfb.h"


! DIMENSION DUMMY COUPLED VARIABLES
      INTEGER(KIND=JWIM), PARAMETER :: NLONW=1
      INTEGER(KIND=JWIM), PARAMETER :: NLATW=1
      INTEGER(KIND=JWIM), PARAMETER :: NWVFIELDS=10
      INTEGER(KIND=JWIM), PARAMETER :: NC=1
      INTEGER(KIND=JWIM), PARAMETER :: NR=1
      INTEGER(KIND=JWIM), PARAMETER :: NGPTOTG=NC*NR
      INTEGER(KIND=JWIM), PARAMETER :: NFIELDS=7

      INTEGER(KIND=JWIM) :: NATMFLX
      INTEGER(KIND=JWIM) :: NGAUSSW, NLON, NLAT
      INTEGER(KIND=JWIM) :: IGRIB_HANDLE_DUM
      INTEGER(KIND=JWIM) :: IGRIB_HANDLE_DUM2
      INTEGER(KIND=JWIM) :: NADV
      INTEGER(KIND=JWIM) :: KSTOP, KSTPW
      INTEGER(KIND=JWIM) :: IDUM
      INTEGER(KIND=JWIM) :: NPR
      INTEGER(KIND=JWIM) :: MAXLEN
      INTEGER(KIND=JWIM) :: KERROR
      INTEGER(KIND=JWIM) :: MASK_IN(NGPTOTG)
      INTEGER(KIND=JWIM) :: MASK_OUT(NLONW,NLATW)

      REAL(KIND=JWRB) :: PSTEP
      REAL(KIND=JWRB) :: RSOUTW, RNORTW
      REAL(KIND=JWRB) :: PRPLRADI, PRPLRG
      REAL(KIND=JWRB) :: RNU_ATM, RNUM_ATM
      REAL(KIND=JWRB) :: WVFLDG(NLONW,NLATW,NWVFIELDS), ZDELATM(NLATW)
      REAL(KIND=JWRB) :: FIELDS(NGPTOTG,NFIELDS)
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: RMISS
      REAL(KIND=JWRB) :: ZRCHAR
      REAL(KIND=JWRU) :: time0, time, timestep_start, timestep0_start

      CHARACTER(LEN=3) :: DBNAME
      CHARACTER(LEN=14) :: ZERO,CBEGDAT

      LOGICAL :: LLSTOP, LLWRRE, LLRESTARTED, LLIRANK, LLNORMWAMOUT_GLOBAL
      LOGICAL :: LLINIT_WVFLDG
      LOGICAL :: LDWCOUNORMS
      LOGICAL :: FRSTIME
      LOGICAL :: LWCUR
      LOGICAL :: LWSTOKES
      LOGICAL :: LLRNL
      LOGICAL :: LFDBIFS

      LOGICAL :: LGSTATS = .TRUE. ! gstats statistics
      INTEGER :: IUGSTATS
      INTEGER :: ISTEP = 0

      DATA LLSTOP, LLWRRE, LLNORMWAMOUT_GLOBAL  / 3*.FALSE. /

! ----------------------------------------------------------------------


      time0=-wam_user_clock()
      TIME_PROPAG=0_JWRB
      TIME_PHYS=0_JWRB
      IU06=6

!     0.1 INITIALISE MESSAGE PASSING PROTOCOL
!         -----------------------------------

      CALL MPL_INIT(KERROR=KERROR)
      IF (KERROR < 0) THEN
        IU06=6
        WRITE (IU06,*) ' ******************************************'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' *      FATAL ERROR PROGRAM RUNWAM         *'
        WRITE (IU06,*) ' *      =========================         *'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' *            PROBLEM WITH                *'
        WRITE (IU06,*) ' *      MESSAGE PASSING INITIALISATION    *'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS       *'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' ******************************************'
        CALL ABORT1
      ENDIF

      IF (LHOOK) CALL DR_HOOK('RUNWAM',0,ZHOOK_HANDLE)

      CALL EC_MEMINFO(6, "runwam:init", KCOMM=MPL_COMM, KBARR=1, KIOTASK=-1, KCALL=0)

!     0.2 GET MODEL PARAMETERS
!         --------------------

      LWCOU=.FALSE.
      LWCOU2W=.FALSE.
      LWCOURNW=.FALSE.
      LWCOUAST=.FALSE.
      LWCOUHMF=.FALSE.
      LWFLUX=.FALSE. ! will be reset to true if ocean fluxes are output.
      LWCUR=.FALSE. ! only used in coupled runs with atmospheric model
      LFDBIFS=.FALSE.
      LLINIT_WVFLDG=.FALSE.


      PRPLRADI=1.0_JWRB
      PRPLRG=1.0_JWRB
      LLRNL=.TRUE.

!     KINEMATIC AIR VISCOSITY AND REDUCED VALUE FOR MOMENTUM TRANSFER
      RNU_ATM=1.5E-5
      RNUM_ATM=0.11_JWRB*RNU_ATM

      CALL WVWAMINIT (LWCOU,IU06,LLRNL,NGAUSSW,NLON,NLAT,RSOUTW,RNORTW)

      CALL WVWAMINIT1 (LWCOU, LWCOU2W, LWCOURNW, LWCOUHMF, LWFLUX, LFDBIFS)

      CALL WAM_GSTATS_SETUP()

      CALL GSTATS(0,0)

!     0.3 DETERMINE GRID DOMAIN DECOMPOSITION
!         -----------------------------------

      NADV=0
      FRSTIME=.TRUE.

      CALL INIWCST(PRPLRADI)

      CALL WVWAMDECOMP


!     1.1  ALLOCATE NECESSARY ARRAYS
!          -------------------------

      CALL WVALLOC

!     1.2 OPEN ODB DATABASE TO RETRIEVE SATELLITE DATA
!         --------------------------------------------
      IF (LODBRALT) THEN
        CALL WAM_ODB_OPEN()
      ENDIF

!     ------------------------------------------------------------------

!*    2. CALLS TO WAVEMDL UNTIL MODEL DATE REACHES END DATE.
!*       EACH CALL INTEGRATES ONE WIND INPUT TIMESTEP, OR ONE
!*       PROPAGATION TIMESTEP, WHAT EVER IS LONGER.
!        ---------------------------------------------------

      ZERO   = ' '
      CDATEE = ZERO
      CDTPRO = ZERO

!* DEFINE DUMMY PARAMETERS TO FILL COUPLED ARRAYS

      CBEGDAT='99999999999999'
      PSTEP=0.0_JWRB ! only used in coupled model
      KSTOP=0  ! only used in coupled model
      KSTPW=0  ! only used in coupled model
      IDUM=0
      IGRIB_HANDLE_DUM=-99 ! only used in coupled model
      IGRIB_HANDLE_DUM2=-99 ! only used in coupled model
      NATMFLX=0
      LWSTOKES=.FALSE.  ! only used in coupled runs with atmospheric model
      RMISS=-999.0_JWRB ! missing data indicator
      ZRCHAR=0.0155_JWRB ! default value for Charnock

      ! WAM-NEMO COUPLING

      IF (LWNEMOCOU) THEN
        IF (IRANK == 1) WRITE (IU06,*)'CALLING NEMOGCMCOUP_INIT'
#ifdef WITH_NEMO
        CALL NEMOGCMCOUP_INIT( IRANK-1, MPL_COMM, NEMOINIDATE, NEMOINITIME,      &
     &                         NEMOITINI, NEMOITEND, NEMOTSTEP,         &
     &                         .TRUE., -1, .FALSE. )
#endif
        IF (IU06 /= 6 .OR. IRANK == 1) THEN
          WRITE(IU06,*)'NEMO INITIAL DATE         : ',NEMOINIDATE
          WRITE(IU06,*)'NEMO INITIAL TIME         : ',NEMOINITIME
          WRITE(IU06,*)'NEMO INITIAL STEP         : ',NEMOITINI
          WRITE(IU06,*)'NEMO FINAL STEP           : ',NEMOITEND
          WRITE(IU06,*)'NEMO TIME STEP            : ',NEMOTSTEP
        ENDIF
        IF (NEMOFRCO==0) THEN
          WRITE(IU06,*) ' NEMOFRCO NEEDS TO BE SET IN WAM-NEMO MODE.'
          CALL ABORT1
        ELSE
          NEMOCSTEP=NEMOITINI
          NEMONSTEP=NINT(NEMOFRCO*IDELPRO/NEMOTSTEP)
          NEMOWSTEP=0
        ENDIF
        IF ((NEMOFRCO*IDELPRO) /= NINT(NEMOTSTEP*NEMONSTEP)) THEN
          WRITE(IU06,*)'Inconsistent NEMO and WAM coupling intervals:'
          WRITE(IU06,*)'WAM coupling interval is : ',NEMOFRCO*IDELPRO
          WRITE(IU06,*)'NEMO coupling interval is: ',NEMOTSTEP*NEMONSTEP
          CALL ABORT1
        ELSE
          IF (IU06 /= 6 .OR. IRANK == 1) THEN
            WRITE(IU06,*)'NEMONSTEP                 : ',NEMONSTEP
            WRITE(IU06,*)'NEMO coupling interval is : ',                &
     &                   NEMOTSTEP*NEMONSTEP
          ENDIF
        ENDIF
      ENDIF

 20   CONTINUE

      IF (IRANK==1) THEN
        WRITE(6,'(4A)') " WAVEMDL STEPPING FROM DATE ", CBPLTDT, " TO DATE ", CEPLTDT
      ENDIF
      timestep0_start = - wam_user_clock()
      DO WHILE (CDTPRO < CDATEE .OR. CDTPRO == ZERO)
        timestep_start = - wam_user_clock()
        CALL WAVEMDL(CBEGDAT, PSTEP, KSTOP, KSTPW,                      &
     &             NFIELDS, NGPTOTG, NC, NR,                            &
     &             IGRIB_HANDLE_DUM, IGRIB_HANDLE_DUM2,                 &
     &             RMISS, ZRCHAR, FIELDS,                               &
     &             NATMFLX,                                             &
     &             LWCUR, LWSTOKES,                                     &
     &             LLINIT_WVFLDG, NWVFIELDS, WVFLDG,                    &
     &             NLONW, NLATW, LLSTOP, LLWRRE,                        &
     &             LLRESTARTED, ZDELATM,                                &
     &             LDWCOUNORMS, LLNORMWAMOUT_GLOBAL,                    &
     &             MASK_IN, MASK_OUT,                                   &
     &             FRSTIME, NADV, PRPLRADI, PRPLRG,                     &
     &             RNU_ATM, RNUM_ATM,                                   &
     &             IDUM,IDUM, .FALSE.)
        time = (timestep_start+wam_user_clock())*1.e-6
        ISTEP = ISTEP+1
        IF (IRANK==1) THEN
          CALL WAM_MEMINFO(6,ISTEP)
          WRITE(6,'(A,I3,3A,F8.2,A)') " WAVEMDL STEP ",ISTEP," TO DATE ", CDTPRO, &
            & "  TOOK ",time," SECONDS"
        ENDIF
        IF (LLSTOP) EXIT
      ENDDO
      IF (IRANK==1) THEN
        time = (timestep0_start+wam_user_clock())*1.e-6
        WRITE(6,'(A,F8.2,A)') " WAVEMDL STEPS TOOK ",time," SECONDS"
      ENDIF


!     3. CLOSE ODB DATABASE OF SATELLITE DATA
!        ------------------------------------
      IF (LODBRALT) THEN
        CALL WAM_ODB_CLOSE()
      ENDIF

!    4.  FINALIZE RUNWAM
!        ---------------

      CALL GSTATS(0,1)

      CALL WAM_GSTATS_FILE_OPEN(IUGSTATS)
      IF (IUGSTATS == -1 ) IUGSTATS = IU06
      CALL WAM_GSTATS_PRINT(IUGSTATS)
      CALL WAM_GSTATS_FILE_CLOSE(IUGSTATS)

#ifdef WITH_NEMO
      IF (LWNEMOCOU) CALL NEMOGCMCOUP_FINAL
#endif

      CALL EC_MEMINFO(6,"runwam:end",MPL_COMM,KBARR=1,KIOTASK=-1,KCALL=1)

      time=time0+wam_user_clock()
      time=time*1E-06
      WRITE (IU06,'(A)') ' ++++++++++++++++++++++++++++++'
      WRITE (IU06,'(A)') ' + TOTAL USER TIME IN SECONDS +'
      WRITE (IU06,'(A,F18.2,A)') ' + ', time, '         +'
      WRITE (IU06,'(A)') ' + WAVE PROPAGATION TIME      +'
      WRITE (IU06,'(A,F18.2,A)') ' + ', TIME_PROPAG, '         +'
      WRITE (IU06,'(A)') ' + SOURCE TERM TIME           +'
      WRITE (IU06,'(A,F18.2,A)') ' + ', TIME_PHYS, '         +'
      WRITE (IU06,'(A)') ' +                            +'
      WRITE (IU06,'(A,I8,A)') ' + ON PE : ', IRANK, '           +'
      WRITE (IU06,'(A)') ' ++++++++++++++++++++++++++++++'
      IF (IRANK==1 .AND. IU06/=6) THEN
        WRITE (6,'(A)') ' ++++++++++++++++++++++++++++++'
        WRITE (6,'(A)') ' + TOTAL USER TIME IN SECONDS +'
        WRITE (6,'(A,F18.2,A)') ' + ', time, '         +'
        WRITE (6,'(A)') ' + WAVE PROPAGATION TIME      +'
        WRITE (6,'(A,F18.2,A)') ' + ', TIME_PROPAG, '         +'
        WRITE (6,'(A)') ' + SOURCE TERM TIME           +'
        WRITE (6,'(A,F18.2,A)') ' + ', TIME_PHYS, '         +'
        WRITE (6,'(A)') ' +                            +'
        WRITE (6,'(A,I8,A)') ' + ON PE : ', IRANK, '           +'
        WRITE (6,'(A)') ' ++++++++++++++++++++++++++++++'
      ENDIF

      CALL MPCLOSE_UNIT

      IF (LHOOK) CALL DR_HOOK('RUNWAM',1,ZHOOK_HANDLE)

CONTAINS

SUBROUTINE WAM_MEMINFO(KOUT,KSTEP)
USE EC_PARKIND, ONLY : JPIM
IMPLICIT NONE
INTEGER(KIND=JPIM), INTENT(IN) :: KOUT, KSTEP
CHARACTER(LEN=32) CLSTEP
INTEGER(KIND=JPIM) :: ICOMM
#include "ec_meminfo.intfb.h"
WRITE(CLSTEP,'(11X,"STEP",I5," :")') KSTEP
ICOMM = -2 ! No headers from EC_MEMINFO by default
IF (KSTEP == 1) ICOMM = -1 ! Print also headers
CALL EC_MEMINFO(KOUT,TRIM(CLSTEP),ICOMM,KBARR=0,KIOTASK=-1,KCALL=-1)
CALL EC_FLUSH(KOUT)
END SUBROUTINE WAM_MEMINFO

      END SUBROUTINE RUNWAM
