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

MODULE CPLNG_EXCHANGE_MOD

   IMPLICIT NONE

   PRIVATE

   PUBLIC CPLNG_EXCHANGE

CONTAINS

   SUBROUTINE CPLNG_EXCHANGE(YDMCC,YDRIP,YDERAD,YDDYNA,KSTAGE)

      USE PARKIND1, ONLY: JPIM, JPRB
      USE YOMCT2,   ONLY: NSTAR2,NSTOP2
      USE YOMMCC,   ONLY: TMCC
      USE YOMRIP,   ONLY: TRIP
      USE YOERAD,   ONLY: TERAD
      USE YOMDYNA , ONLY : TDYNA
#ifdef WITH_OASIS
      USE MOD_OASIS
#endif
      USE CPLNG_LOG_MOD
      USE CPLNG_DATA_MOD

      ! Argument
      TYPE(TMCC), INTENT(INOUT) :: YDMCC
      TYPE(TRIP), INTENT(INOUT) :: YDRIP
      TYPE(TERAD),INTENT(INOUT) :: YDERAD
      TYPE(TDYNA),INTENT(INOUT) :: YDDYNA

      INTEGER(KIND=JPIM), INTENT(IN) :: KSTAGE

      ! Locals
      INTEGER(KIND=JPIM) :: II
      INTEGER(KIND=JPIM) :: ILVL,ICAT
      INTEGER(KIND=JPIM) :: IINFO
      INTEGER(KIND=JPIM) :: ITIME_IN_SECONDS
      CHARACTER(LEN=3)   :: CLERRSTR

      ! Early return if CPLNG not active
      IF (.NOT.YDMCC%CPLNG_ACTIVE) THEN
         CALL CPLNG_WARNING('CPLNG_EXCHANGE called but CPLNG_ACTIVE is false.')
         RETURN
      ENDIF

      ! Early return if KSTAGE is set to zero (which means ignore this field)
      IF (KSTAGE==0) RETURN

      ! If LPERPET is true, this is a perpetual run and RSTATI (needed later)
      ! doesn't represent the time since the model started. Can't handle this.
      ! Note that LPERPET is true for Aqua planet (LAQUA)
      IF (YDERAD%LPERPET) THEN
         CALL ABOR1("CPLNG_EXCHANGE: Coupling doesn't work for perpetual runs.")
      ENDIF

      ! Compute the current time in seconds since the start of the current leg.
      ! By convention, what we want to send to OASIS is the time at the beginning
      ! of the current time step intervall.
      !
      ! NOTE: Since the OASIS interface requires an integer of
      !       SELECTED_INT_KIND(9) for it's time/date argument in OASIS_Put/Get,
      !       the length of a leg is limited to 10**9 seconds (about 31 years)!
      !
      ! To figure out what the current time of the model is, we rely on RSTATI
      ! from module YOMRIP, which is updated early in the time step by UPDTIM.
      ! However, RSTATI does not contain the time at the beginning of the time
      ! step intervall. Instead, it is at time step n:
      !     IF (LTWOTL): RSTATI == (n+1/2)*dt
      !     ELSE       : ???
      ! and, hence, it is corrected accordingly.
      !   We also need to correct for previous restart legs since RSTATI is not
      ! reset during a restart. This is done by substracting NSTAR2*TSTEP, where
      ! NSTAR2 from YOMCT2 has the first time step number of the leg.
      !
      ! Note that NSTEP, which could also be used to compute the current time, is
      ! updated much later (too late) in the time step loop in CNT4!
      !
      ! Moreover, a test is implemented to avoid calling the coupler after the end
      ! of the leg, which is defined by TSTEP*NSTOP2. IFS is doing one additional
      ! time step, probably due to the two-level time stepping scheme.
      !
      IF (YDDYNA%LTWOTL) THEN
         ! Return if time is beyond last step
         IF (YDRIP%RSTATI>YDRIP%TSTEP*NSTOP2) RETURN
         ! Compute time at start of current step
         ITIME_IN_SECONDS = NINT(YDRIP%RSTATI - NSTAR2*YDRIP%TSTEP - 0.5_JPRB*YDRIP%TSTEP,JPIM)
      ELSE
         CALL ABOR1("CPLNG_EXCHANGE: Can't handle LTWOTL==.FALSE. yet.")
      ENDIF

      DO II=1,YDMCC%CPLNG_NUM_FIELDS

         ! Do nothing if this field is not exchanged in the current stage
         IF (YDMCC%CPLNG_FLD(II)%STAGE/=KSTAGE) CYCLE

         ! Decide whether the coupling field is to be sent (put) or
         ! received (get)
         SELECT CASE (YDMCC%CPLNG_FLD(II)%INOUT)

#ifdef WITH_OASIS

         CASE (OASIS_Out) ! Call OASIS_PUT for couple fields that are sent

            DO ICAT=1,YDMCC%CPLNG_FLD(II)%NUM_CAT
               DO ILVL=1,YDMCC%CPLNG_FLD(II)%NUM_LVL
                  CALL OASIS_PUT(YDMCC%CPLNG_FLD(II)%ID(ILVL,ICAT),  &
                     &           ITIME_IN_SECONDS,             &
                     &           YDMCC%CPLNG_FLD(II)%D(:,ILVL,ICAT), &
                     &           IINFO)

                  SELECT CASE (IINFO)

                  CASE (OASIS_Sent,      &
                     &  OASIS_LocTrans,  &
                     &  OASIS_ToRest,    &
                     &  OASIS_Output,    &
                     &  OASIS_SentOut,   &
                     &  OASIS_ToRestOut, &
                     &  OASIS_Waitgroup, &
                     &  OASIS_Ok         )

                  CASE DEFAULT

                     WRITE (CLERRSTR,'(I3)') IINFO
                     CALL ABOR1("CPLNG_EXCHANGE: Error in OASIS_PUT: "//CLERRSTR)

                  END SELECT
               ENDDO
            ENDDO

         CASE (OASIS_In) ! Call OASIS_GET for couple fields that are received

            DO ICAT=1,YDMCC%CPLNG_FLD(II)%NUM_CAT
               DO ILVL=1,YDMCC%CPLNG_FLD(II)%NUM_LVL
                  CALL OASIS_GET(YDMCC%CPLNG_FLD(II)%ID(ILVL,ICAT),  &
                     &           ITIME_IN_SECONDS,             &
                     &           YDMCC%CPLNG_FLD(II)%D(:,ILVL,ICAT), &
                     &           IINFO)

                  SELECT CASE (IINFO)

                  CASE (OASIS_Recvd,       &
                     &  OASIS_FromRest,    &
                     &  OASIS_Input,       &
                     &  OASIS_RecvOut,     &
                     &  OASIS_FromRestOut, &
                     &  OASIS_Ok           )

                  CASE DEFAULT

                     WRITE (CLERRSTR,'(I3)') IINFO
                     CALL ABOR1("CPLNG_EXCHANGE: Error in OASIS_GET: "//CLERRSTR)

                  END SELECT
               ENDDO
            ENDDO

#endif

         CASE DEFAULT ! Everthing else is an error

            WRITE (CLERRSTR,'(I3,A,I5,A)') II,' (INOUT=',YDMCC%CPLNG_FLD(II)%INOUT,')'
            CALL ABOR1("CPLNG_EXCHANGE: Error in definition of field no. "//CLERRSTR)

         END SELECT

      ENDDO

   END SUBROUTINE CPLNG_EXCHANGE

END MODULE CPLNG_EXCHANGE_MOD
