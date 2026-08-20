! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

! ======================================================================

      SUBROUTINE INCDATE (CDATE, ISHIFT)

! ======================================================================

!**** *INCDATE* - TO UPDATE DATE TIME GROUP

!     J. BIDLOT   FEB 2007    RE-INRODUCING THE OLD WAY WITHOUT 
!                             THE NEED FOR ECLIB. ADDING SECONDS

!**   PURPOSE.
!     --------
!       UPDATING DATE TIME GROUP.

!**   INTERFACE.
!     ----------

!       *CALL* *INCDATE (CDATE,ISHIFT)*
!         *CDATE*  CHAR*12 - DATE TIME GROUP (YYYYMMDDHHMM) OR
!                  CHAR*14 - DATE TIME GROUP (YYYYMMDDHHMMSS)
!         *ISHIFT* INTEGER - TIME INCREMENT IN SECONDS, RESTRICTIONS:
!                            SINCE ISHIFT IS AN INTEGER
!                            -MAXINT < ISHIFT <  +MAXINT
!                             WHERE MAXINT=2**(IREP-1)
!                             WHERE IREPR IS THE INTEGER REPRESENTATION
!                             USUALLY 32 BITS --> MAXINT=2147482648
!                             (ABOUT 68 YEARS).

!     METHOD.
!     -------

!       THE DATE AND TIME CAN BE SUPPLIED IN THE Y2K COMPLIANT
!       14 CHARACTER OR 12 CHARACTER FORMAT. IF THE 10 CHARACTER
!       FORMAT IS USED, AN ERROR MESSAGE WILL BE ISSUED. 

!     EXTERNALS.
!     ----------

!     REFERENCES.
!     -----------


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: ISHIFT
      CHARACTER(LEN=*), INTENT(INOUT) :: CDATE

      INTEGER(KIND=JWIM) :: IL, IRET
      INTEGER(KIND=JWIM) :: IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC
      INTEGER(KIND=JWIM) :: MON(12)

      CHARACTER(LEN=80) :: CLFMT

      LOGICAL :: LLND

      DATA MON /31,28,31,30,31,30,31,31,30,31,30,31/

! ----------------------------------------------------------------------

!*    1.0 SPLIT DATE TIME GROUP INTO SECOND, MINUTE, HOUR, DAY, MONTH,
!*        YEAR DEPENDENT ON THE FORMAT OF CDATE.
!         -----------------------------------------------------------

      IL = LEN_TRIM(CDATE)
      IF (IL==10) THEN
        WRITE(CLFMT,'(A49, I2.2, A2)')                                  &
     &     '(" INCDATE : NON-Y2K COMPLIANT  DATE IN USE ", A', IL, ' )'
        WRITE(6, CLFMT ) CDATE
        WRITE(*, CLFMT ) CDATE
        LLND = .FALSE.
        CALL ABORT1
      ELSEIF (IL==12) THEN
        LLND = .TRUE.
        READ(CDATE,'(i4, 4i2)')IYEAR, IMON, IDAY, IHOUR, IMIN
        ISEC = 0
      ELSEIF (IL==14) THEN
        LLND = .TRUE.
        READ(CDATE,'(i4, 5i2)')IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC
      ELSEIF (IL==0) THEN
        WRITE(6, *) ' INCDATE:  FATAL@! DATE IS 0 CHARACTER LONG'
        WRITE(*, *) ' INCDATE:  FATAL@! DATE IS 0 CHARACTER LONG'
        CALL ABORT1
      ELSE
        WRITE(CLFMT,'(A59, I2.2, A2)')                                  &
     &  '(" INCDATE:  FATAL@! DATE IS ",I2," CHARACTERS LONG!!= ",A',   &
     &    IL, ' )'
        WRITE(6, CLFMT) IL, CDATE
        WRITE(*, CLFMT) IL, CDATE
        STOP
      ENDIF

      IRET=0
      MON(2)=MFEB_LENGTH(IYEAR)
      IF (IMON<1 .OR. IMON>12) IRET=-1
      IF (IDAY<1 .OR. IDAY>MON(MIN(MAX(IMON,1),12))) IRET=-2
      IF (IHOUR<0 .OR. IHOUR>23) IRET=-3
      IF (IMIN<0 .OR. IMIN>59) IRET=-4
      IF (ISEC<0 .OR. ISEC>59) IRET=-5
      IF (IRET/=0) THEN
        WRITE(6, '(" INCDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(6, '("           INPUT DATE INCORRECT,")')
        WRITE(6, '("           IYEAR ", I4 )') IYEAR
        WRITE(6, '("           IMON ", I2 )') IMON
        WRITE(6, '("           IDAY ", I2 )') IDAY
        WRITE(6, '("           IHOUR ", I2 )') IHOUR
        WRITE(6, '("           IMIN ", I2 )') IMIN
        WRITE(6, '("           ISEC ", I2 )') ISEC
        WRITE(*, '(" INCDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(*, '("           INPUT DATE INCORRECT,")')
        WRITE(*, '("           IYEAR ", I4 )') IYEAR
        WRITE(*, '("           IMON ", I2 )') IMON
        WRITE(*, '("           IDAY ", I2 )') IDAY
        WRITE(*, '("           IHOUR ", I2 )') IHOUR
        WRITE(*, '("           IMIN ", I2 )') IMIN
        WRITE(*, '("           ISEC ", I2 )') ISEC
        CALL ABORT1
      ENDIF
      IF (IL==12 .AND. MOD(ISHIFT,60) /= 0) THEN
        IRET=-6
        WRITE(6, '(" INCDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(6, '(" WHEN HANDELING CHARACTER*12 DATES,")')
        WRITE(6, '(" ISHIFT MUST BE A MULTIPLE OF 60 SECONDS")')
        WRITE(6, '(" ISHIFT = ", I12 )') ISHIFT
        WRITE(*, '(" INCDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(*, '(" WHEN HANDELING CHARACTER*12 DATES,")')
        WRITE(*, '(" ISHIFT MUST BE A MULTIPLE OF 60 SECONDS")')
        WRITE(*, '(" ISHIFT = ", I12 )') ISHIFT
        CALL ABORT1
      ENDIF

! ----------------------------------------------------------------------

!*    2.0 INCREMENT OR DECREMENT THE DATE AND TIME BY ISHIFT SECONDS.
!         -----------------------------------------------------------


      ISEC=ISEC+MOD(ISHIFT,60)
      IF (ISEC >= 60) THEN
        IMIN = IMIN+ISEC/60
        ISEC = ISEC-(ISEC/60)*60
      ELSE IF (ISEC < 0) THEN
        IMIN = IMIN +(ISEC-59)/60
        ISEC = ISEC-((ISEC-59)/60)*60
      END IF

      IMIN=IMIN+ISHIFT/60

!     2.1 POSITIVE SHIFT GREATER THAN 1 MINUTE.

      IF (IMIN >= 60) THEN
         IHOUR = IHOUR + IMIN/60
         IMIN = IMIN - (IMIN/60)*60
         IF (IHOUR >= 24) THEN
            IDAY = IDAY + IHOUR/24
            IHOUR = IHOUR - (IHOUR/24)*24
            DO WHILE (IDAY > MON(IMON))
               IDAY=IDAY-MON(IMON)
               IMON=IMON+1
               IF (IMON == 13) THEN
                  IMON = 1
                  IYEAR=IYEAR+1
                  MON(2)=MFEB_LENGTH(IYEAR)
               END IF
            ENDDO
         END IF
      ELSE IF (IMIN < 0) THEN
 
!     2.2 NEGATIVE SHIFT.

         IHOUR = IHOUR + (IMIN-59)/60
         IMIN = IMIN - ((IMIN-59)/60)*60
         IF (IHOUR < 0) THEN
            IDAY = IDAY + (IHOUR-23)/24
            IHOUR = IHOUR - ((IHOUR-23)/24)*24
            DO WHILE (IDAY < 1)
               IMON=IMON-1
               IF (IMON == 0) THEN
                  IMON = 12
                  IYEAR=IYEAR-1
                  MON(2)=MFEB_LENGTH(IYEAR)
               END IF
               IDAY=IDAY+MON(IMON)
            ENDDO
         END IF
      END IF


! ----------------------------------------------------------------------

!*    3.0 CREATE OUTPUT STRING OF NEW DATE AND TIME USING THE INPUT
!*        FORMAT.
!         ---------------------------------------------------------

      IRET=0
      MON(2)=MFEB_LENGTH(IYEAR)
      IF (IMON<1 .OR. IMON>12) IRET=-1
      IF (IDAY<1 .OR. IDAY>MON(MIN(MAX(IMON,1),12))) IRET=-2
      IF (IHOUR<0 .OR. IHOUR>23) IRET=-3
      IF (IMIN<0 .OR. IMIN>59) IRET=-4
      IF (ISEC<0 .OR. ISEC>59) IRET=-5
      IF (IRET/=0) THEN
        WRITE(6, '(" INCDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(6, '("           OUTPUT DATE INCORRECT,")')
        WRITE(6, '("           IYEAR ", I4 )') IYEAR
        WRITE(6, '("           IMON ", I2 )') IMON
        WRITE(6, '("           IDAY ", I2 )') IDAY
        WRITE(6, '("           IHOUR ", I2 )') IHOUR
        WRITE(6, '("           IMIN ", I2 )') IMIN
        WRITE(6, '("           ISEC ", I2 )') ISEC
        WRITE(6, '("           ISHIFT ", I12 )') ISHIFT
        CALL ABORT1
      ENDIF

      IF (LLND) THEN
        IF (IL==12) THEN
          WRITE(CDATE,'(I4.4, 4I2.2)')IYEAR, IMON, IDAY, IHOUR,         &
     &                                IMIN
        ELSEIF (IL==14) THEN
          WRITE(CDATE,'(I4.4, 5I2.2)')IYEAR, IMON, IDAY, IHOUR,         &
     &                                IMIN, ISEC
        ELSE
          WRITE(6,'(" THE LENGTH OF THE INPUT CHARACTER CHANGED @!")')
          CALL ABORT1
        ENDIF
      ELSE
        IYEAR = IYEAR - 1900
        IF (IYEAR >= 100) THEN
          WRITE(CLFMT,'(A39, I2.2, A2)')                                &
     &     '(" NON-Y2K COMPLIANT  DATE IN USE ", A', IL, ' )'
          WRITE(6, CLFMT ) CDATE
          WRITE(6,'(" THIS CANNOT BE DONE ANY LONGER")')
          CALL ABORT1
        ENDIF
        WRITE(CDATE,'(5I2.2)')IYEAR, IMON, IDAY, IHOUR, IMIN
      ENDIF

      RETURN

      CONTAINS

      INTEGER(KIND=JWIM) FUNCTION MFEB_LENGTH(IYEAR)
!     LENGTH OF FEBRUARY IN DAYS FOR YEAR IYEAR (YYYY)
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      IMPLICIT NONE
      INTEGER(KIND=JWIM), INTENT(IN) :: IYEAR
      IF (MOD(IYEAR,400) == 0) THEN
        MFEB_LENGTH=29
      ELSEIF (MOD(IYEAR,100) == 0) THEN
        MFEB_LENGTH=28
      ELSEIF (MOD(IYEAR,4) == 0) THEN
        MFEB_LENGTH=29
      ELSE
        MFEB_LENGTH=28
      ENDIF
      END FUNCTION MFEB_LENGTH

      END SUBROUTINE INCDATE
