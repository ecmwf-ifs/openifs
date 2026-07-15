! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

! ======================================================================

      SUBROUTINE DIFDATE (CDATE1, CDATE2, KSHIFT)

! ======================================================================

!****     *DIFDATE* - TO COMPUTE TIME DIFFERENCE.

!     J. BIDLOT     FEB 2007  RE-INRODUCING THE OLD DIFDATE WITHOUT 
!                             THE NEED FOR ECLIB. ADDING SECONDS

!**   PURPOSE.
!     --------

!          COMPUTE THE SECONDS BETWEEN THE INPUT DATES.

!**   INTERFACE.
!     ----------

!         *CALL* *DIFDATE (CDATE1, CDATE2, KSHIFT)*

!     I/   *CDATE1* CHAR*12 - DATE TIME GROUP (YYYYMMDDHHMM) OR
!                   CHAR*14 - DATE TIME GROUP (YYYYMMDDHHMMSS).
!     I/   *CDATE2* CHAR*12 - DATE TIME GROUP (YYYYMMDDHHMM) OR
!                   CHAR*14 - DATE TIME GROUP (YYYYMMDDHHMMSS).
!      /O  *KSHIFT* INTEGER - DIFFERENCE IN SECONDS (CDATE2-CDATE1).
!                             RESTRICTION:
!                             SINCE KSHIFT IS AN INTEGER
!                            -MAXINT < KSHIFT <  +MAXINT
!                             WHERE MAXINT=2**(IREP-1)
!                             WHERE IREP IS THE INTEGER REPRESENTATION
!                             USUALLY 32 BITS --> MAXINT=2147482648,
!                             ABOUT 68 YEARS!


!     METHOD.
!     -------

!       THE DATE AND TIME CAN BE SUPPLIED IN THE Y2K COMPLIANT
!       12 or 14 CHARACTER FORMAT. IF THE 10 CHARACTER FORMAT IS
!       USED AN ERROR MESSAGE WILL BE ISSUED. 


!     EXTERNALS.
!     ----------

!!!!!!!!!!!!!!!!!!!!!!!!!!
!      INCDATE 
!!!!!!!!!!!!!!!!!!!!!!!!!!


!     REFERENCES.
!     -----------


!     MODIFICATIONS
!     -----------

!       NONE.

!----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWABORT, ONLY : WAM_ABORT
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!----------------------------------------------------------------------
      IMPLICIT NONE
#include "abort1.intfb.h"
#include "incdate.intfb.h"

      INTEGER(KIND=JWIM) :: M, IL, ISI, IDUM, IRET, ISHIFT, KSHIFT
      INTEGER(KIND=JWIM) :: IYEAR ,IMON ,IDAY ,IHOUR ,IMIN ,ISEC
      INTEGER(KIND=JWIM) :: IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1,ISEC1
      INTEGER(KIND=JWIM) :: IYEAR2,IMON2,IDAY2,IHOUR2,IMIN2,ISEC2
      INTEGER(KIND=JWIM) :: MON(12)

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      CHARACTER(LEN=*) :: CDATE1, CDATE2
      CHARACTER(LEN=14) :: CDT1, CDT2
      CHARACTER(LEN=80) :: CLFMT

      LOGICAL :: LLND

      DATA MON /31,28,31,30,31,30,31,31,30,31,30,31/


! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('DIFDATE',0,ZHOOK_HANDLE)

!*    1.0 SPLIT DATE TIME GROUPS INTO SECOND, MINUTE, HOUR, DAY, MONTH,
!*        YEAR DEPENDENT ON THE FORMAT OF CDATE.
!         -----------------------------------------------------------

!*        SPLIT CDATE1.

      IL = LEN_TRIM(CDATE1)
      IF (IL==10) THEN
        WRITE(CLFMT,'(A49, I2.2, A2)')                                  &
     &     '(" DIFDATE : NON-Y2K COMPLIANT DATE1 IN USE ", A', IL, ' )'
        WRITE(6, CLFMT ) CDATE1
        WRITE(*, CLFMT ) CDATE1
        LLND = .FALSE.
        CALL WAM_ABORT(__FILENAME__,__LINE__)
      ELSEIF (IL==12) THEN
        LLND = .TRUE.
        READ(CDATE1,'(I4, 5I2)')IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1
        ISEC1 = 0
        CDT1=CDATE1(1:12)//'00'
      ELSEIF (IL==14) THEN
        LLND = .TRUE.
        READ(CDATE1,'(I4, 5I2)')IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1,ISEC1
!       with satellite data it might happen that the seconds were
!       rounded to 60. One needs to increment the date properly.
        IF (ISEC1==60) THEN
          IYEAR=IYEAR1
          IMON=IMON1
          IDAY=IDAY1
          IHOUR=IHOUR1
          IMIN=IMIN1
          ISEC=ISEC1-1
          ISHIFT=1
          WRITE(CDT1,'(I4.4, 5I2.2)')IYEAR, IMON, IDAY, IHOUR,          &
     &                               IMIN, ISEC
          CALL INCDATE(CDT1, ISHIFT)
          READ(CDT1,'(I4, 5I2)')IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1,ISEC1
        ELSE
          CDT1=CDATE1
        ENDIF
      ELSEIF (IL==0) THEN
        WRITE(6, *) ' '
        WRITE(6, *) ' DIFDATE:  FATAL@! DATE1 IS 0 CHARACTER LONG'
        WRITE(6, *) ' '
        CALL WAM_ABORT(__FILENAME__,__LINE__)
      ELSE
        WRITE(6, *) ' '
        WRITE(CLFMT,'(A61, I2.2, A2)')                                  &
     &  '(" DIFDATE :  FATAL@! DATE1 IS ",I2," CHARACTERS LONG!!= ",A', &
     &    IL, ' )'
        WRITE(6, CLFMT) IL, CDATE1
        WRITE(6, *) ' '
        CALL WAM_ABORT(__FILENAME__,__LINE__)
      ENDIF

!*        SPLIT CDATE2.

      IL = LEN_TRIM(CDATE2)
      IF (IL==10) THEN
        WRITE(CLFMT,'(A49, I2.2, A2)')                                  &
     &     '(" DIFDATE : NON-Y2K COMPLIANT DATE2 IN USE ", A', IL, ' )'
        WRITE(6, CLFMT ) CDATE2
        WRITE(*, CLFMT ) CDATE2
        LLND = .FALSE.
        CALL ABORT1
      ELSEIF (IL==12) THEN
        LLND = .TRUE.
        READ(CDATE2,'(I4, 5I2)')IYEAR2,IMON2,IDAY2,IHOUR2,IMIN2
        ISEC2 = 0
        CDT2=CDATE2(1:12)//'00'
      ELSEIF (IL==14) THEN
        LLND = .TRUE.
        READ(CDATE2,'(I4, 5I2)')IYEAR2,IMON2,IDAY2,IHOUR2,IMIN2,ISEC2
!       with satellite data it might happen that the seconds were
!       rounded to 60. One needs to increment the date properly.
        IF (ISEC2==60) THEN
          IYEAR=IYEAR2
          IMON=IMON2
          IDAY=IDAY2
          IHOUR=IHOUR2
          IMIN=IMIN2
          ISEC=ISEC2-1
          ISHIFT=1
          WRITE(CDT2,'(I4.4, 5I2.2)')IYEAR, IMON, IDAY, IHOUR,          &
     &                               IMIN, ISEC
          CALL INCDATE(CDT2, ISHIFT)
          READ(CDT2,'(I4, 5I2)')IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1,ISEC1
        ELSE
          CDT2=CDATE2
        ENDIF
      ELSEIF (IL==0) THEN
        WRITE(6, *) ' '
        WRITE(6, *) ' DIFDATE:  FATAL@! DATE2 IS 0 CHARACTER LONG'
        WRITE(6, *) ' '
        CALL WAM_ABORT(__FILENAME__,__LINE__)
      ELSE
        WRITE(6, *) ' '
        WRITE(CLFMT,'(A61, I2.2, A2)')                                  &
     &  '(" DIFDATE :  FATAL@! DATE2 IS ",I2," CHARACTERS LONG!!= ",A', &
     &    IL, ' )'
        WRITE(6, CLFMT) IL, CDATE2
        WRITE(6, *) ' '
        CALL WAM_ABORT(__FILENAME__,__LINE__)
      ENDIF

      IRET=0
      MON(2)=MFEB_LENGTH(IYEAR1)
      IF (IMON1<1 .OR. IMON1>12) IRET=-1
      IF (IDAY1<1 .OR. IDAY1>MON(MIN(MAX(IMON1,1),12))) IRET=-2
      IF (IHOUR1<0 .OR. IHOUR1>23) IRET=-3
      IF (IMIN1<0 .OR. IMIN1>59) IRET=-4
      IF (ISEC1<0 .OR. ISEC1>59) IRET=-5

      MON(2)=MFEB_LENGTH(IYEAR2)
      IF (IMON2<1 .OR. IMON2>12) IRET=-1
      IF (IDAY2<1 .OR. IDAY2>MON(MIN(MAX(IMON2,1),12))) IRET=-2
      IF (IHOUR2<0 .OR. IHOUR2>23) IRET=-3
      IF (IMIN2<0 .OR. IMIN2>59) IRET=-4
      IF (ISEC2<0 .OR. ISEC2>59) IRET=-5

      IF (IRET/=0) THEN
        WRITE(6, '(" DIFDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(6, '("           INPUT DATE(S) INCORRECT,")')
        WRITE(6, '("           IYEAR1 ", I4 )') IYEAR1 
        WRITE(6, '("           IMON1 ", I2 )') IMON1 
        WRITE(6, '("           IDAY1 ", I2 )') IDAY1 
        WRITE(6, '("           IHOUR1 ", I2 )') IHOUR1 
        WRITE(6, '("           IMIN1 ", I2 )') IMIN1 
        WRITE(6, '("           ISEC1 ", I2 )') ISEC1 
        WRITE(6, '("           IYEAR2 ", I4 )') IYEAR2 
        WRITE(6, '("           IMON2 ", I2 )') IMON2 
        WRITE(6, '("           IDAY2 ", I2 )') IDAY2 
        WRITE(6, '("           IHOUR2 ", I2 )') IHOUR2 
        WRITE(6, '("           IMIN2 ", I2 )') IMIN2 
        WRITE(6, '("           ISEC2 ", I2 )') ISEC2 
        WRITE(*, '(" DIFDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(*, '("           INPUT DATE(S) INCORRECT,")')
        WRITE(*, '("           IYEAR1 ", I4 )') IYEAR1 
        WRITE(*, '("           IMON1 ", I2 )') IMON1 
        WRITE(*, '("           IDAY1 ", I2 )') IDAY1 
        WRITE(*, '("           IHOUR1 ", I2 )') IHOUR1 
        WRITE(*, '("           IMIN1 ", I2 )') IMIN1 
        WRITE(*, '("           ISEC1 ", I2 )') ISEC1 
        WRITE(*, '("           IYEAR2 ", I4 )') IYEAR2 
        WRITE(*, '("           IMON2 ", I2 )') IMON2 
        WRITE(*, '("           IDAY2 ", I2 )') IDAY2 
        WRITE(*, '("           IHOUR2 ", I2 )') IHOUR2 
        WRITE(*, '("           IMIN2 ", I2 )') IMIN2 
        WRITE(*, '("           ISEC2 ", I2 )') ISEC2 
        CALL ABORT1
      ENDIF

!*    CHANGE DATE TIME GROUPS TO ENSURE THAT THE SECOND IS LARGER.
!     ------------------------------------------------------------

      ISI = 1 
      IF (CDT1 > CDT2) THEN
         ISI  = -1
         IDUM=IYEAR1
         IYEAR1=IYEAR2
         IYEAR2=IDUM
         IDUM=IMON1
         IMON1=IMON2
         IMON2=IDUM
         IDUM=IDAY1
         IDAY1=IDAY2
         IDAY2=IDUM
         IDUM=IHOUR1
         IHOUR1=IHOUR2
         IHOUR2=IDUM
         IDUM=IMIN1
         IMIN1=IMIN2
         IMIN2=IDUM
         IDUM=ISEC1
         ISEC1=ISEC2
         ISEC2=IDUM
      ENDIF

! ----------------------------------------------------------------------

!*    2.0 CORRECT DAYS IN FEBRUARY OF FIRST YEAR FOR LEAP-YEAR.
!         -----------------------------------------------------
 
      MON(2)=MFEB_LENGTH(IYEAR1)

! ----------------------------------------------------------------------
 
!*    3.0 COMPUTE TIME DIFFERENCE IN SECONDS.
!         -----------------------------------
 
!     3.1 DIFFERENCE BETWEEN DAY, HOUR ,MINUTE, SECOND.
 
      KSHIFT = (((IDAY2-IDAY1)*24+IHOUR2-IHOUR1)*60+IMIN2-IMIN1)*60+    &
     &         ISEC2-ISEC1
 
!     3.2 ADD DIFFERENCE FROM MONTH.
 
      IF (IYEAR2 > IYEAR1) THEN
 
!     3.2.1 START AND END MONTHS ARE IN DIFFERENT YEARS.
 
         DO WHILE (IYEAR2 > IYEAR1)
            DO M=IMON1,12
               KSHIFT = KSHIFT + INT(MON(M)*86400)
            ENDDO
            IMON1 = 1
            IYEAR1 = IYEAR1 + 1
            MON(2)=MFEB_LENGTH(IYEAR1)

         ENDDO
         IF (IMON2 > 1) THEN
            MON(2)=MFEB_LENGTH(IYEAR2)
            DO M=1,IMON2-1
               KSHIFT = KSHIFT + INT(MON(M)*86400)
            ENDDO
         ENDIF
      ELSE
 
!     3.2.2 START AND END MONTHS ARE IN THE SAME YEAR.
 
         IF (IMON2 > IMON1) THEN
            DO M=IMON1,IMON2-1
               KSHIFT = KSHIFT + INT(MON(M)*86400)
            ENDDO
         ENDIF
      ENDIF
 
! ----------------------------------------------------------------------
 
!*    4.0 CHANGE SIGN OF DIFFERENCE IF DATES HAVE BEEN EXCHANGED
!*        AND CONVERT TO SECONDS.
!         ------------------------------------------------------
 
      KSHIFT = KSHIFT*ISI

      IF (LHOOK) CALL DR_HOOK('DIFDATE',1,ZHOOK_HANDLE)

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

      END SUBROUTINE DIFDATE
