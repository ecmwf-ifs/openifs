! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

! ----------------------------------------------------------------------

      SUBROUTINE GRSTNAME (CDATED, CDATEF, IFCST, FILEID, KCPLEN, CPAD,  &
     &                     FILENAME)

! ----------------------------------------------------------------------

!**** *GRSTNAME*  DETERMINES FILENAME (INCLUDING PATH) FOR RESTART FILES

!     J. BIDLOT     ECMWF         NOVEMBER 1996      PARALLEL RUN 
!     B. HANSEN     ECMWF         APRIL    1996      MOVE CYMPPATH TO
!                                                    COMTEXT.

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *GRSTNAME (CDATED, CDATEF, IFCST, FILEID, KCPLEN, CPAD,
!                         FILENAME)*

!         *CDATED*   CHAR*14    DATE (YYYYMMDDHHmmss) USED TO SPECIFY THE FILE NAME (see below).
!         *CDATEF*   CHAR*14    DATE TIME GROUP OF FORECAST START.
!         *IFCST*               FORECAST STEP IN HOURS. (only used if CDATED < CDATEF)
!         *FILEID*   CHARACTER  FILE ID.
!         *KCPLEN*   INTEGER    LENGTH OF CPAD
!         *CPAD*     CHARACTER  PATH FOR OUTPUT TO DISK  
!         *FILENAME* CHARACTER  FILE NAME INCLUDING PATH

!     EXTERNALS.
!     ----------

!       *DIFDATE*   - COMPUTE TIME DIFFERENCES.

!     METHOD.
!     -------

!        FILE NAME IS BUILT FROM THE 3 CHARACTER FILEID FOLLOWED
!        BY YYYYMMDDHHmmss_ddddddhhmmss WHERE YYYYMMDDHHmmss IS THE LAST ANALYSIS TIME (See CDATEH)
!        AND ddddddhhmmss IS THE FORCECAST RANGE GIVEN AS: 
!        dddddd DAYS, hh HOURS, mi MINUTES AND ss seconds.
!        ddddddhhmmss=000000000000  FOR ANALYSIS FILES.
!        THESE NUMBERS ARE CONSTRUCTED FROM CDATED AND CDATEF.

!        THEN IF KCPLEN .NOT.EQ.0
!        THE PATH NAME IS ADDED. THE PATH IS DETERMINED CPAD
!        IN THE INPUT NAMELIST 

!     REFERENCES.
!     -----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      CHARACTER(LEN=14), INTENT(IN) :: CDATED, CDATEF
      INTEGER(KIND=JWIM), INTENT(IN) :: IFCST
      CHARACTER(LEN=3), INTENT(IN) :: FILEID
      INTEGER(KIND=JWIM), INTENT(IN) :: KCPLEN
      CHARACTER(LEN=*), INTENT(IN) :: CPAD
      CHARACTER(LEN=*), INTENT(OUT) :: FILENAME


      INTEGER(KIND=JWIM) :: IC, IREP, NDAY, ISHIFT, ISHIFTDAY, IYRDIF
      INTEGER(KIND=JWIM) :: IYEARSEC, IMAXYEAR
      INTEGER(KIND=JWIM) :: NCHUNCK, ISFT, LNAME
      INTEGER(KIND=JWIM) :: IYEAR1, IMON1, IDAY1, IHOUR1, IMIN1, ISEC1
      INTEGER(KIND=JWIM) :: IYEAR2, IMON2, IDAY2, IHOUR2, IMIN2, ISEC2
      INTEGER(KIND=JWIM) :: IDDDDDD, IHH, IMIMI, ISS

      CHARACTER(LEN=14) :: CDATEH, CDT
      CHARACTER(LEN=30) :: FILENA

!*    1.0  CONSTRUCT FILE NAME ALONE 
!          -------------------------

      IREP=KIND(IMAXYEAR)
      NDAY=366
      IYEARSEC=NDAY*86400
!      IMAXYEAR=(2**(8*IREP-1)-1)/IYEARSEC  ! integer overflow or invalid flp in 2**x -- depending on compiler impl. (SS/04-May-2018)
      IMAXYEAR=HUGE(IMAXYEAR)/IYEARSEC ! desired coding

      FILENA = FILEID
      IF (CDATED < CDATEF) THEN
        CDATEH = CDATED
        ISHIFT = IFCST 
        ISHIFTDAY=0
      ELSE
        CDATEH = CDATEF

        READ(CDATEF,'(I4, 5I2)')IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1,ISEC1
        READ(CDATED,'(I4, 5I2)')IYEAR2,IMON2,IDAY2,IHOUR2,IMIN2,ISEC2
        ISHIFTDAY=0
        CDT=CDATED
!       IN CASE THE DATE DIFFERENCE IS LARGER THAN WHAT DIFDATE CAN DO.
        IYRDIF=IYEAR2-IYEAR1
        IF (ABS(IYRDIF) >= IMAXYEAR) THEN
          NCHUNCK=ABS(IYRDIF)/IMAXYEAR
          DO IC=1,NCHUNCK
            ISFT=-SIGN(1,IYRDIF)*IMAXYEAR*IYEARSEC
            CALL INCDATE(CDT,ISFT)
            ISHIFTDAY=ISHIFTDAY-ISFT/86400
          ENDDO
        ENDIF

        CALL DIFDATE (CDATEF, CDT, ISHIFT)
      ENDIF
      FILENA(4:17) = CDATEH
      FILENA(18:18) = "_"

      IDDDDDD=ISHIFT/86400
      IHH=(ISHIFT-IDDDDDD*86400)/3600
      IMIMI=(ISHIFT-IDDDDDD*86400-IHH*3600)/60
      ISS=ISHIFT-IDDDDDD*86400-IHH*3600-IMIMI*60

      IDDDDDD=IDDDDDD+ISHIFTDAY

      WRITE (FILENA(19:30),'(I6.6,3I2.2)')IDDDDDD,IHH,IMIMI,ISS
      LNAME = LEN_TRIM(FILENA)

!     2.0   GET PATH AND ADD IT TO FILE NAME
!           --------------------------------

      IF (KCPLEN == 0) THEN
        FILENAME=FILENA(1:LNAME)
      ELSE
        FILENAME=CPAD(1:KCPLEN)//'/'//FILENA(1:LNAME)
      ENDIF

      END SUBROUTINE GRSTNAME
