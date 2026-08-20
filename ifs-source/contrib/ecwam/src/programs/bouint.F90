! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      PROGRAM BOUINT

! ----------------------------------------------------------------------
 
!**** *BOUINT* -  INTERPOLATION OF BOUNDARY VALUE SPECTRA IN TIME.
 
!     SUSANNE HASSELMANN JUNE 1990.
!     H. GUNTHER        GKSS/ECMWF   JAN. 91     MODIFIED FOR CYCLE_4
 
!*    PURPOSE.
!     --------
 
!       INTERPOLATION OF BOUNDARY SPECTRA OF THE COARSE GRID
!       WAVE MODEL OUTPUT IN TIME.
 
!**   INTERFACE.
!     ----------
 
!        *IU01*  INTEGER  INPUT UNIT OF SPECTRA.
!        *IU02*  INTEGER  OUTPUT UNIT OF INTERPOLATED SPECTRA.
!        *IU05*  INTEGER  INPUT UNIT FOR USER INPUT.
!        *IU06*  INTEGER  PRINTER OUTPUT UNIT.
 
!     METHOD.
!     -------
 
!       THE SPECTRA ARE READ FROM THE INPUT UNIT AND INTERPOLATED TO
!       THE OUTPUT TIMES SPECIFIED IN THE USER INPUT. IF THE FIRST
!       DATE OF THE INPUT SPECTRA IS LATER THAN THE FIRST OUTPUT DATE
!       THE FIRST INPUT SPECTRUM IS KEPT FOR ALL DATES BEFORE THE
!       THE FIRST INPUT DATE.
!       THE INTERPOLATION IS DONE IN FOUR STEPS..
!       - ROTATE SPECTRA ACCORDING TO MEAN OF MEAN ANGLES,
!       - TRANSFORM FREQUENCIES ACCORDING TO MEAN OF MEAN FREQUENCIES,
!       - ADJUST ENERGY ACCORDCING TO MEAN OF TOTAL ENERGY AND
!       - INTERPOLATE RESULTING SPECTRA.
 
!       INPUT AND OUTPUT DATA FILES ARE AUTOMATICALLY ASSIGN.
!       FILE NAMES AND PATHES ARE DEFINED BY USER INPUT.
!       THE USER MAY HAVE TO
!       ADOPT THIS SUBS FOR HIS COMPUTER ENVIROMENT.
 
!     EXTERNALS.
!     ----------
 
!       *ABORT*     - TERMINATES PROCESSING.
!       *DIFDATE*   - COMPUTE TIME DIFFERENCE.
!       *GSFILE*    - OPENS A FILE.
!       *INCDATE*   - INCREMENT A DATE.
!       *INTSPEC*   - INTERPOLATE A SPECTRUM.
!       *ROTSPEC*   - ROTATES SPECTRUM.
!       *SFILE*     - SAVES A FILE (FRONT/END COMPUTER ONLY).
!       *STRSPEC*   - TRANSFORM FREQUENCIES.
!       *UIBOU*     - READ USER INPUT.
 
!        FLOWCHART OF BOUINT. CALLS OD SUB ABORT ARE NOT SHOWN.
!        ABORT IS CALLED BY BOUINT, UIBOU, GSFILE.
!         _________        _________
!        |  BOUINT |______|  UIBOU  |
!        |_________|  |   |_________|
!                     |    _________
!                     |___| DIFDATE |
!                     |   |_________|
!                     |    _________
!                     |___| GSFILE  |
!                     |   |_________|
!                     |    _________
!                     |___| INCDATE |   ___| ROTSPEC |
!                     |   |_________|  |   |_________|
!                     |    _________   |    _________
!                     |___| INTSPEC |__|___| STRSPEC |
!                         |_________|      |_________|
 
!      REFERENCES.
!      -----------
 
!        K.HASSELMANN, 1990,
!           INTERPOLATION OF WAVE SPECTRA. WAM NOTE 6/6/90.
 
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWTEXT  , ONLY : PATH

 
! ----------------------------------------------------------------------
      IMPLICIT NONE
 
      INTEGER(KIND=JWIM) :: K, M, IJ 
      INTEGER(KIND=JWIM) :: NANG, NFRE, MBMAX
      INTEGER(KIND=JWIM) :: KL, ML, NBOUNC, IDELPRC
      INTEGER(KIND=JWIM) :: IU01, IU02, IU05, IU06
      INTEGER(KIND=JWIM) :: IDELPRF, IDELFI, NEWOUT
      INTEGER(KIND=JWIM) :: IDEL12, IDEL1L

      REAL(KIND=JWRB) :: FMEAN_PT, EMEAN_PT, THQ_PT
      REAL(KIND=JWRB) :: XANG, XFRE, TH0, FR1, CO, XBOU, XDELC, XDELF
      REAL(KIND=JWRB) :: DEL12, DEL1L
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:) :: XLAT, XLON
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:) :: FR
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:) :: FMEAN, EMEAN, THQ 
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:) :: FL
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:,:,:) :: F

      CHARACTER(LEN=14) :: CDATE1, CDATE2, CDTOUT, CDATEA, CDATEE,      &
     &          CDATES, CDATE
      CHARACTER(LEN=3) :: USERID, RUNDI, FILEDI, RUNDO, FILEDO
      CHARACTER(LEN=60) :: PATHI, PATHO

      LOGICAL :: LREAL, LSWAN
 
!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *USERID*    CHARACTER USER IDENTIFIER.
!      *RUNDI*     CHARACTER RUN  IDENTIFIER OF INPUT DATA.
!      *FILEDI*    CHARACTER FILE IDENTIFIER OF INPUT DATA.
!      *PATHI*     CHARACTER DIRECTORY OF INPUT DATA.
!      *RUNDO*     CHARACTER RUN  IDENTIFIER OF OUTPUT DATA.
!      *FILEDO*    CHARACTER FILE IDENTIFIER OF OUTPUT DATA.
!      *PATHO*     CHARACTER DIRECTORY OF OUTPUT DATA.
!      *LREAL*     LOGICAL, IF TRUE THE OUTPUT WILL ALWAYS BE REAL*4
!      *LSWAN*     LOGICAL, IF TRUE THE OUTPUT WILL FOR SWAN INPUT.
 
! ----------------------------------------------------------------------
 
!*    1. INITIALISATION.
!        ---------------
 
 1000 CONTINUE
 
!*    1.1 UNITS.
!         ------
 
      IU01 = 1
      IU02 = 2
      IU05 = 5
      IU06 = 6
      NANG=-1
      NFRE=-1
      MBMAX=-1
 
!*    1.2 READ USER INPUT.
!         ----------------
 
      CALL UIBOU (IU05, IU06, CDATEA, CDATEE, IDELPRF,                  &
     &                  CDATES, IDELFI, USERID,                         &
     &                  RUNDI, FILEDI, PATHI, RUNDO, FILEDO, PATHO,     &
     &                  LREAL, LSWAN)
      CDTOUT = CDATEA
      CDATE = CDATES
      NEWOUT = 1
 
! ----------------------------------------------------------------------
 
!*    2. FIRST COARSE GRID OUTPUT FILE.
!        ------------------------------
 
 2000 CONTINUE
 
!*    2.1 FETCH FIRST FILE.
!         -----------------
      PATH='' 
      CALL GSFILE (IU06, IU01, 0, CDATE, CDATE, FILEDI, 'G')

!*    2.2 READ BOUNDARY FILE HEADER.
!         --------------------------
 
      READ (IU01, ERR=6000, END=6000)                                   &
     &     XANG, XFRE, TH0, FR1, CO, XBOU, XDELC
      KL = NINT(XANG)
      ML = NINT(XFRE)
      NBOUNC  = NINT(XBOU)
      IDELPRC = NINT(XDELC)
      WRITE(IU06,*) ' '
      WRITE(IU06,*) ' INPUT FILE HEADER:'
      WRITE(IU06,*) ' NO. OF DIRECTIONS IS      KL     = ', KL
      WRITE(IU06,*) ' NO. OF FREQUENCIES IS     ML     = ', ML
      WRITE(IU06,*) ' FIRST DIRECTION IS        TH0    = ', TH0
      WRITE(IU06,*) ' FIRST FREQUENCY IS        FR(1)  = ', FR1
      WRITE(IU06,*) ' FREQUENCY RATIO IS        CO     = ', CO
      WRITE(IU06,*) ' NO. OF BOUNDRAY POINTS IS NBOUNC = ', NBOUNC
      WRITE(IU06,*) ' TIME STEP OF DATA IS      IDELPRC= ', IDELPRC
 
!*    2.3 CHECK DIMENSIONS.
!         -----------------
 
      IF (NANG .EQ. -1 .AND. NFRE .EQ. -1 .AND. MBMAX.EQ.-1) THEN
        NANG=KL
        NFRE=ML
        MBMAX=NBOUNC

        ALLOCATE(F(MBMAX,NANG,NFRE,2))
        ALLOCATE(FL(NANG,NFRE))
        ALLOCATE(FMEAN(MBMAX,2))
        ALLOCATE(EMEAN(MBMAX,2))
        ALLOCATE(THQ(MBMAX,2))
        ALLOCATE(XLAT(MBMAX))
        ALLOCATE(XLON(MBMAX))
        ALLOCATE(FR(NFRE))
      ELSE IF (KL.GT.NANG .OR. ML.GT.NFRE .OR. NBOUNC.GT.MBMAX) THEN
         WRITE(IU06,*) '*******************************************'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '*    FATAL ERROR PROGRAM BOUNINT.         *'
         WRITE(IU06,*) '*    ============================         *'
         WRITE(IU06,*) '* ONE OR MORE DIMENSIONS ARE TO SMALL.    *'
         WRITE(IU06,*) '* NO. OF DIRECTIONS IS      KL     = ', KL
         WRITE(IU06,*) '*    DIMENSION IS           NANG   = ', NANG
         WRITE(IU06,*) '* NO. OF FREQUENCIES IS     ML     = ', ML
         WRITE(IU06,*) '*             DIMENSION IS  NFRE   = ', NFRE
         WRITE(IU06,*) '* NO. OF BOUNDRAY POINTS IS NBOUNC = ', NBOUNC
         WRITE(IU06,*) '*              DIMENSION IS MBMAX  = ', MBMAX
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '*******************************************'
         CALL ABORT1
      ENDIF
 
!*    2.4 GENERATE FREQUENCY ARRAY.
!         -------------------------
 
      FR(1) = FR1 
      DO M=2,ML
         FR(M) = CO*FR(M-1)
      ENDDO
 
!*    2.5 OPEN FILE AND WRITE NEW BOUNDARY HEADER.
!         ----------------------------------------
 
      IF (NEWOUT.EQ.1) THEN

         CALL GSFILE (IU06, IU02, 0, CDATE, CDATE, FILEDO, 'G')

         XDELF = REAL(IDELPRF,JWRB)
         IF(LREAL)THEN
           WRITE(IU02) REAL(XANG,4), REAL(XFRE,4), REAL(TH0,4),         &
     &       REAL(FR(1),4), REAL(CO,4), REAL(XBOU,4), REAL(XDELF,4)
         ELSE
           WRITE(IU02) XANG, XFRE, TH0, FR(1), CO, XBOU, XDELF
         ENDIF
         NEWOUT = 0
      ENDIF
 
! ----------------------------------------------------------------------
 
!*    3. FIRST BOUNDARY VALUES.
!        ----------------------
 
 3000 CONTINUE
 
!*    3.1 READ BOUNDARY VALUES.
!         ---------------------
 
      DO IJ=1,NBOUNC
         READ (IU01, ERR=6001, END=6001)                                &
     &     XLON(IJ), XLAT(IJ), CDATE1, EMEAN(IJ,1),                     &
     &     THQ(IJ,1), FMEAN(IJ,1)
         READ (IU01, ERR=6002, END=6002) ((F(IJ,K,M,1),K=1,KL),M=1,ML)
      ENDDO
      WRITE(IU06,*) ' '
      WRITE(IU06,*) ' DATE OF FIRST INPUT IS CDATE1 = ', CDATE1
 
!*    3.2 ESTIMATE DATE OF NEXT BOUNDARY VALUES.
!         --------------------------------------
 
      CDATE2 = CDATE1
      CALL INCDATE (CDATE2,IDELPRC)
 
!*    3.3 IF SECOND DATE IS BEFORE FIRST OUTPUT DATE.
!         -------------------------------------------
 
      IF (CDATE2.LE.CDATEA) THEN
         WRITE(IU06,*) ' '
         WRITE(IU06,*) ' INPUT IS SKIPPED, BECAUSE ESTIMATED DATE '
         WRITE(IU06,*) ' OF NEXT SPECTRA IS LE THAN START DATE'
         IF (CDATE2.GT.CDATE) THEN
 
!*    3.3.1 CLOSE FILE, UPDATE FILE DATE AND BRANCH BACK TO 2.
!           --------------------------------------------------
 
            CALL INCDATE (CDATE,IDELFI)
            CLOSE (UNIT=IU01)
            CLOSE (UNIT=IU02)
            NEWOUT = 1
            GOTO 2000
         ELSE
 
!*    3.3.2 BRANCH BACK TO 3. TO READ SPECTRA FOR NEXT INPUT TIME.
!           ------------------------------------------------------
 
            GOTO 3000
         ENDIF
      ENDIF
 
! ----------------------------------------------------------------------
 
!*    4. LOOP OVER INPUT VALUES.
!        -----------------------
 
 4000 CONTINUE
 
!*    4.1 WRITE FIRST VALUES TO OUTPUT IF REQUESTED.
!         ------------------------------------------
 
      IF (CDTOUT.LE.CDATE1) THEN
        DO IJ=1,NBOUNC
          IF(LREAL) THEN
            IF(LSWAN)THEN
            WRITE(IU02) REAL(XLON(IJ),4), REAL(XLAT(IJ),4),CDTOUT(3:14), &
     &      REAL(EMEAN(IJ,1),4), REAL(THQ(IJ,1),4), REAL(FMEAN(IJ,1),4)
            WRITE(IU02) ((REAL(F(IJ,K,M,1),4),K=1,KL),M=1,ML)
            ELSE
            WRITE(IU02) REAL(XLON(IJ),4), REAL(XLAT(IJ),4), CDTOUT,      &
     &       REAL(EMEAN(IJ,1),4), REAL(THQ(IJ,1),4), REAL(FMEAN(IJ,1),4)
            WRITE(IU02) ((REAL(F(IJ,K,M,1),4),K=1,KL),M=1,ML)
            ENDIF
          ELSE
            IF(LSWAN)THEN
            WRITE(IU02) XLON(IJ), XLAT(IJ), CDTOUT(3:14),                &
     &        EMEAN(IJ,1), THQ(IJ,1), FMEAN(IJ,1)
            ELSE
            WRITE(IU02) XLON(IJ), XLAT(IJ), CDTOUT,                      &
     &        EMEAN(IJ,1), THQ(IJ,1), FMEAN(IJ,1)
            ENDIF
            WRITE(IU02) ((F(IJ,K,M,1),K=1,KL),M=1,ML)
          ENDIF
        ENDDO
        WRITE(IU06,*) ' INPUT DATA COPIED TO OUTPUT IDTOUT = ', CDTOUT
        CALL INCDATE (CDTOUT, IDELPRF)
        IF (CDTOUT.GE.CDATEE) GOTO 5000
        IF (CDTOUT.LE.CDATE1) THEN
          WRITE(IU06,*) ' '
          WRITE(IU06,*) ' NEXT OUTPUT DATE IS LE FIRST INPUT DATE '
          WRITE(IU06,*) ' SAME DATA WILL BE STORED FOR THIS DATE'
          GOTO 4000
        ENDIF
      ENDIF
 
!*    4.2 NEXT INPUT BOUNDARY VALUES.
!         ---------------------------
 
!*    4.2.1 NEW INPUT AND OUTPUT FILES?
!           ---------------------------
 
      IF (CDATE1.EQ.CDATE) THEN
        CALL INCDATE (CDATE,IDELFI)

        CLOSE (UNIT=IU02)
        CALL GSFILE (IU06, IU02, 0, CDATE, CDATE, FILEDO, 'G')

        CLOSE (UNIT=IU01)
        CALL GSFILE (IU06, IU01, 0, CDATE, CDATE, FILEDI, 'G')

        READ (IU01, ERR=6000, END=6000)                                 &
     &        XANG, XFRE, TH0, FR(1), CO, XBOU, XDELC
        IF(LREAL)THEN
          WRITE(IU02) REAL(XANG,4), REAL(XFRE,4), REAL(TH0,4),          &
     &      REAL(FR(1),4), REAL(CO,4), REAL(XBOU,4), REAL(XDELF,4)
        ELSE
          WRITE(IU02) XANG, XFRE, TH0, FR(1), CO, XBOU, XDELF
        ENDIF
      ENDIF
 
!*    4.2.2 READ NEXT INPUT.
!           ----------------
 
      DO IJ=1,NBOUNC
        READ (IU01, ERR=6003, END=6003)                                 &
     &    XLON(IJ), XLAT(IJ), CDATE2, EMEAN(IJ,2),                      &
     &    THQ(IJ,2), FMEAN(IJ,2)
        READ (IU01, ERR=6004, END=6004)                                 &
     &   ((F(IJ,K,M,2),K=1,KL),M=1,ML)
      ENDDO
      WRITE(IU06,*) ' DATE OF SECOND INPUT IS CDATE2 = ', CDATE2
 
!*    4.3 CHECK CONSISTENCY BEWEEN DATES.
!         -------------------------------
 
      CALL DIFDATE (CDATE1, CDATE2, IDEL12)
      DEL12 = REAL(IDEL12)
      IF (IDEL12.NE.IDELPRC) THEN
         WRITE(IU06,*) '*******************************************'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '*    FATAL ERROR PROGRAM BOUNINT.         *'
         WRITE(IU06,*) '*    ============================         *'
         WRITE(IU06,*) '* DATE INCREMENTS DO NOT MATCH.           *'
         WRITE(IU06,*) '* INCREMENT FROM HEADER IS IDELPRC = ', IDELPRC
         WRITE(IU06,*) '* COMPUTED FROM DATES IS   IDEL12  = ', IDEL12
         WRITE(IU06,*) '* DATE OF FIRST VALUES IS  CDATE1  = ', CDATE1
         WRITE(IU06,*) '* DATE OF SECOND VALUES IS CDATE2  = ', CDATE2
         WRITE(IU06,*) '* IS THE NO. OF BOUNDARY POINTS CORRECT?  *'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '*******************************************'
         CALL ABORT1
      ENDIF
 
!*    4.4 LOOP OVER INTERMEDIATE OUTPUT TIMES.
!         ------------------------------------
 
 4400 CONTINUE
      CALL DIFDATE (CDATE1, CDTOUT, IDEL1L)
      DEL1L = REAL(IDEL1L)
 
!*    4.4.1 LOOP OVER BOUNDARY POINTS.
!           --------------------------
 
      DO IJ=1,NBOUNC
 
!*      4.4.1.1 INTERPOLATE.
!               ------------
 
         CALL INTSPEC (NFRE, NANG, ML, KL, FR, DEL12, DEL1L,            &
     &                 F(IJ,:,:,1), FMEAN(IJ,1),                        &
     &                 EMEAN(IJ,1), THQ(IJ,1),                          &
     &                 F(IJ,:,:,2), FMEAN(IJ,2),                        &
     &                 EMEAN(IJ,2), THQ(IJ,2),                          &
     &                 FL, FMEAN_PT, EMEAN_PT, THQ_PT)
 
!*      4.4.1.2 WRITE TO OUTPUT FILE.
!               ---------------------
 
         IF(LREAL)THEN
           IF(LSWAN)THEN
             WRITE(IU02) REAL(XLON(IJ),4),REAL(XLAT(IJ),4),CDTOUT(3:14), &
     &       REAL(EMEAN_PT,4), REAL(THQ_PT,4), REAL(FMEAN_PT,4)
           ELSE
              WRITE(IU02) REAL(XLON(IJ),4),REAL(XLAT(IJ),4), CDTOUT,     &
     &        REAL(EMEAN_PT,4), REAL(THQ_PT,4), REAL(FMEAN_PT,4)
           ENDIF
           WRITE(IU02) ((REAL(FL(K,M),4),K=1,KL),M=1,ML)
         ELSE
           IF(LSWAN)THEN
             WRITE(IU02) XLON(IJ), XLAT(IJ), CDTOUT(3:14),               &
     &       EMEAN_PT, THQ_PT, FMEAN_PT
           ELSE
             WRITE(IU02) XLON(IJ), XLAT(IJ), CDTOUT,                     &
     &       EMEAN_PT, THQ_PT, FMEAN_PT
           ENDIF
           WRITE(IU02) ((FL(K,M),K=1,KL),M=1,ML)
         ENDIF
      ENDDO
      WRITE(IU06,*) ' INTERPOLATED DATA TO OUTPUT IDTOUT = ', CDTOUT
 
!*    4.4.2 UPDATE AND CHECK DATES.
!           -----------------------
 
      CALL INCDATE (CDTOUT,IDELPRF)
 
!*    4.4.2.1 BRANCH TO 5. IF ALL DONE.
!             -------------------------
 
      IF (CDTOUT.GT.CDATEE) GOTO 5000
 
!*    4.4.2.2 BRANCH BACK TO 4.4 FOR NEXT INTERMEDIATE DATE.
!             ----------------------------------------------
 
      IF (CDTOUT.LT.CDATE2) GOTO 4400
 
!*    4.4.3 COPY SECOND TO FIRST VALUES.
!           ----------------------------
 
      CDATE1 = CDATE2
      DO IJ=1,NBOUNC
        FMEAN(IJ,1) = FMEAN(IJ,2)
        EMEAN(IJ,1) = EMEAN(IJ,2)
        THQ(IJ,1) =  THQ(IJ,2)
      ENDDO
      DO M=1,ML
        DO K=1,KL
          DO IJ=1,NBOUNC
            F(IJ,K,M,1) = F(IJ,K,M,2)
          ENDDO
        ENDDO
      ENDDO
 
!*    4.4.4 BRANCH BACK TO 4. FOR NEXT INPUT.
!           ---------------------------------
 
      GOTO 4000
 
! ----------------------------------------------------------------------
 
!*    5. SAVE LAST OUTPUT FILE AND TERMINATE PROGRAM.
!        --------------------------------------------
 
 5000 CONTINUE
      CLOSE (UNIT=IU02, STATUS='KEEP')

      STOP
 
! ----------------------------------------------------------------------
 
!*    6. ERROR MESSAGES.
!        ---------------
 
 6000 CONTINUE
         WRITE(IU06,*) '*******************************************'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '*    FATAL ERROR PROGRAM BOUNINT.         *'
         WRITE(IU06,*) '*    ============================         *'
         WRITE(IU06,*) '* PROGRAM TRIES TO READ                   *'
         WRITE(IU06,*) '* HEADER OF BOUNDARY VALUES               *'
         WRITE(IU06,*) '* END OF FILE OR READ ERROR.              *'
         WRITE(IU06,*) '* UNIT IS IU01 = ', IU01
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '*******************************************'
         CALL ABORT1
 6001 CONTINUE
         WRITE(IU06,*) '*******************************************'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '*    FATAL ERROR PROGRAM BOUNINT.         *'
         WRITE(IU06,*) '*    ============================         *'
         WRITE(IU06,*) '* PROGRAM TRIES TO READ                   *'
         WRITE(IU06,*) '* A HEADER OF FIRST SET OF SPECTRA        *'
         WRITE(IU06,*) '* SPECTRA COUNTER IS IJ = ', IJ
         WRITE(IU06,*) '* END OF FILE OR READ ERROR.              *'
         WRITE(IU06,*) '* UNIT IS          IU01 = ', IU01
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '*******************************************'
         CALL ABORT1
 6002 CONTINUE
         WRITE(IU06,*) '*******************************************'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '*    FATAL ERROR PROGRAM BOUNINT.         *'
         WRITE(IU06,*) '*    ============================         *'
         WRITE(IU06,*) '* PROGRAM TRIES TO READ                   *'
         WRITE(IU06,*) '* A SPECTRUM OF FIRST SET OF SPECTRA      *'
         WRITE(IU06,*) '* END OF FILE OR READ ERROR.              *'
         WRITE(IU06,*) '* DATE IS        CDATE1 = ', CDATE1
         WRITE(IU06,*) '* SPECTRA COUNTER IS IJ = ', IJ
         WRITE(IU06,*) '* UNIT IS          IU01 = ', IU01
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '*******************************************'
         CALL ABORT1
 6003 CONTINUE
         WRITE(IU06,*) '*******************************************'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '*    FATAL ERROR PROGRAM BOUNINT.         *'
         WRITE(IU06,*) '*    ============================         *'
         WRITE(IU06,*) '* END OF FILE OR READ ERROR.              *'
         WRITE(IU06,*) '* PROGRAM TRIES TO READ                   *'
         WRITE(IU06,*) '* A HEADER OF SECOND SET OF SPECTRA       *'
         WRITE(IU06,*) '* SPECTRA COUNTER IS IJ = ', IJ
         WRITE(IU06,*) '* UNIT IS          IU01 = ', IU01
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '*******************************************'
         CALL ABORT1
 6004 CONTINUE
         WRITE(IU06,*) '*******************************************'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '*    FATAL ERROR PROGRAM BOUNINT.         *'
         WRITE(IU06,*) '*    ============================         *'
         WRITE(IU06,*) '* END OF FILE OR READ ERROR.              *'
         WRITE(IU06,*) '* PROGRAM TRIES TO READ                   *'
         WRITE(IU06,*) '* A SPECTRUM OF SECOND SET OF SPECTRA     *'
         WRITE(IU06,*) '* SPECTRA COUNTER IS IJ = ', IJ
         WRITE(IU06,*) '* DATE IS        CDATE2 = ', CDATE2
         WRITE(IU06,*) '* UNIT IS          IU01 = ', IU01
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
         WRITE(IU06,*) '*                                         *'
         WRITE(IU06,*) '*******************************************'
         CALL ABORT1

      END PROGRAM BOUINT
