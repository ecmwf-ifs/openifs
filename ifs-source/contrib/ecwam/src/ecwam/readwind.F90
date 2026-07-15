! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE READWIND (CDTWIR, FILNM, LLNOTOPENED, IREAD,   &
 &                   NXS, NXE, NYS, NYE, FIELDG)

!***  *READWIND* - PROGRAM TO FORMAT WIND FIELDS FROM MARS 
!                  FROM INPUT FILE sfcwind.

!     PETER JANSSEN      ECMWF        FEBRUARY 1987 


!    MODIFIED BY:                                                 
!    ------------                                                       

!    LIANA ZAMBRESKY     GKSS/ECMWF   JULY 1988

!        1. PASSING OF ARRAYS AS FORMAL PARMETERS TO SAVE MEMORY
!        2. PROCESSING WINDS FOR ONLY ONE WIND TIME STEP

!    PEDRO VITERBO       ECMWF        OCTOBER 1988

!        DEALING WITH MARS PACKED DATA.

!    H. GUNTHER      ECMWF/GKSS       OCTOBER 1991

!        GRIB EDITION 1.

!     J. DOYLE    NRL/ECMWF     OCT. 1996
!                 DECODES GRID WIND FIELDS WHICH ARE PASSED AS ARRAYS
!                 WHEN MODEL IS COUPLED TO ATMOPHERIC MODEL

!    B. HANSEN    ECMWF 1997
!                 RENAME SUBROUTINE TO MATCH THE NAME OF THE FILE.
!                 USE INCLUDE FILES FOR ALL COMMON BLOCKS.

!    J. BIDLOT    ECMWF 1998 RENAMED READWND_MARS AS READWIND

!    S. ABDALLA   ECMWF 2001
!                 GETS ATMOSPHERIC PARAMETERS FOR GUSTINESS
!                 AND AIR DENSITY COMPUTATIONS
!
!    P.A.E.M. JANSSEN   ECMWF 2003
!                 USE NEUTRAL WINDS RATHER THAN REAL WINDS
!                 TRANSFER AIR DENSITY CALCULATION TO IFS

!    J. BIDLOT    ECMWF 2008  REMOVE PART FOR ATMOSPHERIC FIELDS
!                             TO *IFSTOWAM*

!    J. BIDLOR    ECMWF 2010 USE GRIB API

!     PURPOSE                                                       
!     -------                                                      

!     *READWIND*
!         READS TARGET FILES FROM THE MARS ARCHIVE AND INFERS FROM THIS 
!             A) THE DATE OF THE WIND FIELD
!             B) THE DEFINITION OF THE MARS GRID
!             C) THE U AND V COMPONENT OF THE WIND FIELD
!             D) THE TYPE OF WIND FIELD AND
!             E) THE UNPACKED WINDS

!     INTERFACE
!     ---------

!     *CALL* *READWIND (CDTWIR, FILNM, LLNOTOPENED, IREAD,
!                       NXS, NXE, NYS, NYE, FIELDG)
!
!        *CDTWIR* - DATE/TIME OF THE DATA READ
!        *FILNM*    FILENAME OF INPUT FILE
!        *LLNOTOPENED*  TRUE IF THE INPUT FILE HAS TO BE OPENED
!        *IREAD*    PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!                   (IF NEEDED)
!        *NXS:NXE*  FIRST DIMENSION OF FIELDG
!        *NYS:NYE*  SECOND DIMENSION OF FIELDG
!        *FIELDG* - INPUT FORCING FIELDS ON THE WAVE MODEL GRID


!     EXTERNALS
!     ---------

!     *ABORT1*         TERMINATES PROCESSING
!     *INCDATE*        INCREMENTS DATE
!     *KGRIBSIZE*
!     *FLUSH*
!     *GRB2WGRD*       TRANSFORM GRIB FIELDS TO WAM GRID FIELDS

! --------------------------------------------------------------------- 

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWICE   , ONLY : LICERUN  ,IPARAMCI ,LICETH
      USE YOWMAP   , ONLY : IRGG     ,NLONRGG  ,NGY      ,NIBLO    ,CLDOMAIN
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPRECI
      USE YOWPARAM , ONLY : SWAMPWIND,SWAMPWIND2,DTNEWWIND,LTURN90 ,LWDINTS  ,    & 
     &                      SWAMPCIFR,LLUNSTR
      USE YOWSTAT  , ONLY : CDATEA   ,IDELWI   ,LADEN    ,LGUST
      USE YOWTEST  , ONLY : IU06
      USE YOWWNDG  , ONLY : ICODE    ,IWPER    ,ICOORD 
      USE YOWWIND  , ONLY : IUNITW   , NBITW    ,CWDFILE ,LLWSWAVE ,LLWDWAVE 
      USE YOWPCONS , ONLY : RAD      ,ZMISS    ,ROAIR    ,WSTAR0
#ifdef WAM_HAVE_UNWAM
      USE YOWPD, ONLY : MNP => npa
#endif

      USE MPL_MODULE, ONLY : MPL_BARRIER, MPL_BROADCAST
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWGRIB, ONLY : JPKSIZE_T, &
                        & JPGRIB_BUFFER_TOO_SMALL, JPGRIB_END_OF_FILE, &
                        & JPGRIB_SUCCESS, &
                        & IGRIB_NEW_FROM_MESSAGE, &
                        & IGRIB_OPEN_FILE, &
                        & IGRIB_READ_FROM_FILE, &
                        & IGRIB_RELEASE
      USE YOWABORT, ONLY : WAM_ABORT

! --------------------------------------------------------------------- 

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "grib2wgrid.intfb.h"
#include "incdate.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "kgribsize.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD 
      CHARACTER(LEN=14), INTENT(INOUT) :: CDTWIR
      CHARACTER(LEN=24), INTENT(INOUT) :: FILNM
      LOGICAL, INTENT(INOUT) :: LLNOTOPENED
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FIELDG


      INTEGER(KIND=JWIM) :: NFLD  
      INTEGER(KIND=JWIM) :: I, J, IVAR, JSN
      INTEGER(KIND=JWIM) :: ISIZE
      INTEGER(KIND=JWIM) :: IFORP, IPARAM, KZLEV, IDM 
      INTEGER(KIND=JWIM) :: IWTIME, IDTTURN
      INTEGER(KIND=JWIM) :: LNAME 
      INTEGER(KIND=JWIM) :: IRET
      INTEGER(KIND=JWIM) :: KGRIB_HANDLE
      INTEGER(KIND=JWIM) :: IDUM(2)
      INTEGER(KIND=JWIM), ALLOCATABLE :: KGRIB(:)
      INTEGER(KIND=JWIM) :: NLONRGG_LOC(NGY)

      INTEGER(KIND=JPKSIZE_T) :: KBYTES

      REAL(KIND=JWRB) :: ZDUM
      REAL(KIND=JWRB) :: UWIND, VWIND, WSPEED, WTHETA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NXS:NXE, NYS:NYE) :: WORK

      CHARACTER(LEN=14) :: CDTTURN
      CHARACTER(LEN=14), SAVE :: CWDDATE

      LOGICAL :: LLABORT
      LOGICAL :: LLEXIST
      LOGICAL, ALLOCATABLE :: LLNOTREAD(:)
! --------------------------------------------------------------------  

      IF (LHOOK) CALL DR_HOOK('READWIND',0,ZHOOK_HANDLE)

      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
        NLONRGG_LOC(:)=MNP
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ELSE
        NLONRGG_LOC(:) = NLONRGG(:)
      ENDIF

      ICODE=0

      IF (.NOT.(CLDOMAIN == 's' .OR. LWDINTS)) THEN
!       IF IT IS NOT A SWAMP CASE OR A CASE FOR WHICH A TIME SERIES
!       IS PRESCRIBED.
!       
!       WIND FIELDS ARE SUPPLIED FROM INPUT FILE sfcwindin
!       AS IS THE CASE IN THE STAND ALONE VERSION OR AT INITIAL INPUT
!       TIME IN COUPLED MODEL.

!       1.0 INITIALISE VALUES. 
!           ------------------ 

        ICOORD = 1

        FILNM='sfcwindin'
     
        IF (IRANK == IREAD) THEN
          LLEXIST=.FALSE.
          LNAME = LEN_TRIM(FILNM)
          INQUIRE(FILE=FILNM(1:LNAME),EXIST=LLEXIST)
          IF (.NOT. LLEXIST) THEN
            WRITE (IU06,*) '*************************************'
            WRITE (IU06,*) '*                                   *'
            WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
            WRITE (IU06,*) '*  IN READWIND:                     *'
            WRITE (IU06,*) '*  COULD NOT FIND FILE ',FILNM
            WRITE (IU06,*) '*                                   *'
            WRITE (IU06,*) '*************************************'
            WRITE (*,*) '*************************************'
            WRITE (*,*) '*                                   *'
            WRITE (*,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
            WRITE (*,*) '*  IN READWIND:                     *'
            WRITE (*,*) '*  COULD NOT FIND FILE ',FILNM
            WRITE (*,*) '*                                   *'
            WRITE (*,*) '*************************************'
            CALL ABORT1
          ENDIF

          IF (LLNOTOPENED) THEN
            CALL IGRIB_OPEN_FILE(IUNITW,FILNM(1:LNAME),'r') 
            LLNOTOPENED = .FALSE.
          ENDIF
        ENDIF  


! --------------------------------------------------------------------  

!*      2.0 READ MARS U and V WIND COMPONENTS  IN GRIB CODE FORMAT.
!           AND SEA ICE FACTION (OR SST) IF LICERUN=TRUE
!           WIND SPEED FROM PREVIOUS WAVE MODEL RUNS (IF LLWSWAVE)
!           WIND DIRECTION FROM PREVIOUS WAVE MODEL RUNS (IF LLWDWAVE)
!           SURFACE AIR DENSITY (IF LADEN AND NOT COUPLED TO IFS)
!           GUSTINESS (IF LGUST AND NOT COUPLED TO IFS)
!           -------------------------------------------------------

        IWPER = 1

        IF (LICERUN) THEN
          NFLD=3
        ELSE
          NFLD=2
          IPARAMCI=31
        ENDIF

        IF (LICETH) NFLD=NFLD+1
        IF (LLWSWAVE) NFLD=NFLD+1
        IF (LLWDWAVE) NFLD=NFLD+1

!       AIR DENSITY AND GUSTINESS ARE ONLY PROVIDED FROM A FILE IF STAND ALONE RUN
!!      for coupled run we re-initialise witth he value provided by IFS, even though
!!      strickly speaking, it should be provided as part of the restart but it is better than
!!      setting it to a constant !!
        IF (.NOT.LWCOU .AND. LADEN) NFLD=NFLD+1

!!      for coupled run we re-initialise with the value provided by IFS, even though
!!      strickly speaking, it should be provided as part of the restart but it is better than
!!      setting it to a constant !!
        IF (.NOT.LWCOU .AND. LGUST) NFLD=NFLD+1

        ALLOCATE(LLNOTREAD(NFLD))
        LLNOTREAD=.TRUE.

        LLABORT=.FALSE.


!       LOOP OVER INPUT

        WND: DO IVAR=1,NFLD

2002      CONTINUE

          IF (IRANK == IREAD) THEN
1021        ISIZE=NBITW
            KBYTES=ISIZE*NPRECI
            IF (.NOT.ALLOCATED(KGRIB)) ALLOCATE(KGRIB(ISIZE))
            CALL IGRIB_READ_FROM_FILE(IUNITW,KGRIB,KBYTES,IRET)
            IF (IRET == JPGRIB_BUFFER_TOO_SMALL) THEN
!!!           *IGRIB_READ_FROM_FILE* does not read through the file if
!!!            the size is too small, so figure out the size and read again.
              CALL KGRIBSIZE(IU06, KBYTES, NBITW, 'READWIND')
              DEALLOCATE(KGRIB)
              GOTO 1021
            ELSEIF (IRET == JPGRIB_END_OF_FILE) THEN
              WRITE(IU06,*) '**********************************'
              WRITE(IU06,*) '* READWIND: END OF FILE ENCOUNTED'
              WRITE(IU06,*) '**********************************'
              CALL ABORT1
            ELSEIF (IRET /= JPGRIB_SUCCESS) THEN
              WRITE(IU06,*) '**********************************'
              WRITE(IU06,*) '* READWIND: FILE HANDLING ERROR'
              WRITE(IU06,*) '**********************************'
              CALL ABORT1
            ENDIF
          ENDIF

          CALL GSTATS(622,0)
          CALL MPL_BARRIER(CDSTRING='READWIND: KGRIB ')

!         SEND GRIB DATA TO THE OTHER PE'S
!         --------------------------------
          IF (NPROC > 1) THEN
            IF (IRANK == IREAD) THEN
              IDUM(1)=ISIZE
            ENDIF
            CALL MPL_BROADCAST(IDUM(1:1),KROOT=IREAD,KTAG=IVAR,         &
     &                         CDSTRING='READWIND IDUM:')
            IF (IRANK /= IREAD) THEN
              ISIZE=IDUM(1)
              ALLOCATE(KGRIB(ISIZE))
            ENDIF

            CALL MPL_BROADCAST(KGRIB(1:ISIZE),KROOT=IREAD,KTAG=IVAR,    &
     &                         CDSTRING='READWIND KGRIB:')

          ENDIF
          CALL GSTATS(622,1)

! ----------------------------------------------------------------------

!*        3.0 UNPACK MARS FIELDS.
!             -------------------
          KGRIB_HANDLE=-99
          CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,KGRIB)
          ZDUM=0.0_JWRB
          CALL GRIB2WGRID (IU06, NPROMA_WAM,                            &
     &                     KGRIB_HANDLE, KGRIB, ISIZE,                  &
     &                     LLUNSTR,                                     &
     &                     NGY, IRGG, NLONRGG_LOC,                      &
     &                     NXS, NXE, NYS, NYE,                          &
     &                     FIELDG%XLON, FIELDG%YLAT,                    &
     &                     ZMISS, ZDUM, ZDUM,                           &
     &                     CDTWIR, IFORP, IPARAM, KZLEV, IDM, IDM, WORK)

          CALL IGRIB_RELEASE(KGRIB_HANDLE)

          IF (IPARAM == 165 .OR. IPARAM == 33 .OR. IPARAM == 131 .OR.   &
     &        IPARAM == 180 ) THEN

            IF (LLNOTREAD(IVAR)) THEN
              DO J = NYS, NYE
                JSN = NGY-J+1
                DO I = NXS, MIN(NLONRGG_LOC(JSN), NXE)
                  FIELDG%UWND(I,J)=WORK(I,J)
                ENDDO
              ENDDO
              LLNOTREAD(IVAR)=.FALSE.
            ELSE
              LLABORT=.TRUE.
            ENDIF

            IF (IPARAM == 180) THEN
              ICODE = 2
            ELSE
              ICODE = 3
            ENDIF
          ELSEIF (IPARAM == 166 .OR. IPARAM == 34 .OR. IPARAM == 132    &
     &            .OR. IPARAM == 181 ) THEN

            IF (LLNOTREAD(IVAR)) THEN 
              DO J = NYS, NYE
                JSN = NGY-J+1
                DO I = NXS, MIN(NLONRGG_LOC(JSN), NXE)
                  FIELDG%VWND(I,J)=WORK(I,J)
                ENDDO
              ENDDO
              LLNOTREAD(IVAR)=.FALSE.
            ELSE
              LLABORT=.TRUE.
            ENDIF

            IF (IPARAM == 181) THEN
              ICODE = 2
            ELSE
              ICODE = 3
            ENDIF

          ELSEIF (IPARAM == 31 .OR. IPARAM == 139) THEN
             IPARAMCI=IPARAM

            IF (LLNOTREAD(IVAR)) THEN 
              DO J = NYS, NYE
                JSN = NGY-J+1
                DO I = NXS, MIN(NLONRGG_LOC(JSN), NXE)
                  FIELDG%CICOVER(I,J)=WORK(I,J)
                ENDDO
              ENDDO
              LLNOTREAD(IVAR)=.FALSE.
            ELSEIF (.NOT.LICERUN .AND.                                  &
     &              (IPARAM == 31 .OR. IPARAM == 139) ) THEN
!             SKIP SEA ICE MASK INFORMATION AS IT IS NOT NEEDED
              GOTO 2002
            ELSE
              LLABORT=.TRUE.
            ENDIF

          ELSEIF (IPARAM == 92) THEN

            IF (LLNOTREAD(IVAR)) THEN 
              DO J = NYS, NYE
                JSN = NGY-J+1
                DO I = NXS, MIN(NLONRGG_LOC(JSN), NXE)
                  FIELDG%CITHICK(I,J)=WORK(I,J)
                ENDDO
              ENDDO
              LLNOTREAD(IVAR)=.FALSE.
            ELSEIF (.NOT.LICETH .AND.                                   &
     &              (IPARAM == 92) ) THEN
!             SKIP SEA ICE THICKNESS INFORMATION AS IT IS NOT NEEDED
              GOTO 2002
            ELSE
              LLABORT=.TRUE.
            ENDIF

          ELSEIF (IPARAM == 245) THEN

            IF (LLWSWAVE) THEN
              IF (LLNOTREAD(IVAR)) THEN 
                DO J = NYS, NYE
                  JSN = NGY-J+1
                  DO I = NXS, MIN(NLONRGG_LOC(JSN), NXE)
                    IF (WORK(I,J) == ZMISS) THEN
                      FIELDG%WSWAVE(I,J)=0.0_JWRB
                    ELSE
                      FIELDG%WSWAVE(I,J)=WORK(I,J)
                    ENDIF
                  ENDDO
                ENDDO
                LLNOTREAD(IVAR)=.FALSE.
              ELSE
                LLABORT=.TRUE.
              ENDIF
            ENDIF

          ELSEIF (IPARAM == 249) THEN

            IF (LLWDWAVE) THEN
              IF (LLNOTREAD(IVAR)) THEN 
                DO J = NYS, NYE
                  JSN = NGY-J+1
                  DO I = NXS, MIN(NLONRGG_LOC(JSN), NXE)
                    IF (WORK(I,J) == ZMISS) THEN
                      FIELDG%WDWAVE(I,J)=0.0_JWRB
                    ELSE
!                     re-convert to WAM convention
                      FIELDG%WDWAVE(I,J)=RAD*(WORK(I,J)-180.0_JWRB)
                    ENDIF
                  ENDDO
                ENDDO
                LLNOTREAD(IVAR)=.FALSE.
              ELSE
                LLABORT=.TRUE.
              ENDIF
            ENDIF

          ELSEIF (IPARAM == 209) THEN

            IF (LLNOTREAD(IVAR)) THEN 
              DO J = NYS, NYE
                JSN = NGY-J+1
                DO I = NXS, MIN(NLONRGG_LOC(JSN), NXE)
                  IF (WORK(I,J) == ZMISS) THEN
                    FIELDG%AIRD(I,J)=ROAIR
                  ELSE
                    FIELDG%AIRD(I,J)=WORK(I,J)
                  ENDIF
                ENDDO
              ENDDO
              LLNOTREAD(IVAR)=.FALSE.
            ELSE
              LLABORT=.TRUE.
            ENDIF

          ELSEIF (IPARAM == 208) THEN

            IF (LLNOTREAD(IVAR)) THEN 
              DO J = NYS, NYE
                JSN = NGY-J+1
                DO I = NXS, MIN(NLONRGG_LOC(JSN), NXE)
                  IF (WORK(I,J) == ZMISS) THEN
                    FIELDG%WSTAR(I,J) = WSTAR0
                  ELSE
                    FIELDG%WSTAR(I,J)=WORK(I,J)
                  ENDIF
                ENDDO
              ENDDO
              LLNOTREAD(IVAR)=.FALSE.
            ELSE
              LLABORT=.TRUE.
            ENDIF

          ELSE
            WRITE(IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) ' +                                        +'
            WRITE(IU06,*) ' +    WARNING ERROR IN SUB. READWIND      +'
            WRITE(IU06,*) ' +    ==============================      +'
            WRITE(IU06,*) ' + SUSPICIOUS WIND OR SEA ICE FIELD PARAM +'
            WRITE(IU06,*) ' + PARAM IS = ', IPARAM
            WRITE(IU06,*) ' +                                        +'
            WRITE(IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++'
            CALL ABORT1
          ENDIF

          IF (LLABORT) THEN
            WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) ' +                                         +'
            WRITE(IU06,*) ' +    WARNING ERROR IN SUB. READWIND       +'
            WRITE(IU06,*) ' +    ==============================       +'
            WRITE(IU06,*) ' + WIND OR SEA ICE FIELD PARAM READ TWICE  +'
            WRITE(IU06,*) ' + PARAM IS = ', IPARAM
            WRITE(IU06,*) ' + SUSPECT                                 +'
            WRITE(IU06,*) ' + INCOMPLETE LIST OF INPUT PARAMETERS     +'
            WRITE(IU06,*) ' + THE LIST IS :                           +'
            WRITE(IU06,*) ' + 165 or 33 or 131 or 180 U-WIND COMPONENT+'
            WRITE(IU06,*) ' + 166 or 34 or 132 or 181 V-WIND COMPONENT+'
            IF (LICERUN) THEN
            WRITE(IU06,*) ' + 31 or 139 SEA ICE FRACTION OR SST       +'
            ENDIF
            IF (LICETH) THEN
            WRITE(IU06,*) ' + 92 SEA ICE THICKNESS                    +'
            ENDIF
            IF (LLWSWAVE) THEN
            WRITE(IU06,*) ' + 245  WIND SPEED FROM WAVE MODEL         +'
            ENDIF
            IF (LLWDWAVE) THEN
            WRITE(IU06,*) ' + 249  WIND DIRECTION FROM WAVE MODEL     +'
            ENDIF
            IF (LADEN) THEN
            WRITE(IU06,*) ' + 209  AIR DENSITY AT THE SURFACE         +'
            ENDIF
            IF (LGUST) THEN
            WRITE(IU06,*) ' + 208  FREE CONVECTIVE VELOCITY SCALE     +'
            ENDIF
            WRITE(IU06,*) ' +                                         +'
            WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++'
            CALL ABORT1
          ENDIF

          DEALLOCATE(KGRIB)
        ENDDO WND

        DEALLOCATE(LLNOTREAD)

      ELSE
!       SWAMP CASE : CREATE YOUR OWN WIND FIELD
!       =======================================
!       for the swamp case create your own wind fields
!       or if the file windforcing_time_series is present then
!       the wind forcing will be extracted from the wind speed
!       and direction time series specified herewith.

        ICODE = 3

        IF (LLNOTOPENED) THEN
           CDTWIR = CDATEA
           CALL INCDATE(CDTWIR,-IDELWI)
           IWPER = 0
           ICOORD = 1

           IF (LWDINTS) THEN
             IUNITW=IWAM_GET_UNIT(IU06,CWDFILE, 'r', 'f', 0, 'READWRITE')
             OPEN(IUNITW,FILE=CWDFILE,FORM='FORMATTED')
             CWDDATE=CDATEA
           ELSE
             IUNITW=-1
           ENDIF
  
           LLNOTOPENED = .FALSE.
  
        ENDIF

        CALL INCDATE(CDTWIR,IDELWI)

        IF (LWDINTS) THEN

           READ(IUNITW,*,END=110,ERR=120) IWTIME,WSPEED,WTHETA

           CALL INCDATE(CWDDATE,IWTIME)
           WRITE(IU06,'(a28,a14,1x,f6.2,1x,f6.1)')                      &
     &     '  WIND INPUT TIME SERIES AT ', CWDDATE,WSPEED,WTHETA

           WTHETA=RAD*WTHETA
           UWIND=-WSPEED*SIN(WTHETA)
           VWIND=-WSPEED*COS(WTHETA)
        ELSE
!         THE TRADITIONAL ONE GRID POINT FORCING

          CDTTURN = CDATEA
          IDTTURN=DTNEWWIND*3600
          CALL INCDATE(CDTTURN,IDTTURN)
          IF (SWAMPWIND2 <= 0.0_JWRB .OR. CDTWIR <= CDTTURN) THEN
            WRITE(IU06,*) ' READWIND - SWAMP CASE WITH WIND SPEED = ',  &
     &                    SWAMPWIND 
            UWIND=0.0_JWRB
            VWIND=SWAMPWIND
          ELSE
            WRITE(IU06,*) ' READWIND - SWAMP CASE WITH WIND SPEED = ',  &
     &                    SWAMPWIND2 
            IF (LTURN90) THEN
              UWIND=SWAMPWIND2
              VWIND=0.0_JWRB
            ELSE
              UWIND=0.0_JWRB
              VWIND=SWAMPWIND2
            ENDIF
          ENDIF
        ENDIF

        IPARAMCI=31

        DO J = NYS, NYE
          DO I = NXS, NXE
            FIELDG%UWND(I,J)=UWIND
            FIELDG%VWND(I,J)=VWIND
          ENDDO
        ENDDO

!!!! impose the northern part of the swamp domain to be covered with
!!!! a sea ice cover=SWAMPCIFR
        DO J = NYS, NYE/2
          DO I = NXS, NXE
            FIELDG%CICOVER(I,J)=SWAMPCIFR
          ENDDO
        ENDDO
        DO J=NYE/2+1,NYE
          DO I = NXS, NXE
            FIELDG%CICOVER(I,J)=0.0_JWRB
          ENDDO
        ENDDO

      ENDIF

      IF (LHOOK) CALL DR_HOOK('READWIND',1,ZHOOK_HANDLE)

      RETURN

110   CONTINUE
      WRITE(IU06,*) '****************************************'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*  ERROR IN READWIND:                  *'
      WRITE(IU06,*) '*  END OF FILE REACHED FOR ',CWDFILE
      WRITE(IU06,*) '*  THE WIND TIME SERIES SHOULD AT LEAST*'
      WRITE(IU06,*) '*  BE AS LONG AS THE INTENDED RUN !!!! *'
      WRITE(IU06,*) '*  THE RUN HAS ABORTED                 *'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '****************************************'
      CALL ABORT1

120   CONTINUE
      WRITE(IU06,*) '****************************************'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*  ERROR IN READWIND:                  *'
      WRITE(IU06,*) '*  ERROR READING FILE ',CWDFILE
      WRITE(IU06,*) '*  THE RUN HAS ABORTED                 *'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '****************************************'
      CALL ABORT1

      IF (LHOOK) CALL DR_HOOK('READWIND',1,ZHOOK_HANDLE)

END SUBROUTINE READWIND
