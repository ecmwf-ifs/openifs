! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

    SUBROUTINE UIPREP (IFORM, LLGRID)

! ----------------------------------------------------------------------

!**** *UIPREP* - ROUTINE TO READ NAMELIST INPUT FOR PREPROC.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       TO READ USER INPUT OF PROGRAM PREPROC AND CHECKS CONSISTENCY.

!**   INTERFACE.
!     ----------

!       *CALL* *UIPREP (IFORM, LLGRID)*
!          *IFORM*   - OUTPUT FORMAT OPTION = 1 UNFORMATED
!                                           = 2 FORMATED
!                                           OTHERWISE BOTH
!          *LLGRID*  - TRUE IF THE GRID DEFINITION HAS BEEN SPECIFIED
!                      IN INPUT FILE grid_description

!     METHOD.
!     -------
!         NAMELIST READ.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *ADJUST*    - CORRECTS LONGITUDE INPUT.

!     REFERENCE.
!     ----------

!       NONE.

!     MODIFICATIONS.                                                    
!     --------------                                                    

!        B. HANSEN    *ECMWF*      JAN 98
!           REPLACE FORMATED READ OF USER INPUT BY NAMELIST READ.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NFRE     ,NFRE_RED ,LLUNSTR
      USE YOWCPBO  , ONLY : IBOUNC   ,GBOUNC_MAX, GBOUNC ,              &
     &            AMOSOC   ,AMONOC   ,AMOEAC   ,AMOWEC
      USE YOWCINP  , ONLY : NOUT     ,XOUTW    ,XOUTS    ,XOUTE    ,    &
     &            XOUTN    ,NOUTD
      USE YOWFPBO  , ONLY : IBOUNF
      USE YOWFRED  , ONLY : IFRE1    ,FR1      ,FR       ,FRATIO
      USE YOWMAP   , ONLY : NGX      ,NGY      ,IPER     ,IRGG     ,    &
     &            AMOWEP   ,AMOSOP   ,AMOEAP   ,AMONOP   ,              &
     &            NIBLO    ,CLDOMAIN ,                                  &
     &            XDELLA   ,XDELLO   ,NLONRGG  ,LLOBSTRCT,LAQUA
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
#ifdef WAM_HAVE_UNWAM
      USE YOWUNPOOL, ONLY : LPREPROC, LVECTOR, IVECTOR
#endif
      USE YOWABORT, ONLY : WAM_ABORT

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "adjust.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "mfr.intfb.h"

      INTEGER(KIND=JWIM), INTENT(OUT) :: IFORM
      LOGICAL, INTENT(OUT) :: LLGRID

      INTEGER(KIND=JWIM) :: K, M, I, II, KSN
      INTEGER(KIND=JWIM) :: IU05, IU
      INTEGER(KIND=JWIM) :: ISPECTRUNC
      INTEGER(KIND=JWIM) :: IOS, IOUTA, IOUTANEW, IDUM
      INTEGER(KIND=JWIM), ALLOCATABLE :: NDUMP(:)

      REAL(KIND=JWRB) :: ZOUTS, ZOUTN, ZOUTW, ZOUTE, IOUTD
      REAL(KIND=JWRB) :: WEST, EAST, DW, DE, DS, DN
      REAL(KIND=JWRB), ALLOCATABLE :: XDUMP(:)

      CHARACTER(LEN=70) :: CLINE, FILENAME
      LOGICAL :: LLEXISTS

#ifndef WAM_HAVE_UNWAM
      LOGICAL            :: LPREPROC
      LOGICAL            :: LVECTOR
      INTEGER(KIND=JWIM) :: IVECTOR
#endif

! ----------------------------------------------------------------------

      NAMELIST /NALINE/ CLINE, NFRE, NFRE_RED, FR1, IFRE1,              &
     &                  IRGG, XDELLA, XDELLO,                           &
     &                  AMOSOP, AMONOP, AMOWEP, AMOEAP,                 &
     &                  IFORM, ITEST, ITESTB,                           &
     &                  IBOUNC, IBOUNF, AMOSOC, AMONOC, AMOWEC, AMOEAC, &
     &                  NIBLO, CLDOMAIN,LLOBSTRCT,                      &
     &                  LAQUA, LLUNSTR, LPREPROC

      NAMELIST /NACORR/ ZOUTS, ZOUTN, ZOUTW, ZOUTE, IOUTD


! ----------------------------------------------------------------------

!*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
!        ---------------------------------------------

      CLINE  = ' '

      LLUNSTR = .FALSE. 
      LPREPROC = .FALSE. 
      LVECTOR   =.FALSE.
      IVECTOR   = 1
      NFRE     =   0
      NFRE_RED =   0
      FR1    =   0.0_JWRB
      IFRE1 = 1
      IRGG   =  -1
      XDELLA =   0.0_JWRB
      XDELLO =   0.0_JWRB
      AMOSOP =-100.0_JWRB
      AMONOP =-100.0_JWRB
      AMOWEP =   0.0_JWRB
      AMOEAP =   0.0_JWRB
      IFORM  =  -1
      ITEST  =  -1
      ITESTB =  -1
      IBOUNC =  -1
      IBOUNF =  -1
      AMOSOC =-100.0_JWRB
      AMONOC =-100.0_JWRB
      AMOWEC =   0.0_JWRB
      AMOEAC =   0.0_JWRB
      LLOBSTRCT = .TRUE.
      LAQUA =.FALSE. ! to force an aqua planet on irregular grid with xdella

!     THE FOLLOWING DEFAULT VALUES INDICATES THAT THEY ARE DETERMINED
!     INTERNALLY BY PREPROC EXCEPT IF PROVIDED IN THE NAMELIST
      NIBLO  = -1 ! 
      CLDOMAIN= '-'

! ----------------------------------------------------------------------

!*    1. READ NAMELIST NALINE.
!        ---------------------


      INQUIRE(FILE="procin", EXIST=LLEXISTS)
      IF (.NOT. LLEXISTS) THEN
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'+ SUBROUTINE UIPREP :                      +'
        WRITE(IU06,*)'+ READ NAMELIST FAILED                     +'
        WRITE(IU06,*)'+ NAMELIST FILENAME: procin                +'
        WRITE(IU06,*)'+ PROGRAM WILL ABORT                       +'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        CALL WAM_ABORT("Expected namelist file does not exist: 'procin'",__FILENAME__,__LINE__)
      ENDIF

      IU05 =  IWAM_GET_UNIT (IU06, 'procin', 'r', 'f', 0, 'READ')
      READ (IU05, NALINE)

      IF (NFRE_RED > NFRE ) THEN
        WRITE (IU06,*) '**********************************************'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '*       FATAL ERROR IN SUB. UIPREP           *'
        WRITE (IU06,*) '*       ==========================           *'
        WRITE (IU06,*) '* THE REDUCED NUMBER OF FREQUENCIES NFRE_RED *'
        WRITE (IU06,*) '* IS LARGER THAN THE TOTAL NUMNBER NFRE  !!  *'
        WRITE (IU06,*) '* NFRE_RED = ', NFRE_RED
        WRITE (IU06,*) '* NFRE     = ', NFRE
        WRITE (IU06,*) '**********************************************'
        CALL WAM_ABORT(__FILENAME__,__LINE__)
      ENDIF

      IF (LLUNSTR .AND. NFRE_RED /= NFRE ) THEN
        WRITE (IU06,*) '**********************************************'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '*       FATAL ERROR IN SUB. UIPREP           *'
        WRITE (IU06,*) '*       ==========================           *'
        WRITE (IU06,*) '* THE REDUCED NUMBER OF FREQUENCIES NFRE_RED *'
        WRITE (IU06,*) '* IS DIFFERENT THAN THE TOTAL NUMNBER NFRE   *'
        WRITE (IU06,*) '* NFRE_RED = ', NFRE_RED
        WRITE (IU06,*) '* NFRE     = ', NFRE
        WRITE (IU06,*) '* LLUNSTR  = ', LLUNSTR
        WRITE (IU06,*) '* The software has not yet been adapted     *'
        WRITE (IU06,*) '* Use NFRE_RED = NFRE = ',NFRE
        WRITE (IU06,*) '**********************************************'
        CALL WAM_ABORT(__FILENAME__,__LINE__)
      ENDIF

      IF( NFRE <= 0 ) CALL WAM_ABORT( "Expected positive value for NFRE", __FILENAME__, __LINE__ )
      IF( IFRE1 <= 0 ) CALL WAM_ABORT( "Expected positive value for IFRE1",  __FILENAME__, __LINE__ )

      ALLOCATE(FR(NFRE_RED))

      CALL MFR(NFRE_RED, IFRE1, FR1, FRATIO, FR)


      WRITE (IU06,'(2X, A70,/)') CLINE(1:70)


!     1.1 CHECK IF FILE grid_description IS PRESENT
!         IF IT IS THERE IT WILL SUPERSEDE THE NAMLIST DEFINITION

      FILENAME='grid_description'
      INQUIRE(FILE=FILENAME,EXIST=LLGRID)
      IF (LLGRID) THEN
        IU=IWAM_GET_UNIT(IU06,FILENAME,'S','F',0,'READWRITE')
        OPEN(IU,FILE=FILENAME,STATUS='OLD', FORM='FORMATTED')
        READ (IU,*) ISPECTRUNC
        READ (IU,*) AMONOP
        READ (IU,*) AMOSOP
        READ (IU,*) AMOWEP
        READ (IU,*) AMOEAP
        READ (IU,*) IPER
        READ (IU,*) IRGG
        READ (IU,*) NGY
        WRITE(IU06,*) "grid_description read in "
      ENDIF

! ----------------------------------------------------------------------


!*    2. PRINT CONTENTS OF NAMELIST AND PERFORM CONSISTENCY CHECKS.
!        ----------------------------------------------------------

!*    2.1 FREQUENCY AND DIRECTION GRID DEFINITIONS.

      WRITE (IU06,'("   FREQUENCY / DIRECTION GRID"/)')
      WRITE (IU06,'("   NUMBER OF FREQUENCIES IS NFRE = ",I6)') NFRE
      WRITE (IU06,'("   REDUCED NUMBER OF FREQ IS NFRE_RED = ",I6)') NFRE_RED
      IF (IFRE1 /= 1) THEN
        WRITE (IU06,'("!!MINIMUM FREQUENCY WAS RESET TO",F10.6)') FR(1)
      ELSE
        WRITE (IU06,'("  MINIMUM FREQUENCY IS  FR(1) = ",F10.6)') FR(1)
      ENDIF

!     SET DIMENSIONS.

      GBOUNC = 0
      DO II=1,GBOUNC_MAX
        IF (AMOSOC(II) == -100.0_JWRB .AND. AMONOC(II) == -100.0_JWRB .AND. &
     &    AMOWEC(II) == 0.0_JWRB .AND. AMOEAC(II) == 0.0_JWRB)EXIT
        GBOUNC=II
      ENDDO
      IF (IBOUNC == 1 .AND. GBOUNC <= 0 ) THEN
        WRITE (IU06,*) '**********************************************'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '*       FATAL ERROR IN SUB. UIPREP           *'
        WRITE (IU06,*) '*       ==========================           *'
        WRITE (IU06,*) '* THE PRODUCTION OF FINE GRID BOUNDARY WAS   *'
        WRITE (IU06,*) '* REQUESTED BUT NO BOUNDING BOX VALUES WERE  *'
        WRITE (IU06,*) '* SUPPLIED IN THE NAMELIST!                  *'
        WRITE (IU06,*) '**********************************************'
        CALL WAM_ABORT(__FILENAME__,__LINE__)
      ENDIF


! ----------------------------------------------------------------------

!*    2.2 OUTPUT GRID DEFINITIONS.

      WRITE (IU06,'(/"  OUTPUT GRID"/)')
      IF (IRGG == 1) THEN
        WRITE(IU06,'("   USE A REDUCED GAUSSIAN GRID")')
      ELSE
        WRITE(IU06,'("   USE A REGULAR LON LAT GRID")')
      ENDIF      

      IF (LAQUA) THEN
        LLOBSTRCT=.FALSE.
        WRITE(IU06,*) ' '
        WRITE(IU06,*) '   AN AQUA PLANET WAS REQUESTED '
        WRITE(IU06,*) '   THE BATHYMETRY WILL BE DEEP EVERYWHERE !!!'
        WRITE(IU06,*) ' '
      ENDIF

      IF (LLOBSTRCT) THEN
        WRITE(IU06,*) ' '
        WRITE(IU06,*) '   THE NEW FORM OF BATHYMETRY INPUT IS USED IN'
        WRITE(IU06,*) '   CONJUNCTION WITH THE OBSTRUCTION COEFFICIENT'
        WRITE(IU06,*) ' '
      ENDIF

      CALL ADJUST (AMOWEP, AMOEAP)

!*    SET DIMENSIONS.

      IF (LLGRID) THEN
        XDELLA = (AMONOP-AMOSOP)/(NGY-1)
        ALLOCATE(NLONRGG(NGY))

        NGX = 0
        DO K=1,NGY
          KSN=NGY-K+1
          READ(IU,*,IOSTAT=IOS) NLONRGG(KSN)
          IF( IOS < 0 ) THEN
            CALL WAM_ABORT("End of file reached before finishing NLONRGG",__FILENAME__,__LINE__)
          ENDIF
          IF( NLONRGG(KSN) <= 0 ) THEN
            CALL WAM_ABORT("Expected positive value of NLONRGG",__FILENAME__,__LINE__)
          ENDIF
          NGX = MAX(NGX,NLONRGG(KSN))
        ENDDO

        IF (IPER == 1) THEN
          XDELLO  = 360._JWRB/REAL(NGX)
          AMOEAP = AMOWEP + 360._JWRB - XDELLO
        ELSE
          XDELLO = (AMOEAP-AMOWEP)/(NGX-1)
        ENDIF


        CLOSE(IU)

      ELSE
!       RESET FOR AQUA PLANET IF SELECTED
        IF (LAQUA) THEN
          NGY = NINT((AMONOP-AMOSOP)/XDELLA) + 1
          WRITE (IU06,*) ' !! RESETING TO AQUA PLANET CONFIGURATION !!'
          AMONOP=90.0_JWRB
          AMOSOP=-90.0_JWRB
          NGY=INT((AMONOP-AMOSOP)/XDELLA)+1
          IRGG=1
        ENDIF

        IPER = 0
        IF (ABS(AMOEAP-AMOWEP+1.5_JWRB*XDELLO) >= 360.0_JWRB ) IPER = 1
        NGX = NINT((AMOEAP-AMOWEP)/XDELLO) + 1
        NGY = NINT((AMONOP-AMOSOP)/XDELLA) + 1

        ALLOCATE(NLONRGG(NGY))

      ENDIF

      WRITE (IU06,'("   RESOLUTION LAT-LON ",2F12.7)') XDELLA, XDELLO
      WRITE (IU06,'("    SOUTHERN LAT  NORTHERN LAT ",                  &
     &  " WESTERN LONG "," EASTERN LONG",                               &
     &  /,2X,4F14.3)') AMOSOP, AMONOP, AMOWEP, AMOEAP
      IF (IPER == 1) WRITE (IU06,*) '   THE GRID IS EAST-WEST PERIODIC'

      IF (CLDOMAIN == '-' ) THEN
        IF (IPER == 1) THEN
          CLDOMAIN = 'g'
        ELSE
          CLDOMAIN = 'm'
        ENDIF 
      ENDIF 

! ----------------------------------------------------------------------

!*    2.3 OUTPUT GRID CORRECTIONS READ FROM NAMELIST NACORR.
!         --------------------------------------------------

      IF (LLOBSTRCT) THEN
!     when the new bathymetry used, it is assumed that there is not
!     any need for manual correction of the data since it should have
!     been taken care when creating the file. 
        NOUT=0
      ELSE
        REWIND (IU05)
        NOUT=0
        IOS=-9
        IOUTA=500
        ALLOCATE(XOUTW(IOUTA))
        ALLOCATE(XOUTS(IOUTA))
        ALLOCATE(XOUTE(IOUTA))
        ALLOCATE(XOUTN(IOUTA))
        ALLOCATE(NOUTD(IOUTA))
        XOUTS  =-100.0_JWRB
        XOUTN  =-100.0_JWRB
        XOUTW  =   0.0_JWRB
        XOUTE  =   0.0_JWRB
        NOUTD  =   0
        CORR: DO
        READ(IU05, NACORR, ERR=3000, IOSTAT=IOS, END=2300)
        NOUT=NOUT+1
        IF (NOUT > IOUTA) THEN
          ALLOCATE(XDUMP(IOUTA))
          IOUTANEW=IOUTA+100
          WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++'
          WRITE (IU06,*) '+                                      +'
          WRITE (IU06,*) '+     WARNING IN SUB. UIPREP           +'
          WRITE (IU06,*) '+     =============================    +'
          WRITE (IU06,*) '+                                      +'
          WRITE (IU06,*) '+  NUMBER OF AREAS TO BE CORRECTED     +'
          WRITE (IU06,*) '+  EXCEEDS  DIMENSION   IOUTA = ', IOUTA
          WRITE (IU06,*) '+  THE DIMEMSION WAS RESET TO ',IOUTANEW
          WRITE (IU06,*) '+                                      +'
          WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++'

          XDUMP=XOUTW
          DEALLOCATE(XOUTW)
          ALLOCATE(XOUTW(IOUTANEW))
          DO IDUM=1,IOUTA
            XOUTW(IDUM)=XDUMP(IDUM)
          ENDDO

          XDUMP=XOUTS
          DEALLOCATE(XOUTS)
          ALLOCATE(XOUTS(IOUTANEW))
          DO IDUM=1,IOUTA
            XOUTS(IDUM)=XDUMP(IDUM)
          ENDDO

          XDUMP=XOUTE
          DEALLOCATE(XOUTE)
          ALLOCATE(XOUTE(IOUTANEW))
          DO IDUM=1,IOUTA
            XOUTE(IDUM)=XDUMP(IDUM)
          ENDDO

          XDUMP=XOUTN
          DEALLOCATE(XOUTN)
          ALLOCATE(XOUTN(IOUTANEW))
          DO IDUM=1,IOUTA
            XOUTN(IDUM)=XDUMP(IDUM)
          ENDDO

          DEALLOCATE(XDUMP)

          ALLOCATE(NDUMP(IOUTA))
          NDUMP=NOUTD
          DEALLOCATE(NOUTD)
          ALLOCATE(NOUTD(IOUTANEW))
          DO IDUM=1,IOUTA
            NOUTD(IDUM)=NDUMP(IDUM)
          ENDDO
          DEALLOCATE(NDUMP)

          IOUTA=IOUTANEW

          CALL ADJUST (ZOUTW, ZOUTE)
          XOUTS(NOUT) = ZOUTS
          XOUTN(NOUT) = ZOUTN
          XOUTW(NOUT) = ZOUTW
          XOUTE(NOUT) = ZOUTE
          NOUTD(NOUT) = IOUTD
        ELSE
          CALL ADJUST (ZOUTW, ZOUTE)
          XOUTS(NOUT) = ZOUTS
          XOUTN(NOUT) = ZOUTN
          XOUTW(NOUT) = ZOUTW
          XOUTE(NOUT) = ZOUTE
          NOUTD(NOUT) = IOUTD
        ENDIF
      ENDDO CORR
      ENDIF
 2300 IF (NOUT > 0) THEN
        WRITE (IU06,'(/4X," AREAS TO BE CORRECTED IN OUTPUT GRID",      &
     &   /,4X,"  NO.   SOUTHERN LAT ",                                  &
     &   " NORTHERN LAT  WESTERN LONG ",                                &
     &   " EASTERN LONG  DEPTH")')
        DO I=1,NOUT
          CALL ADJUST (XOUTW(I), XOUTE(I))
          WRITE (IU06,'(4X,I5,1X,4F14.3,I7 )') I, XOUTS(I),             &
     &     XOUTN(I), XOUTW(I), XOUTE(I), NOUTD(I)
        ENDDO
      ENDIF


! ----------------------------------------------------------------------

!*    2.5 MODEL OPTIONS.
!         --------------

      WRITE (IU06,'(" OUTPUT OPTION IS       IFORM =",I3,               &
     &         " (1: UNFORM.  2: FORM.  3: BOTH)")') IFORM
      WRITE (IU06,'(" TEST OUTPUT OPTION IS  ITEST =",I3,               &
     &         " (0: NO  >0: UP TO LEVEL ITEST)")') ITEST
      WRITE (IU06,'(" BLOCK TEST OPTION IS  ITESTB =",I3,               &
     &         " (0: NO  >0: UP TO BLOCK ITESTB)")') ITESTB

! ----------------------------------------------------------------------

!*    2.6 NESTED GRID INFORMATION.
!         ------------------------

      IF ((IBOUNC == 1 .OR. IBOUNF == 1) .AND. IRGG == 1) THEN
        WRITE (IU06,*) '**********************************************'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '*       FATAL ERROR IN SUB. UIPREP           *'
        WRITE (IU06,*) '*       ==========================           *'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '* THE NESTING OPTION DOES NOT WORK FOR       *'
        WRITE (IU06,*) '* IRREGULAR ALT-LON GRID                     *' 
        WRITE (IU06,*) '* IRGG = ',IRGG
        WRITE (IU06,*) '* IBOUNC = ',IBOUNC
        WRITE (IU06,*) '* IBOUNF = ',IBOUNF
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '**********************************************'
        CALL WAM_ABORT(__FILENAME__,__LINE__)
      ENDIF


!*    2.6.1 COARSE GRID OPTION.
!           -------------------

      WRITE (IU06,'(" COARSE GRID OPTION IS IBOUNC = ",I3)') IBOUNC

      IF (IBOUNC == 1) THEN

        DO I=1,GBOUNC
          CALL ADJUST (AMOWEC(I), AMOEAC(I))
          WRITE (IU06,*) '   THIS IS A COARSE GRID RUN  INFORMATION',   &
     &     ' FOR A FOLLOW UP FINE GRID WILL BE GENERATED'
          WRITE (IU06,'(/4X," NEST AREA IN COARSE GIRD IS",             &
     &     /,4X,"  SOUTHERN LAT  NORTHERN LAT ",                        &
     &     " WESTERN LONG  EASTERN LONG ")')
          WRITE (IU06,'(4X,4F14.3)') AMOSOC(I), AMONOC(I), AMOWEC(I),   &
     &     AMOEAC(I)

!*    2.6.1.1 ARE ALL CORNER POINTS OF THE NEST GRID POINTS?
!             ----------------------------------------------

          WEST = MOD(AMOWEC(I) - AMOWEP + 720.0_JWRB, 360.0_JWRB)
          EAST = MOD(AMOEAC(I) - AMOEAP + 720.0_JWRB, 360.0_JWRB)
          DW=ABS(NINT(WEST/ XDELLO)-(WEST/ XDELLO)) 
          DE=ABS(NINT(EAST/ XDELLO)-(EAST/ XDELLO))
          DS=ABS(NINT((AMOSOC(I) - AMOSOP)/ XDELLA)-((AMOSOC(I)-AMOSOP)/ XDELLA))
          DN=ABS(NINT((AMONOC(I) - AMONOP)/ XDELLA)-((AMONOC(I)-AMONOP)/ XDELLA))

          IF ((DW > 1.E-10_JWRB) .OR.                                &
     &        (DE > 1.E-10_JWRB) .OR.                                &
     &        (DS > 1.E-10_JWRB) .OR.                                &
     &        (DN > 1.E-10_JWRB)) THEN
            WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '+    WARNING ERROR IN SUB. UIPREP        +'
            WRITE (IU06,*) '+    ============================        +'
            WRITE (IU06,*) '+ ERROR IN NEST SPECIFICATIONS.          +'
            WRITE (IU06,*) '+ ONE OR MORE CORNER POINTS ARE NOT      +'
            WRITE (IU06,*) '+ COARSE GRID POINTS.                    +'
            WRITE (IU06,*) '+ NEST INFORMATION WILL NOT BE GENERATED +'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++'
            WRITE (IU06,*) DW, AMOWEC(I), AMOWEP
            WRITE (IU06,*) DE, AMOWEC(I), AMOEAP
            WRITE (IU06,*) DS, AMOSOC(I), AMOSOP
            WRITE (IU06,*) DN, AMONOC(I), AMONOP
            IBOUNC = 0
            CALL WAM_ABORT(__FILENAME__,__LINE__)
          ENDIF


!*    2.6.1.2 INCLUDES THE COARSE GRID THE NEST GRID?
!             ---------------------------------------

          IF ((IPER /= 1) .AND.                                                   &
     &     (AMOWEP > AMOWEC(I)           .OR. AMOEAC(I) > AMOEAP     ) .AND.      &
     &     (AMOWEP > AMOWEC(I)+360._JWRB .OR. AMOEAC(I)+360._JWRB > AMOEAP) .AND. &
     &     (AMOWEP > AMOWEC(I)-360._JWRB .OR. AMOEAC(I)-360._JWRB > AMOEAP)) THEN

            WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '+     WARNING ERROR IN SUB. UIPREP       +'
            WRITE (IU06,*) '+     ============================       +'
            WRITE (IU06,*) '+ ERROR IN NEST SPECIFICATIONS.          +'
            WRITE (IU06,*) '+ WEST OR EAST BOUNDARY IS NOT IN COARSE +'
            WRITE (IU06,*) '+ GRID AREA.                             +'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '+     AMOWEC AND/OR AMOEAC ARE WRONG     +'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '+ NEST INFORMATION WILL NOT BE GENERATED +'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++'
            IBOUNC = 0
          ENDIF
          IF (AMOSOP > AMOSOC(I) .OR. AMONOC(I) > AMONOP .OR.     &
     &        AMOSOC(I) >= AMONOC(I)) THEN
            WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '+     WARNING ERROR IN SUB. UIPREP       +'
            WRITE (IU06,*) '+     ============================       +'
            WRITE (IU06,*) '+ ERROR IN NEST SPECIFICATIONS.          +'
            WRITE (IU06,*) '+ OF FINE GRID NUMBER : ', I
            WRITE (IU06,*) '+ NORTH OR SOUTH BOUNDARY IS NOT IN      +'
            WRITE (IU06,*) '+ COARSE GRID AREA, OR SOUTH IS GE NORTH +'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '+     AMOSOC AND/OR AMONOC ARE WRONG     +'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '+ NEST INFORMATION WILL NOT BE GENERATED +'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++'
            IBOUNC = 0
          ENDIF
        ENDDO
        WRITE (IU06,*) '   NUMBERS OF FINE GRIDS ',GBOUNC
      ELSE
        WRITE (IU06,*) '   A NEST IS NOT INCLUDED IN THIS GRID'
      ENDIF

!*    2.6.2 FINE GRID OPTION.
!           -----------------

      WRITE (IU06,'(" FINE GRID OPTION IS   IBOUNF = ",I3)') IBOUNF
      IF (IBOUNF == 1) THEN
        WRITE (IU06,*) '   THIS IS A FINE GRID RUN, INPUT FROM',        &
     &   ' A COARSE GRID IS EXPECTED'
      ELSE
        WRITE (IU06,*) '   BOUNDARY VALUES FROM A COARSE GRID',         &
     &   ' ARE NOT EXPECTED'
      ENDIF

      WRITE (IU06,'(" END OF USER INPUT PROTOCOL",/)')

      RETURN

! ----------------------------------------------------------------------


!*    3. ERROR HANDLING.
!        ---------------

 3000 CONTINUE
      WRITE (IU06,*) '**********************************************'
      WRITE (IU06,*) '*                                            *'
      WRITE (IU06,*) '*       FATAL ERROR IN SUB. UIPREP           *'
      WRITE (IU06,*) '*       ==========================           *'
      WRITE (IU06,*) '*                                            *'
      WRITE (IU06,*) '*  ERROR WHILE READING NAMELIST  STATUS=', IOS
      WRITE (IU06,*) '*                                            *'
      WRITE (IU06,*) '**********************************************'
      CALL WAM_ABORT(__FILENAME__,__LINE__)

      END SUBROUTINE UIPREP
