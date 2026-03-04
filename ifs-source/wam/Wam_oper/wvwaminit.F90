! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#define __FILENAME__ "wvwaminit.F90"

      SUBROUTINE WVWAMINIT (LLCOUPLED, IULOG, LLRNL,                    &
     &                      NGAUSSW, NLON, NLAT, RSOUTW, RNORTW)

!****  *WVWAMINIT* - INITIALISES A FEW VARIABLE AND READ WAVE MODEL
!                    NAMELIST (if LLRNL is true) and
!                    CONFIGURATION PARAMETERS AND TABLES.


!     J. BIDLOT    ECMWF   FEBRUARY 2008: MOVE THIS PART OUT OF MPDECOMP 

!*    INTERFACE.
!     ----------
!     CALL *WVWAMINIT* (LLCOUPLED,IULOG,LLRNL,
!                       NGAUSSW, NLON, NLAT, RSOUTW, RNORTW)

!      *LLCOUPLED*  LOGICAL  TRUE IF ATMOSPHERIC COUPLED RUN.
!      *IULOG*      INTEGER  STANDARD OUTPUT UNIT USED BY IFS.
!      *LLRNL*      LOGICAL  IF TRUE READS WAM NAMELIST.
!      *NGAUSSW*    INTEGER  1: GAUSSIAN GRID, 0: LAT-LON GRID
!      *NLON*       INTEGER  MAXIMUM NUMBER OF LONGITUDES USED BY WAM. 
!      *NLAT*       INTEGER  MAXIMUM NUMBER OF LATITUDES USED BY WAM. 
!      *RSOUTW*     REAL     MOST SOUTHERN LATITUDE IN WAM.
!      *RNORTW*     REAL     MOST NORTHERn LATITUDE IN WAM.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMAP   , ONLY : AMOSOP   ,AMONOP   ,IQGAUSS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPRECR   ,NPRECI   ,    &
     &            KTAG
      USE YOWPARAM , ONLY : NGX      ,NGY      ,KWAMVER 
      
      USE YOWTEST  , ONLY : IU06
      USE YOWUNIT  , ONLY : IREADG   ,NPROPAGS ,IU07     ,IU08     ,    &
     &            LWVWAMINIT
      USE YOWSTAT  , ONLY : IPROPAGS
      USE YOWUNPOOL, ONLY : LLUNSTR
      USE YOWABORT, ONLY : WAM_ABORT

      USE MPL_MODULE, ONLY : MPL_MYRANK, MPL_NPROC, MPL_BARRIER
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "expand_string.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "mpuserin.intfb.h"
#include "setwavphys.intfb.h"
#include "readpre.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IULOG
      INTEGER(KIND=JWIM), INTENT(OUT) :: NGAUSSW, NLON, NLAT
      INTEGER(KIND=JWIM) :: IREAD, LFILE
      INTEGER(KIND=JWIM) :: I4(2)

      REAL(KIND=JWRB), INTENT(OUT) :: RSOUTW, RNORTW
      REAL(KIND=JWRB) :: X4(2)
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      CHARACTER(LEN=1) :: C1 
      CHARACTER(LEN=80) :: FILENAME
      CHARACTER(LEN=80) :: LOGFILENAME

      LOGICAL, INTENT(IN) :: LLRNL
      LOGICAL, INTENT(IN) :: LLCOUPLED
      LOGICAL :: LLEXIST
      LOGICAL, SAVE :: LFRST

      DATA LFRST /.TRUE./

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WVWAMINIT',0,ZHOOK_HANDLE)

      LWVWAMINIT=.FALSE.

      IRANK = MPL_MYRANK()
      NPROC = MPL_NPROC()

!     STANDARD OUTPUT UNIT
!     --------------------

      IF (LLCOUPLED) THEN
        IU06=IULOG
      ELSE
        IU06=65
        LLEXIST=.TRUE.
        DO WHILE (LLEXIST)
          IU06=IU06+1
          LOGFILENAME='logfile.%p'
          CALL EXPAND_STRING(IRANK,NPROC,0,0,LOGFILENAME,1)
          INQUIRE(UNIT=IU06, OPENED=LLEXIST)
          IF (.NOT. LLEXIST) THEN
            OPEN(IU06,FILE=LOGFILENAME,STATUS='UNKNOWN')
            WRITE(IU06,*) ' STDOUT OF PE ', IRANK
            WRITE(IU06,*) ' ============ '
            WRITE(IU06,*) ' '
          ENDIF
        ENDDO
      ENDIF

!     GET CONTROLLING FLAG FROM INPUT NAMELIST 
!     ----------------------------------------
      IF (LFRST .AND. LLRNL) THEN
        CALL MPUSERIN
        LFRST=.FALSE.
      ENDIF

      WRITE(IU06,*) ' WAM SOFTWARE VERSION: ', KWAMVER


!     SET WAVE PHYSICS PACKAGE
!     ------------------------

      CALL SETWAVPHYS

!     DETERMINE BYTE STORAGE REPRESENTATION OF REAL NUMBERS
!     -----------------------------------------------------

      X4=1.0_JWRB
      NPRECR = KIND(X4)
      I4=1
      NPRECI = KIND(I4)

!     READ GRID INPUT FROM PREPROC 
!     ----------------------------

      IREAD=IREADG

!     CONSTANT FILES INPUT UNIT
!     -------------------------
      IF (IRANK == IREAD) THEN
        FILENAME='wam_grid_tables'
        LFILE=0
        LLEXIST=.FALSE.
        IF (FILENAME /= ' ') LFILE=LEN_TRIM(FILENAME)
        INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
        IF (.NOT. LLEXIST) THEN
          WRITE(IU06,*) '************************************'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '*  FATAL ERROR IN SUB. WVWAMINIT   *'
          WRITE(IU06,*) '*  =============================   *'
          WRITE(IU06,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE), ' IS MISSING !!!!'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '************************************'
          CALL WAM_ABORT("WAVE MODEL INPUT FILE '"//FILENAME(1:LFILE)//"' IS MISSING !!!",&
                        &__FILENAME__,__LINE__)
        ENDIF
        IU07 = IWAM_GET_UNIT(IU06, FILENAME(1:LFILE) , 'r', 'u',0,'READWRITE')

        IF (IPROPAGS < 0 .OR. IPROPAGS > NPROPAGS) THEN
          WRITE(IU06,*) '************************************'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '*  FATAL ERROR IN SUB. WVWAMINIT   *'
          WRITE(IU06,*) '*  WRONG VALUE FOR IPROPAGS:       *'
          WRITE(IU06,*) '*  IPROPAGS= ',IPROPAGS
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '************************************'
          CALL WAM_ABORT("Wrong value for IPROPAGS",__FILENAME__,__LINE__)
        ENDIF

        WRITE(C1,'(I1)') IPROPAGS 
        FILENAME='wam_subgrid_'//C1

        LFILE=0
        LLEXIST=.FALSE.
        IF (FILENAME /= ' ') LFILE=LEN_TRIM(FILENAME)
        INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
        IF (.NOT. LLEXIST) THEN
          WRITE(IU06,*) '************************************'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '*  FATAL ERROR IN SUB. WVWAMINIT   *'
          WRITE(IU06,*) '*  =============================   *'
          WRITE(IU06,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE), ' IS MISSING !!!!'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '************************************'
          CALL WAM_ABORT("WAVE MODEL INPUT FILE '"//FILENAME(1:LFILE)//"' IS MISSING !!!",&
                        &__FILENAME__,__LINE__)
        ENDIF
        IU08(IPROPAGS) = IWAM_GET_UNIT(IU06, FILENAME(1:LFILE),'r','u',0,'READWRITE')
      ENDIF

      KTAG=1

      CALL MPL_BARRIER(CDSTRING='WVWAMINIT:')

      CALL READPRE (IU07)
      CALL MPL_BARRIER(CDSTRING='WVWAMINIT:')

      WRITE(IU06,*) ' WVWAMINT: WAVE MODEL CONFIGURATION READ IN'
      CALL FLUSH (IU06)

!     RETURN NECESSARY PARAMETERS
      NGAUSSW = IQGAUSS
      NLON=NGX
      NLAT=NGY
      RSOUTW=AMOSOP
      RNORTW=AMONOP

      IF (LHOOK) CALL DR_HOOK('WVWAMINIT',1,ZHOOK_HANDLE)

      END SUBROUTINE WVWAMINIT
