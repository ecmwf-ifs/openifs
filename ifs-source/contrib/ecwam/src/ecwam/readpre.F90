! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE READPRE (IU07)

! ----------------------------------------------------------------------

!**** *READPRE*  READ GRID OUTPUT FROM PREPROC.

!     H. GUNTHER      GKSS/ECMWF     MAY 1990
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     J. BIDLOT            ECMWF       11/2003
!                          IF YOUR ARE RUNNING AT ECMWF:
!                          BE AWARE THAT IF YOU CHANGE ANYTHING TO THE
!                          STRUCTURE OF THE OUTPUT FILE YOU WILL HAVE TO
!                          MAKE SURE THAT IT IS CREATED FOR YOUR RUN,
!                          OTHERWISE IT MIGHT PICK UP THE DEFAULT ONE
!                          THAT IS ALREADY ON DISK.
!                          YOU ALSO HAVE TO CHANGE OUTCOM.
! ALSO CHECK
!            READPREB in Alt
!  and
!            ... /nemo/tools/interpolate/wambingrid.F90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*    PURPOSE.
!     --------

!       INPUT OF PREPROC GRID OUTPUT.

!**   INTERFACE.
!     ----------

!       *CALL* *READPRE (IU07)*
!          *IU07 *  - INPUT UNIT OF PREPROC GRID FILE.

!     METHOD.
!     -------

!       UNFORMATED READ FROM UNIT.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWGRIBHD, ONLY : IMDLGRBID_G
      USE YOWMAP   , ONLY : NGX      ,NGY      ,    &
     &            IPER     ,IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,    &
     &            AMONOP   ,XDELLA   ,XDELLO   ,NLONRGG  ,    &
     &            IQGAUSS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,KTAG
      USE YOWPARAM , ONLY : LLR8TOR4 ,LLUNSTR
      USE YOWSHAL  , ONLY : BATHY
      USE YOWTEST  , ONLY : IU06
      USE YOWABORT , ONLY : WAM_ABORT
#ifdef WAM_HAVE_UNWAM
      USE YOWUNPOOL, ONLY : LPREPROC
      USE UNWAM     ,ONLY : INIT_UNWAM, UNWAM_IN, SET_UNWAM_HANDLES
#endif
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "mpbcastgrid.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU07
      INTEGER(KIND=JWIM) :: IREAD
      INTEGER(KIND=JWIM) :: IP, I, K
      INTEGER(KIND=JWIM) :: KMDLGRDID
      INTEGER(KIND=JWIM) :: NKIND !Precision of file when reading

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('READPRE',0,ZHOOK_HANDLE)

      NKIND=0
      IREAD=1
      KTAG=1

      CALL GSTATS(1771,0)
      IF (IRANK == IREAD) THEN
!       READ MODEL IDENTIFIERS
        CALL READREC(1)
        IF (KMDLGRDID /= IMDLGRBID_G) THEN
          WRITE(IU06,*) '*****************************************'
          WRITE(IU06,*) '*                                       *'
          WRITE(IU06,*) '*  FATAL ERROR(S) IN SUB. READPRE       *'
          WRITE(IU06,*) '*  ==============================       *'
          WRITE(IU06,*) '*                                       *'
          WRITE(IU06,*) '* THE PROGRAM HAS DETECTED DIFFERENT    *'
          WRITE(IU06,*) '* MODEL GRIB IDENTIFIER.                *' 
          WRITE(IU06,*) '* MAKE SURE YOU HAVE RUN PREPROC !!!!   *'
          WRITE(IU06,*)    KMDLGRDID, IMDLGRBID_G 
          WRITE(IU06,*) '*                                       *'
          WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.     *'
          WRITE(IU06,*) '* ---------------   --------------      *'
          WRITE(IU06,*) '*****************************************'
          CALL ABORT1
        ENDIF
        IF (NKIND /= KIND(AMOSOP)) THEN
          WRITE(IU06,*) '*****************************************'
          WRITE(IU06,*) '*                                       *'
          WRITE(IU06,*) '*  FATAL ERROR(S) IN SUB. READPRE       *'
          WRITE(IU06,*) '*  ==============================       *'
          WRITE(IU06,*) '*                                       *'
          WRITE(IU06,*) '* THE PROGRAM HAS DETECTED DIFFERENT    *'
          WRITE(IU06,*) '* PRECISION IN FILE AND MODEL.          *' 
          WRITE(IU06,*)    NKIND, KIND(AMOSOP)
          WRITE(IU06,*) '*                                       *'
          WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.     *'
          WRITE(IU06,*) '* ---------------   --------------      *'
          WRITE(IU06,*) '*****************************************'
          CALL ABORT1
        ENDIF

!*    0. READ YOWPARAM (BLOCK SIZES). 
!        ----------------------------

        CALL READREC(2)


!*    2. READ MODULE YOWGRID (GENERAL GRID ORGANISATION).
!        ------------------------------------------------

        IF (.NOT.ALLOCATED(NLONRGG)) ALLOCATE(NLONRGG(NGY))

        CALL READREC(3)


!*    3. READ MODULE YOWMAP (LONG. AND LAT. INDICES OF GRID POINTS).
!        --------------------------------------------------------

        CALL READREC(4)

!     DETERMINE IF WE ARE USING A QUASI GAUSSIAN GRID OR 
!     LAT-LONG GRID (REGULAR OR IRREGULAR).

        IF (IPER ==1 .AND. AMONOP == ABS(AMOSOP) .AND.                  &
     &     MOD(NGY,2) == 0 .AND. IRGG == 1 ) THEN
          IQGAUSS=1
        ELSE
          IQGAUSS=0
        ENDIF


!*    8. READ MODULE YOWSHAL (DEPTH AND SHALLOW WATER TABLES).
!        ----------------------------------------------------

        IF (ALLOCATED(BATHY)) DEALLOCATE(BATHY)
        ALLOCATE(BATHY(NGX,NGY))

        CALL READREC(5)

!       THE UNSTRUCTURED BITS (if pre-computed by PREPROC)
        IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
          IF (LPREPROC) THEN
            CALL SET_UNWAM_HANDLES
            CALL UNWAM_IN(IU07)
          ENDIF
#else
          CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
        END IF

      ENDIF

      CLOSE (UNIT=IU07)

      CALL GSTATS(1771,1)



!     SEND INFORMATION FROM READPRE TO ALL PE's
      CALL GSTATS(694,0)
      CALL MPBCASTGRID(IU06,IREAD,KTAG)
      CALL GSTATS(694,1)


!     RECALCULATE THE UNSTRUCTURED BITS (if not read in)

      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
        IF (.NOT.LPREPROC) THEN
          CALL INIT_UNWAM
        ENDIF
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ENDIF

      IF (LHOOK) CALL DR_HOOK('READPRE',1,ZHOOK_HANDLE)

      RETURN

      CONTAINS

      SUBROUTINE READREC(KREC)
      IMPLICIT NONE
      INTEGER(KIND=JWIM), INTENT(IN) :: KREC
      INTEGER(KIND=JWIM) :: ISTAT
      REAL(KIND=JWRU) :: R8_AMOEAP,R8_AMONOP,R8_AMOSOP,R8_AMOWEP,R8_XDELLA,R8_XDELLO
      REAL(KIND=JWRU), ALLOCATABLE, DIMENSION(:,:) :: R8_BATHY

!23456789-123456789-123456789-123456789-123456789-123456789-123456789-12
      SELECT CASE(KREC)
      CASE(1)
         LLR8TOR4 = .FALSE.
         READ(IU07,IOSTAT=ISTAT) NKIND, KMDLGRDID
         IF (ISTAT /= 0) GOTO 1000

         IF (KIND(AMOSOP) == 4 .AND. NKIND /= 4) THEN
            LLR8TOR4 = .TRUE.   ! Input REALs are indeed legacy REAL*8, but KIND(AMOSOP) == 4
         ENDIF
         NKIND = KIND(AMOSOP)   ! Pretend NKIND is in "right precision"
         WRITE(IU06,1002) 'READPRE(READREC): NKIND=',NKIND, ', KIND(AMOSOP)=',KIND(AMOSOP),', LLR8TOR4=',LLR8TOR4
 1002    FORMAT(2X,A,I0,A,I0,A,L1)

      CASE(2)
         READ(IU07,IOSTAT=ISTAT) NGX, NGY
         IF (ISTAT /= 0) GOTO 1000

      CASE(3)
          READ(IU07,IOSTAT=ISTAT) NLONRGG
          IF (ISTAT /= 0) GOTO 1000

      CASE(4)
         IF (LLR8TOR4) THEN
            READ(IU07,IOSTAT=ISTAT) IPER, IRGG, &
     &                              R8_AMOWEP, R8_AMOSOP, R8_AMOEAP, R8_AMONOP, &
     &                              R8_XDELLA, R8_XDELLO
            IF (ISTAT /= 0) GOTO 1000

            AMOWEP = R8_AMOWEP
            AMOSOP = R8_AMOSOP
            AMOEAP = R8_AMOEAP
            AMONOP = R8_AMONOP
            XDELLA = R8_XDELLA
            XDELLO = R8_XDELLO

         ELSE
            READ(IU07,IOSTAT=ISTAT) IPER, IRGG, &
     &                              AMOWEP, AMOSOP, AMOEAP, AMONOP, &
     &                              XDELLA, XDELLO
            IF (ISTAT /= 0) GOTO 1000
         ENDIF

      CASE(5)
         IF (LLR8TOR4) THEN
            ALLOCATE(R8_BATHY(NGX,NGY))
            READ(IU07,IOSTAT=ISTAT) R8_BATHY
            IF (ISTAT /= 0) GOTO 1000
            BATHY(:,:) = R8_BATHY(:,:)
            DEALLOCATE(R8_BATHY)
         ELSE
            READ(IU07,IOSTAT=ISTAT) BATHY
            IF (ISTAT /= 0) GOTO 1000
         ENDIF

      CASE DEFAULT
         WRITE(IU06,*)'***ERROR IN READREC: INVALID RECORD NUMBER=',KREC
         CALL FLUSH(IU06)
         CALL ABORT1
      END SELECT
      RETURN

 1000 CONTINUE
      WRITE(IU06,1001)'***ERROR IN READREC(',KREC,') : IOSTAT=',ISTAT,', NKIND=',NKIND
 1001 FORMAT(1X,A,I0,A,I0,A,I0)
      WRITE(IU06,*) '*************************************'
      WRITE(IU06,*) '*                                   *'
      WRITE(IU06,*) '*  READ ERROR IN SUB. READREC       *'
      WRITE(IU06,*) '*  ==========================       *'
      WRITE(IU06,*) '*                                   *'
      WRITE(IU06,'(1X,A,I0)') '*  READ ERROR TO UNIT fort.',IU07
      WRITE(IU06,*) '*  IS THE FILE PRESENT ????         *' 
      WRITE(IU06,*) '*                                   *'
      WRITE(IU06,*) '*************************************'
      CALL FLUSH(IU06)
      CALL ABORT1
      RETURN

      END

END SUBROUTINE READPRE
