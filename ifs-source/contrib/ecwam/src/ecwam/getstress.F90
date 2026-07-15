! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE GETSTRESS(BLK2LOC, WVENVI, FF_NOW, NEMO2WAM,        &
     &                     NBLKS, NBLKE, IREAD)

! ----------------------------------------------------------------------

!*    PURPOSE.
!     --------
!     READS RESTART WIND AND STRESS FIELDS FROM DISK OR USE GRIB WIND
!     AND DRAG COEFFICIENT TO RECONSTRUCT ALL WIND AND STRESS FIELDS.

!**   INTERFACE.
!     ----------
!     *CALL**GETSTRESS(BLK2LOC, WVENVI, FF_NOW, NBLKS, NBLKE, IREAD)
!     *BLK2LOC*   POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!     *WVENVI*    WAVE ENVIRONMENT
!     *FF_NOW*    FORCING FIELDS AT CURRENT TIME
!     *NEMO2WAM*  FIELDS FRON OCEAN MODEL to WAM
!     *NBLKS*     INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NBLKE*     INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *IREAD*     PROCESSOR WHICH WILL ACCESS THE FILE ON DISK

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDLOC, ENVIRONMENT, FORCING_FIELDS, OCEAN2WAVE

      USE YOWCOUT  , ONLY : NREAL    ,LRSTPARALR
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, KIJL4CHNK, IJFROMCHNK
      USE YOWMAP   , ONLY : NIBLO
      USE YOWMESPAS, ONLY : LGRIBIN  ,LWAVEWIND
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NSUP     ,KTAG
      USE YOWREFD  , ONLY : LLUPDTTD
      USE YOWPHYS  , ONLY : PRCHAR 
      USE YOWSTAT  , ONLY : CDATEA   ,CDATEF   ,CDTPRO   ,IREFRA   ,LNSESTART
      USE YOWTEST  , ONLY : IU06
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,LRESTARTED
      USE YOWUBUF  , ONLY : LUPDTWGHT
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_BROADCAST

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "buildstress.intfb.h"
#include "expand_string.intfb.h"
#include "grstname.intfb.h"
#include "mpdistribscfld.intfb.h"
#include "readstress.intfb.h"

      TYPE(WVGRIDLOC), INTENT(IN) :: BLK2LOC
      TYPE(ENVIRONMENT), INTENT(INOUT) :: WVENVI
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NOW
      TYPE(OCEAN2WAVE), INTENT(INOUT) :: NEMO2WAM
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      INTEGER(KIND=JWIM),DIMENSION(NPROC), INTENT(IN) :: NBLKS, NBLKE


      INTEGER(KIND=JWIM) :: IBUFLENGTH
      INTEGER(KIND=JWIM) :: IFLD, IJ, KCOUNT, IC
      INTEGER(KIND=JWIM) :: IJSG, IJLG, ICHNK, KIJS, IJSB, KIJL, IJLB
      INTEGER(KIND=JWIM) :: IJINF, IJSUP
      INTEGER(KIND=JWIM) :: LNAME 
      INTEGER(KIND=JWIM) :: IFCST
      INTEGER(KIND=JWIM),ALLOCATABLE,DIMENSION(:) :: IBUF

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:,:) :: RFIELD

      CHARACTER(LEN= 14) :: ZERO
      CHARACTER(LEN=296) :: FILENAME

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('GETSTRESS',0,ZHOOK_HANDLE)


      ZERO = ' '

      IJSG = IJFROMCHNK(1,1)
      IJLG = IJSG + SUM(KIJL4CHNK) - 1

      DO ICHNK=1,NCHNK
        FF_NOW%Z0B(:,ICHNK) = 0.0_JWRB
      ENDDO

!     READ RESTART FILE FROM PE IREAD

      IF (LGRIBIN .AND. .NOT.LRESTARTED) THEN
!       GRIB RESTART
!       CREATES WIND AND STRESS FIELDS FROM GRIB WINDS AND DRAG COEFFICIENT.
        CALL BUILDSTRESS(BLK2LOC, WVENVI, FF_NOW, NEMO2WAM, IREAD)

      ELSE

!       BINARY RESTART

        IF (LRSTPARALR) THEN
          IJINF = IJSG
          IJSUP = IJLG
        ELSE
          IJINF = 1
          IJSUP = NIBLO
        ENDIF

        ALLOCATE(RFIELD(IJINF:IJSUP,NREAL))

        IFCST = 0
        CALL GRSTNAME(CDATEA, CDATEF, IFCST, 'LAW', ICPLEN, CPATH, FILENAME)
        IF (LRSTPARALR) THEN
!          RESTART FILES FROM ALL PE's
           LNAME = LEN_TRIM(FILENAME)
           FILENAME = FILENAME(1:LNAME)//'.%p_%n'
           CALL EXPAND_STRING(IRANK, NPROC, 0, 0, FILENAME, 1)
        ENDIF

        IF (LRSTPARALR .OR. IRANK == IREAD) THEN
          CALL READSTRESS(IJINF, IJSUP, NREAL, RFIELD, FILENAME, LRSTPARALR)
        ENDIF

        IF (.NOT.LRSTPARALR) THEN
!       DISTRIBUTE THE DIFFERENT CONTRIBUTIONS TO THE OTHER PE'S
          DO IFLD = 1, NREAL
            CALL MPDISTRIBSCFLD(IREAD, KTAG, NBLKS, NBLKE, RFIELD(:,IFLD))
            KTAG=KTAG+1
          ENDDO
        ENDIF

!       KEEP CORRESPONDING CONTRIBUTION 
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, IJSB, KIJL, IJLB)
        DO ICHNK = 1, NCHNK
          KIJS = 1
          IJSB = IJFROMCHNK(KIJS, ICHNK)
          KIJL = KIJL4CHNK(ICHNK)
          IJLB = IJFROMCHNK(KIJL, ICHNK)

          FF_NOW%WSWAVE(KIJS:KIJL, ICHNK)  = RFIELD(IJSB:IJLB,1)
          FF_NOW%WDWAVE(KIJS:KIJL, ICHNK)  = RFIELD(IJSB:IJLB,2)
          FF_NOW%UFRIC(KIJS:KIJL, ICHNK)   = RFIELD(IJSB:IJLB,3)
          FF_NOW%TAUW(KIJS:KIJL, ICHNK)    = RFIELD(IJSB:IJLB,4)
          FF_NOW%TAUWDIR(KIJS:KIJL, ICHNK) = RFIELD(IJSB:IJLB,5)
          FF_NOW%Z0M(KIJS:KIJL, ICHNK)     = RFIELD(IJSB:IJLB,6)
          FF_NOW%Z0B(KIJS:KIJL, ICHNK)     = RFIELD(IJSB:IJLB,7)
          FF_NOW%CHRNCK(KIJS:KIJL, ICHNK)  = RFIELD(IJSB:IJLB,8)
          FF_NOW%AIRD(KIJS:KIJL, ICHNK)    = RFIELD(IJSB:IJLB,9)
          FF_NOW%WSTAR(KIJS:KIJL, ICHNK)   = RFIELD(IJSB:IJLB,10)
          FF_NOW%CICOVER(KIJS:KIJL, ICHNK) = RFIELD(IJSB:IJLB,11)
          FF_NOW%CITHICK(KIJS:KIJL, ICHNK) = RFIELD(IJSB:IJLB,12)
          FF_NOW%USTRA(KIJS:KIJL, ICHNK)   = RFIELD(IJSB:IJLB,13)
          FF_NOW%VSTRA(KIJS:KIJL, ICHNK)   = RFIELD(IJSB:IJLB,14)
          WVENVI%UCUR(KIJS:KIJL, ICHNK)    = RFIELD(IJSB:IJLB,15)
          WVENVI%VCUR(KIJS:KIJL, ICHNK)    = RFIELD(IJSB:IJLB,16)

          IF (KIJL < NPROMA_WAM) THEN
            FF_NOW%WSWAVE(KIJL+1:NPROMA_WAM, ICHNK)  = FF_NOW%WSWAVE(1, ICHNK)
            FF_NOW%WDWAVE(KIJL+1:NPROMA_WAM, ICHNK)  = FF_NOW%WDWAVE(1, ICHNK)
            FF_NOW%UFRIC(KIJL+1:NPROMA_WAM, ICHNK)   = FF_NOW%UFRIC(1, ICHNK)
            FF_NOW%TAUW(KIJL+1:NPROMA_WAM, ICHNK)    = FF_NOW%TAUW(1, ICHNK)
            FF_NOW%TAUWDIR(KIJL+1:NPROMA_WAM, ICHNK) = FF_NOW%TAUWDIR(1, ICHNK)
            FF_NOW%Z0M(KIJL+1:NPROMA_WAM, ICHNK)     = FF_NOW%Z0M(1, ICHNK)
            FF_NOW%Z0B(KIJL+1:NPROMA_WAM, ICHNK)     = FF_NOW%Z0B(1, ICHNK)
            FF_NOW%CHRNCK(KIJL+1:NPROMA_WAM, ICHNK)  = FF_NOW%CHRNCK(1, ICHNK)
            FF_NOW%AIRD(KIJL+1:NPROMA_WAM, ICHNK)    = FF_NOW%AIRD(1, ICHNK)
            FF_NOW%WSTAR(KIJL+1:NPROMA_WAM, ICHNK)   = FF_NOW%WSTAR(1, ICHNK)
            FF_NOW%CICOVER(KIJL+1:NPROMA_WAM, ICHNK) = FF_NOW%CICOVER(1, ICHNK)
            FF_NOW%CITHICK(KIJL+1:NPROMA_WAM, ICHNK) = FF_NOW%CITHICK(1, ICHNK)
            FF_NOW%USTRA(KIJL+1:NPROMA_WAM, ICHNK)   = FF_NOW%USTRA(1, ICHNK)
            FF_NOW%VSTRA(KIJL+1:NPROMA_WAM, ICHNK)   = FF_NOW%VSTRA(1, ICHNK)
            WVENVI%UCUR(KIJL+1:NPROMA_WAM, ICHNK)    = WVENVI%UCUR(1,ICHNK)
            WVENVI%VCUR(KIJL+1:NPROMA_WAM, ICHNK)    = WVENVI%VCUR(1,ICHNK)
          ENDIF
        ENDDO
!$OMP   END PARALLEL DO

        IF (IREFRA /= 0) THEN
          LLUPDTTD = .TRUE.
          LUPDTWGHT=.TRUE.
        ENDIF

        DEALLOCATE(RFIELD)

      ENDIF

!     BROADCAST THE 4 DATES FROM RESTART FILE TO THE OTHER PE'S
!     AS WELL AS LWAVEWIND
      IF (.NOT. (IREAD == 0 .OR. NPROC == 1)) THEN
        IBUFLENGTH = 57+1
        ALLOCATE(IBUF(IBUFLENGTH))
        IF (IRANK == IREAD) THEN
          IF (LWAVEWIND) THEN
            IBUF(1)=1
          ELSE
            IBUF(1)=0
          ENDIF
          KCOUNT=1
          DO IC=1,14
            KCOUNT=KCOUNT+1
            READ(CDTPRO(IC:IC),'(I1)') IBUF(KCOUNT)
          ENDDO
          DO IC=1,14
            KCOUNT=KCOUNT+1
            READ(CDATEWO(IC:IC),'(I1)') IBUF(KCOUNT)
          ENDDO
          DO IC=1,14
            KCOUNT=KCOUNT+1
            READ(CDAWIFL(IC:IC),'(I1)') IBUF(KCOUNT)
          ENDDO
          DO IC=1,14
            KCOUNT=KCOUNT+1
            READ(CDATEFL(IC:IC),'(I1)') IBUF(KCOUNT)
            IF (KCOUNT > IBUFLENGTH) THEN
              WRITE (IU06,*) ' '
              WRITE (IU06,*) ' **************************************'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' *   FATAL ERROR IN SUB. GETSTRESS:   *'
              WRITE (IU06,*) ' *   ==============================   *'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' * IBUFLENGTH TOO SMALL !!!           *'
              WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS   *'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' **************************************'
              CALL ABORT1
            ENDIF
          ENDDO
        ENDIF

        CALL GSTATS(621,0)
        CALL MPL_BROADCAST(IBUF(1:IBUFLENGTH), KROOT=IREAD, KTAG=KTAG, CDSTRING='GETSTRESS:')
        CALL GSTATS(621,1)

        KTAG=KTAG+1

        IF (IRANK /= IREAD) THEN
          IF (IBUF(1) == 1) THEN
            LWAVEWIND=.TRUE. 
          ELSE
            LWAVEWIND=.FALSE. 
          ENDIF
          KCOUNT=1
          DO IC=1,14
            KCOUNT=KCOUNT+1
            WRITE(CDTPRO(IC:IC),'(I1)') IBUF(KCOUNT)
          ENDDO
          DO IC=1,14
            KCOUNT=KCOUNT+1
            WRITE(CDATEWO(IC:IC),'(I1)') IBUF(KCOUNT)
          ENDDO
          DO IC=1,14
            KCOUNT=KCOUNT+1
            WRITE(CDAWIFL(IC:IC),'(I1)') IBUF(KCOUNT)
          ENDDO
          DO IC=1,14
            KCOUNT=KCOUNT+1
            WRITE(CDATEFL(IC:IC),'(I1)') IBUF(KCOUNT)
            IF (KCOUNT > IBUFLENGTH) THEN
              WRITE (IU06,*) ' '
              WRITE (IU06,*) ' **************************************'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' *   FATAL ERROR IN SUB. GETSTRESS:   *'
              WRITE (IU06,*) ' *   ==============================   *'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' * IBUFLENGTH TOO SMALL !!!           *'
              WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS   *'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' **************************************'
              CALL ABORT1
            ENDIF
          ENDDO
        ENDIF

        DEALLOCATE(IBUF)

      ENDIF

      IF (CDTPRO == '00000000000000') CDTPRO = ZERO
      IF (CDATEWO == '00000000000000')CDATEWO = ZERO
      IF (CDAWIFL == '00000000000000')CDAWIFL = ZERO
      IF (CDATEFL == '00000000000000') CDATEFL= ZERO


      IF (LNSESTART .AND. .NOT.LRESTARTED) THEN
!       WHEN INITAL SPECTRA SET TO NOISE LEVEL,
!       RESET WAVE INDUCED STRESS TO ZERO 
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK)
        DO ICHNK = 1, NCHNK
          FF_NOW%CHRNCK(1:NPROMA_WAM,ICHNK) = PRCHAR
          FF_NOW%TAUW(1:NPROMA_WAM,ICHNK) = 0.0_JWRB
        ENDDO
!$OMP   END PARALLEL DO
      ENDIF


      WRITE(IU06,*) ''
      WRITE(IU06,*) ' WIND AND STRESS FILES READ IN....... CDTPRO  = ', CDTPRO
      CALL FLUSH (IU06)

IF (LHOOK) CALL DR_HOOK('GETSTRESS',1,ZHOOK_HANDLE)

END SUBROUTINE GETSTRESS
