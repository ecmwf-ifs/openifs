! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE GETCURR(LWCUR, IREAD, BLK2LOC,            &
 &                 NXS, NXE, NYS, NYE, FIELDG,       &
 &                 NEMO2WAM, WVENVI)

!****  *GETCURR* - READS SURFACE CURRENTS FROM FILE (IF UNCOUPLED)
!                  OR EXTRACT THEM FROM FORCING_FIELDS DATA STRUCTURE (IF COUPLED)

!     PURPOSE.
!     --------

!       GET SURFACE CURRENTS AND COMPUTE THE ASSOCIATED
!       REFRACTION TERMS.

!*    INTERFACE.
!     ----------

!     CALL *GETCURR*(LWCUR, IREAD, IFROMIJ, JFROMIJ,
!                    NXS, NXE, NYS, NYE, FIELDG,
!                    NEMOUCUR, NEMOVCUR, UCUR, VCUR)

!     *LWCUR*  - LOGICAL CONTROLLING THE PRESENCE OF MEANINGFUL
!                SURFACE CURRENTS IN ARRAY FORCING_FIELDS DATA STRUCTURE   
!     *IREAD*  - PROCESSOR WHICH WILL ACCESS THE FILE ON DISK IF NEEDED).
!     *IFROMIJ*  POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!     *JFROMIJ*  POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!     *NXS:NXE*  FIRST DIMENSION OF FIELDG
!     *NYS:NYE*  SECOND DIMENSION OF FIELDG
!     *FIELDG* - INPUT FORCING FIELDS ON THE WAVE MODEL GRID
!     *NEMOUCUR* U-COMPONENT OF CURRENT FROM NEMO (if used)
!     *NEMOVCUR* V-COMPONENT OF CURRENT FROM NEMO (if used)
!     *UCUR*   - U-COMPONENT OF THE SURFACE CURRENT
!     *VCUR*   - V-COMPONENT OF THE SURFACE CURRENT


!     METHOD.
!     -------

!     IF UNCOUPLED
!     READS CURRENTS FROM FILE "currents" WHEN IT IS NEEDED
!     IF COUPLED
!     GET CURRENTS AS PROCESSED BY *IFSTOWAM*
!     THEN 
!     SET LLCHKCFLA=.TRUE. TO CHECK ON THE CFL CRITERIA.

! ------------------------------------------------------------------- 

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS, WVGRIDLOC, ENVIRONMENT, OCEAN2WAVE

      USE YOWCOUP  , ONLY : LWCOU, LWNEMOCOUCUR

      USE YOWCURR  , ONLY : CDTCUR   ,IDELCUR  , LLCHKCFLA, CURRENT_MAX
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, KIJL4CHNK
      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : NANG     ,NFRE_RED
      USE YOWREFD  , ONLY : LLUPDTTD
      USE YOWSTAT  , ONLY : CDTPRO   ,IREFRA
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : LUPDTWGHT
      USE YOWWIND  , ONLY : LLNEWCURR 

      USE YOMHOOK  , ONLY  : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_ALLREDUCE

! --------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "current2wam.intfb.h"
#include "incdate.intfb.h"
#include "wamcur.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      LOGICAL, INTENT(IN) :: LWCUR
      TYPE(WVGRIDLOC), INTENT(IN) :: BLK2LOC
      TYPE(ENVIRONMENT), INTENT(INOUT) :: WVENVI
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(IN) :: FIELDG
      TYPE(OCEAN2WAVE), INTENT(IN) :: NEMO2WAM


      INTEGER(KIND=JWIM) :: LIU
      INTEGER(KIND=JWIM) :: ICHNK, KIJS, KIJL, IJ, IX, JY
      INTEGER(KIND=JWIM) :: KUPDATE

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NCHNK) :: OLDUCUR, OLDVCUR

      CHARACTER(LEN=14) :: CDATEIN, CDTNEWCUR
      CHARACTER(LEN=24) :: FILNM

      LOGICAL :: LLCURRENT
      LOGICAL :: LLUPDATE

! --------------------------------------------------------------------- 

      IF (LHOOK) CALL DR_HOOK('GETCURR',0,ZHOOK_HANDLE)

!     1.0  GET NEW CURRENTS IF IT IS REQUIRED.
!          ----------------------------------

      CALL GSTATS(1984,0)


      IF (LLNEWCURR) THEN
        IF ( (LWCOU .AND. LWCUR ) .OR. LWNEMOCOUCUR .OR. IREFRA == 2 .OR. IREFRA == 3) THEN

          CDTNEWCUR=CDTCUR

          IF (.NOT.LWCOU) CALL INCDATE(CDTNEWCUR,IDELCUR/2)

          IF (CDTPRO >= CDTNEWCUR) THEN

            DO ICHNK=1,NCHNK
              OLDUCUR(:,ICHNK) = WVENVI%UCUR(:,ICHNK)
              OLDVCUR(:,ICHNK) = WVENVI%VCUR(:,ICHNK)
            ENDDO
            LLCURRENT=.FALSE.

            CALL INCDATE(CDTCUR,IDELCUR)

            IF (LWCOU) THEN
!             CURRENTS FROM COUPLING INTERFACE 
!             --------------------------------
              IF (LWCUR .AND. .NOT.LWNEMOCOUCUR) THEN

                CALL GSTATS(1444,0)
!$OMP           PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, KIJL)
                DO ICHNK = 1, NCHNK
                  KIJS=1
                  KIJL=NPROMA_WAM
                  CALL WAMCUR (NXS, NXE, NYS, NYE, FIELDG,                                     &
 &                             KIJS, KIJL, BLK2LOC%IFROMIJ(:,ICHNK), BLK2LOC%JFROMIJ(:,ICHNK), &
 &                             WVENVI%UCUR(:,ICHNK), WVENVI%VCUR(:,ICHNK))
                ENDDO
!$OMP           END PARALLEL DO
                CALL GSTATS(1444,1)

!             CURRENTS FROM NEMO
!             ------------------
              ELSEIF (LWNEMOCOUCUR) THEN
                WRITE(IU06,*)' NEMO CURRENTS OBTAINED'!

                CALL GSTATS(1444,0)
!$OMP           PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, IJ, IX, JY)
                DO ICHNK = 1, NCHNK
                  DO IJ = 1, NPROMA_WAM
                    IX = BLK2LOC%IFROMIJ(IJ,ICHNK)
                    JY = BLK2LOC%JFROMIJ(IJ,ICHNK)
                    IF (FIELDG%LKFR(IX,JY) <=  0.0_JWRB ) THEN
!                     if lake cover = 0, we assume open ocean point, then get currents directly from NEMO 
                      WVENVI%UCUR(IJ,ICHNK) = SIGN(MIN(ABS(NEMO2WAM%NEMOUCUR(IJ,ICHNK)),CURRENT_MAX),NEMO2WAM%NEMOUCUR(IJ,ICHNK))
                      WVENVI%VCUR(IJ,ICHNK) = SIGN(MIN(ABS(NEMO2WAM%NEMOVCUR(IJ,ICHNK)),CURRENT_MAX),NEMO2WAM%NEMOVCUR(IJ,ICHNK))
                    ELSE
!                     no currents over lakes and land
                      WVENVI%UCUR(IJ,ICHNK) = 0.0_JWRB
                      WVENVI%VCUR(IJ,ICHNK) = 0.0_JWRB
                    ENDIF
                  ENDDO
                ENDDO
!$OMP           END PARALLEL DO
                CALL GSTATS(1444,1)

              ELSE
                DO ICHNK=1, NCHNK
                  WVENVI%UCUR(:,ICHNK)=0.0_JWRB
                  WVENVI%VCUR(:,ICHNK)=0.0_JWRB
                ENDDO
                WRITE(IU06,*)' '
                WRITE(IU06,*)'    ****************************'
                WRITE(IU06,*)'     IREFRA = ',IREFRA
                WRITE(IU06,*)'     LWCUR IS FALSE !!!!' 
                WRITE(IU06,*)'     CURRENTS ARE SET TO 0. '
                WRITE(IU06,*)'    ****************************'
                WRITE(IU06,*)' '
                CALL FLUSH(IU06)
              ENDIF
            ELSE
!             CURRENTS FROM INPUT FILE
!             ------------------------
              FILNM = 'currents'
              LIU   = LEN_TRIM(FILNM)
              FILNM=FILNM(1:LIU)
              INQUIRE(FILE=FILNM,EXIST=LLCURRENT)

              IF (LLCURRENT) THEN
                WRITE(IU06,*) ' '
                WRITE(IU06,*) ' SUB. GETCURR: GETTING OCEAN CURRENTS',  &
     &                        ' FOR DATE ',CDTCUR
                CALL FLUSH(IU06)

                CALL CURRENT2WAM (FILNM, IREAD, CDATEIN,        &
     &                            BLK2LOC,                      &
     &                            NXS, NXE, NYS, NYE, FIELDG,   &
     &                            WVENVI)
                

                IF (CDATEIN /= CDTCUR) THEN
                WRITE (IU06,*) ' **************************************'
                WRITE (IU06,*) ' *                                    *'
                WRITE (IU06,*) ' * PROBLEM IN GETCURR :               *'
                WRITE (IU06,*) ' * THE REQUESTED DATE FOR THE CURRENTS*'
                WRITE (IU06,*) ' * DOES NOT CORRESPOND TO THE DECODED *'
                WRITE (IU06,*) ' * DATE !!!!                          *'
                WRITE (IU06,*) ' * CDTCUR =',CDTCUR 
                WRITE (IU06,*) ' * CDATEIN=',CDATEIN
                WRITE (IU06,*) ' *                                    *'
                WRITE (IU06,*) ' **************************************'
                CALL ABORT1
                ENDIF
              ELSE
                DO ICHNK=1, NCHNK
                  WVENVI%UCUR(:,ICHNK)=0.0_JWRB
                  WVENVI%VCUR(:,ICHNK)=0.0_JWRB
                ENDDO
                WRITE(IU06,*)' '
                WRITE(IU06,*)'    ****************************'
                WRITE(IU06,*)'     FILE ',FILNM,' NOT FOUND '
                WRITE(IU06,*)'     CURRENTS ARE SET TO 0. '
                WRITE(IU06,*)'    ****************************'
                WRITE(IU06,*)' '
                CALL FLUSH(IU06)
              ENDIF

            ENDIF

!           COMPUTE REFRACTION TERMS
!           ------------------------
!           CHECK IF UPDATE IS NEEDED
            LLUPDATE = .FALSE.
            KUPDATE = 0
            OUT: DO ICHNK = 1, NCHNK
              DO IJ = 1, KIJL4CHNK(ICHNK)
                IF ( WVENVI%UCUR(IJ,ICHNK) /= OLDUCUR(IJ,ICHNK) .OR. WVENVI%VCUR(IJ,ICHNK) /= OLDVCUR(IJ,ICHNK) ) THEN
                  LLUPDATE=.TRUE.
                  KUPDATE = 1
                  EXIT OUT
                ENDIF 
              ENDDO
            ENDDO OUT

            IF (NPROC > 1) THEN
              CALL MPL_ALLREDUCE(KUPDATE,'MAX',CDSTRING='GETCURR KUPDATE:')
              IF (KUPDATE > 0 ) LLUPDATE = .TRUE.
            ENDIF

            IF (LLUPDATE) THEN
              IF (IREFRA /= 0) LLUPDTTD = .TRUE.

              LLCHKCFLA=.TRUE.

!             SET LOGICAL TO RECOMPUTE THE WEIGHTS IN CTUW.
              LUPDTWGHT=.TRUE.

            ELSE
              LLCHKCFLA=.FALSE.
            ENDIF

          ELSE
            LLCHKCFLA=.FALSE.
          ENDIF
        ENDIF
      ELSE
        LLCHKCFLA=.FALSE.
      ENDIF

      CALL GSTATS(1984,1)

      IF (LHOOK) CALL DR_HOOK('GETCURR',1,ZHOOK_HANDLE)

END SUBROUTINE GETCURR 
