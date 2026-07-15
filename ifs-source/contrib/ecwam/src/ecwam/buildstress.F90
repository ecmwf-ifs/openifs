! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE BUILDSTRESS(BLK2LOC, WVENVI, FF_NOW, NEMO2WAM, IREAD)

! ----------------------------------------------------------------------

!*    PURPOSE.
!     --------
!     CREATES WIND AND STRESS FIELDS FROM GRIB WINDS AND CD.

!**   INTERFACE.
!     ----------
!     CALL *BUILDSTRESS*(BLK2LOC, WVENVI, FF_NOW, NEMO2WAM,  IREAD)* 
!     *BLK2LOC*  - POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!     *WVENVI*   - WAVE ENVIRONMENT.
!     *FF_NOW*   - DATA STRUCTURE WITH THE CURRENT FORCING FIELDS
!     *NEMO2WAM* - FIELDS FRON OCEAN MODEL to WAM
!     *IREAD*    - PROCESSOR WHICH WILL ACCESS THE FILE ON DISK

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDLOC,ENVIRONMENT, FORCING_FIELDS, OCEAN2WAVE

      USE YOWCOUP  , ONLY : LWCOU    , LLCAPCHNK,                     &
     &                      LWNEMOCOUCIC, LWNEMOCOUCIT
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK 
      USE YOWMPP   , ONLY : NPROC    ,IRANK    ,KTAG
      USE YOWMESPAS, ONLY : LNOCDIN  ,LWAVEWIND
      USE YOWPCONS , ONLY : G        ,GM1      ,ROAIR    ,EPSUS    ,EPSU10
      USE YOWPHYS  , ONLY : ALPHA    ,XKAPPA   ,XNLEV    ,RNUM
      USE YOWSTAT  , ONLY : CDATEA   ,CDTPRO
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL  ,             &
     &                      NXFFS    ,NXFFE    ,NYFFS    ,NYFFE
      USE MPL_MODULE, ONLY : MPL_BROADCAST
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "cdustarz0.intfb.h"
#include "chnkmin.intfb.h"
#include "getwnd.intfb.h"
#include "init_fieldg.intfb.h"
#include "readwgrib.intfb.h"

      TYPE(WVGRIDLOC), INTENT(IN) :: BLK2LOC
      TYPE(ENVIRONMENT), INTENT(INOUT) :: WVENVI
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NOW
      TYPE(OCEAN2WAVE), INTENT(INOUT) :: NEMO2WAM
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD


      INTEGER(KIND=JWIM) :: ICODE_WND
      INTEGER(KIND=JWIM) :: ILEN, LIU, IPARAM, KZLEVUWAVE, KZLEVCD
      INTEGER(KIND=JWIM) :: IJ, ICHNK, KIJS, KIJL
      INTEGER(KIND=JWIM) :: NWAVEWIND(1)

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: TEMPXNLEV, CDSQRTINV, Z0TOT, USTAR, CHARNOCKOG
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NCHNK) :: CD
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NCHNK) :: ALPHAOG

!     INPUT FORCING FIELDS ON THE WAVE MODEL GRID:
      TYPE(FORCING_FIELDS) :: FIELDG

      CHARACTER(LEN=24) :: FILNM

      LOGICAL :: LLONLYPOS
      LOGICAL :: LWNDFILE, LCLOSEWND
      LOGICAL :: LCR
      LOGICAL :: LLINIALL, LLOCAL

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('BUILDSTRESS',0,ZHOOK_HANDLE)
      CALL FIELDG%ALLOC(NXFFS, NYFFS, NXFFE, NYFFE)


      CDATEWO = ' '
      CDAWIFL = ' '
      CDATEFL = ' '
      CDTPRO = CDATEA
      LCR=.FALSE.

!     GETWND AND READWGRIB REQUIRES FIELDG
      LLINIALL=.TRUE.
      LLOCAL=.TRUE.
      CALL INIT_FIELDG(BLK2LOC, LLINIALL, LLOCAL,         &
     &                 NXFFS, NXFFE, NYFFS, NYFFE, FIELDG)


!     1.1 GET ATMOSPHERIC MODEL FORCINGS FIELDS 
!         -------------------------------------

      ILEN=1
      LWNDFILE=.TRUE.
      LCLOSEWND=.TRUE.

      CALL GETWND (BLK2LOC,                             &
     &             NXFFS, NXFFE, NYFFS, NYFFE, FIELDG,  &
     &             WVENVI,                              &
     &             FF_NOW,                              &
     &             CDTPRO, LWNDFILE, LCLOSEWND, IREAD,  &
     &             LCR, NEMO2WAM,                       &
     &             ICODE_WND)

!     1.2 USE DATA FROM A FILE CONTAINING WIND SPEED MODIFIED BY
!         ----------------------------------------------------
!     A PREVIOUS WAVE MODEL RUN ON WHICH THIS RESTART IS BASED
!     IF IT EXISTS, OTHERWISE IT WILL USE THE INFORMATION FROM
!     THE ATMOSPHERIC MODEL WIND SPEED ARCHIVE.


      FILNM='uwavein'
      LIU = LEN_TRIM(FILNM)
      FILNM=FILNM(1:LIU)
      IF (IREAD == IRANK) THEN
        INQUIRE(FILE=FILNM, EXIST=LWAVEWIND)
        IF (LWAVEWIND) THEN
          NWAVEWIND(1)=1
        ELSE
          NWAVEWIND(1)=0
        ENDIF
      ENDIF

!     USE MESSAGE PASSING TO SEND FILE STATUS TO THE OTHER PE'S
      CALL GSTATS(696,0)
      IF (NPROC > 1) THEN
        CALL MPL_BROADCAST(NWAVEWIND,KROOT=IREAD,KTAG=KTAG, CDSTRING='BUILDSTRESS :') 
        KTAG=KTAG+1
        IF (NWAVEWIND(1) == 1) THEN
          LWAVEWIND=.TRUE.
        ELSE
          LWAVEWIND=.FALSE.
        ENDIF
      ENDIF
      CALL GSTATS(696,1)

      IF (LWAVEWIND) THEN
        IPARAM=245
        LLONLYPOS=.TRUE.
        CALL READWGRIB(IU06, FILNM, IPARAM, CDTPRO,            &
     &                 BLK2LOC,                                &
     &                 NXFFS, NXFFE, NYFFS, NYFFE, FIELDG,     &
     &                 FF_NOW%WSWAVE , KZLEVUWAVE, LLONLYPOS, IREAD)

        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' A DATA FILE CONTAINING WIND SPEED INFORMATION'
        WRITE(IU06,*) ' AS PROVIDED BY A PREVIOUS WAVE MODEL RUN'
        WRITE(IU06,*) ' WAS USED TO UPDATE THE INPUT ATMOSPHERIC WINDS'
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' THE INPUT WINDS AND DRAG COEFFICIENT ARE FOUND'
        WRITE(IU06,*) ' TO HAVE BEEN DETERMINED FOR HEIGHT AT ',        &
     &                  KZLEVUWAVE,' m'


      ELSE

        KZLEVUWAVE=10
        WRITE(IU06,*) ' '
        WRITE(IU06,*) '          !!!! NOTE !!!!'
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' NO INFORMATION ON WIND SPEEDS FROM THE'
        WRITE(IU06,*) ' WAVE MODEL WAS PROVIDED'
        WRITE(IU06,*) ' NO UPDATE OF THE INPUT ATMOSPHERIC WINDS'
        WRITE(IU06,*) ' WAS POSSIBLE'
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' THE INPUT WIND AND DRAG COEFFICIENT ARE ASSUMED'
        WRITE(IU06,*) ' TO HAVE BEEN DETERMINED FOR HEIGHT AT ',        &
     &                  KZLEVUWAVE,' m'

      ENDIF
 
!     TEST WHETHER THE HEIGHT OF THE INPUT WINDS IS THE SAME AS DEFINED
!     FOR THE REST OF THE RUN (SEE CALL TO STRESS IN PREPROC). IF NOT
!     RECOMPUTE THE TABLE. THE ORIGINAL TABLE WILL BE SWAP BACK AT THE
!     END OF THE ROUTINE.

      IF (KZLEVUWAVE /= NINT(XNLEV)) THEN
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' THE REFERENCE HEIGHT TO BE USED IN WAMODEL'
        WRITE(IU06,*) ' ',XNLEV 
        WRITE(IU06,*) ' IS DIFFERENT THAN THE INPUT FIELDS HEIGHT'
        WRITE(IU06,*) ' ',REAL(KZLEVUWAVE,JWRB)
        WRITE(IU06,*) ' THE NECESSARY ADJUSTMENTS WILL BE MADE'
        WRITE(IU06,*) ' TO DETERMINE THE INITIAL FIELDS.'
        WRITE(IU06,*) ' '

        TEMPXNLEV=KZLEVUWAVE
      ELSE
        TEMPXNLEV=XNLEV
      ENDIF

!     1.3 INITIALISE CD USING THE FRICTION VELOCITY FOR TAUW=0.
!         ----------------------------------------------------

      CALL GSTATS(1444,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK)
      DO ICHNK = 1, NCHNK
        CALL CDUSTARZ0 (1, NPROMA_WAM, FF_NOW%WSWAVE(:,ICHNK), TEMPXNLEV, CD(:,ICHNK), FF_NOW%UFRIC(:,ICHNK), FF_NOW%Z0M(:,ICHNK))
        FF_NOW%TAUW(:,ICHNK) = 0.1_JWRB * FF_NOW%UFRIC(:,ICHNK)**2
        FF_NOW%TAUWDIR(:,ICHNK) = FF_NOW%WDWAVE(:,ICHNK)
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1444,1)

!     1.4  GET DRAG COEFFICIENT
!          --------------------
      IF (.NOT.LNOCDIN) THEN
        IPARAM=233
        LLONLYPOS=.TRUE.
        FILNM='cdwavein'
!       !!!! CD was initialised above !!!!
        CALL READWGRIB(IU06, FILNM, IPARAM, CDTPRO,         &
     &                 BLK2LOC,                             &
     &                 NXFFS, NXFFE, NYFFS, NYFFE, FIELDG,  &
     &                 CD, KZLEVCD, LLONLYPOS, IREAD)

!       TEST REFERENCE LEVEL FOR UWAVE AND CD

        IF (KZLEVUWAVE /= 0 .AND. KZLEVCD /= 0) THEN
          IF (KZLEVUWAVE /= KZLEVCD) THEN
            WRITE(IU06,*)'************************************'
            WRITE(IU06,*)'*                                  *'
            WRITE(IU06,*)'* FATAL ERROR IN SUB BUILDSTRESS   *'
            WRITE(IU06,*)'* ==============================   *'
            WRITE(IU06,*)'* REFERENCE LEVEL FOR CD AND UWAVE *'
            WRITE(IU06,*)'* DIFFER. PROGRAM WILL ABORT       *'
            WRITE(IU06,*)'* REF. LEVEL CD = ',KZLEVCD
            WRITE(IU06,*)'* REF. LEVEL UWAVE = ',KZLEVUWAVE
            WRITE(IU06,*)'*                                  *'
            WRITE(IU06,*)'************************************'
            CALL ABORT1
          ENDIF
        ENDIF

!       1.5 COMPUTE TAUW, UFRIC AND Z0M
!           ----------------------------

        IF (LLCAPCHNK) THEN
          DO ICHNK = 1, NCHNK
            DO IJ = 1, NPROMA_WAM
              ALPHAOG(IJ, ICHNK) = CHNKMIN(FF_NOW%WSWAVE(IJ,ICHNK))*GM1
            ENDDO
          ENDDO
        ELSE
            ALPHAOG(:, :) = ALPHA*GM1
        ENDIF

        CALL GSTATS(1444,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) & 
!$OMP&  PRIVATE(ICHNK, KIJS, KIJL, IJ, USTAR, CDSQRTINV, Z0TOT, CHARNOCKOG)
        DO ICHNK = 1, NCHNK
          KIJS=1
          KIJL=NPROMA_WAM

          DO IJ = KIJS, KIJL
!!          UFRIC WILL FIRST CONTAIN ITS SQUARE
!!          THE NUMERICAL RELATION BETWEEN UFRIC AND WSWAVE SHOULD
!!          ALWAYS BE AS IN OUTBETA
            FF_NOW%UFRIC(IJ,ICHNK) = CD(IJ, ICHNK) * MAX(FF_NOW%WSWAVE(IJ,ICHNK)**2, EPSU10**2)
            FF_NOW%UFRIC(IJ,ICHNK) = MAX(FF_NOW%UFRIC(IJ,ICHNK), EPSUS)
            USTAR = SQRT(FF_NOW%UFRIC(IJ,ICHNK))
            CDSQRTINV = MIN(1._JWRB/SQRT(CD(IJ, ICHNK)), 50.0_JWRB)
            Z0TOT = TEMPXNLEV*EXP(-XKAPPA*CDSQRTINV)
!           Z0M SHOULD ONLY CONTAIN THE CHARNOCK CONTRIBUTION (see taut_z0) TO DETERMINE TAUW
            FF_NOW%Z0M(IJ,ICHNK) = MAX(Z0TOT - RNUM/USTAR, ALPHAOG(IJ, ICHNK)*FF_NOW%UFRIC(IJ,ICHNK))
            CHARNOCKOG = FF_NOW%Z0M(IJ,ICHNK)/FF_NOW%UFRIC(IJ,ICHNK)
            CHARNOCKOG = MAX(CHARNOCKOG, ALPHAOG(IJ, ICHNK))
            FF_NOW%TAUW(IJ,ICHNK) = MAX(FF_NOW%UFRIC(IJ,ICHNK)*(1._JWRB-(ALPHAOG(IJ, ICHNK)/CHARNOCKOG)**2), 0._JWRB)
            FF_NOW%TAUWDIR(IJ,ICHNK) = FF_NOW%WDWAVE(IJ,ICHNK)
            FF_NOW%CHRNCK(IJ,ICHNK) = G * CHARNOCKOG
            FF_NOW%UFRIC(IJ,ICHNK) = USTAR
          ENDDO

        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1444,1)
      ENDIF

      WRITE(IU06,*) ' SUB. BUILDSTRESS: INPUT OF RESTART FILES DONE'

      CALL FIELDG%DEALLOC()

IF (LHOOK) CALL DR_HOOK('BUILDSTRESS',1,ZHOOK_HANDLE)

END SUBROUTINE BUILDSTRESS
