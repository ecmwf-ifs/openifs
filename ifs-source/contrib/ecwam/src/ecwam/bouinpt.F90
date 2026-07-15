! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE BOUINPT (IU02, FL1, NSTART, NEND)

! ----------------------------------------------------------------------

!**** *BOUINPT* - BOUNDARY VALUE INPUT INTO A BLOCK.

!     H. GUNTHER    GKSS/ECMWF         JANUARY 1991

!*    PURPOSE.
!     --------

!       INPUT AND SPACE INTERPOLATION OF BOUNDARY SPECTRA.

!**   INTERFACE.
!     ----------

!       *CALL* *BOUINPT (IU02)*
!          *IU02*    INTEGER      UNIT FOR INPUT OF BOUNDARY VALUES.
!          *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!          *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN

!     EXTERNALS.
!     ----------

!       *ABORT*     - TERMINATES PROCESSING.
!       *GSFILE*    - GETS OR SAVES A FILE.
!       *INTSPEC*   - INTERPOLATE A SPECTRUM.

!     METHOD.
!     -------

!       IN THE FIRST CALL OF THE SUB. THE FILE HEADER IS READ
!       AND THE CONSISTENCY IS CHECKED. THE SUB. READS A COMPLETE
!       SET OF BOUNDARY VALUES EACH PROPAGATION TIMESTEP WHEN
!       IT IS CALLED FOR THE FIRST BLOCK. THE SPECTRA REQUIRED FOR
!       THE ACTUAL BLOCK ARE INTERPOLATED AND STORED IN THE BLOCK.
!       INDICES AND WEIGHTS NECESSARY FOR THE INTERPOLATION AND
!       STORAGE ARE PRECOMPUTED IN PROG. PREPROC.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFPBO  , ONLY : NBOUNF   ,IJARF    ,IBFL     ,              &
     &                      IBFR     ,BFW
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, ICHNKFROMIJ, IPRMFROMIJ
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : CDATEF   ,CDTPRO   ,CDTBC   ,IDELBC    ,    &
     &                      IDELPRO
      USE YOWTEST  , ONLY : IU06
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_BROADCAST, MPL_SEND, MPL_RECV, MPL_WAIT, &
                           & MPL_ABORT, JP_BLOCKING_STANDARD, &
                           & JP_NON_BLOCKING_STANDARD

! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "abort1.intfb.h"
#include "gsfile_new.intfb.h"
#include "intspec.intfb.h"

      INTEGER(KIND=JWIM), INTENT(INOUT) :: IU02
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1
      INTEGER(KIND=JWIM), DIMENSION(NPROC), INTENT(IN) :: NSTART, NEND


      INTEGER(KIND=JWIM) :: I, M, K, IP, IK, IJ, IJF, IC, KL, ML, IDELINP
      INTEGER(KIND=JWIM), SAVE :: NBOINP
      INTEGER(KIND=JWIM) :: ISEND, IREQ
      INTEGER(KIND=JWIM) :: ICOUNT, NZCOMBUF, NINTFLD
      INTEGER(KIND=JWIM) :: KRCOUNT, KRTAG
      INTEGER(KIND=JWIM) :: ILOC, IBCL, IBCR
      INTEGER(KIND=JWIM) :: ISENDREQ(NPROC)
      INTEGER(KIND=JWIM) :: NIJB
      INTEGER(KIND=JWIM), DIMENSION(NBOUNF) :: IJB, IBND

      REAL(KIND=JWRB) :: XANG, XFRE, TH0, FR1, CO, XNBO, XDELIN  
      REAL(KIND=JWRB) :: XLON, XLAT
      REAL(KIND=JWRB) :: DEL12, DEL1L
      REAL(KIND=JWRB) :: FMEAN, EMEAN, THQ
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!      *FMEAN1*    REAL      MEAN FREQUENCIES FROM COARSE GRID.
!      *EMEAN1*    REAL      TOTAL ENERGIES FROM COARSE GRID.
!      *THQ1*      REAL      MEAN DIRECTIONS FROM COARSE GRID (RAD).
!      *NBOINP*    REAL      NUMBER OF INPUT SPECTRA.
!      *FL*        REAL      INTERPOLATED SPECTRUM.
      REAL(KIND=JWRB), ALLOCATABLE :: FMEAN1(:), EMEAN1(:), THQ1(:) 
      REAL(KIND=JWRB), ALLOCATABLE :: ZCOMBUFS(:), ZCOMBUFR(:)
      REAL(KIND=JWRB), ALLOCATABLE :: FL(:,:)
      REAL(KIND=JWRB), ALLOCATABLE :: F1(:,:,:)

      CHARACTER(LEN=14) :: CDTLAST, CDATE1
      CHARACTER(LEN=120) :: FILENAME

      DATA CDTLAST /'00000000000000'/  
!      *CDTLAST*   CHAR*14   DATE OF LAST FETCH FILE.

      LOGICAL :: LLBC(NPROC)

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('BOUINPT',0,ZHOOK_HANDLE)

!*    1. INITIALIZATION.
!        ---------------

!       FIND OUT WHICH PEs NEEDS BOUNDARY CONDITIONS
        DO IP=1,NPROC
          LLBC(IP)=.FALSE.
          DO I=1,NBOUNF
            IJF=IJARF(I)
            IF (IJF >= NSTART(IP) .AND. IJF <= NEND(IP)) THEN
              LLBC(IP)=.TRUE.
              EXIT
            ENDIF
          ENDDO
        ENDDO

!*      GET LOCAL INDICES OF BOUNDARY POINTS
        NIJB=0 
        DO I=1,NBOUNF
          IJF=IJARF(I)
          IF (IJF >= NSTART(IRANK) .AND. IJF <= NEND(IRANK)) THEN
            NIJB=NIJB+1 
            IBND(NIJB)=I
            IJB(NIJB)=IJF
          ENDIF
        ENDDO


!     1. 1 GET BOUNDARY DATA FROM PE ISEND
!          -------------------------------
      ISEND=1

!*    1.2 HAS A NEW BOUNDARY VALUE INPUT FILE TO BE FETCH?
!         ------------------------------------------------
      IF (CDTLAST < CDTBC) THEN

        NZCOMBUF=7
        ALLOCATE(ZCOMBUFS(NZCOMBUF))
        IF (IRANK == ISEND)THEN
!*        FETCH INPUT FILE.
          CALL GSFILE (IU06, IU02, 0, CDTBC, CDTBC, 'FBI', 'G')
!*        READ BOUNDARY FILE HEADER.
          READ (IU02, ERR=5000, END=5000) (ZCOMBUFS(IC), IC=1,NZCOMBUF)
        ENDIF

        CALL MPL_BROADCAST(ZCOMBUFS,KROOT=ISEND,KTAG=KTAG,CDSTRING='BOUINPT:')
        KTAG=KTAG+1

        XANG=ZCOMBUFS(1)
        XFRE=ZCOMBUFS(2)
        TH0=ZCOMBUFS(3)
        FR1=ZCOMBUFS(4)
        CO=ZCOMBUFS(5)
        XNBO=ZCOMBUFS(6)
        XDELIN=ZCOMBUFS(7)
        DEALLOCATE(ZCOMBUFS)

        KL = NINT(XANG)
        ML = NINT(XFRE)
        NBOINP  = NINT(XNBO)
        IDELINP = NINT(XDELIN)

        CDTLAST = CDTBC

!*      CHECK CONSISTENCY.

        IF (KL /= NANG .OR. ML /= NFRE .OR. FR1 /= FR(1) .OR.                               &
     &      TH0 /= TH(1) .OR. MOD(IDELINP,IDELPRO) /= 0 .OR. IDELINP < IDELPRO) THEN
          WRITE (IU06,*) '****************************************'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '*    FATAL ERROR SUB. BOUINP           *'
          WRITE (IU06,*) '*    =======================           *'
          WRITE (IU06,*) '* VALUES IN BOUNDARY FILE HEADER ARE   *'
          WRITE (IU06,*) '* INCONSISTENT WITH MODEL SET-UP.      *'
          WRITE (IU06,*) '* FINE MODEL VALUES ARE:               *'
          WRITE (IU06,*) '* NO. OF DIRECTIONS     NANG   = ', NANG
          WRITE (IU06,*) '* NO. OF FREQUENCIES    NFRE   = ', NFRE
          WRITE (IU06,*) '* FIRST DIRECTION       TH(1)  = ', TH(1)
          WRITE (IU06,*) '* FIRST FREQUENCY       FR(1)  = ', FR(1)
          WRITE (IU06,*) '* PROPAGATION TIMESTEP  IDELPRO= ', IDELPRO
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '* COARSE MODEL VALUES ARE:             *'
          WRITE (IU06,*) '* NO. OF DIRECTIONS     KL     = ', KL 
          WRITE (IU06,*) '* NO. OF FREQUENCIES    ML     = ', ML 
          WRITE (IU06,*) '* FIRST DIRECTION       TH0    = ', TH0
          WRITE (IU06,*) '* FIRST FREQUENCY       FR1    = ', FR1
          WRITE (IU06,*) '* PROPAGATION TIMESTEP  IDELINP= ', IDELINP
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.    *'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '****************************************'
          CALL ABORT1
        ENDIF

      ENDIF

      NINTFLD=0
      IF (LLBC(IRANK) .OR. IRANK == ISEND) THEN
        ALLOCATE(F1(NANG,NFRE,0:NBOINP))
        ALLOCATE(FL(NANG,NFRE))
        NINTFLD=NINTFLD+1
        ALLOCATE(FMEAN1(0:NBOINP))
        NINTFLD=NINTFLD+1
        ALLOCATE(EMEAN1(0:NBOINP))
        NINTFLD=NINTFLD+1
        ALLOCATE(THQ1(0:NBOINP))

!*     INITIALISE FOR LANDPOINTS.
        DO M=1,NFRE
          DO K=1,NANG
            F1(K,M,0) = 0.0_JWRB
          ENDDO
        ENDDO
        FMEAN1(0) = 0.0_JWRB
        EMEAN1(0) = 0.0_JWRB
        THQ1(0) = 0.0_JWRB
      ENDIF

! ----------------------------------------------------------------------

!*    2.1 READ BOUNDARY SPECTRA.
!         ----------------------

      IF (LLBC(IRANK) .OR. IRANK == ISEND) THEN
        NZCOMBUF = NBOINP*(NINTFLD+NANG*NFRE)
        ALLOCATE(ZCOMBUFS(NZCOMBUF))
      ENDIF

      IF (IRANK == ISEND)THEN
 2100   CONTINUE
        DO IJ=1,NBOINP
          READ (IU02, ERR=5001, END=5001) XLON, XLAT, CDATE1, EMEAN1(IJ), THQ1(IJ), FMEAN1(IJ)
          READ (IU02, ERR=5002, END=5002) ((F1(K,M,IJ),K=1,NANG),M=1,NFRE)
        ENDDO

!*    2.2 CHECK DATES.
!         ------------
        IF (CDATE1 < CDTPRO) THEN
          WRITE (IU06,*) '****************************************'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '*    WARNING ERROR SUB. BOUINP         *'
          WRITE (IU06,*) '*    =========================         *'
          WRITE (IU06,*) '* BOUNDARY VALUE INPUT DATE IS BEFORE  *'
          WRITE (IU06,*) '* MODEL DATE.                          *'
          WRITE (IU06,*) '* MODEL DATE IS          CDTPRO = ', CDTPRO
          WRITE (IU06,*) '* BOUNDARY VALUE DATE IS CDATE1 = ', CDATE1
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '* PROGRAM WILL TRY NEXT INPUT.         *'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '****************************************'
          GOTO 2100
        ELSEIF (CDATE1 > CDTPRO) THEN
          WRITE (IU06,*) '****************************************'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '*    FATAL ERROR SUB. BOUINP           *'
          WRITE (IU06,*) '*    =======================           *'
          WRITE (IU06,*) '* DATES DO NOT MATCH.                  *'
          WRITE (IU06,*) '* MODEL DATE IS          CDTPRO = ', CDTPRO
          WRITE (IU06,*) '* BOUNDARY VALUE DATE IS CDATE1 = ', CDATE1
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.    *'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '****************************************'
          CALL ABORT1
        ENDIF
      ENDIF

!     SEND INFORMATION TO OTHER PE'S
      IREQ=0
      IF (IRANK == ISEND)THEN
        ICOUNT=0
        DO IJ=1,NBOINP
          ICOUNT=ICOUNT+1
          ZCOMBUFS(ICOUNT)=EMEAN1(IJ)
          ICOUNT=ICOUNT+1
          ZCOMBUFS(ICOUNT)=THQ1(IJ)
          ICOUNT=ICOUNT+1
          ZCOMBUFS(ICOUNT)=FMEAN1(IJ)
          DO M=1,NFRE
            DO K=1,NANG
            ICOUNT=ICOUNT+1
            ZCOMBUFS(ICOUNT)=F1(K,M,IJ)
            ENDDO
          ENDDO
        ENDDO
        IF (ICOUNT /= NZCOMBUF) THEN
          WRITE (IU06,*) '***************************************'
          WRITE (IU06,*) '*    ERROR SUB. BOUINP                *'
          WRITE (IU06,*) '*    =================                *'
          WRITE (IU06,*) '* ICOUNT SHOULD BE EQUAL TO NZCOMBUF! *'
          WRITE (IU06,*) '* BEFORE SENDING THE DATA             *'
          WRITE (IU06,*) '* ICOUNT = ',ICOUNT
          WRITE (IU06,*) '* NZCOMBUF = ',NZCOMBUF
          WRITE (IU06,*) '***************************************'
          CALL ABORT1
        ENDIF

!       SEND BUFFER
        DO IP=1,NPROC
          IF (LLBC(IP)) THEN
            IREQ=IREQ+1
            CALL GSTATS(673,0)
            CALL MPL_SEND(ZCOMBUFS(1:NZCOMBUF),KDEST=IP,KTAG=KTAG,      &
     &                    KMP_TYPE=JP_NON_BLOCKING_STANDARD,            &
     &                    KREQUEST=ISENDREQ(IREQ),                      &
     &                    CDSTRING='BOUINPT:')
            CALL GSTATS(673,1)
          ENDIF
        ENDDO

      ENDIF
!
      IF (LLBC(IRANK)) THEN
        CALL GSTATS(673,0)
        ALLOCATE(ZCOMBUFR(NZCOMBUF))
        CALL MPL_RECV(ZCOMBUFR(1:NZCOMBUF),                             &
     &                KSOURCE=ISEND, KTAG=KTAG,                         &
     &                KOUNT=KRCOUNT, KRECVTAG=KRTAG,                    &
     &                KMP_TYPE=JP_BLOCKING_STANDARD,                    &
     &                CDSTRING='BOUINPT:')
        IF (KRCOUNT /= NZCOMBUF) CALL MPL_ABORT                          &
     &    ('MPL_RECV ERROR in BOUINPT: WRONG MESSAGE LENGTH')

        CALL GSTATS(673,1)

        ICOUNT=0
        DO IJ=1,NBOINP
          ICOUNT=ICOUNT+1
          EMEAN1(IJ)=ZCOMBUFR(ICOUNT)
          ICOUNT=ICOUNT+1
          THQ1(IJ)=ZCOMBUFR(ICOUNT)
          ICOUNT=ICOUNT+1
          FMEAN1(IJ)=ZCOMBUFR(ICOUNT)
          DO M=1,NFRE
            DO K=1,NANG
            ICOUNT=ICOUNT+1
            F1(K,M,IJ)=ZCOMBUFR(ICOUNT)
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(ZCOMBUFR)
        IF (ICOUNT /= NZCOMBUF) THEN
          WRITE (IU06,*) '***************************************'
          WRITE (IU06,*) '*    ERROR SUB. BOUINP                *'
          WRITE (IU06,*) '*    =================                *'
          WRITE (IU06,*) '* ICOUNT SHOULD BE EQUAL TO NZCOMBUF! *'
          WRITE (IU06,*) '* AFTER RECEIVING THE DATA            *'
          WRITE (IU06,*) '* ICOUNT = ',ICOUNT
          WRITE (IU06,*) '* NZCOMBUF = ',NZCOMBUF
          WRITE (IU06,*) '***************************************'
          CALL ABORT1
        ENDIF
      ENDIF

      IF( IREQ > 0 )THEN
        CALL GSTATS(673,0)
        CALL MPL_WAIT(KREQUEST=ISENDREQ(1:IREQ),                        &
     &                CDSTRING='BOUINPT: WAIT MPL_SEND')
        CALL GSTATS(673,1)
      ENDIF

      KTAG=KTAG+1

! ----------------------------------------------------------------------

!*    3. LOOP OVER BOUNDARY POINTS.
!        --------------------------

      DEL12 = 1.0_JWRB

      DO ILOC=1,NIJB

        IJF=IJB(ILOC)
        IK = ICHNKFROMIJ(IJF)
        IP = IPRMFROMIJ(IJF)

        I=IBND(ILOC)
        IBCL = IBFL(I)
        IBCR = IBFR(I)
        DEL1L = BFW(I)


!*    3.2 CHECK INTERPOLATION WEIGHT.
!         ---------------------------

        IF (DEL1L > 0.0_JWRB ) THEN

!*    3.2.1. WEIGHT IS GT ZERO, INTERPOLATE.
!            -------------------------------
          CALL INTSPEC (NFRE, NANG, NFRE, NANG, FR, DEL12, DEL1L,       &
     &       F1(1,1,IBCL), FMEAN1(IBCL), EMEAN1(IBCL), THQ1(IBCL),      &
     &       F1(1,1,IBCR), FMEAN1(IBCR), EMEAN1(IBCR), THQ1(IBCR),      &
     &       FL1(IP, :, :, IK), FMEAN, EMEAN, THQ)
        ELSE

!*    3.2.2. WEIGHT IS ZERO COPY LEFT POINT.
!            -------------------------------
          DO M=1,NFRE
            DO K=1,NANG
              FL1(IP, K, M, IK) = F1(K,M,IBCL)
            ENDDO
          ENDDO

        ENDIF
      ENDDO

! ----------------------------------------------------------------------

!*    4. RETURN TO CALLING PROG.
!        -----------------------

      IF (ALLOCATED(ZCOMBUFS)) DEALLOCATE(ZCOMBUFS)
      IF (ALLOCATED(F1)) DEALLOCATE(F1)
      IF (ALLOCATED(FL)) DEALLOCATE(FL)
      IF (ALLOCATED(FMEAN1)) DEALLOCATE(FMEAN1)
      IF (ALLOCATED(EMEAN1)) DEALLOCATE(EMEAN1)
      IF (ALLOCATED(THQ1)) DEALLOCATE(THQ1)

      RETURN

! ----------------------------------------------------------------------

!*    5. ERROR MESSAGES.
!        ---------------

 5000 CONTINUE
      WRITE (IU06,*) '*******************************************'
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '*      FATAL ERROR SUB. BOUNINPT.         *'
      WRITE (IU06,*) '*      ==========================         *'
      WRITE (IU06,*) '* END OF FILE OR READ ERROR.              *'
      WRITE (IU06,*) '* PROGRAM TRIES TO READ                   *'
      WRITE (IU06,*) '* HEADER OF BOUNDARY VALUES               *'
      WRITE (IU06,*) '* UNIT IS IU02 = ', IU02
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '*******************************************'
      CALL ABORT1
 5001 CONTINUE
      WRITE (IU06,*) '*******************************************'
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '*      FATAL ERROR SUB. BOUNINPT.         *'
      WRITE (IU06,*) '*      ==========================         *'
      WRITE (IU06,*) '* END OF FILE OR READ ERROR.              *'
      WRITE (IU06,*) '* PROGRAM TRIES TO READ                   *'
      WRITE (IU06,*) '* A SPECTRUM HEADER                       *'
      WRITE (IU06,*) '* SPECTRA COUNTER IS IJ = ', IJ
      WRITE (IU06,*) '* UNIT IS          IU02 = ', IU02
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '*******************************************'
      CALL ABORT1
 5002 CONTINUE
      WRITE (IU06,*) '*******************************************'
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '*      FATAL ERROR SUB. BOUNINPT.         *'
      WRITE (IU06,*) '*      ==========================         *'
      WRITE (IU06,*) '* END OF FILE OR READ ERROR.              *'
      WRITE (IU06,*) '* PROGRAM TRIES TO READ A SPECTRUM        *'
      WRITE (IU06,*) '* DATE IS        CDATE1 = ', CDATE1
      WRITE (IU06,*) '* SPECTRA COUNTER IS IJ = ', IJ
      WRITE (IU06,*) '* UNIT IS          IU02 = ', IU02
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '*******************************************'
      CALL ABORT1

      IF (LHOOK) CALL DR_HOOK('BOUINPT',1,ZHOOK_HANDLE)

END SUBROUTINE BOUINPT
