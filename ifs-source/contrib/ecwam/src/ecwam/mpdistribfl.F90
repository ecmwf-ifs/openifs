! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MPDISTRIBFL(ISEND, ITAG, NBLKS, NBLKE, KINF, KSUP,     &
     &                       MINF, MSUP, FL)

! ----------------------------------------------------------------------

!****  *MPDISTRIBFL* - DISTRIBUTE FL ACROSS PROCESSORS 

!     J. BIDLOT    ECMWF   SEPTEMBER 1997 

!     PURPOSE.
!     --------
!     SEND THEIR RESPECTIVE CONTRIBUTION OF ARRAY FL 
!     FROM PROCESSOR ISEND TO THE OTHER PE's.

!*    INTERFACE.
!     ----------

!      CALL *MPDISTRIBFL(ISEND, ITAG, NBLKS, NBLKE, KINF, KSUP,
!    &                   MINF, MSUP, FL)

!     *ISEND*     RANK OF THE PROCESS ONTO WHICH FIELD IS COLLECTED 
!     *ITAG*      TAG ASSOCIATED WITH AS A PARTICULAR CALL TO SUBROUTINE
!                 THIS IS NECESSARY TO DIFFERENTIATE THE DIFFERENT CALLS 
!     *NBLKS*     INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NBLKE*     INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *KINF*      INDEX OF THE FIRST DIRECTION OF FL TO BE DISTRIBUTED
!     *KSUP*      INDEX OF THE LAST DIRECTION OF FL TO BE DISTRIBUTED
!     *MINF*      INDEX OF THE FIRST FREQUENCY OF FL TO BE DISTRIBUTED
!     *MSUP*      INDEX OF THE LAST FREQUENCY OF FL TO BE DISTRIBUTED
!     *FL*        INPUT/OUTPUT ARRAY CONTAINING THE PART OF THE SPECTRUM 

!     METHOD.
!     -------
!     MPL SEND OF ARRAY FL FROM PROCESSOR CORRESPONDING TO ISEND TO 
!     ALL OTHER PROCESSORS.

!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_SEND
!         MPL_RECV
!         MPL_ABORT 

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMAP   , ONLY : NIBLO
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,MPMAXLENGTH
      USE YOWPARAM , ONLY : LLUNSTR
#ifdef WAM_HAVE_UNWAM
      USE YOWPD,     ONLY : RANK, MNP=>NPA, EXCHANGE
#endif
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_SEND, MPL_RECV, MPL_WAIT, MPL_ABORT, &
                           & JP_NON_BLOCKING_STANDARD
      USE YOWABORT, ONLY : WAM_ABORT

!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: ISEND, ITAG, KINF, KSUP, MINF, MSUP
      INTEGER(KIND=JWIM), DIMENSION(NPROC), INTENT(IN) :: NBLKS, NBLKE

      REAL(KIND=JWRB), DIMENSION(NIBLO,KINF:KSUP,MINF:MSUP), INTENT(INOUT) :: FL


      INTEGER(KIND=JWIM) :: IR, IP, M, K, IJ, MPLENGTH, KCOUNT, KRTAG,  &
     &                      KRCOUNT
      INTEGER(KIND=JWIM),DIMENSION(NPROC) :: ISENDREQ

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), ALLOCATABLE :: ZCOMBUFS(:,:), ZCOMBUFR(:)

      REAL(KIND=JWRU), ALLOCATABLE :: AC(:,:,:)

!----------------------------------------------------------------------

      IF (ISEND.EQ.0 .OR. NPROC.EQ.1) RETURN

      IF (LHOOK) CALL DR_HOOK('MPDISTRIBFL',0,ZHOOK_HANDLE)

      MPLENGTH=MPMAXLENGTH*(KSUP-KINF+1)*(MSUP-MINF+1)
      ALLOCATE(ZCOMBUFR(MPLENGTH))

      IF (IRANK.EQ.ISEND) THEN
!     1.1 SEND NON BLOCKING TO ALL PROCESSORS
!         -----------------------------------
        ALLOCATE(ZCOMBUFS(MPLENGTH,NPROC))

        DO IP=1,NPROC
          KCOUNT=0
          DO M=MINF,MSUP
            DO K=KINF,KSUP
              DO IJ=NBLKS(IP),NBLKE(IP)
                KCOUNT=KCOUNT+1
                ZCOMBUFS(KCOUNT,IP)=FL(IJ,K,M)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        DO IP=1,NPROC
          CALL MPL_SEND(ZCOMBUFS(1:MPLENGTH,IP),KDEST=IP,KTAG=ITAG,     &
     &      KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(IP),    &
     &      CDSTRING='MPDISTRIBFL SEND:')
        ENDDO
      ENDIF

!     1.2 RECEIVE CONTRIBUTION TO THE FIELD FROM PE ISEND 
!         ------------------------------------------------ 

      CALL MPL_RECV(ZCOMBUFR(1:MPLENGTH),KSOURCE=ISEND,KTAG=ITAG,       &
     &   KOUNT=KRCOUNT,KRECVTAG=KRTAG,CDSTRING='MPDISTRIBFL:')
      IF (KRCOUNT.NE.MPLENGTH) CALL MPL_ABORT                           &
     &   ('MPL_RECV ERROR in MPDISTRIBFL:MISMATCHED MESSAGE LENGTH')
      IF (KRTAG.NE.ITAG) CALL MPL_ABORT                                 &
     &   ('MPL_RECV ERROR in MPDISTRIBFL MISMATCHED TAGS' )

      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
        ALLOCATE(AC(MNP,KINF:KSUP,MINF:MSUP))
        DO M=MINF,MSUP
          DO K=KINF,KSUP
            DO IP = 1, RANK(IRANK)%NP
              AC(IP,K,M) = REAL(ZCOMBUFR(IP), KIND=JWRU)  !DBLE(ZCOMBUFR(IP))
            ENDDO
          ENDDO 
        ENDDO
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ELSE
        KCOUNT=0
        DO M=MINF,MSUP
          DO K=KINF,KSUP
            DO IJ=NBLKS(IRANK),NBLKE(IRANK)
              KCOUNT=KCOUNT+1
              FL(IJ,K,M)=ZCOMBUFR(KCOUNT)
            ENDDO
          ENDDO 
        ENDDO 
      ENDIF

!     1.3 WAIT ANY OUTSTANDING SENDS TO COMPLETE
!         --------------------------------------

      IF (IRANK.EQ.ISEND) THEN
        CALL MPL_WAIT(KREQUEST=ISENDREQ, CDSTRING='MPDISTRIBFL:')
        DEALLOCATE(ZCOMBUFS)
      ENDIF

      DEALLOCATE(ZCOMBUFR)

!!! it's not very prety but when unstructured, you also need to have the values on the halo
!!! it should be combined with the previous exchange
      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
        DO M=MINF,MSUP
          DO K=KINF,KSUP
            CALL EXCHANGE(AC(:,K,M))
            DO IP = 1 , MNP
              FL(IP,K,M) = AC(IP,K,M)
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(AC)
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ENDIF

      IF (LHOOK) CALL DR_HOOK('MPDISTRIBFL',1,ZHOOK_HANDLE)

      END SUBROUTINE MPDISTRIBFL
