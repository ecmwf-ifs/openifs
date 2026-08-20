! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MPDISTRIBSCFLD(ISEND, ITAG, NBLKS, NBLKE, FIELD)

! ----------------------------------------------------------------------

!****  *MPDISTRIBSCFLD* - DISTRIBUTE FIELD ACROSS PROCESSORS 

!     J. BIDLOT    ECMWF   SEPTEMBER 1997 

!     PURPOSE.
!     --------
!     SEND THE RESPECTIVE CONTRIBUTION OF ARRAY FIELD FROM PROCESSOR 
!     ISEND TO THE OTHER PE's.

!*    INTERFACE.
!     ----------

!     CALL *MPDISTRIBFL*(ISEND,ITAG,NBLKS,NBLKE,FIELD) 

!     *ISEND*     RANK OF THE PROCESS ONTO WHICH FIELD IS COLLECTED 
!     *ITAG*      TAG ASSOCIATED WITH AS A PARTICULAR CALL TO SUBROUTINE
!                 THIS IS NECESSARY TO DIFFERENTIATE THE DIFFERENT CALLS
!     *NBLKS*     INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NBLKE*     INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *FIELD*     INPUT/OUTPUT ARRAY

!     METHOD.
!     -------
!     MPL SEND OF ARRAY FIELD FROM PROCESSOR CORRESPONDING TO ISEND TO 
!     ALL OTHER PROCESSORS.

!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_SEND
!         MPL_RECV
!         MPL_ABORT
!         MPL_WAIT

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
      USE MPL_MODULE,ONLY : MPL_SEND, MPL_RECV, MPL_WAIT, MPL_ABORT, &
                          & JP_NON_BLOCKING_STANDARD
      USE YOWABORT  ,ONLY : WAM_ABORT

!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: ISEND, ITAG
      INTEGER(KIND=JWIM), DIMENSION(NPROC), INTENT(IN) :: NBLKS, NBLKE

      REAL(KIND=JWRB), DIMENSION(NIBLO), INTENT(INOUT) :: FIELD


      INTEGER(KIND=JWIM) :: IR, IP, IJ, MPLENGTH, KCOUNT, KRTAG, KRCOUNT
      INTEGER(KIND=JWIM), DIMENSION(NPROC) :: ISENDREQ

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), ALLOCATABLE :: ZCOMBUFR(:)
      REAL(KIND=JWRB), ALLOCATABLE :: ZCOMBUFS(:,:)

      REAL(KIND=JWRU), ALLOCATABLE :: AC(:)

!----------------------------------------------------------------------

      IF (ISEND.EQ.0 .OR. NPROC.EQ.1) RETURN 

      IF (LHOOK) CALL DR_HOOK('MPDISTIBSCFLD',0,ZHOOK_HANDLE)

      IF (IRANK.EQ.ISEND) THEN
!     1.1 SEND NON BLOCKING TO ALL PROCESSORS
!         -----------------------------------
        MPLENGTH=MPMAXLENGTH
        ALLOCATE(ZCOMBUFS(MPLENGTH,NPROC))

        DO IP=1,NPROC
          KCOUNT=0
          DO IJ=NBLKS(IP),NBLKE(IP)
            KCOUNT=KCOUNT+1
            ZCOMBUFS(KCOUNT,IP)=FIELD(IJ)
          ENDDO
        ENDDO

        DO IP=1,NPROC
          CALL GSTATS(624,0)
          CALL MPL_SEND(ZCOMBUFS(1:MPLENGTH,IP),KDEST=IP,KTAG=ITAG,     &
     &       KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(IP),   &
     &       CDSTRING='MPDISTRIBSCFLD SEND:')
          CALL GSTATS(624,1)
        ENDDO

      ENDIF

!     1.2 RECEIVE CONTRIBUTION TO THE FIELD FROM PE ISEND 
!         ------------------------------------------------ 

      MPLENGTH=MPMAXLENGTH
      ALLOCATE(ZCOMBUFR(MPLENGTH))

      CALL GSTATS(624,0)
      CALL MPL_RECV(ZCOMBUFR(1:MPLENGTH),KSOURCE=ISEND,KTAG=ITAG,       &
     &   KOUNT=KRCOUNT,KRECVTAG=KRTAG,CDSTRING='MPDISTRIBSCFLD:')
      IF (KRCOUNT.NE.MPLENGTH) CALL MPL_ABORT                           &
     &   ('MPL_RECV ERROR in MPDISTRIBSCFLD:MISMATCHED MSG LENGTH')
      IF (KRTAG.NE.ITAG) CALL MPL_ABORT                                 &
     &   ('MPL_RECV ERROR in MPDISTRIBSCFLD MISMATCHED TAGS' )
      CALL GSTATS(624,1)

      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
        ALLOCATE(AC(MNP))
        DO IP = 1, RANK(IRANK)%NP
          AC(IP) = REAL(ZCOMBUFR(IP), KIND=JWRU)  !DBLE(ZCOMBUFR(IP))
        ENDDO
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ELSE
        KCOUNT=0
        DO IJ=NBLKS(IRANK),NBLKE(IRANK)
          KCOUNT=KCOUNT+1
          FIELD(IJ)=ZCOMBUFR(KCOUNT)
        ENDDO
      ENDIF

!     1.3 WAIT ANY OUTSTANDING SENDS TO COMPLETE
!         --------------------------------------

      IF (IRANK.EQ.ISEND) THEN
        CALL MPL_WAIT(KREQUEST=ISENDREQ, CDSTRING='MPDISTRIBSCFLD:')
        DEALLOCATE(ZCOMBUFS)
      ENDIF

      DEALLOCATE(ZCOMBUFR)


!!! it's not very prety but when unstructured, you also need to have the values on the halo
!!! it should be combined with the previous exchange
      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
        CALL EXCHANGE(AC)
        DO IP = 1 , MNP
          FIELD(IP) = AC(IP)
        ENDDO
        DEALLOCATE(AC)
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ENDIF

      IF (LHOOK) CALL DR_HOOK('MPDISTIBSCFLD',1,ZHOOK_HANDLE)

      END SUBROUTINE MPDISTRIBSCFLD
