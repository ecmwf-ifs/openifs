! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE MPGATHERBC(IRECV, NSCFLD,                   &
 &                    FL1, EMEAN, THQ, FMEAN,          &
 &                    FLPTS, EMPTS, TQPTS, FMPTS)

!***  *MPGATHERBC* - MESSAGE PASSING OF COARSE BOUNDARY CONDITIONS.


!     PURPOSE                                                       
!     -------                                                      

!     *MPGATHERBC* - GATHERS TO PE IRECV ALL COARSE BOUNDARY CONDITIONS

!      *IRECV*   - PE RECEIVING ALL CONTRIBUTIONS.
!      *NSCFLD*  - NUMBER OF SCALAR FIELDS TO BE GATHERED
!      *FL1*     - BLOCK OF SPECTRA.
!      *EMEAN*   - BLOCK TOTAL ENERGY
!      *THQ*     - BLOCK MEAN DIRECTION
!      *FMEAN*   - BLOCK MEAn FREQUENCY
!       THE FOLLOWING ARRAY WILL BE GATHERED TO COARSE BOUNDARY POINTS.
!                  FL1   --> FLPTS
!                  EMEAN --> EMPTS
!                  THQ   --> TQPTS
!                  FMEAN --> FMPTS

!----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCPBO  , ONLY : NBOUNC   ,IJARC
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, ICHNKFROMIJ, IPRMFROMIJ 
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,LL1D
      USE YOWSPEC  , ONLY : IJ2NEWIJ ,NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06
      USE MPL_MODULE, ONLY : MPL_SEND, MPL_RECV, MPL_WAIT, MPL_BARRIER, &
                           & MPL_ABORT, &
                           & JP_NON_BLOCKING_STANDARD, JP_BLOCKING_STANDARD

!----------------------------------------------------------------------
      IMPLICIT NONE

#include "abort1.intfb.h"

      INTEGER(KIND=JWIM) :: IRECV, NSCFLD
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN) :: EMEAN, FMEAN, THQ
      REAL(KIND=JWRB), DIMENSION(NBOUNC), INTENT(INOUT) :: EMPTS, TQPTS, FMPTS
      REAL(KIND=JWRB), DIMENSION(NBOUNC, NANG, NFRE), INTENT(INOUT)   :: FLPTS

      INTEGER(KIND=JWIM) :: ITAG
      INTEGER(KIND=JWIM) :: MAXLENGTH, MESLENGTH
      INTEGER(KIND=JWIM) :: NGOU, IJ, M, K, IP, KCOUNT, IST, IND, IGOU, ICHNK, IK
      INTEGER(KIND=JWIM) :: ISENDREQ
      INTEGER(KIND=JWIM) :: KRFROM, KRCOUNT, KRTAG
      INTEGER(KIND=JWIM), DIMENSION(NPROC) :: NPTS, IBS, IBL
      INTEGER(KIND=JWIM), DIMENSION(NBOUNC) :: IJBC, NGOUG

      REAL(KIND=JWRB), ALLOCATABLE :: ZCOMBUFS(:), ZCOMBUFR(:)

!----------------------------------------------------------------------

      ITAG=1

!     1.0 DEFAULT ACTION IF NO FIELD GATHERING 
!         ------------------------------------
      IF (IRECV == 0 .OR. NPROC == 1) THEN

        DO NGOU=1,NBOUNC
          IJ = IJARC(NGOU)
          IK = ICHNKFROMIJ(IJ)
          IP = IPRMFROMIJ(IJ)

          EMPTS(NGOU) = EMEAN(IP, IK)
          TQPTS(NGOU) = THQ(IP, IK)
          FMPTS(NGOU) = FMEAN(IP, IK)
          DO  M = 1, NFRE
            DO K = 1, NANG
              FLPTS(NGOU,K,M) = FL1(IP, K, M, IK)
            ENDDO
          ENDDO
        ENDDO

      ELSE

!       DETERMINE WHICH PE'S NEED TO COMMUNICATE
!   should be computed only once and moved to a module
        KCOUNT=0
        DO IP=1,NPROC
          NPTS(IP)=0
          IBS(IP)=0
          IBL(IP)=0
          DO NGOU=1,NBOUNC
            IF (LL1D) THEN
              IJ = IJARC(NGOU)
            ELSE
              IJ = IJ2NEWIJ(IJARC(NGOU))
            ENDIF
            IF (IJ >= NSTART(IP) .AND. IJ <= NEND(IP)) THEN
              NPTS(IP)=NPTS(IP)+1
              KCOUNT=KCOUNT+1
              IJBC(KCOUNT)=IJ
              NGOUG(KCOUNT)=NGOU
            ENDIF
            IF (NPTS(IP) == 1) IBS(IP)=KCOUNT
          ENDDO
          IF (NPTS(IP) > 0) THEN
            IBL(IP)=IBS(IP)+NPTS(IP)-1
          ENDIF
        ENDDO

        MESLENGTH=NSCFLD+NANG*NFRE
        MAXLENGTH=MESLENGTH*NBOUNC
        ALLOCATE(ZCOMBUFS(MAXLENGTH))

!     1.1 SEND TO THE PROCESS THAT GATHERS THE WHOLE FIELD
!         ------------------------------------------------

        IF (NPTS(IRANK) > 0) THEN
          KCOUNT=MESLENGTH*(IBS(IRANK)-1)
          IST=MESLENGTH*(IBS(IRANK)-1)+1
          IND=IST+NPTS(IRANK)*MESLENGTH-1

          DO IGOU=IBS(IRANK),IBL(IRANK)
            IJ = IJBC(IGOU)
            IK = ICHNKFROMIJ(IJ)
            IP = IPRMFROMIJ(IJ)

            KCOUNT=KCOUNT+1
            ZCOMBUFS(KCOUNT) = EMEAN(IP, IK)
            KCOUNT=KCOUNT+1
            ZCOMBUFS(KCOUNT) = THQ(IP, IK)
            KCOUNT=KCOUNT+1
            ZCOMBUFS(KCOUNT) = FMEAN(IP, IK)
            DO  M = 1, NFRE
              DO K = 1, NANG
                KCOUNT=KCOUNT+1
                ZCOMBUFS(KCOUNT) = FL1(IP, K, M, IK)
              ENDDO
            ENDDO
          ENDDO
          IF (KCOUNT /= IND) THEN
            WRITE (IU06,*) ' ******************************************'
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' *      FATAL ERROR SUB. MPGATHERBC       *'
            WRITE (IU06,*) ' *      ==========================        *'
            WRITE (IU06,*) ' *      KCOUNT SHOULD BE EQUAL TO IND     *'
            WRITE (IU06,*) ' *      KCOUNT = ',KCOUNT
            WRITE (IU06,*) ' *      IND = ',IND
            WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS       *'
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' ******************************************'
            CALL ABORT1
          ENDIF

!         SEND BUFFER
          CALL GSTATS(673,0)
          CALL MPL_SEND(ZCOMBUFS(IST:IND),KDEST=IRECV,KTAG=ITAG,        &
     &                  KMP_TYPE=JP_NON_BLOCKING_STANDARD,              &
     &                  KREQUEST=ISENDREQ,                              &
     &                  CDSTRING='MPGATHERBC:')
          CALL GSTATS(673,1)
        ENDIF

!       RECEIVE CONTRIBUTIONS
!       ---------------------
        IF (IRANK == IRECV) THEN
          ALLOCATE(ZCOMBUFR(MAXLENGTH))
          CALL GSTATS(673,0)
          DO IP=1,NPROC
            IF (NPTS(IP) > 0) THEN
              IST=MESLENGTH*(IBS(IP)-1)+1
              IND=IST+NPTS(IP)*MESLENGTH-1
              KCOUNT=IND-IST+1
              CALL MPL_RECV(ZCOMBUFR(IST:IND),                          &
     &                      KSOURCE=IP, KTAG=ITAG,                      &
     &                      KOUNT=KRCOUNT, KRECVTAG=KRTAG,              &
     &                      KMP_TYPE=JP_BLOCKING_STANDARD,              &
     &                      CDSTRING='MPGATHERBC:')
              IF (KRTAG /= ITAG) CALL MPL_ABORT                         &
     &            ('MPL_RECV ERROR in MPGATHERBC:  MISMATCHED TAGS' )
              IF (KRCOUNT /= KCOUNT) CALL MPL_ABORT                     &
     &            ('MPL_RECV ERROR in MPGATHERBC: WRONG MESSAGE LENGTH')
            ENDIF
          ENDDO
          CALL GSTATS(673,1)

!         DECODE BUFFER
          KCOUNT=0
          DO IGOU=1,NBOUNC
            NGOU=NGOUG(IGOU)
            KCOUNT=KCOUNT+1
            EMPTS(NGOU)=ZCOMBUFR(KCOUNT)
            KCOUNT=KCOUNT+1
            TQPTS(NGOU)=ZCOMBUFR(KCOUNT)
            KCOUNT=KCOUNT+1
            FMPTS(NGOU)=ZCOMBUFR(KCOUNT)
            DO  M = 1, NFRE
              DO K = 1, NANG
                KCOUNT=KCOUNT+1
                FLPTS(NGOU,K,M) = ZCOMBUFR(KCOUNT)
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(ZCOMBUFR)
        ENDIF

        IF (NPTS(IRANK) > 0) THEN
          CALL GSTATS(673,0)
          CALL MPL_WAIT(KREQUEST=ISENDREQ, CDSTRING='MPGATHERBC:')
          CALL GSTATS(673,1)
        ENDIF

        DEALLOCATE(ZCOMBUFS)
      ENDIF

      CALL MPL_BARRIER(CDSTRING='MPGATHERBC:')

END SUBROUTINE MPGATHERBC
