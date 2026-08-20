! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE MPBCASTGRID(IU06, ISEND, ITAG)

! ----------------------------------------------------------------------
!**** *MPBCASTGRID* - BROADCAST THE CONTENT OF THE GRID FILE READ BY
!                     READPRE TO THE OTHER PE'S. IT WILL ALSO ALLOCATE
!                     THE NECESSARY ARRAYS ON THE RECEIVING PE'S AS THEY
!                     WERE NOT ON THOSE PE'S SINCE READPRE WAS NOT
!                     CALLED FOR THEM.

!     J. BIDLOT    ECMWF   OCTOBER 1997

!     PURPOSE.
!     --------
!     BROADCAST THE CONTENT OF THE GRID FILE READ BY READPRE
!     TO THE OTHER PE'S
!*    INTERFACE.
!     ----------

!     CALL *MPBCASTGRID*(IU06,ISEND,ITAG)

!     *IU06*      UNIT FOR PRINTER MESSAGES.
!     *ISEND*     RANK OF THE PROCESS ONTO WHICH FIELD IS COLLECTED
!     *ITAG*      TAG ASSOCIATED WITH AS A PARTICULAR CALL TO SUBROUTINE
!                 THIS IS NECESSARY TO DIFFERENTIATE THE DIFFERENT CALLS

!     METHOD.
!     -------
!     MPL_BROADCAST FROM PROCESSOR CORRESPONDING TO ISEND TO
!     ALL OTHER PROCESSORS.

!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_BROADCAST

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMAP   , ONLY : NGX      ,NGY      ,              &
     &            IPER     ,IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,    &
     &            AMONOP   ,XDELLA   ,XDELLO   ,NLONRGG  ,    &
     &            IQGAUSS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPRECR   ,NPRECI
      USE YOWSHAL  , ONLY : BATHY

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_BROADCAST

!----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU06, ISEND
      INTEGER(KIND=JWIM), INTENT(INOUT) :: ITAG
      INTEGER(KIND=JWIM), PARAMETER :: MFIRST=2
      INTEGER(KIND=JWIM) :: I, J, IJ, K, K1, K2, M, M1, M2, IC, L, KDEPTH, NGOU 
      INTEGER(KIND=JWIM) :: IKCOUNT, KCOUNT
      INTEGER(KIND=JWIM) :: MIC, MZC 
      INTEGER(KIND=JWIM),ALLOCATABLE :: ICOMBUF(:)

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),ALLOCATABLE :: ZCOMBUF(:)

!----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MPBCASTGRID',0,ZHOOK_HANDLE)

      IF (ISEND == 0 .OR. NPROC == 1) THEN
         WRITE (IU06,*) ''
!     1.1 SEND TO ALL PROCESSORS OTHER THAN ISEND
!         ------------------------------------------------
      ELSE

!       BUFFER SIZE MESSAGE AND THE FEW DIMENSIONS NEEDED TO 
!       ALLOCATE ALL ARRAYS
        ALLOCATE(ICOMBUF(MFIRST))

        IF (IRANK == ISEND) THEN
          IKCOUNT=0
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NGX
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NGY
          IF (IKCOUNT /= MFIRST) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* IKCOUNT /= MFIRST !!!*' 
            WRITE (IU06,*) '* ON IRANK = ',IRANK
            WRITE (IU06,*) '* IKCOUNT = ',IKCOUNT
            WRITE (IU06,*) '* MFIRST  = ',MFIRST
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF
        ENDIF

        CALL MPL_BROADCAST(ICOMBUF,KROOT=ISEND,KTAG=ITAG,CDSTRING='MPBCASTGRID:')
        ITAG=ITAG+1

        IF (IRANK /= ISEND) THEN
          IKCOUNT=0
          IKCOUNT=IKCOUNT+1
          NGX=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NGY=ICOMBUF(IKCOUNT)
          IF (IKCOUNT /= MFIRST) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
            WRITE (IU06,*) '* IKCOUNT .NE. MFIRST !!!*' 
            WRITE (IU06,*) '* ON IRANK = ',IRANK
            WRITE (IU06,*) '* IKCOUNT = ',IKCOUNT
            WRITE (IU06,*) '* MFIRST  - ',MFIRST
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF
        ENDIF
        DEALLOCATE(ICOMBUF)

        MIC=3+NGY
        MZC=6+NGX*NGY

!       ENCODE MAIN MESSAGE BUFFERS (ON PE=ISEND) AND
!       ALLOCATE ALL ARRAYS NEEDED TO KEEP THE BUFFERS ON THE OTHER PE'S

        ALLOCATE(ICOMBUF(MIC))
        ALLOCATE(ZCOMBUF(MZC))

        IF (IRANK /= ISEND) THEN

          IF (.NOT.ALLOCATED(NLONRGG)) ALLOCATE(NLONRGG(NGY))

          IF (ALLOCATED(BATHY)) DEALLOCATE(BATHY)
          ALLOCATE(BATHY(NGX,NGY))

        ELSE 
          KCOUNT=0
          IKCOUNT=0

          DO J=1,NGY
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=NLONRGG(J)
          ENDDO

          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=IPER
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=IRGG
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=IQGAUSS

          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=AMOWEP
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=AMOSOP
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=AMOEAP
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=AMONOP
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=XDELLA
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=XDELLO

          DO K=1,NGY
            DO I=1,NGX
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=BATHY(I,K)
            ENDDO
          ENDDO

          IF (IKCOUNT /= MIC) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
            WRITE (IU06,*) '* IKCOUNT NE MIC PRIOR   *'
            WRITE (IU06,*) '* CALL TO MPL_BROADCAST  *'
            WRITE (IU06,*) '* IKCOUNT =',IKCOUNT
            WRITE (IU06,*) '* MIC =',MIC
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF 
          IF (KCOUNT /= MZC) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
            WRITE (IU06,*) '* KCOUNT NE MZC PRIOR    *'
            WRITE (IU06,*) '* CALL TO MPL_BROADCAST  *'
            WRITE (IU06,*) '* KCOUNT =',KCOUNT
            WRITE (IU06,*) '* MZC =',MZC
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF 
        ENDIF

        CALL MPL_BROADCAST(ICOMBUF,KROOT=ISEND,KTAG=ITAG, CDSTRING='MPBCASTGRID 1:')
        ITAG=ITAG+1

        CALL MPL_BROADCAST(ZCOMBUF,KROOT=ISEND,KTAG=ITAG, CDSTRING='MPBCASTGRID 2:')
        ITAG=ITAG+1

        IF (IRANK /= ISEND) THEN
          KCOUNT=0
          IKCOUNT=0

          DO J=1,NGY
            IKCOUNT=IKCOUNT+1
            NLONRGG(J)=ICOMBUF(IKCOUNT)
          ENDDO

          IKCOUNT=IKCOUNT+1
          IPER=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          IRGG=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          IQGAUSS=ICOMBUF(IKCOUNT)

          KCOUNT=KCOUNT+1
          AMOWEP=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          AMOSOP=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          AMOEAP=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          AMONOP=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          XDELLA=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          XDELLO=ZCOMBUF(KCOUNT)

          DO K=1,NGY
            DO I=1,NGX
              KCOUNT=KCOUNT+1
              BATHY(I,K)=ZCOMBUF(KCOUNT)
            ENDDO
          ENDDO

          IF (IKCOUNT /= MIC) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
            WRITE (IU06,*) '* IKCOUNT NE MIC AFTER   *'
            WRITE (IU06,*) '* CALL TO MPL_BROADCAST  *'
            WRITE (IU06,*) '* IKCOUNT =',IKCOUNT
            WRITE (IU06,*) '* MIC =',MIC
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF 
          IF (KCOUNT /= MZC) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
            WRITE (IU06,*) '* KCOUNT NE MZC AFTER    *'
            WRITE (IU06,*) '* CALL TO MPL_BROADCAST  *'
            WRITE (IU06,*) '* KCOUNT =',KCOUNT
            WRITE (IU06,*) '* MZC =',MZC
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF 
        ENDIF

        DEALLOCATE(ICOMBUF,ZCOMBUF)

      ENDIF

      IF (LHOOK) CALL DR_HOOK('MPBCASTGRID',1,ZHOOK_HANDLE)

END SUBROUTINE MPBCASTGRID
