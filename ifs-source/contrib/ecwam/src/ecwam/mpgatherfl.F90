! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MPGATHERFL(IRECV, NBLKS, NBLKE, KINF, KSUP, MINF, MSUP, FL)

!****  *MPGATHERFL* - GATHER FL ONTO A SINGLE PROCESS 

!     J. BIDLOT    ECMWF   APRIL 1996  MESSAGE PASSING
!     J. BIDLOT    ECMWF   MARCH 1997  add use of MINF AND MSUP
!     J. BIDLOT    ECMWF   JANUARY 2003 use of MPL_GATHERV
!     J. BIDLOT    ECMWF   JANUARY 2009 use of KINF AND KSUP

!     PURPOSE.
!     --------

!     GATHER ARRAY FL DISTRIBUTED ACROSS THE
!     DIFFERENT PROCESSES ONTO THE SINGLE PROCESS IRECV

!*    INTERFACE.
!     ----------

!     CALL *MPGATHERFL*(IRECV,NBLKS,NBLKE,KINF,KSUP,MINF,MSUP,FL)

!     *IRECV*     RANK OF THE PROCESS ONTO WHICH FIELD IS COLLECTED 
!     *NBLKS*     INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NBLKE*     INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *KINF*      INDEX OF THE FIRST DIRECTION OF FL
!     *KSUP*      INDEX OF THE LAST DIRECTION OF FL
!     *MINF*      INDEX OF THE FIRST FREQUENCY OF FL
!     *MSUP*      INDEX OF THE LAST FREQUENCY OF FL
!     *FL*        INPUT/OUTPUT ARRAY CONTAINING THE PART OF THE SPECTRUM 

!     METHOD.
!     -------
!     MPL SEND OF ARRAY FL TO PROCESS CORRESPONDING TO IRECV FOR
!     ALL PROCESS EXCEPT FOR THE PROCESS CORRESPONDING TO IRECV
!     WHERE IT IS RECEIVED.

!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_GATHERV

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMAP   , ONLY : NIBLO
      USE YOWPARAM , ONLY : LLUNSTR
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_GATHERV

!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IRECV, KINF, KSUP, MINF, MSUP
      INTEGER(KIND=JWIM), DIMENSION(NPROC), INTENT(IN) :: NBLKS, NBLKE

      REAL(KIND=JWRB), DIMENSION(NIBLO,KINF:KSUP,MINF:MSUP), INTENT(INOUT) :: FL


      INTEGER(KIND=JWIM) :: N23, ILEN, KCOUNT, M, K, IJ, ILENR, IP

      INTEGER(KIND=JWIM), DIMENSION(NPROC) :: KRECVCOUNTS

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), ALLOCATABLE,DIMENSION(:) :: PSENDBUF, PRECVBUF

!----------------------------------------------------------------------

      IF (IRECV == 0 .OR. NPROC == 1) RETURN

      IF (LHOOK) CALL DR_HOOK('MPGATHERFL',0,ZHOOK_HANDLE)

      N23=(KSUP-KINF+1)*(MSUP-MINF+1)

      ILEN=(NBLKE(IRANK)-NBLKS(IRANK)+1)*N23
      ALLOCATE(PSENDBUF(ILEN))
      KCOUNT=0
      DO M=MINF,MSUP
        DO K=KINF,KSUP
          DO IJ=NBLKS(IRANK),NBLKE(IRANK)
            KCOUNT=KCOUNT+1
            PSENDBUF(KCOUNT)=FL(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

      ILENR=0
      DO IP=1,NPROC
        KRECVCOUNTS(IP)=(NBLKE(IP)-NBLKS(IP)+1)*N23
        ILENR=ILENR+KRECVCOUNTS(IP)
      ENDDO
      ALLOCATE(PRECVBUF(ILENR))

      CALL MPL_GATHERV(PSENDBUF=PSENDBUF,                               &
     &                 KRECVCOUNTS=KRECVCOUNTS,                         &
     &                 PRECVBUF=PRECVBUF,                               &
     &                 KROOT=IRECV,                                     &
     &                 CDSTRING='MPGATHERFL:')

      IF (IRANK == IRECV) THEN
        KCOUNT=0
        DO IP=1,NPROC
          DO M=MINF,MSUP
            DO K=KINF,KSUP
              DO IJ=NBLKS(IP),NBLKE(IP)
                KCOUNT=KCOUNT+1
                FL(IJ,K,M)=PRECVBUF(KCOUNT)
              ENDDO
            ENDDO 
          ENDDO 
        ENDDO 
      ENDIF

      DEALLOCATE(PSENDBUF)
      DEALLOCATE(PRECVBUF)

      IF (LHOOK) CALL DR_HOOK('MPGATHERFL',1,ZHOOK_HANDLE)

      END SUBROUTINE MPGATHERFL
