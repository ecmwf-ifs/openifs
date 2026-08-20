! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MPGATHERSCFLD(IRECV, NBLKS, NBLKE, FIELD, NLEN)

!****  *MPGATHERSCFLD* - GATHER SCALAR BLOCK DATA FIELD ONTO A SINGLE 
!****                    PROCESS 

!     J. BIDLOT    ECMWF   APRIL 1996  MESSAGE PASSING

!     PURPOSE.
!     --------

!     GATHER SCALAR BLOCK DATA CONTAINED IN ARRAY FIELD DISTRIBUTED
!     ACROSS THE DIFFERENT PROCESSES ONTO THE SINGLE PROCESS IRECV 

!*    INTERFACE.
!     ----------

!     CALL *MPGATHERSCFLD*(IRECV,NBLKS,NBLKE,FIELD,NLEN) 

!     *IRECV*     RANK OF THE PROCESS ONTO WHICH FIELD IS COLLECTED 
!     *NBLKS*     INDEX OF THE FIRST LOCAL POINT OF EACH PROCESSOR
!     *NBLKE*     INDEX OF THE LAST LOCAL POINT OF EACH PROCESSOR
!      IF NBLKE(IR) < NBLKS(IR) THEN IT IS ASSUMED THAT DATA VALUES ARE
!      ON PROCESSOR IR
!     *FIELD*     INPUT/OUTPUT ARRAY CONTAINING THE FIELD
!     *NLEN*      DIMENSION OF FIELD

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_GATHERV

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMPP   , ONLY : IRANK     ,NPROC
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_GATHERV, MPL_BARRIER

!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IRECV
      INTEGER(KIND=JWIM), INTENT(IN) :: NLEN
      INTEGER(KIND=JWIM), DIMENSION(NPROC), INTENT(IN) :: NBLKS, NBLKE

      REAL(KIND=JWRB), DIMENSION(NLEN), INTENT(INOUT) :: FIELD

      INTEGER(KIND=JWIM) :: IP, IJ
      INTEGER(KIND=JWIM) :: NST, NND
      INTEGER(KIND=JWIM), DIMENSION(NPROC) :: KRECVCOUNTS

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:) :: ZSENDBUF

!----------------------------------------------------------------------

        IF (IRECV.EQ.0 .OR. NPROC.EQ.1) RETURN

        IF (LHOOK) CALL DR_HOOK('MPGATHERSCFLD',0,ZHOOK_HANDLE)

        CALL GSTATS(674,0)

        DO IP=1,NPROC
          KRECVCOUNTS(IP)=MAX(NBLKE(IP)-NBLKS(IP)+1,0)
        ENDDO
        NST=NBLKS(IRANK)
        NND=MAX(NBLKS(IRANK),NBLKE(IRANK))

        ALLOCATE(ZSENDBUF(NST:NND))

        DO IJ=NBLKS(IRANK),NBLKE(IRANK)
          ZSENDBUF(IJ)=FIELD(IJ)
        ENDDO
        CALL MPL_GATHERV(ZSENDBUF(NST:NND),                             &
     &                   KRECVCOUNTS=KRECVCOUNTS,                       &
     &                   PRECVBUF=FIELD,                                &
     &                   KROOT=IRECV,                                   &
     &                   CDSTRING='MPGATHERSCFLD:')

        DEALLOCATE(ZSENDBUF)

        CALL MPL_BARRIER(CDSTRING='MPGATHERSCFLD:')

        CALL GSTATS(674,1)

        IF (LHOOK) CALL DR_HOOK('MPGATHERSCFLD',1,ZHOOK_HANDLE)

      END SUBROUTINE MPGATHERSCFLD
