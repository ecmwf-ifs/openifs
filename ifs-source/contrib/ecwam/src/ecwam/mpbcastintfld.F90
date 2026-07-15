! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MPBCASTINTFLD(ISEND, ITAG, NXF, NYF, IFIELD) 

!****  *MPBCASTINTFLD* - BROADCAST 2-D INTEGER DATA FIELD FROM 1 PE 
!                        TO ALL OTHERS.

!     J. BIDLOT    ECMWF   FEBRUARY 1999. 

!     PURPOSE.
!     --------

!     BROADCAST DATA CONTAINED IN 2-D ARRAY IFIELD FROM PE ISEND 
!     TO ALL OTHER PE's.

!*    INTERFACE.
!     ----------

!     CALL *MPBCASTINTFLD*(ITAG,NSTART,NEND,IFIELD) 

!     *ISEND*     RANK OF THE PROCESS FROM WHICH THE FIELD IS
!                 BROADCASTED
!     *ITAG*      TAG ASSOCIATED WITH A PARTICULAR CALL TO SUBROUTINE
!     *NXF*       FIRST DIMENSION OF IFIELD
!     *NYF*       SECOND DIMENSION OF IFIELD
!     *IFIELD*    INPUT/OUTPUT ARRAY CONTAINING THE FIELD

!     METHOD.
!     -------
!     MPL BROADCAST OF ARRAY FIELD FROM ISEND
!     AND RECEIVE OF THE CONTRIBUTION ON THE OTHER PE's.

!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_BROADCAST

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
 
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_BROADCAST

!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: ISEND, NXF, NYF
      INTEGER(KIND=JWIM), INTENT(INOUT) :: ITAG
      INTEGER(KIND=JWIM),DIMENSION(NXF,NYF), INTENT(INOUT) :: IFIELD

      INTEGER(KIND=JWIM) :: I, J, KCOUNT
      INTEGER(KIND=JWIM) :: MPLENGTH

      INTEGER(KIND=JWIM),ALLOCATABLE :: ICOMBUF(:)

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MPBCASTINTFLD',0,ZHOOK_HANDLE)

      MPLENGTH=NXF*NYF

      IF (ISEND.NE.0 .AND. NPROC.NE.1) THEN
        ALLOCATE(ICOMBUF(MPLENGTH))

        IF (IRANK.EQ.ISEND) THEN
          KCOUNT=0
          DO J=1,NYF
            DO I=1,NXF
              KCOUNT=KCOUNT+1
              ICOMBUF(KCOUNT)=IFIELD(I,J)
            ENDDO
          ENDDO
        ENDIF

        CALL GSTATS(617,0)
        CALL MPL_BROADCAST(ICOMBUF(1:MPLENGTH),KROOT=ISEND,             &
     &   KTAG=ITAG,CDSTRING='MPBCASTINTFLD:')
        CALL GSTATS(617,1)


        IF (IRANK.NE.ISEND) THEN
          KCOUNT=0
          DO J=1,NYF
            DO I=1,NXF
              KCOUNT=KCOUNT+1
              IFIELD(I,J)=ICOMBUF(KCOUNT)
            ENDDO
          ENDDO
        ENDIF

        DEALLOCATE(ICOMBUF)

      ENDIF

      IF (LHOOK) CALL DR_HOOK('MPBCASTINTFLD',1,ZHOOK_HANDLE)

      END SUBROUTINE MPBCASTINTFLD
