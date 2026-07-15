! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE MPMINMAXAVG(LLGLOBAL, IRECV, LDREPROD, BOUT, WNORM)
! ----------------------------------------------------------------------

!****  *MPMINMAXAVG* - FIND GLOBAL MIN, MAX AND AVERAGE OF BOUT

!           LLGLOBAL -  IF TRUE NORMS FOR SELECTED OUPTUT FIELDS WILL BE GLOBAL
!                       ON PROCESSOR IRECV
!           LDREPROD -  only if LLGLOBAL is false
!                       Reproducibility flag for SUMmation-operator.
!                       Meaningful only for REAL-numbers.
!                       Three modes (applicable for REAL-number only):
!                       1) Not provided at all (the default) ==> MPL_ABORT
!                       2) Provided and .TRUE. ==> Use home-written binary tree
!                          No MPI_ALLREDUCE used.
!                       3) Provided, but .FALSE. ==> let MPI_ALLREDUCE do the summation.
!           *BOUT*    - OUTPUT PARAMETERS BUFFER
!           *WNORM(1,:) : MINIMUM
!           *WNORM(2,:) : MAXIMUM
!           *WNORM(3,:) : AVERAGE
!           *WNORM(4,:) : NUMBER OF VALUES USED TO PRODUCE THE AVERAGE


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT   , ONLY : NFLAG    ,JPPFLAG , ITOBOUT ,NIPRMOUT
      USE YOWMAP    , ONLY : NIBLO
      USE YOWMPP    , ONLY : IRANK    ,NPROC
      USE YOWPARAM  , ONLY : LL1D    ,LLUNSTR
      USE YOWPCONS  , ONLY : ZMISS
      USE YOWGRID   , ONLY : NPROMA_WAM, NCHNK, IJFROMCHNK, KIJL4CHNK
      USE YOWSPEC   , ONLY : NBLKS    ,NBLKE   ,IJ2NEWIJ
      USE YOWTEST   , ONLY : IU06

      USE MPL_MODULE, ONLY : MPL_ALLREDUCE
      USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "abort1.intfb.h"
#include "mpgatherscfld.intfb.h"

      LOGICAL, INTENT(IN) :: LLGLOBAL
      INTEGER(KIND=JWIM), INTENT(IN) :: IRECV
      LOGICAL, INTENT(IN) :: LDREPROD
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NIPRMOUT, NCHNK), INTENT(IN) :: BOUT
      REAL(KIND=JWRB),DIMENSION(4, NIPRMOUT), INTENT(OUT) :: WNORM

      INTEGER(KIND=JWIM) :: IJ, IJOLD, I, IP, ITG, IT, ICOUNT
      INTEGER(KIND=JWIM) :: IJSG, IJLG, IJSB, IJLB, KIJS, KIJL, ICHNK, IPRM

      REAL(KIND=JWRB),DIMENSION(NIPRMOUT) :: ZMIN, ZMAX
      REAL(KIND=JWRB),DIMENSION(2*NIPRMOUT) :: ZSUM
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), ALLOCATABLE :: ZGLOBAL(:)
!     -------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MPMINMAXAVG',0,ZHOOK_HANDLE)

      IJSG = IJFROMCHNK(1,1)
      IJLG = IJSG + SUM(KIJL4CHNK) - 1

      ZSUM(:) = 0.0_JWRB
      ZMIN(:) = HUGE(ZMIN(:))
      ZMAX(:) = -HUGE(ZMAX(:))

      ICOUNT = 0

      IF (LLGLOBAL) THEN

        IF (LLUNSTR) THEN
          WRITE(IU06,*) '************************************'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '*  FATAL ERROR IN SUB. MPMINMAXAVG *'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '* GLOBAL NORM NOT YET CODED        *'
          WRITE(IU06,*) '* FOR UNSTRUCTURED GRID            *'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '************************************'
          CALL ABORT1
        ENDIF


        WNORM(:,:) = 0.0_JWRB

        IF (IJSG /= NBLKS(IRANK) .OR. IJLG /= NBLKE(IRANK) ) THEN
          WRITE(IU06,*) '************************************'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '*  FATAL ERROR IN SUB. MPMINMAXAVG *'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '* IJS /= NBLKS(IRANK) .OR.      *'
          WRITE(IU06,*)  IJSG, NBLKS(IRANK)
          WRITE(IU06,*) '* IJL /= NBLKE(IRANK)           *'
          WRITE(IU06,*)  IJLG, NBLKE(IRANK)
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '* THIS SHOULD NOT HAPPEN  !        *'
          WRITE(IU06,*) '************************************'
          CALL ABORT1
        ENDIF

        ALLOCATE(ZGLOBAL(NIBLO))
        ZGLOBAL(:)=ZMISS

        DO ITG=1,JPPFLAG
          IF (NFLAG(ITG)) THEN
            IT = ITOBOUT(ITG)

!$OMP       PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, IJSB, KIJL, IJLB)
            DO ICHNK = 1, NCHNK
              KIJS = 1
              IJSB = IJFROMCHNK(KIJS,ICHNK)
              KIJL = KIJL4CHNK(ICHNK)
              IJLB = IJFROMCHNK(KIJL,ICHNK)

              ZGLOBAL(IJSB:IJLB) = BOUT(KIJS:KIJL, IT, ICHNK)
            ENDDO
!$OMP       END PARALLEL DO

            CALL MPGATHERSCFLD(IRECV, NBLKS, NBLKE, ZGLOBAL, NIBLO)

            IF (IRANK == IRECV) THEN
              DO IJOLD = 1, NIBLO
                IF (LL1D .OR. LLUNSTR .OR. NPROC == 1) THEN
                  IJ=IJOLD
                ELSE
                  IJ=IJ2NEWIJ(IJOLD)
                ENDIF

                IF (ZGLOBAL(IJ) /= ZMISS) THEN
                  ICOUNT = ICOUNT + 1
                  ZSUM(IT) = ZSUM(IT) + ZGLOBAL(IJ)
                  ZMIN(IT) = MIN(ZMIN(IT), ZGLOBAL(IJ))
                  ZMAX(IT) = MAX(ZMAX(IT), ZGLOBAL(IJ))
                ENDIF
              ENDDO
              WNORM(1,IT)=ZSUM(IT)/MAX(ICOUNT,1)
              WNORM(2,IT)=ZMIN(IT)
              WNORM(3,IT)=ZMAX(IT)
              WNORM(4,IT)=ICOUNT
            ENDIF

          ENDIF
        ENDDO

        DEALLOCATE(ZGLOBAL)

      ELSE

!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(I, ICHNK, IPRM)
        DO I = 1, NIPRMOUT
!!!     outer loop on I !!!!

          DO ICHNK = 1, NCHNK
            DO IPRM = 1, KIJL4CHNK(ICHNK) 
              IF (BOUT(IPRM, I, ICHNK) /= ZMISS) THEN
                ZSUM(I) = ZSUM(I) + BOUT(IPRM, I, ICHNK)
                ZSUM(NIPRMOUT+I) = ZSUM(NIPRMOUT+I) + 1.0_JWRB
                ZMIN(I) = MIN(ZMIN(I), BOUT(IPRM, I, ICHNK))
                ZMAX(I) = MAX(ZMAX(I), BOUT(IPRM, I, ICHNK))
              ENDIF
            ENDDO
          ENDDO

        ENDDO
!$OMP   END PARALLEL DO

        IF (NPROC > 1) THEN
          CALL MPL_ALLREDUCE(ZSUM,'SUM',LDREPROD=LDREPROD, CDSTRING='MPMINMAXAVG:')
          CALL MPL_ALLREDUCE(ZMIN,'MIN',CDSTRING='MPMINMAXAVG VALMIN:')
          CALL MPL_ALLREDUCE(ZMAX,'MAX',CDSTRING='MPMINMAXAVG VALMAX:')
        ENDIF

        DO I=1,NIPRMOUT
          WNORM(2,I) = ZMIN(I)
          WNORM(3,I) = ZMAX(I)
          WNORM(4,I) = ZSUM(NIPRMOUT+I)
          IF (WNORM(4,I) < 1.0_JWRB) THEN
            WNORM(1,I) = -HUGE(WNORM(3,I))
          ELSE
            WNORM(1,I) = ZSUM(I)/WNORM(4,I)
          ENDIF
        ENDDO

      ENDIF
     
IF (LHOOK) CALL DR_HOOK('MPMINMAXAVG',1,ZHOOK_HANDLE)

END SUBROUTINE MPMINMAXAVG
