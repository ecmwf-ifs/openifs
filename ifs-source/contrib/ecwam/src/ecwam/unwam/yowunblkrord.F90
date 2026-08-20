! (C) Copyright 2001- Aron Roland (Roland & Partner, Germany).
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE YOWUNBLKRORD
IMPLICIT NONE
CONTAINS

      SUBROUTINE UNBLKRORD(NLTOG,IJS,IJL,KINF,KSUP,MINF,MSUP,           &
     &                     BLOCK,BLOCK_G)

!****  *UNBLKRORD* - RE-ORDER BLOCK OF DATA FROM WHEN UNSTRUCTURED GRID 
!                    IS USED TO EITHER ITS GLOBAL CONFIGURATION
!                    i.e. independent of the model grid decomposition
!                    OR BACK TO THE LOCAL CONFIGURATION IN USE 
!                    i.e. dependent of the model grid decomposition used

!     PURPOSE.
!     --------

!     RE-ORDER GLOBAL BLOCK STRUCTURE WHEN UNSTRUCTURED GRID. 

!*    INTERFACE.
!     ----------

!     CALL *UNBLKRORD*(NLTOG,IJS,IJL,KINF,KSUP,MINF,MSUP,BLOCK,BLOCK_G)

!     *NLTOG*     SPECIFIES THE LOCAL ORDER IN INPOY OR OUTPUT 
!     *IJS*       INDEX OF THE FIRST FIRST DIMENSION OF BLOCK AND BLOCK_G
!     *IJL*       INDEX OF THE LAST  FIRST DIMENSION OF BLOCK AND BLOCK_G
!     *KINF*      INDEX OF THE FIRST SECOND DIMENSION OF BLOCK AND BLOCK_G
!     *KSUP*      INDEX OF THE LAST  SECOND DIMENSION OF BLOCK AND BLOCK_G
!     *MINF*      INDEX OF THE FIRST THIRD DIMENSION OF BLOCK AND BLOCK_G
!     *MSUP*      INDEX OF THE LAST THIRS DIMENSION OF BLOCK AND BLOCK_G
!     *BLOCK*     ARRAY CONTAINING THE BLOCK DATA WITH THE LOCAL ORDER
!                 INPUT WHEN NLTOG=1, OUTPUT WHEN NLTOG=-1 
!     *BLOCK_G*   ARRAY CONTAINING THE BLOCK DATA WITH THE GLOBAL ORDER
!                 OUTPUT WHEN NLTOG=1, INPUT WHEN NLTOG=-1 

!     METHOD.
!     -------

! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMPP   , ONLY : NPROC
      USE YOWTEST  , ONLY : IU06
      USE YOWPARAM , ONLY : LLUNSTR
      USE YOWPD, ONLY : RANK
      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK

!----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NLTOG,IJS,IJL,KINF,KSUP,MINF,MSUP
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,KINF:KSUP,MINF:MSUP), INTENT(INOUT) :: BLOCK
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,KINF:KSUP,MINF:MSUP), INTENT(INOUT) :: BLOCK_G

      INTEGER(KIND=JWIM) :: IR, IP, M, K, ISTART

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('UNBLKRORD',0,ZHOOK_HANDLE)

      IF(NPROC.EQ.1) THEN
        IF(NLTOG.EQ.1) THEN
          BLOCK_G(:,:,:)=BLOCK(:,:,:)
        ELSE IF(NLTOG.EQ.-1) THEN
          BLOCK(:,:,:)=BLOCK_G(:,:,:)
        ENDIF
      ELSE IF(LLUNSTR) THEN
        IF(SUM(RANK(:)%NP).NE.IJL) THEN
          WRITE(IU06,*)'*******************************************'
          WRITE(IU06,*)'*                                         *'
          WRITE(IU06,*)'* UNBLKRORD : INCONSISTENCY NIBLO AND SUM *' 
          WRITE(IU06,*)'* NOT EQUAL !!!                           *'
          WRITE(IU06,*)'* NIBLO (=IJL) = ',IJL
          WRITE(IU06,*)'* SUM(RANK(:)%NP)= ',SUM(RANK(:)%NP)
          WRITE(IU06,*)'*                                         *'
          WRITE(IU06,*)'*******************************************'
          CALL ABORT1
        ENDIF

!!! we will need to apply OMP on this !!!!!!!
!       LOCAL DISTRIBUTION TO GLOBAL DISTRIBUTION
        IF(NLTOG.EQ.1) THEN
          DO M=MINF,MSUP
            DO K=KINF,KSUP
              DO IR=1,NPROC
                ISTART=RANK(IR)%ISTART
                DO IP=1,RANK(IR)%NP
                  BLOCK_G(RANK(IR)%IPLG(IP),K,M)=BLOCK(ISTART+IP-1,K,M)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ELSE IF(NLTOG.EQ.-1) THEN
!       GLOBAL DISTRIBUTION TO LOCAL DISTRIBUTION
          DO M=MINF,MSUP
            DO K=KINF,KSUP
              DO IR=1,NPROC
                ISTART=RANK(IR)%ISTART
                DO IP=1,RANK(IR)%NP
                  BLOCK(ISTART+IP-1,K,M)=BLOCK_G(RANK(IR)%IPLG(IP),K,M)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ELSE
          WRITE(IU06,*)'**************************************'
          WRITE(IU06,*)'*                                    *'
          WRITE(IU06,*)'* UNBLKRORD : WRONG OPTION FOR NLTOG:*' 
          WRITE(IU06,*)'* NLTOG= ',NLTOG 
          WRITE(IU06,*)'*                                    *'
          WRITE(IU06,*)'**************************************'
          CALL ABORT1
        ENDIF

      ENDIF

      IF (LHOOK) CALL DR_HOOK('UNBLKRORD',1,ZHOOK_HANDLE)

      END SUBROUTINE UNBLKRORD

END MODULE YOWUNBLKRORD
