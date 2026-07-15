! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MFREDIR

! ----------------------------------------------------------------------

!**** *MFREDIR* - ROUTINE TO COMPUTE FREQUENCY DIRECTION CONSTANTS.

!     H. GUNTHER   ECMWF   2/4/90.
!     B. HANSEN    ECMWF   4/11/97. ACTIVATE THE ROTATION OF SPECTRA.

!*    PURPOSE.
!     --------

!       INITIATES THE FREQUENCY AND DIRECTION CONSTANTS WHICH ARE
!       SAVED IN COMMON FREDIR.

!**   INTERFACE.
!     ----------

!       *CALL* *MFREDIR*

!     METHOD.
!     -------

!       STARTING FROM THE FIRST FREQUENCY THE NEXT  ARE INCREMENTED
!       BY THE FACTOR *FRATIO* (SEE PARAMETER STATEMENT IN YOWFRED).
!       THE DIRECTIONS ARE EQUALLY DISTRIBUTED OVER THE CIRCLE
!       THE CIRCLE STARTING FROM 0.5 DELTH. 

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED
      USE YOWFRED  , ONLY : IFRE1    ,FR1      ,FRATIO   ,              & 
     &            FR       ,DFIM     ,GOM      ,C        ,              &
     &            DELTH    ,TH       ,COSTH    ,SINTH
      USE YOWPCONS , ONLY : G        ,PI       ,ZPI      ,DEG
      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "mfr.intfb.h"

      INTEGER(KIND=JWIM) :: M, K

      REAL(KIND=JWRB) :: CO1
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MFREDIR',0,ZHOOK_HANDLE)

      IF (.NOT.ALLOCATED(FR)) ALLOCATE(FR(NFRE))
      IF (.NOT.ALLOCATED(DFIM)) ALLOCATE(DFIM(NFRE))
      IF (.NOT.ALLOCATED(GOM)) ALLOCATE(GOM(NFRE))
      IF (.NOT.ALLOCATED(C)) ALLOCATE(C(NFRE))
      IF (.NOT.ALLOCATED(TH)) ALLOCATE(TH(NANG))
      IF (.NOT.ALLOCATED(COSTH)) ALLOCATE(COSTH(NANG))
      IF (.NOT.ALLOCATED(SINTH)) ALLOCATE(SINTH(NANG))


!*    1. FREQUENCY DEPENDENT CONSTANTS.
!        ------------------------------

!*    1.1 COMPUTE FREQUENCIES.
!         --------------------

      CALL MFR(NFRE, IFRE1, FR1, FRATIO, FR)  


!*    1.2 COMPUTE DEEP WATER GROUP VELOCITIES.
!        ------------------------------------

      DO M=1,NFRE
        GOM(M) = G/(4.0_JWRB*PI*FR(M))
      ENDDO

!*    1.3 COMPUTE PHASE VELOCITY IN DEEP WATER.
!         -------------------------------------

      DO M = 1,NFRE
        C(M) = G/(ZPI*FR(M))
      ENDDO

! ----------------------------------------------------------------------

!*    2. COMPUTATION OF DIRECTIONS, BANDWIDTH, SIN AND COS.
!        --------------------------------------------------

      DELTH = ZPI/REAL(NANG,JWRB)
      DO K=1,NANG
        TH(K) = REAL(K-1,JWRB)*DELTH + 0.5_JWRB*DELTH
        COSTH(K) = COS(TH(K))
        SINTH(K) = SIN(TH(K))
      ENDDO

! ----------------------------------------------------------------------

!*    3. COMPUTATION FREQUENCY DIRECTION AREAS
!        -------------------------------------

      CO1 = 0.5_JWRB*(FRATIO-1.0_JWRB)*DELTH
      DFIM(1)= CO1*FR(1)
      DO M=2,NFRE-1
        DFIM(M)=CO1 * (FR(M)+FR(M-1))
      ENDDO
      DFIM(NFRE)=CO1*FR(NFRE-1)

! ----------------------------------------------------------------------

!*    4. PRINTER PROTOCOL
!         ---------------

      WRITE (IU06,*) ' '
      WRITE (IU06,'(''  FREQUENCY AND DIRECTION GRID'')')
      WRITE (IU06,'(''  NUMBER OF FREQUENCIES IS  NFRE = '',I3)') NFRE
      WRITE (IU06,'(''  REDUCED NUMBER OF FREQUENCIES IS  NFRE_RED = '',I3)') NFRE_RED
      WRITE (IU06,'(''  NUMBER OF DIRECTIONS  IS  NANG = '',I3)') NANG
      WRITE (IU06,'(''  MODEL FREQUENCIES IN HERTZ:'')')
      WRITE (IU06,'(1X,13F10.5)') (FR(M),M=1,NFRE)
      WRITE (IU06,'(''  MODEL FREQUENCY INTERVALS TIMES DIRECTION'',    &
     &              '' INTERVAL IN HERTZ*RADIANS'')')
      WRITE (IU06,'(1X,13F10.5)') (DFIM(M),M=1,NFRE)
      WRITE (IU06,'(''  MODEL DEEP WATER GROUP VELOCITY IN M/S:'')')
      WRITE (IU06,'(1X,13F10.5)') (GOM(M),M=1,NFRE)
      WRITE (IU06,'(''  MODEL DEEP WATER PHASE VELOCITY IN M/S:'')')
      WRITE (IU06,'(1X,13F10.5)') (C(M),M=1,NFRE)
      WRITE (IU06,'(''  MODEL DIRECTIONS IN DEGREE'',                    &
     &              '' (CLOCKWISE FROM NORTH):'')')
      WRITE (IU06,'(1X,13F10.5)') (TH(K)*DEG,K=1,NANG)
      WRITE (IU06,*) ' '

      IF (LHOOK) CALL DR_HOOK('MFREDIR',1,ZHOOK_HANDLE)

      END SUBROUTINE MFREDIR
