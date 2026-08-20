! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE READMDLCONF (IU07)

! ----------------------------------------------------------------------

!**** *READMDLCONF* 

!*    PURPOSE.
!     --------

!     DETERMINES MODEL GRID CONFIGURATION 

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,SINPH    ,COSPH    ,    &
     &            IJS      ,IJL
      USE YOWMAP   , ONLY : BLK2GLO  ,NGX      ,NGY      ,NIBLO    ,    &
     &            IPER     ,IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,    &
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,NLONRGG  ,    &
     &            IQGAUSS
      USE YOWPCONS , ONLY : CIRC     ,RAD
      USE YOWSHAL  , ONLY : BATHY    ,LLOCEANMASK
      USE YOWTEST  , ONLY : IU06
      USE YOWABORT , ONLY : WAM_ABORT

      USE MPL_MODULE,ONLY : MPL_BARRIER
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "readpre.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU07
      INTEGER(KIND=JWIM) :: IREAD
      INTEGER(KIND=JWIM) :: IP, I, K

      REAL(KIND=JWRB), PARAMETER :: XLATMAX=87.5_JWRB
      REAL(KIND=JWRB) :: XLAT, ZLONS, COSPHMIN
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('READMDLCONF',0,ZHOOK_HANDLE)


      CALL MPL_BARRIER(CDSTRING='READMDLCONF:')
      CALL READPRE (IU07)
      CALL MPL_BARRIER(CDSTRING='READMDLCONF:')


!     DETERMINE THE TOTAL NUMBER OF OCEAN POINTS AND THE CONNECTION TO THE GRID (BLK2GLO)
      IF (ALLOCATED(LLOCEANMASK)) DEALLOCATE(LLOCEANMASK)
      ALLOCATE(LLOCEANMASK(NGX,NGY))
      NIBLO = 0
      DO K=1,NGY
        DO I=1,NLONRGG(K)
          IF (BATHY(I,K) > 0.0_JWRB) THEN
            NIBLO = NIBLO + 1
            LLOCEANMASK(I,K) = .TRUE.
          ELSE
            LLOCEANMASK(I,K) = .FALSE.
          ENDIF
        ENDDO
      ENDDO

      IF ( NIBLO <= 0 ) THEN
         WRITE(IU06,*)'***ERROR IN READREC: NIBLO SHOULD > 0, NIBLO = ',NIBLO
         CALL FLUSH(IU06)
         CALL ABORT1
      ENDIF

      IF (ALLOCATED(BLK2GLO%IXLG)) CALL BLK2GLO%DEALLOC
      CALL BLK2GLO%ALLOC(NIBLO)

      IP = 0
      DO K=1,NGY
        DO I=1,NLONRGG(K)
          IF (LLOCEANMASK(I,K)) THEN
            IP = IP+1
            BLK2GLO%IXLG(IP) = I
            BLK2GLO%KXLT(IP) = K
          ENDIF
        ENDDO
      ENDDO



!*    GRID INCREMENTS IN RADIANS AND METRES.
!     --------------------------------------

      IF(ALLOCATED(DELLAM) ) DEALLOCATE(DELLAM)
      ALLOCATE(DELLAM(NGY))
      IF(ALLOCATED(SINPH) ) DEALLOCATE(SINPH) 
      ALLOCATE(SINPH(NGY))
      IF(ALLOCATED(COSPH) ) DEALLOCATE(COSPH)
      ALLOCATE(COSPH(NGY))
      IF(ALLOCATED(ZDELLO) ) DEALLOCATE(ZDELLO)
      ALLOCATE(ZDELLO(NGY))

      DELPHI = XDELLA*CIRC/360.0_JWRB

      DO K=1,NGY
        XLAT       = (AMOSOP + REAL(K-1,JWRB)*XDELLA)*RAD
        SINPH(K)   = SIN(XLAT)
        COSPH(K)   = COS(XLAT)

        IF (NGX == 1 .AND. NGY == 1) THEN
          ZDELLO(K) = XDELLO
          DELLAM(K) = ZDELLO(K)*CIRC/360.0_JWRB
          EXIT
        ENDIF

        IF (IPER == 1) THEN
          ZDELLO(K) = 360.0_JWRB/REAL(NLONRGG(K),JWRB)
        ELSE
          ZDELLO(K) = (AMOEAP-AMOWEP)/REAL(NLONRGG(K)-1,JWRB)
        ENDIF
        DELLAM(K) = ZDELLO(K)*CIRC/360.0_JWRB
      ENDDO

!!!!! IF THE POLES ARE INCLUDED, ARTIFICIALLY REMOVE THE SINGULARITY
      COSPHMIN=COS(XLATMAX*RAD)
      DO K=1,NGY
        IF (COSPH(K) <= COSPHMIN) THEN
          COSPH(K)=COS(XLATMAX*RAD)
          SINPH(K)=SIN(XLATMAX*RAD)
        ENDIF
      ENDDO 

      IJS = 1
      IJL = NIBLO


      IF (LHOOK) CALL DR_HOOK('READMDLCONF',1,ZHOOK_HANDLE)

END SUBROUTINE READMDLCONF
