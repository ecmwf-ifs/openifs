! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MBLOCK (BATHY, KA, KE, IPP)

! ----------------------------------------------------------------------

!**** *MBLOCK* - ROUTINE TO ARRANGE WAMODEL GRID FOR ONE BLOCK.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       *MBLOCK* ARRANGES WAMODEL GRID FOR A BLOCK AND
!                COMPUTES VARIOUS MODEL CONSTANTS

!**   INTERFACE.
!     ----------

!       *CALL* *MBLOCK (BATHY, KA, KE, IPP)*
!          *BATHY*   - BATHYMETRY DATA.
!          *KA*      - NUMBER OF FIRST LAT IN BLOCK.
!          *KE*      - NUMBER OF LAST LAT IN BLOCK.
!          *IPP*     - NUMBER OF SEA POINTS PER LAT.

!     METHOD.
!     -------

!       THE LAND POINTS ARE REMOVED. ALL MODEL CONSTANTS WHICH ARE
!       GRID DEPENDENT ARE YOWPUTED AND STORED IN THE YOWMON BLOCKS.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWGRID  , ONLY : IJS      ,IJL
      USE YOWMAP   , ONLY : BLK2GLO  ,AMOSOP   ,XDELLA, NLONRGG,  &
     &                      NGX      ,NGY      ,NIBLO
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM) :: KA, KE
      INTEGER(KIND=JWIM) :: IC, IJ, IP, K, I 
      INTEGER(KIND=JWIM) :: IPP(NGY)

      REAL(KIND=JWRB) :: BATHY(NGX, NGY)

! ----------------------------------------------------------------------

!*    1. UPDATE BLOCK NUMBER AND INITIALIZES ARRAYS.
!        -------------------------------------------


      CALL BLK2GLO%ALLOC(NIBLO)

      DO IJ=1,NIBLO
        BLK2GLO%IXLG(IJ) = 0
        BLK2GLO%KXLT(IJ) = 0
      ENDDO

! ----------------------------------------------------------------------

!*    2. THE FIRST AND LAST BLOCK MUST CONTAIN MORE THAN 2
!*       ALL OTHER BLOCKS MORE  THAN 3 LATITUDES.
!        -------------------------------------------------

      IF (KA == 1 .AND. KE == 1 .AND. NGY == 1) THEN
        WRITE (IU06,*) '**********************************************'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '* ALLOWS FOR THE 1 GRID POINT MODEL          *'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '**********************************************'
      ELSEIF ((KE == 1) .OR. (KA == NGY) .OR.                            &
     &    ((KA /= 1) .AND. (KE == NGY) .AND. (KE-KA < 2))) THEN
        WRITE (IU06,*) '**********************************************'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '*        FATAL ERROR IN SUB. MBLOCK          *'
        WRITE (IU06,*) '*        ==========================          *'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '* BLOCK LENGTH IS TOO SHORT.                 *'
        WRITE (IU06,*) '* LESS THAN 2 LATITUDES IN FIRST OR LAST, OR *'
        WRITE (IU06,*) '* LESS THAN 3 LATITUDES IN OTHER BLOCKS.     *'
        WRITE (IU06,*) '* BLOCK LENGTH IS             NIBLO = ', NIBLO
        WRITE (IU06,*) '* NUMBER OF FIRST LATITUDE IS    KA = ', KA
        WRITE (IU06,*) '* NUMBER OF LAST  LATITUDE IS    KE = ', KE
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '*  PROGRAM WILL BE ABORTED                   *'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '**********************************************'
        CALL ABORT1
      ENDIF

! ----------------------------------------------------------------------

!*    3. COMPUTE INDICES OF FIRST, SECOND, BEFORE LAST, AND LAST LAT.
!        -----------------------------------------------------------

        IJS = 1
        IJL = NIBLO 

! ----------------------------------------------------------------------

!*    4. REMOVE LAND POINTS AND STORE COS AND SIN OF LAT.
!        ------------------------------------------------

      IP = 0
      DO K=KA,KE
        DO I=1,NLONRGG(K)
          IF (BATHY(I,K) > -990.0_JWRB) THEN
            IP = IP+1
            BLK2GLO%IXLG(IP) = I
            BLK2GLO%KXLT(IP) = K
          ENDIF
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    5. PRINTER PROTOCOL OF BLOCK.
!        --------------------------

        WRITE (IU06,'(1H0,'' BLOCKING INFORMATION:'')')
        WRITE (IU06,'(1H ,''            LATITUDES   '',                 &
     &   ''   SECOND LAT. INDEX '',                                     &
     &   '' SECOND TO LAST LAT  '',                                     &
     &   ''   TOTAL'')')
        WRITE (IU06,'(1H ,''  NO     SOUTH     NORTH'',                 &
     &   ''     START       END'',                                      &
     &   ''     START       END'',                                      &
     &   ''    POINTS'')')
      WRITE (IU06,'(2F10.2,5I10)')                                &
     &        AMOSOP+(KA-1)*XDELLA, AMOSOP+(KE-1)*XDELLA,          &
     &        IJS, IJL

      END SUBROUTINE MBLOCK
