! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MGRID (BATHY)

! ----------------------------------------------------------------------

!**** *MGRID* - ROUTINE TO ARRANGE WAMODEL GRID.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       TO ARRANGE WAMODEL GRID FOR A GIVEN AREA AND COMPUTE VARIOUS
!       MODEL CONSTANTS.

!**   INTERFACE.
!     ----------

!       *CALL* *MGRID (BATHY)*
!          *BATHY*     - BATHYMETRY /DATA OF PART

!     METHOD.
!     -------

!       THE NUMBER OF SEA POINTS PER LATITUDE IS COUNTED AND MODEL
!       BLOCKS OF MAXIMUM LENGTH OF NIBLO ARE CONSTRUCTED.

!     EXTERNALS.
!     ----------

!       *MBLOCK*    - SUB. TO GENERATE A BLOCK.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMAP   , ONLY : NGX      ,NGY      ,NIBLO       ,NLONRGG

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "mblock.intfb.h"

      REAL(KIND=JWRB), INTENT(INOUT) :: BATHY(NGX,NGY)

      INTEGER(KIND=JWIM) :: K, I, IL, KA, KE
      INTEGER(KIND=JWIM), DIMENSION(NGY) :: IPP

! ----------------------------------------------------------------------


!*    1. COUNT NUMBER OF SEA POINTS PER LATITUDE.
!        ----------------------------------------

      DO K=1,NGY
        IPP(K) = 0
        DO I=1,NLONRGG(K)
          IF (BATHY(I,K) > -990.0_JWRB) THEN
            IPP(K) = IPP(K) + 1
          ENDIF
        ENDDO
      ENDDO

      IF (NIBLO <= 0) THEN
        NIBLO=0
        DO K=1,NGY
          NIBLO=NIBLO+IPP(K)
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    2. MAKE BLOCKS.
!        ------------

      IL = 0
      KA = 1
      DO K = 1,NGY
        IL = IL + IPP(K)
        IF (IL > NIBLO) THEN
          KE = K-1
          CALL MBLOCK (BATHY, KA, KE, IPP)
          KA = KE-1
          IL = IPP(KA)+IPP(KE)+IPP(KE+1)
        ENDIF
      ENDDO
      CALL MBLOCK (BATHY, KA, NGY, IPP)

      END SUBROUTINE MGRID
