! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE INITDPTHFLDS(WVENVI, WVPRPT, WVPRPT_LAND)
! ----------------------------------------------------------------------

!**** *INITDPTHFLDS* - 

!*    PURPOSE.
!     --------

!     INITIALISE GRID POINT FIELDS DEPENDENT ON WATER DEPTH:
!     EMAXDPT

!     INITIALISE GRID POINT FIELDS DEPENDENT ON WATER DEPTH AND FREQUENCY:

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY, FREQUENCY_LAND

      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWFRED  , ONLY : ZPIFR
      USE YOWSHAL  , ONLY : GAM_B_J, BATHYMAX 

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "depthprpt.intfb.h"

      TYPE(ENVIRONMENT), INTENT(INOUT) :: WVENVI
      TYPE(FREQUENCY), INTENT(INOUT) :: WVPRPT
      TYPE(FREQUENCY_LAND), INTENT(INOUT) :: WVPRPT_LAND


      INTEGER(KIND=JWIM) :: IJ, ICHNK, KIJS, KIJL

      REAL(KIND=JWRB) :: GAM
      REAL(KIND=JWRB), DIMENSION(1) :: DEEP_DEPTH
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INITDPTHFLDS',0,ZHOOK_HANDLE)

      CALL GSTATS(1504,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK, KIJS, KIJL, IJ, GAM)
      DO ICHNK = 1, NCHNK
        KIJS = 1
        KIJL = NPROMA_WAM

!       COMPUTE THE MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
        DO IJ = KIJS, KIJL
!         REDUCE GAMMA FOR SMALL DEPTH ( < 4m)
!         (might need to be revisted when grid is fine resolution)
          IF (WVENVI%DEPTH(IJ,ICHNK) < 4.0_JWRB) THEN
            GAM=GAM_B_J*WVENVI%DEPTH(IJ,ICHNK)/4.0_JWRB
          ELSE
            GAM=GAM_B_J
          ENDIF
          WVENVI%EMAXDPT(IJ,ICHNK) = 0.0625_JWRB*(GAM*WVENVI%DEPTH(IJ,ICHNK))**2
        ENDDO

        CALL DEPTHPRPT (1, NPROMA_WAM, WVENVI%DEPTH(:,ICHNK),                                            &
 &                      WVPRPT%WAVNUM(:,:,ICHNK), WVPRPT%CINV(:,:,ICHNK), WVPRPT%CGROUP(:,:,ICHNK),      &
 &                      WVPRPT%XK2CG(:,:,ICHNK), WVPRPT%OMOSNH2KD(:,:,ICHNK), WVPRPT%STOKFAC(:,:,ICHNK) )

      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1504,1)


      ! Fictitious values for land point (NSUP+1)

      DEEP_DEPTH(1) = BATHYMAX

      CALL DEPTHPRPT (1, 1, DEEP_DEPTH(:),                                                      &
 &                    WVPRPT_LAND%WAVNUM(:), WVPRPT_LAND%CINV(:), WVPRPT_LAND%CGROUP(:),        &
 &                    WVPRPT_LAND%XK2CG(:), WVPRPT_LAND%OMOSNH2KD(:), WVPRPT_LAND%STOKFAC(:) )

      WVPRPT_LAND%CIWA(:) = 1.0_JWRB

IF (LHOOK) CALL DR_HOOK('INITDPTHFLDS',1,ZHOOK_HANDLE)

END SUBROUTINE INITDPTHFLDS
