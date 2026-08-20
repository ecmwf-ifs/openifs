! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MEANSQS_GC(XKMSS, KIJS, KIJL, HALP, USTAR, XMSSCG, FRGC)

!***  DETERMINE MSS FOR GRAV-CAP WAVES UP TO WAVE NUMBER XKMSS

!     AUTHOR: PETER JANSSEN
!     ------

!     REFERENCES:
!     ----------

!     VIERS PAPER EQ.(29)

!----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : NWAV_GC, XLOGKRATIOM1_GC, XKM_GC,           &
     &                      VG_GC, C2OSQRTVG_GC, DELKCC_GC, DELKCC_GC_NS
      USE YOWPCONS , ONLY : G, ZPI,  SURFT

      USE YOMHOOK  ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!----------------------------------------------------------------------

      IMPLICIT NONE

#include "omegagc.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), INTENT(IN) :: XKMSS ! WAVE NUMBER CUT-OFF
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: HALP  ! 1/2 Phillips parameter
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: USTAR ! friction velocity
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: XMSSCG  ! mean square slope for gravity-capillary waves
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: FRGC  ! Frequency from which the gravity-capillary spectrum is approximated


      INTEGER(KIND=JWIM) :: IJ, I, NE
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL) :: NS
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: XKS, OMS, COEF
   
!     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MEANSQS_GC',0,ZHOOK_HANDLE)

      NE = MIN(MAX(NINT(LOG(XKMSS*XKM_GC(1))*XLOGKRATIOM1_GC ), 1), NWAV_GC)

      CALL OMEGAGC(KIJS, KIJL, USTAR, NS, XKS, OMS)

      DO IJ = KIJS, KIJL
        FRGC(IJ) = OMS(IJ)/ZPI
        IF(XKS(IJ) > XKMSS) THEN
          NS(IJ) = NE
          XMSSCG(IJ) = 0.0_JWRB
        ELSE
          XMSSCG(IJ) = DELKCC_GC_NS(NS(IJ)) * XKM_GC(NS(IJ)) 
        ENDIF
      ENDDO

      DO IJ = KIJS, KIJL
        DO I = NS(IJ)+1, NE 
!         ANALYTICAL FORM INERTIAL SUB RANGE F(k) = k**(-4)*BB
!         BB = COEF(IJ)*SQRT(VG_GC(I))/C_GC(I)**2
!         mss :  integral of k**2 F(k)  k dk
          XMSSCG(IJ) = XMSSCG(IJ) + DELKCC_GC(I) * XKM_GC(I) 
        ENDDO
        COEF(IJ) = C2OSQRTVG_GC(NS(IJ))*HALP(IJ)
        XMSSCG(IJ) = XMSSCG(IJ)*COEF(IJ)
      ENDDO

      IF (LHOOK) CALL DR_HOOK('MEANSQS_GC',1,ZHOOK_HANDLE)
 
      END SUBROUTINE MEANSQS_GC
