! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTEGER FUNCTION NS_GC (USTAR)

! ----------------------------------------------------------------------

!**** *NS_GC* - FUNCTION TO DETERMINE THE CUT-OFF ANGULAR FREQUENCY INDEX 
!               FOR THE GRAVITY-CAPILLARY MODEL
!               !!!! rounded to the closest index of XK_GC  !!!!!

!**   INTERFACE.
!     ----------

!       *FUNCTION* *NS_GC (USTAR)*

!       *USTAR*  - FRICTION VELOCITY. 

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWFRED  , ONLY : XLOGKRATIOM1_GC, NWAV_GC,  XKM_GC
USE YOWPCONS , ONLY : SQRTGOSURFT

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: USTAR

      REAL(KIND=JWRB) :: Y, XKS

! ----------------------------------------------------------------------


!!!Y = 1.0_JWRB/(1.48_JWRB+2.05_JWRB*UST)
!!!Y = (1.0_JWRB + UST**2)/(1.0_JWRB+10.0_JWRB*UST**2)

XKS = SQRTGOSURFT/(1.48_JWRB+2.05_JWRB*USTAR)

NS_GC = MIN( INT(LOG(MAX(XKS*XKM_GC(1), 1.0_JWRB))*XLOGKRATIOM1_GC) + 1, NWAV_GC-1)


END FUNCTION NS_GC
