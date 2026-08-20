! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE INITGC

! ----------------------------------------------------------------------

!**** *INITGC* -

!*    PURPOSE.
!     ---------

!     INITIALISATION FOR THE GRAVITY-CAPILLARY PARAMTERISATION


!**   INTERFACE.
!     ----------

!       *CALL* *INITGC

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : NWAV_GC, KRATIO_GC, XKS_GC, XKL_GC,  &
     &                      XK_GC, XKM_GC, OMEGA_GC, OMXKM3_GC, VG_GC, C_GC, &
     &                      CM_GC, C2OSQRTVG_GC, XKMSQRTVGOC2_GC, OM3GMKM_GC,&
     &                      DELKCC_GC, DELKCC_GC_NS, DELKCC_OMXKM3_GC
      USE YOWPCONS , ONLY : G,  SURFT, SQRTGOSURFT

      USE YOMHOOK  ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: I

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INITGC',0,ZHOOK_HANDLE)

      SQRTGOSURFT = SQRT(G/SURFT)

      NWAV_GC = NINT(LOG(XKL_GC/XKS_GC)/(LOG(KRATIO_GC)))

      IF(.NOT.ALLOCATED(XK_GC)) ALLOCATE(XK_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(XKM_GC)) ALLOCATE(XKM_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(OMEGA_GC)) ALLOCATE(OMEGA_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(OMXKM3_GC)) ALLOCATE(OMXKM3_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(VG_GC)) ALLOCATE(VG_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(C_GC)) ALLOCATE(C_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(CM_GC)) ALLOCATE(CM_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(C2OSQRTVG_GC)) ALLOCATE(C2OSQRTVG_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(XKMSQRTVGOC2_GC)) ALLOCATE(XKMSQRTVGOC2_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(OM3GMKM_GC)) ALLOCATE(OM3GMKM_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(DELKCC_GC)) ALLOCATE(DELKCC_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(DELKCC_GC_NS)) ALLOCATE(DELKCC_GC_NS(NWAV_GC))
      IF(.NOT.ALLOCATED(DELKCC_OMXKM3_GC)) ALLOCATE(DELKCC_OMXKM3_GC(NWAV_GC))

!     COMPUTATION OF WAVENUMBER AND INCREMENT ARRAYS.
!     -------------------------------------------------
      DO I=1,NWAV_GC
        XK_GC(I) = XKS_GC*(KRATIO_GC)**(I-1)
        XKM_GC(I) = 1.0_JWRB/XK_GC(I)
        OMEGA_GC(I) = FOMEG_GC(XK_GC(I))
        OMXKM3_GC(I) = OMEGA_GC(I) * XKM_GC(I)**3
        VG_GC(I) = FVG_GC(XK_GC(I))
        C_GC(I) = FC_GC(XK_GC(I))      
        CM_GC(I) = 1.0_JWRB/C_GC(I)
        C2OSQRTVG_GC(I) = C_GC(I)**2 / SQRT(VG_GC(I))
        XKMSQRTVGOC2_GC(I) = XKM_GC(I) / C2OSQRTVG_GC(I)
        OM3GMKM_GC(I) = OMEGA_GC(I)**3 / (G*XK_GC(I))
      ENDDO

      DELKCC_GC(1) = 0.5*(XK_GC(2)-XK_GC(1))/C2OSQRTVG_GC(1)
      DELKCC_GC_NS(1) = DELKCC_GC(1)
      DO I=2,NWAV_GC-1
        DELKCC_GC(I) = 0.5_JWRB*(XK_GC(I+1)-XK_GC(I-1))/C2OSQRTVG_GC(I)
        DELKCC_GC_NS(I) = 0.5_JWRB*(XK_GC(I+1)-XK_GC(I))/C2OSQRTVG_GC(I)
      ENDDO
      DELKCC_GC(NWAV_GC) = 0.5_JWRB*(XK_GC(NWAV_GC)-XK_GC(NWAV_GC-1))/C2OSQRTVG_GC(NWAV_GC)
      DELKCC_GC_NS(NWAV_GC) = DELKCC_GC(NWAV_GC)

      DO I=1,NWAV_GC
        DELKCC_OMXKM3_GC(I) = DELKCC_GC(I) * OMXKM3_GC(I)
      ENDDO

      IF (LHOOK) CALL DR_HOOK('INITGC',1,ZHOOK_HANDLE)

      END SUBROUTINE INITGC
