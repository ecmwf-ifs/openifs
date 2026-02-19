! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE GPUVS( &
 ! --- INPUT -----------------------------------------------------------------
 & KFLEV,KPROMA,KSTART,KEND,LDER, &
 & PUF,PVF, &
 ! --- OUTPUT ----------------------------------------------------------------
 & PUS,PVS, &
 ! --- OPTIONAL INPUT --------------------------------------------------------
 & PDIVF,PVORF,PUFL,PVFL, &
 ! --- OPTIONAL OUTPUT -------------------------------------------------------
 & PUS_L,PVS_L,PUS_M,PVS_M &
 & )

! GPUVS - Diagnoses "V_surf" and "grad(V_surf)".

! Purpose
! -------
!   Diagnoses "V_surf" and "grad(V_surf)" (surface wind).
!   The current assumptions which are done are:
!   * V_surf=V(l=KFLEV)
!   * grad(V_surf)=grad(V(l=KFLEV))
!   The consistency must be kept with the content of routine GPHLUV.

! Interface
! ---------
!   * INPUT:
!   KFLEV   - number of levels.
!   KPROMA  - horizontal dimension.
!   KSTART  - start of work.
!   KEND    - end of work.
!   LDER    - T: treatment of the derivatives.
!   PUF     - upper air U-wind at full levels.
!   PVF     - upper air V-wind at full levels.

!   * OUTPUT:
!   PUS     - surface U wind.
!   PVS     - surface V wind.

!   * OPTIONAL INPUT:
!   PDIVF   - upper air divergence at full levels.
!   PVORF   - upper air vorticity at full levels.
!   PUFL    - upper air zonal component of grad(U-wind) at full levels.
!   PVFL    - upper air zonal component of grad(V-wind) at full levels.

!   * OPTIONAL OUTPUT:
!   PUS_L   - zonal component of grad(surface U wind).
!   PVS_L   - zonal component of grad(surface V wind).
!   PUS_M   - meridian component of grad(surface U wind).
!   PVS_M   - meridian component of grad(surface V wind).

! Externals
! ---------

! Method
! ------

! Reference
! ---------

! Author
! ------
!   K. Yessad, Dec 2004 (after GNHPDVD, GNHGRP, GPHLUV)

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK

! -----------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)          :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN)          :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)          :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)          :: KEND 
LOGICAL           ,INTENT(IN)          :: LDER
REAL(KIND=JPRB)   ,INTENT(IN)          :: PUF(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)          :: PVF(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)         :: PUS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)         :: PVS(KPROMA) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PDIVF(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PVORF(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PUFL(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PVFL(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)   :: PUS_L(KPROMA) 
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)   :: PVS_L(KPROMA) 
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)   :: PUS_M(KPROMA) 
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)   :: PVS_M(KPROMA) 

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

#include "abor1.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPUVS',0,ZHOOK_HANDLE)

! -----------------------------------------------------------------------------

!*      1. COMPUTES "V_surf" and "grad(V_surf)".
!       ----------------------------------------

DO JROF=KSTART,KEND
  PUS(JROF)=PUF(JROF,KFLEV) 
  PVS(JROF)=PVF(JROF,KFLEV) 
ENDDO

IF (LDER) THEN
  IF (.NOT.(PRESENT(PDIVF).AND.PRESENT(PVORF).AND.PRESENT(PUFL) &
   & .AND.PRESENT(PVFL).AND.PRESENT(PUS_L).AND.PRESENT(PVS_L) &
   & .AND.PRESENT(PUS_M).AND.PRESENT(PVS_M))) &
   & CALL ABOR1(' GPUVS: LDER=T => PDIVF to PVS_M should be present!')
  DO JROF=KSTART,KEND
    PUS_L(JROF)=PUFL(JROF,KFLEV) 
    PVS_L(JROF)=PVFL(JROF,KFLEV) 
    PUS_M(JROF)=PVFL(JROF,KFLEV)-PVORF(JROF,KFLEV)
    PVS_M(JROF)=PDIVF(JROF,KFLEV)-PUFL(JROF,KFLEV)
  ENDDO
ENDIF

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPUVS',1,ZHOOK_HANDLE)

END SUBROUTINE GPUVS

