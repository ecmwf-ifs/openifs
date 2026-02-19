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

SUBROUTINE PRE_GRID_BICONSERV(YDGEOMETRY,KGPOUT,KLATS,KLONS,PLATS,&
                            & PSPIN,PSPOUT,PSPSP,PGPSP,PGPSPG)

!     Purpose.
!     --------
!       Compute pressure fields needed for conserving interpolation.

!     Arguments.
!     ----------
!       KGPOUT  : Number of points at output resolution
!       KLATS   : 
!       KLONS   : 
!       PLATS   : 
!       PSPSP   : LnPs spectral field            (PSPSP or PGPSP or PGPSPG required)
!       PGPSP   : LnPs GP field in NPROMA blocks (PSPSP or PGPSP or PGPSPG required)
!       PGPSPG  : LnPs GP global field NGPTOTG size (PSPSP or PGPSP or PGPSPG required)
!       PSPIN : surface pressure at input  res. (for conserving 3D interp)
!       PSPOUT: surface pressure at output res. (for conserving 3D interp)

!     Author.
!     -------
!       Y. Tremolet (from E.Holm and L.Isaksen code)

!     Modifications.
!     --------------
!       Original : 26-Jan-2005
!       080204 E.Holm    : Use only surface pressure
!       G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
! ----------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1, ONLY: JPIM, JPRB
USE YOMHOOK , ONLY: LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE(GEOMETRY)    , INTENT(IN) :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN) :: KGPOUT
INTEGER(KIND=JPIM), INTENT(IN) :: KLATS
INTEGER(KIND=JPIM), INTENT(IN) :: KLONS(KLATS)
REAL(KIND=JPRB), INTENT(IN)    :: PLATS(KLATS)
REAL(KIND=JPRB), INTENT(OUT)   :: PSPIN(YDGEOMETRY%YRGEM%NGPTOTG)
REAL(KIND=JPRB), INTENT(OUT)   :: PSPOUT(KGPOUT)
REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PSPSP(YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PGPSP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PGPSPG(YDGEOMETRY%YRGEM%NGPTOTG)

INTEGER(KIND=JPIM) :: JJ,JG,JB,IEND
REAL(KIND=JPRB) :: ZPRESSIN(YDGEOMETRY%YRGEM%NGPTOTG,1),ZPRESSLOC(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
#include "gathergpf.intfb.h"
#include "speree.intfb.h"
#include "grid_biconserv.intfb.h"
! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PRE_GRID_BICONSERV',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NDGLG=>YDDIM%NDGLG, NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, &
 & NSPEC2=>YDDIM%NSPEC2, &
 & NGPTOT=>YDGEM%NGPTOT, NGPTOTG=>YDGEM%NGPTOTG, NLOENG=>YDGEM%NLOENG)

IF (.NOT. (PRESENT(PSPSP).OR.PRESENT(PGPSP).OR.PRESENT(PGPSPG)))&
  & CALL ABOR1('PRE_GRID_BICONSERV: PSPSP, PGPSP or PGPSPG required')

! Surface pressure at both input (size=NGPTOTG) and output (size=KGPOUT)
! resolution is required for conserving interpolation method, 
! so these fields are distributed from the io processor
! to all other processors.
! It is assumed that input and output has the same vertical level
! definitions, i.e. VAH and VBH are the same in the two GPPREH calls.

! HR pressure in grid point space required for conserving interpolation of
! 3D fields.


IF (PRESENT(PGPSPG)) THEN
  ZPRESSIN(:,1)=EXP(PGPSPG(:))
ELSE
  ! Transform HR lnPs to local gridpoint space
  ! and convert from ln(ps) to ps.
  IF (PRESENT(PSPSP)) THEN
    CALL SPEREE(YDGEOMETRY,1,1,PSPSP,ZPRESSLOC)
    ZPRESSLOC(:)=EXP(ZPRESSLOC(:))
  ELSEIF (PRESENT(PGPSP)) THEN
    JG=0
    DO JB=1,NGPBLKS
      IEND=MIN(NPROMA,NGPTOT-NPROMA*(JB-1))
      DO JJ=1,IEND
        JG=JG+1
        ZPRESSLOC(JG)=EXP(PGPSP(JJ,JB))
      ENDDO
    ENDDO
  ENDIF
  ! Gather local grid point sections to a global field
  ! and distribute to all processors
  CALL GATHERGPF(YDGEOMETRY,ZPRESSLOC,ZPRESSIN,1,-1)
ENDIF
PSPIN(:)=ZPRESSIN(:,1)

CALL GRID_BICONSERV(YDGEOMETRY%YRVAB%VAH,YDGEOMETRY%YRVAB%VBH, NDGLG,NLOENG(1:NDGLG),YDGEOMETRY%YRCSGLEG%RLATIG(1:NDGLG),&
       & NGPTOTG,KLATS,KLONS,PLATS,KGPOUT,PSPIN,PSPOUT)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('PRE_GRID_BICONSERV',1,ZHOOK_HANDLE)
END SUBROUTINE PRE_GRID_BICONSERV
