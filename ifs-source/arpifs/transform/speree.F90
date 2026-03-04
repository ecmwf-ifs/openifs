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

SUBROUTINE SPEREE(YDGEOMETRY,KFLSUR,KFIELD,PSPEC,PREEL,PREELL,PREELM)

!**** *SPEREE* -Inverse transform to Grid-point space

!     Purpose.
!     --------
!     Transform from Spectral space to grid-point space.
!     Optionally computes horizontal derivatives.

!**   Interface.
!     ----------
!        *CALL* *SPEREE(...)

!        Explicit arguments :
!        --------------------
!           KFLSUR   : first dimension of PSPEC
!           KFIELD   : number of fields to be transformed
!           PSPEC    : spectral input array
!           PREEL    : grid point output array
!           PREELL   : grid point output array: zonal derivative (optional)
!           PREELM   : grid point output array: merid derivative (optional)

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.  INV_TRANS - inverse transform (TRANS lib)
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      MPP Group *ECMWF*
!      Original : 95-10-01

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K.Yessad (Jan 2011): optional calculation of derivatives
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : MYSETV, NPRTRV

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLSUR 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELD 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPEC(KFLSUR,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PREEL(YDGEOMETRY%YRGEM%NGPTOT,KFIELD,1)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PREELL(YDGEOMETRY%YRGEM%NGPTOT,KFIELD,1)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PREELM(YDGEOMETRY%YRGEM%NGPTOT,KFIELD,1)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ILF,ISTG,IFLDPP,IREST,IGF
INTEGER(KIND=JPIM) :: INUMLF(NPRTRV),IVSETSC(KFIELD)
INTEGER(KIND=JPIM) :: JROCB
REAL(KIND=JPRB),ALLOCATABLE :: ZREEL(:,:,:)
LOGICAL :: LLDER
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "inv_trans.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPEREE',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NRESOL=>YDDIM%NRESOL, NSPEC2=>YDDIM%NSPEC2, &
 & NGPTOT=>YDGEM%NGPTOT)
!     ------------------------------------------------------------------

!  INITIALISATION

ILF = 0
IF( KFIELD == 1 )THEN
  IVSETSC(1) = NPRTRV
  IF(NPRTRV == MYSETV) THEN
    ILF = 1
  ENDIF
  ISTG=2
ELSE
  IFLDPP=KFIELD/NPRTRV
  IREST=KFIELD-IFLDPP*NPRTRV
  DO JROCB=1,NPRTRV
    IF( JROCB <= IREST )THEN
      INUMLF(JROCB)=IFLDPP+1
    ELSE
      INUMLF(JROCB)=IFLDPP
    ENDIF
  ENDDO
  ILF = INUMLF(MYSETV)
  ISTG = 1
  DO JROCB=1,NPRTRV
    IVSETSC(ISTG:ISTG+INUMLF(JROCB)-1) = JROCB
    ISTG = ISTG+INUMLF(JROCB)
  ENDDO
ENDIF
IGF = ISTG-1

LLDER=PRESENT(PREELL).AND.PRESENT(PREELM)

IF (LLDER) THEN
  ALLOCATE(ZREEL(NGPTOT,3*KFIELD,1))
  CALL INV_TRANS(PSPSCALAR=PSPEC(1:ILF,:),LDSCDERS=.TRUE.,&
   & KRESOL=NRESOL,KPROMA=NGPTOT,KVSETSC=IVSETSC(1:IGF),PGP=ZREEL)  

  PREEL (1:NGPTOT,1:KFIELD,1)=ZREEL(1:NGPTOT,1:KFIELD,1)
  PREELM(1:NGPTOT,1:KFIELD,1)=ZREEL(1:NGPTOT,KFIELD+1:2*KFIELD,1)
  PREELL(1:NGPTOT,1:KFIELD,1)=ZREEL(1:NGPTOT,2*KFIELD+1:3*KFIELD,1)
  DEALLOCATE(ZREEL)
ELSE
  CALL INV_TRANS(PSPSCALAR=PSPEC(1:ILF,:),&
   & KRESOL=NRESOL,KPROMA=NGPTOT,KVSETSC=IVSETSC(1:IGF),&
   & PGP=PREEL)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPEREE',1,ZHOOK_HANDLE)
END SUBROUTINE SPEREE

