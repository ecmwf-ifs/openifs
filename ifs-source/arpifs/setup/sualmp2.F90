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

SUBROUTINE SUALMP2(YDGEOMETRY)

!**** *SUALMP2 * - Routine to allocate space for some distributed memory arrays

!     Purpose.
!     --------
!           Allocate space for MYLEVS, MYLATS

!**   Interface.
!     ----------
!        *CALL* *SUALMP2*

!     Explicit arguments :  None
!     --------------------

!     Implicit arguments :
!     --------------------

!     Method.
!     -------
!     Externals.
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
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK,DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : NPRINTLEV
USE YOMCT0   , ONLY : LALLOPR
USE YOMLUN   , ONLY : NULOUT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
INTEGER(KIND=JPIM) :: IU
LOGICAL :: LLP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUALMP2',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NFLEVL=>YDDIMV%NFLEVL, &
 & NDGENL=>YDDIM%NDGENL, NDGSAL=>YDDIM%NDGSAL, &
 & NGPTOTG=>YDGEM%NGPTOTG, NGPTOT=>YDGEM%NGPTOT)
!     ------------------------------------------------------------------

LLP = NPRINTLEV >= 1.OR. LALLOPR
IU = NULOUT
ALLOCATE(YDMP%MYLEVS(NFLEVL))
IF(LLP)WRITE(IU,9) 'MYLEVS   ',SIZE(YDMP%MYLEVS   ),SHAPE(YDMP%MYLEVS   )
ALLOCATE(YDMP%MYLATS(NDGSAL:NDGENL))
IF(LLP)WRITE(IU,9) 'MYLATS   ',SIZE(YDMP%MYLATS   ),SHAPE(YDMP%MYLATS   )

ALLOCATE(YDMP%NGLOBALINDEX(NGPTOT))
IF(LLP)WRITE(IU,9) 'NGLOBALIND',SIZE(YDMP%NGLOBALINDEX),SHAPE(YDMP%NGLOBALINDEX)
ALLOCATE(YDMP%NGLOBALAT(NGPTOT))
IF(LLP)WRITE(IU,9) 'NGLOBALAT',SIZE(YDMP%NGLOBALAT),SHAPE(YDMP%NGLOBALAT)
ALLOCATE(YDMP%NGLOBALPROC(NGPTOTG))
IF(LLP)WRITE(IU,9) 'NGLOBALPRO',SIZE(YDMP%NGLOBALPROC),SHAPE(YDMP%NGLOBALPROC)
ALLOCATE(YDMP%NLOCALINDEX(NGPTOTG))
IF(LLP)WRITE(IU,9) 'NLOCALINDE',SIZE(YDMP%NLOCALINDEX),SHAPE(YDMP%NLOCALINDEX)

!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUALMP2',1,ZHOOK_HANDLE)
END SUBROUTINE SUALMP2
