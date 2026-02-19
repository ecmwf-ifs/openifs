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

SUBROUTINE SUALNUD(YDGEOMETRY,YGFL)

!**** *SUALNUD * - Routine to allocate space for nudging variables

!     Purpose.
!     --------
!           Allocate space for nudging variables.

!**   Interface.
!     ----------
!        *CALL* *SUALNUD*

!     Explicit arguments :  None
!     --------------------
!        Called by SU0YOMB.

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
!      Tomas Wilhelmsson *ECMWF*
!      Original : 19-12-2012

!     Modifications.
!     --------------
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : NPRINTLEV
USE YOMCT0   , ONLY : LALLOPR
USE YOMLUN   , ONLY : NULOUT
USE YOMSNU   , ONLY : XPNUDG, TNUDTE, TNUDSH, TNUDDI, TNUDVO, TNUDSV, TNUDLP  
USE YOMNUD   , ONLY : NFNUDG, LNUDG, LNUDDI, LNUDLP, LNUDSH, LNUDSV, LNUDTE, LNUDVO  
USE YOM_YGFL , ONLY : TYPE_GFLD

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(TYPE_GFLD),INTENT(INOUT) :: YGFL
INTEGER(KIND=JPIM) :: IU

LOGICAL ::  LLP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUALNUD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NGFL_EXT=>YGFL%NGFL_EXT, &
 & NSPEC2G=>YDDIM%NSPEC2G, &
 & NFLEVL=>YDDIMV%NFLEVL)
!     ------------------------------------------------------------------

!*       1.    ALLOCATE SPACE FOR ARRAYS.
!              --------------------------

LLP = NPRINTLEV >= 1.OR. LALLOPR
IU = NULOUT
IF (LLP) WRITE(NULOUT,'('' SUALNUD PRINTOUTS '')')

!*       1.1 NUDGING (SPECTRAL FIELDS + TIME WEIGHTS)

IF(LNUDG) THEN
  ALLOCATE(XPNUDG(NFNUDG))
  IF(LLP)WRITE(IU,9) 'XPNUDG   ',SIZE(XPNUDG),SHAPE(XPNUDG)
  IF(LNUDTE) THEN
    ALLOCATE(TNUDTE(NSPEC2G,NFLEVL,NFNUDG))
    IF(LLP)WRITE(IU,9) 'TNUDTE   ',SIZE(TNUDTE),SHAPE(TNUDTE)
  ENDIF
  IF(LNUDSH) THEN
    ALLOCATE(TNUDSH(NSPEC2G,NFLEVL,NFNUDG))
    IF(LLP)WRITE(IU,9) 'TNUDSH   ',SIZE(TNUDSH),SHAPE(TNUDSH)
  ENDIF
  IF(LNUDDI) THEN
    ALLOCATE(TNUDDI(NSPEC2G,NFLEVL,NFNUDG))
    IF(LLP)WRITE(IU,9) 'TNUDDI   ',SIZE(TNUDDI),SHAPE(TNUDDI)
  ENDIF
  IF(LNUDVO) THEN
    ALLOCATE(TNUDVO(NSPEC2G,NFLEVL,NFNUDG))
    IF(LLP)WRITE(IU,9) 'TNUDVO   ',SIZE(TNUDVO),SHAPE(TNUDVO)
  ENDIF
  ! === obsolescent =============================================
  IF(LNUDSV) THEN
    ALLOCATE(TNUDSV(NSPEC2G,NFLEVL,NGFL_EXT,NFNUDG))
    IF(LLP)WRITE(IU,9) 'TNUDSV   ',SIZE(TNUDSV),SHAPE(TNUDSV)
  ENDIF
  ! (TNUDSH+TNUDSV) have to be replaced by TNUDGFL in the future.
  ! =============================================================
  IF(LNUDLP) THEN
    ALLOCATE(TNUDLP(NSPEC2G,NFNUDG))
    IF(LLP)WRITE(IU,9) 'TNUDLP   ',SIZE(TNUDLP),SHAPE(TNUDLP)
  ENDIF
ENDIF

!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUALNUD',1,ZHOOK_HANDLE)
END SUBROUTINE SUALNUD
