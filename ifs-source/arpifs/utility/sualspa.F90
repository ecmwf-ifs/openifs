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

SUBROUTINE SUALSPA(YDGEOMETRY,LDGRADSP)

!**** *SUALSPA * - Routine to allocate space for SPA3

!     Purpose.
!     --------
!           Allocate space for SPA3

!**   Interface.
!     ----------
!        *CALL* *SUALSPA*

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
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-26

!     Modifications.
!     --------------
!      R. El Khatib  01-08-20 make Ozone pointer point somewhere in any case for safety
!      P.Smolikova 02-09-30 SPA3AUX allocation for d4 in NH
!      M.Hamrud   03-08-01 GFL introduction
!      C.Fischer  03-09-30 Aladin mean wind profiles spa1
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      G.Radnoti  04-03-19 nflevg becomes nflevl for Aladin B-level
!      R.Brozkova 05-07-17 "X" part of div.: SPNHX to NFTHER (prev. SPVDAUX)
!      R. El Khatib  04-Jul-2014 LDALLOC_SPEC : temporary optional key to have the input
!       structure allocated with alloc_spec, making the structure fully
!       compatible with spectral_fields_mod tools. Default is .TRUE. if not present.
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0       , ONLY : LALLOPR
USE YOMMP0       , ONLY : NPRINTLEV
USE YOMSP        , ONLY : SPVOR_FLT, SPDIV_FLT, SPUB_FLT, SPVB_FLT
USE YOMLUN       , ONLY : NULOUT
USE DEALLOCATE_IF_ASSOCIATED_MOD, ONLY : DEALLOCATE_IF_ASSOCIATED

!     ------------------------------------------------------------------

IMPLICIT NONE
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
LOGICAL,        INTENT(IN) :: LDGRADSP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUALSPA',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NSPEC2=>YDDIM%NSPEC2, &
 & NFLSUR=>YDDIMV%NFLSUR, NFLEVL=>YDDIMV%NFLEVL)


!     ------------------------------------------------------------------

CALL DEALLOCATE_IF_ASSOCIATED(SPVOR_FLT)
CALL DEALLOCATE_IF_ASSOCIATED(SPDIV_FLT)

IF( LDGRADSP ) THEN
  ALLOCATE(SPVOR_FLT(NFLSUR,NSPEC2))
  SPVOR_FLT(1:NFLSUR,1:NSPEC2)=0._JPRB
  ALLOCATE(SPDIV_FLT(NFLSUR,NSPEC2))
  SPDIV_FLT(1:NFLSUR,1:NSPEC2)=0._JPRB
  ALLOCATE(SPUB_FLT(NFLEVL))
  SPUB_FLT(1:NFLEVL)=0._JPRB
  ALLOCATE(SPVB_FLT(NFLEVL))
  SPVB_FLT(1:NFLEVL)=0._JPRB
  IF(NPRINTLEV >= 1.OR. LALLOPR) THEN
    WRITE(NULOUT,9990) 'SPVOR_FLT     ',SIZE(SPVOR_FLT     ),SHAPE(SPVOR_FLT  )
    WRITE(NULOUT,9990) 'SPDIV_FLT     ',SIZE(SPDIV_FLT     ),SHAPE(SPDIV_FLT  )
  ENDIF
ENDIF

9990 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUALSPA',1,ZHOOK_HANDLE)

END SUBROUTINE SUALSPA
