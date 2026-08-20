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

SUBROUTINE SUVFE_OPER_SETUP( YDCVER,&
 & KTBC, KBBC, KFLEV_IN, &
 & KORDER, KBASIS, KINTERNALS, KOFF, KKNOTS )

!**** *SUVFE_OPER_SETUP*  - Set Up VFE - OPERATOR SETUP

!**   Interface.
!     ----------

!     *CALL* SUVFE_OPER_SETUP
!     Explicit arguments :
!     --------------------
!      * INPUT:
!        KTBC/KBBC    : KTBC/KBBC(1) = 1  - value of basis at top/bottom BC is 0,
!                       KTBC/KBBC(2) = N  - all derivatives up to N-order of basis
!                       at top/bottom are set to 0
!        KFLEV_IN     : size if onput vector data
!      * OUTPUT:
!        KORDER       : order of basis
!        KBASIS       : number of basis function need to satisfy all BCs and to
!                       represent input vector
!        KINTERNALS   : number of internals knots (multiple boundary knots are excluded)
!        KOFF         : offset of function due to TBC
!        KKNOTS       : total size of knot vector

!     Method.
!     -------
! This subroutine is supposed to be called inside suvertfeb.
! This routine is called for IN basis and for OUT basis independently.
!
! An example of order calculation for invertible operators when general order is NVFE_ORDER=4:
! 
! integral is required to be 0 at the surface => KBBC_OUT(1) = 1
! (we integrate from surface to model top). The derivated function must
! therefore fullfill the same condition due to consistency:
! KBBC_IN(1) = 1. The operator is defined for NFLEVG+1 levels. 
!
! The order of integrated function is:
! KBASIS = NFLEVG+1, KTBC_IN+KBBC_IN=0, fixed KINTERNALS = NFLEVG + 2 - 4
! order of in basis is  iorder_in = NFLEVG + 1 - NFLEVG - 2 + 4 = 3
!
! The order of output function is:
! KBASIS = NFLEVG+1+1 (KTBC_OUT+KBBC_OUT=1), fixed KINTERNALS = NFLEVG + 2 - 4
! order of out basis is  iorder_out = NFLEVG + 2 - NFLEVG - 2 + 4 = 4
!
! By default the VFE_SETUP is called for NFLEVG+2 levels with KORDER of splines.
! This is due to fact that in VFE full levels are considered NFLEVG internal levels
! + 2 material boundaries (eta=0 and eta=1). This reflects the field VFE_ETAF(0:NFLEVG+1)

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ALADIN-NH documentation.

!     Author.
!     -------
!        Jozef Vivoda, SHMU/LACE 
!        Original : 2017-09

!     Modifications.
!     --------------
!      P.Smolikova (Sep 2020): VFE pruning.
!     ------------------------------------------------------------------

USE PARKIND1, ONLY : JPRB, JPIM
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN  , ONLY : NULOUT
USE YOMCVER , ONLY : TCVER

IMPLICIT NONE

TYPE(TCVER),INTENT(IN)         :: YDCVER
INTEGER(KIND=JPIM), INTENT(IN) :: KTBC(2)
INTEGER(KIND=JPIM), INTENT(IN) :: KBBC(2)
INTEGER(KIND=JPIM), INTENT(IN) :: KFLEV_IN

INTEGER(KIND=JPIM), INTENT(OUT) :: KORDER
INTEGER(KIND=JPIM), INTENT(OUT) :: KBASIS
INTEGER(KIND=JPIM), INTENT(OUT) :: KINTERNALS
INTEGER(KIND=JPIM), INTENT(OUT) :: KOFF
INTEGER(KIND=JPIM), INTENT(OUT) :: KKNOTS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!------------------------------

! #include "abor1.intfb.h"

!------------------------------

IF (LHOOK) CALL DR_HOOK('SUVFE_OPER_SETUP',0,ZHOOK_HANDLE)

KOFF = KTBC(1) + KTBC(2)

KBASIS = KFLEV_IN + KOFF + KBBC(1) + KBBC(2)

KORDER     = YDCVER%NVFE_ORDER       ! order of splines consistent with namelist
KINTERNALS = KBASIS - KORDER
KKNOTS     = KBASIS + KORDER

WRITE(NULOUT,*) "(SUVFE_OPER_SETUP) NORDER : ", KORDER, " NBASIS : ", KBASIS, &
 & " NINTERNALS : ", KINTERNALS, " NOFF : ", KOFF, " NKNOTS : ", KKNOTS
WRITE(NULOUT,*) "(SUVFE_OPER_SETUP) KTBC   : ", KTBC(1), KTBC(2), " KBBC  :", KBBC(1), KBBC(2)

! check consistency
IF( KORDER < 2 )THEN
    CALL ABOR1("(SUVFE_OPER_SETUP) ORDER OF BASIS LESS THAN LINEAR. INCREASE NVFE_ORDER IN NAMDYNA.")
ENDIF

IF( (KTBC(1) + KTBC(2)) >= KORDER )THEN
    CALL ABOR1("(SUVFE_OPER_SETUP) ORDER OF BASIS NOT SUFFICIENT FOR REQUIRED TBC. INCREASE NVFE_ORDER IN NAMDYNA.")
ENDIF

IF( (KBBC(1) + KBBC(2)) >= KORDER )THEN
    CALL ABOR1("(SUVFE_OPER_SETUP) ORDER OF BASIS NOT SUFFICIENT FOR REQUIRED BBC. INCREASE NVFE_ORDER IN NAMDYNA.")
ENDIF

IF( KINTERNALS < 0 )THEN
    CALL ABOR1("(SUVFE_OPER_SETUP) TOO HIGH ORDER OF SPLINES REQUIRED. INCREASE NFLEVG OR DECREASE NVFE_ORDER IN NAMDYNA")
ENDIF

IF (LHOOK) CALL DR_HOOK('SUVFE_OPER_SETUP',1,ZHOOK_HANDLE)
END SUBROUTINE SUVFE_OPER_SETUP 
  
