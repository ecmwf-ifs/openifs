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

SUBROUTINE GPHLUV(YDDIMV,KPROMA,KSTART,KPROF,PU,PV,PUVH)

!**** *GPHLUV* - wind components calculation in half-levels

!     Purpose.
!     --------
!           Compute wind components in half-levels

!**   Interface.
!     ----------
!        *CALL* *GPHLUV()

!        Explicit arguments :
!        --------------------
!        * INPUT:
!          KPROMA  - length of work
!          KSTART  - start of work
!          KPROF   - end of work
!          PU      - U-wind at full levels
!          PV      - V-wind at full levels

!        * IN/OUT:
!          PUVH    - horizontal wind and weights at half levels
!                    (IN) for weights, (OUT) for half level wind

!        Implicit arguments : none.
!        --------------------

!     Method.
!     -------
!        Interpolation from full-levels using logaritmic pressure profile
!        Then modify values on the top and bottom using boundary condition
!        (at the bottom free-slip condition)
!        Store also the weights of vertical interpolation

!     Externals.
!     ----------
!     Reference.
!     ----------
!        ARPEGE/ALADIN documentation

!     Author.
!     -------
!      Radmila Bubnova & Martin Janousek,  CNRM/GMAP/EXT
!      Original : February 1996

! Modifications
! -------------
!   Modified 02-07-02 by C. Fischer : intents for dummy arrays
!   Modified 08-2002 C. Smith and K. YESSAD:
!    - make optional the vertical averaging.
!    - some cleanings and optimisations.
!   01-Oct-2003 M. Hamrud  CY28 Cleaning
!   09-Jun-2004 J. Masek   NH cleaning (LVSLWBC)
!   21-Jan-2005 K. Yessad  Remove useless dummy arguments.
!   07-Mar-2007 K. Yessad  Remove LVERAVE_HLUV.
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
! End Modifications
!------------------------------------------------------------------

USE YOMDIMV   , ONLY : TDIMV
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE INTDYN_MOD, ONLY : YYTHW

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIMV),       INTENT(IN)    :: YDDIMV
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KPROMA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KPROMA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUVH(KPROMA,0:YDDIMV%NFLEVG,YYTHW%NDIM) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF, JLEV
REAL(KIND=JPRB) :: ZWW

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPHLUV',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG)

!     ------------------------------------------------------------------

!*    1. General case:

DO JLEV=1, NFLEVG - 1
  DO JROF=KSTART,KPROF
    ZWW=PUVH(JROF,JLEV,YYTHW%M_WWI)
    PUVH(JROF,JLEV,YYTHW%M_UH)=ZWW*PU(JROF,JLEV)+(1.0_JPRB-ZWW)*PU(JROF,JLEV+1)  
    PUVH(JROF,JLEV,YYTHW%M_VH)=ZWW*PV(JROF,JLEV)+(1.0_JPRB-ZWW)*PV(JROF,JLEV+1)  
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*    2. Top.

PUVH(KSTART:KPROF,0,YYTHW%M_UH) = PU(KSTART:KPROF,1)
PUVH(KSTART:KPROF,0,YYTHW%M_VH) = PV(KSTART:KPROF,1)

!     ------------------------------------------------------------------

!*    3. Surface.

PUVH(KSTART:KPROF,NFLEVG,YYTHW%M_UH) = PU(KSTART:KPROF,NFLEVG)
PUVH(KSTART:KPROF,NFLEVG,YYTHW%M_VH) = PV(KSTART:KPROF,NFLEVG)

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPHLUV',1,ZHOOK_HANDLE)
END SUBROUTINE GPHLUV
