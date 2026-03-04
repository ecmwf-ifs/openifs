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

SUBROUTINE CAPESHEAR(YDVAB,KPROMA,KST,KND,KFLEVG,&
 & PRESH,PRESF,PU,PV,PCAPE,&
 & PCAPES)

!**** *POS* - Calculate Cape-Shear

!     PURPOSE.
!     -------

!**   INTERFACE.
!     ---------
!        *CALL* *CAPESHEAR(...)

!        EXPLICIT ARGUMENTS :
!        ------------------
!        * INPUT:
!        YDVAB     : vertical geometry
!        KPROMA    : horizontal dimension.
!        KST       : first point to post-process 
!        KND       : last point to post-process 
!        KFLEVG    : number of model vertical levels (NFLEVG)
!        PRESH     : half level pressure
!        PRESF     : full level pressure
!        PU        : zonal wind
!        PV        : meridional wind
!        PCAPE     : CAPE

!        * OUTPUT (post-processed fields):
!        PCAPES    : CAPE-SHEAR

!     EXTERNALS.
!     ---------

! Author
! ------
!   Tomas Wilhelmsson *ECMWF* : 2017-02-16

! Modifications
! -------------
!------------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMVERT  , ONLY : TVAB
USE PARDIM   , ONLY : JPNPPM
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVAB)        ,INTENT(IN)    :: YDVAB
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEVG
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESH(KPROMA,0:KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESF(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KPROMA,KFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KPROMA,KFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCAPE(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAPES(KPROMA)
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPOPLEV=2, JPLEV500=1, JPLEV925=2

!!REAL(KIND=JPRB) :: ZUT0(KPROMA,KFLEVG),ZVT0(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZPRES(KPROMA,JPOPLEV),ZLNPRES(KPROMA,JPOPLEV)
REAL(KIND=JPRB) :: ZPXPD(KPROMA,0:KFLEVG,JPNPPM)
REAL(KIND=JPRB) :: ZPXP(KPROMA,0:KFLEVG,JPNPPM)
REAL(KIND=JPRB) :: ZUDIFF(KPROMA),ZVDIFF(KPROMA)
REAL(KIND=JPRB) :: ZPPU(KPROMA,JPOPLEV),ZPPV(KPROMA,JPOPLEV)
REAL(KIND=JPRB) :: ZPXLEV(JPOPLEV)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: ILOLEV
INTEGER(KIND=JPIM) :: JLEV, JLEVP, JROF
INTEGER(KIND=JPIM) :: ILEVB(KPROMA,JPOPLEV,JPNPPM)

LOGICAL :: LLBELO(KPROMA,JPOPLEV), LLBLOW(JPOPLEV)
LOGICAL :: LLBELS(KPROMA,JPOPLEV), LLBLES(JPOPLEV)

!     ------------------------------------------------------------------

#include "ppflev.intfb.h"
#include "ppinit.intfb.h"
#include "ppuv.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CAPESHEAR',0,ZHOOK_HANDLE)
!*       1.    INITIALIZATIONS.
!              ---------------

ILOLEV=1

!*    1.2  ARRAYS INITIALIZATIONS

!*      1.2.3   Z[X]T0 for GMV and GMVS

!! Ignoring map factor, so this will be wrong on a stretched grid.
!!
!!DO JLEV=1,KFLEVG
!!  DO JROF=KST,KND
!!    ZUT0(JROF,JLEV)=PU(JROF,JLEV)*PGM(JROF)
!!    ZVT0(JROF,JLEV)=PV(JROF,JLEV)*PGM(JROF)
!!  ENDDO
!!ENDDO

!     ------------------------------------------------------------------
!*       2.    PREPARE FOR INTERPOLATIONS.
!              ---------------------------

!*      2.1      COMPUTE POST-PROCESSING LEVELS

ZPXLEV(JPLEV500) = 50000._JPRB ! 500 hPa
ZPXLEV(JPLEV925) = 92500._JPRB ! 925 hPa

DO JLEVP=1,JPOPLEV
  DO JROF=KST,KND
    ZPRES(JROF,JLEVP)=ZPXLEV(JLEVP)
    ZLNPRES(JROF,JLEVP)=LOG(ZPXLEV(JLEVP))
  ENDDO
ENDDO

!*      2.2     SET UP HELP-ARRAYS FOR VERTICAL INTERPOLATION

CALL PPINIT(KPROMA,KST,KND,KFLEVG,JPNPPM,PRESH,PRESF,ZPXP,ZPXPD)

!*      2.3     FIND MODEL LEVEL UNDER SPECIFIED OUTPUT LEVELS

CALL PPFLEV(KPROMA,KST,KND,KFLEVG,JPOPLEV,JPNPPM,ZPRES,PRESH,&
  & PRESF,ILEVB,LLBELO,LLBELS,LLBLOW,LLBLES)  

!     ------------------------------------------------------------------
!*       3.    POST-PROCESSING (INTERPOLATIONS).
!              ---------------------------------

CALL PPUV(KPROMA,KST,KND,KFLEVG,JPOPLEV,ILOLEV,JPNPPM,ILEVB,LLBELO,LLBLOW,&
  & ZLNPRES,ZPXP,ZPXPD,PU,PV,ZPPU,ZPPV)  

DO JROF=KST,KND
  ZUDIFF(JROF)=ZPPU(JROF,JPLEV500)-ZPPU(JROF,JPLEV925)
  ZVDIFF(JROF)=ZPPV(JROF,JPLEV500)-ZPPV(JROF,JPLEV925)
  PCAPES(JROF)=SQRT(PCAPE(JROF)*(ZUDIFF(JROF)*ZUDIFF(JROF)+ZVDIFF(JROF)*ZVDIFF(JROF)))
ENDDO

IF (LHOOK) CALL DR_HOOK('CAPESHEAR',1,ZHOOK_HANDLE)

END SUBROUTINE CAPESHEAR

