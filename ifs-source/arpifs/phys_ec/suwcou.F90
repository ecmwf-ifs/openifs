! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUWCOU(YDRIP,YDEWCOU)

!**** *SETWCOU* - PRESET CONSTANTS *COMWCOU*.

!     PURPOSE.
!     --------

!           PRESET CONSTANTS IN *YOEWCOU*.

!     INTERFACE.
!     ----------

!           *SETWCOU* IS CALLED FROM *INICOM*.

!     EXTERNALS
!     ---------

!           NONE.

!     AUTHOR.
!     -------
!      P. VITERBO      E.C.M.W.F.      07/10/88

!     MODIFICATIONS.
!     --------------
!      Modified 05/01/14 G. Mozdzynski Optimise coupling communications
!      D.Salmond     22-Nov-2005 Mods for coarser/finer physics
!      Modified 06/08/16 J. Bidlot CBEGDAT is now 14 character long 
!      MODIFIED 08/02/01 J. Bidlot introduce WVWAMINIT
!      K. Yessad (July 2014): Move some variables.
! -----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOEWCOU  , ONLY : TEWCOU
USE YOMRIP   , ONLY : TRIP
USE YOMRIP0  , ONLY : NINDAT, NSSSSS

! -----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TEWCOU),INTENT(INOUT):: YDEWCOU
TYPE(TRIP)  ,INTENT(INOUT):: YDRIP
INTEGER(KIND=JPIM) :: IHH, IMN, ISS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------
#include "fcttim.func.h"
! -----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUWCOU',0,ZHOOK_HANDLE)
ASSOCIATE(CBEGDAT=>YDEWCOU%CBEGDAT, &
 & NGRIB_HANDLE_FOR_WAM=>YDEWCOU%NGRIB_HANDLE_FOR_WAM, NLATW=>YDEWCOU%NLATW, &
 & NLONW=>YDEWCOU%NLONW, RMISSW=>YDEWCOU%RMISSW, &
 & NSTOP=>YDRIP%NSTOP, TSTEP=>YDRIP%TSTEP)
! -----------------------------------------------------------------------

RMISSW = -999._JPRB

NLATW=1
NLONW=1

NGRIB_HANDLE_FOR_WAM=-99

!     DEFINE INITIAL DATE

IHH=NCTH(NSSSSS)
IMN=REAL((NSSSSS-IHH*3600),JPRB)/60._JPRB
ISS=NSSSSS-IHH*3600-IMN*60
WRITE(YDEWCOU%CBEGDAT(1:8),'(I8)') NINDAT
WRITE(YDEWCOU%CBEGDAT(9:14),'(I6.6)') IHH*10000+IMN*100+ISS

! -----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUWCOU',1,ZHOOK_HANDLE)
END SUBROUTINE SUWCOU
