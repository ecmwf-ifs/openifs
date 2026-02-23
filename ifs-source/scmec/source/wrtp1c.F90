! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE WRTP1C(YDGEOMETRY,YDMODEL,YDSURF)

!**** *WRTP1C*  - Write prognostic variables of the one-column model

!     Purpose.
!     --------
!     Write out prognostic variables

!**   Interface.
!     ----------
!        *CALL* *WRTP1C

!        Explicit arguments :
!        --------------------


!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------
!        ECMWF Research Department documentation 
!        of the single column model

!     Author.
!     -------
!        Joao Teixeira  *ECMWF*

!     Modifications.
!     --------------
!        Original      94-01-11
!        J.Teixeira   Jan.-1995  new output files.
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE TYPE_MODEL   , ONLY : MODEL
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE PARDIM1C
USE PARDIM
USE YOMPHYDS
USE YOMCT3   , ONLY : NSTEP
USE YOMGT1C0 , ONLY : UT0      ,VT0      ,TT0      ,QT0      ,&
                     &WT0      ,ST0      ,AT0      ,SPT0
USE YOMGP1C0 , ONLY : TSA0     ,WSA0     ,SNS0     ,TL0      ,WL0
USE YOMLOG1C , ONLY : NPOSASC
USE YOETHF
USE YOMCST
USE SURFACE_FIELDS_MIX, ONLY : TSURF

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(MODEL),    INTENT(INOUT) :: YDMODEL
TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF

INTEGER(KIND=JPIM) :: ICSS

INTEGER(KIND=JPIM) :: JSLEV,IST,IEND,JK,JALEV

REAL(KIND=JPRB)    :: ZSNNU0,ZRENU0,ZSCALE,ZQS
REAL(KIND=JPRB)    :: ZLINU0(YDSURF%YSP_SBD%NLEVS),ZPRS(0:YDGEOMETRY%YRDIMV%NFLEVG),ZPRSF(YDGEOMETRY%YRDIMV%NFLEVG), &
 & ZRH(YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB)    :: ZRDAW(YDSURF%YSP_SBD%NLEVS)       ! SOIL LEVEL THICKNESS

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

#include "fcttre.func.h"

!     ------------------------------------------------------------------
#include "surf_inq.h"
#include "gphpre.intfb.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WRTP1C',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM,YDMP=>YDGEOMETRY%YRMP, &
 & YDVAB=>YDGEOMETRY%YRVAB, YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,&
 & YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & YSURF=>YDEPHY%YSURF, &
 & YSP_SBD=>YDSURF%YSP_SBD)
ICSS=YSP_SBD%NLEVS

!*       1.   SCALING SOME VARIABLES AND CALCULATING PRESSURE.
!             ------------------------------------------------

ZSNNU0 = SNS0(1) / RHOH2O
ZRENU0 = WL0(1)  / RHOH2O
CALL SURF_INQ(YSURF,PRDAW=ZRDAW)
ZLINU0(1:ICSS) = WSA0(1:ICSS) / RHOH2O / ZRDAW(1:ICSS)

IST  = 1
IEND = 1

ZPRS(NFLEVG) = EXP(SPT0)

CALL GPHPRE(NPROMA,NFLEVG,IST,IEND,YDVAB,YDGEOMETRY%YRCVER,ZPRS,PRESF=ZPRSF)


!*       2.    WRITE TIME STEP.
!              ----------------


WRITE(NPOSASC,*) 'TIME STEP'
WRITE(NPOSASC,'(I6)') NSTEP

WRITE(NPOSASC,*) 'PROGNOSTIC VARIABLES'

!     -----------------------------------------------------------------
!       3.0  computation of rh.
!            ------------------

DO JK=1,NFLEVG
  ZQS=FOEEWM(TT0(JK))/ZPRSF(JK)
  ZQS=MIN(0.5_JPRB,ZQS)
  ZQS=ZQS/(1.0_JPRB-RETV*ZQS)
  ZRH(JK)=QT0(JK)/ZQS
ENDDO


!*       3.    WRITE ATMOSPHERIC VARIABLES.
!              ---------------------------


WRITE(NPOSASC,3019) '  P  ','  U  ','  V  ','  T  ','  Q  '&
               &,'  A  ','  L  ','  I  ','  RH  '
DO JALEV=1,NFLEVG
  WRITE(NPOSASC,3020) ZPRSF(JALEV),UT0(JALEV),VT0(JALEV)&
   &,TT0(JALEV),QT0(JALEV),AT0(JALEV),WT0(JALEV)&
   &,ST0(JALEV),ZRH(JALEV)
ENDDO

3019 FORMAT(9(2X,A12))
3020 FORMAT(9(2X,E12.6))

WRITE(NPOSASC,*) 'SURFACE PRESSURE  -  LOG SURFACE PRESSURE'
WRITE(NPOSASC,'(2(2x,E12.6))') ZPRS(NFLEVG),SPT0


!       3.1   WRITE IN UPPERAIR FILE.
!       -----------------------------


IF (NSTEP == 0) THEN
  WRITE(51,3119) 'LEVEL ',' STEP ','  P  ','  U  ','  V  '&
    &,'  T  ','  Q  ','  A  ','  L  ','  I  ', '  RH  '
ENDIF
DO JALEV=1,NFLEVG
  WRITE(51,3120) JALEV,NSTEP,ZPRSF(JALEV),UT0(JALEV),VT0(JALEV)&
    &,TT0(JALEV),QT0(JALEV),AT0(JALEV),WT0(JALEV)&
    &,ST0(JALEV),ZRH(JALEV)
ENDDO

3119 FORMAT(2(2X,A8),9(2X,A12))
3120 FORMAT(2(2X,I8),9(2X,E12.6))


!*       4.    WRITE SOIL VARIABLES.
!              ---------------------


WRITE(NPOSASC,*) ' SOIL TEMPERATURE - SOIL MOISTURE '
DO JSLEV=1,ICSS
  WRITE(NPOSASC,4020) TSA0(JSLEV) , ZLINU0(JSLEV)
ENDDO
4020 FORMAT(2(2X,E12.6))


!*       5.    WRITE SKIN VARIABLES.
!              ---------------------


WRITE(NPOSASC,*) ' SKIN TEMP. - SKIN. RES. CONT. '
WRITE(NPOSASC,5020) TL0 , ZRENU0
5020 FORMAT(2(2X,E12.6))


!*       6.    WRITE SNOW DEPTH.
!              -----------------


WRITE(NPOSASC,*) 'SNOW DEPTH'
WRITE(NPOSASC,'(E12.6)') ZSNNU0


!*       7.    WRITE APPARENT SURFACE HUMIDITY.
!              --------------------------------


! WRITE(NPOSASC,*) 'AP. SURF. HUMIDITY'
! WRITE(NPOSASC,'(E12.6)') VFASQ

!        7.1   WRITE IN SURFSOIL FILE.
!        -----------------------------


IF (ICSS == 4) THEN
  IF (NSTEP == 0) THEN
    WRITE(56,7119) 'STEP ',' Tskin ',' SRWC ',' SNOW '&
      &,' Tsoil(1) ',' Tsoil(2) ',' Tsoil(3) '&
      &,' Tsoil(4) ',' Qsoil(1) ',' Qsoil(2) ',' Qsoil(3) '&
      &,' Qsoil(4) '
  ENDIF
  WRITE(56,7120) NSTEP,TL0,ZRENU0,ZSNNU0 &
    &,TSA0(1),TSA0(2),TSA0(3),TSA0(4)&
    &,ZLINU0(1),ZLINU0(2),ZLINU0(3),ZLINU0(4) 

7119 FORMAT(1(2X,A8),11(2X,A12))
7120 FORMAT(1(2X,I8),11(2X,E12.6))
ENDIF

!     --------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('WRTP1C',1,ZHOOK_HANDLE)
END SUBROUTINE WRTP1C
