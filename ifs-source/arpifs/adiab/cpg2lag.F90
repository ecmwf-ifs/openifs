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

SUBROUTINE CPG2LAG(YDGEOMETRY,YDGMV,YDML_GCONF,YDML_DYN,KNUMB,YDSL,KIBL,&
 & PB1,PB2,PGMV,PGMVS,&
 & PUT1,PVT1,PSPT1)  

!**** *CPG2LAG* - Grid point calculations 2D lagged part .

!     Purpose.
!     --------
!           Grid point calculations 2D lagged part.

!**   Interface.
!     ----------
!        *CALL* *CPG2LAG(...)*

!        Explicit arguments :
!        --------------------

!        INPUT and INPUT-OUTPUT:
!         KNUMB     : number of elements of arrays for which computations
!                     are performed (used in message passing version)
!         YDSL      : SL_STRUCT definition.
!         KASLB1    : horizontal dimension for semi-Lagrangian buffer.
!         KIBL      : index into YRGSGEOM/YRCSGEOM instances in YDGEOMETRY
!         PB1       : buffer for quantities to be interpolated.
!         PB2       : buffer for quantities at the arrival grid-point.
!         PGMV      : "t-dt" and "t" upper air GMV variables.
!         PGMVS     : "t-dt" and "t" surface GMV variables.

!        OUTPUT:
!         PUT1      : t+dt U-wind.
!         PVT1      : t+dt V-wind.
!         PSPT1     : t+dt equivalent height.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE documentation vol 2 ch 1 and vol 3 ch 6

!     Author.
!     -------
!      Mats Hamrud  *ECMWF*
!      Original : 92-02-01

!     Modifications.
!     --------------
!      01-Oct-2003 M. Hamrud   CY28 Cleaning
!      30-Jul-2008 J. Masek    Dataflow for new SLHD interpolators.
!      G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      G. Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM and TCSGEOM
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
! End Modifications
!-------------------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD     , ONLY : MODEL_DYNAMICS_TYPE
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE YOMGMV                 , ONLY : TGMV
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT2                 , ONLY : NSTOP2
USE YOMCT3                 , ONLY : NSTEP
USE EINT_MOD               , ONLY : SL_STRUCT

!-------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)               ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)                   ,INTENT(INOUT) :: YDGMV
TYPE(MODEL_DYNAMICS_TYPE)    ,INTENT(INOUT) :: YDML_DYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
INTEGER(KIND=JPIM),INTENT(IN)    :: KNUMB 
TYPE(SL_STRUCT)   ,INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB1(YDSL%NASLB1,YDML_DYN%YRPTRSLB1%NFLDSLB1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,1,YDGMV%NDIMGMV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUT1(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVT1(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPT1(YDGEOMETRY%YRDIM%NPROMA) 
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: IEND, IFLDX, IST, JROF
REAL(KIND=JPRB)    :: ZDUMARR(1,1,1)
LOGICAL   :: LLLSTEP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gptf1.intfb.h"
#include "ladine.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CPG2LAG',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDPTRSLB2=>YDML_DYN%YRPTRSLB2, &
  & YDPTRSLB1=>YDML_DYN%YRPTRSLB1,YDDYNA=>YDML_DYN%YRDYNA)

ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, &
 & MSLB1SP0=>YDPTRSLB1%MSLB1SP0, MSLB1SP9=>YDPTRSLB1%MSLB1SP9, &
 & MSLB1U0=>YDPTRSLB1%MSLB1U0, MSLB1U9=>YDPTRSLB1%MSLB1U9, &
 & MSLB1UR0=>YDPTRSLB1%MSLB1UR0, MSLB1V0=>YDPTRSLB1%MSLB1V0, &
 & MSLB1V9=>YDPTRSLB1%MSLB1V9, MSLB1VR0=>YDPTRSLB1%MSLB1VR0, &
 & NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & MSLB2URL=>YDPTRSLB2%MSLB2URL, MSLB2USI=>YDPTRSLB2%MSLB2USI, &
 & MSLB2VRL=>YDPTRSLB2%MSLB2VRL, MSLB2VSI=>YDPTRSLB2%MSLB2VSI, &
 & NFLDSLB2=>YDPTRSLB2%NFLDSLB2)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY COMPUTATIONS.
!              -------------------------

IST =1
IEND=KNUMB
IF(NSTEP < NSTOP2) THEN
  LLLSTEP=.FALSE.
ELSE
  LLLSTEP=.TRUE.
ENDIF

!     ------------------------------------------------------------------

!*       2.   LAGGED PART OF GRIDPOINT COMPUTATIONS
!             -------------------------------------

!*       2.1  LAGGED SEMI-LAGRANGIAN CALCULATIONS

IF(YDDYNA%LSLAG) THEN

  IFLDX=NFLDSLB1
  CALL LADINE(YDGEOMETRY,YDML_GCONF%YRRIP,YDML_DYN%YRDYN,YDML_DYN%YRDYNA,IST,IEND,YDSL,YDML_DYN%YRSLINT,NPROMA,1,IFLDX,&
   & KIBL,&
   & PB1(1,MSLB1UR0),PB1(1,MSLB1VR0),&
   & PB2(1,MSLB2URL),PB2(1,MSLB2VRL),&
   & PB1(1,MSLB1U0),PB1(1,MSLB1V0),PB1(1,MSLB1SP0),&
   & PB1(1,MSLB1U9),PB1(1,MSLB1V9),PB1(1,MSLB1SP9),&
   & PB2(1,MSLB2USI),PB2(1,MSLB2VSI),&
   & PUT1,PVT1,PSPT1)

ENDIF

CALL GPTF1(YDGEOMETRY,YDGMV,YDML_GCONF,YDML_DYN%YRDYN,YDML_DYN%YRDYNA,.TRUE.,IST,IEND,LLLSTEP,PGMV,PGMVS,ZDUMARR)
 
!*       2.4  DIVIDE BY MAP FACTOR

DO JROF=IST,IEND
  PUT1(JROF)=PUT1(JROF)/YDGSGEOM%GM(JROF)
  PVT1(JROF)=PVT1(JROF)/YDGSGEOM%GM(JROF)
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG2LAG',1,ZHOOK_HANDLE)
END SUBROUTINE CPG2LAG
