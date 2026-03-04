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

SUBROUTINE CPG2(YDGEOMETRY,YDGMV,YDML_GCONF,YDML_DYN,CDCONF,KNUMB,KSTGLO,PB1,PB2,PGMV,PGMVS,&
 & KIBL,&
 & PUT0,PVT0,PUT0L,PVT0L,PDIVT0,PVORT0,&
 & PSPT0,PSPT0L,PSPT0M,&
 & PUT9,PVT9,PDIVT9,PSPNLT9,&
 & PSPT9,PSPT9L,PSPT9M,&
 & PUT1,PVT1,PSPT1,PTRAJ_SLAG)

!**** *CPG2* - Grid point calculations - 2-D model.

!     Purpose.
!     --------
!           Grid point calculations in dynamics - 2-D MODEL.

!**   Interface.
!     ----------
!        *CALL* *CPG2(...)

!        Explicit arguments :
!        --------------------

!        CDCONF    : configuration of work                       (I)
!        KNUMB     : number of elements of arrays for which
!                    computations are performed (MP version)     (I)
!        KSTGLO    : global offset (MP version)                  (I)
!        PB1       : SL buffer for quantities to be interpolated (O)
!        PB2       : Buffer for t+dt quantities                  (O)
!        KIBL      : index into YRGSGEOM/YRCSGEOM instances in YDGEOMETRY (I)
!        PUT0      : zonal wind time t                           (I/O)
!        PVT0      : meridian wind time t                        (I/O)
!        PUT0L     : zonal derivative of zonal wind time t       (I/O)
!        PVT0L     : zonal derivative of meridian wind time t    (I/O)
!        PDIVT0    : divergence  time t                          (I/O)
!        PVORT0    : vorticity  time t                           (I/O)
!        PSPT0     : equivalent height time t                    (I)
!        PSPT0L    : zonal derivative of PSPT0                   (I/O)
!        PSPT0M    : meridian derivative of PSPT0                (I/O)
!        PTRAJ_SLAG : stored trajectory for TL/AD models         (O)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Calls LACDYNSHD if SL scheme.
!        Calls WRSLTRAJ if writing semi-Lagrangian trajectory.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 88-02-04

! Modifications
! -------------
!   Modified 01-07-11 by K. YESSAD: relaxation of thin layer hyp for SL2TL.
!   Modified 01-08-30 by K. YESSAD: pruning and some other cleanings.
!   Modified 03-2002  J.Vivoda iterface to SC2...9 routines
!   Modified 02-09-30 by P. Smolikova: iterface to SC2...9 for d4 in NH
!   01-Oct-2003 M. Hamrud  CY28 Cleaning
!   Modified 04-01-22 by C.Temperton: Bugfixes (use of LUVDER,LVOR,PSPNLT9)
!   09-Jun-2004 J. Masek   NH cleaning (LLFULLIMP)
!   K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   F. Vana  28-Nov-2013 : Redesigned trajectory handling.
!   K. Yessad (July 2014): Move some variables.
! End Modifications
!     ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD     , ONLY : MODEL_DYNAMICS_TYPE
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE YOMGMV                 , ONLY : TGMV
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0                 , ONLY : LRSHW ,LRVEQ, NCONF
USE YOMCT3                 , ONLY : NSTEP
USE YOMVWRK                , ONLY : NTRSLTYPE
USE YOMTRAJ                , ONLY : TRAJ_SLAG_TYPE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)          ,INTENT(INOUT) :: YDGMV
TYPE(MODEL_DYNAMICS_TYPE),INTENT(INOUT):: YDML_DYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
CHARACTER(LEN=1)    ,INTENT(IN)    :: CDCONF
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KNUMB
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KSTGLO
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PB1(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,1,YDGMV%NDIMGMV)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS)
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KIBL
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PUT0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PVT0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PUT0L(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PVT0L(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PDIVT0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PVORT0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(IN)    :: PSPT0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PSPT0L(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PSPT0M(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PUT9(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PVT9(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PDIVT9(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PSPNLT9(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(IN)    :: PSPT9(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PSPT9L(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PSPT9M(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PUT1(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PVT1(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PSPT1(YDGEOMETRY%YRDIM%NPROMA)
TYPE(TRAJ_SLAG_TYPE),INTENT(INOUT) :: PTRAJ_SLAG
!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZBDT(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZDUM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZOROGL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZOROGM(YDGEOMETRY%YRDIM%NPROMA)
INTEGER(KIND=JPIM) :: IEND, IST, JROF
LOGICAL :: LLFSTEP
REAL(KIND=JPRB) :: ZBT, ZDT, ZTE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "gptf2.intfb.h"
#include "lacdynshw.intfb.h"




!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPG2',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL), &
 & YDOROG=>YDGEOMETRY%YROROG(KIBL), YDDYN=>YDML_DYN%YRDYN,YDDYNA=>YDML_DYN%YRDYNA, &
 & YDPTRSLB1=>YDML_DYN%YRPTRSLB1, &
 & YDPTRSLB2=>YDML_DYN%YRPTRSLB2,YDRIP=>YDML_GCONF%YRRIP, &
 & YDDIMF=>YDML_GCONF%YRDIMF)

ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & LUVDER=>YDDIMF%LUVDER, LVOR=>YDDIMF%LVOR, &
 & BETADT=>YDDYN%BETADT, LIMPF=>YDDYN%LIMPF, LSIDG=>YDDYN%LSIDG, &
 & SIVP=>YDDYN%SIVP, &
 & RSTRET=>YDGEM%RSTRET, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, &
 & MSLB1SP0=>YDPTRSLB1%MSLB1SP0, MSLB1SP9=>YDPTRSLB1%MSLB1SP9, &
 & MSLB1U0=>YDPTRSLB1%MSLB1U0, MSLB1U9=>YDPTRSLB1%MSLB1U9, &
 & MSLB1UR0=>YDPTRSLB1%MSLB1UR0, MSLB1V0=>YDPTRSLB1%MSLB1V0, &
 & MSLB1V9=>YDPTRSLB1%MSLB1V9, MSLB1VR0=>YDPTRSLB1%MSLB1VR0, &
 & MSLBUF1=>YDPTRSLB1%MSLBUF1, NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & MSLB2URL=>YDPTRSLB2%MSLB2URL, MSLB2USI=>YDPTRSLB2%MSLB2USI, &
 & MSLB2VRL=>YDPTRSLB2%MSLB2VRL, MSLB2VSI=>YDPTRSLB2%MSLB2VSI, &
 & MSLBUF2=>YDPTRSLB2%MSLBUF2, NFLDSLB2=>YDPTRSLB2%NFLDSLB2, &
 & TDT=>YDRIP%TDT)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

IF    (CDCONF == 'A')THEN
  ZDT=TDT
  ZTE=1.0_JPRB
ELSEIF(CDCONF == 'B')THEN
  ZDT=0.0_JPRB
  ZTE=1.0_JPRB
ELSE
  CALL ABOR1(' CPG2 - UNKNOWN CONFIGURATION ')
ENDIF

IF(NSTEP > 0) THEN
  LLFSTEP=.FALSE.
ELSE
  LLFSTEP=.TRUE.
ENDIF

IST =1
IEND=KNUMB

!     ------------------------------------------------------------------

!*       2.    INTERFACE TO GLOBAL ARRAYS.
!              ---------------------------

!*       2.2   TIME FILTER (PART 2).

CALL GPTF2(YDGEOMETRY,YDGMV,&
 ! --- INPUT ---------------------------------------------------------
 & YDML_GCONF,YDDYN,YDDYNA,IST,IEND,LLFSTEP,&
 ! --- INPUT for P.T0.; INPUT/OUTPUT for P.T9. -----------------------
 & PGMV,PGMVS,ZDUM)

!*       2.3   MAP FACTOR

DO JROF=IST,IEND
  PUT0  (JROF)=PUT0  (JROF)* YDGSGEOM%GM(JROF)
  PVT0  (JROF)=PVT0  (JROF)* YDGSGEOM%GM(JROF)
  IF (LUVDER) THEN
    PUT0L (JROF)=PUT0L (JROF)*(YDGSGEOM%GM(JROF)**2)
    PVT0L (JROF)=PVT0L (JROF)*(YDGSGEOM%GM(JROF)**2)
  ENDIF
  PDIVT0(JROF)=PDIVT0(JROF)*(YDGSGEOM%GM(JROF)**2)
  IF (LVOR) THEN
    PVORT0(JROF)=PVORT0(JROF)*(YDGSGEOM%GM(JROF)**2)
  ENDIF
  PSPT0L(JROF)=PSPT0L(JROF)* YDGSGEOM%GM(JROF)
  PSPT0M(JROF)=PSPT0M(JROF)* YDGSGEOM%GM(JROF)
  IF (.NOT.YDDYNA%LELTRA) THEN
    PUT9  (JROF)=PUT9  (JROF)* YDGSGEOM%GM(JROF)
    PVT9  (JROF)=PVT9  (JROF)* YDGSGEOM%GM(JROF)
  ENDIF
  IF (.NOT.YDDYNA%LTWOTL) THEN
    PDIVT9(JROF)=PDIVT9(JROF)*(YDGSGEOM%GM(JROF)**2)
    PSPT9L(JROF)=PSPT9L(JROF)* YDGSGEOM%GM(JROF)
    PSPT9M(JROF)=PSPT9M(JROF)* YDGSGEOM%GM(JROF)
  ENDIF
  ZOROGL(JROF)=YDOROG%OROGL(JROF)* YDGSGEOM%GM(JROF)
  ZOROGM(JROF)=YDOROG%OROGM(JROF)* YDGSGEOM%GM(JROF)
ENDDO

!     ------------------------------------------------------------------

!*       3.    SHALLOW WATER MODEL.
!              --------------------

IF(LRSHW)THEN

  IF(.NOT.YDDYNA%LSLAG) THEN

!*        3.1  EULERIAN SCHEME.

!*        3.1.1   EXPLICIT STEP FOR EULERIAN SCHEME.

    DO JROF=IST,IEND
      PUT1(JROF)=PUT9(JROF)*ZTE+ZDT*(&
       & PVORT0(JROF)*PVT0(JROF)&
       & - PVT0L(JROF)*PVT0(JROF)&
       & - PUT0L(JROF)*PUT0(JROF)&
       & +(PUT0 (JROF)*PUT0(JROF)+PVT0(JROF)*PVT0(JROF))*0.5_JPRB&
       & *YDCSGEOM%RATATX(JROF)&
       & +YDGSGEOM%RCORI(JROF)*PVT0(JROF)&
       & -PSPT0L(JROF)               )
      PVT1(JROF)=PVT9(JROF)*ZTE+ZDT*(&
       & -PDIVT0(JROF)*PVT0(JROF)&
       & - PVT0L(JROF)*PUT0(JROF)&
       & + PUT0L(JROF)*PVT0(JROF)&
       & -(PUT0 (JROF)*PUT0(JROF)+PVT0(JROF)*PVT0(JROF))*0.5_JPRB&
       & *YDCSGEOM%RATATH(JROF)&
       & -YDGSGEOM%RCORI(JROF)*PUT0(JROF)&
       & -PSPT0M(JROF)               )
      PSPT1(JROF)=PSPT9(JROF)*ZTE-ZDT*(&
       & PUT0(JROF)*(PSPT0L(JROF)-ZOROGL(JROF))&
       & + PVT0(JROF)*(PSPT0M(JROF)-ZOROGM(JROF))&
       & + PDIVT0(JROF)*(PSPT0(JROF)-YDOROG%OROG(JROF)) )
    ENDDO

!*        3.1.2   SEMI-IMPLICIT ADJUSTMENT FOR EULERIAN SCHEME.

    ZBT=BETADT*0.5_JPRB*ZDT*ZTE
    DO JROF=IST,IEND
      IF (LSIDG) THEN
        ZBDT(JROF)=ZBT*SIVP(1)
      ELSE
        ZBDT(JROF)=ZBT*SIVP(1)*RSTRET*RSTRET/(YDGSGEOM%GM(JROF)*YDGSGEOM%GM(JROF))
      ENDIF
    ENDDO

    DO JROF=IST,IEND
      PSPT1(JROF)=PSPT1(JROF)&
       & +ZBDT(JROF)*(2.0_JPRB*PDIVT0(JROF)-PDIVT9(JROF))
      PUT1(JROF)=PUT1(JROF)&
       & +ZBT*(2.0_JPRB*PSPT0L(JROF)-PSPT9L(JROF))
      PVT1(JROF)=PVT1(JROF)&
       & +ZBT*(2.0_JPRB*PSPT0M(JROF)-PSPT9M(JROF))
      IF (LIMPF) THEN
        PUT1(JROF)=PUT1(JROF)&
         & -ZBT*YDGSGEOM%RCORI(JROF)*(2.0_JPRB*PVT0(JROF)-PVT9(JROF))
        PVT1(JROF)=PVT1(JROF)&
         & +ZBT*YDGSGEOM%RCORI(JROF)*(2.0_JPRB*PUT0(JROF)-PUT9(JROF))
      ENDIF
    ENDDO

  ELSE

!*        3.2  SEMI-LAGRANGIAN SCHEME.

!*          3.2.1   WRITE t-dt VARIABLES.

    IF((NCONF == 421.OR.NCONF == 521).AND.(NTRSLTYPE == 1)) THEN
      CALL WRSLTRAJ2(YDDIM,IST,IEND,PTRAJ_SLAG,&
       & PUT9,PVT9,PSPNLT9,ZDUM,ZDUM,&
       & ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM)
    ENDIF

!         * Preliminary fill buffers with zero.

    PB1(IST:IEND,1:NFLDSLB1)=0.0_JPRB
    PB2(IST:IEND,1:NFLDSLB2)=0.0_JPRB

!         * Fill buffers.

    CALL LACDYNSHW(YDGEM,YDDYN,YDDYNA,CDCONF,NPROMA,IST,IEND,ZDT,&
     & PUT9,PVT9,PSPT9,PDIVT9,PSPT9L,PSPT9M,&
     & PUT0,PVT0,PSPT0,PDIVT0,PSPT0L,PSPT0M,&
     & YDGSGEOM,YDOROG%OROG,ZOROGL,ZOROGM,&
     & PSPNLT9,&
     & PB1(1,MSLB1U9),PB1(1,MSLB1V9),PB1(1,MSLB1SP9),&
     & PB1(1,MSLB1U0),PB1(1,MSLB1V0),PB1(1,MSLB1SP0),&
     & PUT1,PVT1,PSPT1,&
     & PB1(1,MSLB1UR0),PB1(1,MSLB1VR0),&
     & PB2(1,MSLB2USI),PB2(1,MSLB2VSI),&
     & PB2(1,MSLB2URL),PB2(1,MSLB2VRL))

!*          3.2.2   WRITE SL BUFFERS (saving trajectory)

    IF((NCONF == 421.OR.NCONF == 521).AND.(NTRSLTYPE == 2)) THEN
      CALL WRSLTRAJ2(YDDIM,IST,IEND,PTRAJ_SLAG,ZDUM,ZDUM,ZDUM,&
       & PB1(1,MSLB1U9),PB1(1,MSLB1V9),&
       & PB1(1,MSLB1SP9),PB1(1,MSLB1U0),&
       & PB1(1,MSLB1V0),PB1(1,MSLB1SP0),&
       & PB1(1,MSLB1UR0),PB1(1,MSLB1VR0),&
       & PB2(1,MSLB2URL),PB2(1,MSLB2VRL),&
       & PB2(1,MSLB2USI),PB2(1,MSLB2VSI),&
       & PUT1,PVT1 )
    ENDIF

  ENDIF

!     ------------------------------------------------------------------

!*       4.    VORTICITY EQUATION MODEL.
!              -------------------------

ELSEIF(LRVEQ)THEN
  DO JROF=IST,IEND
    PUT1(JROF)=PUT9(JROF)*ZTE+ZDT*(&
     & PVORT0(JROF)*PVT0(JROF)&
     & - PVT0L(JROF)*PVT0(JROF)&
     & - PUT0L(JROF)*PUT0(JROF)&
     & +(PUT0 (JROF)*PUT0(JROF)+PVT0(JROF)*PVT0(JROF))*0.5_JPRB&
     & *YDCSGEOM%RATATX(JROF)&
     & +YDGSGEOM%RCORI(JROF)*PVT0(JROF)    )
    PVT1(JROF)=PVT9(JROF)*ZTE+ZDT*(&
     & -PDIVT0(JROF)*PVT0(JROF)&
     & - PVT0L(JROF)*PUT0(JROF)&
     & + PUT0L(JROF)*PVT0(JROF)&
     & -(PUT0 (JROF)*PUT0(JROF)+PVT0(JROF)*PVT0(JROF))*0.5_JPRB&
     & *YDCSGEOM%RATATH(JROF)&
     & -YDGSGEOM%RCORI(JROF)*PUT0(JROF)    )
    PSPT1(JROF)=0.0_JPRB
  ENDDO

ENDIF

!     ------------------------------------------------------------------

!*       6.    MAP FACTOR
!              -----------

DO JROF=IST,IEND
  PUT0  (JROF)=PUT0  (JROF)/ YDGSGEOM%GM(JROF)
  PVT0  (JROF)=PVT0  (JROF)/ YDGSGEOM%GM(JROF)
  IF (LUVDER) THEN
    PUT0L (JROF)=PUT0L (JROF)/(YDGSGEOM%GM(JROF)**2)
    PVT0L (JROF)=PVT0L (JROF)/(YDGSGEOM%GM(JROF)**2)
  ENDIF
  PDIVT0(JROF)=PDIVT0(JROF)/(YDGSGEOM%GM(JROF)**2)
  IF (LVOR) THEN
    PVORT0(JROF)=PVORT0(JROF)/(YDGSGEOM%GM(JROF)**2)
  ENDIF
  PSPT0L(JROF)=PSPT0L(JROF)/ YDGSGEOM%GM(JROF)
  PSPT0M(JROF)=PSPT0M(JROF)/ YDGSGEOM%GM(JROF)
  IF (.NOT.YDDYNA%LELTRA) THEN
    PUT9  (JROF)=PUT9  (JROF)/ YDGSGEOM%GM(JROF)
    PVT9  (JROF)=PVT9  (JROF)/ YDGSGEOM%GM(JROF)
  ENDIF
  IF (.NOT.YDDYNA%LTWOTL) THEN
    PDIVT9(JROF)=PDIVT9(JROF)/(YDGSGEOM%GM(JROF)**2)
    PSPT9L(JROF)=PSPT9L(JROF)/ YDGSGEOM%GM(JROF)
    PSPT9M(JROF)=PSPT9M(JROF)/ YDGSGEOM%GM(JROF)
  ENDIF
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG2',1,ZHOOK_HANDLE)
END SUBROUTINE CPG2
