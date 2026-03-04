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

SUBROUTINE CPG_END(YDGEOMETRY,YDGMV,YDSURF,&
 !---------------------------------------------------------------------
 ! - INPUT .
 & YDEPHY,YDML_GCONF,YDDYN,YDDYNA,YDML_PHY_MF,LDTENC,KST,KEND,KSTC,KENDC,KBL,KSTGLO,LDCONFX,LDLDIAB,&
 & PGM,PQS,PSP_RR,PSP_SG,PSP_CL,&
 !---------------------------------------------------------------------
 ! - INPUT/OUTPUT .
 & PGFL,PGPAR,PGMV,PGMVS,PTRAJ_PHYS)

!**** *CPG_END* - Grid point calculations: end of non lagged part.

!     Purpose.
!     --------
!           Grid point calculations: end of non lagged part.
!           - save quantities at t and t-dt in buffers.
!           - save soil fields in a buffer.
!           - divide p-order derivatives (p=1 or 2) by M**p.

!           Abbreviation "vwv" stands for "vertical wind variable".

!**   Interface.
!     ----------
!        *CALL* *CPG_END(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        KST       : first element of work.
!        KEND      : last element of work.
!        KSTC,KENDC: the same as KST,KEND but including zone C for ALADIN.
!        KBL       : block number.
!        KSTGLO    : global offset.
!        LDCONFX   : (see in CPG)
!        LDLDIAB   : .T. if complete physics is activated, and predictor step.
!        PGM       : mapping factor.
!        PQS       : specific humidity at surface level.
!        PSP_RR    : surface prognostic quantities
!        PSP_CL    : surface diagnostic quantities (cls fields)
!        PSP_SG    : surface snow prognostic quantities

!     INPUT/OUTPUT:
!     -------------
!        PGFL      : GFL variables at time t-dt and t.
!        PGPAR     : surface fields for AROME.
!        PGMV      : upper air GMV variables at time t and t-dt.
!        PGMVS     : surface GMV variables at time t and t-dt.
!        PTRAJ_PHYS : Trajectory for physics (of next timestep)

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
!        K. YESSAD after the old part 9 of CPG. 

! Modifications
! -------------
!   Original : 01-08-14
!   04 JAN 2002: K.YESSAD - pruning; thin layer hyp relaxation for SL.
!   03-2002, J. Vivoda: PC schemes for NH dynamics (LPC_XXXX keys)
!   30-09-02 P.Smolikova : variable d4 in NH
!   12-10-2002, J. Masek : PC bugfix.
!   Modified 08-2002 C. Smith : use "w" as prognostic variable in the
!    semi-lag advection of vertical divergence in the NH model.
!   18 SEP 2002: C.MOUSSY - Surface pressure and its derivatives
!   31-01-03 P.Smolikova : variable d4 in 2tl PC NH
!   25-02-03 J.Vivoda : corrected test for PC-scheme
!   Modified 2003-07-17 C. Fischer - psvdauxt0* come from pgmv
!   01-Oct-2003 M. Hamrud   CY28 Cleaning
!   10-Jan-2004 M. Hamrud   CY28R1 Cleaning
!   Jan-2004 K. YESSAD   split GP_SPV into GP_SPV + GPMPFC_GMVS.
!   20-Jan-2004 JM. Audoin  reactive Tendency surface pressure
!   09-Jun-2004 J. Masek    NH cleaning (LFULLIMP)
!   01-Jul-2004 K. Yessad  Make clearer the tests for PC scheme.
!   2004-11-16 Y. Seity - save surface fields for AROME in GPARBUF 
!   M.Hamrud      01-Jul-2006 Revised surface fields
!   06-Mar-2006 R. Zaaboul  Add LMSE when saving PGPAR for AROME
!   17-Apr-2007 S.Ivatek-S: Over dimensioning of PGPAR to NGPAR+1 and PGPNH to
!                NFGPNH+1, boundary checking problem if NGPAR=0 or NFGPNH=0 bf
!   31-08-2007 Y.Seity : remove LMPA key for GPARBUF reading 
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (Nov 2009): prune lpc_old.
!   K. Yessad (Jan 2011): new architecture for LBC modules and set-up.
!   K. Yessad (Nov 2011): new GPMPFC.
!   K. Yessad (Nov 2012): simplify testings.
!   F. Vana  28-Nov-2013 : Redesigned trajectory handling
! End Modifications
!------------------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF, GPPOPER
USE YOMGMV   , ONLY : TGMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOEPHY   , ONLY : TEPHY
USE YOMDYN   , ONLY : TDYN
USE YOMDYNA  , ONLY : TDYNA
USE YOMGPPB  , ONLY : GPARBUF
USE YOMTRAJ  , ONLY : TRAJ_PHYS_TYPE

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)       ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)           ,INTENT(INOUT) :: YDGMV
TYPE(TSURF)          ,INTENT(INOUT) :: YDSURF
TYPE(TDYN)           ,INTENT(INOUT) :: YDDYN
TYPE(TDYNA)          ,INTENT(INOUT) :: YDDYNA
TYPE(TEPHY)          ,INTENT(INOUT) :: YDEPHY
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT):: YDML_PHY_MF
LOGICAL              ,INTENT(IN)    :: LDTENC
INTEGER(KIND=JPIM)   ,INTENT(IN)    :: KST 
INTEGER(KIND=JPIM)   ,INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM)   ,INTENT(IN)    :: KSTC 
INTEGER(KIND=JPIM)   ,INTENT(IN)    :: KENDC 
INTEGER(KIND=JPIM)   ,INTENT(IN)    :: KBL
INTEGER(KIND=JPIM)   ,INTENT(IN)    :: KSTGLO 
LOGICAL              ,INTENT(IN)    :: LDCONFX
LOGICAL              ,INTENT(IN)    :: LDLDIAB 
REAL(KIND=JPRB)      ,INTENT(IN)    :: PGM(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)      ,INTENT(IN)    :: PQS(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)      ,INTENT(IN)    :: PSP_RR(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSP_RRD%NDIM)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PSP_CL(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSP_CLD%NDIM)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PSP_SG(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSP_SGD%NDIM)
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGPAR(YDGEOMETRY%YRDIM%NPROMA,YDML_PHY_MF%YRPARAR%NGPAR+1) 
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV) 
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS)
TYPE (TRAJ_PHYS_TYPE),INTENT(INOUT) :: PTRAJ_PHYS
!     ------------------------------------------------------------------
LOGICAL    :: LLSTR

INTEGER(KIND=JPIM) :: IFLAG
REAL(KIND=JPRB) :: ZEPS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gpmpfc.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPG_END',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
  & YDSIMPHL=>YDML_PHY_MF%YRSIMPHL,YDARPHY=>YDML_PHY_MF%YRARPHY,YGFL=>YDML_GCONF%YGFL,YDPARAR=>YDML_PHY_MF%YRPARAR)
ASSOCIATE(NDIM=>YGFL%NDIM, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
 & LAGPHY=>YDEPHY%LAGPHY, &
 & RSTRET=>YDGEM%RSTRET, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0, &
 & YT9=>YDGMV%YT9, &
 & YSP_RR=>YDSURF%YSP_RR, YSP_RRD=>YDSURF%YSP_RRD, YSP_SGD=>YDSURF%YSP_SGD, &
 & YSP_CL=>YDSURF%YSP_CL, YSP_CLD=>YDSURF%YSP_CLD, &
 & NGPAR=>YDPARAR%NGPAR, LMSE=>YDARPHY%LMSE, LTRAJPS=>YDSIMPHL%LTRAJPS)
!     ------------------------------------------------------------------

!*    1.    PRELIMINARY INITIALISATIONS.
!           ----------------------------

ZEPS=100.0_JPRB*TINY(1.0_JPRB)
LLSTR=(ABS(RSTRET-1.0_JPRB)>ZEPS)

!     ------------------------------------------------------------------

!*    2.    SURFACE FIELDS TIMESTEPPING.
!           ----------------------------

! NO TIME-STEPPING FOR SURFACE FIELDS IF ADIABATIC CALL OR NHS ITERATION.
IF (NCURRENT_ITER == 0) THEN
  IF ((.NOT.LDCONFX).AND.LDLDIAB) THEN  
    IF (YDDYNA%LTWOTL) THEN
      CALL GPPOPER(YDDYN,'SET0TO1',YDSURF,KBL=KBL)
    ELSE
      CALL GPPOPER(YDDYN,'PHTFILT',YDSURF,KBL=KBL)
    ENDIF
  ENDIF
ENDIF

IF (LTRAJPS) THEN  
  PTRAJ_PHYS%PQSS1MF(KST:KEND)=PQS(KST:KEND)
  PTRAJ_PHYS%PTS1MF(KST:KEND) =PSP_RR(KST:KEND,YSP_RR%YT%MP9)
ENDIF

!     ------------------------------------------------------------------

!*    3.    DIVIDE BY MAPPING FACTOR.
!           -------------------------

IF (LLSTR.AND.(.NOT.LAGPHY)) THEN
  IFLAG=1
  CALL GPMPFC(YDGMV,YDML_GCONF,YDDYN,YDDYNA,NPROMA,NFLEVG,KST,KEND,IFLAG,PGM,PGMV,PGMVS,PGFL)
ENDIF

!     ------------------------------------------------------------------

!*    4.    SURFACE PRESSURE AND ITS DERIVATIVES.
!           -------------------------------------

IF (YDDYNA%LTWOTL.AND.LDTENC) THEN
  PGMVS(KST:KEND,YT9%MSP)=PGMVS(KST:KEND,YT0%MSP)
  PGMVS(KST:KEND,YT9%MSPL)=PGMVS(KST:KEND,YT0%MSPL)
  PGMVS(KST:KEND,YT9%MSPM)=PGMVS(KST:KEND,YT0%MSPM)
ENDIF

!     ------------------------------------------------------------------

!*    5.    "WRITE" TO WORK FILES.
!           ----------------------

! save Arome surface fields in GPARBUF
IF(NCURRENT_ITER == 0)THEN
  IF (LMSE.AND.NGPAR/=0) THEN
    GPARBUF (:, :, 1+(KSTGLO-1)/YDGEOMETRY%YRDIM%NPROMA) = PGPAR(1:NPROMA,1:NGPAR)

  ENDIF
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG_END',1,ZHOOK_HANDLE)
END SUBROUTINE CPG_END
