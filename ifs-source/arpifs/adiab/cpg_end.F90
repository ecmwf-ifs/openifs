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

SUBROUTINE CPG_END(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDGMV, YDSURF, YDEPHY, YDML_GCONF, YDDYN, YDDYNA, &
& YDML_PHY_MF, LDTENC, YDDPHY, PGM, PQS, PRRT9, PGFL, PGPAR, PGMV, PGMVS, PTRAJ_PHYS)


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
USE CPG_OPTS_TYPE_MOD, ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE SURFACE_FIELDS_MIX , ONLY : TSURF, GPPOPER, GPOPER_2
USE YOMGMV   , ONLY : TGMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOEPHY   , ONLY : TEPHY
USE YOMDYN   , ONLY : TDYN
USE YOMDYNA  , ONLY : TDYNA
USE YOMGPPB  , ONLY : GPARBUF
USE YOMTRAJ  , ONLY : TRAJ_PHYS_TYPE
USE YOMDPHY  , ONLY : TDPHY
USE ALGORITHM_STATE_MOD, ONLY : GET_NUPTRA

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)       ,INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_BNDS_TYPE),   INTENT(IN)   :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE),  INTENT(IN)    :: YDCPG_OPTS
TYPE(TGMV)           ,INTENT(INOUT) :: YDGMV
TYPE(TSURF)          ,INTENT(INOUT) :: YDSURF
TYPE(TDYN)           ,INTENT(INOUT) :: YDDYN
TYPE(TDYNA)          ,INTENT(INOUT) :: YDDYNA
TYPE(TEPHY)          ,INTENT(INOUT) :: YDEPHY
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT):: YDML_PHY_MF
TYPE(TDPHY)          ,INTENT(IN) :: YDDPHY
LOGICAL              ,INTENT(IN)    :: LDTENC
REAL(KIND=JPRB)      ,INTENT(IN)    :: PGM(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)      ,INTENT(IN)    :: PQS(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)      ,INTENT(IN)    :: PRRT9(YDGEOMETRY%YRDIM%NPROMA)
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

IF (LHOOK) CALL DR_HOOK('CPG_END', 0, ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDSIMPHL=>YDML_PHY_MF%YRSIMPHL, YDARPHY=>YDML_PHY_MF%YRARPHY, YGFL=>YDML_GCONF%YGFL, YDPARAR=>YDML_PHY_MF%YRPARAR)
ASSOCIATE(NDIM=>YGFL%NDIM, NPROMA=>YDDIM%NPROMA, NFLEVG=>YDDIMV%NFLEVG, NCURRENT_ITER=>YDDYN%NCURRENT_ITER,       &
 & LAGPHY=>YDEPHY%LAGPHY, RSTRET=>YDGEM%RSTRET, NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0,  &
 & YT9=>YDGMV%YT9, YSP_RR=>YDSURF%YSP_RR, YSP_RRD=>YDSURF%YSP_RRD, YSP_SGD=>YDSURF%YSP_SGD, YSP_CL=>YDSURF%YSP_CL, &
 & YSP_CLD=>YDSURF%YSP_CLD, NGPAR=>YDPARAR%NGPAR, LMSE=>YDARPHY%LMSE, LTRAJPS=>YDSIMPHL%LTRAJPS, &
 & LDIRCLSMOD=>YDDPHY%LDIRCLSMOD, &
 & KBL=>YDCPG_BNDS%KBL, KIDIA=>YDCPG_BNDS%KIDIA, KFDIA=>YDCPG_BNDS%KFDIA, KSTGLO=>YDCPG_BNDS%KSTGLO)

!     ------------------------------------------------------------------

!*    1.    PRELIMINARY INITIALISATIONS.
!           ----------------------------

ZEPS=100.0_JPRB*TINY(1.0_JPRB)
LLSTR=(ABS(RSTRET-1.0_JPRB)>ZEPS)

!     ------------------------------------------------------------------

!*    2.    SURFACE FIELDS TIMESTEPPING.
!           ----------------------------

! NO TIME-STEPPING FOR SURFACE FIELDS IF ADIABATIC CALL OR NHS ITERATION.
!  EXCEPT IN CASE OF ASSIMILATION FOR TAKING INTO ACCOUNT THE NEW TRAJECTORY
!  UPDATED WITH THE LAST INCREMENT (CLS FIELDS ONLY)
IF (NCURRENT_ITER == 0) THEN
  IF (YDCPG_OPTS%LDIAB) THEN  
    IF (.NOT.YDCPG_OPTS%LCONFX) THEN  
      IF (YDDYNA%LTWOTL) THEN
        CALL GPPOPER(YDDYN, 'SET0TO1', YDSURF, KBL=KBL)
      ELSE
        CALL GPPOPER(YDDYN, 'PHTFILT', YDSURF, KBL=KBL)
      ENDIF
    ELSEIF (YDDYNA%LTWOTL.AND.LDIRCLSMOD.AND.GET_NUPTRA()>0) THEN
        CALL GPOPER_2(YDDYN,'SET0TO1',YDSURF%SP_CL(:,:,KBL),YDSURF%YSP_CLD,YDSURF%YSP_CL%YCL)
    ENDIF
  ENDIF
ENDIF

IF (LTRAJPS) THEN  
  PTRAJ_PHYS%PQSS1MF(KIDIA:KFDIA)=PQS(KIDIA:KFDIA)
  PTRAJ_PHYS%PTS1MF(KIDIA:KFDIA)=PRRT9(KIDIA:KFDIA)
ENDIF

!     ------------------------------------------------------------------

!*    3.    DIVIDE BY MAPPING FACTOR.
!           -------------------------

IF (LLSTR.AND.(.NOT.LAGPHY)) THEN
  IFLAG=1
  CALL GPMPFC(YDGMV, YDML_GCONF, YDDYN, YDDYNA, NPROMA, NFLEVG, KIDIA, KFDIA, IFLAG, &
  & PGM, PGMV, PGMVS, PGFL)
ENDIF

!     ------------------------------------------------------------------

!*    4.    SURFACE PRESSURE AND ITS DERIVATIVES.
!           -------------------------------------

IF (YDDYNA%LTWOTL.AND.LDTENC) THEN
  PGMVS(KIDIA:KFDIA,YT9%MSP)=PGMVS(KIDIA:KFDIA,YT0%MSP)
  PGMVS(KIDIA:KFDIA,YT9%MSPL)=PGMVS(KIDIA:KFDIA,YT0%MSPL)
  PGMVS(KIDIA:KFDIA,YT9%MSPM)=PGMVS(KIDIA:KFDIA,YT0%MSPM)
ENDIF

!     ------------------------------------------------------------------

!*    5.    "WRITE" TO WORK FILES.
!           ----------------------

! save Arome surface fields in GPARBUF
IF(NCURRENT_ITER == 0)THEN
  IF (LMSE.AND.NGPAR/=0) THEN
    GPARBUF (:, :, 1+(KSTGLO-1)/NPROMA) = PGPAR(1:NPROMA,1:NGPAR)

  ENDIF
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG_END', 1, ZHOOK_HANDLE)
END SUBROUTINE CPG_END
