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

SUBROUTINE SUALDYNB(YDGEOMETRY,YDML_GCONF,YDDYN,YDDYNA)

!**** *SUALDYNB* - Routine to allocate space for dynamics.

!     Purpose.
!     --------
!           Allocate space for the dynamics calc.
!           Arrays which need a separate allocation in LAM model.

!**   Interface.
!     ----------
!        *CALL* *SUALDYNB*

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
!        Karim YESSAD
!        Original : 95-11-20  From SUALDYN

!     Modifications.
!     --------------
!   F. Vana 13-Jan-2009 : Recognized special case LSLHD_STATIC
!   K. Yessad (July 2014): Move some variables.
!   F. Vana and M. Diamantakis (Aug 2016): Simplified LSETTLSVF scheme
!   K. Yessad (June 2017): Introduce NHQE model.
!   K. Yessad (June 2018): Alternate NHEE SI scheme elimination.
!   M. Diamantakis 2020  : Save DP to accelerate convergence at next step
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : NPRINTLEV
USE YOMCT0   , ONLY : LALLOPR
USE YOMLUN   , ONLY : NULOUT
USE YOMDYN   , ONLY : TDYN
USE YOMDYNA  , ONLY : TDYNA
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT) :: YDDYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(TDYNA),INTENT(IN):: YDDYNA
INTEGER(KIND=JPIM) :: IU,IVTHS,IL

LOGICAL :: LLP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUALDYNB',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
 & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDDIMF=>YDML_GCONF%YRDIMF, &
 & YGFL=>YDML_GCONF%YGFL)
ASSOCIATE(NF3D=>YDDIMF%NF3D, &
 & NSPEC=>YDDIM%NSPEC, NSMAX=>YDDIM%NSMAX, &
 & NFLEVL=>YDDIMV%NFLEVL, NFLEVG=>YDDIMV%NFLEVG, &
 & NGPTOT=>YDGEM%NGPTOT, NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, &
 & LSIDG=>YDDYN%LSIDG, LSTRHD=>YDDYN%LSTRHD, &
 & LSLHD_STATIC=>YDDYNA%LSLHD_STATIC,LSLHD_W=>YDDYNA%LSLHD_W,LSLHD_SVD=>YDDYNA%LSLHD_SVD, &
 & LNHX=>YDDYNA%LNHX,LSI_NHEE=>YDDYNA%LSI_NHEE, &
 & LSLDP_XYZ=>YDDYN%LSLDP_XYZ, LSLDP_SAVE=>YDDYN%LSLDP_SAVE)
!     ------------------------------------------------------------------

!*       1.    ALLOCATE SPACE FOR ARRAYS.
!              --------------------------

LLP = NPRINTLEV >= 1.OR. LALLOPR
IU = NULOUT

! * Semi-Lagrangian departure point saving to be used at next timestep
IF (LSLDP_XYZ.AND.LSLDP_SAVE) THEN
  ALLOCATE(YDDYN%RSAVEDP(NPROMA,NFLEVG,3,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'RSAVEDP    ',SIZE(YDDYN%RSAVEDP),SHAPE(YDDYN%RSAVEDP)
ENDIF

! * Semi-implicit scheme for ARPEGE stretched geometry.
IF (LSIDG) THEN
  ALLOCATE(YDDYN%SIHEG(NFLEVG,NSPEC,3))
  IF(LLP)WRITE(IU,9) 'SIHEG    ',SIZE(YDDYN%SIHEG),SHAPE(YDDYN%SIHEG)
  ALLOCATE(YDDYN%SIHEG2(NFLEVG,NSMAX+1,2:3))
  IF(LLP)WRITE(IU,9) 'SIHEG2   ',SIZE(YDDYN%SIHEG2),SHAPE(YDDYN%SIHEG2)
  IF (YDDYNA%LNHEE.AND.(.NOT.LSI_NHEE)) THEN
    ALLOCATE(YDDYN%SIHEGB(NFLEVG,NSPEC,3))
    IF(LLP)WRITE(IU,9) 'SIHEGB   ',SIZE(YDDYN%SIHEGB),SHAPE(YDDYN%SIHEGB)
    ALLOCATE(YDDYN%SIHEGB2(NFLEVG,NSMAX+1,2:3))
    IF(LLP)WRITE(IU,9) 'SIHEGB2  ',SIZE(YDDYN%SIHEGB2),SHAPE(YDDYN%SIHEGB2)
  ENDIF
ENDIF

! * Horizontal diffusion scheme.
IF (LSTRHD) THEN
  ALLOCATE(YDDYN%RDHI  (NFLEVL*NF3D+1,NSPEC,3))
  IF(LLP)WRITE(IU,9) 'RDHI     ',SIZE(YDDYN%RDHI),SHAPE(YDDYN%RDHI)

  ! * Compute the number of prognostic variables on which the HDS
  !   diffusion is applied (currently IVTHS=(IL/NFLEVL) is between 0 and 4).
  !   Value of IVTHS must be the same in SUALDYNB and SUHDU.
  IVTHS=0
  IF ((.NOT.LSLHD_STATIC).AND.LSLHD_W) IVTHS=IVTHS+2
  IF ((.NOT.LSLHD_STATIC).AND.YDDYNA%LNHDYN.AND.LSLHD_SVD) IVTHS=IVTHS+1
  IF ((.NOT.LSLHD_STATIC).AND.LNHX.AND.LSLHD_SVD) IVTHS=IVTHS+1
  IL=IVTHS*NFLEVL
  IF(IVTHS>0) THEN
    ALLOCATE(YDDYN%RDHS  (IL,NSPEC,3))
    IF(LLP)WRITE(IU,9) 'RDHS     ',SIZE(YDDYN%RDHS),SHAPE(YDDYN%RDHS)
  ENDIF
ENDIF
ALLOCATE(YDDYN%RDIVOR(NFLEVL,0:NSMAX))
IF(LLP)WRITE(IU,9) 'RDIVOR   ',SIZE(YDDYN%RDIVOR),SHAPE(YDDYN%RDIVOR)
ALLOCATE(YDDYN%RDIDIV(NFLEVL,0:NSMAX))
IF(LLP)WRITE(IU,9) 'RDIDIV   ',SIZE(YDDYN%RDIDIV),SHAPE(YDDYN%RDIDIV)
ALLOCATE(YDDYN%RDITG (NFLEVG,0:NSMAX))
IF(LLP)WRITE(IU,9) 'RDITG    ',SIZE(YDDYN%RDITG),SHAPE(YDDYN%RDITG)
ALLOCATE(YDDYN%RDIGFL(NFLEVL,0:NSMAX,YGFL%NUMSPFLDS))
IF(LLP)WRITE(IU,9) 'RDIGFL   ',SIZE(YDDYN%RDIGFL),SHAPE(YDDYN%RDIGFL)
ALLOCATE(YDDYN%RDISP (      0:NSMAX))
IF(LLP)WRITE(IU,9) 'RDISP    ',SIZE(YDDYN%RDISP),SHAPE(YDDYN%RDISP)
IF(YDDYNA%LNHDYN) THEN
  ALLOCATE(YDDYN%RDIPD(NFLEVL,0:NSMAX))
  IF(LLP)WRITE(IU,9) 'RDIPD    ',SIZE(YDDYN%RDIPD),SHAPE(YDDYN%RDIPD)
  ALLOCATE(YDDYN%RDIVD(NFLEVL,0:NSMAX))
  IF(LLP)WRITE(IU,9) 'RDIVD    ',SIZE(YDDYN%RDIVD),SHAPE(YDDYN%RDIVD)
ENDIF
IF((.NOT.LSLHD_STATIC).AND.LSLHD_W) THEN
  ALLOCATE(YDDYN%RDSVOR(NFLEVL,0:NSMAX))
  IF(LLP)WRITE(IU,9) 'RDSVOR   ',SIZE(YDDYN%RDSVOR),SHAPE(YDDYN%RDSVOR)
  ALLOCATE(YDDYN%RDSDIV(NFLEVL,0:NSMAX))
  IF(LLP)WRITE(IU,9) 'RDSDIV   ',SIZE(YDDYN%RDSDIV),SHAPE(YDDYN%RDSDIV)
ENDIF
IF((.NOT.LSLHD_STATIC).AND.LSLHD_SVD.AND.YDDYNA%LNHDYN) THEN
  ALLOCATE(YDDYN%RDSVD(NFLEVL,0:NSMAX))
  IF(LLP)WRITE(IU,9) 'RDSVD    ',SIZE(YDDYN%RDSVD),SHAPE(YDDYN%RDSVD)
ENDIF

! * SLHD scheme
! ky: for SLHD.. arrays it is recommended to allocate them even if LSLHD=F,
!     to avoid false alarms in "bound checking" runs,
!     because they are passed to SLag routines via dummy arg.
ALLOCATE(YDDYN%SLHDA(NGPTOT))
ALLOCATE(YDDYN%SLHDD0(NGPTOT))

!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUALDYNB',1,ZHOOK_HANDLE)
END SUBROUTINE SUALDYNB
