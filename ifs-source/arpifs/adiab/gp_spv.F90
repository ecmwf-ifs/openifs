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

SUBROUTINE GP_SPV(&
 ! - INPUT --------------------------------------------------------------
 & YDGEOMETRY, YDDYN,YDDYNA,YDSIMPHL,LDTL,KST,KEND,&
 & PSPT0,PSPT0L,PSPT0M,PSPT9,PSPT9L,PSPT9M,&
 ! - OUTPUT -------------------------------------------------------------
 & PRE0,PRE0L,PRE0M,PRE9,PRE9L,PRE9M)

!     ------------------------------------------------------------------
!**** *GP_SPV* - 

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *GP_SPV(...)*

!        Explicit arguments :
!        --------------------

!        abbreviation "prehyds" means "hydrostatic surface pressure".

!        INPUT:
!        ------
!        LDTL               : true if gp_spv is called from TL/AD models
!        KST                : start of work
!        KEND               : depth of work
!        PSPT0,PSPT0L,PSPT0M: log(prehyds) and derivatives at t.
!        PSPT9,PSPT9L,PSPT9M: log(prehyds) and derivatives at t-dt.

!        OUTPUT:
!        -------
!        PRE0,PRE0L,PRE0M   : prehyds and derivatives at t.
!        PRE9,PRE9L,PRE9M   : prehyds and derivatives at t-dt.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
       
!     Author.
!     -------

! Modifications
! -------------
!   07-Aug-2001 R. El Khatib  Pruning options
!      Mar-2002 J. Vivoda     PC schemes for NH dynamics (LPC_XXXX keys)
!   27-Jun-2002 C. Fischer    cdlock & ldtl
!   01-Oct-2003 M. Hamrud     CY28 Cleaning
!      Dec-2003 K. Yessad     multiplication by GM has moved in GPMPFC_GMVS.
!   09-Jun-2004 J. Masek      NH cleaning (LFULLIMP)
!   01-Jul-2004 K. Yessad     Make clearer the tests for PC scheme.
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (Nov 2012): simplify testings.
!   K. Yessad (July 2014): Move some variables.
! End Modifications
!------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : TDYNA
USE YOMSIMPHL    , ONLY : TSIMPHL

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TDYNA)       ,INTENT(IN)    :: YDDYNA
TYPE(TSIMPHL)     ,INTENT(IN)    :: YDSIMPHL
LOGICAL           ,INTENT(IN)    :: LDTL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0L(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0M(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT9(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT9L(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT9M(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0L(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0M(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE9L(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE9M(YDGEOMETRY%YRDIM%NPROMA) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF
LOGICAL :: LLPRE9
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GP_SPV',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
 & LSIMPH=>YDSIMPHL%LSIMPH)
!     ------------------------------------------------------------------

IF (.NOT.(((LSIMPH.OR.YDDYNA%LNHDYN).AND.LDTL).OR..NOT.LDTL)) THEN
  LLPRE9=.FALSE.
ELSE
  IF (YDDYNA%LTWOTL) THEN
    LLPRE9=(NCURRENT_ITER > 0 .AND. YDDYNA%LPC_FULL)
  ELSE
    LLPRE9=.TRUE.
  ENDIF
ENDIF

DO JROF=KST,KEND
  PRE0 (JROF,NFLEVG)=EXP(PSPT0(JROF))
  PRE0L(JROF       )=PSPT0L(JROF)*PRE0(JROF,NFLEVG)
  PRE0M(JROF       )=PSPT0M(JROF)*PRE0(JROF,NFLEVG)
ENDDO
IF (LLPRE9) THEN
  DO JROF=KST,KEND
    PRE9 (JROF,NFLEVG)=EXP(PSPT9(JROF))
    PRE9L(JROF       )=PSPT9L(JROF)*PRE9(JROF,NFLEVG)
    PRE9M(JROF       )=PSPT9M(JROF)*PRE9(JROF,NFLEVG)
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GP_SPV',1,ZHOOK_HANDLE)
END SUBROUTINE GP_SPV
