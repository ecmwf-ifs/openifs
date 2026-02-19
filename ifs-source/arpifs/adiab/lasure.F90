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

SUBROUTINE LASURE(&
 ! ----- INPUT ---------------------------------------------------------------
 & YDGEOMETRY,YDEPHY,YDDYN,YDEDYN,YDPHY,KSTART,KPROF,PBETADT,PDT,KIBL,&
 ! ----- OUTPUT --------------------------------------------------------------
 & PDTS2,PBT,LD2TLFF1,PBDT,PREDIV,PESGP,PESGM )  

!**** *LASURE*   Semi-Lagrangian scheme.
!                Set-up for other subroutines called by LACDYN.

!     Purpose.
!     --------
!       Computes some intermediate quantities necessary in the other
!       subroutines called by LACDYN.

!**   Interface.
!     ----------
!        *CALL* *LASURE(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KSTART  : first element of work.
!          KPROF   : depth of work.
!          PBETADT : BETADT or 0 according to configuration.
!          PDT     : time step for the first time-integration step of
!                    a leap-frog scheme or all time-integration steps of
!                    a two-time level scheme; 2*time step for the following
!                    time-integration steps of a leap-frog scheme.
!          KIBL    : index into YDGSGEOM in YDGEOMETRY

!        OUTPUT:
!          PDTS2   : 0.5*PDT.
!          PBT     : PDTS2*PBETADT.
!          LD2TLFF1: .T./.F.: Refined treatement of (2*Omega Vec r) at
!                    the origin point when there is t-dt (or t in SL2TL)
!                    physics / Other cases.
!          PBDT    : PBT or PBT*(c**2/GM**2) according to LSIDG.
!          PREDIV  : 1. or c**2/GM**2 according to LSIDG.
!          PESGP   : (1 + uncentering factor).
!          PESGM   : (1 - uncentering factor).

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-lagrangian scheme.

!     Externals.
!     ----------
!           none. 
!           Called by LACDYN, LACDYNTL and LACDYNAD.

!     Reference.
!     ----------
!             Arpege documentation about semi-lagrangian scheme.

!     Author.
!     -------
!      K. YESSAD (METEO FRANCE/CNRM/GMAP) after old part one of LACDYN. 
!      Loops are rewritten according to F90 norms.
!      Original : JULY 1995.

!     Modifications.
!     --------------
!      J.Vivoda 03-2002 PC schemes for NH dynamics (LPC_XXXX keys)
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad Aug 2008: rationalisation of dummy argument interfaces
!      K. Yessad (Dec 2011): use YDGSGEOM.
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOEPHY       , ONLY : TEPHY
USE YOMPHY       , ONLY : TPHY
USE YEMDYN       , ONLY : TEDYN
USE YOMDYN       , ONLY : TDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TEDYN)       ,INTENT(IN)    :: YDEDYN
TYPE(TEPHY)       ,INTENT(IN)    :: YDEPHY
TYPE(TPHY)        ,INTENT(IN)    :: YDPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBETADT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTS2 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBT 
LOGICAL           ,INTENT(OUT)   :: LD2TLFF1 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBDT(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PREDIV(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PESGP 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PESGM 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) ::  JROF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------


IF (LHOOK) CALL DR_HOOK('LASURE',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,&
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL))
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & L2TLFF=>YDDYN%L2TLFF, LSIDG=>YDDYN%LSIDG, VESL=>YDDYN%VESL, &
 & LESIDG=>YDEDYN%LESIDG, &
 & LAGPHY=>YDEPHY%LAGPHY, LEPHYS=>YDEPHY%LEPHYS, &
 & RSTRET=>YDGEM%RSTRET, LMPHYS=>YDPHY%LMPHYS)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS:
!              ----------------------------

!     * Time step.

PDTS2 = 0.5_JPRB*PDT
PBT = PDTS2*PBETADT

IF (LSIDG) THEN
  DO JROF=KSTART,KPROF
    PBDT(JROF)  = PBT
    PREDIV(JROF)= 1.0_JPRB
  ENDDO
ELSEIF (LESIDG) THEN
  DO JROF=KSTART,KPROF
    PBDT(JROF)  = PBT*YDGSGEOM%GMAPPA(JROF)/(YDGSGEOM%GM(JROF)*YDGSGEOM%GM(JROF))
    PREDIV(JROF)= YDGSGEOM%GMAPPA(JROF)/(YDGSGEOM%GM(JROF)*YDGSGEOM%GM(JROF))
  ENDDO
ELSE
  DO JROF=KSTART,KPROF
    PBDT(JROF)  = PBT*RSTRET*RSTRET/(YDGSGEOM%GM(JROF)*YDGSGEOM%GM(JROF))
    PREDIV(JROF)= RSTRET*RSTRET/(YDGSGEOM%GM(JROF)*YDGSGEOM%GM(JROF))
  ENDDO
ENDIF

!     * ( 1 + uncentering factor ) and other quantities linked with
!       uncentering factor.

PESGP = 1.0_JPRB + VESL
PESGM = 1.0_JPRB - VESL

!     * Switch to activate two additional interpolations for refined
!       computation of Coriolis term.

LD2TLFF1=L2TLFF.AND.(LMPHYS.OR.LEPHYS).AND.(.NOT.LAGPHY)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LASURE',1,ZHOOK_HANDLE)
END SUBROUTINE LASURE
