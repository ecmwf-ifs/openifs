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

SUBROUTINE SPCIMPFSOLVE(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY,YDRIP,YDDYN,LDNHDYN,LDNHX,LDONEM,PSDIVP,&
 ! --- OUTPUT ----------------------------------------------------------------
 & PSPDIVP)  

!**** *SPCIMPFSOLVE* - CALL SIMPLICO FROM S-SPACE

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SPCIMPFSOLVE(..)

!        Explicit arguments :  LDONEM  - T if only one m if processed
!                              PSDIVP  - right hand side
!                              PSPDIVP - solution

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Tomas Wilhelmsson  *ECMWF*
!        21-09-2009 Original with extracted code from spcsi.F90

!     Modifications.
!     --------------
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      O. Marsden (May 2016): Remove redundant geometry argument
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST       , ONLY : ROMEGA
USE YOMDYN       , ONLY : TDYN
USE YOMRIP       , ONLY : TRIP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
LOGICAL           ,INTENT(IN)    :: LDNHDYN
LOGICAL           ,INTENT(IN)    :: LDNHX
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
LOGICAL           ,INTENT(IN)    :: LDONEM 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSDIVP (:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSPDIVP(:,:)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZALPHA (0:YDGEOMETRY%YRDIM%NSMAX+1)
REAL(KIND=JPRB) :: ZDENIM (0:YDGEOMETRY%YRDIM%NSMAX+1)
REAL(KIND=JPRB) :: ZEPSI  (0:YDGEOMETRY%YRDIM%NSMAX)
REAL(KIND=JPRB) :: ZFPLUS (0:YDGEOMETRY%YRDIM%NSMAX+1)
REAL(KIND=JPRB) :: ZFMINUS(0:YDGEOMETRY%YRDIM%NSMAX+1)
REAL(KIND=JPRB) :: ZSIVP(YDGEOMETRY%YRDIMV%NFLEVL)

REAL(KIND=JPRB),ALLOCATABLE :: ZSDIVP (:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZSPDIVP(:,:)

INTEGER(KIND=JPIM) :: ILO, IM, ISTA, IEND, JL, JMLOC

REAL(KIND=JPRB) :: ZAL, ZBDT, ZBDT2, ZEM, ZEN, ZF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "simplico.h"
#include "trmtos.intfb.h"
#include "trstom.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPCIMPFSOLVE',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDLAP=>YDGEOMETRY%YRLAP)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, NSPEC2=>YDDIM%NSPEC2, NUMP=>YDDIM%NUMP, &
 & NFLEVL=>YDDIMV%NFLEVL, NFLSUR=>YDDIMV%NFLSUR, &
 & RBT=>YDDYN%RBT, SIVP=>YDDYN%SIVP, &
 & RSTRET=>YDGEM%RSTRET, &
 & MYLEVS=>YDMP%MYLEVS, &
 & TDT=>YDRIP%TDT)
!     ------------------------------------------------------------------

ALLOCATE(ZSDIVP (NFLSUR,NSPEC2))
ALLOCATE(ZSPDIVP(NFLSUR,NSPEC2))

!     ------------------------------------------------------------------

!*       1.    SEMI-IMPLICIT SPECTRAL COMPUTATIONS.
!              ------------------------------------

!*        1.1  Preliminary initialisations.

ZBDT=0.5_JPRB*TDT*RBT
ZBDT2=(ZBDT*RSTRET)**2

DO JL=1,NFLEVL
  ZSIVP(JL) = SIVP(MYLEVS(JL))
ENDDO

!*        1.2  Main calculations.

CALL GSTATS(9,2)
CALL GSTATS(55,0)
CALL TRSTOM(YDGEOMETRY,LDNHDYN,LDNHX,PSPDIVG=PSDIVP,PSPDIV=ZSDIVP,LDFULLM=LDONEM)
CALL GSTATS(55,2)
CALL GSTATS(9,3)

CALL GSTATS(1028,0)
!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JMLOC,IM,ISTA,IEND) &
!$OMP&PRIVATE(ZEM,ZAL,ILO,ZALPHA,ZDENIM,ZEPSI,ZEN,ZFPLUS,ZFMINUS)
DO JMLOC=1,NUMP

  IM=YDLAP%MYMS(JMLOC)
  ISTA=YDLAP%NASM0(IM)
  IEND=ISTA+2*(NSMAX+1-IM)-1

  !*             Set up helper arrays for implicit Coriolis.

  ZEM=REAL(IM,JPRB)
  ZAL=2.0_JPRB*ZBDT*ROMEGA*ZEM
  ILO=IM
  IF (IM == 0) THEN
    ZALPHA(0)=0.0_JPRB
    ZDENIM(0)=0.0_JPRB
    ZEPSI(0)=0.0_JPRB
    ILO=1
  ENDIF
  DO JL=ILO,NSMAX
    ZEN=REAL(JL,JPRB)
    ZALPHA(JL)=ZAL/(ZEN*(ZEN+1.0_JPRB))
    ZDENIM(JL)=1.0_JPRB/(1.0_JPRB+ZALPHA(JL)**2)
    ZEPSI(JL)=SQRT((ZEN*ZEN-ZEM*ZEM)/(4.0_JPRB*ZEN*ZEN-1.0_JPRB))
  ENDDO
  ZALPHA(NSMAX+1)=0.0_JPRB
  ZDENIM(NSMAX+1)=0.0_JPRB

  IF (IM == 0) THEN
    ZFPLUS(0)=0.0_JPRB
    ZFMINUS(0)=0.0_JPRB
  ENDIF
  ZF=2.0_JPRB*ZBDT*ROMEGA
  DO JL=ILO,NSMAX-1
    ZEN=REAL(JL,JPRB)
    ZFPLUS(JL)=ZF*ZEN*ZEPSI(JL+1)/(ZEN+1.0_JPRB)
    ZFMINUS(JL)=ZF*(ZEN+1.0_JPRB)*ZEPSI(JL)/ZEN
  ENDDO
  ZEN=REAL(NSMAX,JPRB)
  ZFPLUS(NSMAX)=0.0_JPRB
  ZFMINUS(NSMAX)=ZF*(ZEN+1.0_JPRB)*ZEPSI(NSMAX)/ZEN
  ZFPLUS(NSMAX+1)=0.0_JPRB
  ZFMINUS(NSMAX+1)=0.0_JPRB

  !           Solve complex pentadiagonal system

  CALL SIMPLICO(IM,NSMAX,NFLEVL,NFLSUR,ZALPHA(IM),&
   & ZDENIM(IM),ZFPLUS(IM),ZFMINUS(IM),ZSIVP,YDLAP%RLAPDI(0:),&
   & ZBDT2,ZSDIVP(1,ISTA),ZSPDIVP(1,ISTA))

ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1028,1)

CALL GSTATS(9,2)
CALL GSTATS(55,3)
CALL TRMTOS(YDGEOMETRY,LDNHDYN,LDNHX,PSPDIV=ZSPDIVP,PSPDIVG=PSPDIVP,LDFULLM=LDONEM)  
CALL GSTATS(55,1)
CALL GSTATS(9,3)

!     ------------------------------------------------------------------

DEALLOCATE(ZSDIVP)
DEALLOCATE(ZSPDIVP)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPCIMPFSOLVE',1,ZHOOK_HANDLE)
END SUBROUTINE SPCIMPFSOLVE
