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

!OCL  NOPREEX
SUBROUTINE SUALGCO(YDGEM,YDDPHY)

!**** *SUALGCO*  - Allocate arrays for coupled fields diagnostics

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SUALGCO

!        Explicit arguments :
!        --------------------
!        None

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
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15 (in SUGEM2)

!     Modifications.
!     --------------
!      K. Yessad (June 2013): move code in SUALGCO
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE YOMGEM   , ONLY : TGEM
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMMP0   , ONLY : NPRINTLEV
USE YOMCT0   , ONLY : LALLOPR
USE YOMLUN   , ONLY : NULOUT
USE YOMGCO   , ONLY : STRSU    ,STRSV    ,FCHAS    ,FRSOS    ,&
 & FHUMS    ,RUIST    ,FCHLL    ,FCHLN    ,FCHSS    ,&
 & FHUML    ,FHUMN    ,FRLDS    ,QWPRO    ,TMAX2    ,&
 & TMIN2    ,HSNOW    ,TSUR2    ,TSTS     ,TSFL     ,&
 & XPUQ     ,XPVQ     ,FCRFTH   ,FCRFSO
USE YOMDPHY  , ONLY : TDPHY

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TGEM) , INTENT(IN) :: YDGEM
TYPE(TDPHY) ,INTENT(INOUT):: YDDPHY
INTEGER(KIND=JPIM) :: IU
LOGICAL :: LLP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUALGCO',0,ZHOOK_HANDLE)
ASSOCIATE(NTSSG=>YDDPHY%NTSSG, &
 & NGPTOT=>YDGEM%NGPTOT)
!     ------------------------------------------------------------------

!*       1.    ALLOCATE ARRAYS FOR COUPLED FIELDS DIAGNOSTICS
!              ----------------------------------------------

LLP = NPRINTLEV >= 1.OR. LALLOPR
IU = NULOUT

ALLOCATE(STRSU(NGPTOT))
STRSU(:)=0._JPRB
IF(LLP)WRITE(IU,9) 'STRSU    ',SIZE(STRSU),SHAPE(STRSU)
ALLOCATE(STRSV(NGPTOT))
STRSV(:)=0._JPRB
IF(LLP)WRITE(IU,9) 'STRSV    ',SIZE(STRSV),SHAPE(STRSV)
ALLOCATE(FCHAS(NGPTOT,NTSSG+1))
FCHAS(:,:)=0._JPRB
IF(LLP)WRITE(IU,9) 'FCHAS    ',SIZE(FCHAS),SHAPE(FCHAS)
ALLOCATE(FRSOS(NGPTOT,NTSSG+1))
FRSOS(:,:)=0._JPRB
IF(LLP)WRITE(IU,9) 'FRSOS    ',SIZE(FRSOS),SHAPE(FRSOS)
ALLOCATE(FHUMS(NGPTOT,NTSSG+1))
FHUMS(:,:)=0._JPRB
IF(LLP)WRITE(IU,9) 'FHUMS    ',SIZE(FHUMS),SHAPE(FHUMS)
ALLOCATE(RUIST(NGPTOT))
RUIST(:)=0._JPRB
IF(LLP)WRITE(IU,9) 'RUIST    ',SIZE(RUIST),SHAPE(RUIST)
ALLOCATE(FCHLL(NGPTOT,NTSSG+1))
FCHLL(:,:)=0._JPRB
IF(LLP)WRITE(IU,9) 'FCHLL    ',SIZE(FCHLL),SHAPE(FCHLL)
ALLOCATE(FCHLN(NGPTOT,NTSSG+1))
FCHLN(:,:)=0._JPRB
IF(LLP)WRITE(IU,9) 'FCHLN    ',SIZE(FCHLN),SHAPE(FCHLN)
ALLOCATE(FCHSS(NGPTOT,NTSSG+1))
FCHSS(:,:)=0._JPRB
IF(LLP)WRITE(IU,9) 'FCHSS    ',SIZE(FCHSS),SHAPE(FCHSS)
ALLOCATE(FHUML(NGPTOT))
FHUML(:)=0._JPRB
IF(LLP)WRITE(IU,9) 'FHUML    ',SIZE(FHUML),SHAPE(FHUML)
ALLOCATE(FHUMN(NGPTOT))
FHUMN(:)=0._JPRB
IF(LLP)WRITE(IU,9) 'FHUMN    ',SIZE(FHUMN),SHAPE(FHUMN)
ALLOCATE(FRLDS(NGPTOT))
FRLDS(:)=0._JPRB
IF(LLP)WRITE(IU,9) 'FRLDS    ',SIZE(FRLDS),SHAPE(FRLDS)
ALLOCATE(QWPRO(NGPTOT))
QWPRO(:)=0._JPRB
IF(LLP)WRITE(IU,9) 'QWPRO    ',SIZE(QWPRO),SHAPE(QWPRO)
ALLOCATE(TMAX2(NGPTOT))
TMAX2(:)=-1000._JPRB
IF(LLP)WRITE(IU,9) 'TMAX2    ',SIZE(TMAX2),SHAPE(TMAX2)
ALLOCATE(TMIN2(NGPTOT))
TMIN2(:)=1000._JPRB
IF(LLP)WRITE(IU,9) 'TMIN2    ',SIZE(TMIN2),SHAPE(TMIN2)
ALLOCATE(HSNOW(NGPTOT))
HSNOW(:)=0._JPRB
IF(LLP)WRITE(IU,9) 'HSNOW    ',SIZE(HSNOW),SHAPE(HSNOW)
ALLOCATE(TSUR2(NGPTOT))
TSUR2(:)=0._JPRB
IF(LLP)WRITE(IU,9) 'TSUR2    ',SIZE(TSUR2),SHAPE(TSUR2)
ALLOCATE(TSTS(NGPTOT))
TSTS(:)=0._JPRB
IF(LLP)WRITE(IU,9) 'TSTS     ',SIZE(TSTS),SHAPE(TSTS)
ALLOCATE(TSFL(NGPTOT))
TSFL(:)=0._JPRB
IF(LLP)WRITE(IU,9) 'TSFL     ',SIZE(TSFL),SHAPE(TSFL)
ALLOCATE(XPUQ(NGPTOT))
XPUQ(:)=0._JPRB
IF(LLP)WRITE(IU,9) 'XPUQ     ',SIZE(XPUQ),SHAPE(XPUQ)
ALLOCATE(XPVQ(NGPTOT))
XPVQ(:)=0._JPRB
IF(LLP)WRITE(IU,9) 'XPVQ     ',SIZE(XPVQ),SHAPE(XPVQ)
ALLOCATE(FCRFTH(NGPTOT,0:1))
FCRFTH(:,:)=0._JPRB
IF(LLP)WRITE(IU,9) 'FCRFTH   ',SIZE(FCRFTH),SHAPE(FCRFTH)
ALLOCATE(FCRFSO(NGPTOT,0:1))
FCRFSO(:,:)=0._JPRB
IF(LLP)WRITE(IU,9) 'FCRFSO   ',SIZE(FCRFSO),SHAPE(FCRFSO)

9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUALGCO',1,ZHOOK_HANDLE)
END SUBROUTINE SUALGCO
