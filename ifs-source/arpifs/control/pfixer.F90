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

SUBROUTINE PFIXER(YDGEOMETRY,YDGMV,YGFL,YDDYN,PSURF0,PSURF1,PGMVT1S)
!
!     Purpose.   Correct pressure field to satisfy global conservation. Scheme based
!                on McGregor's algorithm applied to C-CAM climate model or simple 
!                proportional. Current use intended mainly for runs with CO2, CH4 species.
!                

!     --------   
!
!*   Interface.
!     ----------
!
!        *CALL* *PFIXER(PGMVS,PSURF0,PSURF1,PGMVT1S)
!
!
!!     INPUT:
!     ----------
!        PGMVS       : surface fields at t and t-dt.
!        PSURF0     : surface pressure at time t
!
!     INPUT/OUTPUT:
!     -------------
!        PSURF1     : surface pressure at current time (t+dt)
!                     This array is local to the fixer (other parts of
!                     the model do not see it).
!        PGMVT1S    : ln(PSURF1). May or may not be altered by
!                     the fixer depending on flag LGPMASCOR. This array
!                     is passed to other parts of dynamics/physics.

!        Implicit arguments :  None.
!        --------------------
!
!     Method.
!     -------
!     - McGregor's mass fixer algorithm applied to pressure field
!
!     Externals.   See includes below.
!     ----------

!     Reference.
!     ----------
!       ECMWF Tech Memo 713
!       Diamantakis, J. Flemming : Global mass fixer algorithms for 
!                                  conservative tracer transport in the
!                                  ECMWF model

!     Author.
!     -------
!     Michail Diamantakis   *ECMWF*

!
! End Modifications
!------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM , JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN       , ONLY : NULOUT
USE YOM_YGFL     , ONLY : TYPE_GFLD
USE YOMDYN       , ONLY : TDYN
USE YOMMP0       , ONLY : NPROC
USE MPL_MODULE   , ONLY : MPL_ALLREDUCE

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TDYN)        ,INTENT(INOUT) :: YDDYN
TYPE(TYPE_GFLD)   ,INTENT(INOUT) :: YGFL
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSURF0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSURF1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1S(YDGEOMETRY%YRDIM%NPROMA,YDGMV%YT1%NDIMS,YDGEOMETRY%YRDIM%NGPBLKS) 
!------------------------------------------------------------------------------
INTEGER(KIND=JPIM) :: IST,  IEND, IBL, JKGLO, JROF

LOGICAL :: LLMGFIX

REAL(KIND=JPRB) :: ZPRE1DT(YDGEOMETRY%YRDIM%NPROMA,2,YDGEOMETRY%YRDIM%NGPBLKS),&
 &  ZDIFPOSNEG(YDGEOMETRY%YRDIM%NPROMA,2,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB) :: ZDM, ZMAXD, ZRMSDIF, ZRMSFLD, ZDIFF, ZDMRATIO, ZMAXDIF
REAL(KIND=JPRB) :: ZSUMPOSNEG(3,2), ZPRENORMS(3,2), ZALPHA, ZINVAL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "dist_grid.h"
#include "gath_grid.h"
#include "gpnorm1.intfb.h"

!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PFIXER',0,ZHOOK_HANDLE)
!------------------------------------------------------------------------------

! Initializations
LLMGFIX=YGFL%LTRCMFMG.OR.(YDDYN%NGPMASCOR>0)
IST=1

!----------------------------------------------------------------------------
! Proportional mass fixer
!----------------------------------------------------------------------------

IF (.NOT.LLMGFIX.OR.YGFL%NMFDIAGLEV>1) THEN
! Compute total mass before and after advection
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP& PRIVATE(JKGLO,IEND,IBL,JROF)
  DO JKGLO=1,YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIM%NPROMA
    IEND=MIN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRGEM%NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/YDGEOMETRY%YRDIM%NPROMA+1
    DO JROF=IST,IEND
      ZPRE1DT(JROF,1,IBL)=PSURF0(JROF,IBL)
      ZPRE1DT(JROF,2,IBL)=PSURF1(JROF,IBL)  
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  CALL GPNORM1(YDGEOMETRY,ZPRE1DT,2,.FALSE.,PNORMS=ZPRENORMS)
  IF (.NOT.LLMGFIX) THEN
    ZALPHA=ZPRENORMS(1,1)/ZPRENORMS(1,2)
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP& PRIVATE(JKGLO,IEND,IBL,JROF)
    DO JKGLO=1,YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIM%NPROMA
      IEND=MIN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRGEM%NGPTOT-JKGLO+1)
      IBL=(JKGLO-1)/YDGEOMETRY%YRDIM%NPROMA+1
      DO JROF=IST,IEND
        PSURF1(JROF,IBL)=ZALPHA*PSURF1(JROF,IBL)
        PGMVT1S(JROF,YDGMV%YT1%MSP,IBL)=LOG(PSURF1(JROF,IBL))
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDIF
ENDIF

!----------------------------------------------------------------------------
! Mc Gregor's mass fixer
!----------------------------------------------------------------------------

IF (LLMGFIX) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP& PRIVATE(JKGLO,IEND,IBL,JROF)
  DO JKGLO=1,YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIM%NPROMA
    IEND=MIN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRGEM%NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/YDGEOMETRY%YRDIM%NPROMA+1
    DO JROF=IST,IEND
      ZDIFPOSNEG(JROF,1,IBL)=MAX(0.0_JPRB,PSURF1(JROF,IBL)-PSURF0(JROF,IBL))
      ZDIFPOSNEG(JROF,2,IBL)=MIN(0.0_JPRB,PSURF1(JROF,IBL)-PSURF0(JROF,IBL))
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  CALL GPNORM1(YDGEOMETRY,ZDIFPOSNEG,2,.FALSE.,PNORMS=ZSUMPOSNEG)
  ZALPHA=SQRT(-ZSUMPOSNEG(1,2)/ZSUMPOSNEG(1,1))
  ZINVAL=1.0_JPRB/ZALPHA
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP& PRIVATE(JKGLO,IEND,IBL,JROF)
  DO JKGLO=1,YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIM%NPROMA
    IEND=MIN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRGEM%NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/YDGEOMETRY%YRDIM%NPROMA+1
    DO JROF=IST,IEND
      PSURF1(JROF,IBL)=PSURF0(JROF,IBL)+ZALPHA*ZDIFPOSNEG(JROF,1,IBL)+&
        &                                   ZINVAL*ZDIFPOSNEG(JROF,2,IBL)
      PGMVT1S(JROF,YDGMV%YT1%MSP,IBL)=LOG(PSURF1(JROF,IBL))
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ENDIF

!----------------------------------------------------------------------------
! Compute mass fixer diagnostics
!----------------------------------------------------------------------------
IF (YGFL%NMFDIAGLEV>1) THEN
  ZMAXD=0.0_JPRB
  ZRMSDIF=0.0_JPRB
  ZRMSFLD=0.0_JPRB
  DO JKGLO=1,YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIM%NPROMA
    IEND=MIN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRGEM%NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/YDGEOMETRY%YRDIM%NPROMA+1
  ! calculation of surface pressures  
    DO JROF=IST,IEND
      ZDIFF=ZPRE1DT(JROF,2,IBL)-PSURF1(JROF,IBL)
      ZRMSDIF=ZRMSDIF+ZDIFF*ZDIFF
      ZRMSFLD=ZRMSFLD+ZPRE1DT(JROF,2,IBL)*ZPRE1DT(JROF,2,IBL)
      ZMAXD=MAX(ABS(ZDIFF),ZMAXD)
    ENDDO
  ENDDO
  IF (NPROC>1) THEN
    CALL MPL_ALLREDUCE(ZRMSDIF,'SUM',LDREPROD=.FALSE.,CDSTRING='PFIXER:')
    CALL MPL_ALLREDUCE(ZRMSFLD,'SUM',LDREPROD=.FALSE.,CDSTRING='PFIXER:')
    CALL MPL_ALLREDUCE(ZMAXD,'MAX',LDREPROD=.FALSE.,CDSTRING='PFIXER:')
  ENDIF
  ZRMSDIF=SQRT(ZRMSDIF/YDGEOMETRY%YRGEM%NGPTOTG)
  ZRMSFLD=SQRT(ZRMSFLD/YDGEOMETRY%YRGEM%NGPTOTG)
! write header for fixer output: mass imbalance, imbalance ratio to total mass x 100, max fixer correction
! over RMS norm of corrected field x 100, ratio of rms norm of correction field to norm of field x 100
  WRITE(NULOUT,*) '----------------------------------------------------------------------------'
  WRITE(NULOUT,*) '            ADVECTION PRESSURE MASS FIXER GLOBAL DIAGNOSTICS             '
  WRITE(NULOUT,*) '    FIELD         DM       DM/Mtot %    max(dij)/|field| %   |dij|/|field| %'
  WRITE(NULOUT,*) '----------------------------------------------------------------------------'

  ZDM=ZPRENORMS(1,2)-ZPRENORMS(1,1)
  ZDMRATIO=100.0_JPRB*ZDM/ZPRENORMS(1,1)
  ZMAXDIF=100.0_JPRB*ZMAXD/ZRMSFLD
  ZRMSDIF=100.0_JPRB*ZRMSDIF/ZRMSFLD
  WRITE(NULOUT,'(A13,2X,E10.4,1X,E9.2,6X,E10.2,11X,E9.2)') 'SURF PRESS',&
        &           ZDM,ZDMRATIO,ZMAXDIF,ZRMSDIF

ENDIF

!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PFIXER',1,ZHOOK_HANDLE)
END SUBROUTINE PFIXER
