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

SUBROUTINE GPMASSCOR(YDGEOMETRY,YGFL,YDDYN,KGPMASCOR,PGPF,YDSP,PGMASS0,PGMASSI,PGMASSINC)
!
!     Purpose.   Correct pressure field to satisfy global mass conservation. 
!
!     --------   
!
!*   Interface.
!     ----------
!
!        *CALL* *GPMASSCOR(NGPMASCOR,PGPF,PSPF,PGMASS0,PGMASSI,PGMASSINC)
!
!
!!     INPUT:
!     -------------
!        NGPMASCOR  : mass fixer algorithm: 0 -> simple proportional
!                                           1 -> weighted correction
!        PGMASS0    : initial time total model mass (or dry air mass)
!        PGMASSI    : current time total model mass (or dry air mass)
!
!      INPUT/OUTPUT:
!     ---------------
!        PGPF       : surface grid-point pressure field at current time t
!        PSPF       : surface spectral pressure field at current time t
!        PGMASSINC  : LOG(PGMASS0/PGMASSI)
!
!        Implicit arguments :  None.
!        --------------------
!
!     Method.
!     -------
!     - NGPMASCOR: 0 -> proportional mass fixer 
!     - NGPMASCOR: 1 -> Zerroukat mass fixer algorithm variant (weight proportional to Laplacian)
!
!     Externals.   See includes below.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!     Michail Diamantakis   *ECMWF*

!
! End Modifications
!------------------------------------------------------------------------------

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM , JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN             , ONLY : NULOUT
USE YOM_YGFL           , ONLY : TYPE_GFLD
USE YOMDYN             , ONLY : TDYN
USE YOMSTA             , ONLY : RTSUR
USE YOMCST             , ONLY : RD
USE YOMMP0             , ONLY : NPROC
USE MPL_MODULE         , ONLY : MPL_ALLREDUCE
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)          ,INTENT(IN)    :: YDDYN
TYPE(TYPE_GFLD)     ,INTENT(IN)    :: YGFL
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KGPMASCOR
REAL(KIND=JPRB)     ,INTENT(IN)    :: PGMASS0
REAL(KIND=JPRB)     ,INTENT(IN)    :: PGMASSI
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PGMASSINC
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PGPF(YDGEOMETRY%YRGEM%NGPTOT)
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSP
!------------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IST,  IEND, IBL, JSP, JKGLO, JROF
REAL(KIND=JPRB) :: ZWORK(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB) :: ZDM, ZMAXD, ZRMSDIF, ZRMSFLD, ZDIFF, ZDMRATIO, ZMAXDIF
REAL(KIND=JPRB) :: ZGLBNRM(3,1), ZLAMBDA, ZPS0, ZPS1
REAL(KIND=JPRB) :: ZSP(YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) :: ZGP(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZCMSLP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "gpnorm1.intfb.h"
#include "speree.intfb.h"
#include "reespe.intfb.h"

!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPMASSCOR',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  &
 & YDLAP=>YDGEOMETRY%YRLAP, YDSTA=>YDGEOMETRY%YRSTA, YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA,  &
 & YDVFE=>YDGEOMETRY%YRVFE, YDCSGLEG=>YDGEOMETRY%YRCSGLEG,  &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, &
  & YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB, YDGSGEOM=>YDGEOMETRY%YRGSGEOM, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB,  &
  & YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NMFDIAGLEV=>YGFL%NMFDIAGLEV, &
 & NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, NSPEC2=>YDDIM%NSPEC2, &
 & RCMSLP0=>YDDYN%RCMSLP0, &
 & NGPTOT=>YDGEM%NGPTOT, NGPTOTG=>YDGEM%NGPTOTG, &
 & NVALUE=>YDLAP%NVALUE, RLAPDI=>YDLAP%RLAPDI)
!------------------------------------------------------------------------------

IST=1

IF (KGPMASCOR>0) THEN
  ZCMSLP=RCMSLP0/(RD*RTSUR)
!--------------------------------------------------------------
! Compute Laplacian of smoothened (orography extracted) ln(Ps)
!--------------------------------------------------------------
  DO JSP=1,NSPEC2
    ZSP(JSP)=RLAPDI(NVALUE(JSP))*(ZCMSLP*YDSP%OROG(JSP)+YDSP%SP(JSP))
  ENDDO
  CALL SPEREE(YDGEOMETRY,1,1,ZSP,ZGP) ! to gp space
!----------------------------------------------------------------------------
! calculate weight 
!----------------------------------------------------------------------------
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP& PRIVATE(JKGLO,IEND,IBL,JROF) 
  DO JKGLO=1,NGPTOT,NPROMA
    IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/NPROMA+1
    DO JROF=IST,IEND
      ZWORK(JROF,IBL) = SQRT(ABS(ZGP(JROF+JKGLO-1)))
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  CALL GPNORM1(YDGEOMETRY,ZWORK,1,.FALSE.,PNORMS=ZGLBNRM)
!----------------------------------------------------------------------------
! Correct ln(surface pressure)
!----------------------------------------------------------------------------
  ZDM=PGMASS0-PGMASSI
  ZLAMBDA=ZDM/ZGLBNRM(1,1)
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP& PRIVATE(JKGLO,IEND,IBL,JROF)
  DO JKGLO=1,NGPTOT,NPROMA
    IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/NPROMA+1
    DO JROF=IST,IEND
      ZPS0=EXP(PGPF(JROF+JKGLO-1))
      PGPF(JROF+JKGLO-1)=LOG(ZPS0-ZLAMBDA*ZWORK(JROF,IBL))
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ELSE
! apply simple proportional mass fixer
  DO JKGLO=1,NGPTOT
    PGPF(JKGLO)=PGPF(JKGLO)-PGMASSINC
  ENDDO
ENDIF

! back to spectral space
CALL REESPE(YDGEOMETRY,1,1,YDSP%SP,PGPF)
! Reset log(incr) to 0 to prevent modifying again in spcmascor()
PGMASSINC=0.0_JPRB

!----------------------------------------------------------------------------
! Compute surface grid-point pressure increment diagnostics
!----------------------------------------------------------------------------
IF (NMFDIAGLEV>0) THEN
  ZMAXD=0.0_JPRB
  ZRMSDIF=0.0_JPRB
  ZRMSFLD=0.0_JPRB
  IF (KGPMASCOR>0) THEN
    DO JKGLO=1,NGPTOT,NPROMA
      IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
      IBL=(JKGLO-1)/NPROMA+1
      DO JROF=IST,IEND
        ZDIFF=-ZLAMBDA*ZWORK(JROF,IBL)
        ZRMSDIF=ZRMSDIF+ZDIFF*ZDIFF
        ZPS0=EXP(PGPF(JROF+JKGLO-1))-ZDIFF
        ZRMSFLD=ZRMSFLD+ZPS0*ZPS0
        ZMAXD=MAX(ABS(ZDIFF),ZMAXD)
      ENDDO
    ENDDO
  ELSE
    DO JKGLO=1,NGPTOT,NPROMA
      IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
      IBL=(JKGLO-1)/NPROMA+1
      DO JROF=IST,IEND
        ZPS1=EXP(PGPF(JROF+JKGLO-1))
        ZPS0=ZPS1*(PGMASS0/PGMASSI)
        ZDIFF=ZPS1-ZPS0
        ZRMSDIF=ZRMSDIF+ZDIFF*ZDIFF
        ZRMSFLD=ZRMSFLD+ZPS0*ZPS0
        ZMAXD=MAX(ABS(ZDIFF),ZMAXD)
      ENDDO
    ENDDO
  ENDIF
  IF (NPROC>1) THEN
    CALL MPL_ALLREDUCE(ZRMSDIF,'SUM',LDREPROD=.FALSE.,CDSTRING='GPMASSCOR:')
    CALL MPL_ALLREDUCE(ZRMSFLD,'SUM',LDREPROD=.FALSE.,CDSTRING='GPMASSCOR:')
    CALL MPL_ALLREDUCE(ZMAXD,'MAX',LDREPROD=.FALSE.,CDSTRING='GPMASSCOR:')
  ENDIF
  ZRMSDIF=SQRT(ZRMSDIF/NGPTOTG)
  ZRMSFLD=SQRT(ZRMSFLD/NGPTOTG)
! write header for fixer output: mass imbalance, imbalance ratio to total mass x 100, max fixer correction
! over RMS norm of corrected field x 100, ratio of rms norm of correction field to norm of field x 100
  WRITE(NULOUT,*) '----------------------------------------------------------------------------'
  WRITE(NULOUT,*) '             GPMASSCOR PRESSURE FIXER GLOBAL DIAGNOSTICS                    '
  WRITE(NULOUT,*) '    FIELD         DM       DM/Mtot %    max(dij)/|field| %   |dij|/|field| %'
  WRITE(NULOUT,*) '----------------------------------------------------------------------------'
  ZDM=PGMASS0-PGMASSI
  ZDMRATIO=100.0_JPRB*ZDM/PGMASSI
  ZMAXDIF=100.0_JPRB*ZMAXD/ZRMSFLD
  ZRMSDIF=100.0_JPRB*ZRMSDIF/ZRMSFLD
  WRITE(NULOUT,'(A13,2X,E10.4,1X,E9.2,6X,E10.2,11X,E9.2)') 'SURF PRESS',&
        &           ZDM,ZDMRATIO,ZMAXDIF,ZRMSDIF
ENDIF

!------------------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPMASSCOR',1,ZHOOK_HANDLE)
END SUBROUTINE GPMASSCOR
