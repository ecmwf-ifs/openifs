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

#ifdef RS6K
@PROCESS NOCHECK
#endif

SUBROUTINE TRACMF(YDGEOMETRY,YDGMV,YGFL,PGMVS,PGMVT1S,PGFL,PGFLT1)
!
!     Purpose.   Correct global mass of gfl of tracer after SL Advection 
!     --------   
!
!*   Interface.
!     ----------
!
!        *CALL* *TRACMF (PGMVS,PGMVT1S, PGFLT, PGFLT1)
!
!     Explicit arguments :
!        --------------------
!
!
!!     INPUT:
!     -------------
!        PGMVS       : surface GMV variables at t and t-dt.
!        PGMVS       : surface GMV variables at t +dt 
!        PGFL       : GFL variables buffer t0 t1
!
!     OUTPUT:
!     -------
!        PGFLT1       : GFL variables buffer t + dt
!
!        Implicit arguments :  None.
!        --------------------
!
!     Method.
!     -------
!     - proportional mass fix 
!     - additive mass fix according to different diagnistics 
!        horizontal laplace, vertical laplace and absolute gradient
!     Externals.   See includes below.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Johannes Flemming   *ECMWF*

! Modifications
! -------------
!   M. Diamantakis  01-12-2012 LTRCMFIX_PS moved to module YOM_YGFL
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   M. Diamantakis (Feb 2015) proportional fixer move to qmfixer subroutine
!                             to allow combination with BC fixer
!--------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOM_YGFL     , ONLY : TYPE_GFLD
USE YOMCST       , ONLY : RG
USE MPL_MODULE   , ONLY : MPL_ALLREDUCE


IMPLICIT NONE

! arguments

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TYPE_GFLD)   ,INTENT(INOUT) :: YGFL
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS,YDGEOMETRY%YRDIM%NGPBLKS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMVT1S(YDGEOMETRY%YRDIM%NPROMA,YDGMV%YT1%NDIMS,YDGEOMETRY%YRDIM%NGPBLKS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YGFL%NDIM,YDGEOMETRY%YRDIM%NGPBLKS) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YGFL%NDIM1,YDGEOMETRY%YRDIM%NGPBLKS) 
! Local variables
INTEGER(KIND=JPIM) :: IST, IEND, JLEV , JGFL
REAL(KIND=JPRB) :: ZPRE9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZPRE1(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZSPT1(YDGEOMETRY%YRDIM%NPROMA),ZSPT9(YDGEOMETRY%YRDIM%NPROMA)

REAL(KIND=JPRB) :: ZFIELD3DT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZFIELD3DT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZFIELD2DT1(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZFIELD2DT0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZFIELD1D(YDGEOMETRY%YRDIM%NPROMA,2,YDGEOMETRY%YRDIM%NGPBLKS), ZNORMS(3,2)
REAL(KIND=JPRB) :: ZBETA(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIM%NGPBLKS), ZBETA_SUM(1) 

INTEGER(KIND=JPIM) :: IBL,ICEND,JKGLO,JROF

INTEGER(KIND=JPIM) :: JLEV1, JLEV2, JL, JSP
REAL(KIND=JPRB)    :: ZMASSD, ZMASSCOR
REAL(KIND=JPRB)    :: ZPRES1(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)    :: ZPRES9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIM%NGPBLKS) 
REAL(KIND=JPRB)    :: ZSP(YDGEOMETRY%YRDIMV%NFLSUR,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB)    :: ZGP(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    :: ZPRE1DT(YDGEOMETRY%YRDIM%NPROMA,2,YDGEOMETRY%YRDIM%NGPBLKS), ZPRENORMS(3,2), ZPSCOR


REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "gpnorm1.intfb.h"
#include "gppwc.intfb.h"
#include "gphpre.intfb.h"

#include "speree.intfb.h"
#include "reespe.intfb.h"

IF (LHOOK) CALL DR_HOOK('TRACMF',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDLAP=>YDGEOMETRY%YRLAP,  &
 & YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, YDSTA=>YDGEOMETRY%YRSTA,  &
 & YDCSGLEG=>YDGEOMETRY%YRCSGLEG,YDCVER =>YDGEOMETRY%YRCVER, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,  &
 & YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(LTRCMFA_DIF=>YGFL%LTRCMFA_DIF, LTRCMFA_LAP=>YGFL%LTRCMFA_LAP, &
 & LTRCMFA_VER=>YGFL%LTRCMFA_VER, LTRCMFIX_PS=>YGFL%LTRCMFIX_PS, &
 & NDIM=>YGFL%NDIM, NDIM1=>YGFL%NDIM1, NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, NSPEC2=>YDDIM%NSPEC2, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLSUR=>YDDIMV%NFLSUR, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & NDIMGMVS=>YDGMV%NDIMGMVS, YPH9=>YDGMV%YPH9, YT1=>YDGMV%YT1, &
 & NVALUE=>YDLAP%NVALUE, RLAPDI=>YDLAP%RLAPDI)
!! mass fixer type
!LTRCMFA_LAP=.FALSE.  ! additive based on hor. Laplace
!LTRCMFA_VER=.FALSE.  ! additive based on vert. Lapace
!LTRCMFA_DIF=.FALSE.   ! additiva based on abs gradient

! calculate pressure on interface for t & t + 1
DO JKGLO=1,NGPTOT,NPROMA
     ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
     IBL=(JKGLO-1)/NPROMA+1
     IST=1
     IEND=ICEND
! caluclation of the pressures  
  DO JROF=IST,IEND
        ZSPT9(JROF)=PGMVS(JROF,YPH9%MSP,IBL)
        ZPRE9 (JROF,NFLEVG) = EXP(ZSPT9(JROF))
        ZSPT1(JROF)=PGMVT1S(JROF,YT1%MSP,IBL)
        ZPRE1 (JROF,NFLEVG) = EXP(ZSPT1(JROF))
        ZPRE1DT(JROF,1,IBL)=ZPRE9(JROF,NFLEVG)
        ZPRE1DT(JROF,2,IBL)=ZPRE1(JROF,NFLEVG)
   ENDDO

   CALL GPHPRE(NPROMA,NFLEVG,IST,IEND,YDVAB,YDCVER,ZPRE9)
   CALL GPHPRE(NPROMA,NFLEVG,IST,IEND,YDVAB,YDCVER,ZPRE1)

   DO JLEV=0,NFLEVG
     DO JROF=IST,IEND
       ZPRES9(JROF,JLEV,IBL) = ZPRE9(JROF,JLEV)
       ZPRES1(JROF,JLEV,IBL) = ZPRE1(JROF,JLEV)
     ENDDO
  ENDDO 
ENDDO

! change surface pressure for t+1 to match integral  
IF (LTRCMFIX_PS ) THEN 
! mean surface pressure global fixer
! is it better to use the mass of the t=0 or from the start of the forecast
   CALL GPNORM1(YDGEOMETRY,ZPRE1DT,2,.FALSE.,PNORMS=ZPRENORMS)
   ZPSCOR=ZPRENORMS(1,1)/ZPRENORMS(1,2)

! re-calculate pressure for mass calculation at t + 1 
  DO JKGLO=1,NGPTOT,NPROMA
     ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
     IBL=(JKGLO-1)/NPROMA+1
     IST=1
     IEND=ICEND
      DO JROF=IST,IEND
        ZSPT1(JROF)=PGMVT1S(JROF,YT1%MSP,IBL)  
        ZPRE1 (JROF,NFLEVG) = EXP(ZSPT1(JROF)) * ZPSCOR
      ENDDO

      CALL GPHPRE(NPROMA,NFLEVG,IST,IEND,YDVAB,YDCVER,ZPRE1)
      DO JLEV=0,NFLEVG
        DO JROF=IST,IEND
          ZPRES1(JROF,JLEV,IBL) = ZPRE1(JROF,JLEV)
        ENDDO
      ENDDO
   ENDDO
ENDIF
! end pressure fix         

! FIELD3DT0 AND FIELD3DT1 ARE RESP. THE FIELDS BEFORE AND AFTER ADVECTION
DO JGFL=1,NUMFLDS

! apply correction if flag is set 
  IF ( YCOMP(JGFL)%LMASSFIX .AND. YCOMP(JGFL)%LADV ) THEN
! calculate total mass befoe and after advection
    DO JKGLO=1,NGPTOT,NPROMA
      ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
      IBL=(JKGLO-1)/NPROMA+1
      IST=1
      IEND=ICEND
      DO JLEV=1,NFLEVG
        DO JROF=IST,IEND
          ZFIELD3DT0(JROF,JLEV)=PGFL(JROF,JLEV,YCOMP(JGFL)%MP,IBL)
          ZFIELD3DT1(JROF,JLEV)=PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1,IBL)
          ZPRE9 (JROF,JLEV) =  ZPRES9(JROF,JLEV,IBL) 
          ZPRE1 (JROF,JLEV) =  ZPRES1(JROF,JLEV,IBL)
        ENDDO
      ENDDO

      CALL GPPWC(NPROMA,IST,IEND,NFLEVG,ZFIELD2DT0,ZFIELD3DT0,ZPRE9)
      CALL GPPWC(NPROMA,IST,IEND,NFLEVG,ZFIELD2DT1,ZFIELD3DT1,ZPRE1)

      DO JROF=IST,IEND
        ZFIELD1D(JROF,1,IBL)=ZFIELD2DT0(JROF,NFLEVG)
        ZFIELD1D(JROF,2,IBL)=ZFIELD2DT1(JROF,NFLEVG)
      ENDDO
   ENDDO 

! norms are broadcast already if PNORMS=ZNORMS specified
     CALL GPNORM1(YDGEOMETRY,ZFIELD1D,2,.FALSE.,PNORMS=ZNORMS)
    
     ZMASSCOR = (ZNORMS(1,1)-ZNORMS(1,2))  
    
! Calculate corrected Concentrations PGFLT1    

!    local additive mass fixer based on error measure ZBETA  
!   ( LTRCMFA_DIF .OR. LTRCMFA_LAP .OR. LTRCMFA_VER ) 

! calculation of weight zbeta for the additive MF
     ZBETA_SUM(1)=0.0_JPRB

! MF based on square of laplace
    IF ( LTRCMFA_LAP ) THEN
      DO JKGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
        IBL=(JKGLO-1)/NPROMA+1
        IST=1
        IEND=ICEND
        DO JLEV=1,NFLEVG
          DO JROF=IST,IEND
            ZGP(JROF+JKGLO-1,JLEV) = PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1,IBL)
          ENDDO
        ENDDO
      ENDDO
! calculate horizontal laplace
! for efficency this should be done for all species at once i.e. for NFLEVG*NCHEM
      CALL REESPE(YDGEOMETRY,NFLSUR,NFLEVG,ZSP,ZGP)
! calculate horizontal laplacian in spectral space
      DO JSP=1,NSPEC2
        DO JL=1,NFLSUR
          ZSP(JL,JSP) = RLAPDI( NVALUE(JSP) ) *  ZSP(JL,JSP) 
        ENDDO
      ENDDO
! transfer in gp space
      CALL SPEREE(YDGEOMETRY,NFLSUR,NFLEVG,ZSP,ZGP)

! calculate beta - without weight  
      DO JKGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
        IBL=(JKGLO-1)/NPROMA+1
        IST=1
        IEND=ICEND
        DO JLEV=1,NFLEVG
          DO JROF=IST,IEND
            ZBETA(JROF,JLEV,IBL) = ZGP(JROF+JKGLO-1,JLEV) * ZGP(JROF+JKGLO-1,JLEV) 
          ENDDO
        ENDDO
      ENDDO
! zsed on abs(gradient)
    ELSEIF ( LTRCMFA_DIF ) THEN     
      DO JKGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
        IBL=(JKGLO-1)/NPROMA+1
        IST=1
        IEND=ICEND
! set error measure ZBETA 
        DO JLEV=1,NFLEVG
          DO JROF=IST,IEND
! define ZBETA: difference between t and t+1 is used as a measure for "error"
            ZBETA(JROF,JLEV,IBL) = ABS(PGFL(JROF,JLEV,YCOMP(JGFL)%MP,IBL) - PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1,IBL))
          ENDDO
        ENDDO
      ENDDO
! vertical 
    ELSEIF ( LTRCMFA_VER ) THEN     
      DO JKGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
        IBL=(JKGLO-1)/NPROMA+1
        IST=1
        IEND=ICEND
! set error measure ZBETA 
        DO JLEV=1,NFLEVG
          JLEV1 = MAX(JLEV-1,1)
          JLEV2 = MIN(JLEV+1,NFLEVG)
          DO JROF=IST,IEND
! define ZBETA: vertical laplacian square
            ZBETA(JROF,JLEV,IBL) =  ( PGFLT1(JROF,JLEV1,YCOMP(JGFL)%MP1,IBL) - 2.0_JPRB * PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1,IBL)&
              &      +  PGFLT1(JROF,JLEV2,YCOMP(JGFL)%MP1,IBL)) /&
              &      ( (ZPRES9(JROF,JLEV-1,IBL) -  ZPRES9(JROF,JLEV,IBL)) * (ZPRES9(JROF,JLEV-1,IBL) -  ZPRES9(JROF,JLEV,IBL))) 
            ZBETA(JROF,JLEV,IBL) = ZBETA(JROF,JLEV,IBL) * ZBETA(JROF,JLEV,IBL) 
          ENDDO
        ENDDO
      ENDDO
    ENDIF

! scale zbeta with grid box volume and sum up for norm
    DO JKGLO=1,NGPTOT,NPROMA
      ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
      IBL=(JKGLO-1)/NPROMA+1
      IST=1
      IEND=ICEND
! set error measure ZBETA 
      DO JLEV=1,NFLEVG
        DO JROF=IST,IEND
! account for gird box size
          ZBETA(JROF,JLEV,IBL) = ZBETA(JROF,JLEV,IBL)  * ABS (  YDGSGEOM_NB%GAW(JKGLO-1+JROF)&
            &              * (ZPRES9(JROF,JLEV,IBL) -  ZPRES9(JROF,JLEV-1,IBL)) ) 
          ZBETA_SUM(1) = ZBETA_SUM(1) +  ZBETA(JROF,JLEV,IBL) 
        ENDDO
      ENDDO
    ENDDO
! gather global sum  
!    this might be costly , try to have two loops over NCHEM and stor array of ZBETA_SUM(NCHEM)   
    CALL MPL_ALLREDUCE( ZBETA_SUM,'sum',LDREPROD=.FALSE.,CDSTRING='TRACMF')    
    
    IF (  ZBETA_SUM(1) > 0.0_JPRB ) THEN
      
      DO JKGLO=1,NGPTOT,NPROMA
        ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
        IBL=(JKGLO-1)/NPROMA+1
        IST=1
        IEND=ICEND
        DO JLEV=1,NFLEVG
          DO JROF=IST,IEND
               ! air mass density in kg/m2 
            ZMASSD = (( ZPRES9(JROF,JLEV,IBL) - ZPRES9(JROF,JLEV-1,IBL)) / RG ) *  YDGSGEOM_NB%GAW(JKGLO-1+JROF)  
            PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1,IBL) = PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1,IBL) +&
              &                           ( ZBETA(JROF,JLEV,IBL) / ZBETA_SUM(1)) * ZMASSCOR / ZMASSD
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  ENDIF ! LMASSFIX
 
ENDDO ! loop JGFL

!* END PART FOR TRACER MASS FIXER

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TRACMF',1,ZHOOK_HANDLE)
END SUBROUTINE TRACMF
