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

SUBROUTINE CORMASSDRY(YDGEOMETRY,YGFL,YDDYN,YDSPEC,PGP,KSTEP,LDGPMASCOR)

!**** *CORMASSDRY*  Cheap mass corrector 'MASCOR'.
!                   - ROUTINE FOR CALCULATING GLOBAL MASS
!                     IF REQUIRED (LMASDRY=T) SUBTRACT TOTAL WATER TOTAL MASS

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        *CALL* *CORMASSDRY(...)*

!        INPUT:

!        KRET      - field number whose norm is to be computed
!        PGP       - Input GFL array

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------

!        Nils Wedi + Mats Hamrud *ECMWF*

!     MODIFICATIONS.
!     --------------
!       Original : 08-02-2008
!       K. Yessad (Sep 2008): use array VFE_RDETAH.
!       G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!       M.Diamantakis(Feb 2014): Fix gridpoint version of mass fixer
!       T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!       J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!       S. Malardel Compute dry mass from all active water GFLs
!     -----------------------------------------------------------------

USE YOM_YGFL     , ONLY : TYPE_GFLD
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYN   , ONLY : TDYN
USE INTDYN_MOD,ONLY : YYTXYB
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)

!     -----------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)          ,INTENT(INOUT) :: YDDYN
TYPE(TYPE_GFLD)     ,INTENT(INOUT) :: YGFL
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSPEC
REAL(KIND=JPRB)     ,INTENT(IN)    :: PGP(:,:,:,:) 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KSTEP
LOGICAL             ,INTENT(IN)    :: LDGPMASCOR

!     -----------------------------------------------------------------

REAL(KIND=JPRB) :: ZWORK(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZWORK1(YDGEOMETRY%YRDIM%NPROMA,1,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB) :: ZIN(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZOUT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZNORMS(3,1)
REAL(KIND=JPRB) :: ZXYB(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)

REAL(KIND=JPRB) :: ZPRESH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZGPW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IST, IEND, JF, JROF, JKGLO, IBL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -----------------------------------------------------------------

#include "gpnorm1.intfb.h"
#include "speree.intfb.h"
#include "gphpre.intfb.h"
#include "verdisint.intfb.h"
#include "gpmasscor.intfb.h"

!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CORMASSDRY',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM,&
 &  YDMP=>YDGEOMETRY%YRMP,YDCVER=>YDGEOMETRY%YRCVER,  &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDVAB=>YDGEOMETRY%YRVAB,  &
 & YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, YDSTA=>YDGEOMETRY%YRSTA, &
 & YDLAP=>YDGEOMETRY%YRLAP, YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,  &
  & YDSPGEOM=>YDGEOMETRY%YSPGEOM,YQ=>YGFL%YQ,YL=>YGFL%YL,YI=>YGFL%YI,YR=>YGFL%YR,YS=>YGFL%YS)
ASSOCIATE(NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & GMASS0=>YDDYN%GMASS0, GMASSI=>YDDYN%GMASSI, GMASSINC=>YDDYN%GMASSINC, &
 & LMASDRY=>YDDYN%LMASDRY, NGPMASCOR=>YDDYN%NGPMASCOR, &
 & NGPTOT=>YDGEM%NGPTOT)
!     -----------------------------------------------------------------

!--- LOOP OVER NPROMA PACKETS

! convert surface pressure to gridpoint, store in zwork
CALL SPEREE(YDGEOMETRY,1,1,YDSPEC%SP,ZWORK)

CALL GSTATS(1223,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP& PRIVATE(JKGLO,IST,IEND,JF,JROF,IBL) &
!$OMP& PRIVATE(ZPRESH,ZGPW) &
!$OMP& PRIVATE(ZIN,ZOUT,ZXYB)
DO JKGLO=1,NGPTOT,NPROMA
  IST=1
  IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
  IBL=(JKGLO-1)/NPROMA+1
  
  DO JROF = IST , IEND
    ZPRESH(JROF,NFLEVG) = EXP(ZWORK(JKGLO-1+JROF))
  ENDDO
  ZGPW(:,:)=0.0_JPRB
  IF( LMASDRY .AND. YQ%LACTIVE) THEN
! add at least water vapor in total water
    DO JF=1,NFLEVG
      DO JROF = IST,IEND
        ZGPW(JROF,JF)=PGP(JROF,JF,YQ%MP,IBL) 
      ENDDO
    ENDDO
! add condensates only if they are all active (usual case)
    IF(YL%LACTIVE .AND. YI%LACTIVE .AND. YR%LACTIVE .AND. YS%LACTIVE ) THEN
    ! compute total water
      DO JF=1,NFLEVG
        DO JROF = IST,IEND
          ZGPW(JROF,JF)= ZGPW(JROF,JF) &
 &         + PGP(JROF,JF,YL%MP,IBL) + PGP(JROF,JF,YI%MP,IBL) &
 &         + PGP(JROF,JF,YR%MP,IBL) + PGP(JROF,JF,YS%MP,IBL)
        ENDDO
      ENDDO       
   ENDIF
    ! calculate delta p
    CALL GPHPRE(NPROMA,NFLEVG,IST,IEND,YDVAB,YDCVER,ZPRESH,PXYB=ZXYB)

    ! integrate total water mass at each point
    IF(YDCVER%LVERTFE) THEN
      DO JF=1,NFLEVG
        DO JROF = IST,IEND
          ZIN(JROF,JF) = ZGPW(JROF,JF)&
           & *ZXYB(JROF,JF,YYTXYB%M_DELP)*YDVETA%VFE_RDETAH(JF)
        ENDDO
      ENDDO
      ZIN(IST:IEND,0)=0.0_JPRB
      ZIN(IST:IEND,NFLEVG+1)=0.0_JPRB
      CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',NPROMA,IST,IEND,NFLEVG,ZIN,ZOUT)
    ELSE
      ZOUT(IST:IEND,1)=0.0_JPRB
      DO JF=2,NFLEVG+1
        DO JROF = IST,IEND
          ZOUT(JROF,JF)=ZOUT(JROF,JF-1)&
           & +ZGPW(JROF,JF-1)*ZXYB(JROF,JF-1,YYTXYB%M_DELP)
        ENDDO
      ENDDO
    ENDIF    
    ! store current dry mass in zwork (apply map factor if necessary)
    DO JROF=IST,IEND
      ZWORK1(JROF,1,IBL) = (ZPRESH(JROF,NFLEVG) - ZOUT(JROF,NFLEVG+1))/YDGSGEOM_NB%GM(JKGLO-1+JROF)**2
    ENDDO

  ELSE    
    ! Note Mats Hamrud + Nils Wedi:
    ! Since the exp(integral ln(ps) domega) <> integral(ps domega)
    ! the exact mean conservation of ps does not allow 
    ! (at least without incurring a small error) to simply 
    ! maintain the mean ln(ps) as extracted from the 0-component in spectral space
    DO JROF=IST,IEND
      ZWORK1(JROF,1,IBL) = ZPRESH(JROF,NFLEVG)/YDGSGEOM_NB%GM(JKGLO-1+JROF)**2
    ENDDO
  ENDIF
  
ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1223,1)
  
!--- Now we have all the points in an NGPTOT data structure form
!    and we can calculate grid point norms

CALL GPNORM1(YDGEOMETRY,ZWORK1,1,.FALSE.,PNORMS=ZNORMS)
GMASS0=ZNORMS(1,1)
IF( KSTEP == 0 .AND. GMASSI == 0.0_JPRB ) THEN
  GMASSI = GMASS0
ENDIF
GMASSINC=LOG(GMASS0/GMASSI)  

! possibly spmascor correction is also correct for stretched grid, Nils
!-----------------------------------------------------------------------
! md: apply massfixer in grid-point space instead of spectral
! it costs more but is more accurate (eliminates completely mass error)  
!-----------------------------------------------------------------------
IF (LDGPMASCOR) CALL GPMASSCOR(YDGEOMETRY,YGFL,YDDYN,NGPMASCOR,ZWORK,YDSPEC,GMASS0,GMASSI,GMASSINC)

!     -----------------------------------------------------------------
  
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CORMASSDRY',1,ZHOOK_HANDLE)
END SUBROUTINE CORMASSDRY
