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

SUBROUTINE CORMASS2(YDGEOMETRY,YDSPEC)

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT   
USE DISGRID_MOD, ONLY : DISGRID_SEND, DISGRID_RECV
USE DIWRGRID_MOD, ONLY : DIWRGRID_SEND, DIWRGRID_RECV
USE YOMMP0   , ONLY : MYPROC, MY_REGION_NS, MY_REGION_EW
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)

!**** *CORMASS2*  - COMPUTE THE MASS CORRECTION FOR LONG CLIMATIC RUNS 

!     Purpose.
!     --------
!     This program computes the change in dry surface pressure
!     to maintain a constant global mean

!**   Interface.
!     ----------
!        *CALL* *CORMASS2

!        Explicit arguments :     NONE
!        --------------------

!        Implicit arguments :     NONE
!        --------------------

!     Method.
!     -------
!     - Go to grid point space (SPEREE) ;
!     -   Compute the total vapor water field (Pvap),
!         then the dry-air surface pressure field (Pdry = P - Pvap)
!     -   Uniformally apply the dry-air mass correction to the surface 
!         pressure field.
!     - Go back to spectral  space (REESPE) ;
!     Externals.
!     ----------

!     Reference.
!     ----------
!        CLIMAT-ARPEGE documentation (Meteo-France, CNRM/GMGEC/EAC)

!     Author.
!     -------
!        Jean Philippe Piedelievre and Michel Deque
!          (from CORMASS written by Pascal Marquet)  *Meteo-France*

!     Modifications.
!     --------------
!        Original : 05-02-04
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        A.Alias       26-Oct-2006 MYSETA/MYSETB replaced by
!                                  MY_REGION_NS/MY_REGION_EW
!        Apr 2008  K. Yessad: use DISGRID instead of DISGRID_C + cleanings
!        R. El Khatib : 23-Apr-2010 use disgrid_mod & diwrgrid_mod
!        G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM and TCSGLEG
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSPEC

REAL(KIND=JPRB):: ZSPS (YDGEOMETRY%YRDIM%NSPEC2),  ZSPQ  (YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB):: ZREE  (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB):: ZRW   (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB):: ZPSF  (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB):: ZPVF  (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB):: ZDPF  (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB):: ZPDRYF(YDGEOMETRY%YRGEM%NGPTOT)

REAL(KIND=JPRB):: ZRWG   (YDGEOMETRY%YRGEM%NGPTOTG), ZNLG   (YDGEOMETRY%YRGEM%NGPTOTG), ZGMG (YDGEOMETRY%YRGEM%NGPTOTG)
REAL(KIND=JPRB):: ZPSFG  (YDGEOMETRY%YRGEM%NGPTOTG)
REAL(KIND=JPRB):: ZPDRYFG(YDGEOMETRY%YRGEM%NGPTOTG)

INTEGER(KIND=JPIM)::&
 & IGLG,IOPROC,IFLD,JGL,JLEV,JLON,&
 & JROF, IROF  

REAL(KIND=JPRB)::   ZAPCORM, ZAPDRYF, ZAPDRYI,&
 & ZAPSF, ZDA,&
 & ZDB, ZZNEW, ZWEIG  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "reespe.intfb.h"
#include "speree.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CORMASS2',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDCSGLEG=>YDGEOMETRY%YRCSGLEG,  &
 & YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, YDSTA=>YDGEOMETRY%YRSTA,  &
 & YDLAP=>YDGEOMETRY%YRLAP, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,  &
  & YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NDGLG=>YDDIM%NDGLG, NSPEC2=>YDDIM%NSPEC2, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NGPTOT=>YDGEM%NGPTOT, NGPTOTG=>YDGEM%NGPTOTG, NLOENG=>YDGEM%NLOENG, &
 & MYLATS=>YDMP%MYLATS, NFRSTLOFF=>YDMP%NFRSTLOFF, NLSTLAT=>YDMP%NLSTLAT, &
 & NONL=>YDMP%NONL, NPTRFLOFF=>YDMP%NPTRFLOFF)
!     ------------------------------------------------------------------

!*       1.    START OF ROUTINE.
!              -----------------

!*       1.1   INITIALISISATIONS
!               - - - - - - - - 
WRITE(NULOUT,*) ' '
WRITE(NULOUT,*) '    ====================================== '
WRITE(NULOUT,*) '    =======  ENTRANCE IN CORMASS2 ======== '
WRITE(NULOUT,*) '    ====================================== '
WRITE(NULOUT,*) ' '

! IMPOSED DRY AIR MASS (PA)  
ZAPDRYI=98320.00000000_JPRB

IOPROC=1 

IFLD=0

! One proc for IO

!     ------------------------------------------------------------------
!*       3.    COMPUTE THE TOTAL MASS OF WATER VAPOR.  (see SUSPECA*)
!              ------------------------------------------------------

WRITE(NULOUT,*) '  '
WRITE(NULOUT,*) ' COMPUTE THE TOTAL WATER VAPOR MASS '

ZSPS(1:NSPEC2)=YDSPEC%SP(1:NSPEC2)
CALL SPEREE(YDGEOMETRY,1,1,ZSPS,ZPSF)
ZPSF=EXP(ZPSF)
ZPVF= 0.0_JPRB
DO JLEV=1,NFLEVG

!        *** Compute pressure difference accross layers *** 

  ZDA=YDVAB%VAH(JLEV)-YDVAB%VAH(JLEV-1)
  ZDB=YDVAB%VBH(JLEV)-YDVAB%VBH(JLEV-1)
  DO JROF=1,NGPTOT
    ZDPF(JROF)=ZDA+ZDB*ZPSF(JROF)
  ENDDO

  ZSPQ(1:NSPEC2)=YDSPEC%Q(JLEV,1:NSPEC2)
  CALL SPEREE(YDGEOMETRY,1,1,ZSPQ,ZREE)

  DO JROF=1,NGPTOT
    ZPVF(JROF) = ZPVF(JROF) + ZREE(JROF)*ZDPF(JROF)
  ENDDO

ENDDO !JLEV
DO JROF=1,NGPTOT
  ZPDRYF(JROF) = ZPSF(JROF) -  ZPVF(JROF)
ENDDO
WRITE (NULOUT,*)ZPDRYF(1),ZPSF(1)

!     ------------------------------------------------------------------
!*       5.    COMPUTE AVERAGE CORRECTION OF MASS. (see ZODIGI, ZODIA)
!      -------------------------------------------------------

! Local version of RW

IROF=0
DO JGL=1,NLSTLAT(MY_REGION_NS)-NFRSTLOFF
  IGLG=MYLATS(JGL)
  DO JLON=1,NONL(NPTRFLOFF+JGL,MY_REGION_EW)
    IROF=IROF+1
    ZRW(IROF)=YDCSGLEG%RW(IGLG)
  ENDDO
ENDDO

IF (MYPROC == IOPROC) THEN

  WRITE(NULOUT,*) ' COMPUTE AVERAGE CORRECTION OF MASS '

! Global

  ZAPSF  = 0.0_JPRB
  ZAPDRYF= 0.0_JPRB

! Receive

  IFLD=IFLD+1
  CALL DIWRGRID_RECV(YDGEOMETRY,1,ZPSF,IFLD,ZPSFG)  
  IFLD=IFLD+1
  CALL DIWRGRID_RECV(YDGEOMETRY,1,ZPDRYF,IFLD,ZPDRYFG)  
  IFLD=IFLD+1
  CALL DIWRGRID_RECV(YDGEOMETRY,1,ZRW,IFLD,ZRWG)
  IFLD=IFLD+1
  CALL DIWRGRID_RECV(YDGEOMETRY,1,YDGSGEOM_NB%GM,IFLD,ZGMG)

  IROF=0
  DO JGL=1,NDGLG
    DO JLON=1,NLOENG(JGL)
      IROF=IROF+1
      ZNLG(IROF)=REAL(NLOENG(JGL),JPRB)
    ENDDO
  ENDDO

  DO JROF=1,NGPTOTG
    ZWEIG=ZRWG(JROF)/ZGMG(JROF)**2/ZNLG(JROF)
    ZAPSF  =ZAPSF   +   ZPSFG(JROF)*ZWEIG
    ZAPDRYF=ZAPDRYF + ZPDRYFG(JROF)*ZWEIG
  ENDDO

  ZAPCORM = - ZAPDRYF + ZAPDRYI
  ZZNEW   =   ZAPSF   + ZAPCORM

  WRITE (NULOUT,*) '  '
  WRITE (NULOUT,'(1X,'' MEAN IMPO. DRY PRESSURE '',f19.8)') ZAPDRYI
  WRITE (NULOUT,'(1X,'' MEAN FINAL DRY PRESSURE '',f19.8)') ZAPDRYF
  WRITE (NULOUT,*) '  '
  WRITE (NULOUT,'(1X,'' OLD MEAN FINAL TOT PRES '',f19.8)') ZAPSF
  WRITE (NULOUT,'(1X,'' MEAN CORR. OF  DRY PRES '',f19.8)') ZAPCORM
  WRITE (NULOUT,'(1X,'' NEW MEAN FINAL TOT PRES '',f19.8)') ZZNEW
  WRITE (NULOUT,*) '  '

!     ------------------------------------------------------------------
!*       6.    MAKE THE CORRECTION OF MASS AND WRITE ON *FA* FINAL FILE.
!              ---------------------------------------------------------

  DO JROF=1,NGPTOTG
    ZPSFG(JROF)=LOG(ZPSFG(JROF)+ZAPCORM)
  ENDDO

! Send

  IFLD=IFLD+1
  CALL DISGRID_SEND(YDGEOMETRY,1,ZPSFG,IFLD,ZPSF)
  
ELSE

! Send

  IFLD=IFLD+1
  CALL DIWRGRID_SEND(YDGEOMETRY%YRGEM,IOPROC,1,ZPSF,IFLD)  
  IFLD=IFLD+1
  CALL DIWRGRID_SEND(YDGEOMETRY%YRGEM,IOPROC,1,ZPDRYF,IFLD)
  IFLD=IFLD+1
  CALL DIWRGRID_SEND(YDGEOMETRY%YRGEM,IOPROC,1,ZRW,IFLD)
  IFLD=IFLD+1
  CALL DIWRGRID_SEND(YDGEOMETRY%YRGEM,IOPROC,1,YDGSGEOM_NB%GM,IFLD)

! Receive

  IFLD=IFLD+1
  CALL DISGRID_RECV(YDGEOMETRY,IOPROC,1,ZPSF,IFLD)

ENDIF ! IOPROC

! Go back to spectral space

CALL REESPE(YDGEOMETRY,1,1,ZSPS,ZPSF)
YDSPEC%SP(1:NSPEC2)=ZSPS(1:NSPEC2)

WRITE(NULOUT,*) ' '
WRITE(NULOUT,*) '    ====================================== '
WRITE(NULOUT,*) '    ==========  END OF CORMASS2  ========= '
WRITE(NULOUT,*) '    ====================================== '
WRITE(NULOUT,*) ' '
WRITE(NULOUT,*) ' '

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CORMASS2',1,ZHOOK_HANDLE)
END SUBROUTINE CORMASS2

