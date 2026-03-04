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

SUBROUTINE SUSPECG1(YDGEOMETRY,YGFL,PSPVOR,PSPDIV,PSPT,PSPGFL,&
 & PSPSP,PSPOR,LDINOR,LDSPOR)  

!**** *SUSPECG1*  - Initialize the spectral fields using artificial data

!     Purpose.
!     --------
!           Initialize the spectral fields of the model

!**   Interface.
!     ----------
!        *CALL* *SUSPECG1(.....)*

!        Explicit arguments :
!        --------------------
!               PSPVOR etc. - spectral fields
!               LDINOR - switch for initializing the orography
!               LDSPOR - switch indicating whether spectral orography
!                        field has been read

!     Method.
!     -------
!        See documentation

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      George Mozdzynski *ECMWF*
!      Original : 97-09-18 from SUINIF (CY12R1 MPP version)

!     Modifications.
!     --------------
!      R. El Khatib : 01-08-07 Pruning options
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib : 03 Jul-2013 Limit truncation to 1000 for Vor,Div, Ps, for
!                     the safety of the spectral fit of Ps + ,Bugfix B-level distribution
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST       , ONLY : ROMEGA, RA
USE YOM_YGFL     , ONLY : TYPE_GFLD

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YGFL
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSPVOR(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSPDIV(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSPT(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSPGFL(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSPSP(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSPOR(:) 
LOGICAL           ,INTENT(IN)    :: LDINOR 
LOGICAL           ,INTENT(OUT)   :: LDSPOR 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZREEL(YDGEOMETRY%YRGEM%NGPTOT)

INTEGER(KIND=JPIM) :: IM, IMHW, IN, INST, ISP, J, JLEV,&
 & JMLOC, JN, ISMAX 

REAL(KIND=JPRB) :: ZCC, ZQMAX, ZQMIN, ZSPFACTOR, ZSPVORM1
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "reespe.intfb.h"
#include "speree.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUSPECG1',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  &
 & YDLAP=>YDGEOMETRY%YRLAP, YDSTA=>YDGEOMETRY%YRSTA, YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA,  &
 & YDVFE=>YDGEOMETRY%YRVFE, YDCSGLEG=>YDGEOMETRY%YRCSGLEG,  &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, &
  & YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB, YDGSGEOM=>YDGEOMETRY%YRGSGEOM, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB,  &
  & YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NUMSPFLDS=>YGFL%NUMSPFLDS, YQ=>YGFL%YQ, &
 & NSMAX=>YDDIM%NSMAX, NSPEC2=>YDDIM%NSPEC2, NUMP=>YDDIM%NUMP, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLEVL=>YDDIMV%NFLEVL, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & MYMS=>YDLAP%MYMS, NASM0=>YDLAP%NASM0, &
 & MYLEVS=>YDMP%MYLEVS, NPSP=>YDMP%NPSP, &
 & STTEM=>YDSTA%STTEM)
!     ------------------------------------------------------------------

! CALCULATE 3D INPUT DATA (ROSSBY HAURWITZ WAVE)

! DEFINE WAVE

ISMAX=MIN(NSMAX,1000) ! beyond this value spectral fit apparently does not work
PSPVOR=0._JPRB
PSPDIV=0._JPRB

DO JLEV=1,NFLEVL
  DO JMLOC=1,NUMP
    IM = MYMS(JMLOC)
    DO JN=IM,ISMAX
      ISP = NASM0(IM)+(JN-IM)*2
      PSPVOR(JLEV,ISP)=1.E-8_JPRB * IM * JN * MYLEVS(JLEV)&
       & / ( REAL(ISMAX*ISMAX,JPRB) * REAL(NFLEVG,JPRB) )  
      PSPDIV(JLEV,ISP)=0.0_JPRB
      IF (IM  /=  0) THEN
        PSPVOR(JLEV,ISP+1)=0.0_JPRB
        PSPDIV(JLEV,ISP+1)=0.0_JPRB
      ENDIF
    ENDDO
  ENDDO
ENDDO

IN=4
IMHW=3
DO JLEV=1,NFLEVL
  DO JMLOC=1,NUMP
    IM = MYMS(JMLOC)
    IF (IM == IMHW) THEN
      PSPVOR(JLEV,NASM0(IM)+(IN-IM)*2)=2.E-6_JPRB
      EXIT
    ENDIF
  ENDDO
ENDDO

IF (NPSP == 1) THEN
  DO J=1,NSPEC2
    PSPSP(J)=0.0_JPRB
  ENDDO

! SOLVE LINEAR BALANCE EQUATION.

  ZCC=-2.0_JPRB*ROMEGA*RA**2
  DO JMLOC=1,NUMP
    IM = MYMS(JMLOC)
    IF (IM  <  ISMAX) THEN
      INST=MAX(1,IM)
      DO JN=INST,ISMAX-1
        ISP=NASM0(IM)+(JN-IM)*2
        IF(JN > IM) THEN
          ZSPVORM1=1.E-8_JPRB * IM * (JN-1) / (REAL(ISMAX*ISMAX,JPRB)*&
           & REAL(NFLEVG,JPRB))  
          IF ((IM == IMHW).AND.(JN-1) == IN) ZSPVORM1=2.E-6_JPRB
        ELSE
          ZSPVORM1=0._JPRB
        ENDIF
        IF ((IM == IMHW).AND.(JN+1) == IN) THEN
          ZSPFACTOR=2.E-6_JPRB
        ELSE
          ZSPFACTOR=1.E-8_JPRB*IM*(JN+1)/(REAL(ISMAX*ISMAX,JPRB)*&
           & REAL(NFLEVG,JPRB))  
        ENDIF
        PSPSP(ISP)=ZCC/(JN*(JN+1))*(&
         & SQRT(REAL((2*JN-1)*(JN+IM)*(JN-IM),JPRB)/REAL(2*JN+1,JPRB))*&
         & (1.0_JPRB+1.0_JPRB/REAL(JN,JPRB))*1.0_JPRB/REAL(2*JN-1,JPRB)*ZSPVORM1+&
         & SQRT(REAL((2*JN+3)*(JN-IM+1),JPRB)/&
         & REAL((2*JN+1)*(JN+IM+1),JPRB))*&
         & REAL(JN+IM+1,JPRB)/REAL(2*JN+3,JPRB)*ZSPFACTOR)&
         & *(1.0_JPRB-1.0_JPRB/REAL(JN+1,JPRB))  
      ENDDO
    ENDIF
  ENDDO
  IF (NUMP > 0.AND. MYMS(1) == 0) THEN
    PSPSP(NASM0(0))=1.E5_JPRB
  ENDIF

! CONVERT TO LN(PS)

  ZREEL(:)=0.0_JPRB
  CALL SPEREE(YDGEOMETRY,1,1,PSPSP,ZREEL)
  DO J=1,NGPTOT
    ZREEL(J)=LOG(ZREEL(J))
  ENDDO
  CALL REESPE(YDGEOMETRY,1,1,PSPSP,ZREEL)


ELSE

! PSPSP NOT USED FOR THIS PROCESS, RESET TO UNDEFINED VALUE

  DO J=1,NSPEC2
    PSPSP(J)=-999._JPRB
  ENDDO
ENDIF

! DEFINE TEMPERATURE AND HUMIDITY

IF(NUMSPFLDS > 0) THEN
  PSPGFL(:,:,:) = 0.0_JPRB
ENDIF
DO JLEV=1,NFLEVL
  DO JMLOC=1,NUMP
    IM = MYMS(JMLOC)
    DO JN=IM,NSMAX
      ISP = NASM0(IM)+(JN-IM)*2
      PSPT(JLEV,ISP)=1.E-05_JPRB * IM * JN * MYLEVS(JLEV)&
       & / ( REAL(NSMAX*NSMAX,JPRB) * REAL(NFLEVG,JPRB) )  
      IF( YQ%LSP )THEN
        PSPGFL(JLEV,ISP,YQ%MPSP)=1.E-09_JPRB * IM * JN * MYLEVS(JLEV)&
         & / ( REAL(NSMAX*NSMAX,JPRB) * REAL(NFLEVG,JPRB) )  
      ENDIF
      IF (IM  /=  0) THEN
        PSPT(JLEV,ISP+1)=0.0_JPRB
        IF( YQ%LSP )THEN
          PSPGFL(JLEV,ISP+1,YQ%MPSP)=0.0_JPRB
        ENDIF
      ENDIF
    ENDDO
  ENDDO
ENDDO
ZQMIN= 0.25E-5_JPRB
ZQMAX= 0.01_JPRB
IF (NUMP > 0.AND. MYMS(1) == 0) THEN
  DO JLEV=1,NFLEVL
    PSPT(JLEV,NASM0(0))=STTEM(MYLEVS(JLEV))
    IF( YQ%LSP )THEN
      PSPGFL(JLEV,NASM0(0),YQ%MPSP)=ZQMIN+ (ZQMAX-ZQMIN)*&
       & MAX(0.0_JPRB,REAL(2*MYLEVS(JLEV),JPRB)/REAL(NFLEVG,JPRB)-1.0_JPRB)  
    ENDIF
  ENDDO
ENDIF

! DEFINE OROGRAPHY

IF( LDINOR )THEN
  PSPOR(:)=0.0_JPRB
  DO JMLOC=1,NUMP
    IM = MYMS(JMLOC)
    DO JN=IM,NSMAX
      ISP = NASM0(IM)+(JN-IM)*2
      PSPOR(ISP)=1.E-05_JPRB * IM * JN  /  REAL(NSMAX*NSMAX,JPRB)
      IF (IM  /=  0) THEN
        PSPOR(ISP+1)=0.0_JPRB
      ENDIF
    ENDDO
  ENDDO
  IF (NUMP > 0.AND. MYMS(1) == 0) THEN
    PSPOR(NASM0(0))=1000._JPRB
  ENDIF
  LDSPOR=.TRUE.
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSPECG1',1,ZHOOK_HANDLE)
END SUBROUTINE SUSPECG1
