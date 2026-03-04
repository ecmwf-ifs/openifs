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

SUBROUTINE SPCSI(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY,YDLDDH,YDRIP,YDDYN,YDDYNA,KM,KMLOC,KSTA,KEND,LDONEM,&
 ! --- INOUT -----------------------------------------------------------------
 & PSPVORG,PSPDIVG,PSPTG,PSPSPG,&
 & PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG,&
 ! --- INPUT OPTIONAL --------------------------------------------------------
 & PSPAUXG)

!**** *SPCSI* - SPECTRAL SPACE SEMI-IMPLICIT COMPUTATIONS FOR HYD MODEL.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SPCSI(..)

!        Explicit arguments :
!        --------------------  KM      - Zonal wavenumber 
!                              KMLOC   - Zonal wavenumber (DM-local numbering)
!                              KSTA    - First column processed
!                              KEND    - Last column processed
!                              LDONEM  - T if only one m if processed
!                              PSPVORG - Vorticity columns
!                              PSPDIVG - Divergence columns
!                              PSPTG   - Temperature columns
!                              PSPSPG  - Surface Pressure
!                              PSPTNDSI_VORG - [D vor/Dt]_SI
!                              PSPTNDSI_DIVG - [D div/Dt]_SI
!                              PSPTNDSI_TG   - [D T/Dt]_SI

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-11-24 (before 1997 spcsi.F was part of spc.F)

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      N.Wedi        08-Mar-2005 remove mass correction      
!      K.Yessad 09-Dec-2004: move mass correction in SPCMASCOR + cleanings.
!      K. Yessad 15-May-2006: memory optimisations for stretched geometry
!      N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      K. Yessad (Aug 2009): remove LSITRIC option.
!      F. Voitus: add DDH diagnostics.
!      T.Wilhelmsson 09-09-25: Remove LFULLM requirement for LIMPF
!      K. Yessad (Feb 2012): tests on LL3D, LLDOSI in the caller, simplifications.
!      P. Marguinaud (Nov 2012): Fix unallocated array arguments
!      P. Marguinaud (Sep 2012) : Make PSPAUXG optional
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      O. Marsden (May 2016): Removed redundant geometry arguments
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0       , ONLY : MYSETV
USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : TDYNA
USE YOMLDDH      , ONLY : TLDDH
USE YOMRIP       , ONLY : TRIP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TDYNA)       ,INTENT(IN)    :: YDDYNA
TYPE(TLDDH)       ,INTENT(IN)    :: YDLDDH
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOC
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
LOGICAL           ,INTENT(IN)    :: LDONEM 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPVORG(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDIVG(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTG(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPG(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_VORG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_DIVG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_TG(:,:)
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL :: PSPAUXG(:,:) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB), ALLOCATABLE :: ZSDIVP (:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPDIVP(:,:)

REAL(KIND=JPRB) :: ZSDIV  (YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZHELP  (YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZST    (YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSP    (       KSTA:KEND)

REAL(KIND=JPRB) :: ZSDIVPL (YDGEOMETRY%YRDIMV%NFLEVG,KM:YDGEOMETRY%YRDIM%NSMAX,2)
REAL(KIND=JPRB) :: ZSPDIVPL(YDGEOMETRY%YRDIMV%NFLEVG,KM:YDGEOMETRY%YRDIM%NSMAX,2)

INTEGER(KIND=JPIM) :: II, IN, IOFF, IS0, IS02, ISE, ISPCOL, JLEV, JN, JSP  

REAL(KIND=JPRB) :: ZBDT, ZBDT2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "mxmaop.h"

#include "mxptma.h"
#include "mxture.h"
#include "mxturs.h"
#include "abor1.intfb.h"
#include "sigam.intfb.h"
#include "spcimpfsolve.intfb.h"
#include "sitnu.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPCSI',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  &
 & YDLAP=>YDGEOMETRY%YRLAP, YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & LIMPF=>YDDYN%LIMPF, LSIDG=>YDDYN%LSIDG, RBTS2=>YDDYN%RBTS2, &
 & SIHEG=>YDDYN%SIHEG, SIHEG2=>YDDYN%SIHEG2, SIMI=>YDDYN%SIMI, SIMO=>YDDYN%SIMO, &
 & SIVP=>YDDYN%SIVP, &
 & RSTRET=>YDGEM%RSTRET, &
 & LRSIDDH=>YDLDDH%LRSIDDH, &
 & NPTRSV=>YDMP%NPTRSV, NPTRSVF=>YDMP%NPTRSVF, NSPEC2V=>YDMP%NSPEC2V, &
 & NSPEC2VF=>YDMP%NSPEC2VF, &
 & TDT=>YDRIP%TDT, &
 & SCGMAP=>YDSPGEOM%SCGMAP)
!     ------------------------------------------------------------------

!*       0.    TESTINGS.
!              ---------


IF (LIMPF .AND. .NOT.PRESENT(PSPAUXG)) THEN
  CALL ABOR1(' SPCSI: If LIMPF=T, argument PSPAUXG must be present!')
ENDIF

!     ------------------------------------------------------------------

!*       1.    MEMORY TRANSFER.
!              ----------------

IF (LRSIDDH) THEN
  ! DDH memory transfer
  IF (LIMPF) PSPTNDSI_VORG(1:NFLEVG,KSTA:KEND)=-PSPVORG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_DIVG(1:NFLEVG,KSTA:KEND)=-PSPDIVG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_TG  (1:NFLEVG,KSTA:KEND)=-PSPTG  (1:NFLEVG,KSTA:KEND)
  !the case of surface pressure has not been treated yet
ENDIF

!     ------------------------------------------------------------------

!*       2.    SEMI-IMPLICIT SPECTRAL COMPUTATIONS.
!              ------------------------------------

ALLOCATE(ZSDIVP(NFLEVG,MAX(NSPEC2V,NSPEC2VF)))
ALLOCATE(ZSPDIVP(NFLEVG,MAX(NSPEC2V,NSPEC2VF)))

!*        2.1  Preliminary initialisations.

IF (LDONEM) THEN
  IOFF=NPTRSVF(MYSETV)-1
ELSE
  IOFF=NPTRSV(MYSETV)-1
ENDIF
ISPCOL=KEND-KSTA+1

ZBDT=RBTS2*TDT
ZBDT2=(ZBDT*RSTRET)**2

!*        2.3  Computes right-hand side of Helmholtz equation.

IF( .NOT.LDONEM ) CALL GSTATS(1655,0) ! Main routines and loops in SIGAM chain are parallel
CALL SIGAM(YDGEOMETRY,YDDYN,1,NFLEVG,ZSDIV,PSPTG(:,KSTA:KEND),PSPSPG(KSTA:KEND),ISPCOL,NFLEVG)

IF( .NOT.LDONEM ) CALL GSTATS(1655,1)

IF( .NOT.LDONEM ) CALL GSTATS(1656,0)

IF (LSIDG) THEN
  IF (KM > 0) THEN
!$OMP PARALLEL DO PRIVATE(JSP,JLEV,IN)
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
        IN=YDLAP%NVALUE(JSP+IOFF)
        ZSDIV(JLEV,JSP)=YDLAP%RLAPIN(IN)*PSPDIVG(JLEV,JSP)-ZBDT*ZSDIV(JLEV,JSP)      
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO PRIVATE(JSP,JLEV,IN)
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
        IN=YDLAP%NVALUE(JSP+IOFF)
        ZSDIV(JLEV,JSP)=PSPDIVG(JLEV,JSP)-ZBDT*YDLAP%RLAPDI(IN)*ZSDIV(JLEV,JSP)  
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDIF
ELSE

  ! Case of No Stretching
!$OMP PARALLEL DO PRIVATE(JSP,JLEV,IN)
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      IN=YDLAP%NVALUE(JSP+IOFF)
      ZSDIV(JLEV,JSP)=PSPDIVG(JLEV,JSP)-ZBDT*YDLAP%RLAPDI(IN)*ZSDIV(JLEV,JSP)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ENDIF

!        Add [F] * result to rhs of Helmholtz equation

IF (LIMPF) THEN
!$OMP PARALLEL DO PRIVATE(JSP,JLEV)
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      ZSDIV(JLEV,JSP)=ZSDIV(JLEV,JSP) + PSPAUXG(JLEV,JSP)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ENDIF
IF( .NOT.LDONEM ) CALL GSTATS(1656,1)

!*        2.4  Solve Helmholtz equation

!           Current space --> vertical eigenmodes space.

IF( .NOT.LDONEM ) CALL GSTATS(1660,0) ! MXMAOP Call to SGEMMX Parallelised
CALL MXMAOP(SIMI,1,NFLEVG,ZSDIV,1,NFLEVG,ZSDIVP(:,KSTA:KEND),1,NFLEVG,&
 & NFLEVG,NFLEVG,ISPCOL)  
IF( .NOT.LDONEM ) CALL GSTATS(1660,1)

IF (LSIDG) THEN

  !             Inversion of two tridiagonal systems (Helmholtz equation)
  !                --> (SIMI*DIVprim(t+dt)).

  !             Reorganisation of divergence

  IS0=YDLAP%NSE0L(KMLOC)
  IS02=0
  II=MIN(KM,1)+1
  ZSDIVPL(:,:,:)=0.0_JPRB
  ZSPDIVPL(:,:,:)=0.0_JPRB

  DO JN=KM,NSMAX
    ISE=KSTA+2*(JN-KM)
    ZSDIVPL(:,JN,1:2)=ZSDIVP(:,ISE:ISE+1)
  ENDDO
  IF (KM > 0) THEN

    !               Inversion of a symmetric matrix.

    CALL MXTURS(NSMAX+1-KM,NFLEVG,NFLEVG,II,&
     & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
     & ZSDIVPL,ZSPDIVPL)  
  ELSE

    !               Inversion of a non-symmetric matrix.

    CALL MXTURE(NSMAX+1-KM,NFLEVG,NFLEVG,II,-2,.TRUE.,&
     & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
     & ZSDIVPL,ZSPDIVPL)  
    CALL MXTURE(NSMAX+1-KM,NFLEVG,NFLEVG,II,3,.FALSE.,&
     & SIHEG(1,IS0+1,1),SIHEG2(1,IS02+1,2),&
     & SIHEG2(1,IS02+1,3),ZSDIVPL,ZSPDIVPL)  
  ENDIF

  DO JN=KM,NSMAX
    ISE=KSTA+2*(JN-KM)
    ZSPDIVP(:,ISE:ISE+1)=ZSPDIVPL(:,JN,1:2)
  ENDDO
ELSE

  !             Case with NO Stretching :

  IF (LIMPF) THEN

    !               Solve complex pentadiagonal system

    CALL SPCIMPFSOLVE(YDGEOMETRY,YDRIP,YDDYN,YDDYNA%LNHDYN,YDDYNA%LNHX,LDONEM,ZSDIVP,ZSPDIVP)

  ELSE

    !                 Inversion of a diagonal system (Helmholtz equation)
    !                 --> (SIMI*DIVprim(t+dt)).

    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
        ZSPDIVP(JLEV,JSP)=ZSDIVP(JLEV,JSP)&
         & /(1.0_JPRB-ZBDT2*SIVP(JLEV)*YDLAP%RLAPDI(YDLAP%NVALUE(JSP+IOFF)))  
      ENDDO
    ENDDO
  ENDIF
ENDIF

!           Vertical eigenmodes space --> current space.

IF( .NOT.LDONEM ) CALL GSTATS(1660,0) ! MXMAOP Calls SGEMMX in parallel region
CALL MXMAOP(SIMO,1,NFLEVG,ZSPDIVP(:,KSTA:KEND),1,NFLEVG,PSPDIVG(:,KSTA:KEND),1,&
 & NFLEVG,NFLEVG,NFLEVG,ISPCOL)  
IF( .NOT.LDONEM ) CALL GSTATS(1660,1)

IF (LSIDG) THEN

  !           ZSPDIV=(DIVprim(t+dt)) --> ZSPDIVG=(GM**2 * DIVprim(t+dt)) .

  ZSDIVPL(:,:,:)=0.0_JPRB
  ZSPDIVPL(:,:,:)=0.0_JPRB

  !           Reorganisation of ZSDIVP (Back to the USSR)

  DO JN=KM,NSMAX
    ISE=KSTA+2*(JN-KM)
    ZSDIVPL(:,JN,1:2)=PSPDIVG(:,ISE:ISE+1)
  ENDDO

  !        ZSPDIV=(DIVprim(t+dt)) --> ZPSPDIVG=(GMBAR**2 * DIVprim(t+dt)).

  CALL MXPTMA(NSMAX+1-KM,NFLEVG,NFLEVG,II,SCGMAP(IS0+1,1),&
   & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
   & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
   & ZSDIVPL,ZSPDIVPL)  

  !           Reorganisation of ZSPDIVPL

  DO JN=KM,NSMAX
    ISE=KSTA+2*(JN-KM)
    ZHELP(:,ISE:ISE+1)=ZSPDIVPL(:,JN,1:2)
  ENDDO
ELSE

  !       ZSPDIV=(DIVprim(t+dt)) --> ZSPDIVG=(GMBAR**2 * DIVprim(t+dt)) .

  IF( .NOT.LDONEM ) CALL GSTATS(1656,0)
!$OMP PARALLEL DO PRIVATE(JSP,JLEV)
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      ZHELP(JLEV,JSP)=PSPDIVG(JLEV,JSP)*RSTRET*RSTRET
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  IF( .NOT.LDONEM ) CALL GSTATS(1656,1)

ENDIF

!       If LSIDG:
!         (GM**2 * DIVprim(t+dt)) --> [ tau * (GM**2 * DIVprim(t+dt)) ]
!                                 and [  nu * (GM**2 * DIVprim(t+dt)) ]
!       or if not LSIDG:
!         (GMBAR**2 * DIVprim(t+dt)) --> [ tau * (GMBAR**2 * DIVprim(t+dt)) ]
!                                    and [  nu * (GMBAR**2 * DIVprim(t+dt)) ]

IF( .NOT.LDONEM ) CALL GSTATS(1657,0)  ! Main routines and loops in SITNU chain are parallel
CALL SITNU(YDGEOMETRY,YDDYN,1,NFLEVG,ZHELP,ZST,ZSP,ISPCOL)
IF( .NOT.LDONEM ) CALL GSTATS(1657,1)

!*       2.5  Increment Temperature and surface pressure

IF( .NOT.LDONEM ) CALL GSTATS(1656,0)
!$OMP PARALLEL DO PRIVATE(JSP,JLEV)
DO JSP=KSTA,KEND
  DO JLEV=1,NFLEVG
    PSPTG(JLEV,JSP)=PSPTG(JLEV,JSP)-ZBDT*ZST(JLEV,JSP)
  ENDDO
  PSPSPG(JSP)=PSPSPG(JSP)-ZBDT*ZSP(JSP)
ENDDO
!$OMP END PARALLEL DO
IF( .NOT.LDONEM ) CALL GSTATS(1656,1)

DEALLOCATE(ZSDIVP)
DEALLOCATE(ZSPDIVP)

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF SI TERM AT t+dt FOR DDH.
!              ---------------------------------------

IF (LRSIDDH) THEN
  IF (LIMPF) PSPTNDSI_VORG(1:NFLEVG,KSTA:KEND)=&
   & PSPTNDSI_VORG(1:NFLEVG,KSTA:KEND) + PSPVORG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_DIVG(1:NFLEVG,KSTA:KEND)=PSPTNDSI_DIVG(1:NFLEVG,KSTA:KEND)&
   & + PSPDIVG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_TG(1:NFLEVG,KSTA:KEND)=PSPTNDSI_TG(1:NFLEVG,KSTA:KEND)&
   & + PSPTG(1:NFLEVG,KSTA:KEND)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPCSI',1,ZHOOK_HANDLE)
END SUBROUTINE SPCSI
