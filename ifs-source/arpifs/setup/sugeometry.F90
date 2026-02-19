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

SUBROUTINE SUGEOMETRY(YDGEOMETRY,KSUPERSEDE)

!**** *SUGEOMETRY*  - SETUP MODEL GEOMETRY

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        *CALL* *SUGEOMETRY*

!        EXPLICIT ARGUMENTS
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS
!        --------------------
!        NONE

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      TOMAS WILHELMSSON *ECMWF*
!      Assembled from parts of SU0YOMA and SU0YOMB  2013-08-20

!     MODIFICATIONS.
!     --------------
!     O. Marsden  Aug 2016 Transmit MYLEVS array to trans library, necessary for call to FSPGLH
!     K. Yessad (Feb 2018): remove deep-layer formulations.
!     R. El Khatib 14-May-2018 move allocation of yrorog from suorog/sueorog
!     ------------------------------------------------------------------

USE TYPE_GEOMETRY , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM,JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMCT0   , ONLY : LELAM
USE YOMVV1   , ONLY : DVALH, DVBH
USE YOMVERT  , ONLY : VP00
USE YOMCVER  , ONLY : SUCVER_GEOM,PRT_CVER_GEOM
USE SPGEOM_MOD,ONLY : SUSPGEOM
USE IOSTREAM_MIX, ONLY : INI_IOSTREAM

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), TARGET,INTENT(INOUT) :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KSUPERSEDE
REAL(KIND=JPRB), ALLOCATABLE :: ZGM(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZGNORX(:,:), ZGNORY(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZGELAM(:,:), ZGELAT(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZVALH(:), ZVBH(:)

INTEGER(KIND=JPIM) :: JLEV
INTEGER(KIND=JPIM) :: JKGLO, IEND, IBL
CHARACTER(LEN=35)  ::  CLINE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#ifdef WITH_ATLAS
#include "suatlas_mesh.intfb.h"
#endif
#include "sudim.intfb.h"
#include "suegem_naml.intfb.h"
#include "suegem1a.intfb.h"
#include "suegem1b.intfb.h"
#include "suegem2.intfb.h"
#include "suegeolbc.intfb.h"
#include "suelap.intfb.h"
#include "sugem_naml.intfb.h"
#include "sugem1a.intfb.h"
#include "sugem1b.intfb.h"
#include "sugem2.intfb.h"
#include "sugem3.intfb.h"
#include "suemp.intfb.h"
#include "suetrans.intfb.h"
#include "sulap.intfb.h"
#include "sulega.intfb.h"
#include "suelega.intfb.h"
#include "sump.intfb.h"
#include "susta.intfb.h"
#include "sutrans.intfb.h"
#include "suvert.intfb.h"
#include "suvert2.intfb.h"
#include "suvv1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUGEOMETRY',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  &
 & YDSPGEOM=>YDGEOMETRY%YSPGEOM, YDCSGEOM=>YDGEOMETRY%YRCSGEOM, YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,  &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, &
 & YDLAP=>YDGEOMETRY%YRLAP, YDCSGLEG=>YDGEOMETRY%YRCSGLEG, YDSTA=>YDGEOMETRY%YRSTA)
ASSOCIATE(NDGENG=>YDDIM%NDGENG, NDGLG=>YDDIM%NDGLG, NDGSAG=>YDDIM%NDGSAG, &
 & NDLSUR=>YDDIM%NDLSUR, &
 & NFLEVG=>YDDIMV%NFLEVG, NIOLEVG=>YDDIMV%NIOLEVG)

!     ------------------------------------------------------------------

CLINE='----------------------------------'


YDGEOMETRY%LNONHYD_GEOM=>YDGEOMETRY%YRVERT_GEOM%LNONHYD_GEOM
YDGEOMETRY%YRCVER => YDGEOMETRY%YRVERT_GEOM%YRCVER
!*    Initialize geometry dimensions
WRITE(NULOUT,*) '------ Set up geometry dimensions  ------',CLINE
CALL SUDIM(YDGEOMETRY,KSUPERSEDE)


!*    Initialize geometry (namelist variables); read NAMGEM and NEMGEO
IF (LELAM) THEN
  WRITE(NULOUT,*) '---- Set up geometry (namelist variables) for LAM model -----',CLINE
  ALLOCATE(ZGM(0:NDLSUR,NDGSAG:NDGENG))
  ALLOCATE(ZGNORX(0:NDLSUR,NDGSAG:NDGENG))
  ALLOCATE(ZGNORY(0:NDLSUR,NDGSAG:NDGENG))
  ALLOCATE(ZGELAM(0:NDLSUR,NDGSAG:NDGENG))
  ALLOCATE(ZGELAT(0:NDLSUR,NDGSAG:NDGENG))
  CALL SUEGEM_NAML(YDGEOMETRY,ZGM,ZGNORX,ZGNORY,ZGELAM,ZGELAT,KSUPERSEDE)
ELSE
  WRITE(NULOUT,*) '---- Set up geometry, (namelist variables) for global model -----',CLINE
  CALL SUGEM_NAML(YDGEOMETRY,KSUPERSEDE)
ENDIF
WRITE(NULOUT,*) '---- Set up cver -------',CLINE
CALL SUCVER_GEOM(YDGEOMETRY%YRCVER,YDGEOMETRY%LNONHYD_GEOM)
CALL PRT_CVER_GEOM(YDGEOMETRY%YRCVER)

!*    Allocate and initialize vertical coordinate
WRITE(NULOUT,*) '---- Set up vertical coordinate -------',CLINE
CALL SUVV1(YDGEOMETRY,KSUPERSEDE)
!*    Setup transforms for this resolution
WRITE(NULOUT,*) '------ Set up transforms for this resolution --',CLINE
IF (LELAM) THEN
  CALL SUETRANS(YDGEOMETRY)
ELSE
  CALL SUTRANS(YDGEOMETRY)
ENDIF

!*    Initialize geometry (very first part)
IF (LELAM) THEN
  WRITE(NULOUT,*) '---- Set up geometry, part 1a for LAM model -----',CLINE
  CALL SUEGEM1A(YDGEOMETRY)
ELSE
  WRITE(NULOUT,*) '---- Set up geometry, part 1a for global model -----',CLINE
  CALL SUGEM1A(YDGEOMETRY,KSUPERSEDE)
ENDIF

!*    Initialize distributed memory partitioning parameters
WRITE(NULOUT,*) '-- Set up distributed mem partition -',CLINE
IF (LELAM) THEN
  CALL SUEMP(YDGEOMETRY)
ELSE
  CALL SUMP(YDGEOMETRY)
ENDIF

!*    Initialize geometry variables used for LBC
IF (LELAM) THEN
  WRITE(NULOUT,*) '-- Set up geometry variables used for LBC -',CLINE
  CALL SUEGEOLBC(YDGEOMETRY)
ENDIF

IF (LELAM) THEN
  !*    Initialize YRCSGLEG (YOMLEG) for LAM models
  WRITE(NULOUT,*) '---- Set up YRCSGLEG for LAM models ----',CLINE
  CALL SUELEGA(YDCSGLEG,YDGEOMETRY%YRDIM) 
ELSE
  !*    Initialize Legendre polynomials and YRCSGLEG (YOMLEG) for spherical geometry
  WRITE(NULOUT,*) '---- Set up Legendre polynomials and YRCSGLEG ----',CLINE
  CALL SULEGA(YDCSGLEG,YDGEOMETRY%YRDIM)
ENDIF

!*    Initialize Laplace space constants commons
WRITE(NULOUT,*) '---- Set up Laplace space constants -',CLINE
CALL SULAP(YDLAP,YDGEOMETRY%YRDIM)
IF (LELAM) CALL SUELAP(YDGEOMETRY%YRDIM,YDGEOMETRY%YREDIM,YDGEOMETRY%YRLAP,YDGEOMETRY%YRELAP,YDGEOMETRY%YREGEO)

!*    Initialize vertical geometry
WRITE(NULOUT,*) '---- Set up vertical geometry -------',CLINE
IF (NIOLEVG == 1 .AND. NFLEVG == 1) THEN
  ! for example: 2D models (LR2D=T), conf 923,931,932,933.
  CALL SUVERT2(YDGEOMETRY)

  CALL INI_IOSTREAM(PVAH=YDGEOMETRY%YRVAB%VAH(0:1),PVBH=YDGEOMETRY%YRVAB%VBH(0:1))
ELSE
  CALL SUVERT(YDGEOMETRY)

  IF (NIOLEVG == NFLEVG) THEN
    CALL INI_IOSTREAM(PVAH=YDGEOMETRY%YRVAB%VAH(:),PVBH=YDGEOMETRY%YRVAB%VBH(:))
  ELSE
    ! This case is general ; but the case above is kept for safety
    ALLOCATE(ZVALH(0:NIOLEVG))
    ALLOCATE(ZVBH(0:NIOLEVG))
    DO JLEV=0,NIOLEVG
      ZVALH(JLEV)=REAL(DVALH(JLEV),JPRB)/VP00
      ZVBH (JLEV)=REAL(DVBH(JLEV),JPRB)
    ENDDO
    CALL INI_IOSTREAM(KFLEVG=NIOLEVG,PVAH=ZVALH,PVBH=ZVBH)
    DEALLOCATE(ZVALH,ZVBH)
  ENDIF
ENDIF

!*    Initialize standard atmosphere
WRITE(NULOUT,*) '---- Set up standard atmosphere -----',CLINE
CALL SUSTA(YDGEOMETRY)

!*    Initialize horizontal geometry (part 1b) and Legendre polynomials
IF (LELAM) THEN
  WRITE(NULOUT,*) '---- Set up geometry, part 1b for LAM model -----',CLINE
  CALL SUEGEM1B(YDGEOMETRY) 
ELSE
  WRITE(NULOUT,*) '---- Set up geometry, part 1b for global model -------',CLINE
  WRITE(NULOUT,*) '---- Set up Legendre polynomials ----',CLINE
  CALL SUGEM1B(YDGEOMETRY)
  CALL INI_IOSTREAM(PLATIG=YDCSGLEG%RLATIG(1:NDGLG))
ENDIF

!*    Initialize spectral geometry (GMR,(E)SCGMAP)
WRITE(NULOUT,*) '---- Set up spectral geometry (GMR,(E)SCGMAP) ----',CLINE
CALL SUSPGEOM(YDGEOMETRY)

!*    Initialize geometry (second part)
!     in particuliar, allocate and fill YRGSGEOM and YRCSGEOM
IF (LELAM) THEN
  WRITE(NULOUT,*) '---- Set up geometry, part 2 for LAM (+ set up biper. fit of gm) ----',CLINE
  CALL SUEGEM2(YDGEOMETRY,ZGM,ZGNORX,ZGNORY,ZGELAM,ZGELAT)
  DEALLOCATE(ZGM,ZGNORX,ZGNORY,ZGELAM,ZGELAT)
ELSE
  WRITE(NULOUT,*) '---- Set up geometry, part 2 for global model --------',CLINE
  CALL SUGEM2(YDGEOMETRY)
ENDIF

!*    Initialize geometry (third part)
!     in particuliar, reallocate and fill NLOEN, NMEN, YDCSGLEG%RLATI
WRITE(NULOUT,*) '---- Set up geometry, part 3 --------',CLINE
CALL SUGEM3(YDGEOMETRY)



#ifdef WITH_ATLAS
WRITE(NULOUT,*) '---- Set up ATLAS mesh -------',CLINE
CALL SUATLAS_MESH(YDGEOMETRY)
#endif
!* Allocate orography and derivatives (moved from suorog)
WRITE(NULOUT,*) '---- Allocate orography and derivatives -------',CLINE
ALLOCATE(YDGEOMETRY%YROROG(YDGEOMETRY%YRDIM%NGPBLKS))
DO JKGLO=1,YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIM%NPROMA
  IEND=MIN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRGEM%NGPTOT-JKGLO+1)
  IBL=(JKGLO-1)/YDGEOMETRY%YRDIM%NPROMA+1
  ALLOCATE(YDGEOMETRY%YROROG(IBL)%OROG (1:IEND))
  ALLOCATE(YDGEOMETRY%YROROG(IBL)%OROGL(1:IEND))
  ALLOCATE(YDGEOMETRY%YROROG(IBL)%OROGM(1:IEND))
  IF(YDGEOMETRY%LNONHYD_GEOM) THEN
    ALLOCATE(YDGEOMETRY%YROROG(IBL)%OROGLL(1:IEND))
    ALLOCATE(YDGEOMETRY%YROROG(IBL)%OROGLM(1:IEND))
    ALLOCATE(YDGEOMETRY%YROROG(IBL)%OROGMM(1:IEND))
  ENDIF
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUGEOMETRY',1,ZHOOK_HANDLE)
END SUBROUTINE SUGEOMETRY
