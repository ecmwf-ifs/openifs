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

SUBROUTINE SPECRT(YDGEOMETRY,YDGFL,YGFL,KDIR,YDSP,LDINCR,YDGFL5,PTRAJEC)

!**** *SPECRT* - Computes spectral coeffs of 'TV' (TV defined as RT/Rd) or
!                spectral coeffs of T from TV and Q

!     Purpose.
!     --------
!           COMPUTES SPECTRAL COEFFICIENTS OF 'TV' (i.e., RT/Rd) or T

!**   Interface.
!     ----------
!        *CALL* *SPECRT(...)

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

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Clive Temperton  *ECMWF*
!      Original : 94-05-16

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Modified 04-02-06: Y. Seity : new arguments to GPRCP
!      M.Hamrud  15-Jan-2006  Revised GPRCP
!      J.Hague   01-Aug-2012  Add optional YDSP
!      F. Vana  28-Nov-2013 : Redesigned trajectory handling

!     ------------------------------------------------------------------

USE GEOMETRY_MOD  , ONLY : GEOMETRY
USE YOMGFL        , ONLY : TGFL
USE PARKIND1      , ONLY : JPIM, JPRB
USE YOMHOOK       , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST        , ONLY : RD
USE YOM_YGFL      , ONLY : TYPE_GFLD
USE TRAJECTORY_MOD, ONLY : LTRAJGP
USE YOMLUN        , ONLY : NULOUT
USE YOMTRAJ       , ONLY : TRAJ_TYPE, LPRTTRAJ
!USE SPECTRAL_FIELDS_MOD, ONLY: SPECTRAL_FIELD, ASSIGNMENT(=)
USE SPECTRAL_FIELDS_MOD

IMPLICIT NONE

TYPE(GEOMETRY),            INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL),                INTENT(INOUT) :: YDGFL
TYPE(TYPE_GFLD),           INTENT(IN)    :: YGFL
INTEGER(KIND=JPIM),        INTENT(IN)    :: KDIR 
TYPE(SPECTRAL_FIELD),      INTENT(INOUT) :: YDSP
LOGICAL,         OPTIONAL, INTENT(IN)    :: LDINCR
TYPE(TGFL),      OPTIONAL, INTENT(INOUT) :: YDGFL5
TYPE(TRAJ_TYPE), OPTIONAL, INTENT(IN)    :: PTRAJEC

REAL(KIND=JPRB) :: ZT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB), ALLOCATABLE :: ZR(:,:)
LOGICAL         :: LLINCR
INTEGER(KIND=JPIM) :: IEND, IST, JKGLO, JL, JLEV, IBL  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "inv_trans.h"
#include "dir_trans.h"
#include "gprcp.intfb.h"
#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPECRT',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NDIM5=>YGFL%NDIM5, YI=>YGFL%YI, YL=>YGFL%YL, YQ=>YGFL%YQ, &
 & NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, NRESOL=>YDDIM%NRESOL, &
 & NSPEC2=>YDDIM%NSPEC2, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & NBSETLEV=>YDMP%NBSETLEV)
!     ------------------------------------------------------------------

!*       1.    COMPUTES SPECTRAL COEFFS OF 'TV' OR T
!              -------------------------------------
!   If the fields contain increments we use the trajectory q field to
!   calculate the convertion from Tv to T. The GFL field contains
!   the humidity increment, so it cannot be used to do the convertion
IF(PRESENT(LDINCR)) THEN
  LLINCR=LDINCR
ELSE
  LLINCR=.FALSE.
ENDIF
IF(LLINCR.AND.(.NOT. PRESENT(PTRAJEC))) THEN
  CALL ABOR1('MISSING TRAJECTORY OF MOISTURE TO CONVERT Tv->T')
ENDIF
IF(LLINCR.AND.(.NOT. PRESENT(YDGFL5))) THEN
  CALL ABOR1('MISSING YDGFL5 OF MOISTURE TO CONVERT Tv->T')
ENDIF

IF (LLINCR) THEN
  IF( LTRAJGP )THEN

!   Retrieve grid-point trajectory. We only need humidity

    IF (.NOT.ALLOCATED(YDGFL5%GFL5))  ALLOCATE (YDGFL5%GFL5(NPROMA,NFLEVG,NDIM5,NGPBLKS))
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IEND,IBL)
    DO JKGLO=1,NGPTOT,NPROMA
      IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
      IBL=(JKGLO-1)/NPROMA+1
      IF(.NOT.ASSOCIATED(PTRAJEC%MAIN(IBL)%GFL))&
         & CALL ABOR1('READING TRAJEC%MAIN: TRAJEC%MAIN%GFL NOT ALLOCATED')  
      YDGFL5%GFL5(1:IEND,:,1:NDIM5,IBL)=PTRAJEC%MAIN(IBL)%GFL(1:IEND,:,1:NDIM5)
    ENDDO
!$OMP END PARALLEL DO
    IF (LPRTTRAJ)  WRITE(NULOUT,*)'GREPTRAJ GET_TRAJ_GRID GRID/GFL in SPECRT'
  ELSE
    IF(.NOT.ALLOCATED(YDGFL5%GFL5)) CALL ABOR1('SPECRT:GFL5 NOT ALLOCATED')
  ENDIF
ENDIF

CALL GSTATS(25,0)

CALL INV_TRANS(PSPSCALAR=YDSP%T,KVSETSC=NBSETLEV,KRESOL=NRESOL,KPROMA=NPROMA,PGP=ZT)

IF (.NOT.LLINCR) THEN
  IF(YQ%LSP) THEN
    CALL INV_TRANS(PSPSCALAR=YDSP%Q,KVSETSC=NBSETLEV,KRESOL=NRESOL,KPROMA=NPROMA,&
       & PGP=YDGFL%GFL(:,:,YQ%MP,:))  
  ENDIF
  IF(YL%LSP) THEN
    CALL INV_TRANS(PSPSCALAR=YDSP%L,KVSETSC=NBSETLEV,KRESOL=NRESOL,KPROMA=NPROMA,&
       & PGP=YDGFL%GFL(:,:,YL%MP,:))  
  ENDIF
  IF(YI%LSP) THEN
    CALL INV_TRANS(PSPSCALAR=YDSP%I,KVSETSC=NBSETLEV,KRESOL=NRESOL,KPROMA=NPROMA,&
       & PGP=YDGFL%GFL(:,:,YI%MP,:))  
  ENDIF
ENDIF

CALL GSTATS(1038,0)

!$OMP PARALLEL PRIVATE(JKGLO,IST,IEND,JLEV,JL,IBL,ZR)
ALLOCATE(ZR(NPROMA,NFLEVG))
!$OMP DO SCHEDULE(STATIC)
DO JKGLO=1,NGPTOT,NPROMA
  IST =1
  IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
  IBL=(JKGLO-1)/NPROMA+1
  
  IF (.NOT.LLINCR) THEN
    CALL GPRCP(NPROMA,IST,IEND,NFLEVG,PGFL=YDGFL%GFL(:,:,:,IBL),PR=ZR)
  ELSE
    CALL GPRCP(NPROMA,IST,IEND,NFLEVG,PQ=YDGFL5%GFL5(:,:,YQ%MP5,IBL),PR=ZR)
  ENDIF
  IF (KDIR == 0) THEN
    DO JLEV=1,NFLEVG
      DO JL=IST,IEND
        ZT(JL,JLEV,IBL)=ZR(JL,JLEV)*ZT(JL,JLEV,IBL)/RD
      ENDDO
    ENDDO
  ELSE
    DO JLEV=1,NFLEVG
      DO JL=IST,IEND
        ZT(JL,JLEV,IBL)=RD*ZT(JL,JLEV,IBL)/ZR(JL,JLEV)
      ENDDO
    ENDDO
  ENDIF
ENDDO
!$OMP END DO
DEALLOCATE(ZR)
!$OMP END PARALLEL
CALL GSTATS(1038,1)

CALL DIR_TRANS(PSPSCALAR=YDSP%T,KVSETSC=NBSETLEV,KRESOL=NRESOL,KPROMA=NPROMA,PGP=ZT)

CALL GSTATS(25,1)
!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPECRT',1,ZHOOK_HANDLE)
END SUBROUTINE SPECRT
