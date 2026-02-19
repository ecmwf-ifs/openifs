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

SUBROUTINE DOTPROD3(YDGEOMETRY,KF,PGP1,PGP2,PPROD)

!**** *DOTPROD3*  - ROUTINE FOR CALCULATING GRID POINT NORM OF 3-D FIELDS PGP1 AND PGP2

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        *CALL* *DOTPROD3(...)*

!        INPUT:

!        KF      - Field index
!        PGP1    - Field set 1
!        PGP2    - Field set 2
!        PPROD   - Dot product

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------

!        T.Wilhelmsson *ECMWF*

!     MODIFICATIONS.
!     --------------
!       Original : 25 Oct 2011
!     -----------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMLUN    ,ONLY : NULOUT
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE ORDER_INDEPENDENT_SUMMATION_MOD, ONLY : ORDER_INDEP_GLOBAL_SUM
IMPLICIT NONE

TYPE(GEOMETRY),     INTENT(IN)  :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN)  :: KF
REAL(KIND=JPRB),    INTENT(IN)  :: PGP1(:,:,:,:)
REAL(KIND=JPRB),    INTENT(IN)  :: PGP2(:,:,:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PPROD

INTEGER(KIND=JPIM) :: IEND, IST, JL, JKGLO, IBL, JLEV
INTEGER(KIND=JPIM) :: IPROMA, ILEVS, IGPBLKS
REAL(KIND=JPRB)    :: ZTMP(1)
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


#include "abor1.intfb.h"

!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('DOTPROD3',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, &
 & NGPTOT=>YDGEM%NGPTOT)

IF (ANY(SHAPE(PGP1) /= SHAPE(PGP2))) THEN
  WRITE(NULOUT,*)'DOTPROD3: PGP1 SHAPE ', SHAPE(PGP1)
  WRITE(NULOUT,*)'DOTPROD3: PGP2 SHAPE ', SHAPE(PGP2)
  CALL ABOR1("DOTPROD3 GMV different shapes")
ENDIF

ILEVS = SIZE(PGP1,2)
IF(ILEVS == 0) THEN
  PPROD = 0.0_JPRB
  IF (LHOOK) CALL DR_HOOK('DOTPROD3',1,ZHOOK_HANDLE)  
  RETURN
ENDIF

IPROMA = SIZE(PGP1,1)
IF(IPROMA < NPROMA) THEN
  WRITE(NULOUT,*)'DOTPROD3:FIRST DIM. OF ARGUMENTS TOO SMALL ',IPROMA,NPROMA
  CALL ABOR1('DOTPROD3:FIRST DIMENSION OF ARGUMENTS TOO SMALL ')
ENDIF

IGPBLKS = SIZE(PGP1,4)
IF(IGPBLKS < NGPBLKS) THEN
  WRITE(NULOUT,*)'DOTPROD3:FOUTH DIM. OF ARGUMENTS TOO SMALL ',IGPBLKS,NGPBLKS
  CALL ABOR1('DOTPROD3:FOURTH DIMENSION OF ARGUMENTS TOO SMALL ')
ENDIF

ZTMP(1) = 0.0_JPRB
DO JKGLO=1,NGPTOT,NPROMA
  IST=1
  IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
  IBL=(JKGLO-1)/NPROMA+1
  DO JLEV=1,ILEVS
    DO JL=IST,IEND
        ZTMP(1) = ZTMP(1) + PGP1(JL,JLEV,KF,IBL)*PGP2(JL,JLEV,KF,IBL)
    ENDDO
  ENDDO
ENDDO

PPROD = ORDER_INDEP_GLOBAL_SUM(ZTMP,KNG=1)

!     -----------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('DOTPROD3',1,ZHOOK_HANDLE)
END SUBROUTINE DOTPROD3
