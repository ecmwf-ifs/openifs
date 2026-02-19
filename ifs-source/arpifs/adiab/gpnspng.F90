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

SUBROUTINE GPNSPNG(YGFL,YDSPNG,KPROMA,KFLEVG,KST,KND,PGFL)

! GPNSPNG - gridpoint sponge for purely grid-point GFL.

! Purpose
! -------

! Interface
! ---------
!   KPROMA  - horizontal dimension     (in)
!   KFLEVG  - vertical dimension       (in)
!   KST     - start of work            (in)
!   KND     - end of work              (in)
!   PGFL    - GFL variables            (inout)

! Externals
! ---------

! Method
! ------

! Reference
! ---------
!   ECMWF Research Department documentation of the IFS

! Author
! ------
!   K. Yessad (CNRM/GMAP), after routine GPSPNG.
!   October 2011

! Modifications
! -------------
! End Modifications
! -----------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOM_YGFL , ONLY : TYPE_GFLD
USE SPNG_MOD , ONLY : TSPNG

! -----------------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSPNG)       ,INTENT(IN)    :: YDSPNG
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEVG
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KND
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(KPROMA,KFLEVG,YGFL%NDIM)

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF, JGFL
LOGICAL :: LLDO

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPNSPNG',0,ZHOOK_HANDLE)
ASSOCIATE(NDIM=>YGFL%NDIM, NUMFLDS=>YGFL%NUMFLDS, YA=>YGFL%YA, &
 & YCOMP=>YGFL%YCOMP, YI=>YGFL%YI, YL=>YGFL%YL, YO3=>YGFL%YO3, YQ=>YGFL%YQ)
! -----------------------------------------------------------------------------

! perform sponge calculation for grid-point GFL variables (q,ql,qi,qa,O3 only)
DO JGFL=1,NUMFLDS
  IF(YCOMP(JGFL)%LGP) THEN
    LLDO=( (YCOMP(JGFL)%MP == YQ%MP) .OR.&
     & (YCOMP(JGFL)%MP == YL%MP) .OR. (YCOMP(JGFL)%MP == YI%MP) .OR.&
     & (YCOMP(JGFL)%MP == YA%MP) .OR. (YCOMP(JGFL)%MP == YO3%MP) )
    IF (LLDO) THEN
      DO JLEV=1,KFLEVG
        DO JROF=KST,KND
          PGFL(JROF,JLEV,YCOMP(JGFL)%MP)=YDSPNG%RSPONGF(JLEV,3)*PGFL(JROF,JLEV,YCOMP(JGFL)%MP)
        ENDDO
      ENDDO
    ENDIF
  ENDIF
ENDDO

! -----------------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPNSPNG',1,ZHOOK_HANDLE)
END SUBROUTINE GPNSPNG
