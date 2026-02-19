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

SUBROUTINE GPIAU(YGFL,KPROMA,KFLEVG,KST,KND,PGFL,PIAUGFL)

! GPIAU - gridpoint IAU for purely grid-point GFL.

! Purpose
! -------

! Interface
! ---------
!   KPROMA  - horizontal dimension     (in)
!   KFLEVG  - vertical dimension       (in)
!   KST     - start of work            (in)
!   KND     - end of work              (in)
!   PGFL    - GFL variables            (inout)
!   PIAUGFL - analysis increment fraction for GFL (in)

! Externals
! ---------

! Method
! ------

! Reference
! ---------
!   ECMWF Research Department documentation of the IFS

! Author
! ------
!   P. Brousseau (CNRM/GMAP).
!   February 2014

! Modifications
! -------------
!   Sept. 2018  P. Brousseau   Take care of negative humidity
! End Modifications
! -----------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOM_YGFL , ONLY : TYPE_GFLD

! -----------------------------------------------------------------------------

IMPLICIT NONE

TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEVG
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KND
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(KPROMA,KFLEVG,YGFL%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIAUGFL(KPROMA,KFLEVG,1)

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF, JGFL
LOGICAL :: LLDO

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPIAU',0,ZHOOK_HANDLE)
! -----------------------------------------------------------------------------

! perform sponge calculation for grid-point GFL variables (q,ql,qi,qa,O3 only)
DO JGFL=1,YGFL%NUMFLDS
  IF(YGFL%YCOMP(JGFL)%LGP) THEN
    LLDO=(YGFL%YCOMP(JGFL)%MP == YGFL%YQ%MP)
    IF (LLDO) THEN
      DO JLEV=1,KFLEVG
        DO JROF=KST,KND
          PGFL(JROF,JLEV,YGFL%YCOMP(JGFL)%MP)=MAX(0._JPRB, &
           & PGFL(JROF,JLEV,YGFL%YCOMP(JGFL)%MP) + PIAUGFL(JROF,JLEV,1))
        ENDDO
      ENDDO
    ENDIF
  ENDIF
ENDDO

! -----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPIAU',1,ZHOOK_HANDLE)
END SUBROUTINE GPIAU
