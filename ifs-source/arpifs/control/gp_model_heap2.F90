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

SUBROUTINE GP_MODEL_HEAP2(YDGEOMETRY,YDFIELDS,YDMODEL,CDCONF,LD_DFISTEP,PTRAJEC,PTRAJEC_OOPS)
!****-------------------------------------------------------------------
!**** *GP_MODEL_HEAP2* - Grid-point model
!****-------------------------------------------------------------------
!     Purpose.   gp_model driver using memory on heap for large arrays
!     --------   (allocates memory -- if possible -- only once) 

!**   Interface.
!     ----------
!        *CALL* *GP_MODEL_HEAP2 (..) -- for use by NOPT_MEMORY=2

!        Explicit arguments :  CDCONF - configuration of work (see doc.)
!        --------------------

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   See includes below.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        S. Saarinen, ECMWF, 30-Oct-2017

! Modifications
! -------------
! End Modifications
!-------------------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE FIELDS_MOD   , ONLY : FIELDS
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMTRAJ      , ONLY : TRAJ_TYPE
USE YOMTRAJ_OOPS , ONLY : TRAJ_TYPE_OOPS

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(FIELDS)      ,INTENT(INOUT) :: YDFIELDS
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
CHARACTER(LEN=9)  ,INTENT(IN)    :: CDCONF
LOGICAL           ,INTENT(IN)    :: LD_DFISTEP
TYPE(TRAJ_TYPE)        ,OPTIONAL,INTENT(INOUT) :: PTRAJEC
TYPE(TRAJ_TYPE_OOPS)   ,OPTIONAL,INTENT(INOUT) :: PTRAJEC_OOPS

!     ------------------------------------------------------------------
LOGICAL, SAVE :: LLINITALLOC = .FALSE.
LOGICAL :: LLREALLOC
INTEGER(KIND=JPIM), SAVE :: IDX_IL0(3) = (/0,0,0/)
INTEGER(KIND=JPIM), SAVE :: IDX_ZLSCAW(4) = (/0,0,0,0/)
INTEGER(KIND=JPIM), SAVE :: IDX_ZSLBUF1(2) = (/0,0/)
INTEGER(KIND=JPIM), SAVE :: IDX_ZSLBUF2(3) = (/0,0,0/)
INTEGER(KIND=JPIM), ALLOCATABLE, SAVE :: IL0(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE, SAVE :: ZLSCAW(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE, SAVE :: ZSLBUF1(:,:), ZSLBUF2(:,:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gp_model.intfb.h"

IF (LHOOK) CALL DR_HOOK('GP_MODEL_HEAP2',0,ZHOOK_HANDLE)

LLREALLOC = .FALSE.
IF (LLINITALLOC) THEN
   IF (.not.ALL(IDX_IL0 == &
        & (/YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIM%NGPBLKS/))) LLREALLOC = .TRUE.
   IF (.not.ALL(IDX_ZLSCAW == &
        & (/YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_DYN%YYTLSCAW%NDIM,YDGEOMETRY%YRDIM%NGPBLKS/))) LLREALLOC = .TRUE.
   IF (.not.ALL(IDX_ZSLBUF1 == &
        & (/YDMODEL%YRML_DYN%YRSL%NASLB1,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1/))) LLREALLOC = .TRUE.
   IF (.not.ALL(IDX_ZSLBUF2 == &
        & (/YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB2%NFLDSLB2,YDGEOMETRY%YRDIM%NGPBLKS/))) LLREALLOC = .TRUE.
   IF (LLREALLOC) THEN
      DEALLOCATE(ZSLBUF2)
      DEALLOCATE(ZSLBUF1)
      DEALLOCATE(ZLSCAW)
      DEALLOCATE(IL0)
   ENDIF
ELSE ! (.not.LLINITALLOC) -- occurs only once
   LLINITALLOC = .TRUE.
   LLREALLOC = .TRUE.
ENDIF

IF (LLREALLOC) THEN
   IDX_IL0 = (/YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIM%NGPBLKS/)
   IDX_ZLSCAW = (/YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_DYN%YYTLSCAW%NDIM,YDGEOMETRY%YRDIM%NGPBLKS/)
   IDX_ZSLBUF1 = (/YDMODEL%YRML_DYN%YRSL%NASLB1,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1/)
   IDX_ZSLBUF2 = (/YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB2%NFLDSLB2,YDGEOMETRY%YRDIM%NGPBLKS/)

   ALLOCATE(IL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3,YDGEOMETRY%YRDIM%NGPBLKS))
   ALLOCATE(ZLSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_DYN%YYTLSCAW%NDIM,YDGEOMETRY%YRDIM%NGPBLKS))
   ALLOCATE(ZSLBUF1(YDMODEL%YRML_DYN%YRSL%NASLB1,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1))
   ALLOCATE(ZSLBUF2(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB2%NFLDSLB2,YDGEOMETRY%YRDIM%NGPBLKS))
ENDIF

CALL GP_MODEL(YDGEOMETRY,YDFIELDS,YDMODEL,CDCONF,LD_DFISTEP,IL0,ZLSCAW,ZSLBUF1,ZSLBUF2,PTRAJEC=PTRAJEC,PTRAJEC_OOPS=PTRAJEC_OOPS)

IF (LHOOK) CALL DR_HOOK('GP_MODEL_HEAP2',1,ZHOOK_HANDLE)

END SUBROUTINE GP_MODEL_HEAP2
