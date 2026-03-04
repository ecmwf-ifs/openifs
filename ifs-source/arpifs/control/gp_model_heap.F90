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

SUBROUTINE GP_MODEL_HEAP(YDGEOMETRY,YDFIELDS,YDMODEL,CDCONF,LD_DFISTEP,PTRAJEC,PTRAJEC_OOPS)
!****-------------------------------------------------------------------
!**** *GP_MODEL_HEAP* - Grid-point model
!****-------------------------------------------------------------------
!     Purpose.   gp_model driver using memory on heap for large arrays
!     --------

!**   Interface.
!     ----------
!        *CALL* *GP_MODEL_HEAP (..)

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
!        R. El Khatib *Meteo-France* 02-Oct-2014

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
INTEGER(KIND=JPIM), ALLOCATABLE :: IL0(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZLSCAW(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSLBUF1(:,:), ZSLBUF2(:,:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gp_model.intfb.h"

IF (LHOOK) CALL DR_HOOK('GP_MODEL_HEAP',0,ZHOOK_HANDLE)

ALLOCATE(IL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3,YDGEOMETRY%YRDIM%NGPBLKS))
ALLOCATE(ZLSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_DYN%YYTLSCAW%NDIM,YDGEOMETRY%YRDIM%NGPBLKS))
ALLOCATE(ZSLBUF1(YDMODEL%YRML_DYN%YRSL%NASLB1,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1))
ALLOCATE(ZSLBUF2(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB2%NFLDSLB2,YDGEOMETRY%YRDIM%NGPBLKS))
CALL GP_MODEL(YDGEOMETRY,YDFIELDS,YDMODEL,CDCONF,LD_DFISTEP,IL0,ZLSCAW,ZSLBUF1,ZSLBUF2,PTRAJEC=PTRAJEC,PTRAJEC_OOPS=PTRAJEC_OOPS)
DEALLOCATE(ZSLBUF2)
DEALLOCATE(ZSLBUF1)
DEALLOCATE(ZLSCAW)
DEALLOCATE(IL0)

IF (LHOOK) CALL DR_HOOK('GP_MODEL_HEAP',1,ZHOOK_HANDLE)

END SUBROUTINE GP_MODEL_HEAP
