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

SUBROUTINE SUSPGPG(YDGEOMETRY,YDGFL,YGFL,PFPDAT,KPARAM)

!**** *SUSPGPG*  - Initialize some upper air gridpoint fields from
!                  spectral GRIB.

!     Purpose.
!     --------
!           Initialize the some upper air gridpoint fields of the model
!           from spectral GRIB file(s)

!**   Interface.
!     ----------
!        *CALL* *SUSPGPG(.....)*

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.  SPEREE - spectral to grid-point transform
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud *ECMWF*
!      Original : 94-06-13

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad (Jan 2010): remove useless variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGFL       , ONLY : TGFL
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOM_YGFL     , ONLY : TYPE_GFLD

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL)        ,INTENT(INOUT) :: YDGFL
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YGFL
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPDAT(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
INTEGER(KIND=JPIM),INTENT(IN)    :: KPARAM 
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) ::  IPT,JGFL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "inv_trans.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUSPGPG',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & NPROMA=>YDDIM%NPROMA, NRESOL=>YDDIM%NRESOL, NSPEC2=>YDDIM%NSPEC2, &
 & NFLEVL=>YDDIMV%NFLEVL, &
 & GFL=>YDGFL%GFL, &
 & NBSETLEV=>YDMP%NBSETLEV)
!     ------------------------------------------------------------------

!*       1.  CONVERT TO GRIDPOINT SPACE.
!            ---------------------------

DO JGFL=1,NUMFLDS
  IF(KPARAM == YCOMP(JGFL)%IGRBCODE) THEN
    IPT=YCOMP(JGFL)%MP
    EXIT
  ENDIF
ENDDO
CALL INV_TRANS(PSPSCALAR=PFPDAT,KRESOL=NRESOL,KPROMA=NPROMA,KVSETSC=NBSETLEV,&
 & PGP=GFL(:,:,IPT,:))  

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSPGPG',1,ZHOOK_HANDLE)
END SUBROUTINE SUSPGPG

