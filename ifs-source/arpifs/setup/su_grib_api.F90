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

SUBROUTINE SU_GRIB_API(YDGEOMETRY,YDVAB,LDMCC04,YDGBH,YDFPUSERGEO)

!**** *SU_GRIB_API* - Routine to intitialize GRIB API usage

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SU_GRIB_API*

!        Explicit arguments :  None.
!        --------------------

!        Implicit arguments :  
!        --------------------

!     Method.  
!     -------  

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 15 May 2007

!     ------------------------------------------------------------------

USE TYPE_FPUSERGEO, ONLY : TFPUSERGEO
USE GRIB_HANDLES_MOD  , ONLY : TYPE_GRIB_HANDLES
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMVERT  , ONLY : TVAB
USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE(GEOMETRY),          INTENT(IN), TARGET :: YDGEOMETRY
TYPE(TVAB),              INTENT(IN)         :: YDVAB
LOGICAL,                 INTENT(IN)         :: LDMCC04
TYPE(TYPE_GRIB_HANDLES), INTENT(INOUT)      :: YDGBH
TYPE (TFPUSERGEO),       INTENT(IN), OPTIONAL :: YDFPUSERGEO

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "preset_grib_template.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SU_GRIB_API',0,ZHOOK_HANDLE)

CALL GSTATS(1977,0)

!!YDGBH%GEOM => YDGEOMETRY
IF(PRESENT(YDFPUSERGEO)) THEN
  CALL PRESET_GRIB_TEMPLATE(YDGEOMETRY,YDVAB,LDMCC04,'GG',YDGBH,YDFPUSERGEO=YDFPUSERGEO)
  CALL PRESET_GRIB_TEMPLATE(YDGEOMETRY,YDVAB,LDMCC04,'GG2',YDGBH,YDFPUSERGEO=YDFPUSERGEO)
  CALL PRESET_GRIB_TEMPLATE(YDGEOMETRY,YDVAB,LDMCC04,'GG_ML',YDGBH,YDFPUSERGEO=YDFPUSERGEO)
  CALL PRESET_GRIB_TEMPLATE(YDGEOMETRY,YDVAB,LDMCC04,'SH',YDGBH,YDFPUSERGEO=YDFPUSERGEO)
  CALL PRESET_GRIB_TEMPLATE(YDGEOMETRY,YDVAB,LDMCC04,'SH_ML',YDGBH,YDFPUSERGEO=YDFPUSERGEO)
  !CALL PRESET_GRIB_TEMPLATE(YDGEOMETRY,YDVAB,LDMCC04,'DIAG',YDGBH,YDFPUSERGEO=YDFPUSERGEO)
ELSE  
  CALL PRESET_GRIB_TEMPLATE(YDGEOMETRY,YDVAB,LDMCC04,'GG',YDGBH)
  CALL PRESET_GRIB_TEMPLATE(YDGEOMETRY,YDVAB,LDMCC04,'GG2',YDGBH)
  CALL PRESET_GRIB_TEMPLATE(YDGEOMETRY,YDVAB,LDMCC04,'GG_ML',YDGBH)
  CALL PRESET_GRIB_TEMPLATE(YDGEOMETRY,YDVAB,LDMCC04,'SH',YDGBH)
  CALL PRESET_GRIB_TEMPLATE(YDGEOMETRY,YDVAB,LDMCC04,'SH_ML',YDGBH)
  !CALL PRESET_GRIB_TEMPLATE(YDGEOMETRY,YDVAB,LDMCC04,'DIAG',YDGBH)
ENDIF  
CALL GSTATS(1977,1)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SU_GRIB_API',1,ZHOOK_HANDLE)
END SUBROUTINE SU_GRIB_API
