! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUINIF_FP(YDFPS,YDGEOMETRY,YDSPEC,YDGFL,YDSURF,YDCFU,YDXFU,YDMODEL,KSTEP)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE TYPE_MODEL , ONLY : MODEL
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMGFL , ONLY : TGFL
USE YOMCFU , ONLY : TCFU
USE YOMXFU , ONLY : TXFU
use parkind1 , only:&
 & jpim
USE YOMFP_SERV , ONLY : FP_SERV
use spectral_fields_mod, only:&
 & spectral_field
TYPE (FP_SERV) , INTENT (INOUT) :: YDFPS
TYPE (GEOMETRY) , INTENT (INOUT) :: YDGEOMETRY
TYPE(SPECTRAL_FIELD), INTENT (INOUT) :: YDSPEC
TYPE (TGFL) , INTENT (INOUT) :: YDGFL
TYPE (TSURF) , INTENT (INOUT) :: YDSURF
TYPE (TCFU) , INTENT (INOUT) :: YDCFU
TYPE (TXFU) , INTENT (INOUT) :: YDXFU
TYPE(MODEL) , INTENT (INOUT) :: YDMODEL
INTEGER (KIND=JPIM) , INTENT (INOUT) :: KSTEP
call abor1("suinif_fp.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUINIF_FP
