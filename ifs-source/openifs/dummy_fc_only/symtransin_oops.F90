! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SYMTRANSIN_OOPS(YDGEOM,YDMODEL,YDGMV,YDGFL,YGFL,YDGFLT5Q)
USE GEOMETRY_MOD , ONLY : GEOMETRY
use parkind1 , only:&
 & jprb
USE YOMGFL , ONLY : TGFL
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOMGMV , ONLY : TGMV
USE TYPE_MODEL , ONLY : MODEL
USE FJBCHVAR_MOD
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOM
TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
TYPE(TGMV) ,INTENT(IN) :: YDGMV
TYPE(TGFL) ,INTENT(INOUT) :: YDGFL
TYPE(TYPE_GFLD) ,INTENT(IN) :: YGFL
REAL(KIND=JPRB) ,INTENT(IN) :: YDGFLT5Q(YDGEOM%YRDIM%NPROMA,YDGEOM%YRDIMV%NFLEVG,YDGEOM%YRDIM%NGPBLKS)
call abor1("SYMTRANSIN_OOPS should never be called with this build configuration - EXITING")
END SUBROUTINE SYMTRANSIN_OOPS
