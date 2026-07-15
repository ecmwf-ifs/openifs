! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUPERTPAR(YDML_PHY_MF,YDECUMF,YDERAD)
USE RANDOM_NUMBERS_MIX
USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE YOECUMF , ONLY : TECUMF
USE YOERAD , ONLY : TERAD
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT), TARGET :: YDML_PHY_MF
TYPE(TECUMF) ,INTENT(INOUT) :: YDECUMF
TYPE(TERAD) ,INTENT(INOUT) :: YDERAD
call abor1("SUPERTPAR should never be called with this build configuration - EXITING")
END SUBROUTINE SUPERTPAR
