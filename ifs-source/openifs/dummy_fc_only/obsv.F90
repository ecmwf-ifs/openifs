! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE OBSV(YDEPHY,YDML_PHY_MF,YDJOT,YDVARBC,YDTCV5,YDGOM5,YDODB,KSLOT,CDOBJOB)
USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE YOEPHY , ONLY : TEPHY
use parkind1 , only:&
 & jpim
USE JO_TABLE_MOD , ONLY : JO_TABLE
USE TOVSCV_MOD , ONLY : TOVSCV
USE VARBC_CLASS , ONLY : CLASS_VARBC
USE SUPERGOM_CLASS , ONLY : CLASS_SUPERGOM
USE DBASE_MOD , ONLY : DBASE
TYPE(TEPHY) , INTENT(INOUT) :: YDEPHY
TYPE(MODEL_PHYSICS_MF_TYPE), INTENT(INOUT) :: YDML_PHY_MF
TYPE(JO_TABLE) , INTENT(INOUT) :: YDJOT
TYPE(CLASS_VARBC) , INTENT(INOUT) :: YDVARBC
TYPE(TOVSCV) , INTENT(IN) :: YDTCV5
TYPE(CLASS_SUPERGOM) , INTENT(IN) :: YDGOM5
CLASS(DBASE) , INTENT(INOUT) :: YDODB
INTEGER(KIND=JPIM) , INTENT(IN) :: KSLOT
CHARACTER(LEN=2) , INTENT(IN) :: CDOBJOB
call abor1("OBSV should never be called with this build configuration - EXITING")
END SUBROUTINE OBSV
