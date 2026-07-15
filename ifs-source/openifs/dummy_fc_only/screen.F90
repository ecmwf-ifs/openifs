! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SCREEN(YDMODEL,YDODB,YDGOM5,YDJOT)
USE TYPE_MODEL, ONLY : MODEL
USE JO_TABLE_MOD, ONLY : JO_TABLE
USE SUPERGOM_CLASS, ONLY : CLASS_SUPERGOM
USE DBASE_MOD, ONLY : DBASE
TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
CLASS(DBASE) ,INTENT(INOUT) :: YDODB
TYPE(CLASS_SUPERGOM) ,INTENT(IN) :: YDGOM5
TYPE(JO_TABLE) ,INTENT(INOUT) :: YDJOT
call abor1("SCREEN should never be called with this build configuration - EXITING")
END SUBROUTINE SCREEN
