! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUENSCOV(YDGEOMETRY,YDFIELDS,YDMODEL)
 USE TYPE_MODEL , ONLY : MODEL
 USE FIELDS_MOD , ONLY : FIELDS
 USE YOMENSCOV
 USE GEOMETRY_MOD , ONLY : GEOMETRY
 TYPE(GEOMETRY) ,INTENT(INOUT) :: YDGEOMETRY
 TYPE(FIELDS) ,INTENT(INOUT) :: YDFIELDS
 TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
call abor1("suenscov.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUENSCOV
