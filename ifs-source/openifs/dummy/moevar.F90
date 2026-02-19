! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE MOEVAR(YDRIP,YDML_LBC,KLSGTS)
use parkind1 , only:&
 & jpim
USE YOMRIP , ONLY : TRIP
USE YEMLBC_MODEL,ONLY : TELBC_MODEL
TYPE(TRIP) ,INTENT(INOUT) :: YDRIP
TYPE(TELBC_MODEL) ,INTENT(IN) :: YDML_LBC
INTEGER(KIND=JPIM),INTENT(OUT) :: KLSGTS(0:YDRIP%NSTOP/YDML_LBC%NFRLSG)
call abor1("moevar.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE MOEVAR
