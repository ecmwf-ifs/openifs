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

MODULE YOMFPOBJ

USE PARKIND1  ,ONLY : JPIM
USE PARFPOS  , ONLY : JPOSOBJ

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
! NFPOBJ   : number of fullpos objects
! NFPCONF  : configuration of each object
! CNAMELIST : namelist file for each object

! LFPDISPLAY_PARAMETERS : to display hard-coded parameters

INTEGER(KIND=JPIM) :: NFPOBJ=1
INTEGER(KIND=JPIM) :: NFPCONF(JPOSOBJ)
CHARACTER(LEN=256) :: CNAMELIST(JPOSOBJ)
LOGICAL :: LFPDISPLAY_PARAMETERS=.FALSE.

!     ------------------------------------------------------------------
END MODULE YOMFPOBJ
