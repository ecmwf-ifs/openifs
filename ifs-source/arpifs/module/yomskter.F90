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

MODULE YOMSKTER

USE PARKIND1,ONLY : JPIM  , JPRB
USE PARDIM,  ONLY : JPMXLE, JPMXGL
IMPLICIT NONE

SAVE

!*     YOMSKTER - SKT FIRST GUESS ERRORS AND RELATED PARAMETERS

INTEGER(KIND=JPIM) :: JPMAXVAR
PARAMETER (JPMAXVAR=200)

TYPE SKTER
  INTEGER(KIND=JPIM)                           :: MNROWS,MNVARS,MNFLDS,MNPTE,MNLEVMX
  INTEGER(KIND=JPIM), DIMENSION(JPMXGL)        :: MLONE
  INTEGER(KIND=JPIM), DIMENSION(JPMAXVAR)      :: MFIELDS,MNLEVS
  REAL(KIND=JPRB), DIMENSION(JPMXGL)           :: RZLATE,RZLON0E,RZDLONE
  REAL(KIND=JPRB), DIMENSION(JPMXLE,JPMAXVAR)  :: RZA,RZB
  REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE :: RZESKTGRID
END TYPE

TYPE(SKTER) :: YGSKTER

!-----------------------------------------------------------------------

END MODULE YOMSKTER
