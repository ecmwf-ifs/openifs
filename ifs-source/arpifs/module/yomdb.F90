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

MODULE YOMDB

USE PARKIND1  ,ONLY : JPIM, JPRD

USE PARCMA, ONLY : JPMXUP, JPMXENKF, JPMXENDA, JPMXFCDIAG, JPXTSL, JPMX_LIMB_TAN, &
  & JPMX_RADAR_NIV, JPMX_AK, JP_NUMEV
USE OML_MOD , ONLY : OML_MY_THREAD
USE YOMDB_CONSTANTS

IMPLICIT NONE
SAVE
PUBLIC

PRIVATE :: JPIM, JPRD

LOGICAL :: LODB = .TRUE.
INTEGER(KIND=JPIM), PARAMETER :: JPMAXTSLDB = JPXTSL
INTEGER(KIND=JPIM), PRIVATE   :: I

!*** Disable UNMAPDB-function
LOGICAL :: LUNMAPDB_DO = .TRUE.

!***  Database handles
INTEGER(KIND=JPIM) :: NUMSIMDB  =  0   ! No. of presently opened DBs
INTEGER(KIND=JPIM) :: NHANDLEDB(JPMAXSIMDB) = (/( 0 , I=1,JPMAXSIMDB)/)
INTEGER(KIND=JPIM) :: NPOOLSDB(JPMAXSIMDB)  = (/( 0 , I=1,JPMAXSIMDB)/)
CHARACTER(LEN= 20) :: COPENDB(JPMAXSIMDB)   = (/(' ', i=1,JPMAXSIMDB)/)
CHARACTER(LEN= 20) :: CSTATUSDB(JPMAXSIMDB) = (/(' ', i=1,JPMAXSIMDB)/)

!***  For each thread
INTEGER(KIND=JPIM), ALLOCATABLE :: NACTIVEDB(:) ! Index to NHANDLEDB() for active DB
CHARACTER(LEN= 20), ALLOCATABLE :: CACTIVEDB(:) ! Say 'ECMA' or 'CCMA'

!***  Miscellaneous
LOGICAL :: L_DEBUG   = .FALSE.
LOGICAL :: L_SWAPOUT = .FALSE.
LOGICAL :: L_CANARI  = .FALSE.
LOGICAL :: L_CANALT  = .FALSE.
LOGICAL :: L_SETACTIVE   = .FALSE.   ! when true it will set data to active if true_obsvalue!=obsvalue (for averaging) 
                                     ! This allows to use obsort SQLs request to shuffle from ECMA to CCMA

TYPE ROBAUX_T
SEQUENCE
REAL(KIND=JPRD), POINTER :: ROBAUX(:,:) => NULL()
END TYPE ROBAUX_T

TYPE YOMDB_T
SEQUENCE
#include "yomdb_vars.h"
END TYPE YOMDB_T

TYPE(YOMDB_T), ALLOCATABLE :: O_(:)

END MODULE YOMDB
