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

MODULE YOMRES

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Control variables for the restart control

! NFRRES  : frequency of restart write_ups (time-steps)
! NRESTS  : array containing restart write-up steps

! NFLREST :  flag set non-zero by signal to indicate that
!            restart files should be written

! NFLSTOP :  flag set non-zero by signal to indicate that
!            execution should stop after writing restart files

! LSDHM   :  .TRUE. if time stamp is written as 'ddddhhmm' ;
!            otherwise time stamp is controled by LINC

! LDELRES :  .TRUE. if old restart files are deleted

INTEGER(KIND=JPIM), PARAMETER :: JPNWST=40
INTEGER(KIND=JPIM), PARAMETER :: JPNRF=4
INTEGER(KIND=JPIM) :: NRESTS(0:JPNWST)

CHARACTER (LEN = 8)  :: CRFTIME(JPNRF)
CHARACTER (LEN = 14) :: CTIME
CHARACTER (LEN = 8)  :: CSTEP
INTEGER(KIND=JPIM) :: N1RFS
INTEGER(KIND=JPIM) :: N2RFS
INTEGER(KIND=JPIM) :: NFRRES
INTEGER(KIND=JPIM) :: NFLSTOP
INTEGER(KIND=JPIM) :: NFLREST
LOGICAL :: LSDHM
LOGICAL :: LDELRES

END MODULE YOMRES
