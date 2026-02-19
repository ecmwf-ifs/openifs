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

MODULE YOMMSE

USE PARKIND1, ONLY : JPIM, JPRB

!**** *YOMMSE*  - Store MSE parameters

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE*
!      Original : 01-01-2011

!      Modifications :
!      P.Marguinaud : 11-09-2012 : New parameters


IMPLICIT NONE

TYPE TMSE

LOGICAL :: LPGDFWR   = .FALSE.       ! Default is do not write PGD 
LOGICAL :: LHISFWR   = .TRUE.        ! Default is write SURFEX historic data 
LOGICAL :: LFTZERO   = .TRUE.        ! Force SURFEX fields to 0 in extension zone

INTEGER(KIND=JPIM) :: NSURFEXCTLMAX = 3_JPIM
INTEGER(KIND=JPIM) :: NSURFEXCTL = 3_JPIM

CHARACTER(LEN=4) :: CPREFIX_SFX1 = 'SFX.', &
                  & CPREFIX_SFX2 = 'X'

REAL (KIND=JPRB) :: XZSEPS = 1.E-2_JPRB ! Admitted relative difference between model and SURFEX orography

END TYPE TMSE

!! TYPE(TMSE), POINTER :: YRMSE => NULL()

SAVE

END MODULE YOMMSE

