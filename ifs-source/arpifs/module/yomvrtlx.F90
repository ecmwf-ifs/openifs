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

MODULE YOMVRTLX

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Switches for variational assimilation: use of tangent linear model
! L801TL  : .T. = the sensitivity  is to be run with the tangent model.
!                 (otherwise a finite-difference approximation is used)
LOGICAL :: L801TL=.TRUE.
LOGICAL :: LMINI=.TRUE.

!     ------------------------------------------------------------------
END MODULE YOMVRTLX
