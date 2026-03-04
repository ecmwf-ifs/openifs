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

MODULE YOETLDIAG

! Module characterizing parameters for different tangent-linear (TL) 
! and adjoint (AD) diagnostics

!   NSIM4DWR   : value of NSIM4D to write out increments from 4D-Var
!   LINCRPR    : .TRUE. when increments should be stored unpacked and with
!                64-bit precision to fully re-construct TL/AD run from 
!                the saved files
!   LINCTV     : .TRUE. when T increments to be saved directly as Tv used
!                with grid-point q
!   LEDIAGTL   : .TRUE. when writting/setting some fields or writting norms
!                for TL diagnostics
!   LTL4DREP   : .TRUE. when the model integration of 4D-Var should be
!                exactly reproduced

! -------------------------------------------------------------------------------

USE PARKIND1, ONLY: JPIM, JPRB

IMPLICIT NONE

SAVE

INTEGER (KIND=JPIM) :: NSIM4DWR

LOGICAL :: LINCRPR
LOGICAL :: LINCTV
LOGICAL :: LEDIAGTL
LOGICAL :: LTL4DREP

! ---------------------------------------------------------------------

END MODULE YOETLDIAG
