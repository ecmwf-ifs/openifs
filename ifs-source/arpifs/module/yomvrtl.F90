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

MODULE YOMVRTL

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Switches for variational assimilation: use of tangent linear model
! L131TL  : .T. = the incremental is to be run with the tangent model.
! LTLINT  : .T. = Flag telling if we are in the tangent integration
!                 of the model (incremental)
! LOBSTL  : .T. = we use the tangent linear of the observation operators
!                 (otherwise a finite-difference approximation is used)
LOGICAL :: L131TL
LOGICAL :: LTLINT
LOGICAL :: LOBSTL

!*    Alterations of the TL and AD models
LOGICAL :: LDRYTL       ! .T. = remove the coupling of q and T in the dynamics

!*    Global arrays and variables for the truncated Newton algorithm

LOGICAL :: LIDMODEL     ! .T. = Replace tl-model by identity operator

!     ------------------------------------------------------------------
END MODULE YOMVRTL
