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

MODULE YOMMTS

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

! LMTS : Switch d'activation de la production de temperatures de brillance
! IF LMTS and LMTSCL production of clear sky TB
! NHLIM : max hour  (production for hour of forecast < -NHLIM)
! LUBIQUITAIRE : To have the satellite at the vetical of all the grid points
! LISP_HYBRID : MSG simulated data in 'real' conditions over the MSG domain, 
!               ubiquitaire computation elsewhere

INTEGER(KIND=JPIM) :: NTYPE

LOGICAL :: LMTS,LMTSCL,LUBIQUITAIRE,LISP_HYBRID


END MODULE YOMMTS
