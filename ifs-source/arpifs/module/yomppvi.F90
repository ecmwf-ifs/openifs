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

MODULE YOMPPVI

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!     Variables for vertical interpolator used for example in FULL-POS or
!     in the observation vertical interpolator.

!     LESCALE and LESCALE_[X]: use of ESCALE system in APACHE.
!     LRPPUV_CSTEXT: use of constant top extrapolation in PPUV
!     LRPPUV_CALLITPQ: call to PPITPQ required in PPUV
!     LPPVIVX: if true, limitation on maximum wind velocity 
!     RPPVIVX: maximum wind velocity (m/s)
!     RPPVIVP: pressure threshold (Pa) to limit maximum wind velocity
!     LNOTS_T: alternate lower vertical extrapolation for temperature, which does not use T_star (output of CTSTAR).
!     ------------------------------------------------------------------

LOGICAL :: LESCALE
LOGICAL :: LESCALE_T
LOGICAL :: LESCALE_Q
LOGICAL :: LESCALE_U
LOGICAL :: LESCALE_PD
LOGICAL :: LESCALE_GFL
LOGICAL :: LRPPUV_CSTEXT
LOGICAL :: LRPPUV_CALLITPQ
LOGICAL :: LPPVIVX
LOGICAL :: LNOTS_T
REAL(KIND=JPRB) :: RPPVIVX
REAL(KIND=JPRB) :: RPPVIVP

!     ------------------------------------------------------------------
END MODULE YOMPPVI
