! (C) Copyright 2014- ECMWF.
! (C) Copyright 2014- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOMRIP0

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ---------------------------------------------------------------------------
!*    Real time related variables: variables, the final value is set-up in SURIP0
!     Values are identical for all models run under the OOPS layer.
!     ---------------------------------------------------------------------------

!     NINDAT : run initial date in the form AAAAMMDD
!     NSSSSS : initial time in seconds (e.g. for 12h, 43200)
!     RTIMST : ABSOLUTE TIME OF THE MODEL AT START

INTEGER(KIND=JPIM) :: NINDAT
INTEGER(KIND=JPIM) :: NSSSSS
REAL(KIND=JPRB) :: RTIMST

!     ------------------------------------------------------------------
!     LASTRF : Keep insolation as for years around 2000, used to prevent any drift
!              due to the formulation of RET which is not correct unless over the period 1980-2020

LOGICAL         :: LASTRF

END MODULE YOMRIP0
