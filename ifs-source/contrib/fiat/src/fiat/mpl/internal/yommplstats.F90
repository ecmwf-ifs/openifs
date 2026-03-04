! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE YOMMPLSTATS

USE EC_PARKIND  ,ONLY : JPRD, JPIM

IMPLICIT NONE

SAVE

PRIVATE :: JPRD, JPIM
PUBLIC

!     ------------------------------------------------------------------
! Module for communications statistics.
! Module is internal to the MPLSTATS package -
! routines MPL_SENDSTATS, MPL_RECVSTATS

! LMPLSTATS - TRUE for gathering communications statistics


LOGICAL :: LMPLSTATS = .FALSE.
REAL(KIND=JPRD), ALLOCATABLE    :: MPLSENDBYTES(:), MPLRECVBYTES(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: MPLSENDNUM(:), MPLRECVNUM(:)

END MODULE YOMMPLSTATS




