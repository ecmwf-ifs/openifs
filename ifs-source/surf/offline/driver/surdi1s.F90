SUBROUTINE SURDI1S
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOERDI1S , ONLY : REMISS

#ifdef DOC
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SURDI1s* - Initialize common YOERDI1S
!     ------------------------------------------------------------------
#endif
IMPLICIT NONE

REMISS = 0.99_JPRB

RETURN
END SUBROUTINE SURDI1S
