MODULE YOERIP
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

USE PARKIND1  ,ONLY : JPIM     ,JPRB,JPRD

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Real time related variables for ECMWF radiation: updated UPDTIER

!     RIP0M  : I0 WEIGHTED BY THE DISTANCE EARTH-SUN

!     RCODECM: COSINE OF THE DECLINATION
!     RSIDECM:   SINE OF THE DECLINATION

!     RCOVSRM: COSINE OF TRUE SOLAR TIME
!     RSIVSRM:   SINE OF TRUE SOLAR TIME

REAL(KIND=JPRB) :: RIP0M
REAL(KIND=JPRD) :: RCODECM
REAL(KIND=JPRD) :: RSIDECM
REAL(KIND=JPRD) :: RCOVSRM
REAL(KIND=JPRD) :: RSIVSRM
!     ------------------------------------------------------------------
END MODULE YOERIP
