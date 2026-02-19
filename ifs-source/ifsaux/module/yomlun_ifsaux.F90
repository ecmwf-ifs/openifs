! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOMLUN_IFSAUX

USE EC_LUN  ,ONLY : NULOUT, NULERR

IMPLICIT NONE

SAVE
PRIVATE
PUBLIC :: NULOUT, NULERR

!     ------------------------------------------------------------------

!*    Logical units used by code

!     NULOUT :   output unit
!     NULERR :   unit number for comparison with reference run

!     ------------------------------------------------------------------
END MODULE YOMLUN_IFSAUX
