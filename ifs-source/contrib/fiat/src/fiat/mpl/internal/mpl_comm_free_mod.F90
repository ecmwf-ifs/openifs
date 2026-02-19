! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_COMM_FREE_MOD

!**** *MPL_COMM_FREE_MOD*  - Release ressources used by a communicator

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE*
!      Original : 11-09-2012

USE EC_PARKIND, ONLY : JPIM

IMPLICIT NONE

PRIVATE
PUBLIC :: MPL_COMM_FREE

CONTAINS

SUBROUTINE MPL_COMM_FREE (KCOMM, KERR, CDSTRING)
INTEGER (KIND=JPIM), INTENT (IN)  :: KCOMM
INTEGER (KIND=JPIM), INTENT (OUT) :: KERR
CHARACTER (LEN=*),   INTENT (IN), OPTIONAL :: CDSTRING

CALL MPI_COMM_FREE (KCOMM, KERR)

END SUBROUTINE MPL_COMM_FREE

END MODULE MPL_COMM_FREE_MOD
