! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_COMM_SPLIT_MOD

!**** *MPL_COMM_SPLIT* - Split a communicator

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE*
!      Original : 11-09-2012


USE EC_PARKIND, ONLY : JPIM 
USE MPL_MPIF, ONLY : MPI_UNDEFINED

IMPLICIT NONE

PRIVATE
PUBLIC MPL_COMM_SPLIT

CONTAINS 

SUBROUTINE MPL_COMM_SPLIT (KCOMM, KCOLOR, KKEY, KNEWCOMM, KERROR, CDSTRING)
INTEGER (KIND=JPIM), INTENT (IN)  :: KCOMM
INTEGER (KIND=JPIM), INTENT (IN)  :: KCOLOR
INTEGER (KIND=JPIM), INTENT (IN)  :: KKEY
INTEGER (KIND=JPIM), INTENT (OUT) :: KNEWCOMM
INTEGER (KIND=JPIM), INTENT (OUT) :: KERROR
CHARACTER (LEN=*), OPTIONAL, INTENT (IN) :: CDSTRING
INTEGER (KIND=JPIM) :: ICOLOR

ICOLOR=KCOLOR
IF(ICOLOR<0) ICOLOR=MPI_UNDEFINED

CALL MPI_COMM_SPLIT (KCOMM, ICOLOR, KKEY, KNEWCOMM, KERROR)

END SUBROUTINE MPL_COMM_SPLIT

END MODULE MPL_COMM_SPLIT_MOD
