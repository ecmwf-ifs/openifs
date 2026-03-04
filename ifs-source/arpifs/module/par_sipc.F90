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

MODULE PAR_SIPC

! -- jpbyteint : number of bytes per integer

! -- jpbyterea : number of bytes per real

! -- jpbytecha : number of bytes per character

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE
SAVE
INTEGER(KIND=JPIM), PARAMETER :: JPBYTEINT=4
INTEGER(KIND=JPIM), PARAMETER :: JPBYTEREA=8
INTEGER(KIND=JPIM), PARAMETER :: JPBYTECHA=1

! -- jptest :  The models will test during 2*jptest seconds if the file
!              DUMMY_SIPC has been created by OASIS, signaling that the
!              SHM pools are opened. After, they will abort.

INTEGER(KIND=JPIM), PARAMETER ::  JPTEST=100

END MODULE PAR_SIPC

