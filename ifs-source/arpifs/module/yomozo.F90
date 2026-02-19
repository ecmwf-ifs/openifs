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

MODULE YOMOZO

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

TYPE :: TOZO

!     GRID POINT ARRAYS FOR OZONE: DM-GLOBAL VERSIONS

!     REAL TOZ2DG(NFLEVG*NVCLIS*NTOZ2D,NDGSAG:NDGENG)
REAL(KIND=JPRB),ALLOCATABLE:: TOZ2DG(:,:)

!     GRID POINT ARRAYS FOR OZONE: DM-LOCAL VERSIONS

REAL(KIND=JPRB),ALLOCATABLE:: TOZ2DL(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: TOZ3DBL(:,:)

!     TEECL2: EQUIVALENT CHLORINE CONTENT**2
REAL(KIND=JPRB) ::  TEECL2

!     THE FIELDS FOR OZONE PARAMETERIZATION MAY BE STORED IN A FULL 3D BUFFER
!     OR ONLY IN A LATITUDExLEVEL ARRAY (VALID ONLY IF NSTTYP = 1 )

END TYPE TOZO

!!TYPE(TOZO), POINTER :: YROZO => NULL()

END MODULE YOMOZO
