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

MODULE TYPE_FAOPH

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!========== HANDLING OF OUTPUT FIELDS IN FA FILES  ======

! CFPCA : names of output FA frames
! NHEADMAX : maximum size of the header of a record
! NGPSIZPK : size of a gridpoint field in FA file
! NSPSIZPK : size of a spectral field in FA file

TYPE TFAOPH

CHARACTER(LEN=16)  :: CFPCA
INTEGER(KIND=JPIM) :: NHEADMAX
INTEGER(KIND=JPIM) :: NGPSIZPK
INTEGER(KIND=JPIM) :: NSPSIZPK

END TYPE TFAOPH

!     ------------------------------------------------------------------
END MODULE TYPE_FAOPH
