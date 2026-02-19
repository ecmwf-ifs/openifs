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

MODULE YOMFP4L

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!     Descriptors of a post-processing request

! === ALL FIELDS ===

!     NFIELDG : number of fields
!     ICOD    : internal field codes
!     ZLEV    : level values
!     LLSURF  : Surface / upper-air 
!     CLPREF  : Kind of level / FA prefix name of each field
!     CLNAME  : FA suffix name of each field
!     IGRIB   : GRIB code of each field

! === GRID POINT FIELDS ONLY ===

!     IDOM : number of subdomains for each field
!     IDMP : indexes of subdomains for each field

! === SPECTRAL FIELDS ONLY ===

!     NFIELDL : local number of fields in the "Vertical" distribution 
!     IVSETG  : V-set for each field (global)
!     IVLOCG  : Local field adress for each field


TYPE TRQFP

! All fields
INTEGER(KIND=JPIM) :: NFIELDG = 0
INTEGER(KIND=JPIM), ALLOCATABLE :: ICOD(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IGRIB(:)
LOGICAL,            ALLOCATABLE :: LLSURF(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IVEC(:)
REAL(KIND=JPRB),    ALLOCATABLE :: ZLEV(:)
CHARACTER(LEN=12),  ALLOCATABLE :: CLNAME(:)
CHARACTER(LEN=8),   ALLOCATABLE :: CLPREF(:)

! Gridpoint fields only (however ...)
INTEGER(KIND=JPIM), ALLOCATABLE :: IDOM(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IDMP(:,:)

! Spectral fields only

INTEGER(KIND=JPIM) :: NFIELDL = 0
INTEGER(KIND=JPIM), ALLOCATABLE :: IVSETG(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IVLOCG(:)

END TYPE TRQFP


CONTAINS

INTEGER(KIND=JPIM) FUNCTION IFPSEARCH(YDTYPE,KCOD)

! Search for KCOD in YDTYPE%ICOD. Return 0 if not found.

TYPE(TRQFP), INTENT(IN) :: YDTYPE
INTEGER(KIND=JPIM), INTENT(IN) :: KCOD

INTEGER(KIND=JPIM) :: J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMFP4L:IFPSEARCH',0,ZHOOK_HANDLE)

IFPSEARCH=0

DO J=1,YDTYPE%NFIELDG
  IF (YDTYPE%ICOD(J)==KCOD) THEN
    IFPSEARCH=J
    EXIT
  ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('YOMFP4L:IFPSEARCH',1,ZHOOK_HANDLE)

END FUNCTION IFPSEARCH

!     ------------------------------------------------------------------
END MODULE YOMFP4L
