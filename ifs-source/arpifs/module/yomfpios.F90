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

MODULE YOMFPIOS

USE PARKIND1  ,ONLY : JPIM
USE YOMCT0    ,ONLY : JPNPST

IMPLICIT NONE

SAVE

! General post-processing I/O parameters

! NFPGRIB   : level of GRIB encoding
! NFPWRITE  : write (1 or 2) or not (0 or -1) output files. If (-1 or 2) then the norms are not written.
! NFPXFLD   : maximum number of fields to be extracted from a pp buffer at a time
! NFPDIGITS : number of digits for time stamp extension
!             If negative the NFPDIGITS is used as abs(NFPDIGITS) without
!             the character "+" preceeding it. This facility enables to
!             make file names like boundary files in some sense.
! NFPDITS    : array containing postprocessing steps
! NFPDITSMIN : array containing postprocessing steps in minutes for sub-hour outputs
! NFRFPDI    : frequency of post-processing events
! LPGDFWR   : write Surfex PGD file or not
! LHISFWR   : write SURFEX historic data or not

TYPE TNAMFPIOS

INTEGER(KIND=JPIM) :: NFPGRIB
INTEGER(KIND=JPIM) :: NFPWRITE = 1
INTEGER(KIND=JPIM) :: NFPXFLD = -999
INTEGER(KIND=JPIM) :: NFPDIGITS
INTEGER(KIND=JPIM) :: NFPDITS(0:JPNPST) = 0
INTEGER(KIND=JPIM) :: NFPDITSMIN(0:JPNPST) = 0
INTEGER(KIND=JPIM) :: NFRFPDI = 1
LOGICAL            :: LFPPGDFWR = .FALSE.
LOGICAL            :: LFPHISFWR = .FALSE.

END TYPE TNAMFPIOS

!     ------------------------------------------------------------------
END MODULE YOMFPIOS
