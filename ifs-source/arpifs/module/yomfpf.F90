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

MODULE YOMFPF

USE PARKIND1, ONLY : JPIM, JPRB
USE PARFPOS,  ONLY : JPOSDOM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Full-POS spectral filter

! NFMAX    : maximum truncation of the output subdomains.
!            If negative : means no filter at all.

! LFPBED   : .TRUE. to use the low-pass filter ("THX" filter) on the homogenous resolution
!            space (arpege only); otherwise a gaussian filter is used.
! RFPBED   : coefficient of the exponential function in the "THX" filter.
! RFPLTF   : coefficient of the exponential function in the gaussian filter
! NOTICE : that parameters applies only if all fields should be smoothed (%ISF=3)

! RFPSEL   : controls the selectivity of the THX filter; always > 0.
!            increasing (resp. decreasing) RFPSEL: less (resp. more) selective THX filter.
! RFPMAXDEV: maximum deviation against identity of dilatation/contraction matrixes composition
! CFPDILA  : filename of dilatation matrix
! CFPCONT  : filename of contraction matrix
! CFPMATRD : filename of filtering matrixes to be read
! CFPMATWR : filename of filtering matrixes to be computed

! LFPWRFIL : Write out filtering matrixes on files
! LFPRDFIL : Read filtering matrixes from files
! NFPREADALL : 1 = all MPI tasks read filtering matrixes (1)
!              0 = matrixes reading is distributed among the V-set,
!                  then V-sets communicate to each other
! NFPGLFI  : granularity factor for LFI filtering matrixes
! NFPCMAX  : maximum truncation of the stretched model (> NSMAX*RSTRET)

TYPE TNAMFPF

LOGICAL :: LFPBED
LOGICAL :: LFPWRFIL
LOGICAL :: LFPRDFIL

INTEGER(KIND=JPIM) :: NFPREADALL

INTEGER(KIND=JPIM) :: NFPGLFI
INTEGER(KIND=JPIM) :: NFPCMAX
INTEGER(KIND=JPIM) :: NFMAX(JPOSDOM)

REAL(KIND=JPRB) :: RFPMAXDEV
REAL(KIND=JPRB) :: RFPLTF
REAL(KIND=JPRB) :: RFPBED
REAL(KIND=JPRB) :: RFPSEL(JPOSDOM)

CHARACTER(LEN=256) :: CFPDILA
CHARACTER(LEN=256) :: CFPCONT
CHARACTER(LEN=256) :: CFPMATRD(JPOSDOM)
CHARACTER(LEN=256) :: CFPMATWR(JPOSDOM)

END TYPE TNAMFPF

!     ------------------------------------------------------------------
END MODULE YOMFPF
