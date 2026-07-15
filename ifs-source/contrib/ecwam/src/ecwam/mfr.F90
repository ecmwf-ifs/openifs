! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE MFR(NFRE, IFRE1, FR1, FRATIO, FR)

! ----------------------------------------------------------------------

!**** *MFR* - ROUTINE TO COMPUTE FREQUENCY DIRECTION CONSTANTS.


!*    PURPOSE.
!     --------

!       INITIATES THE FREQUENCY ARRAY

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JWIM), INTENT(IN) :: NFRE  ! total number of frequencies
INTEGER(KIND=JWIM), INTENT(IN) :: IFRE1 ! index of the reference frequency
REAL(KIND=JWRB), INTENT(IN) :: FR1 ! reference frequency (Hz)
REAL(KIND=JWRB), INTENT(IN) :: FRATIO ! ratio between 2 concecutive frequncies
REAL(KIND=JWRB), DIMENSION(NFRE), INTENT(OUT) :: FR ! frequency array

INTEGER(KIND=JWIM) :: M

! ----------------------------------------------------------------------

FR(IFRE1) = FR1
DO M=IFRE1-1,1,-1
  FR(M) = FR(M+1)/FRATIO
ENDDO
DO M=IFRE1+1,NFRE
  FR(M) = FRATIO*FR(M-1)
ENDDO

END SUBROUTINE MFR
