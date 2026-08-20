! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTEGER(KIND=JWIM) FUNCTION JAFU (CL, J, IAN)

! ----------------------------------------------------------------------

!**** *JAFU* - FUNCTION TO COMPUTE THE INDEX ARRAY FOR THE
!              ANGLES OF THE INTERACTING WAVENUMBERS.

!     S. HASSELMANN        MPIFM        01/12/1985.

!*    PURPOSE.
!     --------

!       INDICES DEFINING BINS IN FREQUENCY AND DIRECTION PLANE INTO
!       WHICH NONLINEAR ENERGY TRANSFER INCREMENTS ARE STORED. NEEDED
!       FOR COMPUTATION OF THE NONLINEAR ENERGY TRANSFER.

!**   INTERFACE.
!     ----------

!       *FUNCTION* *JAFU (CL, J, IAN)*
!          *CL*  - WEIGHTS.
!          *J*   - INDEX IN ANGULAR ARRAY.
!          *IAN* - NUMBER OF ANGLES IN ARRAY.

!     METHOD.
!     -------

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!        S. HASSELMANN AND K. HASSELMANN,JPO, 1985 B.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE
      REAL(KIND=JWRB), INTENT(IN):: CL
      INTEGER(KIND=JWIM), INTENT(IN) :: J, IAN

      INTEGER(KIND=JWIM) :: IDPH, JA

      IDPH = CL
      JA = J+IDPH
      IF (JA <= 0)   JA = IAN+JA-1
      IF (JA >= IAN) JA = JA-IAN+1
      JAFU = JA

END FUNCTION JAFU
