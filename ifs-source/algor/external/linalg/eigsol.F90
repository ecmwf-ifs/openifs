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

SUBROUTINE EIGSOL(KFLEVG,KNFLEVG,PA,PFR,PFI,K,PMO,KWO,PWO,KER)

!**** EIGSOL - Routine to avoid calling a routine named RG (compiler problem)
!              Calculation of eigenvalues and eigenvectors of a matrix.

!     Purpose.
!     --------
!      Computes eigenvalues and eigenvectors of matrix PA.

!**   Interface.
!     ----------
!        *CALL* *EIGSOL(...)

!        Explicit arguments :
!        --------------------
!         KFLEVG    - row dimension of the 2D array PA                    (in)
!         KNFLEVG   - order of the matrix PA                              (in)
!         PA        - input matrix, the eigenvalues of which are computed (in)
!         PFR       - real part of eigenvalues                            (out)
!         PFI       - imaginary part of eigenvalues                       (out)
!         K         - 0: only eigenvalues are desired                     (in)
!                     otherwise: eigenvalues and eigenvectors are desired
!         PMO       - eigenvectors                                        (out)
!         KWO       - temporary storage array                             (out)
!         PWO       - temporary storage array                             (out)
!         KER       - error completion code                               (out)

!        Implicit arguments :
!        --------------------
!        none.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      ???

!     Modifications.
!     --------------
!      K. Yessad (Apr 2009): add missing comments, restore INTENT.
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)   :: KFLEVG
INTEGER(KIND=JPIM),INTENT(IN)   :: KNFLEVG
REAL(KIND=JPRB),INTENT(IN)      :: PA(*)
REAL(KIND=JPRB),INTENT(OUT)     :: PFR(*)
REAL(KIND=JPRB),INTENT(OUT)     :: PFI(*)
INTEGER(KIND=JPIM),INTENT(IN)   :: K
REAL(KIND=JPRB),INTENT(OUT)     :: PMO(*)
INTEGER(KIND=JPIM),INTENT(OUT)  :: KWO(*)
REAL(KIND=JPRB),INTENT(OUT)     :: PWO(*)
INTEGER(KIND=JPIM),INTENT(OUT)  :: KER

!     ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('EIGSOL',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

CALL RG(KFLEVG,KNFLEVG,PA,PFR,PFI,K,PMO,KWO,PWO,KER)

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('EIGSOL',1,ZHOOK_HANDLE)
END SUBROUTINE EIGSOL
