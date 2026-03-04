! (C) Copyright 2001- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE ASRE1AD_MOD
CONTAINS
SUBROUTINE ASRE1AD(KM,KMLOC,KF_OUT_LT,PAOA1,PSOA1)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!USE TPM_TRANS

USE ASRE1BAD_MOD    ,ONLY : ASRE1BAD


!**** *ASRE1AD* - Recombine antisymmetric and symmetric parts - adjoint

!     Purpose.
!     --------
!        To recombine the antisymmetric and symmetric parts of the
!        Fourier arrays and update the correct parts of the state
!        variables.

!**   Interface.
!     ----------
!       *CALL* *ASRE1AD(...)

!        Explicit arguments :
!        --------------------
!                              KM - zonal wavenumber
!                              KMLOC - local zonal wavenumber
!                              PAOA1 - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM (basic
!                              variables and N-S derivatives)
!                              PSOA1 - symmetric part of Fourier
!                              fields for zonal wavenumber KM (basic
!                              variables and N-S derivatives)

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------

!     Externals.   ASRE1BAD - basic recombination routine
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From ASRE1AD in IFS CY22R1

!     ------------------------------------------------------------------


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM) , INTENT(IN)  :: KM
INTEGER(KIND=JPIM) , INTENT(IN)  :: KMLOC
INTEGER(KIND=JPIM) , INTENT(IN)  :: KF_OUT_LT

REAL(KIND=JPRB)    , INTENT(OUT) :: PSOA1(:,:),       PAOA1(:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IFLDS


!     ------------------------------------------------------------------

IFLDS = KF_OUT_LT

CALL ASRE1BAD(IFLDS,KM,KMLOC,PAOA1,PSOA1)

!     ------------------------------------------------------------------

END SUBROUTINE ASRE1AD
END MODULE ASRE1AD_MOD
