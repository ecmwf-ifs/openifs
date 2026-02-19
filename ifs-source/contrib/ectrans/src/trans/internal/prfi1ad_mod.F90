! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PRFI1AD_MOD
CONTAINS
SUBROUTINE PRFI1AD(KM,KF_UV,KF_SCALARS,PIA,PSPVOR,PSPDIV,PSPSCALAR,&
 & KFLDPTRUV,KFLDPTRSC)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!USE TPM_DISTR
!USE TPM_TRANS

USE PRFI1BAD_MOD   ,ONLY : PRFI1BAD


!**** *PRFI1AD* - Prepare spectral fields for inverse Legendre transform

!     Purpose.
!     --------
!        To extract the spectral fields for a specific zonal wavenumber
!        and put them in an order suitable for the inverse Legendre           .
!        tranforms.The ordering is from NSMAX to KM for better conditioning.
!        Elements 1,2 and NLCM(KM)+1 are zeroed in preparation for computing
!        u,v and derivatives in spectral space.

!**   Interface.
!     ----------
!        *CALL* *PRFI1AD(KM,PIA,PSPVOR,PSPDIV,PSPSCALAR

!        Explicit arguments :  KM     - zonal wavenumber
!        ------------------    PIA    - spectral components for transform
!                              PSPVOR - vorticity
!                              PSPDIV - divergence
!                              PSPSCALAR - scalar variables

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From PRFI1AD in IFS CY22R1

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KM,KF_UV,KF_SCALARS
REAL(KIND=JPRB) ,OPTIONAL, INTENT(INOUT)  :: PSPVOR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(INOUT)  :: PSPDIV(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(INOUT)  :: PSPSCALAR(:,:)
REAL(KIND=JPRB) ,          INTENT(IN)     :: PIA(:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KFLDPTRSC(:)


!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IDIV, IFIRST, ILAST, IVOR


!     ------------------------------------------------------------------

!*       1.    EXTRACT FIELDS FROM SPECTRAL ARRAYS.
!              ------------------------------------

IFIRST = 1
ILAST  = 4*KF_UV

!*       1.1      VORTICITY AND DIVERGENCE.

IF(KF_UV > 0)THEN
  IVOR = 1
  IDIV = 2*KF_UV+1
  CALL PRFI1BAD(KM,PIA(:,IVOR:IDIV-1),PSPVOR,KF_UV,KFLDPTRUV)
  CALL PRFI1BAD(KM,PIA(:,IDIV:ILAST) ,PSPDIV,KF_UV,KFLDPTRUV)
  ILAST = ILAST+4*KF_UV
ENDIF

!*       1.2    SCALAR VARIABLES.

IF(KF_SCALARS > 0)THEN
  IFIRST = ILAST+1
  ILAST  = IFIRST - 1 + 2*KF_SCALARS
  CALL PRFI1BAD(KM,PIA(:,IFIRST:ILAST),PSPSCALAR(:,:),KF_SCALARS,KFLDPTRSC)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE PRFI1AD
END MODULE PRFI1AD_MOD


