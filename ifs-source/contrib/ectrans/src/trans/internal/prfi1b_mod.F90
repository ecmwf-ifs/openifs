! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PRFI1B_MOD
CONTAINS
SUBROUTINE PRFI1B(KM,PIA,PSPEC,KFIELDS,KFLDPTR)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D


!**** *PRFI1* - Prepare spectral fields for inverse Legendre transform

!     Purpose.
!     --------
!        To extract the spectral fields for a specific zonal wavenumber
!        and put them in an order suitable for the inverse Legendre           .
!        tranforms.The ordering is from NSMAX to KM for better conditioning.
!        Elements 1,2 and NLCM(KM)+1 are zeroed in preparation for computing
!        u,v and derivatives in spectral space.

!**   Interface.
!     ----------
!        *CALL* *PRFI1B(...)*

!        Explicit arguments :  KM     - zonal wavenumber
!        ------------------    PIA    - spectral components for transform
!                              PSPEC  - spectral array
!                              KFIELDS  - number of fields


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
!        Original : 00-02-01 From PRFI1B in IFS CY22R1

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)   :: KM,KFIELDS
REAL(KIND=JPRB)   ,INTENT(IN)   :: PSPEC(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)  :: PIA(:,:)
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: II, INM, IR, J, JFLD, ILCM, IOFF,IFLD


!     ------------------------------------------------------------------

!*       1.    EXTRACT FIELDS FROM SPECTRAL ARRAYS.
!              --------------------------------------------------


ILCM = R%NSMAX+1-KM
IOFF = D%NASM0(KM)

IF(PRESENT(KFLDPTR)) THEN
  DO JFLD=1,KFIELDS
    IR = 2*(JFLD-1)+1
    II = IR+1
    IFLD = KFLDPTR(JFLD)
    DO J=1,ILCM
      INM = IOFF+(ILCM-J)*2
      PIA(J+2,IR) = PSPEC(IFLD,INM  )
      PIA(J+2,II) = PSPEC(IFLD,INM+1)
    ENDDO
  ENDDO

ELSE
  DO J=1,ILCM
    INM = IOFF+(ILCM-J)*2
    !DIR$ IVDEP
    !OCL NOVREC
    DO JFLD=1,KFIELDS
      IR = 2*(JFLD-1)+1
      II = IR+1
      PIA(J+2,IR) = PSPEC(JFLD,INM  )
      PIA(J+2,II) = PSPEC(JFLD,INM+1)
    ENDDO
  ENDDO

ENDIF

DO JFLD=1,2*KFIELDS
  PIA(1,JFLD) = 0.0_JPRB
  PIA(2,JFLD) = 0.0_JPRB
  PIA(ILCM+3,JFLD) = 0.0_JPRB
ENDDO


!     ------------------------------------------------------------------

END SUBROUTINE PRFI1B
END MODULE PRFI1B_MOD
