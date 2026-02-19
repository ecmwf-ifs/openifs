! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE ASRE1BAD_MOD
CONTAINS
SUBROUTINE ASRE1BAD(KFIELD,KM,KMLOC,PAOA,PSOA)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_TRANS       ,ONLY : FOUBUF_IN
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_DISTR       ,ONLY : D


!**** *ASRE1BAD* - Recombine antisymmetric and symmetric parts - adjoint

!     Purpose.
!     --------
!        To recombine the antisymmetric and symmetric parts of the
!        Fourier arrays and update the correct parts of the state
!        variables.

!**   Interface.
!     ----------
!        *CALL* *ASRE1BAD(..)

!        Explicit arguments :
!        -------------------   KFIELD - number of fields (input-c)
!                              KM - zonal wavenumber(input-c)
!                              KMLOC - local version of KM (input-c)
!                              PAOA - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM (input)
!                              PSOA - symmetric part of Fourier
!                              fields for zonal wavenumber KM (input)

!        Implicit arguments :  FOUBUF_IN - output buffer (output)
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
!        Original : 00-02-01 From ASRE1BAD in IFS CY22R1

!     ------------------------------------------------------------------


IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KFIELD,KM,KMLOC
REAL(KIND=JPRB), INTENT(OUT)    :: PSOA(:,:)
REAL(KIND=JPRB), INTENT(OUT)    :: PAOA(:,:)

!     LOCAL INTEGERS
INTEGER(KIND=JPIM) :: ISL, IGLS, JFLD, JGL ,IPROC, IPROCS, IDGNH
INTEGER(KIND=JPIM) :: ISTAN(R%NDGNH),ISTAS(R%NDGNH)

!     ------------------------------------------------------------------

!*       1.    RECOMBINATION  OF SYMMETRIC AND ANTSYMMETRIC PARTS.
!              ---------------------------------------------------

ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)
IDGNH = R%NDGNH

!*       1.2      RECOMBINE

DO JGL=ISL,IDGNH
  IPROC = D%NPROCL(JGL)
  ISTAN(JGL) = (D%NSTAGT0B(IPROC) + D%NPNTGTB1(KMLOC,JGL))*2*KFIELD
  IGLS = R%NDGL+1-JGL
  IPROCS = D%NPROCL(IGLS)
  ISTAS(JGL) = (D%NSTAGT0B(IPROCS) + D%NPNTGTB1(KMLOC,IGLS))*2*KFIELD
ENDDO

DO JGL=ISL,IDGNH
!OCL      NOVREC
  DO JFLD=1,2*KFIELD
    PSOA(JFLD,JGL) = FOUBUF_IN(ISTAN(JGL)+JFLD)+FOUBUF_IN(ISTAS(JGL)+JFLD)
    PAOA(JFLD,JGL) = FOUBUF_IN(ISTAN(JGL)+JFLD)-FOUBUF_IN(ISTAS(JGL)+JFLD)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE ASRE1BAD
END MODULE ASRE1BAD_MOD

