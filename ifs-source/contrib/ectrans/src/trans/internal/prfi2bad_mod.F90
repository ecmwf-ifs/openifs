! (C) Copyright 1990- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PRFI2BAD_MOD
CONTAINS
SUBROUTINE PRFI2BAD(KFIELD,KM,KMLOC,PAIA,PSIA)

!**** *PRFI2BAD* - Prepare input work arrays for direct transform

!     Purpose.
!     --------
!        To extract the Fourier fields for a specific zonal wavenumber
!        and put them in an order suitable for the direct Legendre
!        tranforms, i.e. split into symmetric and anti-symmetric part.

!**   Interface.
!     ----------
!     *CALL* *PRFI2BAD(..)

!        Explicit arguments :
!        -------------------   KFIELD - number of fields
!                              KM - zonal wavenumber
!                              KMLOC - local zonal wavenumber
!                              PAOA - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              PSOA - symmetric part of Fourier
!                              fields for zonal wavenumber KM

!        Implicit arguments :  FOUBUF in TPM_TRANS
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
!        Original : 90-07-01
!        MPP Group: 95-10-01 Support for Distributed Memory version
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_TRANS       ,ONLY : FOUBUF
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_DISTR       ,ONLY : D
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD,KM,KMLOC
REAL(KIND=JPRB)  , INTENT(IN)  :: PSIA(:,:),   PAIA(:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IGLS,  ISL, ISTAN, ISTAS, JF, JGL


!     ------------------------------------------------------------------

!*       1.    EXTRACT SYM./ANTISYM. FIELDS FROM FOURIER ARRAY.
!              ------------------------------------------------

ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)

DO JGL=ISL,R%NDGNH
  IGLS = R%NDGL+1-JGL
  ISTAN = (D%NSTAGT1B(D%NPROCL(JGL) )+D%NPNTGTB1(KMLOC,JGL ))*2*KFIELD
  ISTAS = (D%NSTAGT1B(D%NPROCL(IGLS))+D%NPNTGTB1(KMLOC,IGLS))*2*KFIELD
!DIR$ IVDEP
!OCL    NOVREC
  DO JF=1,KFIELD*2
    FOUBUF(ISTAN+JF) = PSIA(JF,JGL)+PAIA(JF,JGL)
    FOUBUF(ISTAS+JF) = PSIA(JF,JGL)-PAIA(JF,JGL)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE PRFI2BAD
END MODULE PRFI2BAD_MOD
