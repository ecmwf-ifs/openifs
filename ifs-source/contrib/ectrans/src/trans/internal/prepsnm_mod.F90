! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PREPSNM_MOD
CONTAINS
SUBROUTINE PREPSNM(KM,KMLOC,PEPSNM)


!**** *PREPSNM* - Prepare REPSNM for wavenumber KM

!     Purpose.
!     --------
!        Copy the REPSNM values for specific zonal wavenumber M
!        to work array

!**   Interface.
!     ----------
!        CALL PREPSNM(...)

!        Explicit arguments :  KM - zonal wavenumber
!        -------------------   KMLOC - local zonal wavenumber
!                              PEPSNM - REPSNM for zonal
!                                      wavenumber KM

!        Implicit arguments :
!        --------------------

!     Method.
!     -------


!     Reference.
!     ----------


!     Author.
!     -------
!        Lars Isaksen *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From LTINV in IFS CY22R1

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_FIELDS      ,ONLY : F
USE TPM_DISTR       ,ONLY : D
!

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KM,KMLOC
REAL(KIND=JPRB),    INTENT(OUT) :: PEPSNM(0:R%NTMAX+2)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: JN

!     ------------------------------------------------------------------

!*       1.       COPY REPSNM.
!                 ------------


IF (KM > 0) THEN
  PEPSNM(0:KM-1) = 0.0_JPRB
ENDIF

DO JN=KM,R%NTMAX+2
  PEPSNM(JN) = F%REPSNM(D%NPMT(KM)+KMLOC-KM+JN)
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE PREPSNM
END MODULE PREPSNM_MOD

