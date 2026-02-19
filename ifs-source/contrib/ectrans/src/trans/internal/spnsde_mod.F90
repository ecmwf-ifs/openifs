! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SPNSDE_MOD
CONTAINS
SUBROUTINE SPNSDE(KM,KF_SCALARS,PEPSNM,PF,PNSD)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!USE TPM_GEN
USE TPM_DIM         ,ONLY : R
USE TPM_FIELDS      ,ONLY : F
!USE TPM_TRANS


!**** *SPNSDE* - Compute North-South derivative in spectral space

!     Purpose.
!     --------
!        In Laplace space compute the the North-south derivative

!**   Interface.
!     ----------
!        CALL SPNSDE(...)

!        Explicit arguments :
!        --------------------
!        KM -zonal wavenumber (input-c)
!        PEPSNM - REPSNM for wavenumber KM (input-c)
!        PF  (NLEI1,2*KF_SCALARS) - input field (input)
!        PNSD(NLEI1,2*KF_SCALARS) - N-S derivative (output)

!        Organisation within NLEI1:
!        NLEI1 = NSMAX+4+mod(NSMAX+4+1,2)
!                        overdimensioning
!        1        : n=NSMAX+2
!        2        : n=NSMAX+1
!        3        : n=NSMAX
!        .        :
!        .        :
!        NSMAX+3  : n=0
!        NSMAX+4  : n=-1

!        Implicit arguments :  YOMLAP
!        --------------------

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Temperton, 1991, MWR 119 p1303

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From SPNSDE in IFS CY22R1

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KM
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_SCALARS
REAL(KIND=JPRB),    INTENT(IN)  :: PEPSNM(0:R%NTMAX+2)
REAL(KIND=JPRB),    INTENT(IN)  :: PF(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PNSD(:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IJ, ISKIP, J, JN,JI,ISMAX
REAL(KIND=JPRB) :: ZEPSNM(-1:R%NSMAX+4)
REAL(KIND=JPRB) :: ZN(-1:R%NTMAX+4)


!     ------------------------------------------------------------------

!*       1.    COMPUTE NORTH SOUTH DERIVATIVE.
!              -------------------------------


!*       1.1      COMPUTE

ISMAX = R%NSMAX
DO JN=KM-1,ISMAX+2
  IJ = ISMAX+3-JN
  ZN(IJ) = F%RN(JN)
  IF( JN >= 0 ) ZEPSNM(IJ) = PEPSNM(JN)
ENDDO
ZN(0) = F%RN(ISMAX+3)

IF(KM == 0) THEN
  ISKIP = 2
ELSE
  ISKIP = 1
ENDIF

DO J=1,2*KF_SCALARS,ISKIP
  DO JI=2,ISMAX+3-KM
    PNSD(JI,J) = -ZN(JI+1)*ZEPSNM(JI)*PF(JI+1,J)+&
     &ZN(JI-2)*ZEPSNM(JI-1)*PF(JI-1,J)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE SPNSDE
END MODULE SPNSDE_MOD
