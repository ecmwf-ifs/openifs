! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FSCAD_MOD
CONTAINS
SUBROUTINE FSCAD(KGL,KF_UV,KF_SCALARS,KF_SCDERS,&
 & PUV,PSCALAR,PNSDERS,PEWDERS,PUVDERS)

!**** *FSCAD - Division by a*cos(theta), east-west derivatives - adjoint

!     Purpose.
!     --------
!        In Fourier space divide u and v and all north-south
!        derivatives by a*cos(theta). Also compute east-west derivatives
!        of u,v,thermodynamic, passiv scalar variables and surface
!        pressure.

!**   Interface.
!     ----------
!        CALL FSCAD(..)
!        Explicit arguments :  PUV     - u and v
!        --------------------  PSCALAR - scalar valued varaibles
!                              PNSDERS - N-S derivative of S.V.V.
!                              PEWDERS - E-W derivative of S.V.V.
!                              PUVDERS - E-W derivative of u and v
!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03 (From SC2FSC)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_TRANS       ,ONLY : LUVDER
USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_FIELDS      ,ONLY : F
USE TPM_GEOMETRY    ,ONLY : G
!

IMPLICIT NONE

INTEGER(KIND=JPIM) , INTENT(IN) :: KGL,KF_UV,KF_SCALARS,KF_SCDERS
REAL(KIND=JPRB) , INTENT(INOUT) :: PUV(:,:)
REAL(KIND=JPRB) , INTENT(INOUT) :: PSCALAR(:,:)
REAL(KIND=JPRB) , INTENT(INOUT) :: PNSDERS(:,:)
REAL(KIND=JPRB) , INTENT(INOUT) :: PEWDERS(:,:)
REAL(KIND=JPRB) , INTENT(INOUT) :: PUVDERS(:,:)

REAL(KIND=JPRB) :: ZACHTE,ZMUL
INTEGER(KIND=JPIM) :: IMEN,ISTAGTF


INTEGER(KIND=JPIM) :: JLON,JF,IGLG,II,IR,JM

!     ------------------------------------------------------------------

IGLG = D%NPTRLS(MYSETW)+KGL-1
ZACHTE  = F%RACTHE(IGLG)
IMEN    = G%NMEN(IGLG)
ISTAGTF = D%NSTAGTF(KGL)


!     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF(LUVDER)THEN
  DO JM=0,IMEN
    IR = ISTAGTF+2*JM+1
    II = IR+1
    ZMUL = ZACHTE*JM
    DO JF=1,2*KF_UV
      PUV(JF,II) = PUV(JF,II) - PUVDERS(JF,IR)*ZMUL
      PUV(JF,IR) = PUV(JF,IR) + PUVDERS(JF,II)*ZMUL
!      PUVDERS(JF,IR) = _ZERO_
!      PUVDERS(JF,II) = _ZERO_
    ENDDO
  ENDDO
ENDIF

!*       2.2     SCALAR VARIABLES

IF(KF_SCDERS > 0)THEN
  DO JM=0,IMEN
    IR = ISTAGTF+2*JM+1
    II = IR+1
    ZMUL = ZACHTE*JM
    DO JF=1,KF_SCALARS
      PSCALAR(JF,II) = PSCALAR(JF,II) - PEWDERS(JF,IR)*ZMUL
      PSCALAR(JF,IR) = PSCALAR(JF,IR) + PEWDERS(JF,II)*ZMUL
!      PEWDERS(JF,IR) = _ZERO_
!      PEWDERS(JF,II) = _ZERO_
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!*       1.    DIVIDE U V AND N-S DERIVATIVES BY A*COS(THETA)
!              ----------------------------------------------

  
!*       1.1      U AND V.

IF(KF_UV > 0) THEN
  DO JLON=ISTAGTF+1,ISTAGTF+2*(IMEN+1)
    DO JF=1,2*KF_UV
      PUV(JF,JLON) = PUV(JF,JLON)*ZACHTE
    ENDDO
  ENDDO
ENDIF

!*      1.2      N-S DERIVATIVES

IF(KF_SCDERS > 0)THEN
  DO JLON=ISTAGTF+1,ISTAGTF+2*(IMEN+1)
    DO JF=1,KF_SCALARS
      PNSDERS(JF,JLON) = PNSDERS(JF,JLON)*ZACHTE
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FSCAD
END MODULE FSCAD_MOD
