! (C) Copyright 2015- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE VD2UV_MOD
CONTAINS
SUBROUTINE VD2UV(KM,KMLOC,KF_UV,KLEI2,PSPVOR,PSPDIV,PU,PV)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_CONSTANTS
USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D

USE PREPSNM_MOD     ,ONLY : PREPSNM
USE PRFI1B_MOD      ,ONLY : PRFI1B
USE VDTUV_MOD       ,ONLY : VDTUV


!**** *VD2UV* - U and V from Vor/div
!
!     Purpose.
!     --------
!       
!**   Interface.
!     ----------
!        *CALL* *VD2UV(...)

!        Explicit arguments :
!        --------------------
!          KM        - zonal wavenumber
!          KMLOC     - local zonal wavenumber
!          PSPVOR    - spectral vorticity
!          PSPDIV    - spectral divergence
!          PU(:,:)   - spectral U (out)
!          PV(:,:)   - spectral V (out)


!        Implicit arguments :  

!     Method.
!     -------

!     Externals.
!     ----------

!         PREPSNM - prepare REPSNM for wavenumber KM
!         PRFI1B  - prepares the spectral fields
!         VDTUV   - compute u and v from vorticity and divergence

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Temperton, 1991, MWR 119 p1303

!     Author.
!     -------
!        Mats Hamrud  *ECMWF*

!     Modifications.
!     --------------
!        Original : July 2015
!       
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KM
INTEGER(KIND=JPIM), INTENT(IN) :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM), INTENT(IN) :: KLEI2

REAL(KIND=JPRB)   , INTENT(IN)  :: PSPVOR(:,:)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSPDIV(:,:)
REAL(KIND=JPRB)   , INTENT(OUT) :: PU(:,:)
REAL(KIND=JPRB)   , INTENT(OUT) :: PV(:,:)

REAL(KIND=JPRB) :: ZIA(R%NLEI1,KLEI2)
REAL(KIND=JPRB) :: ZEPSNM(0:R%NTMAX+2),ZA_R

INTEGER(KIND=JPIM) :: IFC, ISTA, IIFC, IDGLU, JGL, JFLD,ILCM
INTEGER(KIND=JPIM) :: IVORL,IVORU,IDIVL,IDIVU,IUL,IUU,IVL,IVU,ISL,II,IR,INM,J
INTEGER(KIND=JPIM) :: IFIRST, ILAST, IOFF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

IF (LHOOK) CALL DR_HOOK('VD2UV_MOD',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    PREPARE ZEPSNM.
!              ---------------

CALL PREPSNM(KM,KMLOC,ZEPSNM)

!     ------------------------------------------------------------------


!*       3.    SPECTRAL COMPUTATIONS FOR U,V AND DERIVATIVES.
!              ----------------------------------------------

IFIRST = 1
ILAST  = 4*KF_UV

IF (KF_UV > 0) THEN
  IVORL = 1
  IVORU = 2*KF_UV
  IDIVL = 2*KF_UV+1
  IDIVU = 4*KF_UV
  IUL   = 4*KF_UV+1
  IUU   = 6*KF_UV
  IVL   = 6*KF_UV+1
  IVU   = 8*KF_UV
  CALL PRFI1B(KM,ZIA(:,IVORL:IVORU),PSPVOR,KF_UV)
  CALL PRFI1B(KM,ZIA(:,IDIVL:IDIVU),PSPDIV,KF_UV)

  CALL VDTUV(KM,KF_UV,ZEPSNM,ZIA(:,IVORL:IVORU),ZIA(:,IDIVL:IDIVU),&
           & ZIA(:,IUL:IUU),ZIA(:,IVL:IVU))
  ILCM = R%NSMAX+1-KM
  IOFF = D%NASM0(KM)
  ZA_R = 1.0_JPRB/RA
  DO J=1,ILCM
    INM = IOFF+(ILCM-J)*2
    DO JFLD=1,KF_UV
      IR = 2*(JFLD-1)+1
      II = IR+1
      PU(JFLD,INM  ) = ZIA(J+2,IR+IUL-1)*ZA_R
      PU(JFLD,INM+1) = ZIA(J+2,II+IUL-1)*ZA_R
      PV(JFLD,INM  ) = ZIA(J+2,IR+IVL-1)*ZA_R
      PV(JFLD,INM+1) = ZIA(J+2,II+IVL-1)*ZA_R
    ENDDO
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('VD2UV_MOD',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE VD2UV
END MODULE VD2UV_MOD




