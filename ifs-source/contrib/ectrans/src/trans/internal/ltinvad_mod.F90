! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTINVAD_MOD
CONTAINS
SUBROUTINE LTINVAD(KM,KMLOC,KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,KLEI2,KDIM1,&
 & PSPVOR,PSPDIV,PSPSCALAR,&
 & PSPSC3A,PSPSC3B,PSPSC2 , &
 & KFLDPTRUV,KFLDPTRSC,FSPGL_PROC)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_TRANS       ,ONLY : LDIVGP, LVORGP, NF_SC2, NF_SC3A, NF_SC3B
USE TPM_GEOMETRY

!USE PRLE1AD_MOD
USE PREPSNM_MOD     ,ONLY : PREPSNM
USE PRFI1BAD_MOD    ,ONLY : PRFI1BAD
USE VDTUVAD_MOD     ,ONLY : VDTUVAD
USE SPNSDEAD_MOD    ,ONLY : SPNSDEAD
USE LEINVAD_MOD     ,ONLY : LEINVAD
USE ASRE1BAD_MOD    ,ONLY : ASRE1BAD
!USE FSPGL_INT_MOD


!**** *LTINVAD* - Inverse Legendre transform

!     Purpose.
!     --------
!        Tranform from Laplace space to Fourier space, compute U and V
!        and north/south derivatives of state variables.

!**   Interface.
!     ----------
!        *CALL* *LTINVAD(...)

!        Explicit arguments :
!        --------------------
!          KM        - zonal wavenumber
!          KMLOC     - local zonal wavenumber
!          PSPVOR    - spectral vorticity
!          PSPDIV    - spectral divergence
!          PSPSCALAR - spectral scalar variables

!        Implicit arguments :  The Laplace arrays of the model.
!        --------------------  The values of the Legendre polynomials
!                              The grid point arrays of the model
!     Method.
!     -------

!     Externals.
!     ----------
!         PRLE1AD - prepares the Legendre polonymials
!         PREPSNM - prepare REPSNM for wavenumber KM
!         PRFI1AD - prepares the spectral fields
!         VDTUVAD - compute u and v from vorticity and divergence
!         SPNSDEAD- compute north-south derivatives
!         LEINVAD - Inverse Legendre transform
!         ASRE1AD - recombination of symmetric/antisymmetric part

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Temperton, 1991, MWR 119 p1303

!     Author.
!     -------
!        Mats Hamrud  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From LTINVAD in IFS CY22R1
!     ------------------------------------------------------------------

IMPLICIT NONE


INTEGER(KIND=JPIM), INTENT(IN) :: KM
INTEGER(KIND=JPIM), INTENT(IN) :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN) :: KF_OUT_LT
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCDERS
INTEGER(KIND=JPIM), INTENT(IN) :: KLEI2
INTEGER(KIND=JPIM), INTENT(IN) :: KDIM1

REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(INOUT) :: PSPSC2(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(INOUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(INOUT) :: PSPSC3B(:,:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KFLDPTRSC(:)
EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC

REAL(KIND=JPRB) :: ZIA(R%NLEI1,KLEI2)
REAL(KIND=JPRB) :: ZEPSNM(0:R%NTMAX+2)
REAL(KIND=JPRB) :: ZSOA1(KDIM1,R%NLEI3),ZAOA1(KDIM1,R%NLEI3)


!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IFC, ISTA, IIFC, IDGLU
INTEGER(KIND=JPIM) :: IVORL,IVORU,IDIVL,IDIVU,IUL,IUU,IVL,IVU,ISL,ISU,IDL,IDU
INTEGER(KIND=JPIM) :: ILAST,IFIRST,IDIM1,IDIM3,J3

!     LOCAL LOGICAL SCALARS

!     LOCAL REAL SCALARS

!     ------------------------------------------------------------------

!*       1.    PREPARE AND ZEPSNM.
!              -------------------

CALL PREPSNM(KM,KMLOC,ZEPSNM)

!     ------------------------------------------------------------------
!     6. OPTIONAL COMPUTATIONS IN FOURIER SPACE

!IF(PRESENT(FSPGL_PROC)) THEN
!  CALL FSPGL_INT(KM,KMLOC,FSPGL_PROC)
!ENDIF

!     ------------------------------------------------------------------

!*       5.    RECOMBINATION SYMMETRIC/ANTISYMMETRIC PART.
!              --------------------------------------------

CALL ASRE1BAD(KF_OUT_LT,KM,KMLOC,ZAOA1,ZSOA1)

!     ------------------------------------------------------------------

!*       4.    INVERSE LEGENDRE TRANSFORM.
!              ---------------------------


ISTA = 1
IFC  = 2*KF_OUT_LT
IF(KF_UV > 0 .AND. .NOT. LVORGP) THEN
  ISTA = ISTA+2*KF_UV
ENDIF
IF(KF_UV > 0 .AND. .NOT. LDIVGP) THEN
  ISTA = ISTA+2*KF_UV
ENDIF

ZIA(:,ISTA:ISTA+IFC-1) = 0.0_JPRB

IDGLU = MIN(R%NDGNH,G%NDGLU(KM))
IIFC=IFC
IF(KM == 0)THEN
  IIFC=IFC/2
ENDIF
CALL LEINVAD(KM,KMLOC,IFC,IIFC,KF_OUT_LT,IDGLU,ZIA(:,ISTA:ISTA+IFC-1),ZAOA1,ZSOA1)

!     ------------------------------------------------------------------

!*       3.    SPECTRAL COMPUTATIONS FOR U,V AND DERIVATIVES.
!              ----------------------------------------------

ZIA(:,1:ISTA-1) = 0.0_JPRB

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
  CALL VDTUVAD(KM,KF_UV,ZEPSNM,ZIA(:,IVORL:IVORU),ZIA(:,IDIVL:IDIVU),&
            & ZIA(:,IUL:IUU),ZIA(:,IVL:IVU))
  CALL PRFI1BAD(KM,ZIA(:,IVORL:IVORU),PSPVOR,KF_UV,KFLDPTRUV)
  CALL PRFI1BAD(KM,ZIA(:,IDIVL:IDIVU),PSPDIV,KF_UV,KFLDPTRUV)
  ILAST = ILAST+4*KF_UV
ENDIF

IF (KF_SCDERS > 0) THEN
  ISL = 2*(4*KF_UV)+1
  ISU = ISL+2*KF_SCALARS-1
  IDL = 2*(4*KF_UV+KF_SCALARS)+1
  IDU = IDL+2*KF_SCDERS-1
  CALL SPNSDEAD(KM,KF_SCALARS,ZEPSNM,ZIA(:,ISL:ISU),ZIA(:,IDL:IDU))
ENDIF

IF(KF_SCALARS > 0)THEN
  IF(PRESENT(PSPSCALAR)) THEN
    IFIRST = ILAST+1
    ILAST  = IFIRST - 1 + 2*KF_SCALARS
    CALL PRFI1BAD(KM,ZIA(:,IFIRST:ILAST),PSPSCALAR(:,:),KF_SCALARS,KFLDPTRSC)
  ELSE
    IF(PRESENT(PSPSC2) .AND. NF_SC2 > 0) THEN
      IFIRST = ILAST+1
      ILAST  = IFIRST-1+2*NF_SC2
      CALL PRFI1BAD(KM,ZIA(:,IFIRST:ILAST),PSPSC2(:,:),NF_SC2)
    ENDIF
    IF(PRESENT(PSPSC3A) .AND. NF_SC3A > 0) THEN
      IDIM1=NF_SC3A
      IDIM3=UBOUND(PSPSC3A,3)
      DO J3=1,IDIM3
        IFIRST = ILAST+1
        ILAST  = IFIRST-1+2*IDIM1
        CALL PRFI1BAD(KM,ZIA(:,IFIRST:ILAST),PSPSC3A(:,:,J3),IDIM1)
      ENDDO
    ENDIF
    IF(PRESENT(PSPSC3B) .AND. NF_SC3B > 0) THEN
      IDIM1=NF_SC3B
      IDIM3=UBOUND(PSPSC3B,3)
      DO J3=1,IDIM3
        IFIRST = ILAST+1
        ILAST  = IFIRST-1+2*IDIM1
        CALL PRFI1BAD(KM,ZIA(:,IFIRST:ILAST),PSPSC3B(:,:,J3),IDIM1)
      ENDDO
    ENDIF
  ENDIF
ENDIF


!     ------------------------------------------------------------------


END SUBROUTINE LTINVAD
END MODULE LTINVAD_MOD




