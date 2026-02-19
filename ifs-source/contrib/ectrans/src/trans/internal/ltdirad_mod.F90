! (C) Copyright 1987- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTDIRAD_MOD
CONTAINS
SUBROUTINE LTDIRAD(KM,KMLOC,KF_FS,KF_UV,KF_SCALARS,KLED2,&
 & PSPVOR,PSPDIV,PSPSCALAR,&
 & PSPSC3A,PSPSC3B,PSPSC2 , &
 & KFLDPTRUV,KFLDPTRSC)


USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_GEOMETRY

USE PREPSNM_MOD     ,ONLY : PREPSNM
USE PRFI2AD_MOD     ,ONLY : PRFI2AD
USE LDFOU2AD_MOD    ,ONLY : LDFOU2AD
USE LEDIRAD_MOD     ,ONLY : LEDIRAD
USE UVTVDAD_MOD
USE UPDSPAD_MOD     ,ONLY : UPDSPAD
 

!**** *LTDIRAD* - Control of Direct Legendre transform step - adjoint

!     Purpose.
!     --------
!        Tranform from Fourier space to spectral space, compute
!        vorticity and divergence.

!**   Interface.
!     ----------
!        *CALL* *LTDIRAD(...)*

!        Explicit arguments :
!        --------------------  KM     - zonal wavenumber
!                              KMLOC  - local zonal wavenumber
!                              PSPVOR - spectral vorticity
!                              PSPDIV - spectral divergence
!                              PSPSCALAR - spectral scalar variables

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!         PREPSNM  - prepare REPSNM for wavenumber KM
!         PRFI2AD  - prepares the Fourier work arrays for model variables.
!         LDFOU2AD - computations in Fourier space
!         LEDIRAD  - direct Legendre transform
!         UVTVDAD  -
!         UPDSPAD  - updating of spectral arrays (fields)

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-11-24
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified 93-03-19 D. Giard - CDCONF='T' for tendencies
!        Modified 93-11-18 M. Hamrud - use only one Fourier buffer
!        Modified 94-04-06 R. El khatib Full-POS implementation
!        M.Hamrud  : 94-11-01 New conf 'G' - vor,div->vor,div
!                             instead of u,v->vor,div
!        MPP Group : 95-10-01 Support for Distributed Memory version
!        K. YESSAD (AUGUST 1996):
!               - Legendre transforms for transmission coefficients.
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!        R. El Khatib 12-Jul-2012 LDSPC2AD replaced by UVTVDAD
!     ------------------------------------------------------------------

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM), INTENT(IN)  :: KM
INTEGER(KIND=JPIM), INTENT(IN)  :: KMLOC
INTEGER(KIND=JPIM),INTENT(IN)   :: KF_FS,KF_UV,KF_SCALARS,KLED2

REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(INOUT) :: PSPSC2(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(INOUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(INOUT) :: PSPSC3B(:,:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KFLDPTRSC(:)


!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IFC, IIFC, IDGLU
INTEGER(KIND=JPIM) :: IUS, IUE, IVS, IVE, IVORS, IVORE, IDIVS, IDIVE

!     LOCAL REALS
REAL(KIND=JPRB) :: ZSIA(KLED2,R%NDGNH),       ZAIA(KLED2,R%NDGNH)
REAL(KIND=JPRB) :: ZEPSNM(0:R%NTMAX+2)
REAL(KIND=JPRB) :: ZOA1(R%NLED4,KLED2),         ZOA2(R%NLED4,MAX(4*KF_UV,1))


!     ------------------------------------------------------------------

!*       1.    PREPARE LEGENDRE POLONOMIALS AND EPSNM
!              --------------------------------------




!     ------------------------------------------------------------------

!*       6.    UPDATE SPECTRAL ARRAYS.
!              -----------------------

CALL UPDSPAD(KM,KF_UV,KF_SCALARS,ZOA1,ZOA2, &
 & PSPVOR,PSPDIV,PSPSCALAR,&
 & PSPSC3A,PSPSC3B,PSPSC2 , &
 & KFLDPTRUV,KFLDPTRSC)

!     ------------------------------------------------------------------

!*       5.    COMPUTE VORTICITY AND DIVERGENCE.
!              ---------------------------------

IF( KF_UV > 0 ) THEN
  CALL PREPSNM(KM,KMLOC,ZEPSNM)
  IUS = 1
  IUE = 2*KF_UV
  IVS = 2*KF_UV+1
  IVE = 4*KF_UV
  IVORS = 1
  IVORE = 2*KF_UV
  IDIVS = 2*KF_UV+1
  IDIVE = 4*KF_UV
! SET PART OF ZOA1 CONTAINING U AND V TO 0.
  ZOA1(:,IUS:IVE) = 0.0_JPRB
  CALL UVTVDAD(KM,KF_UV,ZEPSNM,ZOA1(:,IUS:IUE),ZOA1(:,IVS:IVE),&
   & ZOA2(:,IVORS:IVORE),ZOA2(:,IDIVS:IDIVE))
ENDIF

!     ------------------------------------------------------------------

!*       4.    DIRECT LEGENDRE TRANSFORM.
!              --------------------------
IFC = 2*KF_FS
IDGLU = MIN(R%NDGNH,G%NDGLU(KM))
IIFC = IFC
IF(KM == 0)THEN
  IIFC = IFC/2
ENDIF
CALL LEDIRAD(KM,KMLOC,IFC,IIFC,IDGLU,KLED2,ZAIA,ZSIA,ZOA1)

!     ------------------------------------------------------------------

!*       3.    FOURIER SPACE COMPUTATIONS.
!              ---------------------------

CALL LDFOU2AD(KM,KF_UV,ZAIA,ZSIA)

!     ------------------------------------------------------------------

!*       2.    PREPARE WORK ARRAYS.
!              --------------------

CALL PRFI2AD(KM,KMLOC,KF_FS,ZAIA,ZSIA)


!     ------------------------------------------------------------------

END SUBROUTINE LTDIRAD
END MODULE LTDIRAD_MOD

