! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTINV_CTLAD_MOD
CONTAINS
SUBROUTINE FTINV_CTLAD(KF_UV_G,KF_SCALARS_G,&
 & KF_UV,KF_SCALARS,KF_SCDERS,KF_GP,KF_FS,KF_OUT_LT,KVSETUV,KVSETSC,KPTRGP, &
 & KVSETSC3A,KVSETSC3B,KVSETSC2,&
 & PGP,PGPUV,PGP3A,PGP3B,PGP2)


!**** *FTINV_CTLAD - Inverse Fourier transform control - adjoint

!     Purpose. Control routine for Fourier to Gridpoint transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINV_CTLAD(..)

!        Explicit arguments :
!        --------------------
!        PGP     -  gridpoint array
!        KF_UV_G      - global number of spectral u-v fields
!        KF_SCALARS_G - global number of scalar spectral fields
!        KF_UV        - local number of spectral u-v fields
!        KF_SCALARS   - local number of scalar spectral fields
!        KF_SCDERS    - local number of derivatives of scalar spectral fields
!        KF_GP        - total number of output gridpoint fields
!        KF_FS        - total number of fields in fourier space
!        KF_OUT_LT    - total number of fields coming out from inverse LT
!        KVSETUV - "B"  set in spectral/fourier space for
!                   u and v variables
!        KVSETSC - "B" set in spectral/fourier space for
!                  scalar variables
!        KPTRGP - pointer array to fi3elds in gridpoint space


!     Method.
!     -------

!     Externals.  TRLTOG      - transposition routine
!     ----------  FOURIER_IN  - copy fourier data from Fourier buffer
!                 FTINV       - fourier transform
!                 FSC         - Fourier space computations

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        R. El Khatib 09-Sep-2020 NSTACK_MEMORY_TR

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : NERR    ,NSTACK_MEMORY_TR
USE TPM_TRANS       ,ONLY : FOUBUF, LDIVGP, LSCDERS, LUVDER, LVORGP
USE TPM_DISTR       ,ONLY : D, MYPROC, NPROC

USE FOURIER_INAD_MOD ,ONLY : FOURIER_INAD
USE FSCAD_MOD       ,ONLY : FSCAD
USE FTINVAD_MOD     ,ONLY : FTINVAD
USE TRGTOL_MOD      ,ONLY : TRGTOL
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_UV_G
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_SCALARS
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_SCDERS
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_GP
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_FS
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_OUT_LT
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC2(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP2(:,:,:)

!     ------------------------------------------------------------------

REAL(KIND=JPRB),TARGET  :: ZGTF_STACK(KF_FS*MIN(1,MAX(0,NSTACK_MEMORY_TR)),D%NLENGTF)
REAL(KIND=JPRB),TARGET, ALLOCATABLE :: ZGTF_HEAP(:,:)
REAL(KIND=JPRB),POINTER  :: ZGTF(:,:)
REAL(KIND=JPRB),TARGET  :: ZDUM(1,D%NLENGTF)
REAL(KIND=JPRB),POINTER :: ZUV(:,:)
REAL(KIND=JPRB),POINTER :: ZSCALAR(:,:)
REAL(KIND=JPRB),POINTER :: ZNSDERS(:,:)
REAL(KIND=JPRB),POINTER :: ZEWDERS(:,:)
REAL(KIND=JPRB),POINTER :: ZUVDERS(:,:)

INTEGER(KIND=JPIM) :: IST,IBLEN
INTEGER(KIND=JPIM) :: IVSETUV(KF_UV_G)
INTEGER(KIND=JPIM) :: IVSETSC(KF_SCALARS_G)
INTEGER(KIND=JPIM) :: IVSET(KF_GP)
INTEGER(KIND=JPIM) :: J3,JGL,IGL,IOFF,IFGP2,IFGP3A,IFGP3B,IGP3APAR,IGP3BPAR
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC

!     ------------------------------------------------------------------

!   4.  Transposition

IF(PRESENT(KVSETUV)) THEN
  IVSETUV(:) = KVSETUV(:)
ELSE
  IVSETUV(:) = -1
ENDIF

IVSETSC(:)=-1
IF(PRESENT(KVSETSC)) THEN
  IVSETSC(:) = KVSETSC(:)
ELSE
  IOFF=0
  IF(PRESENT(KVSETSC2)) THEN
    IFGP2=UBOUND(KVSETSC2,1)
    IVSETSC(1:IFGP2)=KVSETSC2(:)
    IOFF=IOFF+IFGP2
  ENDIF
  IF(PRESENT(KVSETSC3A)) THEN
    IFGP3A=UBOUND(KVSETSC3A,1)
    IGP3APAR=UBOUND(PGP3A,3)
    IF(LSCDERS) IGP3APAR=IGP3APAR/3
    DO J3=1,IGP3APAR
      IVSETSC(IOFF+1:IOFF+IFGP3A)=KVSETSC3A(:)
      IOFF=IOFF+IFGP3A
    ENDDO
  ENDIF
  IF(PRESENT(KVSETSC3B)) THEN
    IFGP3B=UBOUND(KVSETSC3B,1)
    IGP3BPAR=UBOUND(PGP3B,3)
    IF(LSCDERS) IGP3BPAR=IGP3BPAR/3
    DO J3=1,IGP3BPAR
      IVSETSC(IOFF+1:IOFF+IFGP3B)=KVSETSC3B(:)
      IOFF=IOFF+IFGP3B
    ENDDO
  ENDIF
  IF(IOFF > 0 .AND. IOFF /= KF_SCALARS_G ) THEN
    WRITE(NERR,*)'FTINV_CTLAD:IOFF,KF_SCALARS_G ',IOFF,KF_SCALARS_G
    CALL ABORT_TRANS('FTINV_CTLAD_MOD:IOFF /= KF_SCALARS_G')
  ENDIF
ENDIF

IST = 1
IF(KF_UV_G > 0) THEN
  IF( LVORGP) THEN
    IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
    IST = IST+KF_UV_G
  ENDIF
  IF( LDIVGP) THEN
    IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
    IST = IST+KF_UV_G
  ENDIF
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
ENDIF
IF(KF_SCALARS_G > 0) THEN
  IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
  IST = IST+KF_SCALARS_G
  IF(LSCDERS) THEN
    IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
    IST = IST+KF_SCALARS_G
  ENDIF
ENDIF
IF(KF_UV_G > 0 .AND. LUVDER) THEN
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
ENDIF
IF(KF_SCALARS_G > 0) THEN
  IF(LSCDERS) THEN
    IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
    IST = IST+KF_SCALARS_G
  ENDIF
ENDIF

IF (NSTACK_MEMORY_TR == 1) THEN
  ZGTF => ZGTF_STACK(:,:)
ELSE
  ALLOCATE(ZGTF_HEAP(KF_FS,D%NLENGTF))
! Now, force the OS to allocate this shared array right now, not when it starts
! to be used which is an OPEN-MP loop, that would cause a threads synchronization lock :
  IF (KF_FS > 0 .AND. D%NLENGTF > 0) THEN
    ZGTF_HEAP(1,1)=HUGE(1._JPRB)
  ENDIF
  ZGTF => ZGTF_HEAP(:,:)
ENDIF

CALL GSTATS(182,0)
CALL TRGTOL(ZGTF,KF_FS,KF_GP,KF_SCALARS_G,IVSET,KPTRGP,&
 &PGP,PGPUV,PGP3A,PGP3B,PGP2)
CALL GSTATS(182,1)

!   3.  Fourier transform

IF(KF_UV > 0 .OR. KF_SCDERS > 0) THEN
  IST = 1
  IF(LVORGP) THEN
    IST = IST+KF_UV
  ENDIF
  IF(LDIVGP) THEN
    IST = IST+KF_UV
  ENDIF
  IF(KF_UV>0) THEN
    ZUV => ZGTF(IST:IST+2*KF_UV-1,:)
  ELSE
    ZUV => ZDUM(1:1,:)
  ENDIF
  IST = IST+2*KF_UV
  IF(KF_SCALARS>0) THEN
    ZSCALAR => ZGTF(IST:IST+KF_SCALARS-1,:)
  ELSE
    ZSCALAR => ZDUM(1:1,:)
  ENDIF
  IST = IST+KF_SCALARS
  IF(KF_SCDERS>0) THEN
    ZNSDERS => ZGTF(IST:IST+KF_SCDERS-1,:)
  ELSE
    ZNSDERS => ZDUM(1:1,:)
  ENDIF
  IST = IST+KF_SCDERS
  IF(LUVDER) THEN
    ZUVDERS => ZGTF(IST:IST+2*KF_UV-1,:)
    IST = IST+2*KF_UV
  ELSE
    ZUVDERS => ZDUM(1:1,:)
  ENDIF
  IF(KF_SCDERS > 0) THEN
    ZEWDERS => ZGTF(IST:IST+KF_SCDERS-1,:)
  ELSE
    ZEWDERS => ZDUM(1:1,:)
  ENDIF
ENDIF

IBLEN = D%NLENGT0B*2*KF_OUT_LT
IF (ALLOCATED(FOUBUF)) THEN
  IF (MAX(1,IBLEN) > SIZE(FOUBUF)) THEN
    DEALLOCATE(FOUBUF)
    ALLOCATE(FOUBUF(MAX(1,IBLEN)))
  ENDIF
ELSE
  ALLOCATE(FOUBUF(MAX(1,IBLEN)))
ENDIF

CALL GSTATS(132,0)

IF(MYPROC > NPROC/2)THEN
  IBEG=1
  IEND=D%NDGL_FS
  IINC=1
ELSE
  IBEG=D%NDGL_FS
  IEND=1
  IINC=-1
ENDIF

CALL GSTATS(1641,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JGL,IGL)
DO JGL=IBEG,IEND,IINC
  IGL = JGL
  IF(KF_FS > 0) THEN
    CALL FTINVAD(ZGTF,KF_FS,IGL)
  ENDIF

!   2.  Fourier space computations

  IF(KF_UV > 0 .OR. KF_SCDERS > 0) THEN
    CALL FSCAD(IGL,KF_UV,KF_SCALARS,KF_SCDERS,&
     & ZUV,ZSCALAR,ZNSDERS,ZEWDERS,ZUVDERS)
  ENDIF

!   1. Copy Fourier data to local array

  CALL FOURIER_INAD(ZGTF,KF_OUT_LT,IGL)

ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1641,1)

IF(KF_UV > 0 .OR. KF_SCDERS > 0) THEN
  NULLIFY(ZUV)
  NULLIFY(ZSCALAR)
  NULLIFY(ZNSDERS)
  NULLIFY(ZUVDERS)
  NULLIFY(ZEWDERS)
ENDIF

CALL GSTATS(132,1)

!     ------------------------------------------------------------------

END SUBROUTINE FTINV_CTLAD
END MODULE FTINV_CTLAD_MOD



