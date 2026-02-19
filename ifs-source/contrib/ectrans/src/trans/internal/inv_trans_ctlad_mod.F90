! (C) Copyright 2001- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE INV_TRANS_CTLAD_MOD
CONTAINS
SUBROUTINE INV_TRANS_CTLAD(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_OUT_LT,&
 & KF_UV,KF_SCALARS,KF_SCDERS,&
 & PSPVOR,PSPDIV,PSPSCALAR,KVSETUV,KVSETSC,PGP,&
 & PSPSC3A,PSPSC3B,PSPSC2,KVSETSC3A,KVSETSC3B,KVSETSC2,PGPUV,PGP3A,PGP3B,PGP2)

!**** *INV_TRANS_CTLAD* - Control routine for inverse spectral transform adj.

!     Purpose.
!     --------
!        Control routine for the inverse spectral transform

!**   Interface.
!     ----------
!     CALL INV_TRANS_CTLAD(...)

!     Explicit arguments :
!     --------------------
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     KF_OUT_LT    - total number of fields coming out from inverse LT
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     KF_SCDERS    - local number of derivatives of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity (input)
!     PSPDIV(:,:)  - spectral divergence (input)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
!     KVSETUV(:)  - indicating which 'b-set' in spectral space owns a
!                   vor/div field. Equivalant to NBSETLEV in the IFS.
!                   The length of KVSETUV should be the GLOBAL number
!                   of u/v fields which is the dimension of u and v releated
!                   fields in grid-point space.
!     KVESETSC(:) - indicating which 'b-set' in spectral space owns a
!                   scalar field. As for KVSETUV this argument is required
!                   if the total number of processors is greater than
!                   the number of processors used for distribution in
!                   spectral wave space.
!     PGP(:,:,:)  - gridpoint fields (output)

!                  The ordering of the output fields is as follows (all
!                  parts are optional depending on the input switches):
!
!       vorticity     : KF_UV_G fields
!       divergence    : KF_UV_G fields
!       u             : KF_UV_G fields
!       v             : KF_UV_G fields
!       scalar fields : KF_SCALARS_G fields
!       N-S derivative of scalar fields : KF_SCALARS_G fields
!       E-W derivative of u : KF_UV_G fields
!       E-W derivative of v : KF_UV_G fields
!       E-W derivative of scalar fields : KF_SCALARS_G fields
!
!     Method.
!     -------

!     Externals.  SHUFFLE     - reshuffle fields for load balancing
!     ----------  FIELD_SPLIT - split fields in NPROMATR packets
!                 LTINV_CTLAD   - control of Legendre transform
!                 FTINV_CTLAD   - control of Fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 01-01-03

!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : NPROMATR
USE TPM_TRANS       ,ONLY : LDIVGP, LSCDERS, LUVDER, LVORGP
!USE TPM_DISTR

USE SHUFFLE_MOD     ,ONLY : SHUFFLE
USE FIELD_SPLIT_MOD ,ONLY : FIELD_SPLIT
USE LTINV_CTLAD_MOD ,ONLY : LTINV_CTLAD
USE FTINV_CTLAD_MOD ,ONLY : FTINV_CTLAD
!

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_GP
INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_OUT_LT
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCDERS
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC2(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC2(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGPUV(:,:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP3A(:,:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP3B(:,:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP2(:,:,:)

! Local variables

INTEGER(KIND=JPIM) :: IPTRGP(KF_GP),IPTRSPUV(NPROMATR),IPTRSPSC(NPROMATR)
INTEGER(KIND=JPIM) :: ISHFUV_G(KF_GP),ISHFSC_G(KF_GP)
INTEGER(KIND=JPIM) :: IVSETUV(KF_GP),IVSETSC(KF_GP)
INTEGER(KIND=JPIM) :: IBLKS,JBLK,ISTUV_G,IENUV_G
INTEGER(KIND=JPIM) :: IF_UV_G,IF_UV,ISTUV,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP
INTEGER(KIND=JPIM) :: JFLD,ISTSC_G,IENSC_G,ISTSC,IENSC,IENUV,IF_SCDERS,IF_OUT_LT
INTEGER(KIND=JPIM) :: IOFFD,IOFFU,IOFFV,IOFFUVD,IOFFSC,IOFFSCNS,IOFFSCEW,IOFF,IF_GPB

!     ------------------------------------------------------------------

! Perform transform


IF_GPB = 2*KF_UV_G+KF_SCALARS_G
IF(NPROMATR > 0 .AND. IF_GPB > NPROMATR) THEN

  ! Fields to be split into packets

  CALL SHUFFLE(KF_UV_G,KF_SCALARS_G,ISHFUV_G,IVSETUV,ISHFSC_G,IVSETSC,&
 & KVSETUV,KVSETSC)

  IBLKS=(IF_GPB-1)/NPROMATR+1

  DO JBLK=1,IBLKS
  
    CALL FIELD_SPLIT(JBLK,IF_GPB,KF_UV_G,IVSETUV,IVSETSC,&
     & ISTUV_G,IENUV_G,IF_UV_G,ISTSC_G,IENSC_G,IF_SCALARS_G,&
     & ISTUV,IENUV,IF_UV,ISTSC,IENSC,IF_SCALARS)

    IF(LSCDERS) THEN
      IF_SCDERS = IF_SCALARS
    ELSE
      IF_SCDERS = 0
    ENDIF
        
    IF_OUT_LT = 2*IF_UV + IF_SCALARS+IF_SCDERS
    IF(LVORGP) THEN
      IF_OUT_LT = IF_OUT_LT+IF_UV
    ENDIF
    IF(LDIVGP) THEN
      IF_OUT_LT = IF_OUT_LT+IF_UV
    ENDIF
    IF_FS = IF_OUT_LT+IF_SCDERS
    IF(LUVDER) THEN
      IF_FS = IF_FS+2*IF_UV
    ENDIF

    IF_GP = 2*IF_UV_G+IF_SCALARS_G
    IOFFD = 0
    IOFFU = 0
    IOFFV = KF_UV_G
    IOFFUVD = 2*KF_UV_G+KF_SCALARS_G
    IOFFSC = 2*KF_UV_G
    IF(LVORGP) THEN
      IF_GP = IF_GP+IF_UV_G
      IOFFD = KF_UV_G
      IOFFU = IOFFU+KF_UV_G
      IOFFV = IOFFV+KF_UV_G
      IOFFUVD =IOFFUVD+KF_UV_G
      IOFFSC = IOFFSC+KF_UV_G
    ENDIF
    IF(LDIVGP) THEN
      IF_GP = IF_GP+IF_UV_G
      IOFFU = IOFFU+KF_UV_G
      IOFFV = IOFFV+KF_UV_G
      IOFFUVD =IOFFUVD+KF_UV_G
      IOFFSC = IOFFSC+KF_UV_G
    ENDIF
    IF(LSCDERS) THEN
      IF_GP  = IF_GP+2*IF_SCALARS_G
      IOFFUVD =IOFFUVD+KF_SCALARS_G
      IOFFSCNS = IOFFSC+KF_SCALARS_G
      IOFFSCEW = IOFFSC+2*KF_SCALARS_G
    ENDIF
    IF(LUVDER) THEN
      IF_GP = IF_GP+2*IF_UV_G
      IOFFSCEW = IOFFSCEW+2*KF_UV_G
    ENDIF

    DO JFLD=1,IF_UV_G
      IOFF = 0
      IF(LVORGP) THEN
        IPTRGP(JFLD+IOFF) = ISHFUV_G(ISTUV_G+JFLD-1)
        IOFF = IOFF+IF_UV_G
      ENDIF
      IF(LDIVGP) THEN
        IPTRGP(JFLD+IOFF) = IOFFD+ISHFUV_G(ISTUV_G+JFLD-1)
        IOFF = IOFF+IF_UV_G
      ENDIF
      IPTRGP(JFLD+IOFF) = IOFFU+ISHFUV_G(ISTUV_G+JFLD-1)
      IOFF = IOFF+IF_UV_G
      IPTRGP(JFLD+IOFF) = IOFFV+ISHFUV_G(ISTUV_G+JFLD-1)
      IOFF = IOFF+IF_UV_G+IF_SCALARS_G
      IF(LSCDERS) THEN
        IOFF  = IOFF+IF_SCALARS_G
      ENDIF
      IF(LUVDER) THEN
        IPTRGP(JFLD+IOFF) = IOFFUVD+ISHFUV_G(ISTUV_G+JFLD-1)
        IOFF = IOFF+IF_UV_G
        IPTRGP(JFLD+IOFF) = IOFFUVD+KF_UV_G+ISHFUV_G(ISTUV_G+JFLD-1)
      ENDIF
    ENDDO

    DO JFLD=1,IF_SCALARS_G
      IOFF = 2*IF_UV_G
      IF (LVORGP) IOFF = IOFF+IF_UV_G
      IF (LDIVGP) IOFF = IOFF+IF_UV_G
      IPTRGP(JFLD+IOFF) = IOFFSC+ISHFSC_G(ISTSC_G+JFLD-1)
      IOFF = IOFF+IF_SCALARS_G
      IF(LSCDERS) THEN
        IPTRGP(JFLD+IOFF) = IOFFSCNS+ISHFSC_G(ISTSC_G+JFLD-1)
        IOFF  = IOFF+IF_SCALARS_G
        IF(LUVDER) THEN
          IOFF  = IOFF+2*IF_UV_G
        ENDIF
        IPTRGP(JFLD+IOFF) = IOFFSCEW+ISHFSC_G(ISTSC_G+JFLD-1)
      ENDIF
    ENDDO
    DO JFLD=1,IF_UV
      IPTRSPUV(JFLD) = ISTUV+JFLD-1
    ENDDO
    DO JFLD=1,IF_SCALARS
      IPTRSPSC(JFLD) = ISTSC+JFLD-1
    ENDDO

    IF(IF_UV_G > 0 .AND. IF_SCALARS_G > 0) THEN
      CALL FTINV_CTLAD(IF_UV_G,IF_SCALARS_G,&
       & IF_UV,IF_SCALARS,IF_SCDERS,IF_GP,IF_FS,IF_OUT_LT,&
       & KVSETUV=IVSETUV(ISTUV_G:IENUV_G),&
       & KVSETSC=IVSETSC(ISTSC_G:IENSC_G),KPTRGP=IPTRGP,&
       & PGP=PGP)
   ELSEIF(IF_UV_G > 0) THEN
      CALL FTINV_CTLAD(IF_UV_G,IF_SCALARS_G,&
       & IF_UV,IF_SCALARS,IF_SCDERS,IF_GP,IF_FS,IF_OUT_LT,&
       & KVSETUV=IVSETUV(ISTUV_G:IENUV_G),&
       & KPTRGP=IPTRGP,&
       & PGP=PGP)
    ELSEIF(IF_SCALARS_G > 0) THEN
      CALL FTINV_CTLAD(IF_UV_G,IF_SCALARS_G,&
       & IF_UV,IF_SCALARS,IF_SCDERS,IF_GP,IF_FS,IF_OUT_LT,&
       & KVSETSC=IVSETSC(ISTSC_G:IENSC_G),KPTRGP=IPTRGP,&
       & PGP=PGP)
    ENDIF
    CALL LTINV_CTLAD(IF_OUT_LT,IF_UV,IF_SCALARS,IF_SCDERS, &
     & PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPSCALAR=PSPSCALAR,&
     & KFLDPTRUV=IPTRSPUV,KFLDPTRSC=IPTRSPSC)
    
  ENDDO

ELSE

  ! No splitting of fields, transform done in one go

  CALL FTINV_CTLAD(KF_UV_G,KF_SCALARS_G,&
   & KF_UV,KF_SCALARS,KF_SCDERS,KF_GP,KF_FS,KF_OUT_LT,&
   & KVSETUV=KVSETUV,KVSETSC=KVSETSC,&
   & KVSETSC3A=KVSETSC3A,KVSETSC3B=KVSETSC3B,KVSETSC2=KVSETSC2,&
   & PGP=PGP,PGPUV=PGPUV,PGP3A=PGP3A,PGP3B=PGP3B,PGP2=PGP2)

  CALL LTINV_CTLAD(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS, &
   &PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPSCALAR=PSPSCALAR,&
   &PSPSC3A=PSPSC3A,PSPSC3B=PSPSC3B,PSPSC2=PSPSC2)
 
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE INV_TRANS_CTLAD
END MODULE INV_TRANS_CTLAD_MOD
