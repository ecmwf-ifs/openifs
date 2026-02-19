! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE SPECTRAL_FIELDS_PARA_MOD
USE PARKIND1, ONLY : JPRB,JPIM
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK

!   Purpose.
!   --------
!     This module contains the distributed routines associated
!     to the SPECTRAL_FIELD type.
!
!   Author.
!   -------
!     Yannick Tremolet
!
!   Modifications.
!   --------------
!     Original   12-Mar-2004
! ------------------------------------------------------------------

USE SPECTRAL_FIELDS_DATA, ONLY: SPECTRAL_FIELD, NPRTRV, NPRTRW, LELAM, NULOUT
USE MPL_MODULE

IMPLICIT NONE
PRIVATE
PUBLIC DOT_PRODUCT, SPECTRAL_NORM, SPECTRAL_NORM_LEVS, SPECTRAL_MAX_LEVS

INTERFACE DOT_PRODUCT
MODULE PROCEDURE DOT_PRODUCT_SPEC
END INTERFACE

INTERFACE SPECTRAL_NORM
MODULE PROCEDURE SPECNORM_PRT, SPECNORM_ARRAY, SPECNORM_ARRPRT
END INTERFACE

INTERFACE SPECTRAL_NORM_LEVS
MODULE PROCEDURE SPECNORM_PRT_LEV, SPECNORM_ARRAY_LEV
END INTERFACE

INTERFACE SPECTRAL_MAX_LEVS
MODULE PROCEDURE SPECMAX_PRT_LEV, SPECMAX_ARRAY_LEV
END INTERFACE

! ------------------------------------------------------------------
CONTAINS
! ------------------------------------------------------------------

REAL(KIND=JPRB) FUNCTION DOT_PRODUCT_SPEC(YDSP1,YDSP2,LDWEIGHT)

!   Purpose.
!   --------
!     Scalar product for distributed spectral variable
!                lweight : .T. => spectral inner product (weighted with 2's and 1's)
!                lweight : .F. (or missing) => simple, unweighted dot product
!
!   Author.
!   -------
!     Yannick Tremolet
!
!   Modifications.
!   --------------
!     Original   02-03-19
!     C.Fischer  03-03-05 Aladin spectral treatment added.
!     C. Soci    04-10-07 CY29 cleaning
! ------------------------------------------------------------------

IMPLICIT NONE

TYPE(SPECTRAL_FIELD), INTENT(IN) :: YDSP1, YDSP2
LOGICAL, OPTIONAL, INTENT(IN)    :: LDWEIGHT

INTEGER(KIND=JPIM) :: II,JB,JL,IL,IFLDL,IFLDG,INUMLEV(NPRTRV),JJ
REAL(KIND=JPRB), ALLOCATABLE :: ZBUFL(:),ZBUFG(:),ZPROD3(:)
REAL(KIND=JPRB) :: ZPROD2, ZZPROD
LOGICAL :: LLWEIGHT

! ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:DOT_PRODUCT_SPEC',0,ZHOOK_HANDLE)
CALL GSTATS(671,0)

IF (PRESENT(LDWEIGHT)) THEN
  LLWEIGHT=LDWEIGHT
ELSE
  LLWEIGHT=.FALSE.
ENDIF

!yt IF (ydsp1.NEQV.ydsp2) CALL ABOR1('Spectral fields are not compatible')
IF ((YDSP1%NSMAX/=YDSP2%NSMAX).OR.((YDSP1%NMSMAX/=YDSP2%NMSMAX).AND.LELAM)) &
    & CALL ABOR1('Scalar product requires fields at the same resolution')

IFLDL=YDSP1%NFLEVL
IF (YDSP1%NS2D>0) IFLDL=IFLDL+1
ALLOCATE(ZBUFL(IFLDL))
IFLDG=YDSP1%NFLEVG
IF (YDSP1%NS2G>0) IFLDG=IFLDG+1
ALLOCATE(ZBUFG(IFLDG))
ALLOCATE(ZPROD3(YDSP1%NFLEVG))
INUMLEV(:)=YDSP1%NUMLL(:)
WHERE (YDSP1%NS2L(:)>0) INUMLEV(:)=INUMLEV(:)+1

! Scalar product by level
! -----------------------
CALL DOT_PROD_WAVES(YDSP1,YDSP2,ZBUFL,LLWEIGHT)

! Gather values for all levels
! ----------------------------
CALL MPL_ALLGATHERV(PSENDBUF=ZBUFL(1:IFLDL),PRECVBUF=ZBUFG,&
                  & KRECVCOUNTS=INUMLEV,KCOMM=MPL_ALL_LEVS_COMM, &
                  & CDSTRING='SPECTRAL_FIELDS_PARA_MOD:DOT_PRODUCT_SPEC')

! Reorder levels for reproducibility
! ----------------------------------
II=0
DO JB=1,NPRTRV
  DO JL=1,YDSP1%NUMLL(JB)
    IL=YDSP1%NPTRLL(JB)+JL-1
    II=II+1
    ZPROD3(IL)=ZBUFG(II)
  ENDDO
  IF (YDSP1%NS2L(JB)>0) THEN
    II=II+1
    ZPROD2=ZBUFG(II)
  ENDIF
ENDDO

! Sum all levels for 3D fields
! ----------------------------
ZZPROD=0.0_JPRB
DO JL=1,YDSP1%NFLEVG
  ZZPROD=ZZPROD+ZPROD3(JL)
ENDDO
IF (YDSP1%NS2G>0) THEN
  ZZPROD=ZZPROD+ZPROD2
ENDIF

DEALLOCATE(ZBUFL)
DEALLOCATE(ZBUFG)
DEALLOCATE(ZPROD3)

DOT_PRODUCT_SPEC = ZZPROD

CALL GSTATS(671,1)
IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:DOT_PRODUCT_SPEC',1,ZHOOK_HANDLE)
END FUNCTION DOT_PRODUCT_SPEC

! ------------------------------------------------------------------

SUBROUTINE DOT_PROD_WAVES(YDSP1,YDSP2,PRODSP,LDWEIGHT)

TYPE(SPECTRAL_FIELD), INTENT(IN) :: YDSP1, YDSP2
REAL(KIND=JPRB), INTENT(OUT) :: PRODSP(:)
LOGICAL, INTENT(IN) :: LDWEIGHT

INTEGER(KIND=JPIM) :: JF, JM, JN, JL, JA, JS, IFLDS, IM, ILEN, ILENS, ILENG, JJ, J3D
INTEGER(KIND=JPIM) :: ILENR(NPRTRW), IIISMAX
REAL(KIND=JPRB), ALLOCATABLE :: ZBUFS(:), ZBUFR(:)
REAL(KIND=JPRB), ALLOCATABLE :: ZALLMS(:,:)

! ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:DOT_PROD_WAVES',0,ZHOOK_HANDLE)
IFLDS=YDSP1%NFLEVL
IF (YDSP1%NS2D>0) IFLDS=IFLDS+1
IF (SIZE(PRODSP)/=IFLDS) CALL ABOR1('prodsp wrong size')

ILENS=IFLDS*YDSP1%NUMP
ILENR(:)=IFLDS*YDSP1%NUMPP(:)
ILENG=SUM(ILENR)
ALLOCATE(ZBUFS(ILENS))
ALLOCATE(ZBUFR(ILENG))
IF (LELAM) THEN
  IIISMAX=YDSP1%NMSMAX
ELSE
  IIISMAX=YDSP1%NSMAX
ENDIF
ALLOCATE(ZALLMS(IFLDS,0:IIISMAX))
ZBUFS(:)=0.0_JPRB

! Partial products for local Ms in send buffer
! --------------------------------------------
JL=YDSP1%NFLEVL
IF (.NOT.LELAM) THEN
  DO JF=1,YDSP1%NS3D
    JJ=0
    DO JM=1,YDSP1%NUMP
      IM=YDSP1%MYMS(JM)
      JS=YDSP1%NASM0(IM)-1
      ILEN=2*(YDSP1%NSMAX-IM+1)
      DO JN=1,ILEN
        ZBUFS(JJ+1:JJ+JL) = ZBUFS(JJ+1:JJ+JL) &
            & + YDSP1%SP3D(:,JS+JN,JF) * YDSP2%SP3D(:,JS+JN,JF)
      ENDDO
      JJ=JJ+JL
    ENDDO
  ENDDO
ELSE
  DO JF=1,YDSP1%NS3D
    JJ=0
    DO JM=1,YDSP1%NUMP
      IM=YDSP1%MYMS(JM)
      JS=YDSP1%NESM0(IM)-1
      ILEN=YDSP1%NCPL4M(IM)
      DO JN=1,ILEN
        ZBUFS(JJ+1:JJ+JL) = ZBUFS(JJ+1:JJ+JL) &
            & + YDSP1%SP3D(:,JS+JN,JF) * YDSP2%SP3D(:,JS+JN,JF)
      ENDDO
      JJ=JJ+JL
    ENDDO
  ENDDO
ENDIF

J3D=YDSP1%NFLEVL*YDSP1%NUMP
IF (.NOT.LELAM) THEN
  DO JF=1,YDSP1%NS2D
    DO JM=1,YDSP1%NUMP
      IM=YDSP1%MYMS(JM)
      JS=YDSP1%NASM0(IM)-1
      ILEN=2*(YDSP1%NSMAX-IM+1)
      DO JN=1,ILEN
        ZBUFS(J3D+JM) = ZBUFS(J3D+JM) &
            & + YDSP1%SP2D(JS+JN,JF) * YDSP2%SP2D(JS+JN,JF)
      ENDDO
    ENDDO
  ENDDO
ELSE
  DO JF=1,YDSP1%NS2D
    DO JM=1,YDSP1%NUMP
      IM=YDSP1%MYMS(JM)
      JS=YDSP1%NESM0(IM)-1
      ILEN=YDSP1%NCPL4M(IM)
      DO JN=1,ILEN
        ZBUFS(J3D+JM) = ZBUFS(J3D+JM) &
            & + YDSP1%SP2D(JS+JN,JF) * YDSP2%SP2D(JS+JN,JF)
      ENDDO
    ENDDO
  ENDDO
ENDIF
! Communicate to get all the Ms
! -----------------------------
CALL MPL_ALLGATHERV(PSENDBUF=ZBUFS(1:ILENS),PRECVBUF=ZBUFR,&
                  & KRECVCOUNTS=ILENR,KCOMM=MPL_ALL_MS_COMM, &
                  & CDSTRING='SPECTRAL_FIELDS_PARA_MOD:DOT_PROD_WAVES')

! Reorder Ms for reproducibility
! ------------------------------
JJ=0
DO JA=1,NPRTRW
  DO JM=1,YDSP1%NUMPP(JA)
    IM=YDSP1%MS_PROC(JM,JA)
    DO JL=1,YDSP1%NFLEVL
      JJ=JJ+1
      ZALLMS(JL,IM)=ZBUFR(JJ)
    ENDDO
  ENDDO
  IF (YDSP1%NS2D>0) THEN
    JL=YDSP1%NFLEVL+1
    DO JM=1,YDSP1%NUMPP(JA)
      IM=YDSP1%MS_PROC(JM,JA)
      JJ=JJ+1
      ZALLMS(JL,IM)=ZBUFR(JJ)
    ENDDO
  ENDIF
ENDDO

! Total scalar products
! ---------------------

PRODSP(:) = 0.0_JPRB

IF (LDWEIGHT) THEN
  DO JM=1,IIISMAX
    PRODSP(:) = PRODSP(:) + ZALLMS(:,JM)
  ENDDO
  
  PRODSP(:) = 2.0_JPRB*PRODSP(:) + ZALLMS(:,0)
ELSE
  DO JM=0,IIISMAX
    PRODSP(:) = PRODSP(:) + ZALLMS(:,JM)
  ENDDO
ENDIF

DEALLOCATE(ZBUFS)
DEALLOCATE(ZBUFR)
DEALLOCATE(ZALLMS)

IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:DOT_PROD_WAVES',1,ZHOOK_HANDLE)
END SUBROUTINE DOT_PROD_WAVES

! ------------------------------------------------------------------
SUBROUTINE SPECNORM_PRT(YDSPEC,CDGREP)
! ------------------------------------------------------------------
! Computes norms of distributed spectral fields and prints
! one value per field (3D or 2D). Adds string to output for
! easier grep. Equivalent of SPNORM for type SPECTRAL_FIELD.
! ------------------------------------------------------------------

IMPLICIT NONE
TYPE(SPECTRAL_FIELD), INTENT(IN) :: YDSPEC
CHARACTER(LEN=*), INTENT(IN) :: CDGREP

REAL(KIND=JPRB), ALLOCATABLE :: ZNOR2(:),ZNOR3(:)
CHARACTER(LEN=32) :: CLFMT

! ----------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:SPECNORM_PRT',0,ZHOOK_HANDLE)
ALLOCATE(ZNOR3(YDSPEC%NS3D))
ALLOCATE(ZNOR2(YDSPEC%NS2G))

CALL SPECNORM(YDSPEC,ZNOR3,ZNOR2)

CLFMT='(2(X,A),5(X,ES18.12))'
IF (YDSPEC%NS2G>0) WRITE(NULOUT,CLFMT)CDGREP,'SP NORMS 2D:',ZNOR2
IF (YDSPEC%NS3D>0) WRITE(NULOUT,CLFMT)CDGREP,'SP NORMS 3D:',ZNOR3

DEALLOCATE(ZNOR2)
DEALLOCATE(ZNOR3)

IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:SPECNORM_PRT',1,ZHOOK_HANDLE)
END SUBROUTINE SPECNORM_PRT

! ------------------------------------------------------------------
SUBROUTINE SPECNORM_ARRAY(YDSPEC,PNORMS)
! ------------------------------------------------------------------
! Computes norms of distributed spectral fields and returns
! one value per field (3D or 2D) in array pnorms.
! ------------------------------------------------------------------

IMPLICIT NONE
TYPE(SPECTRAL_FIELD), INTENT(IN) :: YDSPEC
REAL(KIND=JPRB), INTENT(OUT) :: PNORMS(:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:SPECNORM_ARRAY',0,ZHOOK_HANDLE)
IF (SIZE(PNORMS)/=YDSPEC%NS3D+YDSPEC%NS2G) &
  & CALL ABOR1('SPECNORM_ARRAY: pnorms wrong size')

CALL SPECNORM(YDSPEC,PNORMS(1:YDSPEC%NS3D), &
                  & PNORMS(YDSPEC%NS3D+1:YDSPEC%NS3D+YDSPEC%NS2G))

IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:SPECNORM_ARRAY',1,ZHOOK_HANDLE)
END SUBROUTINE SPECNORM_ARRAY

! ------------------------------------------------------------------
SUBROUTINE SPECNORM_ARRPRT(YDSPEC,CDGREP,PNORMS)
! ------------------------------------------------------------------
! Computes norms of distributed spectral fields and prints one
! value per field (3D or 2D) and returns values in array pnorms.
! ------------------------------------------------------------------

IMPLICIT NONE
TYPE(SPECTRAL_FIELD), INTENT(IN) :: YDSPEC
CHARACTER(LEN=*), INTENT(IN) :: CDGREP
REAL(KIND=JPRB), INTENT(OUT) :: PNORMS(:)

CHARACTER(LEN=32) :: CLFMT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:SPECNORM_ARRPRT',0,ZHOOK_HANDLE)
IF (SIZE(PNORMS)/=YDSPEC%NS3D+YDSPEC%NS2G) &
  & CALL ABOR1('SPECNORM_ARRPRT: pnorms wrong size')

CALL SPECNORM(YDSPEC,PNORMS(1:YDSPEC%NS3D), &
                  & PNORMS(YDSPEC%NS3D+1:YDSPEC%NS3D+YDSPEC%NS2G))

CLFMT='(2(X,A),5(X,ES18.12))'
IF (YDSPEC%NS2G>0) WRITE(NULOUT,CLFMT)CDGREP,'SP NORMS 2D:', &
                       & PNORMS(YDSPEC%NS3D+1:YDSPEC%NS3D+YDSPEC%NS2G)
IF (YDSPEC%NS3D>0) WRITE(NULOUT,CLFMT)CDGREP,'SP NORMS 3D:',PNORMS(1:YDSPEC%NS3D)

IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:SPECNORM_ARRPRT',1,ZHOOK_HANDLE)
END SUBROUTINE SPECNORM_ARRPRT

! ------------------------------------------------------------------
SUBROUTINE SPECNORM_PRT_LEV(YDSPEC,CDGREP)
! ------------------------------------------------------------------
! Computes norms of distributed spectral fields and prints
! one value per 2D field and one value per level for 3D fields.
! Adds string to output for easier grep.
! ------------------------------------------------------------------

IMPLICIT NONE
TYPE(SPECTRAL_FIELD), INTENT(IN) :: YDSPEC
CHARACTER(LEN=*), INTENT(IN) :: CDGREP

INTEGER(KIND=JPIM) :: JL
REAL(KIND=JPRB), ALLOCATABLE :: ZNOR2(:),ZNOR3(:,:)
CHARACTER(LEN=32) :: CLFMT

! ----------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:SPECNORM_PRT_LEV',0,ZHOOK_HANDLE)
ALLOCATE(ZNOR3(YDSPEC%NFLEVG,YDSPEC%NS3D))
ALLOCATE(ZNOR2(YDSPEC%NS2G))

CALL SPECNORM_ARRAY_LEV(YDSPEC,ZNOR3,ZNOR2)

IF (YDSPEC%NS2G>0) THEN
  CLFMT='(2(X,A),5(X,ES18.12))'
  WRITE(NULOUT,CLFMT)CDGREP,'SP NORMS 2D:',ZNOR2
ENDIF

IF (YDSPEC%NS3D>0) THEN
  CLFMT='(2(X,A),X,I3,5(X,ES18.12))'
  DO JL=1,YDSPEC%NFLEVG
    WRITE(NULOUT,CLFMT)CDGREP,'SP NORMS 3D:',JL,ZNOR3(JL,:)
  ENDDO
ENDIF

DEALLOCATE(ZNOR2)
DEALLOCATE(ZNOR3)

IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:SPECNORM_PRT_LEV',1,ZHOOK_HANDLE)
END SUBROUTINE SPECNORM_PRT_LEV

! ------------------------------------------------------------------
SUBROUTINE SPECNORM(YDSPEC,PROD3,PROD2)
! ------------------------------------------------------------------
! Performs norms computation
! ------------------------------------------------------------------

IMPLICIT NONE

TYPE(SPECTRAL_FIELD), INTENT(IN) :: YDSPEC
REAL(KIND=JPRB), INTENT(OUT) :: PROD2(:),PROD3(:)

INTEGER(KIND=JPIM) :: JF,JL
REAL(KIND=JPRB), ALLOCATABLE :: ZPROD3(:,:)
REAL(KIND=JPRB) :: ZLEV

! ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:SPECNORM',0,ZHOOK_HANDLE)
IF (SIZE(PROD2)/=YDSPEC%NS2G) CALL ABOR1('SPECNORM: prod2 wrong size')
IF (SIZE(PROD3)/=YDSPEC%NS3D) CALL ABOR1('SPECNORM: prod3 wrong size')

ALLOCATE(ZPROD3(YDSPEC%NFLEVG,YDSPEC%NS3D))

CALL SPECNORM_ARRAY_LEV(YDSPEC,ZPROD3,PROD2)

! Average levels for 3D fields
! ----------------------------
ZLEV=REAL(YDSPEC%NFLEVG,JPRB)
DO JF=1,YDSPEC%NS3D
  PROD3(JF)=0.0_JPRB
  DO JL=1,YDSPEC%NFLEVG
    PROD3(JF)=PROD3(JF)+ZPROD3(JL,JF)
  ENDDO
  PROD3(JF)=PROD3(JF)/ZLEV
ENDDO

DEALLOCATE(ZPROD3)

IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:SPECNORM',1,ZHOOK_HANDLE)
END SUBROUTINE SPECNORM

! ------------------------------------------------------------------
SUBROUTINE SPECNORM_ARRAY_LEV(YDSPEC,PROD3,PROD2)
! ------------------------------------------------------------------
!   Purpose.
!   --------
!     Actual computation of norms of distributed spectral fields.
!     Returns one value per field in prod2 for 2D fields, one value
!     per fields and per level in prod3 for 3D norms.
!
!   Author.
!   -------
!     Yannick Tremolet
!
!   Modifications.
!   --------------
!     Original   11-Mar-2004
!     C. Soci    07-Oct-2004 Aladin spectral treatment added 
! ------------------------------------------------------------------

IMPLICIT NONE

TYPE(SPECTRAL_FIELD), INTENT(IN) :: YDSPEC
REAL(KIND=JPRB), INTENT(OUT) :: PROD2(:)
REAL(KIND=JPRB), INTENT(OUT) :: PROD3(:,:)

INTEGER(KIND=JPIM) :: JA,JB,JF,JJ,JL,JM,JN,JS
INTEGER(KIND=JPIM) :: II,IL,IM,IFLDL,IFLDG,ILEN,ILENS,ILENG,IIISMAX
INTEGER(KIND=JPIM) :: INUMLEV(NPRTRV),ILENR(NPRTRW)
REAL(KIND=JPRB), ALLOCATABLE :: ZBUFL(:),ZBUFG(:), ZBUFS(:), ZBUFR(:)
REAL(KIND=JPRB), ALLOCATABLE :: ZALLMS(:,:)

! ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:SPECNORM_ARRAY_LEV',0,ZHOOK_HANDLE)
IF (SIZE(PROD2)/=YDSPEC%NS2G) &
 & CALL ABOR1('SPECNORM_ARRAY_LEV: prod2 wrong size')
IF (SIZE(PROD3,1)/=YDSPEC%NFLEVG) &
 & CALL ABOR1('SPECNORM_ARRAY_LEV: prod3 wrong size 1')
IF (SIZE(PROD3,2)/=YDSPEC%NS3D)   &
 & CALL ABOR1('SPECNORM_ARRAY_LEV: prod3 wrong size 2')

IFLDL=YDSPEC%NFLEVL*YDSPEC%NS3D+YDSPEC%NS2D

INUMLEV(:)=YDSPEC%NUMLL(:)*YDSPEC%NS3D+YDSPEC%NS2L(:)

! -----------------------------
! Scalar product level by level
! -----------------------------

ILENS=IFLDL*YDSPEC%NUMP
ILENR(:)=IFLDL*YDSPEC%NUMPP(:)
ILENG=SUM(ILENR)
ALLOCATE(ZBUFS(ILENS))
ALLOCATE(ZBUFR(ILENG))
IF (LELAM) THEN
  IIISMAX=YDSPEC%NMSMAX
ELSE
  IIISMAX=YDSPEC%NSMAX
ENDIF
ZBUFS(:)=0.0_JPRB

! Partial products for local Ms in send buffer
! --------------------------------------------
JJ=0
IF (.NOT.LELAM) THEN
  DO JM=1,YDSPEC%NUMP
    IM=YDSPEC%MYMS(JM)
    JS=YDSPEC%NASM0(IM)-1
    ILEN=2*(YDSPEC%NSMAX-IM+1)
    DO JF=1,YDSPEC%NS3D
      DO JN=1,ILEN
        DO JL=1,YDSPEC%NFLEVL
          ZBUFS(JJ+JL) = ZBUFS(JJ+JL) &
              & + YDSPEC%SP3D(JL,JS+JN,JF) * YDSPEC%SP3D(JL,JS+JN,JF)
        ENDDO
      ENDDO
      JJ=JJ+YDSPEC%NFLEVL
    ENDDO
  ENDDO
ELSE
  DO JM=1,YDSPEC%NUMP
    IM=YDSPEC%MYMS(JM)
    JS=YDSPEC%NESM0(IM)-1
    ILEN=YDSPEC%NCPL4M(IM)
    DO JF=1,YDSPEC%NS3D
      DO JN=1,ILEN 
        DO JL=1,YDSPEC%NFLEVL
          ZBUFS(JJ+JL) = ZBUFS(JJ+JL) &
              & + YDSPEC%SP3D(JL,JS+JN,JF) * YDSPEC%SP3D(JL,JS+JN,JF)      
        ENDDO 
      ENDDO
      JJ=JJ+YDSPEC%NFLEVL
    ENDDO   
  ENDDO
ENDIF

IF (.NOT.LELAM) THEN
  DO JM=1,YDSPEC%NUMP
    IM=YDSPEC%MYMS(JM)
    JS=YDSPEC%NASM0(IM)-1
    ILEN=2*(YDSPEC%NSMAX-IM+1)
    DO JF=1,YDSPEC%NS2D
      JJ=JJ+1
      DO JN=1,ILEN
        ZBUFS(JJ) = ZBUFS(JJ) &
            & + YDSPEC%SP2D(JS+JN,JF) * YDSPEC%SP2D(JS+JN,JF)
      ENDDO
    ENDDO
  ENDDO
ELSE
  DO JM=1,YDSPEC%NUMP
    IM=YDSPEC%MYMS(JM) 
    JS=YDSPEC%NESM0(IM)-1
    ILEN=YDSPEC%NCPL4M(IM)
    DO JF=1,YDSPEC%NS2D
      JJ=JJ+1
      DO JN=1,ILEN
        ZBUFS(JJ) = ZBUFS(JJ) &
            & + YDSPEC%SP2D(JS+JN,JF) * YDSPEC%SP2D(JS+JN,JF)
      ENDDO
    ENDDO
  ENDDO 
ENDIF

IF (JJ/=ILENS) CALL ABOR1('SPECNORM_ARRAY_LEV: wrong ilens')

! Communicate to get all the Ms
! -----------------------------
CALL MPL_ALLGATHERV(PSENDBUF=ZBUFS,PRECVBUF=ZBUFR,&
                  & KRECVCOUNTS=ILENR,KCOMM=MPL_ALL_MS_COMM, &
                  & CDSTRING='SPECTRAL_FIELDS_PARA_MOD:SPECNORM_ARRAY_LEV')

ALLOCATE(ZALLMS(IFLDL,0:IIISMAX))

! Reorder Ms for reproducibility
! ------------------------------
JJ=0
DO JA=1,NPRTRW
  DO JM=1,YDSPEC%NUMPP(JA)
    IM=YDSPEC%MS_PROC(JM,JA)
    DO JL=1,YDSPEC%NS3D*YDSPEC%NFLEVL
      JJ=JJ+1
      ZALLMS(JL,IM)=ZBUFR(JJ)
    ENDDO
  ENDDO
  DO JM=1,YDSPEC%NUMPP(JA)
    IM=YDSPEC%MS_PROC(JM,JA)
    DO JF=1,YDSPEC%NS2D
      JL=YDSPEC%NS3D*YDSPEC%NFLEVL+JF
      JJ=JJ+1
      ZALLMS(JL,IM)=ZBUFR(JJ)
    ENDDO
  ENDDO
ENDDO

DEALLOCATE(ZBUFS)
DEALLOCATE(ZBUFR)

! Total scalar products for each level
! ------------------------------------
ALLOCATE(ZBUFL(IFLDL))
ZBUFL(:) = 0.0_JPRB
DO JM=1,IIISMAX  ! Start from M=1
  ZBUFL(:) = ZBUFL(:) + ZALLMS(:,JM)
ENDDO

! Apply weight for computing spectral norm
IF (.NOT.LELAM) ZBUFL(:) = 2.0_JPRB * ZBUFL(:)

! Weights do not apply to M=0
ZBUFL(:) = ZBUFL(:) + ZALLMS(:,0)

! Square root
ZBUFL(:) = SQRT( ZBUFL(:) )

DEALLOCATE(ZALLMS)
! ------------------------------------------------------------------

! Gather values for all levels
! ----------------------------
IFLDG=YDSPEC%NFLEVG*YDSPEC%NS3D+YDSPEC%NS2G
ALLOCATE(ZBUFG(IFLDG))

CALL MPL_ALLGATHERV(PSENDBUF=ZBUFL,PRECVBUF=ZBUFG,&
                  & KRECVCOUNTS=INUMLEV,KCOMM=MPL_ALL_LEVS_COMM, &
                  & CDSTRING='SPECTRAL_FIELDS_PARA_MOD:SPECNORM_ARRAY_LEV')

! Reorder levels for reproducibility
! ----------------------------------
II=0
DO JB=1,NPRTRV
  DO JF=1,YDSPEC%NS3D
    DO JL=1,YDSPEC%NUMLL(JB)
      IL=YDSPEC%NPTRLL(JB)+JL-1
      II=II+1
      PROD3(IL,JF)=ZBUFG(II)
    ENDDO
  ENDDO
  DO JF=1,YDSPEC%NS2L(JB)
    II=II+1
    PROD2(JF)=ZBUFG(II)
  ENDDO
ENDDO

! Finish
! ------
DEALLOCATE(ZBUFL)
DEALLOCATE(ZBUFG)

IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:SPECNORM_ARRAY_LEV',1,ZHOOK_HANDLE)
END SUBROUTINE SPECNORM_ARRAY_LEV

! ------------------------------------------------------------------
SUBROUTINE SPECMAX_PRT_LEV(YDSPEC,CDGREP)
! ------------------------------------------------------------------

IMPLICIT NONE
TYPE(SPECTRAL_FIELD), INTENT(IN) :: YDSPEC
CHARACTER(LEN=*), INTENT(IN) :: CDGREP

INTEGER(KIND=JPIM) :: JL
REAL(KIND=JPRB), ALLOCATABLE :: ZNOR2(:),ZNOR3(:,:)
CHARACTER(LEN=32) :: CLFMT

! ----------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:SPECMAX_PRT_LEV',0,ZHOOK_HANDLE)
ALLOCATE(ZNOR3(YDSPEC%NFLEVG,YDSPEC%NS3D))
ALLOCATE(ZNOR2(YDSPEC%NS2G))

CALL SPECMAX_ARRAY_LEV(YDSPEC,ZNOR3,ZNOR2)

IF (YDSPEC%NS2G>0) THEN
  CLFMT='(2(X,A),5(X,ES18.12))'
  WRITE(NULOUT,CLFMT)CDGREP,'SP MAX 2D:',ZNOR2
ENDIF

IF (YDSPEC%NS3D>0) THEN
  CLFMT='(2(X,A),X,I3,5(X,ES18.12))'
  DO JL=1,YDSPEC%NFLEVG
    WRITE(NULOUT,CLFMT)CDGREP,'SP MAX 3D:',JL,ZNOR3(JL,:)
  ENDDO
ENDIF

DEALLOCATE(ZNOR2)
DEALLOCATE(ZNOR3)

IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:SPECMAX_PRT_LEV',1,ZHOOK_HANDLE)
END SUBROUTINE SPECMAX_PRT_LEV

! ------------------------------------------------------------------
SUBROUTINE SPECMAX_ARRAY_LEV(YDSPEC,PROD3,PROD2)
! ------------------------------------------------------------------
!   Purpose.
!   --------
!     Actual computation of max of distributed spectral fields.
!     Returns one value per field in prod2 for 2D fields, one value
!     per fields and per level in prod3 for 3D norms.
!
!   Author.
!   -------
!     Yannick Tremolet
!
!   Modifications.
!   --------------
!     Original   02-Nov-2004
! ------------------------------------------------------------------

IMPLICIT NONE

TYPE(SPECTRAL_FIELD), INTENT(IN) :: YDSPEC
REAL(KIND=JPRB), INTENT(OUT) :: PROD2(:)
REAL(KIND=JPRB), INTENT(OUT) :: PROD3(:,:)

INTEGER(KIND=JPIM) :: JF,JJ,JL,JS,JA,JB,II,IL,IFLDL,IFLDG
INTEGER(KIND=JPIM) :: INUMLEV(NPRTRV)
REAL(KIND=JPRB), ALLOCATABLE :: ZBUFL(:),ZBUFG(:)

! ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:SPECMAX_ARRAY_LEV',0,ZHOOK_HANDLE)
IF (SIZE(PROD2)/=YDSPEC%NS2G) &
 & CALL ABOR1('SPECTRAL_FIELDS_PARA_MOD:SPECMAX_ARRAY_LEV: prod2 wrong size')
IF (SIZE(PROD3,1)/=YDSPEC%NFLEVG) &
 & CALL ABOR1('SPECTRAL_FIELDS_PARA_MOD:SPECMAX_ARRAY_LEV: prod3 wrong size 1')
IF (SIZE(PROD3,2)/=YDSPEC%NS3D)   &
 & CALL ABOR1('SPECTRAL_FIELDS_PARA_MOD:SPECMAX_ARRAY_LEV: prod3 wrong size 2')

! -----------------------------
! Max level by level
! -----------------------------
IFLDL=YDSPEC%NFLEVL*YDSPEC%NS3D+YDSPEC%NS2D
ALLOCATE(ZBUFL(IFLDL))
ZBUFL(:)=0.0_JPRB

JJ=0
DO JF=1,YDSPEC%NS3D
  DO JS=1,YDSPEC%NSPEC2
    DO JL=1,YDSPEC%NFLEVL
      ZBUFL(JJ+JL) = MAX(ZBUFL(JJ+JL), ABS(YDSPEC%SP3D(JL,JS,JF)))
    ENDDO
  ENDDO
  JJ=JJ+YDSPEC%NFLEVL
ENDDO

DO JF=1,YDSPEC%NS2D
  JJ=JJ+1
  DO JS=1,YDSPEC%NSPEC2
    ZBUFL(JJ) = MAX(ZBUFL(JJ), ABS(YDSPEC%SP2D(JS,JF)))
  ENDDO
ENDDO

! Communicate to get global max
! -----------------------------
CALL MPL_ALLREDUCE(PSENDBUF=ZBUFL,CDOPER='MAX',KCOMM=MPL_ALL_MS_COMM, &
                 & CDSTRING='SPECTRAL_FIELDS_PARA_MOD:SPECMAX_ARRAY_LEV')

! Gather values for all levels
! ----------------------------
IFLDG=YDSPEC%NFLEVG*YDSPEC%NS3D+YDSPEC%NS2G
ALLOCATE(ZBUFG(IFLDG))
INUMLEV(:)=YDSPEC%NUMLL(:)*YDSPEC%NS3D+YDSPEC%NS2L(:)

CALL MPL_ALLGATHERV(PSENDBUF=ZBUFL,PRECVBUF=ZBUFG,&
                  & KRECVCOUNTS=INUMLEV,KCOMM=MPL_ALL_LEVS_COMM, &
                  & CDSTRING='SPECTRAL_FIELDS_PARA_MOD:SPECMAX_ARRAY_LEV')

! Reorder levels
! --------------
II=0
DO JB=1,NPRTRV
  DO JF=1,YDSPEC%NS3D
    DO JL=1,YDSPEC%NUMLL(JB)
      IL=YDSPEC%NPTRLL(JB)+JL-1
      II=II+1
      PROD3(IL,JF)=ZBUFG(II)
    ENDDO
  ENDDO
  DO JF=1,YDSPEC%NS2L(JB)
    II=II+1
    PROD2(JF)=ZBUFG(II)
  ENDDO
ENDDO

! Finish
! ------
DEALLOCATE(ZBUFL)
DEALLOCATE(ZBUFG)

IF (LHOOK) CALL DR_HOOK('SPECTRAL_FIELDS_PARA_MOD:SPECMAX_ARRAY_LEV',1,ZHOOK_HANDLE)
END SUBROUTINE SPECMAX_ARRAY_LEV
! ------------------------------------------------------------------

END MODULE SPECTRAL_FIELDS_PARA_MOD
