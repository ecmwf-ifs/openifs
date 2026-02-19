! (C) Copyright 2015- ECMWF.
! (C) Copyright 2015- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE INTERPOL_DECOMP_MOD

! Compute Interpolative Decomposions (ID)

! See Cheng,H., Gimbutas,Z., Martinsson,P.G. and Rokhlin,V. (2005) 
! "On the compression of low rank matrices", SIAM.J.Sci.Comput.,
!  Vol. 26, No. 4, pp1389-1404

! Also lecture notes "Mulilevel compression of Linear Operators:
! Descendents of Fast Multiple Methods and Calderon-Zygmund Theory"
! P.G.Martinsson and Mark Tygert, 2011. Chapter 7.

! Author: Mats Hamrud


USE PARKIND1, ONLY : JPRB, JPIM, JPRD, JPIB
IMPLICIT NONE
CONTAINS
!===========================================================================
SUBROUTINE COMPUTE_ID(PEPS,KM,KN,PMAT,KRANK,KBCLIST,PNONIM)
IMPLICIT NONE

! Compute ID

REAL(KIND=JPRD),INTENT(IN)     :: PEPS  ! Precision for computation
                                        ! of numerical rank
INTEGER(KIND=JPIM),INTENT(IN)  :: KM    ! Number of rows in matrix pmat
INTEGER(KIND=JPIM),INTENT(IN)  :: KN    ! Number of columns in matrix pmat
REAL(KIND=JPRD)   ,INTENT(IN)  :: PMAT(:,:)  ! Original matrix
INTEGER(KIND=JPIM),INTENT(OUT) :: KRANK      ! Numerical rank
INTEGER(KIND=JPIM),INTENT(OUT) :: KBCLIST(:) ! List of  columns
REAL(KIND=JPRD)   ,INTENT(OUT) :: PNONIM(:,:)  ! Non-identity part of projection
                                               ! matrix

INTEGER(KIND=JPIM) :: JM,JN
REAL(KIND=JPRD) :: ZR(KM,KN)
REAL(KIND=JPRD),ALLOCATABLE :: ZS(:,:),ZT(:,:)
!----------------------------------------------------------------------------
!Avoid destroying input matrix
ZR(:,:) = PMAT(1:KM,1:KN)
! Householder QR
CALL ALG541(PEPS,KM,KN,ZR,KRANK,KBCLIST)

DO JN=1,KN
  DO JM=JN+1,KM
    ZR(JM,JN) = 0.0_JPRD
  ENDDO
ENDDO
! S leftmost kxk block of R
ALLOCATE(ZS(KRANK,KRANK))
DO JN=1,KRANK
  DO JM=1,KRANK
    IF(JM <= KM ) THEN
      ZS(JM,JN) = ZR(JM,JN) 
    ELSE
      ZS(JM,JN) = 0.0_JPRD
    ENDIF
  ENDDO
ENDDO
! T Rightmost kx(k-n) block of R 
ALLOCATE(ZT(KRANK,KN-KRANK))
DO JN=1,KN-KRANK
  DO JM=1,KRANK
    IF(JM <= KM ) THEN
      ZT(JM,JN) = ZR(JM,JN+KRANK) 
    ELSE
      ZT(JM,JN) = 0.0_JPRD
    ENDIF
  ENDDO
ENDDO
!Solve linear equation (BLAS level 3 routine)
IF( KRANK <= 0 ) THEN
  write(0,*) 'warning: KRANK DTRSM ', KRANK, KM, KN
  CALL ABOR1('DTRSM : KRANK <=0 not allowed')
ENDIF

!  IF (JPRB == JPRD) THEN
CALL DTRSM('Left','Upper','No transpose','Non-unit',KRANK,KN-KRANK,1.0_JPRD, &
 & ZS,KRANK,ZT,KRANK)
!  ELSE
!    CALL STRSM('Left','Upper','No transpose','Non-unit',KRANK,KN-KRANK,1.0_JPRD, &
!     & ZS,KRANK,ZT,KRANK)
!  ENDIF

DO JM=1,KRANK
  DO JN=1,KN-KRANK
    PNONIM(JM,JN) = ZT(JM,JN)
  ENDDO
ENDDO
DEALLOCATE(ZS,ZT)
!!$IF(KRANK < KN) THEN
!!$  PRINT *,'MAXVAL PNONIM ',KM,KM,KRANK,MAXVAL( PNONIM(1:KRANK,1:KN-KRANK))
!!$ENDIF
END SUBROUTINE COMPUTE_ID

!==============================================================================
SUBROUTINE ALG541(PEPS,KM,KN,PA,KRANK,KLIST)
IMPLICIT NONE

! Householder QR with Column Pivoting
! Algorithm 5.4.1 from Matrix Computations, G.H.Golub & C.F van Loen, third ed.

! Algorithm modified to terminate at numerical precision "peps"

REAL(KIND=JPRD),INTENT(IN)     :: PEPS     ! Precision
INTEGER(KIND=JPIM),INTENT(IN)  :: KM       ! Number of rows in matrix pa
INTEGER(KIND=JPIM),INTENT(IN)  :: KN       ! Number of columns in matrix pa
REAL(KIND=JPRD),INTENT(INOUT)  :: PA(:,:)  ! On input : original matrix
                                           ! on output : R in upper triangle etc
                                           ! see Golub&Van Loen
INTEGER(KIND=JPIM),INTENT(OUT) :: KRANK    ! Numerical rank of matrix
INTEGER(KIND=JPIM),INTENT(OUT) :: KLIST(:) ! List of columns (pivots)

INTEGER(KIND=JPIM)           :: JM,JN,ISWAP,JK,IK,IIK,IM,IN,IMIN,ILIST(KN)
REAL(KIND=JPRD) :: ZC(KN),ZTAU,ZSWAPA(KM),ZSWAP,ZV(KM),ZBETA,ZWORK(KN),ZTAU_IN
REAL(KIND=JPRD) :: ZTAU_REC,ZEPS
!-------------------------------------------------------------------------------
ZEPS = 10000.0_JPRD*EPSILON(ZEPS)
IMIN=MIN(KM,KN)
! Compute initial column norms,its max and the first column where c=tau
IK = 0
ZTAU = 0._JPRD
DO JN=1,KN
  ZC(JN) = DOT_PRODUCT(PA(1:KM,JN),PA(1:KM,JN))
  IF(ZC(JN) > ZTAU) THEN
    IK = JN
    ZTAU = ZC(JN)
  ENDIF
ENDDO
ZTAU_IN = ZTAU
ZTAU_REC= ZTAU
KRANK = 0
DO WHILE (ZTAU > PEPS**2*ZTAU_IN)
  KRANK = KRANK+1
  IF( KRANK <= IMIN ) THEN
    ILIST(KRANK) = IK
    ! Column swap KRANK with IK 
    ZSWAPA(:) = PA(:,KRANK)
    PA(:,KRANK) = PA(:,IK)
    PA(:,IK) = ZSWAPA(:)
    ZSWAP = ZC(KRANK)
    ZC(KRANK) = ZC(IK)
    ZC(IK) = ZSWAP
    ! Compute Householder vector
    ZBETA=0.0_JPRD
    IF( KM-KRANK >= 0 ) THEN
      CALL ALG511(ZEPS,KM-KRANK+1,PA(KRANK:KM,KRANK),ZV,ZBETA)
    ENDIF
    ! Apply Householder matrix
    IM = KM-KRANK+1
    IN = KN-KRANK+1
    ! LAPACK
    CALL DLARF('Left',IM,IN,ZV,1,ZBETA,PA(KRANK,KRANK),KM,ZWORK)
  ENDIF

! Update column norms
  ZTAU = 0.0_JPRD
  IF(KRANK < IMIN) THEN
    PA(KRANK+1:KM,KRANK) = ZV(2:IM)
    DO JN=KRANK+1,KN
      ZC(JN) = ZC(JN)-PA(KRANK,JN)**2
      IF(ZC(JN) > ZTAU) THEN
        IK = JN
        ZTAU = ZC(JN)
      ENDIF
    ENDDO
! Re-compute column norms due to round-off error
    IF(ZTAU < ZEPS*ZTAU_REC .OR. ZTAU < 0._JPRD .or. (KN-KRANK) > 100 ) THEN
      DO JN=KRANK+1,KN
        ZC(JN) = DOT_PRODUCT(PA(KRANK+1:,JN),PA(KRANK+1:,JN))
        IF(ZC(JN) > ZTAU) THEN
          IK = JN
          ZTAU = ZC(JN)
        ENDIF
      ENDDO
      !write(0,*) 'RECOMPUTE TAU ',KRANK,ZTAU_REC,ZTAU
      ZTAU_REC = ZTAU
    ENDIF
  ENDIF
ENDDO
! Make sure klist is filled also beyond krank
DO JN=1,KN
  KLIST(JN) = JN
ENDDO
DO JN=1,KRANK
  ISWAP = KLIST(JN)
  KLIST(JN) = KLIST(ILIST(JN))
  KLIST(ILIST(JN)) = ISWAP
ENDDO
  
END SUBROUTINE ALG541

!==============================================================================
SUBROUTINE ALG511(PEPS,KSIZE,PX,PV,PBETA)
IMPLICIT NONE
! Compute Householder vector
! Algorithm 5.1.1 from Matrix Computations, G.H.Golub & C.F van Loen, third ed.
REAL(KIND=JPRD),INTENT(IN)     :: PEPS     ! Precision
REAL(KIND=JPRD),INTENT(IN)     :: PX(:)
INTEGER(KIND=JPIM), INTENT(IN)     :: KSIZE
REAL(KIND=JPRD),INTENT(OUT)    :: PV(:)
REAL(KIND=JPRD),INTENT(OUT)    :: PBETA
INTEGER(KIND=JPIM) :: IL
REAL(KIND=JPRD) :: ZSIGMA,ZMU, ZNORM
REAL(KIND=JPRD) :: ZX(KSIZE)
!-------------------------------------------------------------------------------
! normalize
ZNORM=0._JPRD
DO IL=1,KSIZE
  ZNORM = ZNORM + PX(IL)*PX(IL)
ENDDO
ZNORM=SQRT(ZNORM)
ZX(:)=PX(1:KSIZE)
IF( ZNORM > PEPS ) ZX(:)=PX(1:KSIZE)/ZNORM

ZSIGMA=0._JPRD
IF( KSIZE > 1 ) ZSIGMA = DOT_PRODUCT(ZX(2:KSIZE),ZX(2:KSIZE))
PV(1) = 1.0_JPRD
IF( KSIZE > 1 ) PV(2:KSIZE) = ZX(2:KSIZE)
IF(ABS(ZSIGMA) <  PEPS**2) THEN
  PBETA = 0.0_JPRD
ELSE
  ZMU = SQRT(ZX(1)**2+ZSIGMA)
  IF(ZX(1) <= 0.0_JPRD) THEN
    PV(1) = ZX(1)-ZMU
  ELSE
    PV(1) = -ZSIGMA/(ZX(1)+ZMU)
  ENDIF
  PBETA = 2.0_JPRD*PV(1)**2/(ZSIGMA+PV(1)**2)
  PV(:) = PV(:)/(PV(1))
ENDIF

END SUBROUTINE ALG511
!================================================================================

END MODULE INTERPOL_DECOMP_MOD
