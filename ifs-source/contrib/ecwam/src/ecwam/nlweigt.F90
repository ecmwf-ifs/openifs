! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE NLWEIGT

! ----------------------------------------------------------------------

!**** *NLWEIGT* - COMPUTATION OF INDEX ARRAYS AND WEIGHTS FOR THE
!                 COMPUTATION OF THE NONLINEAR TRANSFER RATE.

!     SUSANNE HASSELMANN JUNE 86.

!     H. GUNTHER   ECMWF/GKSS  DECEMBER 90 - CYCLE_4 MODIFICATIONS.
!                                            4 FREQUENCIES ADDED.

!*    PURPOSE.
!     --------

!       COMPUTATION OF PARAMETERS USED IN DISCRETE INTERACTION
!       PARAMETERIZATION OF NONLINEAR TRANSFER.

!**   INTERFACE.
!     ----------

!       *CALL* *NLWEIGT*

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *JAFU*      - FUNCTION FOR COMPUTATION OF ANGULAR INDICES
!                     OF K(F,THET).

!     REFERENCE.
!     ----------
!       S. HASSELMANN AND K. HASSELMANN, JPO, 1985 B.


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DELTH    ,FRATIO
      USE YOWINDN  , ONLY : IKP      ,IKP1     ,IKM      ,IKM1     ,    &
     &            K1W      ,K2W      ,K11W     ,K21W     ,AF11     ,    &
     &            FKLAP    ,FKLAP1   ,FKLAM    ,FKLAM1   ,ACL1     ,    &
     &            ACL2     ,CL11     ,CL21     ,DAL1     ,DAL2     ,    &
     &            FRH      ,KFRH     ,MFRSTLW  ,MLSTHG
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : PI       ,DEG
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "jafu.intfb.h"
!
!*    *PARAMETER*  FOR DISCRETE APPROXIMATION OF NONLINEAR TRANSFER

      REAL(KIND=JWRB), PARAMETER :: ALAMD=0.25_JWRB
      REAL(KIND=JWRB), PARAMETER :: CON=3000.0_JWRB
!
!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *ALAMD*     REAL      LAMBDA
!      *CON*       REAL      WEIGHT FOR DISCRETE APPROXIMATION OF
!                            NONLINEAR TRANSFER
! ----------------------------------------------------------------------

      INTEGER(KIND=JWIM) :: I, ISP, ISM, M, K
      INTEGER(KIND=JWIM) :: KLP1, IC, KH, KLH, KS, ISG, K1, K11, K2, K21, IKN
      INTEGER(KIND=JWIM), ALLOCATABLE :: JA1(:,:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: JA2(:,:)

      REAL(KIND=JWRB) :: F1P1, XF, COSTH3, DELPHI1, COSTH4, DELPHI2, CL1, CL2, CH
      REAL(KIND=JWRB) :: CL1H, CL2H, FRG, FLP, FLM, FKP, FKM, DELTHA, AL11, AL12 
      REAL(KIND=JWRB), ALLOCATABLE :: FRLON(:)

! ----------------------------------------------------------------------

!     0. ALLOCATE ARRAYS
!        ---------------

      F1P1 = LOG10(FRATIO)
      ISP = INT(LOG10(1.0_JWRB+ALAMD)/F1P1+.000001_JWRB)
      ISM = FLOOR(LOG10(1.0_JWRB-ALAMD)/F1P1+.0000001_JWRB)

      MFRSTLW=1+ISM

      MLSTHG=NFRE-ISM

      KFRH=-ISM+ISP+2

      ALLOCATE(JA1(NANG,2))
      ALLOCATE(JA2(NANG,2))
      ALLOCATE(FRLON(MFRSTLW:NFRE+KFRH))

      IF (.NOT.ALLOCATED(IKP)) ALLOCATE(IKP(MFRSTLW:MLSTHG))
      IF (.NOT.ALLOCATED(IKP1)) ALLOCATE(IKP1(MFRSTLW:MLSTHG))
      IF (.NOT.ALLOCATED(IKM)) ALLOCATE(IKM(MFRSTLW:MLSTHG))
      IF (.NOT.ALLOCATED(IKM1)) ALLOCATE(IKM1(MFRSTLW:MLSTHG))
      IF (.NOT.ALLOCATED(K1W)) ALLOCATE(K1W(NANG,2))
      IF (.NOT.ALLOCATED(K2W)) ALLOCATE(K2W(NANG,2))
      IF (.NOT.ALLOCATED(K11W)) ALLOCATE(K11W(NANG,2))
      IF (.NOT.ALLOCATED(K21W)) ALLOCATE(K21W(NANG,2))
      IF (.NOT.ALLOCATED(AF11)) ALLOCATE(AF11(MFRSTLW:MLSTHG))
      IF (.NOT.ALLOCATED(FKLAP)) ALLOCATE(FKLAP(MFRSTLW:MLSTHG))
      IF (.NOT.ALLOCATED(FKLAP1)) ALLOCATE(FKLAP1(MFRSTLW:MLSTHG))
      IF (.NOT.ALLOCATED(FKLAM)) ALLOCATE(FKLAM(MFRSTLW:MLSTHG))
      IF (.NOT.ALLOCATED(FKLAM1)) ALLOCATE(FKLAM1(MFRSTLW:MLSTHG))
      IF (.NOT.ALLOCATED(FRH)) ALLOCATE(FRH(KFRH))


!*    1. COMPUTATION FOR ANGULAR GRID.
!        -----------------------------
!
!*    1.1 DETERMINE ANGLES DELPHI USING RESONANCE CONDITION.
!         --------------------------------------------------
!     
      XF      = ((1.0_JWRB+ALAMD)/(1.0_JWRB-ALAMD))**4
      COSTH3  = (1.0_JWRB+2.0_JWRB*ALAMD+2.0_JWRB*ALAMD**3)/(1.0_JWRB+ALAMD)**2
      DELPHI1 = -180.0_JWRB/PI*ACOS(COSTH3)
      COSTH4  = SQRT(1.0_JWRB-XF+XF*COSTH3**2)
      DELPHI2 = 180.0_JWRB/PI*ACOS(COSTH4)

      DELTHA = DELTH*DEG
      CL1 = DELPHI1/DELTHA
      CL2 = DELPHI2/DELTHA

!*    1.1 COMPUTATION OF INDICES OF ANGULAR CELL.
!         ---------------------------------------

      KLP1 = NANG+1
      IC = 1
      DO KH=1,2
        KLH = NANG 
        IF (KH == 2) KLH=KLP1
        DO K=1,KLH
          KS = K
          IF (KH > 1) KS=KLP1-K+1
          IF (KS > NANG) GO TO 1002
          CH = IC*CL1
          JA1(KS,KH) = JAFU(CH,K,KLP1)
          CH = IC*CL2
          JA2(KS,KH) = JAFU(CH,K,KLP1)
 1002     CONTINUE
        ENDDO
        IC = -1
      ENDDO

!*    1.2 COMPUTATION OF ANGULAR WEIGHTS.
!         -------------------------------

      CL1  = CL1-INT(CL1)
      CL2  = CL2-INT(CL2)
      ACL1 = ABS(CL1)
      ACL2 = ABS(CL2)
      CL11 = 1.0_JWRB-ACL1
      CL21 = 1.0_JWRB-ACL2
      AL11 = (1.0_JWRB+ALAMD)**4
      AL12 = (1.0_JWRB-ALAMD)**4
      DAL1 = 1.0_JWRB/AL11
      DAL2 = 1.0_JWRB/AL12

!*    1.3 COMPUTATION OF ANGULAR INDICES.
!         -------------------------------

      ISG = 1
      DO KH=1,2
        CL1H = ISG*CL1
        CL2H = ISG*CL2
        DO K=1,NANG
          KS = K
          IF (KH == 2) KS = NANG-K+2
          IF (K == 1) KS = 1
          K1 = JA1(K,KH)
          K1W(KS,KH) = K1
          IF (CL1H < 0.0_JWRB) THEN
            K11 = K1-1
            IF (K11 < 1) K11 = NANG 
          ELSE
            K11 = K1+1
            IF (K11 > NANG) K11 = 1
          ENDIF
          K11W(KS,KH) = K11
          K2 = JA2(K,KH)
          K2W(KS,KH) = K2
          IF (CL2H < 0) THEN
            K21 = K2-1
            IF (K21 < 1) K21 = NANG 
          ELSE
            K21 = K2+1
            IF (K21 > NANG) K21 = 1
          ENDIF
          K21W(KS,KH) = K21
        ENDDO
        ISG = -1
      ENDDO

!*    2. COMPUTATION FOR FREQUENCY GRID.
!        -------------------------------

      DO M=1,NFRE
        FRLON(M) = FR(M)
      ENDDO
      DO M=0,MFRSTLW,-1
        FRLON(M)=FRLON(M+1)/FRATIO
      ENDDO
      DO M=NFRE+1,NFRE+KFRH
        FRLON(M) = FRATIO*FRLON(M-1)
      ENDDO
      DO M=MFRSTLW,MLSTHG
        FRG = FRLON(M)
        AF11(M) = CON * FRG**11
        FLP = FRG*(1.0_JWRB+ALAMD)
        FLM = FRG*(1.0_JWRB-ALAMD)
        IKN = M+ISP
        IKP(M) = IKN
        FKP = FRLON(IKP(M))
        IKP1(M) = IKP(M)+1
        FKLAP(M) = (FLP-FKP)/(FRLON(IKP1(M))-FKP)
        FKLAP1(M) = 1.0_JWRB-FKLAP(M)

        IKN = M+ISM
        IF (IKN >= MFRSTLW) THEN
          IKM(M) = IKN
          FKM = FRLON(IKM(M))
          IKM1(M) = IKM(M)+1
          FKLAM(M) = (FLM-FKM)/(FRLON(IKM1(M))-FKM)
          FKLAM1(M) = 1.0_JWRB-FKLAM(M)
        ELSEIF (IKN+1 == MFRSTLW) THEN
          IKM(M) = 1 
          IKM1(M) = MFRSTLW 
          FKM = FRLON(IKM1(M))/FRATIO
          FKLAM(M) = (FLM-FKM)/(FRLON(IKM1(M))-FKM)
          FKLAM1(M) = 0.0_JWRB
        ELSE
          IKM(M) = 1
          FKLAM(M) = 0.0_JWRB
          IKM1(M) = 1
          FKLAM1(M) = 0.0_JWRB
        ENDIF

      ENDDO

!*    3. COMPUTE TAIL FREQUENCY RATIOS.
!        ------------------------------

      DO I=1,KFRH
        M = NFRE+I-1
        FRH(I) = (FRLON(NFRE)/FRLON(M))**5
      ENDDO

!*    4. PRINTER PROTOCOL.
!        -----------------

      WRITE(IU06,'(1H1,'' NON LINEAR INTERACTION PARAMETERS:'')')
      WRITE(IU06,'(1H0,'' COMMON INDNL: CONSTANTS'')')
      WRITE(IU06,*)'    ALAMD = ', ALAMD
      WRITE(IU06,*)'      CON = ', CON
      WRITE(IU06,*)'  DELPHI1 = ',DELPHI1
      WRITE(IU06,*)'  DELPHI2 = ',DELPHI2
      WRITE(IU06,'(1X,''    ACL1       ACL2   '',                       &
     &             ''    CL11       CL21   '',                          &
     &             ''    DAL1       DAL2'')')
      WRITE(IU06,'(1X,6F11.8)') ACL1, ACL2, CL11, CL21, DAL1, DAL2

      WRITE(IU06,'(1H0,'' COMMON INDNL: FREQUENCY ARRAYS'')')
      WRITE(IU06,'(1X,'' M   IKP IKP1  IKM IKM1'',                      &
     &          ''   FKLAP       FKLAP1 '',                             &
     &          ''   FKLAM       FKLAM1     AF11'')')
      DO M=MFRSTLW,MLSTHG
        WRITE(IU06,'(1X,I2,4I5,4F11.8,E11.3)')                          &
     &   M, IKP(M), IKP1(M), IKM(M), IKM1(M),                           &
     &   FKLAP(M), FKLAP1(M), FKLAM(M), FKLAM1(M), AF11(M)
      ENDDO

      WRITE(IU06,'(1H0,'' COMMON INDNL: ANGULAR ARRAYS'')')
      WRITE(IU06,'(1X,''  |--------KH = 1----------|'',                 &
     &              ''|--------KH = 2----------|'')')
      WRITE(IU06,'(1X,'' K   K1W   K2W  K11W  K21W'',                   &
     &              ''   K1W   K2W  K11W  K21W'')')
      DO K=1,NANG
        WRITE(IU06,'(1X,I2,8I6)') K,(K1W(K,KH), K2W(K,KH), K11W(K,KH),  &
     &   K21W(K,KH),KH=1,2)
      ENDDO
      WRITE(IU06,'(1H0,'' COMMON INDNL: TAIL ARRAY FRH'')')
      WRITE(IU06,'(1X,8F10.7)') (FRH(M),M=1,KFRH)

!     5. DEALLOCATE LOCAL ARRAYS
!        -----------------------

      DEALLOCATE(JA1)
      DEALLOCATE(JA2)
      DEALLOCATE(FRLON)

      END SUBROUTINE NLWEIGT
