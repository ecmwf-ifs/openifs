! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE INISNONLIN

! ----------------------------------------------------------------------

!**** *INISNONLIN* - INITIALISE ALL FREQUENCY DEPENDENT ARRAYS USED BY
!                    SNONLIN

!     J. BIDLOT   ECMWF  MAY 2012

!*    PURPOSE.
!     --------

!       USED TO BE IN SNONLIN BUT NOW IT IS ONLY COMPUTED ONCE. 

!**   INTERFACE.
!     ----------

!       *CALL* *INISNONLIN*

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FRATIO
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWINDN  , ONLY : IKP      ,IKP1     ,IKM      ,IKM1     ,    &
     &            K1W      ,K2W      ,K11W     ,K21W     ,AF11     ,    &
     &            FKLAP    ,FKLAP1   ,FKLAM    ,FKLAM1   ,ACL1     ,    &
     &            ACL2     ,CL11     ,CL21     ,DAL1     ,DAL2     ,    &
     &            FRH      ,FTRF     ,MFRSTLW  ,MLSTHG   ,              &
     &            KFRH     ,NINL     ,NRNL     ,INLCOEF  ,RNLCOEF
      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
include "nlweigt.intfb.h"

      INTEGER(KIND=JWIM) :: ICOUNT, IRCOUNT
      INTEGER(KIND=JWIM) :: MC, MP, MP1, MM, MM1, IC, IP, IP1, IM , IM1, ITEMP

      REAL(KIND=JWRB) :: ALPH, FRR
      REAL(KIND=JWRB) :: FFACP, FFACP1, FFACM, FFACM1, FTAIL, FKLAMP, FKLAMP1
      REAL(KIND=JWRB) :: FKLAMPA, FKLAMPB, FKLAMP2, FKLAPA2, FKLAPB2
      REAL(KIND=JWRB) :: FKLAP12, FKLAP22, FKLAMM, FKLAMM1, FKLAMMA, FKLAMMB 
      REAL(KIND=JWRB) :: FKLAMM2, FKLAMA2, FKLAMB2, FKLAM12, FKLAM22
      REAL(KIND=JWRB) :: GW1, GW2, GW3, GW4, GW5, GW6, GW7, GW8

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     INLINE FUNCTION (PIERSON-MOSKOWITZ SMOOTH CUT-OFF)
!     X == FR(1)/FREQUENCY
      REAL(KIND=JWRB) :: EPMMA, X
      EPMMA(X)= EXP(-MIN(1.25_JWRB*X**4,50.0_JWRB))*(X**5)

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INISNONLIN',0,ZHOOK_HANDLE)

!*    MODULE INDNL (WEIGHT OF NON-LINEAR INTERACTION).
!     ------------------------------------------------

      CALL NLWEIGT


!     1. FRONT SPECTRAL TAIL REDUCTION COEFFICIENTS

      IF (.NOT.ALLOCATED(FTRF)) ALLOCATE(FTRF(MFRSTLW:1))
      ALPH=1.0_JWRB/EPMMA(1.0_JWRB)
      FRR=1.0_JWRB
      DO MC=1,MFRSTLW,-1
         FTRF(MC)=ALPH*EPMMA(FRR)
         FRR=FRR*FRATIO
      ENDDO

!     2. WORK ARRAYS STORING THE DIFFERENT INDICES AND COEFFICIENTS

      IF (.NOT.ALLOCATED(INLCOEF)) ALLOCATE(INLCOEF(NINL,1:MLSTHG))
      IF (.NOT.ALLOCATED(RNLCOEF)) ALLOCATE(RNLCOEF(NRNL,1:MLSTHG))

!*    3. FREQUENCY LOOP.
!        ---------------
      DO MC=1,MLSTHG
        MP  = IKP (MC)
        MP1 = IKP1(MC)
        MM  = IKM (MC)
        MM1 = IKM1(MC)
        FFACP  = 1.0_JWRB
        FFACP1 = 1.0_JWRB
        FFACM  = 1.0_JWRB
        FFACM1 = 1.0_JWRB
        FTAIL  = 1.0_JWRB
        IC  = MC
!       front tail protection (keep IC >=1)
        IF (IC < 1 ) IC = 1
        IP  = MP
        IP1 = MP1
        IM  = MM
        IM1 = MM1
!       LOW FREQUENCY FRONT TAIL
        IF (IP < 1) THEN
          FFACP = FTRF(IP)
          IP=1
        ENDIF
        IF (IP1 < 1) THEN
          FFACP1 = FTRF(IP1)
          IP1=1
        ENDIF
        IF (IM < MFRSTLW) THEN
          FFACM = 0.0_JWRB
          IM = 1
        ELSE IF (IM < 1) THEN
          FFACM = FTRF(IM)
          IM = 1
        ENDIF
        IF (IM1 < MFRSTLW) THEN
          FFACM1 = 0.0_JWRB
          IM1 = 1
        ELSE IF (IM1 < 1) THEN
          FFACM1 = FTRF(IM1)
          IM1 = 1
        ENDIF

!       HIGH FREQUENCY TAIL
        IF (IP1 > NFRE) THEN
! Quick fix from Deborah
          ITEMP=IP1-NFRE+1
          IF (ITEMP > SIZE(FRH))THEN
            ITEMP=SIZE(FRH)
          ENDIF
!         FFACP1 = FRH(IP1-NFRE+1)
          FFACP1 = FRH(ITEMP)

          IP1 = NFRE
          IF (IP > NFRE) THEN
            FFACP  = FRH(IP -NFRE+1)
            IP  = NFRE
            IF (IC > NFRE) THEN
              FTAIL  = FRH(IC -NFRE+1)
              IC  = NFRE
              IF (IM1 > NFRE) THEN
                FFACM1 = FRH(IM1-NFRE+1)
                IM1 = NFRE
              ENDIF
            ENDIF
          ENDIF
        ENDIF

        ICOUNT=1
        INLCOEF(ICOUNT,MC) = IC
        ICOUNT=ICOUNT+1
        INLCOEF(ICOUNT,MC) = IP
        ICOUNT=ICOUNT+1
        INLCOEF(ICOUNT,MC) = IP1
        ICOUNT=ICOUNT+1
        INLCOEF(ICOUNT,MC) = IM
        ICOUNT=ICOUNT+1
        INLCOEF(ICOUNT,MC) = IM1

        FKLAMP  = FKLAP(MC)
        FKLAMP1 = FKLAP1(MC)
        GW2 = FKLAMP1*FFACP*DAL1
        GW1 = GW2*CL11
        GW2 = GW2*ACL1
        GW4 = FKLAMP*FFACP1*DAL1
        GW3 = GW4*CL11
        GW4 = GW4*ACL1
        FKLAMPA = FKLAMP*CL11
        FKLAMPB = FKLAMP*ACL1
        FKLAMP2 = FKLAMP1*ACL1
        FKLAMP1 = FKLAMP1*CL11
        FKLAPA2 = FKLAMPA**2
        FKLAPB2 = FKLAMPB**2
        FKLAP12 = FKLAMP1**2
        FKLAP22 = FKLAMP2**2
        IRCOUNT=1
        RNLCOEF(IRCOUNT,MC) = FTAIL
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = GW1
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = GW2
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = GW3
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = GW4
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMPA 
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMPB
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMP2
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMP1
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAPA2
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAPB2
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAP12
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAP22

        FKLAMM  = FKLAM(MC)
        FKLAMM1 = FKLAM1(MC)
        GW6 = FKLAMM1*FFACM*DAL2
        GW5 = GW6*CL21
        GW6 = GW6*ACL2
        GW8 = FKLAMM*FFACM1*DAL2
        GW7 = GW8*CL21
        GW8 = GW8*ACL2
        FKLAMMA = FKLAMM*CL21
        FKLAMMB = FKLAMM*ACL2
        FKLAMM2 = FKLAMM1*ACL2
        FKLAMM1 = FKLAMM1*CL21
        FKLAMA2 = FKLAMMA**2
        FKLAMB2 = FKLAMMB**2
        FKLAM12 = FKLAMM1**2
        FKLAM22 = FKLAMM2**2
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = GW5 
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = GW6 
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = GW7 
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = GW8 
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMMA
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMMB
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMM2
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMM1
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMA2
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMB2
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAM12
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAM22

      ENDDO

      IF (ICOUNT /= NINL) THEN
        WRITE(IU06,*) '*************************************'
        WRITE(IU06,*) 'ERROR IN INISNONLIN : ICOUNT NE NINL'
        WRITE(IU06,*) 'ICOUNT= ',ICOUNT
        WRITE(IU06,*) 'NINL= ',NINL
        WRITE(IU06,*) '*************************************'
        CALL ABORT1
      ENDIF
      IF (IRCOUNT /= NRNL) THEN
        WRITE(IU06,*) '*************************************'
        WRITE(IU06,*) 'ERROR IN INISNONLIN : IRCOUNT NE NRNL'
        WRITE(IU06,*) 'IRCOUNT= ',IRCOUNT
        WRITE(IU06,*) 'NRNL= ',NRNL
        WRITE(IU06,*) '*************************************'
        CALL ABORT1
      ENDIF

      IF (LHOOK) CALL DR_HOOK('INISNONLIN',1,ZHOOK_HANDLE)

      END SUBROUTINE INISNONLIN
