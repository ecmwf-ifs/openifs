! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE MPBCASTGRID(IU06, ISEND, ITAG)

! ----------------------------------------------------------------------
!**** *MPBCASTGRID* - BROADCAST THE CONTENT OF THE GRID FILE READ BY
!                     READPRE TO THE OTHER PE'S. IT WILL ALSO ALLOCATE
!                     THE NECESSARY ARRAYS ON THE RECEIVING PE'S AS THEY
!                     WERE NOT ON THOSE PE'S SINCE READPRE WAS NOT
!                     CALLED FOR THEM.

!     J. BIDLOT    ECMWF   OCTOBER 1997

!     PURPOSE.
!     --------
!     BROADCAST THE CONTENT OF THE GRID FILE READ BY READPRE
!     TO THE OTHER PE'S
!*    INTERFACE.
!     ----------

!     CALL *MPBCASTGRID*(IU06,ISEND,ITAG)

!     *IU06*      UNIT FOR PRINTER MESSAGES.
!     *ISEND*     RANK OF THE PROCESS ONTO WHICH FIELD IS COLLECTED
!     *ITAG*      TAG ASSOCIATED WITH AS A PARTICULAR CALL TO SUBROUTINE
!                 THIS IS NECESSARY TO DIFFERENTIATE THE DIFFERENT CALLS

!     METHOD.
!     -------
!     MPL_BROADCAST FROM PROCESSOR CORRESPONDING TO ISEND TO
!     ALL OTHER PROCESSORS.

!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_BROADCAST

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : NGOUT    ,IJAR
      USE YOWFRED  , ONLY : FR       ,DFIM     ,GOM      ,C        ,    &
     &            DELTH    ,DELTR    ,TH       ,COSTH    ,SINTH
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,SINPH    ,COSPH    ,    &
     &            IJS      ,IJL
      USE YOWINDN  , ONLY : KFRH     ,IKP      ,IKP1     ,IKM      ,    &
     &            IKM1     ,K1W      ,K2W      ,K11W     ,K21W     ,    &
     &            AF11     ,FKLAP    ,FKLAP1   ,FKLAM    ,FKLAM1   ,    &
     &            ACL1     ,ACL2     ,CL11     ,CL21     ,DAL1     ,    &
     &            DAL2     ,FRH      ,MFRSTLW  ,MLSTHG
      USE YOWMAP   , ONLY : BLK2GLO  ,NX       ,NY       ,              &
     &            IPER     ,IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,    &
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,NLONRGG  ,    &
     &            IQGAUSS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPRECR   ,NPRECI
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED ,              &
     &            NGX      ,NGY      ,                                  &
     &            NIBLO    ,NOVER    ,NIBL1    ,CLDOMAIN
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTH_INPUT,DEPTHA   ,DEPTHD   ,  &
     &            TCGOND   ,TFAK     ,TSIHKD   ,TFAC_ST
      USE YOWTABL  , ONLY : FAC0     ,FAC1     ,FAC2     ,FAC3     ,    &
     &            FAK      ,FRHF     ,DFIMHF   ,NFREHF   ,              &
     &            MR       ,XMR      ,MA       ,XMA      ,NFREH    ,    &
     &            NANGH    ,NMAX     ,OMEGA    ,DFDTH    ,THH      ,    &
     &            DELTHH   ,IM_P     ,IM_M     ,TA       ,TB       ,    &
     &            TC_QL    ,TT_4M    ,TT_4P    ,TFAKH

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_BROADCAST

!----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU06, ISEND
      INTEGER(KIND=JWIM), INTENT(INOUT) :: ITAG
      INTEGER(KIND=JWIM), PARAMETER :: MFIRST=19
      INTEGER(KIND=JWIM) :: I, J, IJ, K, K1, K2, M, M1, M2, IC, L,      &
     &                      KDEPTH, NGOU
      INTEGER(KIND=JWIM) :: IKCOUNT, KCOUNT
      INTEGER(KIND=JWIM) :: MIC, MZC 
      INTEGER(KIND=JWIM),ALLOCATABLE :: ICOMBUF(:)

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),ALLOCATABLE :: ZCOMBUF(:)

!----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MPBCASTGRID',0,ZHOOK_HANDLE)

      IF (ISEND == 0 .OR. NPROC == 1) THEN
         WRITE (IU06,*) ''
!     1.1 SEND TO ALL PROCESSORS OTHER THAN ISEND
!         ------------------------------------------------
      ELSE

!       BUFFER SIZE MESSAGE AND THE FEW DIMENSIONS NEEDED TO 
!       ALLOCATE ALL ARRAYS
        ALLOCATE(ICOMBUF(MFIRST))

        IF (IRANK == ISEND) THEN
          IKCOUNT=0
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NANG
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NFRE
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NFRE_RED
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NGX
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NGY
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NIBLO
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NOVER
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=KFRH
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=MFRSTLW
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=MLSTHG
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NIBL1
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=ICHAR(CLDOMAIN)
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NGOUT
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NDEPTH
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=MR
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=MA
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NFREH
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NANGH
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NMAX
          IF (IKCOUNT /= MFIRST) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* IKCOUNT .NE. MFIRST !!!*' 
            WRITE (IU06,*) '* ON IRANK = ',IRANK
            WRITE (IU06,*) '* IKCOUNT = ',IKCOUNT
            WRITE (IU06,*) '* MFIRST  = ',MFIRST
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF
        ENDIF

        CALL MPL_BROADCAST(ICOMBUF,KROOT=ISEND,KTAG=ITAG,CDSTRING='MPBCASTGRID:')
        ITAG=ITAG+1

        IF (IRANK /= ISEND) THEN
          IKCOUNT=0
          IKCOUNT=IKCOUNT+1
          NANG=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NFRE=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NFRE_RED=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NGX=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NGY=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NIBLO=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NOVER=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          KFRH=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          MFRSTLW=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          MLSTHG=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NIBL1=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          CLDOMAIN=CHAR(ICOMBUF(IKCOUNT))
          IKCOUNT=IKCOUNT+1
          NGOUT=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NDEPTH=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          MR=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          MA=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NFREH=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NANGH=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NMAX=ICOMBUF(IKCOUNT)
          IF (IKCOUNT /= MFIRST) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* IKCOUNT .NE. MFIRST !!!*' 
            WRITE (IU06,*) '* ON IRANK = ',IRANK
            WRITE (IU06,*) '* IKCOUNT = ',IKCOUNT
            WRITE (IU06,*) '* MFIRST  - ',MFIRST
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF
        ENDIF
        DEALLOCATE(ICOMBUF)

        MIC=7+NGY+2*NIBLO+4*(MLSTHG-MFRSTLW+1)+             &
     &      8*NANG+NGOUT+2*NFREH*NFREH
        MZC=17+(4+4*NDEPTH)*NFRE+5*(MLSTHG-MFRSTLW+1)+3*NANG+4*NGY+     &
     &      KFRH+                                                       &
     &      NIBLO+4*NANG*NANG*NFREHF*NFREHF+3*NFREHF+              &
     &      2+2*NFREH+NANGH+NFREH*NDEPTH+5*NANGH*NDEPTH*NFREH*NFREH

!       ENCODE MAIN MESSAGE BUFFERS (ON PE=ISEND) AND
!       ALLOCATE ALL ARRAYS NEEDED TO KEEP THE BUFFERS ON THE OTHER PE'S

        ALLOCATE(ICOMBUF(MIC))
        ALLOCATE(ZCOMBUF(MZC))

        IF (IRANK /= ISEND) THEN

          IF (.NOT.ALLOCATED(FR)) ALLOCATE(FR(NFRE))
          IF (.NOT.ALLOCATED(DFIM)) ALLOCATE(DFIM(NFRE))
          IF (.NOT.ALLOCATED(GOM)) ALLOCATE(GOM(NFRE))
          IF (.NOT.ALLOCATED(C)) ALLOCATE(C(NFRE))
          IF (.NOT.ALLOCATED(TH)) ALLOCATE(TH(NANG))
          IF (.NOT.ALLOCATED(COSTH)) ALLOCATE(COSTH(NANG))
          IF (.NOT.ALLOCATED(SINTH)) ALLOCATE(SINTH(NANG))
          IF (.NOT.ALLOCATED(DELLAM)) ALLOCATE(DELLAM(NGY))
          IF (.NOT.ALLOCATED(NLONRGG)) ALLOCATE(NLONRGG(NGY))
          IF (.NOT.ALLOCATED(SINPH)) ALLOCATE(SINPH(NGY))
          IF (.NOT.ALLOCATED(COSPH)) ALLOCATE(COSPH(NGY))
          IF (.NOT.ALLOCATED(BLK2GLO)) ALLOCATE(BLK2GLO(NIBLO))
          IF (.NOT.ALLOCATED(ZDELLO)) ALLOCATE(ZDELLO(NGY))
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

          IF (NGOUT > 0) THEN
            IF (.NOT.ALLOCATED(IJAR)) ALLOCATE(IJAR(NGOUT))
          ENDIF

          IF (ALLOCATED(DEPTH_INPUT)) DEALLOCATE(DEPTH_INPUT)
          ALLOCATE(DEPTH_INPUT(NIBLO))

          IF (.NOT.ALLOCATED(TCGOND)) ALLOCATE(TCGOND(NDEPTH,NFRE))
          IF (.NOT.ALLOCATED(TFAK)) ALLOCATE(TFAK(NDEPTH,NFRE))
          IF (.NOT.ALLOCATED(TSIHKD)) ALLOCATE(TSIHKD(NDEPTH,NFRE))
          IF (.NOT.ALLOCATED(TFAC_ST)) ALLOCATE(TFAC_ST(NDEPTH,NFRE))

          IF (.NOT.ALLOCATED(FAC0))                                     &
     &       ALLOCATE(FAC0(NANG,NANG,NFREHF,NFREHF))
          IF (.NOT.ALLOCATED(FAC1))                                     &
     &       ALLOCATE(FAC1(NANG,NANG,NFREHF,NFREHF))
          IF (.NOT.ALLOCATED(FAC2))                                     &
     &       ALLOCATE(FAC2(NANG,NANG,NFREHF,NFREHF))
          IF (.NOT.ALLOCATED(FAC3))                                     &
     &       ALLOCATE(FAC3(NANG,NANG,NFREHF,NFREHF))
          IF (.NOT.ALLOCATED(FAK)) ALLOCATE(FAK(NFREHF))
          IF (.NOT.ALLOCATED(FRHF)) ALLOCATE(FRHF(NFREHF))
          IF (.NOT.ALLOCATED(DFIMHF)) ALLOCATE(DFIMHF(NFREHF))

          IF (.NOT.ALLOCATED(OMEGA)) ALLOCATE(OMEGA(NFREH))
          IF (.NOT.ALLOCATED(THH))   ALLOCATE(THH(NANGH))
          IF (.NOT.ALLOCATED(DFDTH)) ALLOCATE(DFDTH(NFREH))
          IF (.NOT.ALLOCATED(TA)) ALLOCATE(TA(NDEPTH,NANGH,NFREH,NFREH))
          IF (.NOT.ALLOCATED(TB)) ALLOCATE(TB(NDEPTH,NANGH,NFREH,NFREH))
          IF (.NOT.ALLOCATED(TC_QL))                                    &
     &                    ALLOCATE(TC_QL(NDEPTH,NANGH,NFREH,NFREH))
          IF (.NOT.ALLOCATED(TT_4M))                                    &
     &                    ALLOCATE(TT_4M(NDEPTH,NANGH,NFREH,NFREH))
          IF (.NOT.ALLOCATED(TT_4P))                                    &
     &                    ALLOCATE(TT_4P(NDEPTH,NANGH,NFREH,NFREH))
          IF (.NOT.ALLOCATED(IM_P))                                     &
     &                    ALLOCATE(IM_P(NFREH,NFREH))
          IF (.NOT.ALLOCATED(IM_M))                                     &
     &                    ALLOCATE(IM_M(NFREH,NFREH))
          IF (.NOT.ALLOCATED(TFAKH)) ALLOCATE(TFAKH(NFREH,NDEPTH))

        ELSE 
          KCOUNT=0
          IKCOUNT=0
          DO M=1,NFRE
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FR(M)
          ENDDO
          DO M=1,NFRE
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=DFIM(M)
          ENDDO
          DO M=1,NFRE
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=GOM(M)
          ENDDO
          DO M=1,NFRE
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=C(M)
          ENDDO
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=DELTH
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=DELTR
          DO K=1,NANG
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=TH(K)
          ENDDO
          DO K=1,NANG
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=COSTH(K)
          ENDDO
          DO K=1,NANG
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=SINTH(K)
          ENDDO

          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=DELPHI
          DO J=1,NGY
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=DELLAM(J)
          ENDDO
          DO J=1,NGY
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=NLONRGG(J)
          ENDDO
          DO J=1,NGY
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=SINPH(J)
          ENDDO
          DO J=1,NGY
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=COSPH(J)
          ENDDO
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=IJS
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=IJL

          DO IJ=1,NIBLO
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=BLK2GLO(IJ)%IXLG
          ENDDO
          DO IJ=1,NIBLO
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=BLK2GLO(IJ)%KXLT
          ENDDO
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NX
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NY
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=IPER
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=AMOWEP
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=AMOSOP
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=AMOEAP
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=AMONOP
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=XDELLA
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=XDELLO
          DO J=1,NGY
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=ZDELLO(J)
          ENDDO
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=IRGG
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=IQGAUSS

          DO M=MFRSTLW,MLSTHG
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=IKP(M)
          ENDDO
          DO M=MFRSTLW,MLSTHG
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=IKP1(M)
          ENDDO
          DO M=MFRSTLW,MLSTHG
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=IKM(M)
          ENDDO
          DO M=MFRSTLW,MLSTHG
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=IKM1(M)
          ENDDO
          DO IC=1,2
            DO K=1,NANG
              IKCOUNT=IKCOUNT+1
              ICOMBUF(IKCOUNT)=K1W(K,IC)
            ENDDO
          ENDDO
          DO IC=1,2
            DO K=1,NANG
              IKCOUNT=IKCOUNT+1
              ICOMBUF(IKCOUNT)=K2W(K,IC)
            ENDDO
          ENDDO
          DO IC=1,2
            DO K=1,NANG
              IKCOUNT=IKCOUNT+1
              ICOMBUF(IKCOUNT)=K11W(K,IC)
            ENDDO
          ENDDO
          DO IC=1,2
            DO K=1,NANG
              IKCOUNT=IKCOUNT+1
              ICOMBUF(IKCOUNT)=K21W(K,IC)
            ENDDO
          ENDDO
          DO M=MFRSTLW,MLSTHG
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=AF11(M)
          ENDDO
          DO M=MFRSTLW,MLSTHG
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FKLAP(M)
          ENDDO
          DO M=MFRSTLW,MLSTHG
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FKLAP1(M)
          ENDDO
          DO M=MFRSTLW,MLSTHG
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FKLAM(M)
          ENDDO
          DO M=MFRSTLW,MLSTHG
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FKLAM1(M)
          ENDDO
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=ACL1
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=ACL2
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=CL11
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=CL21
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=DAL1
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=DAL2
          DO IC=1,KFRH
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FRH(IC)
          ENDDO

          DO NGOU=1,NGOUT
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=IJAR(NGOU)
          ENDDO

          DO IJ=1,NIBLO
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=DEPTH_INPUT(IJ)
          ENDDO
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=DEPTHA
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=DEPTHD
          DO M=1,NFRE
            DO KDEPTH=1,NDEPTH
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=TCGOND(KDEPTH,M)
            ENDDO
          ENDDO
          DO M=1,NFRE
            DO KDEPTH=1,NDEPTH
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=TFAK(KDEPTH,M)
            ENDDO
          ENDDO
          DO M=1,NFRE
            DO KDEPTH=1,NDEPTH
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=TSIHKD(KDEPTH,M)
            ENDDO
          ENDDO
          DO M=1,NFRE
            DO KDEPTH=1,NDEPTH
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=TFAC_ST(KDEPTH,M)
            ENDDO
          ENDDO

          DO M2=1,NFREHF
            DO M1=1,NFREHF
              DO K2=1,NANG
                DO K1=1,NANG
                  KCOUNT=KCOUNT+1
                  ZCOMBUF(KCOUNT)=FAC0(K1,K2,M1,M2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO M2=1,NFREHF
            DO M1=1,NFREHF
              DO K2=1,NANG
                DO K1=1,NANG
                  KCOUNT=KCOUNT+1
                  ZCOMBUF(KCOUNT)=FAC1(K1,K2,M1,M2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO M2=1,NFREHF
            DO M1=1,NFREHF
              DO K2=1,NANG
                DO K1=1,NANG
                  KCOUNT=KCOUNT+1
                  ZCOMBUF(KCOUNT)=FAC2(K1,K2,M1,M2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO M2=1,NFREHF
            DO M1=1,NFREHF
              DO K2=1,NANG
                DO K1=1,NANG
                  KCOUNT=KCOUNT+1
                  ZCOMBUF(KCOUNT)=FAC3(K1,K2,M1,M2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO M=1,NFREHF
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FAK(M)
          ENDDO

          DO M=1,NFREHF
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FRHF(M)
          ENDDO

          DO M=1,NFREHF
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=DFIMHF(M)
          ENDDO

          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=XMR
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=XMA

          DO M=1,NFREH
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=OMEGA(M)
          ENDDO

          DO K=1,NANGH
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=THH(K)
          ENDDO

          DO M=1,NFREH
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=DFDTH(M)
          ENDDO

          DO M1=1,NFREH
            DO M2=1,NFREH
              DO K1=1,NANGH
                DO KDEPTH=1,NDEPTH
                   KCOUNT=KCOUNT+1
                   ZCOMBUF(KCOUNT)=TA(KDEPTH,K1,M2,M1)
                   KCOUNT=KCOUNT+1
                   ZCOMBUF(KCOUNT)=TB(KDEPTH,K1,M2,M1)
                   KCOUNT=KCOUNT+1
                   ZCOMBUF(KCOUNT)=TC_QL(KDEPTH,K1,M2,M1)
                   KCOUNT=KCOUNT+1
                   ZCOMBUF(KCOUNT)=TT_4M(KDEPTH,K1,M2,M1)
                   KCOUNT=KCOUNT+1
                   ZCOMBUF(KCOUNT)=TT_4P(KDEPTH,K1,M2,M1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO M1=1,NFREH
            DO M2=1,NFREH
               IKCOUNT=IKCOUNT+1
               ICOMBUF(IKCOUNT)=IM_P(M2,M1)
               IKCOUNT=IKCOUNT+1
               ICOMBUF(IKCOUNT)=IM_M(M2,M1)
            ENDDO
          ENDDO
 
          DO KDEPTH=1,NDEPTH
            DO M=1,NFREH
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=TFAKH(M,KDEPTH)
            ENDDO
          ENDDO

          IF (IKCOUNT /= MIC) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
            WRITE (IU06,*) '* IKCOUNT NE MIC PRIOR   *'
            WRITE (IU06,*) '* CALL TO MPL_BROADCAST  *'
            WRITE (IU06,*) '* IKCOUNT =',IKCOUNT
            WRITE (IU06,*) '* MIC =',MIC
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF 
          IF (KCOUNT /= MZC) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
            WRITE (IU06,*) '* KCOUNT NE MZC PRIOR    *'
            WRITE (IU06,*) '* CALL TO MPL_BROADCAST  *'
            WRITE (IU06,*) '* KCOUNT =',KCOUNT
            WRITE (IU06,*) '* MZC =',MZC
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF 
        ENDIF

        CALL MPL_BROADCAST(ICOMBUF,KROOT=ISEND,KTAG=ITAG,               &
     &                     CDSTRING='MPBCASTGRID 1:')
        ITAG=ITAG+1

        CALL MPL_BROADCAST(ZCOMBUF,KROOT=ISEND,KTAG=ITAG,               &
     &                     CDSTRING='MPBCASTGRID 2:')
        ITAG=ITAG+1

        IF (IRANK /= ISEND) THEN
          KCOUNT=0
          IKCOUNT=0
          DO M=1,NFRE
            KCOUNT=KCOUNT+1
            FR(M)=ZCOMBUF(KCOUNT)
          ENDDO
          DO M=1,NFRE
            KCOUNT=KCOUNT+1
            DFIM(M)=ZCOMBUF(KCOUNT)
          ENDDO
          DO M=1,NFRE
            KCOUNT=KCOUNT+1
            GOM(M)=ZCOMBUF(KCOUNT)
          ENDDO
          DO M=1,NFRE
            KCOUNT=KCOUNT+1
            C(M)=ZCOMBUF(KCOUNT)
          ENDDO
          KCOUNT=KCOUNT+1
          DELTH=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          DELTR=ZCOMBUF(KCOUNT)
          DO K=1,NANG
            KCOUNT=KCOUNT+1
            TH(K)=ZCOMBUF(KCOUNT)
          ENDDO
          DO K=1,NANG
            KCOUNT=KCOUNT+1
            COSTH(K)=ZCOMBUF(KCOUNT)
          ENDDO
          DO K=1,NANG
            KCOUNT=KCOUNT+1
            SINTH(K)=ZCOMBUF(KCOUNT)
          ENDDO

          KCOUNT=KCOUNT+1
          DELPHI=ZCOMBUF(KCOUNT)
          DO J=1,NGY
            KCOUNT=KCOUNT+1
            DELLAM(J)=ZCOMBUF(KCOUNT)
          ENDDO
          DO J=1,NGY
            IKCOUNT=IKCOUNT+1
            NLONRGG(J)=ICOMBUF(IKCOUNT)
          ENDDO
          DO J=1,NGY
            KCOUNT=KCOUNT+1
            SINPH(J)=ZCOMBUF(KCOUNT)
          ENDDO
          DO J=1,NGY
            KCOUNT=KCOUNT+1
            COSPH(J)=ZCOMBUF(KCOUNT)
          ENDDO
          IKCOUNT=IKCOUNT+1
          IJS=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          IJL=ICOMBUF(IKCOUNT)

          DO IJ=1,NIBLO
            IKCOUNT=IKCOUNT+1
            BLK2GLO(IJ)%IXLG=ICOMBUF(IKCOUNT)
          ENDDO
          DO IJ=1,NIBLO
            IKCOUNT=IKCOUNT+1
            BLK2GLO(IJ)%KXLT=ICOMBUF(IKCOUNT)
          ENDDO
          IKCOUNT=IKCOUNT+1
          NX=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NY=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          IPER=ICOMBUF(IKCOUNT)
          KCOUNT=KCOUNT+1
          AMOWEP=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          AMOSOP=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          AMOEAP=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          AMONOP=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          XDELLA=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          XDELLO=ZCOMBUF(KCOUNT)
          DO J=1,NGY
            KCOUNT=KCOUNT+1
            ZDELLO(J)=ZCOMBUF(KCOUNT)
          ENDDO
          IKCOUNT=IKCOUNT+1
          IRGG=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          IQGAUSS=ICOMBUF(IKCOUNT)


          DO M=MFRSTLW,MLSTHG
            IKCOUNT=IKCOUNT+1
            IKP(M)=ICOMBUF(IKCOUNT)
          ENDDO
          DO M=MFRSTLW,MLSTHG
            IKCOUNT=IKCOUNT+1
            IKP1(M)=ICOMBUF(IKCOUNT)
          ENDDO
          DO M=MFRSTLW,MLSTHG
            IKCOUNT=IKCOUNT+1
            IKM(M)=ICOMBUF(IKCOUNT)
          ENDDO
          DO M=MFRSTLW,MLSTHG
            IKCOUNT=IKCOUNT+1
            IKM1(M)=ICOMBUF(IKCOUNT)
          ENDDO
          DO IC=1,2
            DO K=1,NANG
              IKCOUNT=IKCOUNT+1
              K1W(K,IC)=ICOMBUF(IKCOUNT)
            ENDDO
          ENDDO
          DO IC=1,2
            DO K=1,NANG
              IKCOUNT=IKCOUNT+1
              K2W(K,IC)=ICOMBUF(IKCOUNT)
            ENDDO
          ENDDO
          DO IC=1,2
            DO K=1,NANG
              IKCOUNT=IKCOUNT+1
              K11W(K,IC)=ICOMBUF(IKCOUNT)
            ENDDO
          ENDDO
          DO IC=1,2
            DO K=1,NANG
              IKCOUNT=IKCOUNT+1
              K21W(K,IC)=ICOMBUF(IKCOUNT)
            ENDDO
          ENDDO
          DO M=MFRSTLW,MLSTHG
            KCOUNT=KCOUNT+1
            AF11(M)=ZCOMBUF(KCOUNT)
          ENDDO
          DO M=MFRSTLW,MLSTHG
            KCOUNT=KCOUNT+1
            FKLAP(M)=ZCOMBUF(KCOUNT)
          ENDDO
          DO M=MFRSTLW,MLSTHG
            KCOUNT=KCOUNT+1
            FKLAP1(M)=ZCOMBUF(KCOUNT)
          ENDDO
          DO M=MFRSTLW,MLSTHG
            KCOUNT=KCOUNT+1
            FKLAM(M)=ZCOMBUF(KCOUNT)
          ENDDO
          DO M=MFRSTLW,MLSTHG
            KCOUNT=KCOUNT+1
            FKLAM1(M)=ZCOMBUF(KCOUNT)
          ENDDO
          KCOUNT=KCOUNT+1
          ACL1=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          ACL2=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          CL11=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          CL21=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          DAL1=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          DAL2=ZCOMBUF(KCOUNT)
          DO IC=1,KFRH
            KCOUNT=KCOUNT+1
            FRH(IC)=ZCOMBUF(KCOUNT)
          ENDDO

          DO NGOU=1,NGOUT
            IKCOUNT=IKCOUNT+1
            IJAR(NGOU)=ICOMBUF(IKCOUNT)
          ENDDO

          DO IJ=1,NIBLO
            KCOUNT=KCOUNT+1
            DEPTH_INPUT(IJ)=ZCOMBUF(KCOUNT)
          ENDDO
          KCOUNT=KCOUNT+1
          DEPTHA=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          DEPTHD=ZCOMBUF(KCOUNT)
          DO M=1,NFRE
            DO KDEPTH=1,NDEPTH
              KCOUNT=KCOUNT+1
              TCGOND(KDEPTH,M)=ZCOMBUF(KCOUNT)
            ENDDO
          ENDDO
          DO M=1,NFRE
            DO KDEPTH=1,NDEPTH
              KCOUNT=KCOUNT+1
              TFAK(KDEPTH,M)=ZCOMBUF(KCOUNT)
            ENDDO
          ENDDO
          DO M=1,NFRE
            DO KDEPTH=1,NDEPTH
              KCOUNT=KCOUNT+1
              TSIHKD(KDEPTH,M)=ZCOMBUF(KCOUNT)
            ENDDO
          ENDDO
          DO M=1,NFRE
            DO KDEPTH=1,NDEPTH
              KCOUNT=KCOUNT+1
              TFAC_ST(KDEPTH,M)=ZCOMBUF(KCOUNT)
            ENDDO
          ENDDO

          DO M2=1,NFREHF
            DO M1=1,NFREHF
              DO K2=1,NANG
                DO K1=1,NANG
                    KCOUNT=KCOUNT+1
                    FAC0(K1,K2,M1,M2)=ZCOMBUF(KCOUNT)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO M2=1,NFREHF
            DO M1=1,NFREHF
              DO K2=1,NANG
                DO K1=1,NANG
                    KCOUNT=KCOUNT+1
                    FAC1(K1,K2,M1,M2)=ZCOMBUF(KCOUNT)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO M2=1,NFREHF
            DO M1=1,NFREHF
              DO K2=1,NANG
                DO K1=1,NANG
                    KCOUNT=KCOUNT+1
                    FAC2(K1,K2,M1,M2)=ZCOMBUF(KCOUNT)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO M2=1,NFREHF
            DO M1=1,NFREHF
              DO K2=1,NANG
                DO K1=1,NANG
                    KCOUNT=KCOUNT+1
                    FAC3(K1,K2,M1,M2)=ZCOMBUF(KCOUNT)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO M=1,NFREHF
            KCOUNT=KCOUNT+1
            FAK(M)=ZCOMBUF(KCOUNT)
          ENDDO

          DO M=1,NFREHF
            KCOUNT=KCOUNT+1
            FRHF(M)=ZCOMBUF(KCOUNT)
          ENDDO

          DO M=1,NFREHF
            KCOUNT=KCOUNT+1
            DFIMHF(M)=ZCOMBUF(KCOUNT)
          ENDDO

          KCOUNT=KCOUNT+1
          XMR=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          XMA=ZCOMBUF(KCOUNT)

          DO M=1,NFREH
            KCOUNT=KCOUNT+1
            OMEGA(M)=ZCOMBUF(KCOUNT)
          ENDDO

          DO K=1,NANGH
            KCOUNT=KCOUNT+1
            THH(K)=ZCOMBUF(KCOUNT)
          ENDDO

          DO M=1,NFREH
            KCOUNT=KCOUNT+1
            DFDTH(M)=ZCOMBUF(KCOUNT)
          ENDDO

          DO M1=1,NFREH
            DO M2=1,NFREH
              DO K1=1,NANGH
                DO KDEPTH=1,NDEPTH
                   KCOUNT=KCOUNT+1
                   TA(KDEPTH,K1,M2,M1)=ZCOMBUF(KCOUNT)
                   KCOUNT=KCOUNT+1
                   TB(KDEPTH,K1,M2,M1)=ZCOMBUF(KCOUNT)
                   KCOUNT=KCOUNT+1
                   TC_QL(KDEPTH,K1,M2,M1)=ZCOMBUF(KCOUNT)
                   KCOUNT=KCOUNT+1
                   TT_4M(KDEPTH,K1,M2,M1)=ZCOMBUF(KCOUNT)
                   KCOUNT=KCOUNT+1
                   TT_4P(KDEPTH,K1,M2,M1)=ZCOMBUF(KCOUNT)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO M1=1,NFREH
            DO M2=1,NFREH
               IKCOUNT=IKCOUNT+1
               IM_P(M2,M1)=ICOMBUF(IKCOUNT)
               IKCOUNT=IKCOUNT+1
               IM_M(M2,M1)=ICOMBUF(IKCOUNT)
            ENDDO
          ENDDO

          DO KDEPTH=1,NDEPTH
            DO M=1,NFREH
              KCOUNT=KCOUNT+1
              TFAKH(M,KDEPTH)=ZCOMBUF(KCOUNT)
            ENDDO
          ENDDO

          IF (IKCOUNT /= MIC) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
            WRITE (IU06,*) '* IKCOUNT NE MIC AFTER   *'
            WRITE (IU06,*) '* CALL TO MPL_BROADCAST  *'
            WRITE (IU06,*) '* IKCOUNT =',IKCOUNT
            WRITE (IU06,*) '* MIC =',MIC
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF 
          IF (KCOUNT /= MZC) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
            WRITE (IU06,*) '* KCOUNT NE MZC AFTER    *'
            WRITE (IU06,*) '* CALL TO MPL_BROADCAST  *'
            WRITE (IU06,*) '* KCOUNT =',KCOUNT
            WRITE (IU06,*) '* MZC =',MZC
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF 
        ENDIF

        DEALLOCATE(ICOMBUF,ZCOMBUF)

      ENDIF

      IF (LHOOK) CALL DR_HOOK('MPBCASTGRID',1,ZHOOK_HANDLE)

END SUBROUTINE MPBCASTGRID
