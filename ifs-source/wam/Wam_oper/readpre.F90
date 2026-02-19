! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE READPRE (IU07)

! ----------------------------------------------------------------------

!**** *READPRE*  READ GRID OUTPUT FROM PREPROC.

!     H. GUNTHER      GKSS/ECMWF     MAY 1990
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     J. BIDLOT            ECMWF       11/2003
!                          IF YOUR ARE RUNNING AT ECMWF:
!                          BE AWARE THAT IF YOU CHANGE ANYTHING TO THE
!                          STRUCTURE OF THE OUTPUT FILE YOU WILL HAVE TO
!                          MAKE SURE THAT IT IS CREATED FOR YOUR RUN,
!                          OTHERWISE IT MIGHT PICK UP THE DEFAULT ONE
!                          THAT IS ALREADY ON DISK.
!                          YOU ALSO HAVE TO CHANGE OUTCOM.
! ALSO CHECK
!            READPREB in Alt
!  and
!            ... /nemo/tools/interpolate/wambingrid.F90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*    PURPOSE.
!     --------

!       INPUT OF PREPROC GRID OUTPUT.

!**   INTERFACE.
!     ----------

!       *CALL* *READPRE (IU07)*
!          *IU07 *  - INPUT UNIT OF PREPROC GRID FILE.

!     METHOD.
!     -------

!       UNFORMATED READ FROM UNIT.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

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
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED ,              &
     &            NGX      ,NGY      ,                                  &
     &            NIBLO    ,NOVER    ,NIBL1    ,CLDOMAIN ,IMDLGRDID
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTH_INPUT,DEPTHA   ,DEPTHD   ,  &
     &            TCGOND   ,TFAK     ,TSIHKD   ,TFAC_ST  ,TOOSHALLOW
      USE YOWTABL  , ONLY : FAC0     ,FAC1     ,FAC2     ,FAC3     ,    &
     &            FAK      ,FRHF     ,DFIMHF   ,NFREHF   ,              &
     &            MR       ,XMR      ,MA       ,XMA      ,NFREH    ,    &
     &            NANGH    ,NMAX     ,OMEGA    ,DFDTH    ,THH      ,    &
     &            DELTHH   ,IM_P     ,IM_M     ,TA       ,TB       ,    &
     &            TC_QL    ,TT_4M    ,TT_4P    ,TFAKH
      USE YOWTEST  , ONLY : IU06
      USE YOWUNPOOL, ONLY : LLUNSTR  ,LPREPROC, LLR8TOR4
      USE UNWAM     ,ONLY : INIT_UNWAM, UNWAM_IN, SET_UNWAM_HANDLES
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "mpbcastgrid.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU07
      INTEGER(KIND=JWIM) :: IREAD
      INTEGER(KIND=JWIM) :: IDUM, KIBLD, KIBLC
      INTEGER(KIND=JWIM) :: KMDLGRDID, KMDLGRBID_G, KMDLGRBID_M
      INTEGER(KIND=JWIM) :: NKIND !Precision of file when reading

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('READPRE',0,ZHOOK_HANDLE)

      NKIND=0
      IREAD=1
      KTAG=1

      CALL GSTATS(1771,0)
      IF (IRANK == IREAD) THEN
!       READ MODEL IDENTIFIERS
        CALL READREC(1)
        IF (KMDLGRDID /= IMDLGRDID) THEN
          WRITE(IU06,*) '*****************************************'
          WRITE(IU06,*) '*                                       *'
          WRITE(IU06,*) '*  FATAL ERROR(S) IN SUB. READPRE       *'
          WRITE(IU06,*) '*  ==============================       *'
          WRITE(IU06,*) '*                                       *'
          WRITE(IU06,*) '* THE PROGRAM HAS DETECTED DIFFERENT    *'
          WRITE(IU06,*) '* MODEL GRID IDENTIFIER.                *' 
          WRITE(IU06,*) '* MAKE SURE YOU HAVE RUN PREPROC !!!!   *'
          WRITE(IU06,*)    KMDLGRDID, IMDLGRDID
          WRITE(IU06,*) '*                                       *'
          WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.     *'
          WRITE(IU06,*) '* ---------------   --------------      *'
          WRITE(IU06,*) '*****************************************'
          CALL ABORT1
        ENDIF
        IF (NKIND /= KIND(DELPHI)) THEN
          WRITE(IU06,*) '*****************************************'
          WRITE(IU06,*) '*                                       *'
          WRITE(IU06,*) '*  FATAL ERROR(S) IN SUB. READPRE       *'
          WRITE(IU06,*) '*  ==============================       *'
          WRITE(IU06,*) '*                                       *'
          WRITE(IU06,*) '* THE PROGRAM HAS DETECTED DIFFERENT    *'
          WRITE(IU06,*) '* PRECISION IN FILE AND MODEL.          *' 
          WRITE(IU06,*)    NKIND, KIND(DELPHI)
          WRITE(IU06,*) '*                                       *'
          WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.     *'
          WRITE(IU06,*) '* ---------------   --------------      *'
          WRITE(IU06,*) '*****************************************'
          CALL ABORT1
        ENDIF

!*    0. READ YOWPARAM (BLOCK SIZES). 
!        ----------------------------

        CALL READREC(2)

!*    1. READ MODULE YOWFRED (FREQUENCY DIRECTION GRID).
!        ----------------------------------------------

        IF (.NOT.ALLOCATED(FR)) ALLOCATE(FR(NFRE))
        IF (.NOT.ALLOCATED(DFIM)) ALLOCATE(DFIM(NFRE))
        IF (.NOT.ALLOCATED(GOM)) ALLOCATE(GOM(NFRE))
        IF (.NOT.ALLOCATED(C)) ALLOCATE(C(NFRE))
        IF (.NOT.ALLOCATED(TH)) ALLOCATE(TH(NANG))
        IF (.NOT.ALLOCATED(COSTH)) ALLOCATE(COSTH(NANG))
        IF (.NOT.ALLOCATED(SINTH)) ALLOCATE(SINTH(NANG))

        CALL READREC(3)

!*    2. READ MODULE YOWGRID (GENERAL GRID ORGANISATION).
!        ------------------------------------------------

        IF (.NOT.ALLOCATED(DELLAM)) ALLOCATE(DELLAM(NGY))
        IF (.NOT.ALLOCATED(NLONRGG)) ALLOCATE(NLONRGG(NGY))
        IF (.NOT.ALLOCATED(SINPH)) ALLOCATE(SINPH(NGY))
        IF (.NOT.ALLOCATED(COSPH)) ALLOCATE(COSPH(NGY))

        CALL READREC(4)

!*    3. READ MODULE YOWMAP (LONG. AND LAT. INDICES OF GRID POINTS).
!        --------------------------------------------------------

        IF (ALLOCATED(BLK2GLO)) DEALLOCATE(BLK2GLO)
        ALLOCATE(BLK2GLO(NIBLO))

        IF (.NOT.ALLOCATED(ZDELLO)) ALLOCATE(ZDELLO(NGY))

        CALL READREC(5)

!     DETERMINE IF WE ARE USING A QUASI GAUSSIAN GRID OR 
!     LAT-LONG GRID (REGULAR OR IRREGULAR).

        IF (IPER ==1 .AND. AMONOP == ABS(AMOSOP) .AND.                  &
     &     MOD(NY,2) == 0 .AND. IRGG == 1 ) THEN
          IQGAUSS=1
        ELSE
          IQGAUSS=0
        ENDIF

!*    4. READ MODULE YOWINDNL (NON-LINEAR INTERACTION).
!       -----------------------------------------------

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

        CALL READREC(6)


!*    7. READ MODULE YOWCOUT (INDICES OF OUTPUT POINTS).
!        --------------------------------------------

        CALL READREC(12)
        IF (NGOUT > 0) THEN
          IF (.NOT.ALLOCATED(IJAR)) ALLOCATE(IJAR(NGOUT))
          CALL READREC(13)
        ENDIF

!*    8. READ MODULE YOWSHAL (DEPTH AND SHALLOW WATER TABLES).
!        ----------------------------------------------------

        CALL READREC(14)
!!!   note that the size of DEPTH_INPUT will be readjusted in mpdecomp !!!
        IF (ALLOCATED(DEPTH_INPUT)) DEALLOCATE(DEPTH_INPUT)
        ALLOCATE(DEPTH_INPUT(NIBLO))
        IF (.NOT.ALLOCATED(TCGOND)) ALLOCATE(TCGOND(NDEPTH,NFRE))
        IF (.NOT.ALLOCATED(TFAK)) ALLOCATE(TFAK(NDEPTH,NFRE))
        IF (.NOT.ALLOCATED(TSIHKD)) ALLOCATE(TSIHKD(NDEPTH,NFRE))
        IF (.NOT.ALLOCATED(TFAC_ST)) ALLOCATE(TFAC_ST(NDEPTH,NFRE))

        CALL READREC(15)

!*    9. READ MODULE YOWTABL (2ND AND 3RD part).
!        -------------------

        IF (.NOT.ALLOCATED(FAC0)) ALLOCATE(FAC0(NANG,NANG,NFREHF,NFREHF))
        IF (.NOT.ALLOCATED(FAC1)) ALLOCATE(FAC1(NANG,NANG,NFREHF,NFREHF))
        IF (.NOT.ALLOCATED(FAC2)) ALLOCATE(FAC2(NANG,NANG,NFREHF,NFREHF))
        IF (.NOT.ALLOCATED(FAC3)) ALLOCATE(FAC3(NANG,NANG,NFREHF,NFREHF))
        IF (.NOT.ALLOCATED(FAK)) ALLOCATE(FAK(NFREHF))
        IF (.NOT.ALLOCATED(FRHF)) ALLOCATE(FRHF(NFREHF))
        IF (.NOT.ALLOCATED(DFIMHF)) ALLOCATE(DFIMHF(NFREHF))

        CALL READREC(16)

        CALL READREC(17)

        IF (.NOT.ALLOCATED(OMEGA)) ALLOCATE(OMEGA(NFREH))
        IF (.NOT.ALLOCATED(THH))   ALLOCATE(THH(NANGH))
        IF (.NOT.ALLOCATED(DFDTH)) ALLOCATE(DFDTH(NFREH))
        IF (.NOT.ALLOCATED(TA)) ALLOCATE(TA(NDEPTH,NANGH,NFREH,NFREH))
        IF (.NOT.ALLOCATED(TB)) ALLOCATE(TB(NDEPTH,NANGH,NFREH,NFREH))
        IF (.NOT.ALLOCATED(TC_QL))                                       &
     &         ALLOCATE(TC_QL(NDEPTH,NANGH,NFREH,NFREH))
        IF (.NOT.ALLOCATED(TT_4M))                                       &
     &         ALLOCATE(TT_4M(NDEPTH,NANGH,NFREH,NFREH))
        IF (.NOT.ALLOCATED(TT_4P))                                       &
     &         ALLOCATE(TT_4P(NDEPTH,NANGH,NFREH,NFREH))
        IF (.NOT.ALLOCATED(IM_P))                                        &
     &         ALLOCATE(IM_P(NFREH,NFREH))
        IF (.NOT.ALLOCATED(IM_M))                                        &
     &         ALLOCATE(IM_M(NFREH,NFREH))
        IF (.NOT.ALLOCATED(TFAKH)) ALLOCATE(TFAKH(NFREH,NDEPTH))

        CALL READREC(18)

!       THE UNSTRUCTURED BITS (if pre-computed by PREPROC)
        IF (LLUNSTR .AND. LPREPROC) THEN
          CALL SET_UNWAM_HANDLES
          CALL UNWAM_IN(IU07)
        END IF

      ENDIF

      CLOSE (UNIT=IU07)

      CALL GSTATS(1771,1)

!     SEND INFORMATION FROM READPRE TO ALL PE's
      CALL GSTATS(694,0)
      CALL MPBCASTGRID(IU06,IREAD,KTAG)
      CALL GSTATS(694,1)

      TOOSHALLOW=0.1_JWRB*DEPTHA

!     RECALCULATE THE UNSTRUCTURED BITS (if not read in)

      IF (LLUNSTR .AND. .NOT. LPREPROC) THEN
        CALL INIT_UNWAM
      ENDIF

      IF (LHOOK) CALL DR_HOOK('READPRE',1,ZHOOK_HANDLE)

      RETURN

      CONTAINS

      SUBROUTINE READREC(KREC)
      IMPLICIT NONE
      INTEGER(KIND=JWIM), INTENT(IN) :: KREC
      INTEGER(KIND=JWIM) :: ISTAT
      REAL(KIND=JWRU) ::                                                &
     & R8_ACL1,R8_ACL2,                                                 &
     & R8_AMOEAP,R8_AMONOP,R8_AMOSOP,R8_AMOWEP,                         &
     & R8_CL11,R8_CL21,R8_DAL1,R8_DAL2,R8_DELPHI,R8_DELTH,R8_DELTHH,    &
     & R8_DELTR,R8_DEPTHA,R8_DEPTHD,                                    &
     & R8_XDELLA,R8_XDELLO,R8_XMA,R8_XMR
      REAL(KIND=JWRU), ALLOCATABLE, DIMENSION(:) ::                     &
     &     R8_FR,R8_DFIM,R8_GOM,R8_C,                                   &
     &     R8_TH,R8_COSTH,R8_SINTH,                                     &
     &     R8_DELLAM,R8_SINPH,R8_COSPH,                                 &
     &     R8_ZDELLO,                                                   &
     &     R8_AF11,                                                     &
     &     R8_FKLAP,                                                    &
     &     R8_FKLAP1,                                                   &
     &     R8_FKLAM,                                                    &
     &     R8_FKLAM1,                                                   &
     &     R8_FRH,                                                      &
     &     R8_FRHF,                                                     &
     &     R8_DFIMHF,                                                   &
     &     R8_OMEGA,                                                    &
     &     R8_THH,                                                      &
     &     R8_DFDTH,                                                    &
     &     R8_DEPTH_INPUT,                                              &
     &     R8_FAK                
      REAL(KIND=JWRU), ALLOCATABLE, DIMENSION(:,:) ::                   &
     &     R8_TCGOND,                                                   &
     &     R8_TFAK,                                                     &
     &     R8_TSIHKD,                                                   &
     &     R8_TFAC_ST,                                                  &
     &     R8_TFAKH
      REAL(KIND=JWRU), ALLOCATABLE, DIMENSION(:,:,:,:) ::               &
     &     R8_FAC0,                                                     &
     &     R8_FAC1,                                                     &
     &     R8_FAC2,                                                     &
     &     R8_FAC3,                                                     &
     &     R8_TA,                                                       &
     &     R8_TB,                                                       &
     &     R8_TC_QL,                                                    &
     &     R8_TT_4M,                                                    &
     &     R8_TT_4P

!23456789-123456789-123456789-123456789-123456789-123456789-123456789-12
      SELECT CASE(KREC)
      CASE(1)
         LLR8TOR4 = .FALSE.
         READ(IU07,IOSTAT=ISTAT) NKIND,KMDLGRDID,KMDLGRBID_G,KMDLGRBID_M
         IF (ISTAT /= 0) THEN
            ! Is this legacy code i.e. NKIND missing from the 1st record ?
            BACKSPACE(IU07)
            READ(IU07,IOSTAT=ISTAT) KMDLGRDID,KMDLGRBID_G,KMDLGRBID_M
            IF (ISTAT /= 0) GOTO 1000 ! Abort: Don't support this file format at all
         ENDIF
         IF (KIND(DELPHI) == 4 .AND. NKIND /= 4) THEN
            LLR8TOR4 = .TRUE.   ! Input REALs are indeed legacy REAL*8, but KIND(DELPHI) == 4
         ENDIF
         NKIND = KIND(DELPHI)   ! Pretend NKIND is in "right precision"
         WRITE(IU06,1002) 'READPRE(READREC): NKIND=',NKIND,             &
     &        ', KIND(DELPHI)=',KIND(DELPHI),', LLR8TOR4=',LLR8TOR4
 1002    FORMAT(2X,A,I0,A,I0,A,L1)
      CASE(2)
         READ(IU07,IOSTAT=ISTAT)                                        &
     &        NANG, NFRE, NFRE_RED, NGX, NGY, NIBLO, NOVER,             &
     &        KFRH, MFRSTLW, MLSTHG,                                    &
     &        NIBL1, IDUM, KIBLD, KIBLC, CLDOMAIN
         IF (ISTAT /= 0) GOTO 1000
      CASE(3)
         IF (LLR8TOR4) THEN
            ALLOCATE(R8_FR(NFRE),R8_DFIM(NFRE),R8_GOM(NFRE),R8_C(NFRE))
            ALLOCATE(R8_TH(NANG),R8_COSTH(NANG),R8_SINTH(NANG))

            READ(IU07,IOSTAT=ISTAT) R8_FR, R8_DFIM, R8_GOM, R8_C,       &
     &           R8_DELTH, R8_DELTR, R8_TH, R8_COSTH, R8_SINTH
            IF (ISTAT /= 0) GOTO 1000

            FR = R8_FR
            DFIM = R8_DFIM
            GOM = R8_GOM
            C = R8_C
            DELTH = R8_DELTH
            DELTR = R8_DELTR
            TH = R8_TH
            COSTH = R8_COSTH
            SINTH = R8_SINTH
            DEALLOCATE(R8_FR,R8_DFIM,R8_GOM, R8_C)
            DEALLOCATE(R8_TH,R8_COSTH,R8_SINTH)
         ELSE
            READ(IU07,IOSTAT=ISTAT) FR, DFIM, GOM, C,                   &
     &           DELTH, DELTR, TH, COSTH, SINTH
            IF (ISTAT /= 0) GOTO 1000
         ENDIF
      CASE(4)
         IF (LLR8TOR4) THEN
            ALLOCATE(R8_DELLAM(NGY),R8_SINPH(NGY),R8_COSPH(NGY))

            READ(IU07,IOSTAT=ISTAT) R8_DELPHI, R8_DELLAM, NLONRGG,      &
     &           R8_SINPH, R8_COSPH,                                    &
     &           IJS, IJL
            IF (ISTAT /= 0) GOTO 1000

            DELPHI = R8_DELPHI
            DELLAM = R8_DELLAM
            SINPH = R8_SINPH
            COSPH = R8_COSPH
            DEALLOCATE(R8_DELLAM,R8_SINPH,R8_COSPH)
         ELSE
            READ(IU07,IOSTAT=ISTAT) DELPHI, DELLAM, NLONRGG,            &
     &           SINPH, COSPH,                                          &
     &           IJS, IJL
            IF (ISTAT /= 0) GOTO 1000
         ENDIF
      CASE(5)
         IF (LLR8TOR4) THEN
            ALLOCATE(R8_ZDELLO(NGY))

            READ(IU07,IOSTAT=ISTAT) BLK2GLO%IXLG, BLK2GLO%KXLT, NX, NY, IPER, &
     &           R8_AMOWEP, R8_AMOSOP, R8_AMOEAP, R8_AMONOP,                  &
     &           R8_XDELLA, R8_XDELLO,                                        &
     &           R8_ZDELLO, IRGG
            IF (ISTAT /= 0) GOTO 1000

            AMOWEP = R8_AMOWEP
            AMOSOP = R8_AMOSOP
            AMOEAP = R8_AMOEAP
            AMONOP = R8_AMONOP
            XDELLA = R8_XDELLA
            XDELLO = R8_XDELLO
            ZDELLO = R8_ZDELLO

            DEALLOCATE(R8_ZDELLO)
         ELSE
            READ(IU07,IOSTAT=ISTAT) BLK2GLO%IXLG, BLK2GLO%KXLT, NX, NY, IPER, &
     &           AMOWEP, AMOSOP, AMOEAP, AMONOP,                              &
     &           XDELLA, XDELLO,                                              &
     &           ZDELLO, IRGG
            IF (ISTAT /= 0) GOTO 1000
         ENDIF
      CASE(6)
         IF (LLR8TOR4) THEN
            ALLOCATE(R8_AF11(MFRSTLW:MLSTHG))
            ALLOCATE(R8_FKLAP(MFRSTLW:MLSTHG))
            ALLOCATE(R8_FKLAP1(MFRSTLW:MLSTHG))
            ALLOCATE(R8_FKLAM(MFRSTLW:MLSTHG))
            ALLOCATE(R8_FKLAM1(MFRSTLW:MLSTHG))
            ALLOCATE(R8_FRH(KFRH))

            READ(IU07,IOSTAT=ISTAT)                                     &
     &           IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W,            &
     &           R8_AF11, R8_FKLAP, R8_FKLAP1, R8_FKLAM, R8_FKLAM1,     &
     &           R8_ACL1, R8_ACL2,  R8_CL11, R8_CL21,                   &
     &           R8_DAL1, R8_DAL2, R8_FRH
            IF (ISTAT /= 0) GOTO 1000

            AF11 = R8_AF11
            FKLAP = R8_FKLAP
            FKLAP1 = R8_FKLAP1
            FKLAM = R8_FKLAM
            FKLAM1 = R8_FKLAM1
            ACL1 = R8_ACL1
            ACL2 = R8_ACL2
            CL11 = R8_CL11
            CL21 = R8_CL21
            DAL1 = R8_DAL1
            DAL2 = R8_DAL2
            FRH = R8_FRH

            DEALLOCATE(R8_AF11)
            DEALLOCATE(R8_FKLAP)
            DEALLOCATE(R8_FKLAP1)
            DEALLOCATE(R8_FKLAM)
            DEALLOCATE(R8_FKLAM1)
            DEALLOCATE(R8_FRH)
         ELSE
            READ(IU07,IOSTAT=ISTAT)                                     &
     &           IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W,            &
     &           AF11, FKLAP, FKLAP1, FKLAM, FKLAM1,                    &
     &           ACL1, ACL2,  CL11, CL21,                               &
     &           DAL1, DAL2, FRH
            IF (ISTAT /= 0) GOTO 1000
         ENDIF
      CASE(12)
         READ(IU07,IOSTAT=ISTAT) NGOUT
         IF (ISTAT /= 0) GOTO 1000
      CASE(13)
         READ(IU07,IOSTAT=ISTAT) IJAR
         IF (ISTAT /= 0) GOTO 1000
      CASE(14)
         IF (LLR8TOR4) THEN
            READ(IU07,IOSTAT=ISTAT) NDEPTH, R8_DEPTHA, R8_DEPTHD
            IF (ISTAT /= 0) GOTO 1000
            DEPTHA = R8_DEPTHA
            DEPTHD = R8_DEPTHD
         ELSE
            READ(IU07,IOSTAT=ISTAT) NDEPTH, DEPTHA, DEPTHD
            IF (ISTAT /= 0) GOTO 1000
         ENDIF
      CASE(15)
         IF (LLR8TOR4) THEN
            ALLOCATE(R8_DEPTH_INPUT(NIBLO))
            ALLOCATE(R8_TCGOND(NDEPTH,NFRE))
            ALLOCATE(R8_TFAK(NDEPTH,NFRE))
            ALLOCATE(R8_TSIHKD(NDEPTH,NFRE))
            ALLOCATE(R8_TFAC_ST(NDEPTH,NFRE))
            READ(IU07,IOSTAT=ISTAT) R8_DEPTH_INPUT, R8_TCGOND,          &
     &           R8_TFAK, R8_TSIHKD, R8_TFAC_ST
            IF (ISTAT /= 0) GOTO 1000
            DEPTH_INPUT(:) = R8_DEPTH_INPUT(:)
            TCGOND = R8_TCGOND
            TFAK = R8_TFAK
            TSIHKD = R8_TSIHKD
            TFAC_ST = R8_TFAC_ST

            DEALLOCATE(R8_DEPTH_INPUT)
            DEALLOCATE(R8_TCGOND)
            DEALLOCATE(R8_TFAK)
            DEALLOCATE(R8_TSIHKD)
            DEALLOCATE(R8_TFAC_ST)
         ELSE
            READ(IU07,IOSTAT=ISTAT) DEPTH_INPUT, TCGOND,                &
     &           TFAK, TSIHKD, TFAC_ST
            IF (ISTAT /= 0) GOTO 1000
         ENDIF
      CASE(16)
         IF (LLR8TOR4) THEN
            ALLOCATE(R8_FAC0(NANG,NANG,NFREHF,NFREHF))
            ALLOCATE(R8_FAC1(NANG,NANG,NFREHF,NFREHF))
            ALLOCATE(R8_FAC2(NANG,NANG,NFREHF,NFREHF))
            ALLOCATE(R8_FAC3(NANG,NANG,NFREHF,NFREHF))
            ALLOCATE(R8_FAK(NFREHF))
            ALLOCATE(R8_FRHF(NFREHF))
            ALLOCATE(R8_DFIMHF(NFREHF))
            READ(IU07,IOSTAT=ISTAT) R8_FAC0,R8_FAC1,R8_FAC2,R8_FAC3,    &
     &           R8_FAK,R8_FRHF,R8_DFIMHF
            IF (ISTAT /= 0) GOTO 1000
            FAC0 = R8_FAC0
            FAC1 = R8_FAC1
            FAC2 = R8_FAC2
            FAC3 = R8_FAC3
            FAK = R8_FAK
            FRHF = R8_FRHF
            DFIMHF = R8_DFIMHF

            DEALLOCATE(R8_FAC0)
            DEALLOCATE(R8_FAC1)
            DEALLOCATE(R8_FAC2)
            DEALLOCATE(R8_FAC3)
            DEALLOCATE(R8_FAK)
            DEALLOCATE(R8_FRHF)
            DEALLOCATE(R8_DFIMHF)
         ELSE
            READ(IU07,IOSTAT=ISTAT) FAC0,FAC1,FAC2,FAC3,                &
     &           FAK,FRHF,DFIMHF
            IF (ISTAT /= 0) GOTO 1000
         ENDIF
      CASE(17)
         IF (LLR8TOR4) THEN
            READ(IU07,IOSTAT=ISTAT) MR, R8_XMR, MA, R8_XMA,             &
     &           NFREH, NANGH, NMAX
            IF (ISTAT /= 0) GOTO 1000
            XMR = R8_XMR
            XMA = R8_XMA
         ELSE
            READ(IU07,IOSTAT=ISTAT) MR, XMR, MA, XMA,                   &
     &           NFREH, NANGH, NMAX
            IF (ISTAT /= 0) GOTO 1000
         ENDIF
      CASE(18)
         IF (LLR8TOR4) THEN
            ALLOCATE(R8_OMEGA(NFREH))
            ALLOCATE(R8_THH(NANGH))
            ALLOCATE(R8_DFDTH(NFREH))
            ALLOCATE(R8_TA(NDEPTH,NANGH,NFREH,NFREH))
            ALLOCATE(R8_TB(NDEPTH,NANGH,NFREH,NFREH))
            ALLOCATE(R8_TC_QL(NDEPTH,NANGH,NFREH,NFREH))
            ALLOCATE(R8_TT_4M(NDEPTH,NANGH,NFREH,NFREH))
            ALLOCATE(R8_TT_4P(NDEPTH,NANGH,NFREH,NFREH))
            ALLOCATE(R8_TFAKH(NFREH,NDEPTH))
            READ(IU07,IOSTAT=ISTAT) R8_OMEGA, R8_DFDTH, R8_THH,         &
     &           R8_DELTHH, IM_P, IM_M,                                 &
     &           R8_TA, R8_TB, R8_TC_QL, R8_TT_4M, R8_TT_4P, R8_TFAKH
            IF (ISTAT /= 0) GOTO 1000
            OMEGA = R8_OMEGA
            THH = R8_THH
            DFDTH = R8_DFDTH
            TA = R8_TA
            TB = R8_TB
            TC_QL = R8_TC_QL
            TT_4M = R8_TT_4M
            TT_4P = R8_TT_4P
            TFAKH = R8_TFAKH

            DEALLOCATE(R8_OMEGA)
            DEALLOCATE(R8_THH)
            DEALLOCATE(R8_DFDTH)
            DEALLOCATE(R8_TA)
            DEALLOCATE(R8_TB)
            DEALLOCATE(R8_TC_QL)
            DEALLOCATE(R8_TT_4M)
            DEALLOCATE(R8_TT_4P)
            DEALLOCATE(R8_TFAKH)
         ELSE
            READ(IU07,IOSTAT=ISTAT) OMEGA, DFDTH, THH,                  &
     &           DELTHH, IM_P, IM_M,                                    &
     &           TA, TB, TC_QL, TT_4M, TT_4P, TFAKH
            IF (ISTAT /= 0) GOTO 1000
         ENDIF
      CASE DEFAULT
         WRITE(IU06,*)'***ERROR IN READREC: INVALID RECORD NUMBER=',KREC
         CALL FLUSH(IU06)
         CALL ABORT1
      END SELECT
      RETURN

 1000 CONTINUE
      WRITE(IU06,1001)'***ERROR IN READREC(',KREC,') : IOSTAT=',        &
     &     ISTAT,', NKIND=',NKIND
 1001 FORMAT(1X,A,I0,A,I0,A,I0)
      WRITE(IU06,*) '*************************************'
      WRITE(IU06,*) '*                                   *'
      WRITE(IU06,*) '*  READ ERROR IN SUB. READREC       *'
      WRITE(IU06,*) '*  ==========================       *'
      WRITE(IU06,*) '*                                   *'
      WRITE(IU06,'(1X,A,I0)') '*  READ ERROR TO UNIT fort.',IU07
      WRITE(IU06,*) '*  IS THE FILE PRESENT ????         *' 
      WRITE(IU06,*) '*                                   *'
      WRITE(IU06,*) '*************************************'
      CALL FLUSH(IU06)
      CALL ABORT1
      RETURN

      END

END SUBROUTINE READPRE
