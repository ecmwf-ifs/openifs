! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      SUBROUTINE OUTCOM (IU07, IU17, IFORM)

! ----------------------------------------------------------------------

!**** *OUTCOM* - ROUTINE TO WRITE MODULES TO DISK

!     H.GUNTHER            ECMWF       04/04/1990

!     J. BIDLOT            ECMWF       10/1998.
!!!!!!                     COMMON BLOCKS HAVE BEEN CONVERTED TO MODULES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     J. BIDLOT            ECMWF       11/2003
!                          IF YOUR ARE RUNNING AT ECMWF:
!                          BE AWARE THAT IF YOU CHANGE ANYTHING TO THE
!                          STRUCTURE OF THE OUTPUT FILE YOU WILL HAVE TO
!                          MAKE SURE THAT IT IS CREATED FOR YOUR RUN, 
!                          OTHERWISE IT MIGHT PICK UP THE DEFAULT ONE
!                          THAT IS ALREADY ON DISK.
!                          YOU ALSO HAVE TO CHANGE READPRE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*    PURPOSE.
!     -------

!       TO WRITE OUT THE COMPUTED MODULES
!       (MODULE UBUF IS WRITTEN IN MUBUF)

!**   INTERFACE.
!     ----------

!       *CALL* *OUTCOM (IU07, IU17, IFORM)*
!          *IU07*   - LOGICAL UNIT FOR  UNFORMATED WRITE.
!          *IU17*   - LOGICAL UNIT FOR    FORMATED WRITE.
!          *IFORM*   - FORMAT OPTION  = 1  UNFORMATED WRITE.
!                                     = 2  FORMATED WRITE.
!                                     OTHERWISE BOTH.

!     METHOD.
!     -------

!       MODULES YOWPARAM, YOWCOUPL, YOWCURR, YOWFRED, YOWINDNL, YOWGRID,
!       YOWMAP, YOWCOUT, YOWTABL, AND YOWSHAL ARE WRITTEN TO UNIT.
!       ALL FREQUENCY AND DIRECTION DEPENDENT ARRAYS
!       ARE WRITTEN FROM 1 TO THE USED NUMBER OF FREQUENCIES (NFRE, NFRE_RED),
!       AND THE USED NUMBER OF DIRECTIONS (NANG). OTHER ARRAYS ARE
!       WRITTEN ACCORDING TO THEIR DIMENSIONS.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWGRIBHD, ONLY : IMDLGRBID_G,IMDLGRBID_M 
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED  ,             &
     &            NGX      ,NGY      ,                                  &
     &            NIBLO    ,NOVER    ,NIBL1    ,NIBLD    ,              &
     &            NIBLC    ,CLDOMAIN ,IMDLGRDID
      USE YOWCPBO  , ONLY : NBOUNC
      USE YOWFRED  , ONLY : FR       ,DFIM     ,GOM      ,C        ,    &
     &            DELTH    ,DELTR    ,TH       ,COSTH    ,SINTH
      USE YOWFPBO  , ONLY : NBOUNF
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,SINPH    ,COSPH    ,    &
     &            IJS      ,IJL
      USE YOWINDN  , ONLY : IKP      ,IKP1     ,IKM      ,IKM1     ,    &
     &            K1W      ,K2W      ,K11W     ,K21W     ,AF11     ,    &
     &            FKLAP    ,FKLAP1   ,FKLAM    ,FKLAM1   ,ACL1     ,    &
     &            ACL2     ,CL11     ,CL21     ,DAL1     ,DAL2     ,    &
     &            FRH      ,KFRH     ,MFRSTLW  ,MLSTHG
      USE YOWMAP   , ONLY : BLK2GLO  ,NX       ,NY       ,    &
     &            IPER     ,IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,    &
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,NLONRGG
      USE YOWCOUT  , ONLY : NGOUT    ,IJAR
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTH_INPUT,DEPTHA   ,DEPTHD   ,  &
     &            TCGOND   ,TFAK     ,TSIHKD   ,TFAC_ST
      USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,IUSTAR   ,IALPHA   ,    &
     &            FAC0     ,FAC1     ,FAC2     ,FAC3     ,              &
     &            FAK      ,FRHF     ,DFIMHF   ,                        &
     &            MR       ,XMR      ,MA       ,XMA      ,NFREH    ,    &
     &            NANGH    ,NMAX     ,OMEGA    ,DFDTH    ,THH      ,    &
     &            DELTHH   ,IM_P     ,IM_M     ,TA       ,TB       ,    &
     &            TC_QL    ,TT_4M    ,TT_4P    ,TFAKH
      USE YOWUNPOOL ,ONLY : LLUNSTR  ,LPREPROC

      USE UNWAM , ONLY : UNWAM_OUT
! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "outnam.intfb.h"
 
      INTEGER(KIND=JWIM), INTENT(IN) :: IU07, IU17, IFORM
      INTEGER(KIND=JWIM) :: IDUM, K, M, L
      INTEGER(KIND=JWIM) :: NBINP, NOUTT
      INTEGER(KIND=JWIM) :: NKIND !Precision used when writing

! ----------------------------------------------------------------------

  995 FORMAT(5I8)
  996 FORMAT(I8,1X,2E16.7)
  997 FORMAT(14I8,1X,A1)
  998 FORMAT(10I8)
  999 FORMAT(5E16.7)

! ----------------------------------------------------------------------
!        WRITE IDENTIFIERS (MAKE SURE TO UPDATE THE VALUES IF YOU CHANGE
!        ANYTHING TO THE MODEL). 
      
      NKIND = KIND(DELPHI)

      IF (IFORM /= 2) THEN
        WRITE(IU07) NKIND, IMDLGRDID, IMDLGRBID_G, IMDLGRBID_M
      ENDIF
      IF (IFORM /= 1) THEN
        WRITE(IU17,998) NKIND, IMDLGRDID, IMDLGRBID_G, IMDLGRBID_M
      ENDIF
!*    0. WRITE YOWPARAM (BLOCK SIZES).
!        ----------------------------

!     IDUM REPLACES IREFRA WHICH IS NO LONGER USED IN PREPROC.
      IDUM=0

      IF (IFORM /= 2) THEN
        WRITE(IU07) NANG, NFRE, NFRE_RED, NGX, NGY, NIBLO, NOVER,       &
     &              KFRH, MFRSTLW, MLSTHG,                              &
     &              NIBL1, IDUM, NIBLD, NIBLC, CLDOMAIN
      ENDIF
      IF (IFORM /= 1) THEN
        WRITE(IU17,997) NANG, NFRE, NFRE_RED, NGX, NGY, NIBLO, NOVER,   &
     &              KFRH, MFRSTLW, MLSTHG,                              &
     &              NIBL1, IDUM, NIBLD, NIBLC, CLDOMAIN
      ENDIF


!*    1. WRITE MODULE YOWFRED 
!        --------------------

      IF (IFORM /= 2) THEN
        WRITE (IU07) (FR(M),M=1,NFRE), (DFIM(M),M=1,NFRE),              &
     &   (GOM(M),M=1,NFRE), (C(M),M=1,NFRE),                            &
     &   DELTH, DELTR, (TH(K),K=1,NANG),                                &
     &   (COSTH(K),K=1,NANG), (SINTH(K),K=1,NANG)
      ENDIF
      IF (IFORM /= 1) THEN
        WRITE (IU17,999) (FR(M),M=1,NFRE), (DFIM(M),M=1,NFRE),          &
     &   (GOM(M),M=1,NFRE), (C(M),M=1,NFRE),                            &
     &   DELTH, DELTR, (TH(K),K=1,NANG),                                &
     &   (COSTH(K),K=1,NANG), (SINTH(K),K=1,NANG)
      ENDIF

! ----------------------------------------------------------------------

!*    2. WRITE MODULE YOWGRID.
!        ---------------------

      IF (IFORM /= 2) THEN
        WRITE (IU07) DELPHI, (DELLAM(L),L=1,NY), (NLONRGG(L),L=1,NY),   &
     &   (SINPH(L),L=1,NY), (COSPH(L),L=1,NY),                          &
     &   IJS, IJL
      ENDIF
      IF (IFORM /= 1) THEN
        WRITE (IU17,999) DELPHI,(DELLAM(L),L=1,NY),(REAL(NLONRGG(L),JWRB),L=1,NY), &
     &   (SINPH(L),L=1,NY), (COSPH(L),L=1,NY)
        WRITE (IU17,998) IJS, IJL
      ENDIF

! ----------------------------------------------------------------------

!*    3. WRITE MODULE YOWMAP.
!        --------------------

      IF (IFORM /= 2) THEN
        WRITE (IU07) BLK2GLO%IXLG, BLK2GLO%KXLT, NX, NY, IPER,          &
     &   AMOWEP, AMOSOP, AMOEAP, AMONOP, XDELLA, XDELLO,                &
     &   ZDELLO, IRGG
      ENDIF
      IF (IFORM /= 1) THEN
        WRITE (IU17,998) BLK2GLO%IXLG, BLK2GLO%KXLT, NX, NY, IPER
        WRITE (IU17,999) AMOWEP, AMOSOP, AMOEAP, AMONOP,                &
     &   XDELLA, XDELLO,ZDELLO
      ENDIF

! ----------------------------------------------------------------------

!*    4. WRITE MODULE YOWINDNL.
!        ---------------------

      IF (IFORM /= 2) THEN
        WRITE(IU07)(IKP(M),M=MFRSTLW,MLSTHG),                           &
     &   (IKP1(M),M=MFRSTLW,MLSTHG),                                    &
     &   (IKM(M),M=MFRSTLW,MLSTHG), (IKM1(M),M=MFRSTLW,MLSTHG),         &
     &   ((K1W(K,L),K=1,NANG),L=1,2),                                   &
     &   ((K2W(K,L),K=1,NANG),L=1,2),                                   &
     &   ((K11W(K,L),K=1,NANG),L=1,2),                                  &
     &   ((K21W(K,L),K=1,NANG),L=1,2),                                  &
     &   (AF11(M),M=MFRSTLW,MLSTHG), (FKLAP(M),M=MFRSTLW,MLSTHG),       &
     &   (FKLAP1(M),M=MFRSTLW,MLSTHG), (FKLAM(M),M=MFRSTLW,MLSTHG),     &
     &   (FKLAM1(M),M=MFRSTLW,MLSTHG),                                  &
     &   ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH
      ENDIF
      IF (IFORM /= 1) THEN
        WRITE(IU17,998)(IKP(M),M=MFRSTLW,MLSTHG),                       &
     &   (IKP1(M),M=MFRSTLW,MLSTHG),                                    &
     &   (IKM(M),M=MFRSTLW,MLSTHG), (IKM1(M),M=MFRSTLW,MLSTHG),         &
     &   ((K1W(K,L),K=1,NANG),L=1,2),                                   &
     &   ((K2W(K,L),K=1,NANG),L=1,2),                                   &
     &   ((K11W(K,L),K=1,NANG),L=1,2),                                  &
     &   ((K21W(K,L),K=1,NANG),L=1,2)
        WRITE(IU17,999)(AF11(M),M=MFRSTLW,MLSTHG),                      &
     &    (FKLAP(M),M=MFRSTLW,MLSTHG),                                  &
     &   (FKLAP1(M),M=MFRSTLW,MLSTHG), (FKLAM(M),M=MFRSTLW,MLSTHG),     &
     &   (FKLAM1(M),M=MFRSTLW,MLSTHG),                                  &
     &   ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH
      ENDIF

! ----------------------------------------------------------------------

!*    7. WRITE MODULE YOWCOUT.
!        ---------------------

      IF (IFORM /= 2) THEN
        WRITE (IU07)  NGOUT
        IF(NGOUT.GT.0) WRITE (IU07)  IJAR
      ENDIF
      IF (IFORM /= 1) THEN
        WRITE (IU17,998)  NGOUT
        IF(NGOUT.GT.0) WRITE (IU17,998)  IJAR
      ENDIF

! ----------------------------------------------------------------------

!*    8. WRITE MODULE YOWSHAL.
!        --------------------

      IF (IFORM /= 2) THEN
        WRITE (IU07) NDEPTH, DEPTHA, DEPTHD
        WRITE (IU07) DEPTH_INPUT,                                       &
     &   ((TCGOND(L,M),L=1,NDEPTH),M=1,NFRE),                           &
     &   ((TFAK(L,M),L=1,NDEPTH),M=1,NFRE),                             &
     &   ((TSIHKD(L,M),L=1,NDEPTH),M=1,NFRE),                           &
     &   ((TFAC_ST(L,M),L=1,NDEPTH),M=1,NFRE)
      ENDIF
      IF (IFORM /= 1) THEN
        WRITE (IU17,996) NDEPTH, DEPTHA, DEPTHD
        WRITE (IU17,999) DEPTH_INPUT,                                   &
     &   ((TCGOND(L,M),L=1,NDEPTH),M=1,NFRE),                           &
     &   ((TFAK(L,M),L=1,NDEPTH),M=1,NFRE),                             &
     &   ((TSIHKD(L,M),L=1,NDEPTH),M=1,NFRE),                           &
     &   ((TFAC_ST(L,M),L=1,NDEPTH),M=1,NFRE)
      ENDIF

! ----------------------------------------------------------------------

!*    9. WRITE MODULE YOWTABL (2ND AND 3RD PART).
!        ---------------------

      IF (IFORM /= 2) THEN
        WRITE (IU07) FAC0,FAC1,FAC2,FAC3,FAK,FRHF,DFIMHF 
      ENDIF
      IF (IFORM /= 1) THEN
        WRITE (IU17,999) FAC0,FAC1,FAC2,FAC3,FAK,FRHF,DFIMHF 
      ENDIF

      IF (IFORM /= 2) THEN
        WRITE (IU07) MR, XMR, MA, XMA, NFREH, NANGH, NMAX
 
        WRITE (IU07) OMEGA, DFDTH, THH, DELTHH, IM_P, IM_M,             &
     &               TA, TB, TC_QL, TT_4M, TT_4P, TFAKH
      ENDIF
      IF (IFORM /= 1) THEN
        WRITE (IU17,995) MR, MA, NFREH, NANGH, NMAX 
        WRITE (IU17,998) IM_P, IM_M 
        WRITE (IU17,999) XMR, XMA, OMEGA, DFDTH, THH, DELTHH,           &
     &                   TA, TB, TC_QL, TT_4M, TT_4P, TFAKH 
      ENDIF
! ----------------------------------------------------------------------

!*    10. WRITE MODULE YOWCURR.
!         ---------------------

!      THE CURRENTS ARE NO LONGER PART OF THE CONSTANT FILES.

! ----------------------------------------------------------------------

!*   11. WRITE NAMELIST PARWAM.
!        ----------------------
!     NBINP NOUTT have no meaning here since the size of the arrays
!     will be determined at input.
      NBINP=-1
      NOUTT=-1
      
 
      CALL OUTNAM                                                       &
     & (NANG, NFRE,                                                     &
     &  NGX, NGY, NIBLO, NOVER, NGOUT, NOUTT,                           &
     &  KFRH, MFRSTLW, MLSTHG,                                          &
     &  NBOUNC, NBOUNF, NBINP, NIBL1, NIBLD, NIBLC,                     &
     &  ITAUMAX, JUMAX, IUSTAR, IALPHA, NDEPTH, IDUM, IPER)

      IF (LLUNSTR .AND. LPREPROC) THEN
        CALL UNWAM_OUT(IU07)
      END IF

      END SUBROUTINE OUTCOM
