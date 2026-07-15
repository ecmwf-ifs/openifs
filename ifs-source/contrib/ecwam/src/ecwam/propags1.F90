! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE PROPAGS1 (F1, F3, NINF, NSUP, KIJS, KIJL,           &
 &                   BLK2GLO,                                  &
 &                   DEPTH_EXT,                                &
 &                   CGROUP_EXT, OMOSNH2KD_EXT,                &
 &                   DELLAM1_EXT, COSPHM1_EXT,                 &
 &                   U_EXT, V_EXT,                             &
 &                   L1STCALL)

! ----------------------------------------------------------------------

!**** *PROPAGS1* - COMPUTATION OF A PROPAGATION TIME STEP.

!     S.D. HASSELMANN.
!     OPTIMIZED BY: L. ZAMBRESKY AND H. GUENTHER

!     MODIFIED BY   H. GUNTHER   01/06/90    -   LAND POINTS ARE TAKEN
!                             OUT OF BLOCKS AND REFRACTION INTEGRATION
!                             CORRECTED FOR N-S AND S-N PROPAGATION.

!     K.P. HUBBERT                /07/89    -   DEPTH AND CURRENT
!     S. HASSELMANN   MPIFM       /04/90        REFRACTION SHALLOW

!     H. GUNTHER   GKSS/ECMWF   17/01/91    -   MODIFIED FOR CYCLE_4

!     J. BIDLOT    ECMWF   APRIL 1997       - MODIFIED ADVECTION SCHEME 
!     J. BIDLOT    ECMWF   2003             - OBSTRUCTION COEFFICIENTS  

!     J. BIDLOT    ECMWF   2004        - FIRST ORDER PROPAGATION SCHEME
!                                        BASED ON 2 ROTATED QUADRANTS

!!!!!!!!!!!!!!!!!!!!!!    WARNING    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! SO FAR ONLY MODIFIED FOR PROPAGATION WITHOUT REFRACTION ON SPHERICAL GRID !!!!
!!!!! (i.e the 3000 block)
!!!!! THE OTHERS ARE AS IN PROPAGS.
!!!!!!!!!!!!!!!!!!!!!!    WARNING    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*    PURPOSE.
!     --------

!       COMPUTATION OF A PROPAGATION TIME STEP.

!**   INTERFACE.
!     ----------

!       *CALL* *PROPAGS1(F1, F3, NINF, NSUP, KIJS, KIJL,
!                        BLK2GLO,
!                        DEPTH_EXT,
!                        CGROUP_EXT, OMOSNH2KD_EXT, 
!                        U_EXT, V_EXT, 
!                        L1STCALL)
!          *F1*          - SPECTRUM AT TIME T (with exchange halo).
!          *F3*          - SPECTRUM AT TIME T+DELT (without halo).
!          *NINF:NSUP+1* - 1st DIMENSION OF F1 and F3
!          *KIJS*        - ACTIVE INDEX OF FIRST POINT
!          *KIJL*        - ACTIVE INDEX OF LAST POINT
!          *BLK2GLO*     - BLOCK TO GRID TRANSFORMATION
!          *DEPTH_EXT*   - WATER DEPTH 
!          *CGROUP_EXT*  - GROUP VELOCITY
!          *OMOSNH2KD_EXT- OMEGA / SINH(2KD)
!          *U_EXT        - U-COMPONENT OF SURFACE CURRENT
!          *V_EXT        - V-COMPONENT OF SURFACE CURRENT
!          *L1STCALL*    - LOGICAL SHOULD BE FALSE AFTER THE FIRST CALL.

!     METHOD.
!     -------

!       FIRST ORDER FLUX SCHEME.

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO

      USE YOWCURR  , ONLY : LLCHKCFL
      USE YOWFRED  , ONLY : FR       ,GOM      ,DELTH    ,FRATIO   ,    &
     &            COSTH    ,SINTH
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,SINPH    ,              &
     &            COSPH    ,CDR      ,SDR      ,PRQRT
      USE YOWMAP   , ONLY : IRGG     ,NIBLO
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED
      USE YOWPCONS , ONLY : PI       ,ZPI      ,R
      USE YOWREFD  , ONLY : THDD     ,THDC     ,SDOT
      USE YOWSTAT  , ONLY : IDELPRO  ,ICASE    ,IREFRA
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KRLAT     ,             &
     &            KRLON    ,WLAT     ,WRLAT     ,                       &
     &            WRLON    ,OBSLAT   ,OBSLON   ,OBSRLAT  ,OBSRLON
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "checkcfl.intfb.h"

      REAL(KIND=JWRB),DIMENSION(NINF:NSUP+1, NANG, NFRE_RED), INTENT(IN) :: F1
      REAL(KIND=JWRB),DIMENSION(NINF:NSUP+1, NANG, NFRE_RED), INTENT(OUT) :: F3
      INTEGER(KIND=JWIM), INTENT(IN) :: NINF, NSUP
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN):: DEPTH_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: DELLAM1_EXT 
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: COSPHM1_EXT 
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED), INTENT(IN) :: CGROUP_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED), INTENT(IN) :: OMOSNH2KD_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: U_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: V_EXT
      LOGICAL, INTENT(IN) :: L1STCALL


      INTEGER(KIND=JWIM) :: K, M, IJ, JH
      INTEGER(KIND=JWIM) :: NLAND
      INTEGER(KIND=JWIM) :: IC, IJLA, IJPH, KP1, KM1, MP1, MM1, KRLA,   &
     &                      KRLA2
      INTEGER(KIND=JWIM) :: KRLO, KRLO2, IJLAR, IJLAR1, IJPHR, IJPHR1,  &
     &                      IJLA1, IJPH1
      INTEGER(KIND=JWIM),ALLOCATABLE,DIMENSION(:) :: IJRLA, IJRPH

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: DELPRO, DELPH0, DELTH0, DELFR0
      REAL(KIND=JWRB) :: SD, CD, SDA, CDA, DNO, DTT, SDD, CDD, DTHP,    &
     &                   DTHM, DFP, DFM
      REAL(KIND=JWRB) :: CGS, CGC, DLWE, DLEA, DPSO, DPNO, SDA2, SP,    &
     &                   SM, TANPH  
      REAL(KIND=JWRB) :: DPNO2, DPSO2, XX, YY, XR, YR, FF
      REAL(KIND=JWRB) :: SQRT2, COFDELPH0R, CGOMFULL

      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: DELLA0,DCO,DP1,DP2
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: DCOM1
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: DPN,DPN2,DPS,DPS2
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: DEA,DLE,DLW,DPH,DLA
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: RLE,RLE2,RPN,RPN2
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: DTP,DTM,DTC,DRGP,DRGM
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: DOP,DOM,DRCP,DRCM
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: DRDP,DRDM
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: DLADCO
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: DELPH0_CDA
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: DELLA0R,DELPH0R
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: SDDL
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:,:) :: WLATM1, WRLATM1,    &
     &                                              WRLONM1
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:,:,:) :: CGKLON, CGKLAT
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: CFLEA, CFLWE ,CFLNO,  &
     &                                            CFLSO,CFLNO2 
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: CFLSO2, CFLTP, CFLTM, &
     &                                            CFLOP, CFLOM
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:,:) :: ACGKRLON, ACGKRLAT
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: CFLEAR, CFLWER
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: CFLNOR, CFLSOR,CFLNO2R


! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PROPAGS1',0,ZHOOK_HANDLE)

      ALLOCATE(DELLA0(NINF:NSUP+1))
      ALLOCATE(DPH(NINF:NSUP+1))
      ALLOCATE(DLA(NINF:NSUP+1))
      ALLOCATE(DTP(KIJS:KIJL))
      ALLOCATE(DTM(KIJS:KIJL))
      ALLOCATE(DTC(KIJS:KIJL))

      NLAND=NSUP+1

      DELPRO = REAL(IDELPRO)   

!*    0.2 SPHERICAL OR CARTESIAN GRID?
!         ----------------------------

      IF (ICASE == 1) THEN

!*    0.2.1 SPHERICAL GRID.
!           ---------------
        ALLOCATE(DCO(NINF:NSUP+1))
        ALLOCATE(DCOM1(NINF:NSUP+1))
        ALLOCATE(DP1(KIJS:KIJL))
        ALLOCATE(DP2(KIJS:KIJL))

!*    0.2.1.1 COSINE OF LATITUDE.
!             -------------------

        DELLA0(NLAND) = 0.0_JWRB
        DO IJ = NINF,NSUP
          DCO(IJ) = COSPHM1_EXT(IJ)
          DELLA0(IJ) = DELPRO*DELLAM1_EXT(IJ)
        ENDDO
        DO IJ = NINF,NSUP
          DCOM1(IJ) = 1.0_JWRB/COSPHM1_EXT(IJ)
        ENDDO

!*    0.2.1.2 COMPUTE COS PHI FACTOR FOR ADJOINING GRID POINT.
!             ------------------------------------------------

        DO IJ = KIJS,KIJL
          JH = KLAT(IJ,1,1)
          IF (JH == NLAND) THEN
            DP1(IJ) = 1.0_JWRB
          ELSE
            DP1(IJ) = DCO(IJ)*DCOM1(JH)
          ENDIF
          JH = KLAT(IJ,2,1)
          IF (JH == NLAND) THEN
            DP2(IJ) = 1.0_JWRB
          ELSE
            DP2(IJ) = DCO(IJ)*DCOM1(JH)
          ENDIF
        ENDDO

        IF (IREFRA /= 2 .AND. IREFRA /= 3 ) THEN

!*       BRANCH TO 3. IF WITHOUT REFRACTION OR DEPTH.
!        --------------------------------------------

          GOTO 3000
        ELSE

!*       BRANCH TO 4. IF DEPTH AND CURRENT REFRACTION.
!        ---------------------------------------------

          GOTO 4000
        ENDIF
      ELSE

!*    0.2.2 CARTESIAN GRID.
!           ---------------

        DELLA0(NLAND) = 0.0_JWRB
        DO IJ = NINF,NSUP
          DELLA0(IJ) = DELPRO*DELLAM1_EXT(IJ)
        ENDDO

!*    0.2.2.1 BRANCH TO 2. IF DEPTH AND CURRENT REFRACTION.
!             ---------------------------------------------

        IF (IREFRA == 2 .OR. IREFRA == 3 ) GOTO 2000

      ENDIF

! ----------------------------------------------------------------------

!*    1. PROPAGATION FOR CARTESIAN GRID
!*       WITHOUT REFRACTION OR DEPTH REFRATION.
!        --------------------------------------

      DELPH0 = DELPRO/DELPHI
      DELTH0 = 0.25_JWRB*DELPRO/DELTH

      IF (IREFRA == 1) THEN
        ALLOCATE(DRDP(KIJS:KIJL))
        ALLOCATE(DRDM(KIJS:KIJL))
      ENDIF

!*    1.1 LOOP OVER DIRECTIONS.
!         ---------------------

      DO K=1,NANG
        SD = SINTH(K)
        CD = COSTH(K)

!*    1.1.1 INDEX FOR ADJOINING POINTS.
!           ---------------------------

        IF (SD < 0) THEN
          IJLA = 2
        ELSE
          IJLA = 1
        ENDIF
        IF (CD < 0) THEN
          IJPH = 2
        ELSE
          IJPH = 1
        ENDIF

          SD = 0.5_JWRB*SD
          CD = 0.5_JWRB*CD

!*    1.1.3.1 DEPTH REFRACTION.
!             -----------------

          IF (IREFRA == 1) THEN
            KP1 = K+1
            IF (KP1 > NANG) KP1 = 1
            KM1 = K-1
            IF (KM1 < 1) KM1 = NANG
            DO IJ = KIJS,KIJL
              DRDP(IJ) = (THDD(IJ,K) + THDD(IJ,KP1))*DELTH0
              DRDM(IJ) = (THDD(IJ,K) + THDD(IJ,KM1))*DELTH0
            ENDDO
          ENDIF

!*    1.1.3.2 LOOP OVER FREQUENCIES.
!             ----------------------

          DO M=1,NFRE_RED

!*    1.1.3.2.2 WEIGHTS IN INTEGRATION SCHEME.
!               ------------------------------

            IF (SD >= 0.0_JWRB) THEN
              DO IJ=KIJS,KIJL
                SDD = SD*DELLA0(IJ)
                DLA(IJ) = SDD*(CGROUP_EXT(KLON(IJ,1),M) + CGROUP_EXT(IJ,M))
                DTC(IJ) = SDD*(CGROUP_EXT(KLON(IJ,2),M) + CGROUP_EXT(IJ,M))
              ENDDO
            ELSE
              DO IJ=KIJS,KIJL
                SDD = SD*DELLA0(IJ)
                DLA(IJ) =-SDD*(CGROUP_EXT(KLON(IJ,2),M) + CGROUP_EXT(IJ,M))
                DTC(IJ) =-SDD*(CGROUP_EXT(KLON(IJ,1),M) + CGROUP_EXT(IJ,M))
              ENDDO
            ENDIF

            IF (CD >= 0.0_JWRB) THEN
              DO IJ=KIJS,KIJL
                CDD = CD*DELPH0
                DPH(IJ) = CDD*(CGROUP_EXT(KLAT(IJ,1,1),M) + CGROUP_EXT(IJ,M))
                DTC(IJ) = DTC(IJ)                                       &
     &           + CDD*(CGROUP_EXT(KLAT(IJ,2,1),M) + CGROUP_EXT(IJ,M))
              ENDDO
            ELSE
              DO IJ=KIJS,KIJL
                CDD = CD*DELPH0
                DPH(IJ) =-CDD*(CGROUP_EXT(KLAT(IJ,2,1),M) + CGROUP_EXT(IJ,M))
                DTC(IJ) = DTC(IJ)                                       &
     &           -CDD*(CGROUP_EXT(KLAT(IJ,1,1),M) + CGROUP_EXT(IJ,M))
              ENDDO
            ENDIF
            IF (IREFRA == 1) THEN
              DO IJ = KIJS,KIJL
                DTHP = OMOSNH2KD_EXT(IJ,M)*DRDP(IJ)
                DTHM = OMOSNH2KD_EXT(IJ,M)*DRDM(IJ)
                DTC(IJ) = DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
                DTP(IJ) = -DTHP+ABS(DTHP)
                DTM(IJ) =  DTHM+ABS(DTHM)
              ENDDO
            ENDIF

!*    1.1.3.2.3 LOOP OVER GRIDPOINTS.
!               ---------------------

            DO IJ = KIJS,KIJL
              F3(IJ,K,M) = (1.0_JWRB-DTC(IJ))*F1(IJ,K,M )               &
     &         + DPH(IJ) * F1(KLAT(IJ,IJPH,1),K  ,M)                    &
     &         + DLA(IJ) * F1(KLON(IJ,IJLA),K  ,M)
            ENDDO
            IF (IREFRA == 1) THEN
              DO IJ = KIJS,KIJL
                F3(IJ,K,M) = F3(IJ,K,M )                                &
     &           + DTP(IJ) * F1(IJ,KP1,M)                               &
     &           + DTM(IJ) * F1(IJ,KM1,M)
              ENDDO
            ENDIF

!*    BRANCH BACK TO 1.1.3.2 FOR NEXT FREQUENCY.

          ENDDO

!*    BRANCH BACK TO 1.1 FOR NEXT DIRECTION.

      ENDDO

!*    1.2 END OF PROPAGATION FOR CARTESIAN GRID
!*        WITHOUT REFRACTION OR DEPTH REFRACTION, RETURN.
!         -----------------------------------------------

      DEALLOCATE(DELLA0)
      DEALLOCATE(DPH,DLA)
      DEALLOCATE(DTP,DTM,DTC)
      IF (IREFRA == 1) THEN
        DEALLOCATE(DRDP)
        DEALLOCATE(DRDM)
      ENDIF

      IF (LHOOK) CALL DR_HOOK('PROPAGS1',1,ZHOOK_HANDLE)

      RETURN

! ----------------------------------------------------------------------

!*    2. PROPAGATION FOR CARTESIAN GRID
!*       WITH DEPTH AND CURRENT REFRACTION.
!        ----------------------------------

 2000 CONTINUE

      ALLOCATE(DPN(KIJS:KIJL))
      ALLOCATE(DPS(KIJS:KIJL))
      ALLOCATE(DLE(KIJS:KIJL))
      ALLOCATE(DLW(KIJS:KIJL))

      ALLOCATE(DOP(KIJS:KIJL))
      ALLOCATE(DOM(KIJS:KIJL))
      ALLOCATE(DRCP(KIJS:KIJL))
      ALLOCATE(DRCM(KIJS:KIJL))
      ALLOCATE(DRDP(KIJS:KIJL))
      ALLOCATE(DRDM(KIJS:KIJL))

      DELPH0 = 0.25_JWRB*DELPRO/DELPHI
      DELTH0 = 0.25_JWRB*DELPRO/DELTH
      DELLA0 = 0.25_JWRB*DELLA0
      DELFR0 = 0.25_JWRB*DELPRO/((FRATIO-1.0_JWRB)*ZPI)

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

      DO K=1,NANG
        KP1 = K+1
        IF (KP1 > NANG) KP1 = 1
        KM1 = K-1
        IF (KM1 < 1) KM1 = NANG
        SD = SINTH(K)
        CD = COSTH(K)

!*    2.1.1 DEPTH REFRACTION IF SHALLOW WATER.
!           ----------------------------------

        DO IJ = KIJS,KIJL
          DRDP(IJ) = (THDD(IJ,K) + THDD(IJ,KP1))*DELTH0
          DRDM(IJ) = (THDD(IJ,K) + THDD(IJ,KM1))*DELTH0
        ENDDO

!*    2.1.2 CURRENT REFRACTION.
!           -------------------

        DO IJ = KIJS,KIJL
          DRCP(IJ) = (THDC(IJ,K) + THDC(IJ,KP1))*DELTH0
          DRCM(IJ) = (THDC(IJ,K) + THDC(IJ,KM1))*DELTH0
        ENDDO

!*    2.1.3 LOOP OVER FREQUENCIES.
!           ----------------------

        DO M=1,NFRE_RED

            MP1 = MIN(NFRE_RED,M+1)
            MM1 = MAX(1,M-1)
            DFP = DELFR0/FR(M)
            DFM = DELFR0/FR(MM1)


!*    2.1.3.2.2 WEIGHTS IN INTEGRATION SCHEME.
!               ------------------------------

            DLA(NLAND) = SD*CGROUP_EXT(NLAND,M)*DELLA0(NINF)
            DPH(NLAND) = CD*CGROUP_EXT(NLAND,M)*DELPH0
            DO IJ=NINF,NSUP
              DLA(IJ) = (U_EXT(IJ) + SD*CGROUP_EXT(IJ,M))*DELLA0(IJ)
              DPH(IJ) = (V_EXT(IJ) + CD*CGROUP_EXT(IJ,M))*DELPH0
            ENDDO
            DO IJ=KIJS,KIJL
              DLWE = DLA(IJ) + DLA(KLON(IJ,1))
              DLEA = DLA(IJ) + DLA(KLON(IJ,2))
              DLE(IJ) = -DLEA+ABS(DLEA)
              DLW(IJ) =  DLWE+ABS(DLWE)
              DTC(IJ) = DLEA+ABS(DLEA)-DLWE+ABS(DLWE)

              DPSO = DPH(IJ) + DPH(KLAT(IJ,1,1))
              DPNO = DPH(IJ) + DPH(KLAT(IJ,2,1))
              DPN(IJ) = -DPNO+ABS(DPNO)
              DPS(IJ) =  DPSO+ABS(DPSO)
              DTC(IJ) = DTC(IJ) + DPNO+ABS(DPNO)-DPSO+ABS(DPSO)

              DTHP = OMOSNH2KD_EXT(IJ,M)*DRDP(IJ) + DRCP(IJ)
              DTHM = OMOSNH2KD_EXT(IJ,M)*DRDM(IJ) + DRCM(IJ)
              DTC(IJ) = DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
              DTP(IJ) = -DTHP+ABS(DTHP)
              DTM(IJ) =  DTHM+ABS(DTHM)

              DTHP = (SDOT(IJ,K,M) + SDOT(IJ,K,MP1))*DFP
              DTHM = (SDOT(IJ,K,M) + SDOT(IJ,K,MM1))*DFM
              DTC(IJ) =  DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
              DOP(IJ) = (-DTHP+ABS(DTHP))/FRATIO
              DOM(IJ) = ( DTHM+ABS(DTHM))*FRATIO
            ENDDO

!*    2.1.3.3 LOOP OVER GRIDPOINTS.
!             ---------------------

          DO IJ = KIJS,KIJL
            F3(IJ,K,M) = (1.0_JWRB-DTC(IJ))*F1(IJ,K,M )                 &
     &       + DPN(IJ) * F1(KLAT(IJ,2,1),K  ,M)                         &
     &       + DPS(IJ) * F1(KLAT(IJ,1,1),K  ,M)                         &
     &       + DLE(IJ) * F1(KLON(IJ,2),K  ,M)                           &
     &       + DLW(IJ) * F1(KLON(IJ,1),K  ,M)                           &
     &       + DTP(IJ) * F1(IJ        ,KP1,M)                           &
     &       + DTM(IJ) * F1(IJ        ,KM1,M)                           &
     &       + DOP(IJ) * F1(IJ        ,K  ,MP1)                         &
     &       + DOM(IJ) * F1(IJ        ,K  ,MM1)
          ENDDO

!*    BRANCH BACK TO 2.1.3 FOR NEXT FREQUENCY.

        ENDDO

!*    BRANCH BACK TO 2.1 FOR NEXT DIRECTION.

      ENDDO

!*    2.2 END OF PROPAGATION FOR CARTESIAN GRID
!*        WITH DEPTH AND CURRENT REFRACTION, RETURN.
!         ------------------------------------------

      DEALLOCATE(DELLA0,DPN,DPS)
      DEALLOCATE(DLE,DLW,DPH,DLA)
      DEALLOCATE(DTP,DTM,DTC)

      DEALLOCATE(DOP,DOM,DRCP,DRCM)
      DEALLOCATE(DRDP)
      DEALLOCATE(DRDM)

      IF (LHOOK) CALL DR_HOOK('PROPAGS1',1,ZHOOK_HANDLE)

      RETURN

! ----------------------------------------------------------------------

!*    3. PROPAGATION FOR SPHERICAL LATITUDE/LONGITUDE GRID
!*       WITHOUT CURRENT OR DEPTH REFRACTION.
!        -------------------------------------------------

 3000 CONTINUE

      ALLOCATE(DPN(KIJS:KIJL))
      ALLOCATE(DPN2(KIJS:KIJL))
      ALLOCATE(DLE(KIJS:KIJL))
      ALLOCATE(RPN(KIJS:KIJL))
      ALLOCATE(RPN2(KIJS:KIJL))
      ALLOCATE(RLE(KIJS:KIJL))
      ALLOCATE(RLE2(KIJS:KIJL))
      ALLOCATE(DRGP(KIJS:KIJL))
      ALLOCATE(DRGM(KIJS:KIJL))
      ALLOCATE(DLADCO(KIJS:KIJL))
      ALLOCATE(IJRLA(KIJS:KIJL))
      ALLOCATE(IJRPH(KIJS:KIJL))
      ALLOCATE(DELPH0_CDA(KIJS:KIJL))
      ALLOCATE(DELLA0R(KIJS:KIJL))
      ALLOCATE(DELPH0R(KIJS:KIJL))
      ALLOCATE(SDDL(KIJS:KIJL))
      ALLOCATE(CFLEA(KIJS:KIJL))
      ALLOCATE(CFLNO(KIJS:KIJL))
      ALLOCATE(CFLTP(KIJS:KIJL))
      ALLOCATE(CFLTM(KIJS:KIJL))
      ALLOCATE(CFLEAR(KIJS:KIJL))
      ALLOCATE(CFLNOR(KIJS:KIJL))

      IF (IREFRA == 1) THEN
        ALLOCATE(DRDP(KIJS:KIJL))
        ALLOCATE(DRDM(KIJS:KIJL))
      ENDIF

      DELTH0 = 0.25_JWRB*DELPRO/DELTH
      DELPH0 = 0.5_JWRB*DELPRO/DELPHI

      SQRT2=2.0_JWRB*SIN(0.25_JWRB*PI)
!     in rotated grid there is always an average of the Cg's
      COFDELPH0R = 0.5_JWRB*DELPRO/(SQRT2*DELPHI)

      DO IJ=KIJS,KIJL
        DELPH0R(IJ) = PRQRT(IJ)*COFDELPH0R
!       by construct DELLA0R=DELPH0R
        DELLA0R(IJ)=DELPH0R(IJ)
      ENDDO

      DO IJ=KIJS,KIJL
        DLADCO(IJ) = (1.0_JWRB-PRQRT(IJ))*(DCO(IJ)*DELLA0(IJ))
      ENDDO

      ALLOCATE(WLATM1(KIJS:KIJL,2))
      ALLOCATE(WRLATM1(KIJS:KIJL,2))
      ALLOCATE(WRLONM1(KIJS:KIJL,2))
      DO IC=1,2
        DO IJ = KIJS,KIJL
          WLATM1(IJ,IC) = 1.0_JWRB - WLAT(IJ,IC)
          WRLATM1(IJ,IC) = 1.0_JWRB - WRLAT(IJ,IC)
          WRLONM1(IJ,IC) = 1.0_JWRB - WRLON(IJ,IC)
        ENDDO
      ENDDO

!     ACGKRLON AND ACGKRLAT WILL RECOMPUTED FOR EACH DIRECTION AND FREQUENCY
      ALLOCATE(ACGKRLON(KIJS:KIJL,2))
      ALLOCATE(ACGKRLAT(KIJS:KIJL,2))

        ALLOCATE(CGKLON(KIJS:KIJL,NFRE_RED,2))
        ALLOCATE(CGKLAT(KIJS:KIJL,NFRE_RED,2))

        DO IC=1,2
          DO M=1,NFRE_RED
            DO IJ=KIJS,KIJL
              IF (KLON(IJ,IC) /= NLAND) THEN
                CGKLON(IJ,M,IC) = CGROUP_EXT(KLON(IJ,IC),M) + CGROUP_EXT(IJ,M)
              ELSE
                CGKLON(IJ,M,IC) = 2.0_JWRB*CGROUP_EXT(IJ,M)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        IC=1
          DO M=1,NFRE_RED
            DO IJ=KIJS,KIJL
              IF (KLAT(IJ,IC,1) == NLAND .AND.                          &
     &           KLAT(IJ,IC,2) == NLAND) THEN
                CGKLAT(IJ,M,IC) = CGROUP_EXT(IJ,M)*(DP1(IJ)+1.0_JWRB)
              ELSE IF (KLAT(IJ,IC,1) == NLAND) THEN
                CGKLAT(IJ,M,IC) = CGROUP_EXT(IJ,M) + DP1(IJ)*               &
     &                          (WLAT(IJ,IC)*CGROUP_EXT(KLAT(IJ,IC,2),M) +  &
     &                           WLATM1(IJ,IC)*CGROUP_EXT(KLAT(IJ,IC,2),M))
              ELSE IF (KLAT(IJ,IC,2) == NLAND) THEN
                CGKLAT(IJ,M,IC) = CGROUP_EXT(IJ,M) + DP1(IJ)*               &
     &                          (WLAT(IJ,IC)*CGROUP_EXT(KLAT(IJ,IC,1),M) +  &
     &                           WLATM1(IJ,IC)*CGROUP_EXT(KLAT(IJ,IC,1),M))
              ELSE
                CGKLAT(IJ,M,IC) = CGROUP_EXT(IJ,M) + DP1(IJ)*               &
     &                          (WLAT(IJ,IC)*CGROUP_EXT(KLAT(IJ,IC,1),M) +  &
     &                           WLATM1(IJ,IC)*CGROUP_EXT(KLAT(IJ,IC,2),M))
              ENDIF
            ENDDO
          ENDDO
        IC=2
          DO M=1,NFRE_RED
            DO IJ=KIJS,KIJL
              IF (KLAT(IJ,IC,1) == NLAND .AND.                          &
     &           KLAT(IJ,IC,2) == NLAND) THEN
                CGKLAT(IJ,M,IC) = CGROUP_EXT(IJ,M)*(DP2(IJ)+1.0_JWRB)
              ELSE IF (KLAT(IJ,IC,1) == NLAND) THEN
                CGKLAT(IJ,M,IC) = CGROUP_EXT(IJ,M) + DP2(IJ)*               &
     &                          (WLAT(IJ,IC)*CGROUP_EXT(KLAT(IJ,IC,2),M) +  &
     &                           WLATM1(IJ,IC)*CGROUP_EXT(KLAT(IJ,IC,2),M))
              ELSE IF (KLAT(IJ,IC,2) == NLAND) THEN
                CGKLAT(IJ,M,IC) = CGROUP_EXT(IJ,M) + DP2(IJ)*               &
     &                          (WLAT(IJ,IC)*CGROUP_EXT(KLAT(IJ,IC,1),M) +  &
     &                           WLATM1(IJ,IC)*CGROUP_EXT(KLAT(IJ,IC,1),M))
              ELSE
                CGKLAT(IJ,M,IC) = CGROUP_EXT(IJ,M) + DP2(IJ)*               &
     &                         (WLAT(IJ,IC)*CGROUP_EXT(KLAT(IJ,IC,1),M) +   &
     &                          WLATM1(IJ,IC)*CGROUP_EXT(KLAT(IJ,IC,2),M))
              ENDIF
            ENDDO
          ENDDO
!!!!!!!!!
        IF (L1STCALL) THEN
!         THE OBSTRUCTION COEFFICIENTS ARE RESET TO BE THE PRODUCT
!         OF THE OBSTRUCTION BY THE GROUP VELOCITY AT MIDPOINT
!!! only for the non rotated grid
          DO IC=1,2
            DO M=1,NFRE_RED
              DO IJ=KIJS,KIJL
                OBSLON(IJ,M,IC) = OBSLON(IJ,M,IC) * CGKLON(IJ,M,IC)
                OBSLAT(IJ,M,IC) = OBSLAT(IJ,M,IC) * CGKLAT(IJ,M,IC)
              ENDDO
            ENDDO
          ENDDO
        ENDIF


!*    3.1 LOOP OVER DIRECTIONS.
!         ---------------------

      DO K=1,NANG
        KP1 = K+1
        IF (KP1 > NANG) KP1 = 1
        KM1 = K-1
        IF (KM1 < 1) KM1 = NANG

        SD = SINTH(K)
        CD = COSTH(K)
        IF (SD < 0) THEN
          IJLA = 2
        ELSE
          IJLA = 1
        ENDIF
        IF (CD < 0) THEN
          IJPH = 2
        ELSE
          IJPH = 1
        ENDIF
        SDA = ABS(SD)
        SDA2 = 0.5_JWRB*SDA
        CDA = ABS(CD)

        DO IJ=KIJS,KIJL
          DELPH0_CDA(IJ) = (1.0_JWRB-PRQRT(IJ))*DELPH0*CDA
          SDDL(IJ)=SDA2*DLADCO(IJ)
        ENDDO

        DO IJ = KIJS,KIJL
          IF (SDR(IJ,K) < 0) THEN
            IJRLA(IJ) = 2
          ELSE
            IJRLA(IJ) = 1
          ENDIF
          IF (CDR(IJ,K) < 0) THEN
            IJRPH(IJ) = 2
          ELSE
            IJRPH(IJ) = 1
          ENDIF
        ENDDO

!*    3.1.1 COMPUTE GRID REFRACTION.
!           ------------------------

        SP  = DELTH0*(SINTH(K)+SINTH(KP1))/R
        SM  = DELTH0*(SINTH(K)+SINTH(KM1))/R
        DO IJ = KIJS,KIJL
          JH = BLK2GLO%KXLT(IJ)
          TANPH = SINPH(JH)/COSPH(JH)
          DRGP(IJ) = TANPH*SP
          DRGM(IJ) = TANPH*SM
        ENDDO

!*    3.1.4.1 COMPUTE DEPTH REFRACTION.
!             -------------------------

          IF (IREFRA == 1) THEN
            DO IJ = KIJS,KIJL
              DRDP(IJ) = (THDD(IJ,K) + THDD(IJ,KP1))*DELTH0
              DRDM(IJ) = (THDD(IJ,K) + THDD(IJ,KM1))*DELTH0
            ENDDO
          ENDIF

!*    3.1.4.2 LOOP OVER FREQUENCIES.
!             ----------------------

          DO M=1,NFRE_RED


!*    3.1.4.3.2 LAT / LONG WEIGHTS IN INTEGRATION SCHEME.
!               -----------------------------------------

            IJLA1=3-IJLA
            DO IJ=KIJS,KIJL
              CFLEA(IJ) = SDDL(IJ)*CGKLON(IJ,M,IJLA1)
              DTC(IJ) = CFLEA(IJ)  
              DLE(IJ) = SDDL(IJ)*OBSLON(IJ,M,IJLA)
            ENDDO

            IJPH1=3-IJPH
            DO IJ=KIJS,KIJL
              CFLNO(IJ) = DELPH0_CDA(IJ)*CGKLAT(IJ,M,IJPH1)
              DTC(IJ) = DTC(IJ) + CFLNO(IJ)
              XX = DELPH0_CDA(IJ)*OBSLAT(IJ,M,IJPH)
              DPN(IJ) = WLAT(IJ,IJPH)*XX
              DPN2(IJ) = WLATM1(IJ,IJPH)*XX
            ENDDO


!           for rotated grid the advection velocities are function of 
!           frequency and direction
            DO IC=1,2
              DO IJ=KIJS,KIJL
                IF (KRLAT(IJ,IC,1) /= NLAND) THEN
                  KRLA=KRLAT(IJ,IC,1)
                ELSE
                  KRLA=IJ
                ENDIF
                IF (KRLAT(IJ,IC,2) /= NLAND) THEN
                  KRLA2=KRLAT(IJ,IC,2)
                ELSE
                  KRLA2=IJ
                ENDIF
                IF (KRLON(IJ,IC,1) /= NLAND) THEN
                  KRLO=KRLON(IJ,IC,1)
                ELSE
                  KRLO=IJ
                ENDIF
                IF (KRLON(IJ,IC,2) /= NLAND) THEN
                  KRLO2=KRLON(IJ,IC,2)
                ELSE
                  KRLO2=IJ
                ENDIF
                ACGKRLAT(IJ,IC) = CGROUP_EXT(IJ,M)*CDR(IJ,K) +             &
     &             WRLAT(IJ,IC)*CGROUP_EXT(KRLA,M)*CDR(KRLA,K) +           &
     &             WRLATM1(IJ,IC)*CGROUP_EXT(KRLA2,M)*CDR(KRLA2,K)
                ACGKRLON(IJ,IC) =  CGROUP_EXT(IJ,M)*SDR(IJ,K) +            &
     &             WRLON(IJ,IC)*CGROUP_EXT(KRLO,M)*SDR(KRLO,K) +           &
     &             WRLONM1(IJ,IC)*CGROUP_EXT(KRLO2,M)*SDR(KRLO2,K)
                ACGKRLAT(IJ,IC) = ABS(ACGKRLAT(IJ,IC)) 
                ACGKRLON(IJ,IC) = ABS(ACGKRLON(IJ,IC))
              ENDDO
            ENDDO

            DO IJ=KIJS,KIJL
              IJLAR=IJRLA(IJ)
              IJLAR1=3-IJLAR
              XX= DCO(IJ)*DELLA0R(IJ)
              CFLEAR(IJ)=XX*ACGKRLON(IJ,IJLAR1)
              DTC(IJ) = DTC(IJ) + CFLEAR(IJ) 
              FF=XX*ACGKRLON(IJ,IJLAR)
              FF = FF*OBSRLON(IJ,M,IJLAR)
              RLE(IJ) = WRLON(IJ,IJLAR)*FF
              RLE2(IJ) = WRLONM1(IJ,IJLAR)*FF
            ENDDO

            DO IJ=KIJS,KIJL
              IJPHR=IJRPH(IJ)
              IJPHR1=3-IJPHR
              XX= DCO(IJ)*DELPH0R(IJ)
              CFLNOR(IJ)=XX*ACGKRLAT(IJ,IJPHR1)
              DTC(IJ) = DTC(IJ) + CFLNOR(IJ) 
              FF=XX*ACGKRLAT(IJ,IJPHR)
              FF= FF*OBSRLAT(IJ,M,IJPHR)
              RPN(IJ) = WRLAT(IJ,IJPHR)*FF
              RPN2(IJ) = WRLATM1(IJ,IJPHR)*FF
            ENDDO

!*    3.1.4.2.3 REFRACTION WEIGHTS IN INTEGRATION SCHEME.
!               -----------------------------------------

            IF (IREFRA == 0) THEN
              DO IJ=KIJS,KIJL
                DTHP = DRGP(IJ)*CGROUP_EXT(IJ,M)
                DTHM = DRGM(IJ)*CGROUP_EXT(IJ,M)
                CFLTP(IJ) = DTHP+ABS(DTHP)
                CFLTM(IJ) = -DTHM+ABS(DTHM) 
                DTC(IJ) =  DTC(IJ) + CFLTP(IJ) + CFLTM(IJ)
                DTP(IJ) = -DTHP+ABS(DTHP)
                DTM(IJ) =  DTHM+ABS(DTHM)
              ENDDO
            ELSE
              DO IJ=KIJS,KIJL
                DTHP = DRGP(IJ)*CGROUP_EXT(IJ,M)+OMOSNH2KD_EXT(IJ,M)*DRDP(IJ)
                DTHM = DRGM(IJ)*CGROUP_EXT(IJ,M)+OMOSNH2KD_EXT(IJ,M)*DRDM(IJ)
                CFLTP(IJ) = DTHP+ABS(DTHP)
                CFLTM(IJ) = -DTHM+ABS(DTHM) 
                DTC(IJ) =  DTC(IJ) + CFLTP(IJ) + CFLTM(IJ)
                DTP(IJ) = -DTHP+ABS(DTHP)
                DTM(IJ) =  DTHM+ABS(DTHM)
              ENDDO
            ENDIF


!*    3.1.4.2.4 LOOP OVER GRIDPOINTS.
!               ---------------------

            DO IJ = KIJS,KIJL
              IJPHR=IJRPH(IJ)
              IJLAR=IJRLA(IJ)
              F3(IJ,K,M) = (1.0_JWRB-DTC(IJ))*F1(IJ,K,M )               &
     &         + DPN(IJ) * F1(KLAT(IJ,IJPH,1),K  ,M)                    &
     &         + DPN2(IJ)* F1(KLAT(IJ,IJPH,2),K  ,M)                    &
     &         + DLE(IJ) * F1(KLON(IJ,IJLA),K  ,M)                      &
     &         + RPN(IJ) * F1(KRLAT(IJ,IJPHR,1),K  ,M)                  &
     &         + RPN2(IJ)* F1(KRLAT(IJ,IJPHR,2),K  ,M)                  &
     &         + RLE(IJ) * F1(KRLON(IJ,IJLAR,1),K  ,M)                  &
     &         + RLE2(IJ)* F1(KRLON(IJ,IJLAR,2),K  ,M)                  &
     &         + DTP(IJ) * F1(IJ           ,KP1,M)                      &
     &         + DTM(IJ) * F1(IJ           ,KM1,M)
            ENDDO

!           TEST THE STABILITY OF THE ADVECTION SCHEME
!           ------------------------------------------
            IF (LLCHKCFL .AND. M == 1) THEN
!             the non rotated quadrant has the strictest CFL criteria
              DO IJ = KIJS,KIJL
                CGOMFULL=1.0_JWRB/(1.0_JWRB-PRQRT(IJ))
                CFLEA(IJ) = CFLEA(IJ)*CGOMFULL
                CFLNO(IJ) = CFLNO(IJ)*CGOMFULL
                DTC(IJ) = CFLEA(IJ)+CFLNO(IJ)+CFLTP(IJ)+CFLTM(IJ)
              ENDDO
              CALL CHECKCFL(KIJS, KIJL, DEPTH_EXT(KIJS), DTC,     &
     &                      CFLEA,CFLEA,CFLNO,CFLNO,CFLNO,        &
     &                      CFLNO,CFLTP,CFLTM,CFLTP,CFLTM)
            ENDIF


!*    BRANCH BACK TO 3.1.4.2 FOR NEXT FREQUENCY.
          ENDDO

!*    BRANCH BACK TO 3.1 FOR NEXT DIRECTION.

      ENDDO

!*    3.2 END OF PROPAGATION FOR SPHERICAL GRID
!*        WITHOUT REFRACTION OR DEPTH REFRACTION, RETURN.
!         -----------------------------------------------

      DEALLOCATE(DELLA0,DCO,DLADCO,DP1,DP2,DPN,DPN2)
      DEALLOCATE(DCOM1)
      DEALLOCATE(DLE,DPH,DLA)
      DEALLOCATE(RLE,RLE2,RPN,RPN2)
      DEALLOCATE(DTP,DTM,DRGP,DRGM,DTC)
      DEALLOCATE(WLATM1)
      DEALLOCATE(WRLATM1)
      DEALLOCATE(WRLONM1)
      DEALLOCATE(IJRLA)
      DEALLOCATE(IJRPH)
      DEALLOCATE(DELPH0_CDA)
      DEALLOCATE(DELPH0R)
      DEALLOCATE(DELLA0R)
      DEALLOCATE(SDDL)
      DEALLOCATE(CGKLON)
      DEALLOCATE(CGKLAT)
      DEALLOCATE(ACGKRLON)
      DEALLOCATE(ACGKRLAT)
      IF (IREFRA == 1) THEN
        DEALLOCATE(DRDP)
        DEALLOCATE(DRDM)
      ENDIF

  
      DEALLOCATE(CFLEA)
      DEALLOCATE(CFLNO)
      DEALLOCATE(CFLTP)
      DEALLOCATE(CFLTM)
      DEALLOCATE(CFLEAR)
      DEALLOCATE(CFLNOR)

      IF (LHOOK) CALL DR_HOOK('PROPAGS1',1,ZHOOK_HANDLE)

      RETURN

! ----------------------------------------------------------------------

!*    4. PROPAGATION FOR SPHERICAL LATITUDE/LONGITUDE GRID
!*       WITH DEPTH AND CURRENT REFRACTION.
!        -------------------------------------------------

 4000 CONTINUE

      ALLOCATE(DPN(KIJS:KIJL))
      ALLOCATE(DPN2(KIJS:KIJL))
      ALLOCATE(DPS(KIJS:KIJL))
      ALLOCATE(DPS2(KIJS:KIJL))
      ALLOCATE(DLE(KIJS:KIJL))
      ALLOCATE(DLW(KIJS:KIJL))

      ALLOCATE(DOP(KIJS:KIJL))
      ALLOCATE(DOM(KIJS:KIJL))
      ALLOCATE(DRCP(KIJS:KIJL))
      ALLOCATE(DRCM(KIJS:KIJL))
      ALLOCATE(DRDP(KIJS:KIJL))
      ALLOCATE(DRDM(KIJS:KIJL))

      ALLOCATE(DRGP(KIJS:KIJL))
      ALLOCATE(DRGM(KIJS:KIJL))

      ALLOCATE(CFLEA(KIJS:KIJL))
      ALLOCATE(CFLWE(KIJS:KIJL))
      ALLOCATE(CFLNO(KIJS:KIJL))
      ALLOCATE(CFLSO(KIJS:KIJL))
      ALLOCATE(CFLNO2(KIJS:KIJL))
      ALLOCATE(CFLSO2(KIJS:KIJL))
      ALLOCATE(CFLTP(KIJS:KIJL))
      ALLOCATE(CFLTM(KIJS:KIJL))
      ALLOCATE(CFLOP(KIJS:KIJL))
      ALLOCATE(CFLOM(KIJS:KIJL))

      DELPH0 = 0.25_JWRB*DELPRO/DELPHI
      DELTH0 = 0.25_JWRB*DELPRO/DELTH
      DELLA0 = 0.25_JWRB*DELLA0
      DELFR0 = 0.25_JWRB*DELPRO/((FRATIO-1.0_JWRB)*ZPI)

      ALLOCATE(WLATM1(KIJS:KIJL,2))
      DO IC=1,2
        DO IJ = KIJS,KIJL
          WLATM1(IJ,IC) = 1.0_JWRB - WLAT(IJ,IC)
        ENDDO
      ENDDO



!*    4.1 LOOP OVER DIRECTIONS.
!         ---------------------

      DO K=1,NANG
        KP1 = K+1
        IF (KP1 > NANG) KP1 = 1
        KM1 = K-1
        IF (KM1 < 1) KM1 = NANG
        SD = SINTH(K)
        CD = COSTH(K)

!*    4.1.1 COMPUTE GRID REFRACTION.
!           ------------------------

        SP = DELTH0*(SINTH(K)+SINTH(KP1))/R
        SM = DELTH0*(SINTH(K)+SINTH(KM1))/R
        DO IJ = KIJS,KIJL
          JH = BLK2GLO%KXLT(IJ)
          TANPH = SINPH(JH)/COSPH(JH)
          DRGP(IJ) = TANPH*SP
          DRGM(IJ) = TANPH*SM
        ENDDO

!*    4.1.2 COMPUTE DEPTH REFRACTION.
!           -------------------------

        DO IJ = KIJS,KIJL
          DRDP(IJ) = (THDD(IJ,K) + THDD(IJ,KP1))*DELTH0
          DRDM(IJ) = (THDD(IJ,K) + THDD(IJ,KM1))*DELTH0
        ENDDO

!*    4.1.3 COMPUTE CURRENT REFRACTION.
!           ---------------------------

        DO IJ = KIJS,KIJL
          DRCP(IJ) = (THDC(IJ,K) + THDC(IJ,KP1))*DELTH0
          DRCM(IJ) = (THDC(IJ,K) + THDC(IJ,KM1))*DELTH0
        ENDDO

!*    4.1.4 LOOP OVER FREQUENCIES.
!           ----------------------

        DO M=1,NFRE_RED
          MP1 = MIN(NFRE_RED,M+1)
          MM1 = MAX(1,M-1)

            DFP = DELFR0/FR(M)
            DFM = DELFR0/FR(MM1)

!*    4.1.4.2.2 LON/LAT/DIR WEIGHTS IN INTEGRATION SCHEME.
!               ------------------------------------------

            DLA(NLAND) = SD*CGROUP_EXT(NLAND,M)
            DPH(NLAND) = CD*CGROUP_EXT(NLAND,M)*DELPH0
            DO IJ=NINF,NSUP
              DLA(IJ) = U_EXT(IJ)+SD*CGROUP_EXT(IJ,M)
              DPH(IJ) =(V_EXT(IJ)+CD*CGROUP_EXT(IJ,M))*DELPH0
            ENDDO
            DO IJ=KIJS,KIJL
              DLWE = (DLA(IJ) + DLA(KLON(IJ,1)))*DELLA0(IJ)*DCO(IJ)
              DLEA = (DLA(IJ) + DLA(KLON(IJ,2)))*DELLA0(IJ)*DCO(IJ)
              DLE(IJ) = (-DLEA+ABS(DLEA))*OBSLON(IJ,M,2)
              DLW(IJ) = ( DLWE+ABS(DLWE))*OBSLON(IJ,M,1)
              CFLEA(IJ) =  DLEA+ABS(DLEA)
              CFLWE(IJ) =  -DLWE+ABS(DLWE)
              DTC(IJ) =  CFLEA(IJ) + CFLWE(IJ)

!             IRREGULAR GRID
              IF (IRGG == 1) THEN
                DPNO = (DPH(IJ)+ DPH(KLAT(IJ,2,1))*DP2(IJ))*WLAT(IJ,2)
                DPN(IJ) = (-DPNO+ABS(DPNO))*OBSLAT(IJ,M,2)
                DPNO2= (DPH(IJ)+ DPH(KLAT(IJ,2,2))*DP2(IJ))*WLATM1(IJ,2)
                DPN2(IJ)= (-DPNO2+ABS(DPNO2))*OBSLAT(IJ,M,2)

                DPSO = (DPH(IJ)+ DPH(KLAT(IJ,1,1))*DP1(IJ))*WLAT(IJ,1)
                DPS(IJ) = ( DPSO+ABS(DPSO))*OBSLAT(IJ,M,1)
                DPSO2= (DPH(IJ)+ DPH(KLAT(IJ,1,2))*DP1(IJ))*WLATM1(IJ,1)
                DPS2(IJ)= ( DPSO2+ABS(DPSO2))*OBSLAT(IJ,M,1)
                CFLNO(IJ) =  DPNO+ABS(DPNO)
                CFLSO(IJ) =  -DPSO+ABS(DPSO)
                CFLNO2(IJ) =  DPNO2+ABS(DPNO2)
                CFLSO2(IJ) = -DPSO2+ABS(DPSO2)
                DTC(IJ) = DTC(IJ) +  CFLNO(IJ) + CFLSO(IJ)+             &
     &                               CFLNO2(IJ) + CFLSO2(IJ)
              ELSE
!             REGULAR GRID
                DPNO = DPH(IJ) + DPH(KLAT(IJ,2,1))*DP2(IJ)
                DPN(IJ) = (-DPNO+ABS(DPNO))*OBSLAT(IJ,M,2)
                DPSO = DPH(IJ) + DPH(KLAT(IJ,1,1))*DP1(IJ)
                DPS(IJ) = ( DPSO+ABS(DPSO))*OBSLAT(IJ,M,1)
                CFLNO(IJ) = DPNO+ABS(DPNO)
                CFLSO(IJ) = -DPSO+ABS(DPSO)
                CFLNO2(IJ) = 0.0_JWRB
                CFLSO2(IJ) = 0.0_JWRB
                DTC(IJ) = DTC(IJ) +  CFLNO(IJ) + CFLSO(IJ)
              ENDIF

              DTHP=DRGP(IJ)*CGROUP_EXT(IJ,M)                                 &
     &         +OMOSNH2KD_EXT(IJ,M)*DRDP(IJ)+DRCP(IJ)
              DTHM=DRGM(IJ)*CGROUP_EXT(IJ,M)                                 &
     &         +OMOSNH2KD_EXT(IJ,M)*DRDM(IJ)+DRCM(IJ)
              CFLTP(IJ) = DTHP+ABS(DTHP)
              CFLTM(IJ) = -DTHM+ABS(DTHM) 
              DTC(IJ) =  DTC(IJ) + CFLTP(IJ) + CFLTM(IJ)

              DTP(IJ) = -DTHP+ABS(DTHP)
              DTM(IJ) =  DTHM+ABS(DTHM)

              DTHP = (SDOT(IJ,K,M) + SDOT(IJ,K,MP1))*DFP
              DTHM = (SDOT(IJ,K,M) + SDOT(IJ,K,MM1))*DFM
              CFLOP(IJ) = DTHP+ABS(DTHP)
              CFLOM(IJ) = -DTHM+ABS(DTHM) 
              DTC(IJ) =  DTC(IJ) + CFLOP(IJ) + CFLOM(IJ)
              DOP(IJ) = (-DTHP+ABS(DTHP))/FRATIO
              DOM(IJ) = ( DTHM+ABS(DTHM))*FRATIO
            ENDDO

!         TEST THE STABILITY OF THE ADVECTION SCHEME
!         ------------------------------------------
          IF (LLCHKCFL .AND. M == 1) THEN
            CALL CHECKCFL(KIJS, KIJL, DEPTH_EXT(KIJS), DTC,      &
     &                    CFLEA,CFLWE,CFLNO,CFLSO,CFLNO2,        &
     &                    CFLSO2,CFLTP,CFLTM,CFLOP,CFLOM)
          ENDIF

!*    4.1.4.3 LOOP OVER GRIDPOINTS.
!             ---------------------

!         IRREGULAR GRID
          IF (IRGG == 1) THEN
            DO IJ = KIJS,KIJL
              F3(IJ,K,M) = (1.0_JWRB-DTC(IJ))*F1(IJ,K,M )               &
     &         + DPN(IJ) * F1(KLAT(IJ,2,1),K  ,M)                       &
     &         + DPN2(IJ)* F1(KLAT(IJ,2,2),K  ,M)                       &
     &         + DPS(IJ) * F1(KLAT(IJ,1,1),K  ,M)                       &
     &         + DPS2(IJ)* F1(KLAT(IJ,1,2),K  ,M)                       &
     &         + DLE(IJ) * F1(KLON(IJ,2),K  ,M)                         &
     &         + DLW(IJ) * F1(KLON(IJ,1),K  ,M)                         &
     &         + DTP(IJ) * F1(IJ        ,KP1,M)                         &
     &         + DTM(IJ) * F1(IJ        ,KM1,M)                         &
     &         + DOP(IJ) * F1(IJ        ,K  ,MP1)                       &
     &         + DOM(IJ) * F1(IJ        ,K  ,MM1)
            ENDDO
          ELSE
!           REGULAR GRID
            DO IJ = KIJS,KIJL
              F3(IJ,K,M) = (1.0_JWRB-DTC(IJ))*F1(IJ,K,M )               &
     &         + DPN(IJ) * F1(KLAT(IJ,2,1),K  ,M)                       &
     &         + DPS(IJ) * F1(KLAT(IJ,1,1),K  ,M)                       &
     &         + DLE(IJ) * F1(KLON(IJ,2),K  ,M)                         &
     &         + DLW(IJ) * F1(KLON(IJ,1),K  ,M)                         &
     &         + DTP(IJ) * F1(IJ        ,KP1,M)                         &
     &         + DTM(IJ) * F1(IJ        ,KM1,M)                         &
     &         + DOP(IJ) * F1(IJ        ,K  ,MP1)                       &
     &         + DOM(IJ) * F1(IJ        ,K  ,MM1)
            ENDDO
          ENDIF

!*    BRANCH BACK TO 4.1.4 FOR NEXT FREQUENCY.

        ENDDO

!*    BRANCH BACK TO 4.2 FOR NEXT DIRECTION.

      ENDDO

!*    4.4 END OF PROPAGATION FOR SPHERICAL GRID
!*        WITH DEPTH AND CURRENT REFRACTION, RETURN.
!         ------------------------------------------

      DEALLOCATE(DELLA0,DCO,DP1,DP2,DPN,DPN2,DPS,DPS2)
      DEALLOCATE(DCOM1)
      DEALLOCATE(DLE,DLW,DPH,DLA)
      DEALLOCATE(DTP,DTM,DRGP,DRGM,DTC)
      DEALLOCATE(WLATM1)

      DEALLOCATE(DOP,DOM,DRCP,DRCM)

      DEALLOCATE(DRDP)
      DEALLOCATE(DRDM)

      DEALLOCATE(CFLEA)
      DEALLOCATE(CFLWE)
      DEALLOCATE(CFLNO)
      DEALLOCATE(CFLSO)
      DEALLOCATE(CFLNO2)
      DEALLOCATE(CFLSO2)
      DEALLOCATE(CFLTP)
      DEALLOCATE(CFLTM)
      DEALLOCATE(CFLOP)
      DEALLOCATE(CFLOM)

IF (LHOOK) CALL DR_HOOK('PROPAGS1',1,ZHOOK_HANDLE)

END SUBROUTINE PROPAGS1
