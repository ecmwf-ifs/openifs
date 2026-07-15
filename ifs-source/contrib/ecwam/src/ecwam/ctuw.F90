! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE CTUW (DELPRO, MSTART, MEND,                    & 
 &               KIJS, KIJL, NINF, NSUP, LCFLFAIL, ICALL, &
 &               BLK2GLO,                                 &
 &               WLATM1, WCORM1, DP,                      &
 &               CGROUP_EXT, OMOSNH2KD_EXT,               &
 &               COSPHM1_EXT, DEPTH_EXT, U_EXT, V_EXT )
! ----------------------------------------------------------------------

!**** *CTUW* - COMPUTATION OF THE CONER TRANSPORT SCHEME WEIGHTS.


!*    PURPOSE.
!     --------

!       COMPUTATION OF THE CORNER TRANSPORT UPSTREAM WEIGHTS
!       USED IN THE PROPAGATION FOR A GIVEN TIME STEP.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO

      USE YOWCURR  , ONLY : LLCFLCUROFF
      USE YOWFRED  , ONLY : FR       ,DELTH    ,FRATIO   ,COSTH    ,SINTH
      USE YOWGRID  , ONLY : SINPH    ,COSPH
      USE YOWMAP   , ONLY : IRGG    ,IPER     ,                                   &
     &                      XDELLA   ,ZDELLO   ,AMOWEP   ,AMOSOP 
      USE YOWPARAM , ONLY : NANG     ,NFRE_RED
      USE YOWPCONS , ONLY : ZPI      ,R        ,CIRC
      USE YOWREFD  , ONLY : THDD     ,THDC     ,SDOT
      USE YOWSTAT  , ONLY : ICASE    ,IREFRA
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,WLAT     ,KCOR     ,WCOR     ,    &
     &                      SUMWN    ,WLATN    ,WLONN    ,WCORN    ,              &
     &                      WKPMN    ,WMPMN    ,OBSLAT   ,OBSLON   ,OBSCOR   ,    &
     &                      JXO      ,JYO      ,KCR

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE EC_LUN   , ONLY : NULERR

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: DELPRO ! ADVECTION TIME STEP
      INTEGER(KIND=JWIM), INTENT(IN) :: MSTART, MEND ! FREQUENCY START AND END INDEXES OVER WHICH THE WEIGHT WILL BE COMPUTED
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL  ! GRID POINT INDEXES
      INTEGER(KIND=JWIM), INTENT(IN) :: NINF, NSUP  ! HALO EXTEND NINF:NSUP+1
      INTEGER(KIND=JWIM), INTENT(IN) :: ICALL ! INDICATES IF IT IS THE FIRST OR SECOND CALL
!                                               THE SECOND CALL WILL TRY TO TEST WHETHER OR NOT
!                                               CFL IS VIOLATED WHEN THE CURRENT REFRACTION TERMS ARE SET TO 0
!                                               FOR THOSE POINTS WHERE IT WAS VIOLATED AT THE FIRST CALL  
      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO  ! BLOCK TO GRID TRANSFORMATION
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP,2), INTENT(IN) :: WLATM1  ! 1 - WLAT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP,4), INTENT(IN) :: WCORM1  ! 1 - WCOR
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP,2), INTENT(IN) :: DP      ! COS PHI FACTOR
      LOGICAL, DIMENSION(KIJS:KIJL), INTENT(INOUT) :: LCFLFAIL ! TRUE IF CFL CRITERION WAS VIOLATED.
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED), INTENT(IN) :: CGROUP_EXT  ! GROUP VELOCITY
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED), INTENT(IN) :: OMOSNH2KD_EXT ! OMEGA / SINH(2KD)
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: COSPHM1_EXT ! 1/COSPH
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: DEPTH_EXT ! WATER DEPTH
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: U_EXT ! U-COMPONENT OF SURFACE CURRENT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: V_EXT ! V-COMPONENT OF SURFACE CURRENT


      INTEGER(KIND=JWIM) :: NLAND
      INTEGER(KIND=JWIM) :: IP,IJP
      INTEGER(KIND=JWIM) :: K,M,IJ,IC,IX,KY,KK,KKM
      INTEGER(KIND=JWIM) :: KP1,KM1,JH
      INTEGER(KIND=JWIM) :: ICL,ICR,ICC, ICRM, JCR
      INTEGER(KIND=JWIM) :: MP1, MM1
      INTEGER(KIND=JWIM) :: ISSU(2),ISSV(2)

      REAL(KIND=JWRB) :: DTNEW
      REAL(KIND=JWRB) :: DXP, DYP, XM 
      REAL(KIND=JWRB) :: CMTODEG
      REAL(KIND=JWRB) :: GRIDAREAM1
      REAL(KIND=JWRB) :: CGYP, CGTH
      REAL(KIND=JWRB) :: DELTH0, DELFR0, SP, SM, DTHP, DTHM, DFP, DFM 
      REAL(KIND=JWRB) :: TANPH
      REAL(KIND=JWRB) :: DXX, DYY
      REAL(KIND=JWRB) :: UU, VV, UREL, VREL
      REAL(KIND=JWRB) :: XLAT, XLON
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(2) :: ADXP, ADYP
      REAL(KIND=JWRB), DIMENSION(2) :: DXUP, DXDW, DYUP, DYDW
      REAL(KIND=JWRB), DIMENSION(4) :: WEIGHT
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: DRGP,DRGM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: DRDP,DRDM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: DRCP,DRCM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CURMASK
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,2) :: CGX, CGY


! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CTUW',0,ZHOOK_HANDLE)


      NLAND = NSUP+1

      CMTODEG = 360.0_JWRB/CIRC

      IF (ICALL == 1) THEN
        LCFLFAIL(:) = .FALSE.
        CURMASK(:) = 1.0_JWRB
      ELSE
!       determine mask to turn off current refraction where it had failed before
        DO IJ=KIJS,KIJL
          IF (LCFLFAIL(IJ)) THEN
            CURMASK(IJ) = 0.0_JWRB
          ELSE
            CURMASK(IJ) = 1.0_JWRB
          ENDIF
        ENDDO
!       reset
        LCFLFAIL(:) = .FALSE.
      ENDIF

!*    ADVECTION IN PHYSICAL SPACE
!     =========================== 

!*    SPHERICAL OR CARTESIAN GRID?
!     ----------------------------
      IF (ICASE == 1) THEN

!*      SPHERICAL GRID.
!       ---------------

!         FIND THE RELATIVE WEIGHT IN
!         THE CONER TRANSPORT UPSTREAM SCHEME.

!*        LOOP OVER FREQUENCIES.
!         ----------------------
          DO M = MSTART, MEND

!*          LOOP OVER DIRECTIONS.
!           ---------------------
            DO K=1,NANG

!             FIND MEAN GROUP VELOCITY COMPONENTS FOR DIRECTION TH(K)+180
!             -----------------------------------------------------------
                DO IC=1,2
                  DO IJ=KIJS,KIJL
                    CGX(IJ,IC)=                                              &
     &                 0.5_JWRB*(CGROUP_EXT(IJ,M)+CGROUP_EXT(KLON(IJ,IC),M)) &
     &                    *SINTH(K)*COSPHM1_EXT(IJ)
!                   IRREGULAR GRID
                    IF (IRGG == 1) THEN
                      CGYP=WLAT(IJ,IC)*CGROUP_EXT(KLAT(IJ,IC,1),M)+          &
     &                 (1.0_JWRB-WLAT(IJ,IC))*CGROUP_EXT(KLAT(IJ,IC,2),M)
                    ELSE
!                   REGULAR GRID
                      CGYP=CGROUP_EXT(KLAT(IJ,IC,1),M)
                    ENDIF
                    CGY(IJ,IC)=0.5_JWRB*(CGROUP_EXT(IJ,M)+DP(IJ,IC)*CGYP)*COSTH(K)
                  ENDDO
                ENDDO


!             LOOP OVER GRID POINTS
!             ---------------------
              DO IJ=KIJS,KIJL
                IX=BLK2GLO%IXLG(IJ)
                KY=BLK2GLO%KXLT(IJ)

!               FLUX VELOCITIES AT THE GRID BOX INTERFACE 

                DO IC=1,2

                  IF (IREFRA == 2 .OR. IREFRA == 3 ) THEN
                    UU=U_EXT(IJ)*COSPHM1_EXT(IJ)
                    UREL=CGX(IJ,IC)+UU
                    ISSU(IC)=ISAMESIGN(UREL,CGX(IJ,IC))
                    VV=V_EXT(IJ)*0.5_JWRB*(1.0_JWRB+DP(IJ,IC))
                    VREL=CGY(IJ,IC)+VV
                    ISSV(IC)=ISAMESIGN(VREL,CGY(IJ,IC))
                  ELSE
                    UREL=CGX(IJ,IC)
                    ISSU(IC)=1
                    VREL=CGY(IJ,IC)
                    ISSV(IC)=1
                  ENDIF
                  DXP=-DELPRO*UREL*CMTODEG
                  DYP=-DELPRO*VREL*CMTODEG 

                  ADXP(IC)=ABS(DXP)
                  ADYP(IC)=ABS(DYP)

!                 BASIC CFL CHECKS (IN EACH DIRECTION)
!                 ----------------
                  IF (ADXP(IC) > ZDELLO(KY))THEN
                    WRITE (IU06,*) '********************************'
                    WRITE (IU06,*) '* CTUW:                        *'
                    WRITE (IU06,*) '* CFL VIOLATED IN X DIRECTION. *'
                    WRITE (IU06,*) '* ADXP SHOULD BE < ZDELLO, BUT *'
                    WRITE (IU06,*) '* ADXP = ',ADXP(IC),IC
                    WRITE (IU06,*) '* ZDELLO = ',ZDELLO(KY)
                    DTNEW=ZDELLO(KY)*DELPRO/ADXP(IC)
                    WRITE (IU06,*) '* TIME STEP ',DELPRO
                    WRITE (IU06,*) '* SHOULD BE REDUCED TO ', DTNEW
                    WRITE (IU06,*) '*                              *'
                    WRITE (IU06,*) '********************************'
                    LCFLFAIL(IJ)=.TRUE.
                  ENDIF
                  IF (ADYP(IC) > XDELLA)THEN
                    XLON=AMOWEP+(IX-1)*ZDELLO(KY)
                    XLAT=AMOSOP+(KY-1)*XDELLA
                    DTNEW=XDELLA*DELPRO/ADYP(IC)
                    WRITE (IU06,*) '********************************'
                    WRITE (IU06,*) '* CTUW:                        *'
                    WRITE (IU06,*) '* CFL VIOLATED IN Y DIRECTION. *'
                    WRITE (IU06,*) '* ADYP SHOULD BE < XDELLA, BUT *'
                    WRITE (IU06,*) '* ADYP = ',ADYP(IC),IC
                    WRITE (IU06,*) '* XDELLA = ',XDELLA
                    WRITE (IU06,*) '* XLAT= ',XLAT,' XLON= ',XLON 
                    WRITE (IU06,*) '* DEPTH= ',DEPTH_EXT(IJ)
                    WRITE (IU06,*) '* TIME STEP ',DELPRO
                    WRITE (IU06,*) '* SHOULD BE REDUCED TO ', DTNEW
                    WRITE (IU06,*) '*                              *'
                    WRITE (IU06,*) '********************************'
                    LCFLFAIL(IJ)=.TRUE.
                  ENDIF

                  DXUP(IC)=ADXP(IC)*ISSU(IC)
                  DXDW(IC)=ADXP(IC)*(1-ISSU(IC))
                  DYUP(IC)=ADYP(IC)*ISSV(IC)
                  DYDW(IC)=ADYP(IC)*(1-ISSV(IC))
    
                ENDDO

!                 GET ADVECTION WEIGHT FOR ALL NEIGHBOURING GRID POINTS

                  DXX=ZDELLO(KY)-DXUP(JXO(K,2))-DXDW(JXO(K,1))
                  DYY=XDELLA-DYUP(JYO(K,2))-DYDW(JYO(K,1))

                  GRIDAREAM1 =  1.0_JWRB/(ZDELLO(KY)*XDELLA)

!                 WEIGHTED CONTRIBUTION FROM NORTH-SOUTH DIRECTION (WLATN)

                  WEIGHT(JYO(K,1))=DXX*DYUP(JYO(K,1))*GRIDAREAM1
                  WEIGHT(JYO(K,2))=DXX*DYDW(JYO(K,2))*GRIDAREAM1
                  DO IC=1,2
                    WLATN(IJ,K,M,IC,1)=WLAT(IJ,IC)*WEIGHT(IC)
                    WLATN(IJ,K,M,IC,2)=WLATM1(IJ,IC)*WEIGHT(IC)
                  ENDDO

!                 WEIGHTED CONTRIBUTION FROM EAST-WEST DIRECTION (WLONN)

                  WLONN(IJ,K,M,JXO(K,1))=DYY*DXUP(JXO(K,1))*GRIDAREAM1
                  WLONN(IJ,K,M,JXO(K,2))=DYY*DXDW(JXO(K,2))*GRIDAREAM1


!                 CONTRIBUTION FROM CORNERS (KCOR)
                  WEIGHT(1)=DXUP(JXO(K,1))*DYUP(JYO(K,1))*GRIDAREAM1
                  WEIGHT(2)=DXDW(JXO(K,2))*DYUP(JYO(K,1))*GRIDAREAM1
                  WEIGHT(3)=DXUP(JXO(K,1))*DYDW(JYO(K,2))*GRIDAREAM1
                  WEIGHT(4)=DXDW(JXO(K,2))*DYDW(JYO(K,2))*GRIDAREAM1
                  DO ICR=1,4
                    WCORN(IJ,K,M,ICR,1)=WCOR(IJ,KCR(K,ICR))*WEIGHT(ICR)
                    WCORN(IJ,K,M,ICR,2)=WCORM1(IJ,KCR(K,ICR))*WEIGHT(ICR)
                  ENDDO

!                 CONTRIBUTIONS FOR IJ
                  SUMWN(IJ,K,M)=(ZDELLO(KY)*                            &
     &                            (DYDW(JYO(K,1))+DYUP(JYO(K,2))) +     &
     &                           XDELLA*                                &
     &                            (DXUP(JXO(K,2))+DXDW(JXO(K,1))) -     &
     &                           (DXDW(JXO(K,1))+DXUP(JXO(K,2)))*       &
     &                           (DYDW(JYO(K,1))+DYUP(JYO(K,2)))  )     &
     &                           *GRIDAREAM1

                ENDDO  ! END LOOP OVER GRID POINTS


            ENDDO  ! END LOOP OVER DIRECTIONS

          ENDDO  ! END LOOP OVER FREQUENCIES


      ELSE
!*    CARTESIAN GRID.
!     ---------------
        IF (IREFRA == 2 .OR. IREFRA == 3 ) THEN
!*      WITHOUT DEPTH OR/AND CURRENT REFRACTION.
!       ----------------------------------------
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '* CTUW:                                  *'
          WRITE (IU06,*) '* CORNER TRANSPORT SCHEME NOT YET READY  *' 
          WRITE (IU06,*) '* FOR  CARTESIAN GRID !                  *'
          WRITE (IU06,*) '* FOR DEPTH OR/AND CURRENT REFRACTION !  *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          IF (LHOOK) CALL DR_HOOK('CTUW',1,ZHOOK_HANDLE)
          CALL ABORT1
        ELSE
!*      DEPTH AND CURRENT REFRACTION.
!       ----------------------------
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '* CTUW:                                  *'
          WRITE (IU06,*) '* CORNER TRANSPORT SCHEME NOT YET READY  *' 
          WRITE (IU06,*) '* FOR  CARTESIAN GRID !                  *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          IF (LHOOK) CALL DR_HOOK('CTUW',1,ZHOOK_HANDLE)
          CALL ABORT1
        ENDIF

      ENDIF


!*    REFRACTIONS TERMS
!     ================= 

      DELTH0 = 0.25*DELPRO/DELTH

!*    LOOP OVER DIRECTIONS.
!     ---------------------

      DO K=1,NANG
        KP1 = K+1
        IF (KP1 > NANG) KP1 = 1
        KM1 = K-1
        IF (KM1 < 1) KM1 = NANG

!*      COMPUTE GRID REFRACTION.
!       ------------------------
        SP  = DELTH0*(SINTH(K)+SINTH(KP1))/R
        SM  = DELTH0*(SINTH(K)+SINTH(KM1))/R

        DO IJ = KIJS,KIJL
          JH=BLK2GLO%KXLT(IJ)
          TANPH = SINPH(JH)/COSPH(JH)
          DRGP(IJ) = TANPH*SP
          DRGM(IJ) = TANPH*SM
        ENDDO

!*      COMPUTE DEPTH REFRACTION.
!       -------------------------
        IF (IREFRA == 1) THEN
          DO IJ = KIJS,KIJL
            DRDP(IJ) = (THDD(IJ,K) + THDD(IJ,KP1))*DELTH0
            DRDM(IJ) = (THDD(IJ,K) + THDD(IJ,KM1))*DELTH0
          ENDDO
        ELSE
          DO IJ = KIJS,KIJL
            DRDP(IJ) =  0.0_JWRB
            DRDM(IJ) =  0.0_JWRB
          ENDDO
        ENDIF

!*      COMPUTE CURRENT REFRACTION.
!       ---------------------------

        IF (IREFRA == 2 .OR. IREFRA == 3 ) THEN
          DO IJ = KIJS,KIJL
            DRCP(IJ) = CURMASK(IJ)*(THDC(IJ,K) + THDC(IJ,KP1))*DELTH0
            DRCM(IJ) = CURMASK(IJ)*(THDC(IJ,K) + THDC(IJ,KM1))*DELTH0
          ENDDO
        ELSE
          DO IJ = KIJS,KIJL
            DRCP(IJ) = 0.0_JWRB 
            DRCM(IJ) = 0.0_JWRB
          ENDDO
        ENDIF


!*      REFRACTION WEIGHTS IN INTEGRATION SCHEME.
!       -----------------------------------------

!*      NO DEPTH REFRACTION.
!       -------------------
        IF (IREFRA == 0) THEN
          DO M = MSTART, MEND
            DO IJ=KIJS,KIJL
              DTHP = DRGP(IJ)*CGROUP_EXT(IJ,M) + DRCP(IJ)
              DTHM = DRGM(IJ)*CGROUP_EXT(IJ,M) + DRCM(IJ)
              WKPMN(IJ,K,M,0)=(DTHP+ABS(DTHP))+(ABS(DTHM)-DTHM)
              WKPMN(IJ,K,M,1)=-DTHP+ABS(DTHP)
              WKPMN(IJ,K,M,-1)=DTHM+ABS(DTHM)
            ENDDO
          ENDDO
        ELSE
!*      SHALLOW WATER AND DEPTH REFRACTION.
!       -----------------------------------
          DO M = MSTART, MEND
            DO IJ=KIJS,KIJL
              DTHP = DRGP(IJ)*CGROUP_EXT(IJ,M)+OMOSNH2KD_EXT(IJ,M)*DRDP(IJ)+DRCP(IJ)
              DTHM = DRGM(IJ)*CGROUP_EXT(IJ,M)+OMOSNH2KD_EXT(IJ,M)*DRDM(IJ)+DRCM(IJ)
              WKPMN(IJ,K,M,0)=(DTHP+ABS(DTHP))+(ABS(DTHM)-DTHM)
              WKPMN(IJ,K,M,1)=-DTHP+ABS(DTHP)
              WKPMN(IJ,K,M,-1)=DTHM+ABS(DTHM)
            ENDDO
          ENDDO
        ENDIF

!*      COMPUTE FREQUENCY SHIFTING DUE TO CURRENTS.
!       -------------------------------------------

        IF (IREFRA == 2 .OR. IREFRA == 3 ) THEN

          DELFR0 = 0.25_JWRB*DELPRO/((FRATIO-1)*ZPI)

            DO M = MSTART, MEND
              MP1 = MIN(NFRE_RED,M+1)
              MM1 = MAX(1,M-1)
              DFP = DELFR0/FR(M)
              DFM = DELFR0/FR(MM1)

              DO IJ=KIJS,KIJL
                DTHP = CURMASK(IJ) * (SDOT(IJ,K,M) + SDOT(IJ,K,MP1))*DFP
                DTHM = CURMASK(IJ) * (SDOT(IJ,K,M) + SDOT(IJ,K,MM1))*DFM
                WMPMN(IJ,K,M,0) =(DTHP+ABS(DTHP))+(ABS(DTHM)-DTHM)
                WMPMN(IJ,K,M,1) =(-DTHP+ABS(DTHP))/FRATIO
                WMPMN(IJ,K,M,-1)=(DTHM+ABS(DTHM))*FRATIO
              ENDDO
            ENDDO
        ENDIF

      ENDDO  ! END LOOP ON DIRECTIONS


!     CHECK THAT WEIGHTS ARE LESS THAN 1
!     AND COMPUTE THEIR SUM AND CHECK IT IS LESS THAN 1 AS WELL
!!!   THE SUM IS NEEDED LATER ON !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO K=1,NANG
        DO M = MSTART, MEND
          DO IJ=KIJS,KIJL
            DO IC=1,2
              DO ICL=1,2
              IF (WLATN(IJ,K,M,IC,ICL) > 1.0_JWRB .OR.                       &
     &            WLATN(IJ,K,M,IC,ICL) < 0.0_JWRB     ) THEN
                WRITE(IU06,*) '***********************************'
                WRITE(IU06,*) '* CTUW:                           *'
                WRITE(IU06,*) '* CFL VIOLATED IN Y DIRECTION     *'
                WRITE(IU06,*) '* ICALL = ', ICALL
                WRITE(IU06,*) '* WLATN SHOULD BE < 1 AND > 0, BUT*'
                WRITE(IU06,*) '* WLATN(IJ,K,M,IC,ICL)=',IJ,K,M,IC,ICL, &
     &                            WLATN(IJ,K,M,IC,ICL)
                WRITE(IU06,*) '*                                 *'
                WRITE(IU06,*) '***********************************'
                LCFLFAIL(IJ)=.TRUE.
                CALl FLUSH(IU06)
              ENDIF
              ENDDO
            ENDDO

            DO IC=1,2
              IF (WLONN(IJ,K,M,IC) > 1.0_JWRB .OR.                           &
     &            WLONN(IJ,K,M,IC) < 0.0_JWRB      ) THEN
                WRITE (IU06,*) '***********************************'
                WRITE (IU06,*) '* CTUW:                           *'
                WRITE (IU06,*) '* CFL VIOLATED IN X DIRECTION     *'
                WRITE (IU06,*) '* WLONN SHOULD BE < 1 AND > 0, BUT*'
                WRITE (IU06,*) '* WLONN(IJ,K,M,IC)= ',IJ,K,M,IC,        &
     &                          WLONN(IJ,K,M,IC)
                WRITE (IU06,*) '*                                 *'
                WRITE (IU06,*) '***********************************'
                LCFLFAIL(IJ)=.TRUE.
                CALl FLUSH(IU06)
              ENDIF
            ENDDO

            DO ICR=1,4
              DO ICL=1,2
              IF (WCORN(IJ,K,M,ICR,ICL) > 1.0_JWRB .OR.                      &
     &            WCORN(IJ,K,M,ICR,ICL) < 0.0_JWRB     ) THEN
                WRITE (IU06,*) '***********************************'
                WRITE (IU06,*) '* CTUW:                           *'
                WRITE (IU06,*) '* CFL VIOLATED IN CORNER DIRECTION*'
                WRITE (IU06,*) '* WCORN SHOULD BE < 1 AND > 0, BUT*'
                WRITE (IU06,*) '*WCORN(IJ,K,M,ICR,ICL)=',IJ,K,M,ICR,ICL, &
     &                            WCORN(IJ,K,M,ICR,ICL)
                WRITE (IU06,*) '*                                 *'
                WRITE (IU06,*) '***********************************'
                LCFLFAIL(IJ)=.TRUE.
                CALl FLUSH(IU06)
              ENDIF
              ENDDO
            ENDDO

            DO IC=-1,1
              IF (WKPMN(IJ,K,M,IC) > 1.0_JWRB .OR.                           &
     &            WKPMN(IJ,K,M,IC) < 0.0_JWRB      ) THEN
                WRITE (IU06,*) '***********************************'
                WRITE (IU06,*) '* CTUW:                           *'
                WRITE (IU06,*) '* CFL VIOLATED IN DIRECTION SPACE *'
                WRITE (IU06,*) '* WKPMN SHOULD BE < 1 AND > 0, BUT*'
                WRITE (IU06,*) '* WKPMN(IJ,K,M,IC)= ',IJ,K,M,IC,        &
     &                            WKPMN(IJ,K,M,IC)
                WRITE (IU06,*) '*                                 *'
                WRITE (IU06,*) '***********************************'
                LCFLFAIL(IJ)=.TRUE.
                CALl FLUSH(IU06)
              ENDIF
            ENDDO

            SUMWN(IJ,K,M)=SUMWN(IJ,K,M)+WKPMN(IJ,K,M,0)


            IF (IREFRA == 2 .OR. IREFRA == 3 ) THEN

              DO IC=-1,1
                IF (WMPMN(IJ,K,M,IC) > 1.0_JWRB .OR.                         &
     &              WMPMN(IJ,K,M,IC) < 0.0_JWRB      ) THEN
                  WRITE(IU06,*) '***********************************'
                  WRITE(IU06,*) '* CTUW:                           *'
                  WRITE(IU06,*) '* CFL VIOLATED IN FREQUENCY SPACE *'
                  WRITE(IU06,*) '* ICALL = ', ICALL
                  WRITE(IU06,*) '* WMPMN SHOULD BE < 1 AND > 0, BUT*'
                  WRITE(IU06,*) '* WMPMN(IJ,K,M,IC)= ',IJ,K,M,IC,      &
     &                              WMPMN(IJ,K,M,IC)
                  WRITE(IU06,*) '*                                 *'
                  WRITE(IU06,*) '***********************************'

                  WRITE(NULERR,*) '* CTUW: CFL VIOLATED IN FREQUENCY*',ICALL,IJ,K,M,IC,WMPMN(IJ,K,M,IC),U_EXT(IJ),V_EXT(IJ)

                  LCFLFAIL(IJ)=.TRUE.
                  CALl FLUSH(IU06)
                ENDIF
              ENDDO

              SUMWN(IJ,K,M)=SUMWN(IJ,K,M)+WMPMN(IJ,K,M,0)
            ENDIF

!           SUM < 1  ?
            IF (SUMWN(IJ,K,M) > 1.0_JWRB .OR. SUMWN(IJ,K,M) < 0.0_JWRB) THEN
              IX=BLK2GLO%IXLG(IJ)
              KY=BLK2GLO%KXLT(IJ)
              XLON=AMOWEP+(IX-1)*ZDELLO(KY)
              XLAT=AMOSOP+(KY-1)*XDELLA
              WRITE(IU06,*) '***********************************'
              WRITE(IU06,*) '* CTUW:                           *'
              WRITE(IU06,*) '* CFL VIOLATED                    *'
              WRITE(IU06,*) '* ICALL = ', ICALL
              WRITE(IU06,*) '* SUMW SHOULD BE < 1 AND > 0, BUT*'
              WRITE(IU06,*) '* IJ, SUMWN(IJ) = ',IJ,SUMWN(IJ,K,M)
              WRITE(IU06,*) '* XLAT= ',XLAT,' XLON= ',XLON 
              WRITE(IU06,*) '* DEPTH= ',DEPTH_EXT(IJ)
              IF (.NOT. LLCFLCUROFF .OR. ICALL > 1 ) THEN 
                WRITE(NULERR,*) '* CTUW: SUMW IS NOT <1 AND >0, BUT*',ICALL,IJ,K,M,SUMWN(IJ,K,M),XLAT,&
     &                     XLON,DEPTH_EXT(IJ),U_EXT(IJ),V_EXT(IJ)
              ENDIF
              DO IP=1,2
              DO IC=1,2
                IJP = KLAT(IJ,IC,IP)
                IF (IJP /= NLAND) THEN
                  WRITE(IU06,*) '* DEPTH= ',IJP,IP,IC,DEPTH_EXT(IJP)
                  WRITE(IU06,*) '*     U= ',U_EXT(IJP)
                  WRITE(IU06,*) '*     V= ',V_EXT(IJP)
                ELSE
                  WRITE(IU06,*) '* DEPTH= ',IJP,IP,IC,'LAND'
                ENDIF
              ENDDO
              ENDDO
              DO IC=1,2
                IJP = KLON(IJ,IC)
                IF (IJP /= NLAND) THEN
                  WRITE(IU06,*) '* DEPTH= ',IJP,IC,DEPTH_EXT(IJP)
                  WRITE(IU06,*) '*     U= ',U_EXT(IJP)
                  WRITE(IU06,*) '*     V= ',V_EXT(IJP)
                ELSE
                  WRITE(IU06,*) '* DEPTH= ',IJP,IC,'LAND'
                ENDIF
              ENDDO

              IF (IREFRA == 2 .OR. IREFRA == 3 ) THEN
              WRITE(IU06,*) '* U = ',U_EXT(IJ),' V = ',V_EXT(IJ)
              ENDIF
              WRITE(IU06,*) '*                                 *'
              WRITE(IU06,*) '***********************************'
              LCFLFAIL(IJ)=.TRUE.

              CALl FLUSH(IU06)
            ENDIF

          ENDDO  ! END LOOP OVER GRID POINTS
        ENDDO  ! END LOOP OVER FREQUENCIES
      ENDDO  ! END LOOP OVER DIRECTIONS

      DO IJ=KIJS,KIJL
        IF (LCFLFAIL(IJ)) THEN
          IF (LHOOK) CALL DR_HOOK('CTUW',1,ZHOOK_HANDLE)
          RETURN
        ENDIF
      ENDDO


!!!!!!INCLUDE THE BLOCKING COEFFICIENTS INTO THE WEIGHTS OF THE
!     SURROUNDING POINTS.

      DO K=1,NANG
        DO M = MSTART, MEND
          DO IJ=KIJS,KIJL

!           POINTS ON SURROUNDING LATITUDES 
            DO IC=1,2
              DO ICL=1,2
                WLATN(IJ,K,M,IC,ICL) = WLATN(IJ,K,M,IC,ICL)*OBSLAT(IJ,M,IC) 
              ENDDO
            ENDDO

!           POINTS ON SURROUNDING LONGITUDE
            DO IC=1,2
              WLONN(IJ,K,M,IC) = WLONN(IJ,K,M,IC)*OBSLON(IJ,M,IC)
            ENDDO

!           SURROUNDING CORNER POINTS
            DO ICR=1,4
              DO ICL=1,2
                WCORN(IJ,K,M,ICR,ICL) = WCORN(IJ,K,M,ICR,ICL)*OBSCOR(IJ,M,KCR(K,ICR))
              ENDDO
            ENDDO

          ENDDO  ! END LOOP OVER GRID POINTS
        ENDDO  ! END LOOP ON FREQUENCIES
      ENDDO  ! END LOOP OVER DIRECTIONS

IF (LHOOK) CALL DR_HOOK('CTUW',1,ZHOOK_HANDLE)

      RETURN

      CONTAINS

! ----------------------------------------------------------------------

      INTEGER(KIND=JWIM) FUNCTION ISAMESIGN(A,B)
!       =1 IF A AND B HAVE THE SAME SIGN, 0 OTHERWISE
        USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
        IMPLICIT NONE

        REAL(KIND=JWRB), INTENT(IN) :: A,B

        IF (SIGN(1.0_JWRB,A) == SIGN(1.0_JWRB,B)) THEN
          ISAMESIGN=1
        ELSE
          ISAMESIGN=0
        ENDIF
      END FUNCTION ISAMESIGN

END SUBROUTINE CTUW
