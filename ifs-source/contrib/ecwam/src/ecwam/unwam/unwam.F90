! (C) Copyright 2001- Aron Roland (Roland & Partner, Germany).
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE UNWAM

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWPARAM    , ONLY : NANG, NFRE, NFRE_RED
      USE YOWSTAT     , ONLY : IREFRA
      USE YOWFRED     , ONLY : COSTH, SINTH
      USE yowpd, only: MNE=>ne, INE, MNP=>npa, NP_RES => np
      USE yowpd, only: XP=>x, YP=>y, DEP=>z
      USE yowpd, only: exchange
      USE yowMpiModule
      USE UNSTRUCT_BOUND , ONLY : SPSIG, IWBMNP, IWBNDLC, LBCWA
      USE UNSTRUCT_BOUND , ONLY : IOBP, IOBPD, IOBWB
      USE UNSTRUCT_BOUND,  ONLY : WBAC
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal
      USE UNSTRUCT_CURR  , ONLY : CURTXY
! ------------------------------------------------------------------------------
      IMPLICIT NONE

      SAVE
      PRIVATE

      LOGICAL :: APPLY_DXP_CORR = .TRUE.
!!!! JB:
!!!! USE_EXACT_FORMULA_SPHERICAL_AREA = .TRUE. DOES NOT SEEM TO WORK !!!
!!! I need to investigate further
      LOGICAL :: USE_EXACT_FORMULA_SPHERICAL_AREA = .FALSE.
      LOGICAL :: LNANINFCHK = .TRUE.
      REAL(KIND=JWRU) :: DT4A, DT4D, DT4F
      REAL(KIND=JWRU) :: DDIR
      REAL(KIND=JWRU) :: KDMAX = 300.0_JWRU
      REAL(KIND=JWRU) :: ONE = 1.0_JWRU
      REAL(KIND=JWRU) :: ZERO = 0.0_JWRU
      REAL(KIND=JWRU) :: DMIN = 0.0_JWRU
      REAL(KIND=JWRU), ALLOCATABLE    :: CAD_THE(:,:,:)
      REAL(KIND=JWRU), ALLOCATABLE    :: CAS_SIG(:,:,:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: ID_PREV(:), ID_NEXT(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: Jstart(:)
      INTEGER(KIND=JWIM) :: recTime1, recTime2

      INTEGER(KIND=JWIM), ALLOCATABLE :: JA_IE(:,:,:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: IA(:), JA(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: POSI(:,:), I_DIAG(:)
      INTEGER(KIND=JWIM) :: NNZ
      INTEGER(KIND=JWIM) :: REFRA_METHOD = 2;

      REAL(KIND=JWRU) :: DGRTH, DGRTHM1
      REAL(KIND=JWRU), PARAMETER :: YPMAX=87.5_JWRU

      REAL(KIND=JWRB), ALLOCATABLE :: DEPTH_INPUT(:)

      REAL(KIND=JWRU), ALLOCATABLE :: COS2TH(:)
      REAL(KIND=JWRU), ALLOCATABLE :: SIN2TH(:)
      REAL(KIND=JWRU), ALLOCATABLE :: SINCOSTH(:)
      REAL(KIND=JWRU), ALLOCATABLE :: DS_BAND(:)
      REAL(KIND=JWRU), ALLOCATABLE :: ASPAR_JAC(:,:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: B_JAC(:,:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: DTP_I(:,:,:), DTM_I(:,:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: DOP_I(:,:,:), DOM_I(:,:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: U_JACOBI(:,:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: DS_INCR(:)
      REAL(KIND=JWRU), ALLOCATABLE :: DDEP(:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: DCUX(:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: DCUY(:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: DEPDT(:)
      !
      ! Solver thresholds
      !
      INTEGER(KIND=JWIM) :: maxiter = 100
      REAL(KIND=JWRU) :: PMIN
      REAL(KIND=JWRU) :: PTAIL5
      REAL(KIND=JWRU), ALLOCATABLE :: INVTRANS1(:)


PUBLIC :: PROPAG_UNWAM, &
     &    INIT_UNWAM, &
     &    UNWAM_IN, &
     &    SET_UNWAM_HANDLES, &
     &    UNWAM_OUT, &
     &    EXCHANGE_FOR_FL1

!     THE FOLLOWING FLAG ARE PART OF THE INPUT NAMELIST
!     SEE *MPUSERIN*
      REAL(KIND=JWRU), PUBLIC :: WAE_SOLVERTHR
      REAL(KIND=JWRU), PUBLIC :: JGS_DIFF_SOLVERTHR
      LOGICAL, PUBLIC :: LIMPLICIT
      LOGICAL, PUBLIC :: SOURCE_IMPL
      LOGICAL, PUBLIC :: LNONL
      LOGICAL, PUBLIC :: BLOCK_GAUSS_SEIDEL
      LOGICAL, PUBLIC :: LLIMT
      LOGICAL, PUBLIC :: L_SOLVER_NORM
      LOGICAL, PUBLIC :: LCHKCONV

INTERFACE PROPAG_UNWAM 
  MODULE PROCEDURE PROPAG_UNWAM
END INTERFACE

INTERFACE INIT_UNWAM
  MODULE PROCEDURE INIT_UNWAM
END INTERFACE

INTERFACE UNWAM_IN  
  MODULE PROCEDURE UNWAM_IN
END INTERFACE

INTERFACE SET_UNWAM_HANDLES 
  MODULE PROCEDURE SET_UNWAM_HANDLES
END INTERFACE

INTERFACE UNWAM_OUT
  MODULE PROCEDURE UNWAM_OUT
END INTERFACE

INTERFACE EXCHANGE_FOR_FL1
  MODULE PROCEDURE EXCHANGE_FOR_FL1
END INTERFACE

      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PROPAG_UNWAM(FL1, FLNEW)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        ------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!       Driver routine for the fluctuation splitting schemes 
!       Parallelization over freq. and directions unsing OpenMP untill dom. decomp. is not ready 
!       References: Roland, 2008, Roland et al. 2006, Zanke et al. 2006, Csik et al. 2002 
!                   Hubbard & Roe, 2000, Struis et al. 1998 

!     Externals.  EXPLICIT_N_SCHEME, EXPLICIT_LFPSI_SCHEME, EXPLICIT_PSI_SCHEME
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF*

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      USE YOWUNPOOL
!      USE yowpd
      USE YOWTEST  , ONLY : IU06
      USE yownodepool, ONLY : NP, iplg, ipgl
      USE yowdatapool, ONLY: myrank
      USE YOWMPP   , ONLY : NINF, NSUP
      USE UNSTRUCT_BOUND, ONLY : SET_UP_WBAC, APPLY_BOUNDARY_CONDITION

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(INOUT)  :: FL1(NINF-1:NSUP,NANG,NFRE_RED)
      REAL(KIND=JWRB), INTENT(INOUT)  :: FLNEW(1:MNP,NANG,NFRE)
 
      INTEGER(KIND=JWIM) :: IS, ID, IP
      REAL(KIND=JWRU) :: FLsing(NANG,NFRE)
      !
      ! Compute differentials if needed
      !

      IF (LCALC) THEN
        IF (IREFRA .ne. 0) THEN
          CALL DIFFERENTIATE_XYDIR_LSPHE(DEP, DDEP)
          CALL DIFFERENTIATE_XYDIR_LSPHE(CURTXY(1,:), DCUX)
          CALL DIFFERENTIATE_XYDIR_LSPHE(CURTXY(2,:), DCUY)
          DEPDT(:) = 0.0_JWRU
        END IF
      END IF
      IF (LBCWA) THEN
        CALL SET_UP_WBAC
        CALL APPLY_BOUNDARY_CONDITION(FL1)
      END IF
      IF (LIMPLICIT) THEN
        CALL IMPLICIT_N_SCHEME_BLOCK(FL1, FLNEW)
      ELSE
        IF (LVECTOR) THEN 
          CALL EXPLICIT_N_SCHEME_VECTOR(FL1,FLNEW)
        ELSE
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ID,IS)
          DO ID = 1, NANG
            DO IS = 1, NFRE
              CALL EXPLICIT_N_SCHEME(IS,ID,FL1(:,ID,IS),FLNEW(:,ID,IS))
             !CALL EXPLICIT_PSI_SCHEME(IS,ID,FL1(:,ID,IS),FLNEW(:,ID,IS))
             !CALL EXPLICIT_LF_SCHEME(IS,ID,FL1(:,ID,IS),FLNEW(:,ID,IS))
            END DO
          END DO
!$OMP     END PARALLEL DO
          DO IP=1,MNP
            FLsing=REAL(FLNEW(IP,:,:), JWRU)
            CALL REFRACTION_FREQSHIFT_EXPLICIT_SINGLE(FLsing, IP)
            FLNEW(IP,:,:)=REAL(FLsing,JWRB)
          END DO
        ENDIF
      ENDIF
      LCALC = .FALSE.
      END SUBROUTINE PROPAG_UNWAM
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DIFFERENTIATE_XYDIR_LSPHE(VAR, GRAD)
      USE YOWUNPOOL, ONLY : LSPHE, DEGRAD, REARTH
      USE YOWUNPOOL, ONLY : IEN, TRIA
      IMPLICIT NONE
      REAL(KIND=JWRU), INTENT(IN)  :: VAR(MNP)
      REAL(KIND=JWRU), INTENT(OUT) :: GRAD(2,MNP)
      REAL(KIND=JWRU) :: DVDX(MNP), DVDY(MNP)
      INTEGER(KIND=JWIM) :: NI(3)
      INTEGER(KIND=JWIM) :: IE, I1, I2, I3, IP
      REAL(KIND=JWRU) :: DEDY(3),DEDX(3)
      REAL(KIND=JWRU) :: DVDXIE, DVDYIE
      REAL(KIND=JWRU) :: WEI(MNP)
      WEI(:)  = 0._JWRU
      DVDX(:) = 0._JWRU
      DVDY(:) = 0._JWRU
!      WRITE(740+MyRankGlobal,*) 'DIFFERENTIATE_XYDIR'
!      WRITE(740+MyRankGlobal,*) 'sum(VAR )=', sum(VAR)
!      WRITE(740+MyRankGlobal,*) 'sum(IEN)=', sum(IEN)
!      WRITE(740+MyRankGlobal,*) 'sum(TRIA)=', sum(TRIA)
      DO IE = 1, MNE
        NI = INE(:,IE)
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        WEI(NI) = WEI(NI) + 2.*TRIA(IE)
        DEDX(1) = IEN(1,IE)
        DEDX(2) = IEN(3,IE)
        DEDX(3) = IEN(5,IE)
        DEDY(1) = IEN(2,IE)
        DEDY(2) = IEN(4,IE)
        DEDY(3) = IEN(6,IE)
        DVDXIE  = DOT_PRODUCT( VAR(NI),DEDX)
        DVDYIE  = DOT_PRODUCT( VAR(NI),DEDY)
        DVDX(NI) = DVDX(NI) + DVDXIE
        DVDY(NI) = DVDY(NI) + DVDYIE
      END DO
      WRITE(740+MyRankGlobal,*) 'sum(WEI)=', sum(WEI)
      IF (LSPHE) THEN
        DO IP=1,MNP
          DVDX(IP) = DVDX(IP) * INVTRANS1(IP) / WEI(IP)
          DVDY(IP) = DVDY(IP) * DGRTHM1 / WEI(IP)
        END DO
      ELSE
        DO IP=1,MNP
          DVDX(IP) = DVDX(IP) / WEI(IP)
          DVDY(IP) = DVDY(IP) / WEI(IP)
        END DO
      END IF
      CALL exchange(DVDX)
      CALL exchange(DVDY)
!      WRITE(740+MyRankGlobal,*) 'sum(DVDX)=', sum(DVDX)
!      WRITE(740+MyRankGlobal,*) 'sum(DVDY)=', sum(DVDY)
      DO IP=1,MNP
        GRAD(1,IP) = DVDX(IP)
        GRAD(2,IP) = DVDY(IP)
      END DO
      END SUBROUTINE DIFFERENTIATE_XYDIR_LSPHE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DIFFERENTIATE_XYDIR(VAR, GRAD)
      USE YOWUNPOOL, ONLY : IEN, TRIA
      IMPLICIT NONE
      REAL(KIND=JWRU), INTENT(IN)  :: VAR(MNP)
      REAL(KIND=JWRU), intent(out) :: GRAD(2,MNP)
      REAL(KIND=JWRU) :: DVDX(MNP), DVDY(MNP)
      INTEGER(KIND=JWIM) :: NI(3)
      INTEGER(KIND=JWIM) :: IE, I1, I2, I3, IP
      REAL(KIND=JWRU) :: DEDY(3),DEDX(3)
      REAL(KIND=JWRU) :: DVDXIE, DVDYIE
      REAL(KIND=JWRU) :: WEI(MNP)
      WEI(:)  = 0._JWRU
      DVDX(:) = 0._JWRU
      DVDY(:) = 0._JWRU
!      WRITE(740+MyRankGlobal,*) 'DIFFERENTIATE_XYDIR'
!      WRITE(740+MyRankGlobal,*) 'sum(VAR )=', sum(VAR)
!      WRITE(740+MyRankGlobal,*) 'sum(IEN)=', sum(IEN)
!      WRITE(740+MyRankGlobal,*) 'sum(TRIA)=', sum(TRIA)
      DO IE = 1, MNE
        NI = INE(:,IE)
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        WEI(NI) = WEI(NI) + 2.*TRIA(IE)
        DEDX(1) = IEN(1,IE)
        DEDX(2) = IEN(3,IE)
        DEDX(3) = IEN(5,IE)
        DEDY(1) = IEN(2,IE)
        DEDY(2) = IEN(4,IE)
        DEDY(3) = IEN(6,IE)
        DVDXIE  = DOT_PRODUCT( VAR(NI),DEDX)
        DVDYIE  = DOT_PRODUCT( VAR(NI),DEDY)
        DVDX(NI) = DVDX(NI) + DVDXIE
        DVDY(NI) = DVDY(NI) + DVDYIE
      END DO
      WRITE(740+MyRankGlobal,*) 'sum(WEI)=', sum(WEI)
      DVDX(:) = DVDX(:)/WEI(:)
      DVDY(:) = DVDY(:)/WEI(:)
      CALL exchange(DVDX)
      CALL exchange(DVDY)
      WRITE(740+MyRankGlobal,*) 'sum(DVDX)=', sum(DVDX)
      WRITE(740+MyRankGlobal,*) 'sum(DVDY)=', sum(DVDY)
      DO IP=1,MNP
        GRAD(1,IP) = DVDX(IP)
        GRAD(2,IP) = DVDY(IP)
      END DO
      END SUBROUTINE DIFFERENTIATE_XYDIR
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_CURTXY
      USE YOWUNPOOL, ONLY : LCALC
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, KIJL4CHNK, IJFROMCHNK
      USE YOWSHAL  , ONLY : WVENVI
      IMPLICIT NONE
      INTEGER(KIND=JWIM) :: IP, ICHNK, IPRM
      DO ICHNK = 1, NCHNK
        DO IPRM = 1, KIJL4CHNK(ICHNK)
          IP = IJFROMCHNK(IPRM, ICHNK)
          CURTXY(1,IP)=REAL(WVENVI%UCUR(IPRM, ICHNK),JWRU)
          CURTXY(2,IP)=REAL(WVENVI%VCUR(IPRM, ICHNK),JWRU)
        ENDDO
      ENDDO
      LCALC=.TRUE.
      END SUBROUTINE SET_CURTXY
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PROPTHETA(IP, CAD)
      USE YOWUNPOOL, ONLY : LCUR, LSPHE, DEGRAD, REARTH
      USE YOWUNPOOL, ONLY : WK, CG
      IMPLICIT NONE
      INTEGER(KIND=JWIM), INTENT(IN)        :: IP
      REAL(KIND=JWRU), INTENT(OUT)   :: CAD(NANG,NFRE)
      INTEGER(KIND=JWIM) :: IS, ID
      REAL(KIND=JWRU) :: WKDEP, DWDH, CFL
#ifdef DEBUG
!      WRITE(740+MyRankGlobal,*) ' maxval(DDEP)', maxval(DDEP)
!      WRITE(740+MyRankGlobal,*) 'DEP   = ', DEP(IP)
!      WRITE(740+MyRankGlobal,*) 'DMIN  = ', DMIN
!      WRITE(740+MyRankGlobal,*) 'LSPHE = ', LSPHE
!      WRITE(740+MyRankGlobal,*) 'LCUR  = ', LCUR
!      WRITE(740+MyRankGlobal,*) 'maxval(CG)=', maxval(CG)
#endif
      IF (LSPHE) THEN
        DO IS = 1, NFRE
          DO ID = 1, NANG
!        it is SINTH(ID) because when theta is 0, it means propagation towards the north,
!        hence along a meridian, which is a great circle, hence no need to change the
!        wave propagation direction.
            CAD(ID,IS) = CG(IS,IP)*SINTH(ID)*TAN(YP(IP)*DEGRAD)/REARTH
          END DO
        END DO
      ELSE
        CAD(:,:) = 0.0_JWRU
      END IF

!     Assume a meaningful water depth before considering bottom and/or current refraction
      IF (DEP(IP) .GT. DMIN) THEN
        IF ((IREFRA.eq.1).or.(IREFRA.eq.3)) THEN
          DO IS = 1, NFRE
            WKDEP = WK(IS,IP) * DEP(IP)
            IF (WKDEP .LT. 13._JWRU) THEN
              DWDH = SPSIG(IS)/SINH(MIN(KDMAX,2._JWRU*WKDEP))
              DO ID = 1, NANG
                CAD(ID,IS) = CAD(ID,IS) + DWDH * ( SINTH(ID)*DDEP(1,IP)-COSTH(ID)*DDEP(2,IP) )
              END DO
            ENDIF
          END DO
        END IF
        IF (LCUR.and.((IREFRA.eq.2).or.(IREFRA.eq.3))) THEN
          DO IS = 1, NFRE
            DO ID = 1, NANG
              CAD(ID,IS) = CAD(ID,IS) + SIN2TH(ID)*DCUY(1,IP)-COS2TH(ID)*DCUX(2,IP)+SINCOSTH(ID)*( DCUX(1,IP)-DCUY(2,IP) )
            END DO
          END DO
        END IF
      END IF
      END SUBROUTINE PROPTHETA
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PROPSIGMA(IP, CAS)
      USE YOWUNPOOL, ONLY : CG, WK
      IMPLICIT NONE
      INTEGER(KIND=JWIM), INTENT(IN) :: IP
      REAL(KIND=JWRU), INTENT(OUT) :: CAS(NFRE,NANG)
      REAL(KIND=JWRU) :: DWDH, WKDEP
      INTEGER(KIND=JWIM) :: ID, IS
      IF (DEP(IP) .GT. DMIN) THEN
        DO IS = 1, NFRE
          WKDEP = WK(IS,IP) * DEP(IP)
          IF (WKDEP .LT. 13._JWRU) THEN
            DWDH = SPSIG(IS)/SINH(MIN(KDMAX,2._JWRU*WKDEP))
          ELSE
            DWDH = 0._JWRU
          END IF
          DO ID = 1, NANG
             CAS(IS,ID) = DWDH * WK(IS,IP) * ( DEPDT(IP) + CURTXY(1,IP)*DDEP(IP,1) + &
&             CURTXY(2,IP)*DDEP(IP,2) ) - CG(IS,IP) * WK(IS,IP) * ( COS2TH(ID)*DCUX(IP,1) + SIN2TH(ID)*DCUY(IP,2) + &
&             SINCOSTH(ID)*( DCUY(IP,1) + DCUX(IP,2) ) )
          END DO
        END DO
      END IF
      END SUBROUTINE PROPSIGMA
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE QUICKEST_FREQ(MX, Q, CS, DT, DX1, DX2)
      IMPLICIT NONE
      INTEGER(KIND=JWIM), INTENT(IN)        :: MX
      REAL(KIND=JWRU), INTENT(IN)    :: DT
      REAL(KIND=JWRU), INTENT(INOUT) :: Q(0:MX+1)
      REAL(KIND=JWRU), INTENT(IN)    :: CS(0:MX+1)
      REAL(KIND=JWRU), INTENT(IN)    :: DX1(0:MX+1), DX2(0:MX+1)
      REAL(KIND=JWRU)              :: FLA(0:MX+1)
      INTEGER(KIND=JWIM)           :: IXY, IXYC, IXYU, IXYD
      REAL(KIND=JWRU)              :: CSM, CFL
      REAL(KIND=JWRU)              :: DQ, DQNZ, QCN, QBN, QBR, QB
      REAL(KIND=JWRU)              :: CFAC
!
!     Fluxes for central points (Fi,+)
!
      DO IXY = 1, MX-1
        CSM    = 0.5_JWRU * ( CS(IXY) + CS(IXY+1) )
        CFL    = DT *  CSM / DX2(IXY)
        IXYC   = IXY - 1 * INT(MIN(ZERO , SIGN(1.1_JWRU,CFL) ) )
        QB     = 0.5_JWRU *((1.-CFL)*Q(IXY+1)+(1._JWRU+CFL)*Q(IXY)) -                   &
     &       DX2(IXY)**2/DX1(IXYC) * (1._JWRU-CFL**2) / 6._JWRU *                            &
     &       ( (Q(IXYC+1)-Q(IXYC))/DX2(IXYC)-(Q(IXYC)-Q(IXYC-1))/DX2(IXYC-1) )
        IXYU   = IXYC - 1 * INT ( SIGN (1.1_JWRU,CFL) )
        IXYD   = 2*IXYC - IXYU
        DQ     = Q(IXYD) - Q(IXYU)
        DQNZ   = SIGN ( MAX(1.E-15_JWRU,ABS(DQ)) , DQ )
        QCN    = ( Q(IXYC) - Q(IXYU) ) / DQNZ
        QCN    = MIN ( 1.1_JWRU, MAX ( -0.1_JWRU , QCN ) )
        QBN    = MAX ( (QB-Q(IXYU))/DQNZ , QCN )
        QBN    = MIN ( QBN , ONE , QCN/MAX(1.E-10_JWRU,ABS(CFL)) )
        QBR    = Q(IXYU) + QBN*DQ
        CFAC   = REAL(INT( 2._JWRU * ABS(QCN-0.5_JWRU) ), JWRU)
        QB     = (1.-CFAC)*QBR + CFAC*Q(IXYC)
        FLA(IXY) = CSM * QB
      END DO
!
!     Fluxes for points with boundary above
!
      IXY  = 0
      IXYC = IXY - 1 * INT( MIN( 0.0_JWRU, SIGN(1.1_JWRU,CFL) ) )
      FLA(0) = CS(0) * Q(0) ! The flux at the uper boundary of the spectrum is given in the calling routine ...
!
!     Fluxes for points with boundary below
!
      IXY  = MX
      IXYC = IXY - 1 * INT( MIN( 0.0_JWRU, SIGN(1.1_JWRU,CFL) ) )
      FLA(MX) = CS(MX+1) * Q(MX+1) ! The flux at the uper boundary of the spectrum is given in the calling routine
!
!     Fluxes for points with boundary above
!
      DO IXY = 1, MX
        Q(IXY) = MAX(ZERO, Q(IXY) - (FLA(IXY)-FLA(IXY-1)) * DT/DX1(IXY))
      END DO
      END SUBROUTINE QUICKEST_FREQ
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE QUICKEST_DIR(MX, Q, CS, DT, DX)
      IMPLICIT NONE
      INTEGER(KIND=JWIM), INTENT(IN) :: MX
      REAL(KIND=JWRU), INTENT(INOUT) :: Q(-1:MX+2), CS(-1:MX+2)
      REAL(KIND=JWRU), INTENT(IN)    :: DT, DX
      REAL(KIND=JWRU) :: CFLL(0:MX+1), FLA(0:MX+1)
      INTEGER(KIND=JWIM), PARAMETER :: INC = 1
      INTEGER(KIND=JWIM)  :: IXY, IXYC, IXYU, IXYD

      REAL(KIND=JWRU) :: CFL
      REAL(KIND=JWRU) :: DQ, DQNZ, QCN, QBN, QBR, QB
      REAL(KIND=JWRU) :: CFAC

      CS(-1)   = CS(MX-1)
      CS(0)    = CS(MX)
      CS(MX+1) = CS(1)
      CS(MX+2) = CS(2)
      Q(-1)    = Q(MX-1)
      Q(0)     = Q(MX)
      Q(MX+1)  = Q(1)
      Q(MX+2)  = Q(2)
      
      DO IXY = 0, MX+1
        CFLL(IXY) = CS(IXY) * DT / DX
      END DO
      !
      !     Fluxes for central points (Fi,+)
      !
      DO IXY = 0, MX
        CFL = 0.5_JWRU * ( CFLL(IXY) + CFLL(IXY+INC) )
        IXYC = IXY - INC * INT(MIN(ZERO, SIGN(1.1_JWRU,CFL) ) )
        QB   = 0.5_JWRU*((1.0-CFL)*Q(IXY+INC)+(1.0+CFL)*Q(IXY))-(1.0-CFL**2.0)/6.0*(Q(IXYC-INC)-2.0*Q(IXYC)+Q(IXYC+INC))
        IXYU = IXYC - INC * INT( SIGN(1.1_JWRU,CFL) )
        IXYD = 2*IXYC - IXYU
        DQ   = Q(IXYD) - Q(IXYU)   ! DEL = Qd - Qu
        DQNZ = SIGN( MAX(10E-15_JWRU, ABS(DQ)), DQ )
        QCN  = ( Q(IXYC) - Q(IXYU) ) / DQNZ   ! ~Qc = ( Qc - Qu ) / DEL
        QCN  = MIN( 1.1_JWRU, MAX( -0.1_JWRU, QCN ) )   ! -0.1 < ~Qc < 1.1_JWRU
        QBN  = MAX( (QB - Q(IXYU)) / DQNZ, QCN )   ! ~Qf = ( Qf - Qu ) / DEL
        QBN  = MIN( QBN, ONE, QCN/MAX( 1.0E-10_JWRU, ABS(CFL) ) )
        QBR  = Q(IXYU) + QBN * DQ   ! Qf = ~Qf * DEL + Qu
        CFAC = REAL(INT(2.0*ABS(QCN-0.5_JWRU)),JWRU)  ! if 0 < ~Qc < 1, then CFAC = 1.0, Qf = Qc
        QB   = (1.0 - CFAC) * QBR + CFAC * Q(IXYC)
        FLA(IXY) = CFL * QB
      END DO
      !
      !     Propagation
      !
      DO IXY = 1, MX
        Q(IXY) = MAX(ZERO,Q(IXY) + FLA(IXY-INC) - FLA(IXY))
      END DO
      END SUBROUTINE QUICKEST_DIR
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE REFRACTION_FREQSHIFT_EXPLICIT_SINGLE(FLsing, IP)
      USE YOWUNPOOL
      IMPLICIT NONE
      REAL(KIND=JWRU), intent(inout) :: FLsing(NANG,NFRE)
      INTEGER(KIND=JWIM), intent(in) :: IP
      REAL(KIND=JWRU) ACQ_F(0:NFRE+1)
      REAL(KIND=JWRU) ACQ_D(-1:NANG+2)
      REAL(KIND=JWRU) CASS(0:NFRE+1)
      REAL(KIND=JWRU) CAS(NFRE,NANG)
      REAL(KIND=JWRU) CAD(NANG,NFRE)
      REAL(KIND=JWRU) TheFL(NFRE)
      REAL(KIND=JWRU) REST, DT4FI, CFLCAS
      REAL(KIND=JWRU) DT4DI, CFLCAD
      REAL(KIND=JWRU) CADS(-1:NANG+2)
      INTEGER(KIND=JWIM) ID, ITER, IT, IS
      !
      IF (ABS(IOBP(IP)) .EQ. 1) RETURN
      IF (DEP(IP) .LT. DMIN) RETURN
      IF ((IREFRA .eq. 2).or.(IREFRA .eq. 3)) THEN
        CALL PROPSIGMA(IP,CAS)
!#ifdef DEBUG
!        WRITE(740+MyRankGlobal,*) 'NFRE=', NFRE
!        WRITE(740+MyRankGlobal,*) ' NANG=', NANG
!        FLUSH(740+MyRankGlobal)
!#endif
        DO ID = 1, NANG
!#ifdef DEBUG
!          WRITE(740+MyRankGlobal,*) 'ID=', ID
!          FLUSH(740+MyRankGlobal)
!#endif
          ACQ_F(1:NFRE)  = FLsing(ID,:) / SPSIG
!#ifdef DEBUG
!          WRITE(740+MyRankGlobal,*) 'after ACQ_F initialization'
!          FLUSH(740+MyRankGlobal)
!#endif
          CFLCAS  = MAXVAL(ABS(CAS(:,ID))*DT4F/DS_BAND(1:NFRE))
          REST  = ABS(MOD(CFLCAS,ONE))
          IF (REST .GT. THR .AND. REST .LT. 0.5_JWRU) THEN
            ITER = ABS(NINT(CFLCAS)) + 1
          ELSE
            ITER = ABS(NINT(CFLCAS))
          END IF
          DT4FI = DT4F / REAL(ITER,JWRU)
          CASS(1:NFRE) = CAS(:,ID)
          CASS(0)     = 0._JWRU
          CASS(NFRE+1) = CASS(NFRE)
          DO IT = 1, ITER ! Iteration
            ACQ_F(0)      = ACQ_F(1)
            ACQ_F(NFRE+1)  = ACQ_F(NFRE) * PTAIL5
            CALL QUICKEST_FREQ(NFRE,ACQ_F,CASS,DT4FI,DS_BAND,DS_INCR)
          END DO
          TheFL = ACQ_F(1:NFRE) * SPSIG
          FLsing(ID,:) = MAX(ZERO,TheFL)
        END DO
      END IF
      IF ((IREFRA .eq. 1).or.(IREFRA .eq. 3).or.(LSPHE.eqv..TRUE.)) THEN
        CALL PROPTHETA(IP,CAD)
#ifdef DEBUG
!        WRITE(740+MyRankGlobal,*) 'NFRE=', NFRE
!        WRITE(740+MyRankGlobal,*) ' NANG=', NANG
!        FLUSH(740+MyRankGlobal)
!        WRITE(740+MyRankGlobal,*) 'size(FLsing,1)=', size(FLsing,1)
!        WRITE(740+MyRankGlobal,*) 'size(FLsing,2)=', size(FLsing,2)
!        FLUSH(740+MyRankGlobal)
#endif
        DO IS = 1, NFRE
!#ifdef DEBUG
!          WRITE(740+MyRankGlobal,*) 'IS=', IS
!          FLUSH(740+MyRankGlobal)
!#endif
          ACQ_D(1:NANG) = FLsing(:,IS)
!#ifdef DEBUG
!          WRITE(740+MyRankGlobal,*) 'after ACQ_D initialization'
!          FLUSH(740+MyRankGlobal)
!#endif
          CADS(1:NANG) = CAD(:,IS)
          CFLCAD = MAXVAL(ABS(CAD(:,IS)))*DT4D/DDIR
#ifdef DEBUG
!          WRITE(740+MyRankGlobal,*) 'DDIR=', DDIR
!          WRITE(740+MyRankGlobal,*) 'DT4D=', DT4D
!          DO ID=1,NANG
!            WRITE(740+MyRankGlobal,*) 'ID=', ID, ' CAD=', CAD(ID,IS)
!          END DO
!          FLUSH(740+MyRankGlobal)
#endif
          IF (CFLCAD .LT. THR) CYCLE
          REST  = ABS(MOD(CFLCAD,ONE))
          IF (REST .GT. THR .AND. REST .LT. 0.5_JWRU) THEN
            ITER = ABS(NINT(CFLCAD)) + 1
          ELSE
            ITER = ABS(NINT(CFLCAD))
          END IF
          ITER = MAX(1,ITER)
#ifdef DEBUG
!          WRITE(740+MyRankGlobal,*) 'ITER=', ITER
!          FLUSH(740+MyRankGlobal)
#endif
          DT4DI = DT4D / REAL(ITER,JWRU)
          DO IT = 1, ITER ! Iteration
            CALL QUICKEST_DIR(NANG,ACQ_D,CADS,DT4DI,DDIR)
          END DO          ! end Interation
          FLsing(:,IS) = MAX(ZERO,ACQ_D(1:NANG))
        END DO
      END IF
      END SUBROUTINE REFRACTION_FREQSHIFT_EXPLICIT_SINGLE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE FLUCTCFL(IS, ID, DTMAX)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  EXPLICIT_N_SCHEME, EXPLICIT_LFPSI_SCHEME, EXPLICIT_PSI_SCHEME
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF*

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      USE YOWUNPOOL
!     --------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=JWRU) :: K(3,MNE)
      REAL(KIND=JWRB), INTENT(OUT) :: DTMAX
      INTEGER(KIND=JWIM) :: IS, ID
      INTEGER(KIND=JWIM) :: I, J, I1, I2, I3
      INTEGER(KIND=JWIM) :: IP, IE, POS
      REAL(KIND=JWRU) :: KSUM, KMAX, LAMBDA(2)
      REAL(KIND=JWRU) :: DTMAX_EXP, DTMAX_GLOBAL_EXP
      REAL(KIND=JWRU) :: REST, C(2,MNP)
      REAL(KIND=JWRB) :: DIFRU, USOC
      DTMAX_GLOBAL_EXP = 10.E14_JWRU
      CALL CADVXY(IS,ID,C)
      DO IE = 1, MNE
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        LAMBDA(1) = ONESIXTH * (C(1,I1)+C(1,I2)+C(1,I3))
        LAMBDA(2) = ONESIXTH * (C(2,I1)+C(2,I2)+C(2,I3))
        K(1,IE)  = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
        K(2,IE)  = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
        K(3,IE)  = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
      END DO
      J = 0
      DO IP = 1, MNP
        KSUM = 0.0_JWRU
        KMAX = 0.0_JWRU
        Jstart(IP) = J
        DO I = 1, CCON(IP)
          J = J + 1
          IE    = IE_CELL(J)
          POS   = POS_CELL(J)
          KSUM  = KSUM + MAX(K(POS,IE),0.0_JWRU)
          IF ( ABS(K(POS,IE)) > KMAX ) KMAX = ABS(K(POS,IE))
        END DO
        IF (KSUM > 0.0_JWRU) THEN
          DTMAX_EXP = SI(IP)/KSUM
        ELSE
          DTMAX_EXP = 10.E14_JWRU
        END IF
        IF (DTMAX_GLOBAL_EXP>DTMAX_EXP) DTMAX_GLOBAL_EXP = DTMAX_EXP
      END DO
      DTMAX = REAL(DTMAX_GLOBAL_EXP)
      REST  = ABS(MOD(DT4A/DTMAX_GLOBAL_EXP,1.0_JWRU))
      IF (REST > THR .AND. REST < 0.5_JWRU) THEN
        ITER_EXP(ID,IS) = ABS(NINT(DT4A/DTMAX_GLOBAL_EXP)) + 1
      ELSE
        ITER_EXP(ID,IS) = ABS(NINT(DT4A/DTMAX_GLOBAL_EXP))
      END IF
      END SUBROUTINE FLUCTCFL
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPLICIT_N_SCHEME(IS, ID, UOLD, UNEW)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Narrow Stencil (N) scheme using contour integration method (Csik et al.) for flux conservation 
!                 1st order time and space, monot1.d0, conservative, quasi-positive
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF*

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------
         USE YOWUNPOOL
         use yowdatapool, only: myrank
         USE yowpd, only: comm
         USE YOWMPP   , ONLY : NINF, NSUP, IRANK
         USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     -------------------------------------------------------------------------

         IMPLICIT NONE

         INTEGER(KIND=JWIM), INTENT(IN) :: IS,ID
         REAL(KIND=JWRB), INTENT(IN) :: UOLD(NINF-1:NSUP)
         REAL(KIND=JWRB), INTENT(OUT) :: UNEW(NINF-1:NSUP)
!
! local INTEGER
!
         INTEGER(KIND=JWIM) :: IP, IE, IT
         INTEGER(KIND=JWIM) :: I1, I2, I3
         INTEGER(KIND=JWIM) :: NI(3),K
!
! local double
!
         REAL(KIND=JWRU) :: FT
         REAL(KIND=JWRU) :: UTILDE
         REAL(KIND=JWRU) :: DTMAX_GLOBAL_EXP(1), DTMAX_EXP, DTMAX_GLOBAL_EXP_LOC(1)
         REAL(KIND=JWRU) :: REST
         REAL(KIND=JWRU) :: LAMBDA(2), DT4AI
         REAL(KIND=JWRU) :: FL11,FL12,FL21,FL22,FL31,FL32
         REAL(KIND=JWRU) :: KTMP(3)
         REAL(KIND=JWRU) :: KKSUM(MNP), ST(MNP), N(MNE), U3(3), ST3(3)
         REAL(KIND=JWRU) :: C(2,MNP), DTSI(MNP), CFLXY
         REAL(KIND=JWRU) :: FLALL(3,MNE), U(MNP)
         REAL(KIND=JWRU) :: FL111, FL112, FL211, FL212, FL311, FL312
         REAL(KIND=JWRU) :: KELEM(3,MNE)
!
! local parameter
!
         REAL(KIND=JWRU) :: TMP
         REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

         INTEGER(KIND=JWIM) :: ierr

#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) '--------------------------------------------'
      FLUSH(740+MyRankGlobal)
#endif

!     -------------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('EXPLICIT_N_SCHEME',0,ZHOOK_HANDLE)
#endif
!
!        Calculate phase speeds for the certain spectral comp1.d0nt ...
!
         CALL CADVXY(IS,ID,C)

!???JB is the call needed ???
!!!!!         call EXCHANGE(C(1,:))
!???JB is the call needed ???
!!!!         call EXCHANGE(C(2,:))

!         if(IS==1 .and. ID==1) then
!           write(DBG%FHNDL,*) "C"
!           do ip=1, np_global
!             if(ipgl(ip) /= 0) then
!               write(DBG%FHNDL,*) ip, C(:,ipgl(ip))
!             endif
!             if(ghostgl(ip) /= 0) then
!               write(DBG%FHNDL,*) ip, C(:,np+ghostgl(ip))
!             endif
!           end do
!         endif

!         DO IP = 1, MNP
!           C(:,IP) = C(:,IP) * IOBPD(ID,IP)
!         END DO
!
!        Calculate K-Values and contour based quantities ...
!

         DO IE = 1, MNE
            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)

            LAMBDA(1) = ONESIXTH *(C(1,I1)+C(1,I2)+C(1,I3))
            LAMBDA(2) = ONESIXTH *(C(2,I1)+C(2,I2)+C(2,I3))
            KELEM(1,IE) = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
            KELEM(2,IE) = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
            KELEM(3,IE) = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
            KTMP  = KELEM(:,IE)
            TMP   = SUM(MIN(0.0_JWRU,KTMP))
            N(IE) = - 1.0_JWRU/MIN(-THR,TMP)
            KELEM(:,IE) = MAX(0.0_JWRU,KTMP)
            FL11  = C(1,I2) * IEN(1,IE) + C(2,I2) * IEN(2,IE)
            FL12  = C(1,I3) * IEN(1,IE) + C(2,I3) * IEN(2,IE)
            FL21  = C(1,I3) * IEN(3,IE) + C(2,I3) * IEN(4,IE)
            FL22  = C(1,I1) * IEN(3,IE) + C(2,I1) * IEN(4,IE)
            FL31  = C(1,I1) * IEN(5,IE) + C(2,I1) * IEN(6,IE)
            FL32  = C(1,I2) * IEN(5,IE) + C(2,I2) * IEN(6,IE)
            FL111 = 2.0_JWRU*FL11+FL12
            FL112 = 2.0_JWRU*FL12+FL11
            FL211 = 2.0_JWRU*FL21+FL22
            FL212 = 2.0_JWRU*FL22+FL21
            FL311 = 2.0_JWRU*FL31+FL32
            FL312 = 2.0_JWRU*FL32+FL31
            FLALL(1,IE) = (FL311 + FL212) * ONESIXTH + KELEM(1,IE)
            FLALL(2,IE) = (FL111 + FL312) * ONESIXTH + KELEM(2,IE)
            FLALL(3,IE) = (FL211 + FL112) * ONESIXTH + KELEM(3,IE)
         END DO

         IF (LCALC) THEN ! If the current field or water level --> new CFL 
           KKSUM = 0.0_JWRU
           DO IE = 1, MNE
             NI = INE(:,IE)
             KKSUM(NI) = KKSUM(NI) + KELEM(:,IE)
           END DO

! THOMAS find global dt max/min
           DTMAX_GLOBAL_EXP(1) = LARGE
           DTMAX_GLOBAL_EXP_LOC(1) = LARGE
           DO IP = 1, MNP

             DTMAX_EXP = SI(IP)/MAX(SMALL, KKSUM(IP))
             IF (LCFL) THEN
               CFLCXY(1,IP) = MAX(CFLCXY(1,IP), C(1,IP))
               CFLCXY(2,IP) = MAX(CFLCXY(2,IP), C(2,IP))
               CFLCXY(3,IP) = MAX(CFLCXY(3,IP), DT4A/DTMAX_EXP)
             END IF
             DTMAX_GLOBAL_EXP_LOC(1) = MIN(DTMAX_GLOBAL_EXP_LOC(1),DTMAX_EXP)

           END DO


           CALL MPI_ALLREDUCE(DTMAX_GLOBAL_EXP_LOC, DTMAX_GLOBAL_EXP,    &
     &                        1, MPI_REAL8, MPI_MIN, COMM, IERR)

           CFLXY = DT4A/DTMAX_GLOBAL_EXP(1)
           REST  = ABS(MOD(CFLXY,1.0_JWRU))
           IF (REST .LT. THR) THEN
             ITER_EXP(ID,IS) = ABS(NINT(CFLXY))
           ELSE IF (REST .GT. THR .AND. REST .LT. 0.5_JWRU) THEN
             ITER_EXP(ID,IS) = ABS(NINT(CFLXY)) + 1
           ELSE
             ITER_EXP(ID,IS) = ABS(NINT(CFLXY))
           END IF
         END IF
         DT4AI    = DT4A/REAL(ITER_EXP(ID,IS),JWRU)
         DTSI(:)  = DT4AI/SI(:)

         U(1:MNP) = REAL(UOLD(1:MNP),JWRU)
         CALL EXCHANGE(U)
!
!  Loop over all sub time steps, all quantities in this loop depend on the solution U itself !!!
!
         DO IT = 1, ITER_EXP(ID,IS)
           ST = 0.0_JWRU ! Init. ... only used over the residual nodes see IP loop
           DO IE = 1, MNE
             NI     = INE(:,IE)
             U3     = U(NI) 
             UTILDE = N(IE)*(DOT_PRODUCT(FLALL(:,IE),U3*IOBPD(ID,NI)))
             !UTILDE = N(IE)*(DOT_PRODUCT(FLALL(:,IE),U3))
             ST(NI) = ST(NI)*IOBPD(ID,NI)+KELEM(:,IE)*(U3 - UTILDE) 
             !ST(NI) = ST(NI)+KELEM(:,IE)*(U3 - UTILDE)
           END DO
           U = MAX(0.0_JWRU,U-DTSI*ST*REAL(IOBWB,JWRU))!*REAL(IOBPD(ID,:),JWRU)
           CALL EXCHANGE(U)
         END DO  ! ----> End Iteration

         UNEW(NINF-1)=0.0_JWRB
         UNEW(1:MNP) = REAL(U(1:MNP),JWRB)


#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('EXPLICIT_N_SCHEME',1,ZHOOK_HANDLE)
#endif
      END SUBROUTINE EXPLICIT_N_SCHEME
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXCHANGE_FOR_FL1 (AC)
      USE YOWMPP   , ONLY : NINF, NSUP
      IMPLICIT NONE
      REAL(KIND=JWRB), intent(inout) :: AC(NINF-1:NSUP,NANG,NFRE)
      INTEGER(KIND=JWIM) :: IS, ID, IP
      REAL(KIND=JWRU) :: ACexch(MNP)

      DO is=1, NFRE
        DO id=1, NANG
          DO ip=1, mnp
            ACexch(ip) = REAL(AC(ip,id,is),JWRU)
          ENDDO
          CALL exchange(ACexch)
          DO ip = 1 , mnp
            AC(ip,id,is) = REAL(ACexch(ip),JWRB)
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE EXCHANGE_FOR_FL1
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPLICIT_PSI_SCHEME  ( IS, ID, UOLD, UNEW )

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  PSI (positive streamline invariant) CRD (contour approach) scheme accoridng to Struijs et al. Csik et al. and Roland, 2008
!                 2nd order in space on regular meshes in steady state, better than 1st order in space on irregular meshes. Best scheme according 
!                 to me ... 
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF*

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      USE YOWUNPOOL
      USE YOWMPP   , ONLY : NINF, NSUP
      use yowexchangemodule
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN)    :: IS,ID
      REAL(KIND=JWRB), INTENT(IN)   :: UOLD(NINF-1:NSUP)
      REAL(KIND=JWRB), INTENT(OUT)  :: UNEW(NINF-1:NSUP)
!
! local INTEGER
!
      INTEGER(KIND=JWIM) :: IP, IE, IT
      INTEGER(KIND=JWIM) :: I1, I2, I3
      INTEGER(KIND=JWIM) :: NI(3),K
!
! local double
!
      REAL(KIND=JWRU) :: FT
      REAL(KIND=JWRU) :: UTILDE
      REAL(KIND=JWRU) :: DTMAX_GLOBAL_EXP(1), DTMAX_EXP
      REAL(KIND=JWRU) :: REST
      REAL(KIND=JWRU) :: LAMBDA(2), DT4AI
      REAL(KIND=JWRU) :: FL11,FL12,FL21,FL22,FL31,FL32
      REAL(KIND=JWRU) :: KTMP(3)
      REAL(KIND=JWRU) :: KKSUM(MNP), ST(MNP), N(MNE), U3(3), ST3(3)
      REAL(KIND=JWRU) :: C(2,MNP), DTSI(MNP), CFLXY
      REAL(KIND=JWRU) :: FLALL(3,MNE), U(MNP)
      REAL(KIND=JWRU) :: FL111, FL112, FL211, FL212, FL311, FL312
      REAL(KIND=JWRU) :: KELEM(3,MNE)
      REAL(KIND=JWRU) :: THETA_L(3)
      REAL(KIND=JWRU) :: BET1(3), BETAHAT(3)
!
! local parameter
!
      REAL(KIND=JWRU) :: TMP
!
!        Calculate phase speeds for the certain spectral comp1.d0nt ...
!
      CALL CADVXY(IS,ID,C)
!
!        Calculate K-Values and contour based quantities ...
!
      DO IE = 1, MNE
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        LAMBDA(1) = ONESIXTH *(C(1,I1)+C(1,I2)+C(1,I3))
        LAMBDA(2) = ONESIXTH *(C(2,I1)+C(2,I2)+C(2,I3))
        KELEM(1,IE) = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
        KELEM(2,IE) = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
        KELEM(3,IE) = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
        KTMP  = KELEM(:,IE)
        TMP   = SUM(MIN(0.0_JWRU,KTMP))
        N(IE) = - 1.0_JWRU/MIN(-THR,TMP)
        KELEM(:,IE) = MAX(0.0_JWRU,KTMP)
        FL11  = C(1,I2) * IEN(1,IE) + C(2,I2) * IEN(2,IE)
        FL12  = C(1,I3) * IEN(1,IE) + C(2,I3) * IEN(2,IE)
        FL21  = C(1,I3) * IEN(3,IE) + C(2,I3) * IEN(4,IE)
        FL22  = C(1,I1) * IEN(3,IE) + C(2,I1) * IEN(4,IE)
        FL31  = C(1,I1) * IEN(5,IE) + C(2,I1) * IEN(6,IE)
        FL32  = C(1,I2) * IEN(5,IE) + C(2,I2) * IEN(6,IE)
        FL111 = 2.0_JWRU*FL11+FL12
        FL112 = 2.0_JWRU*FL12+FL11
        FL211 = 2.0_JWRU*FL21+FL22
        FL212 = 2.0_JWRU*FL22+FL21
        FL311 = 2.0_JWRU*FL31+FL32
        FL312 = 2.0_JWRU*FL32+FL31
        FLALL(1,IE) = FL311 + FL212
        FLALL(2,IE) = FL111 + FL312
        FLALL(3,IE) = FL211 + FL112
      END DO

      IF (LCALC) THEN ! If the current field or water level --> new CFL 
        KKSUM = 0.0_JWRU
        DO IE = 1, MNE
          NI = INE(:,IE)
          KKSUM(NI) = KKSUM(NI) + KELEM(:,IE)
        END DO
        DTMAX_GLOBAL_EXP(1) = LARGE
        DO IP = 1, MNP
          DTMAX_EXP = SI(IP)/MAX(SMALL,KKSUM(IP))
          DTMAX_GLOBAL_EXP(1) = MIN ( DTMAX_GLOBAL_EXP(1), DTMAX_EXP)
        END DO
        CFLXY = DT4A/DTMAX_GLOBAL_EXP(1)
        REST  = ABS(MOD(CFLXY,1.0_JWRU))
        IF (REST .LT. THR) THEN
          ITER_EXP(ID,IS) = ABS(NINT(CFLXY))
        ELSE IF (REST .GT. THR .AND. REST .LT. 0.5_JWRU) THEN
          ITER_EXP(ID,IS) = ABS(NINT(CFLXY)) + 1
        ELSE
          ITER_EXP(ID,IS) = ABS(NINT(CFLXY))
        END IF
      END IF

      DT4AI    = DT4A/REAL(ITER_EXP(ID,IS),JWRU)
      DTSI(:)  = DT4AI/SI(:)

      U(1:MNP) = REAL(UOLD(1:MNP),JWRU)
!
!  Loop over all sub time steps, all quantities in this loop depend on the solution U itself !!!
!
      DO IT = 1, ITER_EXP(ID,IS)
        ST = 0.0_JWRU
        DO IE = 1, MNE
          NI   = INE(:,IE)
          FT   = -ONESIXTH*DOT_PRODUCT(U(NI)*IOBPD(ID,NI),FLALL(:,IE))
          UTILDE = N(IE)*(DOT_PRODUCT(KELEM(:,IE),U(NI)*IOBPD(ID,NI))-FT)
          THETA_L(:) = KELEM(:,IE) * (U(NI) - UTILDE)
          IF (ABS(FT) .GT. 0.0_JWRU) THEN
            BET1(:) = THETA_L(:)/FT
            IF (ANY( BET1 .LT. 0.0_JWRU) ) THEN
              BETAHAT(1)= BET1(1) + 0.5_JWRU * BET1(2)
              BETAHAT(2)= BET1(2) + 0.5_JWRU * BET1(3)
              BETAHAT(3)= BET1(3) + 0.5_JWRU * BET1(1)
              BET1(1)= MAX(0.0_JWRU,MIN(BETAHAT(1),1.0_JWRU-BETAHAT(2),1.0_JWRU))
              BET1(2)= MAX(0.0_JWRU,MIN(BETAHAT(2),1.0_JWRU-BETAHAT(3),1.0_JWRU))
              BET1(3)= MAX(0.0_JWRU,MIN(BETAHAT(3),1.0_JWRU-BETAHAT(1),1.0_JWRU))
              THETA_L(:) = FT * BET1
            END IF
          ELSE
            THETA_L(:) = 0.0_JWRU
          END IF
          ST(NI) = ST(NI)*IOBPD(ID,NI) + THETA_L ! the 2nd term are the theta values of each node ...
        END DO
        U = MAX(0.0_JWRU,U-DTSI*ST*REAL(IOBWB,JWRU))!*REAL(IOBPD(ID,:),JWRU)
      END DO  ! ----> End Iteration
      CALL EXCHANGE(U)
      UNEW(NINF-1)=0.0_JWRB
      UNEW(1:MNP) = REAL(U(1:MNP),JWRB)

      END SUBROUTINE EXPLICIT_PSI_SCHEME
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPLICIT_LFPSI_SCHEME(IS,ID,UOLD,UNEW)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Law-Friedrich Flux Corrected Transport blended PSI schemes based on Hubbard & Roe (2000), Roland (2008)
!                 True 2nd order space/time scheme on narrow stencil. A bit expensive, needs GSE correction if mesh fine and directioan resolution 
!                 coarse 
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF*

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

         USE YOWUNPOOL
         USE YOWMPP   , ONLY : NINF, NSUP
         IMPLICIT NONE
         INTEGER(KIND=JWIM), INTENT(IN)    :: IS,ID
         REAL(KIND=JWRB), INTENT(IN)   :: UOLD(NINF-1:NSUP)
         REAL(KIND=JWRB), INTENT(OUT)  :: UNEW(NINF-1:NSUP)
!
! local INTEGER
!
         INTEGER(KIND=JWIM) :: IP, IE, IT
         INTEGER(KIND=JWIM) :: I1, I2, I3, K
         INTEGER(KIND=JWIM) :: NI(3)
!
! local double
!
         REAL(KIND=JWRU) :: FT
         REAL(KIND=JWRU) :: UTILDE

         REAL(KIND=JWRU) :: DTMAX_GLOBAL_EXP, DTMAX_EXP

         REAL(KIND=JWRU) :: REST
         REAL(KIND=JWRU) :: TMP(3), TMP1

         REAL(KIND=JWRU) :: LAMBDA(2), DT4AI
         REAL(KIND=JWRU) :: BET1(3), BETAHAT(3), BL

         REAL(KIND=JWRU) :: FL11,FL12,FL21,FL22,FL31,FL32

         REAL(KIND=JWRU) :: THETA_L(3,MNE),THETA_H(3),THETA_ACE(3,MNE),UTMP(3)
         REAL(KIND=JWRU) :: WII(2,MNP), UL(MNP), USTARI(2,MNP), KTMP(3)

         REAL(KIND=JWRU) :: KKSUM(MNP), ST(MNP)
         REAL(KIND=JWRU) :: PM(MNP), PP(MNP), UIM(MNP), UIP(MNP)

         REAL(KIND=JWRU) :: C(2,MNP), U(MNP), DTSI(MNP), CFLXY, N(MNE)
         REAL(KIND=JWRU) :: FL111, FL112, FL211, FL212, FL311, FL312
         REAL(KIND=JWRU) :: KELEM(3,MNE), FLALL(3,MNE)
!
! local parameter
!
         BL = 0.0_JWRU
!
!        Calculate phase speeds for the certain spectral comp1.d0nt ...
!
         CALL CADVXY(IS,ID,C)
!
!        Calculate K-Values and contour based quantities ...
!
         DO IE = 1, MNE
            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)
            LAMBDA(1) = ONESIXTH *(C(1,I1)+C(1,I2)+C(1,I3))
            LAMBDA(2) = ONESIXTH *(C(2,I1)+C(2,I2)+C(2,I3))
!            WRITE(DBG%FHNDL,*) 'LAMBDA', LAMBDA
            KELEM(1,IE) = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
            KELEM(2,IE) = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
            KELEM(3,IE) = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
!            WRITE(DBG%FHNDL,'(A10,6D15.4)') 'TEST IEN', IEN(:,IE)
            KTMP  = KELEM(:,IE)
            TMP1   = SUM(MIN(0.0_JWRU,KTMP))
            N(IE) = - 1.0_JWRU/MIN(-THR,TMP1)
            KELEM(:,IE) = MAX(0.0_JWRU,KTMP)
!            WRITE(DBG%FHNDL,'(A20,3I10,3F15.10)') 'TESTING KELEM', 
!     &                                    IS, ID, IE, KELEM(:,IE)
            FL11  = C(1,I2) * IEN(1,IE) + C(2,I2) * IEN(2,IE)
            FL12  = C(1,I3) * IEN(1,IE) + C(2,I3) * IEN(2,IE)
            FL21  = C(1,I3) * IEN(3,IE) + C(2,I3) * IEN(4,IE)
            FL22  = C(1,I1) * IEN(3,IE) + C(2,I1) * IEN(4,IE)
            FL31  = C(1,I1) * IEN(5,IE) + C(2,I1) * IEN(6,IE)
            FL32  = C(1,I2) * IEN(5,IE) + C(2,I2) * IEN(6,IE)
            FL111 = 2.0_JWRU*FL11+FL12
            FL112 = 2.0_JWRU*FL12+FL11
            FL211 = 2.0_JWRU*FL21+FL22
            FL212 = 2.0_JWRU*FL22+FL21
            FL311 = 2.0_JWRU*FL31+FL32
            FL312 = 2.0_JWRU*FL32+FL31
            FLALL(1,IE) = FL311 + FL212
            FLALL(2,IE) = FL111 + FL312
            FLALL(3,IE) = FL211 + FL112
         END DO

         IF (LCALC) THEN ! If the current field or water level --> new CFL 
           KKSUM = 0.0_JWRU
           DO IE = 1, MNE
             NI = INE(:,IE)
             KKSUM(NI) = KKSUM(NI) + KELEM(:,IE)
           END DO
           DTMAX_GLOBAL_EXP = LARGE
           DO IP = 1, MNP
             DTMAX_EXP = SI(IP)/MAX(SMALL,KKSUM(IP))
             DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
             !WRITE(DBG%FHNDL,*)'LOCAL DT',IP,SI(IP),KKSUM(IP),DTMAX_EXP
           END DO
           CFLXY = DT4A/DTMAX_GLOBAL_EXP
           REST  = ABS(MOD(CFLXY,1.0_JWRU))
           IF (REST .LT. THR) THEN
             ITER_EXP(ID,IS) = ABS(NINT(CFLXY))
           ELSE IF (REST .GT. THR .AND. REST .LT. 0.5_JWRU) THEN
             ITER_EXP(ID,IS) = ABS(NINT(CFLXY)) + 1
           ELSE
             ITER_EXP(ID,IS) = ABS(NINT(CFLXY))
           END IF
!           WRITE(DBG%FHNDL,'(2I5,3F15.8,I10)') IS, ID, DT4A, CFLXY, 
!     &                            DTMAX_GLOBAL_EXP, ITER_EXP(IS,ID)
         END IF

         DT4AI    = REAL(DT4A,JWRU)/REAL(ITER_EXP(ID,IS),JWRU)
         DTSI(:)  = DT4AI/SI(:)

         U(1:MNP)  = REAL(UOLD(1:MNP),JWRU)
!
!  Loop over all sub time steps, all quantities in this loop depend on the solution U itself !!!
!
         DO IT = 1, ITER_EXP(ID,IS)
!
! Element loop
!
           ST = 0.0_JWRU
           PM = 0.0_JWRU
           PP = 0.0_JWRU
           DO IE = 1, MNE
             NI   = INE(:,IE)
             FT   = -ONESIXTH*DOT_PRODUCT(U(NI)*IOBPD(ID,NI),FLALL(:,IE))
             UTILDE = N(IE)*(DOT_PRODUCT(KELEM(:,IE),U(NI)*IOBPD(ID,NI))-FT)
             THETA_L(:,IE) = KELEM(:,IE) * (U(NI) - UTILDE)
             IF (ABS(FT) .GT. 0.0_JWRU) THEN
               BET1(:) = THETA_L(:,IE)/FT
               IF (ANY( BET1 .LT. 0.0_JWRU) ) THEN
                 BETAHAT(1)= BET1(1) + 0.5_JWRU * BET1(2)
                 BETAHAT(2)= BET1(2) + 0.5_JWRU * BET1(3)
                 BETAHAT(3)= BET1(3) + 0.5_JWRU * BET1(1)
                 BET1(1)= MAX(0.0_JWRU,MIN(BETAHAT(1),1.0_JWRU-BETAHAT(2),1.0_JWRU))
                 BET1(2)= MAX(0.0_JWRU,MIN(BETAHAT(2),1.0_JWRU-BETAHAT(3),1.0_JWRU))
                 BET1(3)= MAX(0.0_JWRU,MIN(BETAHAT(3),1.0_JWRU-BETAHAT(1),1.0_JWRU))
                 THETA_L(:,IE) = FT * BET1
               END IF
             ELSE
               THETA_L(:,IE) = 0.0_JWRU
             END IF
              ST(NI)  = ST(NI) + THETA_L(:,IE)
              !THETA_H = (ONETHIRD+DT4AI/(2.*TRIA(IE)) * KELEM(:,IE) )*FT ! LAX
              THETA_H = (ONETHIRD+TWOTHIRD*KELEM(:,IE)/SUM(MAX(0.0_JWRU,KELEM(:,IE))))*FT  ! CENTRAL
              THETA_ACE(:,IE) = THETA_H-THETA_L(:,IE)
              PP(NI) =  PP(NI) + MAX( 0.0_JWRU, -THETA_ACE(:,IE)) * DTSI(NI)
              PM(NI) =  PM(NI) + MIN( 0.0_JWRU, -THETA_ACE(:,IE)) * DTSI(NI)
            END DO
!
            !U = MAX(0.0_JWRU,U-DTSI*ST)*REAL(IOBPD(ID,:),JWRU)
            UL = MAX(0.0_JWRU,U-DTSI*ST*REAL(IOBWB,JWRU))

            USTARI(1,:) = MAX(UL,U) 
            USTARI(2,:) = MIN(UL,U) 

            UIP = 0._JWRU
            UIM = 0._JWRU
            DO IE = 1, MNE
              NI = INE(:,IE)
              UIP(NI) = MAX (UIP(NI), MAXVAL( USTARI(1,NI) ))
              UIM(NI) = MIN (UIM(NI), MINVAL( USTARI(2,NI) ))
            END DO

            WII(1,:) = MIN(1.0_JWRU,(UIP-UL)/MAX( REAL(SMALL,JWRU),PP))
            WII(2,:) = MIN(1.0_JWRU,(UIM-UL)/MIN(-REAL(SMALL,JWRU),PM))

            ST = 0.0_JWRU
            DO IE = 1, MNE
               I1 = INE(1,IE)
               I2 = INE(2,IE)
               I3 = INE(3,IE)
               IF (THETA_ACE(1,IE) .LT. 0.0_JWRU) THEN
                 TMP(1) = WII(1,I1)
               ELSE
                 TMP(1) = WII(2,I1)
               END IF
               IF (THETA_ACE(2,IE) .LT. 0.0_JWRU) THEN
                 TMP(2) = WII(1,I2)
               ELSE
                 TMP(2) = WII(2,I2)
               END IF
               IF (THETA_ACE(3,IE) .LT. 0.0_JWRU) THEN
                 TMP(3) = WII(1,I3)
               ELSE
                 TMP(3) = WII(2,I3)
               END IF
               TMP1 = MINVAL(TMP)
               ST(I1) = ST(I1) + THETA_ACE(1,IE) * TMP1! * (1.d0 - BL) + BL * THETA_L(1,IE)
               ST(I2) = ST(I2) + THETA_ACE(2,IE) * TMP1! * (1.d0 - BL) + BL * THETA_L(2,IE)
               ST(I3) = ST(I3) + THETA_ACE(3,IE) * TMP1! * (1.d0 - BL) + BL * THETA_L(3,IE)
            END DO
            U = MAX(0.0_JWRU,UL-DTSI*ST*REAL(IOBWB,JWRU))!*REAL(IOBPD(ID,:),JWRU)
         END DO  ! ----> End Iteration
         CALL EXCHANGE(U)
         UNEW(NINF-1)=0.0_JWRB
         UNEW(1:MNP) = REAL(U(1:MNP),JWRB)

      END SUBROUTINE EXPLICIT_LFPSI_SCHEME
!**********************************************************************
!*
!**********************************************************************
      SUBROUTINE CADVXY(IS,ID,C)
      USE YOWUNPOOL, ONLY : LADVTEST, LCUR, LSPHE, CG, DEGRAD, REARTH
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN)  :: IS, ID
      REAL(KIND=JWRU), INTENT(OUT)  :: C(2,MNP)
      INTEGER(KIND=JWIM) :: IP
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRU) :: eCX, eCY

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('CADVXY',0,ZHOOK_HANDLE)
#endif
      IF (LCUR) THEN
        IF (LSPHE) THEN
          DO IP=1,MNP
            eCX = CG(IS,IP)*SINTH(ID)+CURTXY(1,IP)
            eCY = CG(IS,IP)*COSTH(ID)+CURTXY(2,IP)
            C(1,IP) = eCX * INVTRANS1(IP)
            C(2,IP) = eCY * DGRTHM1
          END DO
        ELSE
          DO IP=1,MNP
            C(1,IP) = CG(IS,IP)*SINTH(ID)+CURTXY(1,IP)
            C(2,IP) = CG(IS,IP)*COSTH(ID)+CURTXY(2,IP)
          END DO
        END IF
      ELSE
        IF (LSPHE) THEN
          DO IP=1,MNP
            eCX = CG(IS,IP)*SINTH(ID)
            eCY = CG(IS,IP)*COSTH(ID)
            C(1,IP) = eCX * INVTRANS1(IP)
            C(2,IP) = eCY * DGRTHM1
          END DO
        ELSE
          DO IP=1,MNP
            C(1,IP) = CG(IS,IP)*SINTH(ID)
            C(2,IP) = CG(IS,IP)*COSTH(ID)
          END DO
        END IF
      END IF
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('CADVXY',1,ZHOOK_HANDLE)
#endif
      END SUBROUTINE CADVXY
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_UNWAM_ARRAYS

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Intialize UnWam arrays 
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF*

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      USE YOWUNPOOL 
      IMPLICIT NONE
      INTEGER(KIND=JWIM) :: istat
      IF(ALLOCATED(DEPTH_INPUT)) DEALLOCATE(DEPTH_INPUT)
      ALLOCATE( DEPTH_INPUT(MNP))
      ALLOCATE( CCON(MNP) ); CCON = 0
      ALLOCATE( Jstart(MNP) ); Jstart = 0
      ALLOCATE( SI(MNP) ); SI = 0.0_JWRU
      ALLOCATE( ITER_EXP(NANG,NFRE) ); ITER_EXP = 0
      ALLOCATE( TRIA(MNE) ); TRIA = 0.0_JWRU
      IF (LCFL) THEN
        ALLOCATE(CFLCXY(3,MNP)); CFLCXY = 0.0_JWRU
      END IF
      ALLOCATE( IEN(6,MNE) ); IEN = 0
      ALLOCATE( CG(NFRE,MNP) ); CG = 0.0_JWRU
      ALLOCATE( WK(NFRE,MNP) ); CG = 0.0_JWRU
      IF (IREFRA .ne. 0) THEN
        allocate(DCUX(2,MNP), DCUY(2,MNP), DDEP(2,MNP), DEPDT(MNP), stat=istat)
      END IF
      ALLOCATE(CURTXY(2,MNP), stat=istat) ; CURTXY=0.0_JWRU
      END SUBROUTINE INIT_UNWAM_ARRAYS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_UNWAM

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Initalize Unwam ... 
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF*

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------
      USE yowMpiModule
      USE YOWUNPOOL, ONLY : DEGRAD, REARTH
      USE YOWUNPOOL, ONLY : BND, DBG, GRID, FILEDEF
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, KIJL4CHNK, IJFROMCHNK
      USE YOWSHAL  , ONLY : WVENVI, DEPTHA    ,TOOSHALLOW
      USE YOWMAP   , ONLY : NIBLO
      USE yowpd, only: z, comm, np_global, initPD, setDimSize
      USE UNSTRUCT_BOUND, ONLY : INIT_BOUNDARY
      USE YOWFRED, ONLY : DELTH
      USE YOWSTAT  , ONLY : IDELPRO
      IMPLICIT NONE
      INTEGER(KIND=JWIM) :: IERR, IP, ICHNK, IPRM
      INTEGER(KIND=JWIM) :: eNext, ePrev, ID, istat
      CALL setDimSize(NFRE, NANG)
      CALL SET_UNWAM_HANDLES ! set file handles
      CALL initPD("system.dat")
      CALL mpi_barrier(comm, ierr)
      NIBLO = np_global
      CALL INIT_UNWAM_ARRAYS
!JB do not allow negative depth !!!!
      DO IP = 1, MNP 
        DEPTH_INPUT(IP) = REAL(MAX(DEP(IP),REAL(TOOSHALLOW,JWRU)),JWRB)
      ENDDO


!!    Transfer the water depth to array used by the wave physics
      DO ICHNK = 1, NCHNK
        DO IPRM = 1, KIJL4CHNK(ICHNK)
          IP = IJFROMCHNK(IPRM, ICHNK)
          WVENVI%DEPTH(IPRM, ICHNK) = DEPTH_INPUT(IP)
        ENDDO
      ENDDO

      DDIR = DELTH
!
! set time step
!
      allocate(INVTRANS1(MNP))
      DGRTH = DEGRAD*REARTH
      DGRTHM1 = 1.0_JWRU/(DGRTH)
      DO IP=1,MNP
        INVTRANS1(IP) = 1.0_JWRU/(DGRTH*COS(MIN(ABS(YP(IP)),YPMAX)*DEGRAD))
      END DO
!
! set time step
!
      DT4A = REAL(IDELPRO,JWRU)
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'DT4A=', DT4A
      FLUSH(740+MyRankGlobal)
#endif

      DT4D = DT4A
      DT4F = DT4A
!
! set ID_PREVNEXT
!
      allocate(ID_PREV(NANG), ID_NEXT(NANG), stat=istat)
      DO ID=1,NANG
        IF (ID.eq.1) THEN
          ePrev=NANG
        ELSE
          ePrev=ID-1
        END IF
        IF (ID.eq.NANG) THEN
          eNext=1
        ELSE
          eNext=ID+1
        END IF
        ID_PREV(ID)=ePrev
        ID_NEXT(ID)=eNext
      END DO
      CALL INIT_FLUCT        ! init fluctuation splitting stuff
      CALL INIT_BOUNDARY
      END SUBROUTINE INIT_UNWAM
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CORRECT_SINGLE_DXP(DXP)
      IMPLICIT NONE
      REAL(KIND=JWRU), INTENT(INOUT) :: DXP
      IF (DXP .LE. -180.0_JWRU) THEN
        DXP=DXP + 360.0_JWRU
      END IF
      IF (DXP .GE. 180.0_JWRU) THEN
        DXP=DXP - 360.0_JWRU
      END IF
      END SUBROUTINE CORRECT_SINGLE_DXP
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_FLUCT
      USE YOWUNPOOL
      USE YOWPD,     ONLY: COMM, RANK, MYRANK
      use YOWERROR,  ONLY: PARALLEL_ABORT
      USE YOWFRED  , ONLY : FR
      USE YOWTEST  , ONLY : IU06
      USE YOWSPHERE, ONLY : SPHERICAL_COORDINATE_DISTANCE, &
                          & SPHERICAL_COORDINATE_AREA
      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: IERR, IPGL
      INTEGER(KIND=JWIM) :: I, J, K
      INTEGER(KIND=JWIM) :: IP, IE, IS, ISTAT
      INTEGER(KIND=JWIM) :: POS, POS_J, POS_K
      INTEGER(KIND=JWIM) :: IP_I, IP_J, IP_K
      INTEGER(KIND=JWIM) :: I1, I2, I3, NI(3)
      INTEGER(KIND=JWIM) :: CHILF(MNP), COUNT_MAX
      INTEGER(KIND=JWIM) :: ITMP(MNP)
      INTEGER(KIND=JWIM) :: TMPINE
      INTEGER(KIND=JWIM) :: INEXT, IP_NEXT, IE_ADJ, nbMatch, ICON
      INTEGER(KIND=JWIM) :: IADJ, POS_PREV, IP_ADJ_PREV
      INTEGER(KIND=JWIM) :: IE2, POS_NEXT, IP_ADJ_NEXT
      INTEGER(KIND=JWIM) :: MAX_DEG
      INTEGER(KIND=JWIM) :: POS_TRICK(3,2)
      INTEGER(KIND=JWIM), ALLOCATABLE :: PTABLE(:,:)

      REAL(KIND=JWRU) :: P1(2), P2(2), P3(2)
      REAL(KIND=JWRU) :: R1(2), R2(2), R3(2)
      REAL(KIND=JWRU) :: N1(2), N2(2), N3(2)
      REAL(KIND=JWRU) :: TRIA03
      REAL(KIND=JWRU) :: DXP1, DXP2, DXP3 
      REAL(KIND=JWRU) :: DYP1, DYP2, DYP3
      REAL(KIND=JWRU) :: DBLTMP, WN, WVC, WVK, WVCG
      REAL(KIND=JWRU) :: TL1, TL2, TL3
      REAL(KIND=JWRU) :: TLMIN,TMPTLMIN
      REAL(KIND=JWRU) :: TLMAX,TMPTLMAX
      REAL(KIND=JWRU) :: AVETA
      REAL(KIND=JWRU) :: PROV1, PROV2, PROV3
      REAL(KIND=JWRU) :: AREA, AREA_RAD

      LOGICAL   :: LWRONG = .FALSE.
      LOGICAL   :: LSINEG = .FALSE.
      LOGICAL   :: CHECK_COMBIN_ORIENT

      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

      TLMIN = LARGE
      TLMAX = SMALL 
      AVETA = 0.0_JWRU
!
!     CALCULATE NORMAL VECTORS
! 
      DO IE = 1, MNE
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        DXP1 = XP(I2) - XP(I1) ! dx dy for each element and node ...
        DYP1 = YP(I2) - YP(I1)
        DXP2 = XP(I3) - XP(I2)
        DYP2 = YP(I3) - YP(I2)
        DXP3 = XP(I1) - XP(I3)
        DYP3 = YP(I1) - YP(I3)
        IF (APPLY_DXP_CORR) THEN
          CALL CORRECT_SINGLE_DXP(DXP1)
          CALL CORRECT_SINGLE_DXP(DXP2)
          CALL CORRECT_SINGLE_DXP(DXP3)
        END IF
        IEN(1,IE) = - DYP2
        IEN(2,IE) =   DXP2
        IEN(3,IE) = - DYP3
        IEN(4,IE) =   DXP3
        IEN(5,IE) = - DYP1
        IEN(6,IE) =   DXP1
        DBLTMP = (DXP3*DYP1 - DYP3*DXP1)*0.5_JWRU
        IF (LSPHE .and. USE_EXACT_FORMULA_SPHERICAL_AREA) THEN
          CALL SPHERICAL_COORDINATE_AREA(XP(I1), XP(I2), XP(I3),    &
     &           YP(I1), YP(I2), YP(I3), AREA_RAD)
          AREA=AREA_RAD*RADDEG*RADDEG
        ELSE
          AREA=DBLTMP
        END IF
        TRIA(IE) = AREA
        IF (TRIA(IE) .LT. SMALL) THEN
          TMPINE = INE(2,IE)
          INE(2,IE) = INE(3,IE)
          INE(3,IE) = TMPINE
          I2 = INE(2,IE)
          I3 = INE(3,IE)
          TRIA(IE) = -TRIA(IE)
          PROV1=IEN(6,IE) ! DXP1
          PROV2=IEN(2,IE) ! DXP2
          PROV3=IEN(4,IE) ! DXP3
          IEN(6,IE)=-PROV3
          IEN(2,IE)=-PROV2
          IEN(4,IE)=-PROV1
          PROV1= - IEN(5,IE)
          PROV2= - IEN(1,IE)
          PROV3= - IEN(3,IE)
          IEN(1,IE) = PROV2
          IEN(3,IE) = PROV1
          IEN(5,IE) = PROV3
          LWRONG = .TRUE.     ! check element orientation and swap if needed ...
        END IF
        TL1 = SQRT(IEN(5,IE)**2 + IEN(6,IE)**2)  ! Edge length's 
        TL2 = SQRT(IEN(3,IE)**2 + IEN(4,IE)**2)
        TL3 = SQRT(IEN(1,IE)**2 + IEN(2,IE)**2)
        TMPTLMIN = MIN(TL1, TL2, TL3)                ! max. min. edge lengts of hte whole mesh 
        TMPTLMAX = MAX(TL1, TL2, TL3)
        IF (TLMIN > TMPTLMIN) THEN
          TLMIN = TMPTLMIN
        END IF
        IF (TLMAX < TMPTLMAX) THEN
          TLMAX = TMPTLMAX
        END IF
        AVETA = AVETA+TRIA(IE)
      END DO
      AVETA = AVETA / MNE ! Average edge area ...
!
! Calculate the max. number of connected elements count_max
!
      SI(:)   = 0.0_JWRU ! Median Dual Patch Area of each Node
      CCON(:) = 0        ! Number of connected Elements
      DO IE = 1 , MNE
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        CCON(I1) = CCON(I1) + 1
        CCON(I2) = CCON(I2) + 1
        CCON(I3) = CCON(I3) + 1
        TRIA03 = ONETHIRD * TRIA(IE)
        SI(I1) = SI(I1) + TRIA03
        SI(I2) = SI(I2) + TRIA03
        SI(I3) = SI(I3) + TRIA03
      ENDDO

      CALL exchange(SI)

      DO IP =1, MNP
        IF(SI(IP) .LT. SMALL) THEN
          write(IU06,*) 'INIT_FLUCT : SI(IP) NEG. OR LESS THAN THR. IP, SI(IP) =', IP, SI(IP)
          IPGL=RANK(MYRANK+1)%IPLG(IP)
          write(IU06,*) 'INIT_FLUCT IP, IPGL, XP, YP, MYRANK',IP, IPGL, XP(IP), YP(IP), MYRANK
          write(IU06,*) 'INIT_FLUCT IP, IPGL, MYRANK',IP, IPGL, MYRANK
          write(DBG%FHNDL,*) 'SI(IP) NEG. OR LESS THAN THR. IP=', IP
          LSINEG=.TRUE.
        ENDIF
      END DO
      CALL MPI_BARRIER(COMM, IERR)
      IF(LSINEG) THEN
        CALL FLUSH(IU06)
        call parallel_abort('INIT_FLUCY: SI NEGATIVE OR LESS THAN THRESHOLD')
      ENDIF

      MAXMNECON  = MAXVAL(CCON)
      ALLOCATE(CELLVERTEX(MNP,MAXMNECON,2), stat=istat)
      IF (istat/=0) stop 'wwm_fluctsplit, allocate error 4'
      CELLVERTEX(:,:,:) = 0 ! Stores for each node the Elementnumbers of the connected Elements
      CHILF             = 0 ! and the Position of the position of the Node in the Element Index
      DO IE = 1, MNE
        DO J=1,3
          I = INE(J,IE)
          CHILF(I) = CHILF(I)+1
          CELLVERTEX(I,CHILF(I),1) = IE
          CELLVERTEX(I,CHILF(I),2) = J
        END DO
      ENDDO
!
!        Emulates loop structure and counts max. entries in the different pointers that have to be designed
!
      J = 0
      DO IP = 1, MNP
        DO I = 1, CCON(IP)
          J = J + 1
        END DO
      END DO
      COUNT_MAX = J ! Max. Number of entries in the pointers used in the calculations
#ifdef DEBUG
      IF (COUNT_MAX.ne.3*MNE) THEN
        Print *, 'COUNT_MAX=', COUNT_MAX
        Print *, 'MNE=', MNE
        STOP 'Do Not Sleep Before solving the problem'
      ENDIF
#endif
      
      ALLOCATE (IE_CELL(COUNT_MAX), POS_CELL(COUNT_MAX), stat=istat)
      IF (istat/=0) STOP 'MEMORY ERROR AT unwam.F 1210'
      ALLOCATE (IE_CELL2(MNP,MAXMNECON), POS_CELL2(MNP,MAXMNECON), stat=istat)
      IF (istat/=0) STOP 'MEMORY ERROR AT unwam.F 1212'

      IE_CELL  = 0
      POS_CELL = 0
      IE_CELL2  = 0
      POS_CELL2 = 0

      J = 0
      DO IP = 1, MNP
        DO I = 1, CCON(IP)
          J = J + 1
          IE_CELL(J)      = CELLVERTEX(IP,I,1)
          POS_CELL(J)     = CELLVERTEX(IP,I,2)
          IE_CELL2(IP,I)  = CELLVERTEX(IP,I,1)
          POS_CELL2(IP,I) = CELLVERTEX(IP,I,2)
        END DO
      END DO
      DEALLOCATE(CELLVERTEX)

      CHECK_COMBIN_ORIENT=.TRUE.
      IF (CHECK_COMBIN_ORIENT) THEN
        DO IE=1,MNE
          DO I=1,3
            INEXT=POS_TRICK(I,1)
            IP=INE(I, IE)
            IP_NEXT=INE(I, IE)
            nbMatch=0
            IE_ADJ=-1
            DO ICON=1,CCON(IP)
              IE2=IE_CELL2(IP,ICON)
              IF (IE .ne. IE2) THEN
                POS=POS_CELL2(IP, ICON)
                POS_NEXT=POS_TRICK(POS,1)
                IP_ADJ_NEXT=INE(POS_NEXT,IE2)
#ifdef DEBUG
                IF (IP_ADJ_NEXT .eq. IP_NEXT) THEN
                  Print *, 'Combinatorial orientability problem'
                  Print *, 'IE=', IE, ' IE2=', IE2
                  Print *, 'IP=', IP, ' IP_NEXT=', IP_NEXT
                  STOP
                END IF
#endif
                POS_PREV=POS_TRICK(POS,2)
                IP_ADJ_PREV=INE(POS_PREV,IE2)
                IF (IP_ADJ_PREV .eq. IP_NEXT) THEN
                  nbMatch=nbMatch+1
                  IE_ADJ=IE2
                END IF
              END IF
            END DO
#ifdef DEBUG
            IF (nbMatch .gt. 1) THEN
              Print *, 'nbMatch is too large.'
              Print *, 'Should be 0 for boundary edge'
              Print *, 'Should be 1 for interior edges'
              Print *, 'nbMatch=', nbMatch
              STOP
            END IF
#endif           
          END DO
        END DO
      END IF


      IF (LIMPLICIT) THEN
        ALLOCATE(PTABLE(COUNT_MAX,7), JA_IE(3,3,MNE), stat=istat)
        IF (istat/=0) stop 'wwm_fluctsplit, allocate error 6'
        PTABLE(:,:) = 0 ! Table storing some other values needed to design the sparse matrix pointers.
        J = 0
        DO IP = 1, MNP
          DO I = 1, CCON(IP)
            J = J + 1
            IE    = IE_CELL(J)
            POS   = POS_CELL(J)
            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)
            IF (POS == 1) THEN
              POS_J = 2
              POS_K = 3
            ELSE IF (POS == 2) THEN
              POS_J = 3
              POS_K = 1
            ELSE
              POS_J = 1
              POS_K = 2
            END IF
            IP_I = IP
            IP_J = INE(POS_J,IE)
            IP_K = INE(POS_K,IE)
            PTABLE(J,1) = IP_I ! Node numbers of the connected elements
            PTABLE(J,2) = IP_J
            PTABLE(J,3) = IP_K
            PTABLE(J,4) = POS  ! Position of the nodes in the element index
            PTABLE(J,5) = POS_J
            PTABLE(J,6) = POS_K
            PTABLE(J,7) = IE   ! Element numbers same as IE_CELL
          END DO
        END DO
!
! Count number of nonzero entries in the matrix ...
! Basically, each connected element may have two off-diagonal
! contribution and one diagonal related to the connected vertex itself ...
!
        J = 0
        NNZ = 0
        DO IP = 1, MNP
          ITMP(:) = 0
          DO I = 1, CCON(IP)
            J = J + 1
            IP_J  = PTABLE(J,2)
            IP_K  = PTABLE(J,3)
            ITMP(IP)   = 1
            ITMP(IP_J) = 1
            ITMP(IP_K) = 1
          END DO
          NNZ = NNZ + SUM(ITMP)
        END DO
!
! Allocate sparse matrix pointers using the Compressed Sparse Row Format CSR ... this is now done only of MNP nodes
! The next step is to do it for the whole Matrix MNP * MSC * MDC
! see ...:x
!
! JA Pointer according to the convention in my thesis see p. 123
! IA Pointer according to the convention in my thesis see p. 123
        ALLOCATE (JA(NNZ), IA(MNP+1), POSI(3,COUNT_MAX), I_DIAG(MNP), ASPAR_JAC(NANG,NFRE,NNZ), stat=istat)
        IF (istat/=0) stop 'unwam, allocate error 6'
        IF (.NOT. BLOCK_GAUSS_SEIDEL) THEN
          ALLOCATE (U_JACOBI(NANG,NFRE,MNP), stat=istat)
          IF (istat/=0) stop 'unwam, allocate error for U_JACOBI'
        END IF
        IF ((.NOT. LNONL) .AND. SOURCE_IMPL) THEN
          ALLOCATE(B_JAC(NANG, NFRE, MNP), stat=istat)
          IF (istat/=0) stop 'unwam, allocate B_JAC'
        END IF
        IF (REFRA_METHOD .eq. 2) THEN
          allocate(CAD_THE(MDC,MSC,MNP), CAS_SIG(MSC,MDC,MNP), stat=istat)
          IF (istat/=0) stop 'unwam, allocate CAD_THE and CAS_SIG'
        END IF
        JA = 0
        IA = 0
        POSI = 0
! Points to the position of the matrix entry in the mass matrix
! according to the CSR matrix format see p. 124
        J = 0
        K = 0
        IA(1) = 1
        MAX_DEG=0
        DO IP = 1, MNP ! Run through all rows
          ITMP=0
          DO I = 1, CCON(IP) ! Check how many entries there are ...
            J = J + 1
            IP_J  = PTABLE(J,2)
            IP_K  = PTABLE(J,3)
            ITMP(IP)   = 1
            ITMP(IP_J) = 1
            ITMP(IP_K) = 1
          END DO
          IADJ=0
          DO I = 1, MNP ! Run through all columns
            IF (ITMP(I) .GT. 0) THEN
              K = K + 1
              IF (I .ne. IP) THEN
                IADJ=IADJ + 1
              END IF
              JA(K) = I
            END IF
          END DO
          IF (IADJ .gt. MAX_DEG) THEN
            MAX_DEG=IADJ
          END IF
          IA(IP + 1) = K + 1
        END DO


        J = 0
        DO IP = 1, MNP
          DO I = 1, CCON(IP)
            J = J + 1
            IP_J  = PTABLE(J,2)
            IP_K  = PTABLE(J,3)
            DO K = IA(IP), IA(IP+1) - 1
              IF (IP   == JA(K)) POSI(1,J)  = K
              IF (IP   == JA(K)) I_DIAG(IP) = K
              IF (IP_J == JA(K)) POSI(2,J)  = K
              IF (IP_K == JA(K)) POSI(3,J)  = K
            END DO
          END DO
        END DO

        J=0
        DO IP=1,MNP
          DO I = 1, CCON(IP)
            J=J+1
            IE    =  IE_CELL(J)
            POS   =  POS_CELL(J)
            I1    =  POSI(1,J)
            I2    =  POSI(2,J)
            I3    =  POSI(3,J)
            JA_IE(POS,1,IE) = I1
            JA_IE(POS,2,IE) = I2
            JA_IE(POS,3,IE) = I3
          END DO
        END DO
        !
        ! Arrays for Jacobi-Gauss-Seidel solver
        !
        DEALLOCATE(PTABLE)
     ENDIF
!
!     Group velocities
!
      DO IP = 1, MNP 
        DO IS = 1, NFRE
          CALL WAVEKCG8(DEP(IP), FR(IS)*PI2, WN, WVC, WVK, WVCG)
          CG(IS,IP) = WVCG
          WK(IS,IP) = WVK
        END DO 
      END DO
      CALL INIT_SPECTRAL_ARR
!
!     Field used for refractions and/or freq-shift
!
      IF (IREFRA .ne. 0) THEN
        allocate(DDEP(2,MNP), DCUX(2,MNP), DCUY(2,MNP), DEPDT(MNP), stat=istat)
      END IF
      END SUBROUTINE INIT_FLUCT
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_SPECTRAL_ARR
      USE YOWPCONS , ONLY : ZPI, PI
      USE YOWFRED  , ONLY : FR, WETAIL, FRATIO
      IMPLICIT NONE
      INTEGER(KIND=JWIM) :: IS, ID
      INTEGER(KIND=JWIM) :: istat
      allocate(DS_BAND(0:NFRE+1), SPSIG(NFRE), stat=istat)
      SPSIG = FR * ZPI
      PTAIL5 = ONE / (REAL(FRATIO,JWRU))**5
      DS_BAND(0)   = SPSIG(2) - SPSIG(1)
      DS_BAND(1)   = DS_BAND(0)
      DS_BAND(NFRE) = SPSIG(NFRE) - SPSIG(NFRE-1)
      DS_BAND(NFRE+1) = DS_BAND(NFRE)
      DO IS=2,NFRE-1
        DS_BAND(IS) = (SPSIG(IS)-SPSIG(IS-1))/2. + (SPSIG(IS+1)-SPSIG(IS))/2.
      END DO
      !
      allocate(COS2TH(NANG), SIN2TH(NANG), SINCOSTH(NANG), stat=istat)
      DO ID=1,NANG
        COS2TH(ID) = COSTH(ID) * COSTH(ID)
        SIN2TH(ID) = SINTH(ID) * SINTH(ID)
        SINCOSTH(ID) = SINTH(ID) * COSTH(ID)
      END DO
      END SUBROUTINE INIT_SPECTRAL_ARR
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CHECK_SYSTEM

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Check system size and allocate 
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

        USE YOWUNPOOL, ONLY : GRID, FILEDEF
        USE YOWMAP   , ONLY : NIBLO
        USE UNSTRUCT_BOUND, ONLY : RHEADER_NODE, RHEADER_ELEMENT
        IMPLICIT NONE

        INTEGER(KIND=JWIM) :: NBOUND, NDOMAIN
        INTEGER(KIND=JWIM) :: I, ITMP, NGESAMT

        OPEN(GRID%FHNDL, FILE = GRID%FNAME, STATUS='old')
        CALL RHEADER_NODE(GRID%FHNDL,NBOUND,NDOMAIN)
        MNP = NBOUND + NDOMAIN 
        NIBLO = MNP
        DO I = 1, MNP
          READ(GRID%FHNDL, *) 
        END DO 
        CALL RHEADER_ELEMENT(GRID%FHNDL,MNE)
        CLOSE(GRID%FHNDL)  
      END SUBROUTINE CHECK_SYSTEM
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE UNWAM_OUT(IHANDLE)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Write everything out to disk ...

!     Externals.  
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

        USE YOWUNPOOL
        IMPLICIT NONE
        INTEGER(KIND=JWIM), INTENT(IN) :: IHANDLE
        WRITE(IHANDLE) MNP, MNE, NFRE, NANG
        WRITE(IHANDLE) XP
        WRITE(IHANDLE) YP 
        WRITE(IHANDLE) DEPTH_INPUT
        WRITE(IHANDLE) CCON
        WRITE(IHANDLE) SI
        WRITE(IHANDLE) TRIA
        WRITE(IHANDLE) INE
        WRITE(IHANDLE) IEN
        WRITE(IHANDLE) CG
        WRITE(IHANDLE) IOBP
        WRITE(IHANDLE) IOBPD
        WRITE(IHANDLE) IOBWB
      END SUBROUTINE UNWAM_OUT
!**********************************************************************
!!*                                                                   *
!**********************************************************************
      SUBROUTINE UNWAM_IN(IHANDLE)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Read everything from disk ...

!     Externals.  
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

        USE YOWUNPOOL
        USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, KIJL4CHNK, IJFROMCHNK
        USE YOWSHAL  , ONLY : WVENVI
        IMPLICIT NONE

        INTEGER(KIND=JWIM), INTENT(IN) :: IHANDLE

        INTEGER(KIND=JWIM) :: IP, ICHNK, IPRM
! 
!         READ(IHANDLE) MNP, MNE, MSC, MDC

        CALL INIT_UNWAM_ARRAYS
!         READ(IHANDLE) XP
!         READ(IHANDLE) YP 
        READ(IHANDLE) DEPTH_INPUT
        READ(IHANDLE) CCON
        READ(IHANDLE) SI
        READ(IHANDLE) TRIA
!         READ(IHANDLE) INE
        READ(IHANDLE) IEN
        READ(IHANDLE) CG
        READ(IHANDLE) IOBP
        READ(IHANDLE) IOBPD
        READ(IHANDLE) IOBWB

!!      Transfer the water depth to array it is used by the wave physics
        DO ICHNK = 1, NCHNK
          DO IPRM = 1, KIJL4CHNK(ICHNK)
            IP = IJFROMCHNK(IPRM, ICHNK)
            WVENVI%DEPTH(IPRM, ICHNK) = DEPTH_INPUT(IP)
          ENDDO
        ENDDO

      END SUBROUTINE UNWAM_IN
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_UNWAM_HANDLES

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Set file handles using ECMWF tool ....

!     Externals.  
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

        USE YOWUNPOOL
        USE YOWTEST  , ONLY : IU06
        USE YOWMPP   , ONLY : IRANK, NPROC
        IMPLICIT NONE

        INTEGER(KIND=JWIM) :: IWAM_GET_UNIT

!        GRID%FNAME = 'system.dat'
!        GRID%FHNDL = IWAM_GET_UNIT(IU06, GRID%FNAME, 'r', 'f', 0) 

        BND%FNAME = 'sysbnd.dat'
!        BND%FHNDL = IWAM_GET_UNIT(IU06, BND%FNAME, 'w', 'f', 0)

        DBG%FNAME = 'unwamdbg.%p.dat'
        CALL EXPAND_STRING(IRANK,NPROC,0,0,DBG%FNAME,1)
        DBG%FHNDL = IWAM_GET_UNIT(IU06, DBG%FNAME, 'w', 'f', 0, 'READWRITE')

        IF(LLUNBINOUT) THEN
          XFN_HS%FNAME = 'erghs.bin'
          XFN_HS%FHNDL = IWAM_GET_UNIT(IU06, XFN_HS%FNAME, 'w', 'u', 0, 'READWRITE')

          XFN_TM%FNAME = 'ergtm.bin'
          XFN_TM%FHNDL = IWAM_GET_UNIT(IU06, XFN_TM%FNAME, 'w', 'u', 0, 'READWRITE')

          XFN_TEST%FNAME = 'ergtest.bin'
          XFN_TEST%FHNDL = IWAM_GET_UNIT(IU06, XFN_TEST%FNAME, 'w', 'u', 0, 'READWRITE')
        ENDIF

      END SUBROUTINE SET_UNWAM_HANDLES
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVEKCG8(DEP, SIGIN, WN, WVC, WVK, WVCG)
!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------
!        Explicit arguments :  
!        --------------------   
!        Implicit arguments :     N1.d0
!        --------------------
!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Estimate group vel. 
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------
         USE YOWUNPOOL, ONLY : SMALL
         USE YOWPCONS , ONLY : G

         IMPLICIT NONE

         REAL(KIND=JWRU), INTENT(IN)  :: DEP
         REAL(KIND=JWRU), INTENT(IN)  :: SIGIN
         REAL(KIND=JWRU), INTENT(OUT) :: WVC, WVK, WVCG, WN
         REAL(KIND=JWRU) :: SGDLS , AUX1, AUX2
         REAL(KIND=JWRU) :: WKDEP, DMIN
! 
!AR: I put the mindepth for kcg to 0.1d0
!
         DMIN = 0.1_JWRU

         IF (SIGIN .LT. SMALL) THEN
            WN = 0.0_JWRU
            WVK=10.0_JWRU
            WVCG=0.0_JWRU
            RETURN
         END IF

         IF (DEP > DMIN) THEN
            SGDLS = SIGIN*SIGIN*DEP/G
            AUX1 = 1.0_JWRU+0.6522_JWRU*SGDLS+0.4622_JWRU*               &
     &             (SGDLS**2)+0.0864_JWRU*(SGDLS**4)+                    &
     &              0.0675_JWRU*(SGDLS**5)
            AUX2 = 1.0_JWRU/(SGDLS+1.0_JWRU/AUX1)
            WVC = SQRT(AUX2*G*REAL(DEP,JWRU))
            WVK = SIGIN/WVC
            WKDEP = WVK*REAL(DEP,JWRU)
            IF (WKDEP > 13.0_JWRU) THEN
               WN = 0.5_JWRU
            ELSE
               WN = 0.5_JWRU*(1.0_JWRU+2.0_JWRU*WKDEP/                    &
     &              SINH(MIN(300.0_JWRU,2.0_JWRU*WKDEP)))
            END IF
            WVCG = WN*WVC
         ELSE
            WVC  = 0.0_JWRU
            WVK  = 10.0_JWRU
            WVCG = 0.0_JWRU
         END IF

!!!!!!debile: make it all deep
!            WVC = G/SIGIN
!            WN = 0.5_JWRU
!            WVK = SIGIN/WVC
!            WVCG = WN*WVC

      END SUBROUTINE WAVEKCG8
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPLICIT_N_SCHEME_VECTOR(FL1,FLNEW)
      USE yowMpiModule
      USE YOWUNPOOL
      USE yowpd, only: comm
      USE YOWMPP   , ONLY : NINF, NSUP
!*--********************************************************************
! calls       AC2      CCON     CG       COSTH    CURTXY   DIFRM
!             EXCHANGE_P2D      EXCHANGE_P3D_WWM  EXCHANGE_P4D_WWM
!             IEN      IE_CELL2 INE      INVSPHTRANS
!             MPI_ALLREDUCE     MY_WTIME POS_CELL2         SI
!             SINTH    SPSIG    WK
! called by   ** NOTHING **
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  CFLXY    COMM     CX       CY       DIFRU    DT4A
!             DT4AI    DTMAX_EXP         DTMAX_GLOBAL_EXP  FL11
!             FL111    FL112    FL12     FL21     FL211    FL212
!             FL22     FL31     FL311    FL312    FL32     FLALL    I
!             I1       I2       I3       ID       IDIFFR   IE
!             IERR     IOBPD    IOBWB    IP       IPOS     IS       IT
!             ITER_MAX IVECTOR  KELEM    KKSUM    KTMP     LAMBDA
!             LCALC    LDIFR    LSECU    LSPHE    LSTCU    MDC
!             MNE      MNP      MPI_MIN  MSC      N        NI
!             NP_RES   ONE      ONEHALF  ONESIXTH REST     RKIND
!             RTYPE    ST       STAT     THR      TIME1    TIME2
!             TMP      TWO      U        U3       USOC     UTILDE
!             UTILDE3  VERYLARGE         VERYSMALL         WILD
!             WILD2D   WVC      ZERO
! uses PARAMs *** NONE ****
!*++********************************************************************
 
      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: FL1(NINF-1:NSUP,NANG,NFRE)
      REAL(KIND=JWRB), INTENT(INOUT) :: FLNEW(1:MNP,NANG,NFRE)
!
! local INTEGER
!
      INTEGER(KIND=JWIM) :: IP, IE, IT, IS, ID
      INTEGER(KIND=JWIM) :: I1, I2, I3, IERR
      INTEGER(KIND=JWIM) :: NI(3), I, IPOS
!
! local double
!

      REAL(KIND=JWRU) :: UTILDE
      REAL(KIND=JWRU) :: DTMAX_GLOBAL_EXP, DTMAX_EXP
      REAL(KIND=JWRU) :: WILD(MNP), WILD2D(NANG,MNP)
      REAL(KIND=JWRU) :: REST, CFLXY
      REAL(KIND=JWRU) :: LAMBDA(2,NANG,NFRE), DT4AI
      REAL(KIND=JWRU) :: FL11(NANG,NFRE),FL12(NANG,NFRE),FL21(NANG,NFRE)
      REAL(KIND=JWRU) :: FL22(NANG,NFRE),FL31(NANG,NFRE),FL32(NANG,NFRE)
      REAL(KIND=JWRU) :: KTMP(3,NANG,NFRE)
      REAL(KIND=JWRU) :: U3(3)
      REAL(KIND=JWRU) :: KKSUM(NANG,NFRE,MNP)
      REAL(KIND=JWRU) :: ST(NANG,NFRE,MNP)
      REAL(KIND=JWRU) :: N(NANG,NFRE,MNE)
      REAL(KIND=JWRU) :: CX(NANG,NFRE,MNP), CY(NANG,NFRE,MNP)
      REAL(KIND=JWRU) :: U(NANG,NFRE,MNP)
      REAL(KIND=JWRU) :: FLALL(3,NANG,NFRE,MNE)
      REAL(KIND=JWRU) :: KELEM(3,NANG,NFRE,MNE)
      REAL(KIND=JWRU) :: FL111(NANG,NFRE), FL112(NANG,NFRE), FL211(NANG,NFRE)
      REAL(KIND=JWRU) :: FL212(NANG,NFRE), FL311(NANG,NFRE), FL312(NANG,NFRE)
      REAL(KIND=JWRU) :: UTILDE3(MNE)
      REAL(KIND=JWRU) :: USOC, WVC, DIFRU
      REAL(KIND=JWRU) :: ONEOSIX 
!
! local parameter
!
      REAL(KIND=JWRU) :: TMP(NANG,NFRE)
!
!        Calculate phase speeds for the certain spectral comp1.d0nt ...
!
      flall = 0.0_JWRU
      kelem = 0.0_JWRU
      kksum = 0.0_JWRU
      st = 0.0_JWRU
      n = 0.0_JWRU
 
      DO IP = 1, MNP
        DO IS = 1, NFRE
          DO ID = 1, NANG
            CX(ID,IS,IP) = CG(IS,IP)*SINTH(ID)
            CY(ID,IS,IP) = CG(IS,IP)*COSTH(ID)
            IF (LSPHE) THEN
              CX(ID,IS,IP) = CX(ID,IS,IP) * INVTRANS1(IP)
              CY(ID,IS,IP) = CY(ID,IS,IP) * DGRTHM1
            END IF
          END DO
        END DO
      END DO
 
!        Calculate K-Values and contour based quantities ...

      ONEOSIX = 1.0_JWRU/6.0_JWRU
 
!!$OMP DO PRIVATE(IE,I1,I2,I3,LAMBDA,KTMP,TMP,FL11,FL12,FL21,FL22,FL31,FL32,FL111,FL112,FL211,FL212,FL311,FL312)
      DO ie = 1 , mne
        i1 = INE(1,ie)
        i2 = INE(2,ie)
        i3 = INE(3,ie)
        lambda(1,:,:) = ONEOSIX*(cx(:,:,i1)+cx(:,:,i2)+cx(:,:,i3))
        lambda(2,:,:) = ONEOSIX*(cy(:,:,i1)+cy(:,:,i2)+cy(:,:,i3))
        kelem(1,:,:,ie) = lambda(1,:,:)*IEN(1,ie) + lambda(2,:,:)*IEN(2,ie)
        kelem(2,:,:,ie) = lambda(1,:,:)*IEN(3,ie) + lambda(2,:,:)*IEN(4,ie)
        kelem(3,:,:,ie) = lambda(1,:,:)*IEN(5,ie) + lambda(2,:,:)*IEN(6,ie)
        ktmp(:,:,:) = kelem(:,:,:,ie)
        tmp(:,:) = SUM(MIN(0.0_JWRU,ktmp(:,:,:)),DIM=1)
        n(:,:,ie) = -1.0_JWRU/MIN(-thr,tmp(:,:))
        kelem(:,:,:,ie) = MAX(0.0_JWRU,ktmp(:,:,:))
!            WRITE(DBG%FHNDL,'(3I10,3F15.4)') IS, ID, IE, KELEM(:,IE)
        fl11 = cx(:,:,i2)*IEN(1,ie) + cy(:,:,i2)*IEN(2,ie)
        fl12 = cx(:,:,i3)*IEN(1,ie) + cy(:,:,i3)*IEN(2,ie)
        fl21 = cx(:,:,i3)*IEN(3,ie) + cy(:,:,i3)*IEN(4,ie)
        fl22 = cx(:,:,i1)*IEN(3,ie) + cy(:,:,i1)*IEN(4,ie)
        fl31 = cx(:,:,i1)*IEN(5,ie) + cy(:,:,i1)*IEN(6,ie)
        fl32 = cx(:,:,i2)*IEN(5,ie) + cy(:,:,i2)*IEN(6,ie)
        fl111 = 2.0_JWRU*fl11 + fl12
        fl112 = 2.0_JWRU*fl12 + fl11
        fl211 = 2.0_JWRU*fl21 + fl22
        fl212 = 2.0_JWRU*fl22 + fl21
        fl311 = 2.0_JWRU*fl31 + fl32
        fl312 = 2.0_JWRU*fl32 + fl31
        flall(1,:,:,ie) = (fl311+fl212)*ONEOSIX + kelem(1,:,:,ie)
        flall(2,:,:,ie) = (fl111+fl312)*ONEOSIX + kelem(2,:,:,ie)
        flall(3,:,:,ie) = (fl211+fl112)*ONEOSIX + kelem(3,:,:,ie)
      ENDDO
 
      IF ( lcalc ) THEN

#ifdef ebug_adv
         write(dbg%fhndl,*) 'starting the ks buisness'
#endif
         kksum = 0.0_JWRU
         DO ie = 1 , mne
            ni = INE(:,ie)
            kksum(:,:,ni(1)) = kksum(:,:,ni(1)) + kelem(1,:,:,ie)
            kksum(:,:,ni(2)) = kksum(:,:,ni(2)) + kelem(2,:,:,ie)
            kksum(:,:,ni(3)) = kksum(:,:,ni(3)) + kelem(3,:,:,ie)
         ENDDO
#ifdef ebug_adv
         write(dbg%fhndl,*) maxval(kksum), minval(kksum)
         write(dbg%fhndl,*) 'finished the ks buisness'
#endif

         IF ( ivector.EQ.1 ) THEN
#ifdef ebug_adv
         write(dbg%fhndl,*) 'estimate iter_max', ivector
#endif
            DO id = 1 , mdc
               DO is = 1 , msc
                  dtmax_global_exp = 10.E10_JWRU 
                  DO ip = 1 , mnp
                     dtmax_exp = SI(ip)/MAX(thr,kksum(id,is,ip))
                     dtmax_global_exp = MIN(dtmax_global_exp,dtmax_exp)
                  ENDDO
                  dtmax_exp = dtmax_global_exp
                  CALL MPI_ALLREDUCE(dtmax_exp,dtmax_global_exp,1,mpi_real8,mpi_min,comm,ierr)
                  cflxy = dt4a/dtmax_global_exp
                  rest = ABS(MOD(cflxy,1.0_JWRU))
                  IF ( rest.LT.thr ) THEN
                     ITER_EXP(id,is) = ABS(NINT(cflxy))
                  ELSEIF ( rest.GT.thr .AND. rest.LT.0.5_JWRU ) THEN
                     ITER_EXP(id,is) = ABS(NINT(cflxy)) + 1
                  ELSE
                     ITER_EXP(id,is) = ABS(NINT(cflxy))
                  ENDIF
               ENDDO
            ENDDO
#ifdef ebug_adv
         write(dbg%fhndl,*) 'iter_max = ', maxval(ITER_EXP)
#endif
         ELSEIF ( ivector.EQ.2 ) THEN
#ifdef ebug_adv
         write(dbg%fhndl,*) 'estimate iter_max', ivector
#endif
            dtmax_global_exp = 10.E10
            DO ip = 1 , mnp 
               dtmax_exp = SI(ip)/MAX(thr,MAXVAL(kksum(:,:,ip)))
               dtmax_global_exp = MIN(dtmax_global_exp,dtmax_exp)
            ENDDO
            dtmax_exp = dtmax_global_exp
            CALL MPI_ALLREDUCE(dtmax_exp,dtmax_global_exp,1,mpi_real8,mpi_min,comm,ierr)
            cflxy = dt4a/dtmax_global_exp
            rest = ABS(MOD(cflxy,1.0_JWRU))
            IF ( rest.LT.thr ) THEN
               iter_max = ABS(NINT(cflxy))
            ELSEIF ( rest.GT.thr .AND. rest.LT.0.5_JWRU ) THEN
               iter_max = ABS(NINT(cflxy)) + 1
            ELSE
               iter_max = ABS(NINT(cflxy))
            ENDIF
#ifdef ebug_adv
         write(dbg%fhndl,*) 'iter_max = ', iter_max
#endif
         ELSEIF ( ivector.EQ.3 ) THEN
#ifdef ebug_adv
         write(dbg%fhndl,*) 'estimate iter_max', ivector
#endif
            DO is = 1 , msc
               dtmax_global_exp = 10.E10_JWRU
               DO ip = 1 , mnp 
                  dtmax_exp = SI(ip)/MAX(thr,MAXVAL(kksum(:,is,ip)))
                  dtmax_global_exp = MIN(dtmax_global_exp,dtmax_exp)
               ENDDO
               dtmax_exp = dtmax_global_exp
               CALL MPI_ALLREDUCE(dtmax_exp,dtmax_global_exp,1,mpi_real8,mpi_min,comm,ierr)
               cflxy = dt4a/dtmax_global_exp
               rest = ABS(MOD(cflxy,1.0_JWRU))
               IF ( rest.LT.thr ) THEN
                  ITER_EXPD(is) = ABS(NINT(cflxy))
               ELSEIF ( rest.GT.thr .AND. rest.LT.0.5_JWRU ) THEN
                  ITER_EXPD(is) = ABS(NINT(cflxy)) + 1
               ELSE
                  ITER_EXPD(is) = ABS(NINT(cflxy))
               ENDIF
            ENDDO
#ifdef ebug_adv
         write(dbg%fhndl,*) 'iter_max = ', maxval(iter_expd)
#endif
         ELSEIF ( ivector.EQ.4 ) THEN
#ifdef ebug_adv
         write(dbg%fhndl,*) 'estimate iter_max', ivector
#endif
            DO is = 1 , msc
               dtmax_global_exp = 10.E10
               DO ip = 1 , mnp 
                  dtmax_exp = SI(ip)/MAX(thr,MAXVAL(kksum(:,is,ip)))
                  dtmax_global_exp = MIN(dtmax_global_exp,dtmax_exp)
               ENDDO
               dtmax_exp = dtmax_global_exp
               CALL MPI_ALLREDUCE(dtmax_exp,dtmax_global_exp,1,mpi_real8,mpi_min,comm,ierr)
               cflxy = dt4a/dtmax_global_exp
               rest = ABS(MOD(cflxy,1.0_JWRU))
               IF ( rest.LT.thr ) THEN
                  ITER_EXPD(is) = ABS(NINT(cflxy))
               ELSEIF ( rest.GT.thr .AND. rest.LT.0.5_JWRU ) THEN
                  ITER_EXPD(is) = ABS(NINT(cflxy)) + 1
               ELSE
                  ITER_EXPD(is) = ABS(NINT(cflxy))
               ENDIF
            ENDDO
#ifdef ebug_adv
         write(dbg%fhndl,*) 'iter_max = ', maxval(iter_expd) 
#endif
         ELSEIF ( ivector.EQ.5 ) THEN
#ifdef ebug_adv
         write(dbg%fhndl,*) 'estimate iter_exp', lvector 
#endif
            dtmax_global_exp = 10.E10_JWRU
            DO ip = 1 , mnp 
               dtmax_exp = SI(ip)/MAX(thr,MAXVAL(kksum(:,:,ip)))
               dtmax_global_exp = MIN(dtmax_global_exp,dtmax_exp)
            ENDDO
            dtmax_exp = dtmax_global_exp
            CALL MPI_ALLREDUCE(dtmax_exp,dtmax_global_exp,1,MPI_REAL8,mpi_min,comm,ierr)
            cflxy = dt4a/dtmax_global_exp
            rest = ABS(MOD(cflxy,1.0_JWRU))
            IF ( rest.LT.thr ) THEN
               iter_max = ABS(NINT(cflxy))
            ELSEIF ( rest.GT.thr .AND. rest.LT.0.5_JWRU ) THEN
               iter_max = ABS(NINT(cflxy)) + 1
            ELSE
               iter_max = ABS(NINT(cflxy))
            ENDIF
#ifdef ebug_adv
         write(dbg%fhndl,*) 'iter_max = ', iter_max
#endif
         ENDIF    !IVECTOR
         CALL FLUSH(dbg%fhndl)
      ENDIF     !LCALC


      DO ip = 1 , mnp
        u(:,:,ip) = FL1(ip,:,:)
      ENDDO

#ifdef ebug_adv
      write(dbg%fhndl,*) 'doing time integration'
      call cpu_time(time1)
#endif

      IF ( ivector.EQ.1 ) THEN
         DO id = 1 , mdc
            DO is = 1 , msc
               dt4ai = dt4a/ITER_EXP(id,is)
               DO it = 1 , ITER_EXP(id,is)
                  st(id,is,:) = 0.0_JWRU
                  DO ie = 1 , mne
                     ni = INE(:,ie)
                     u3(:) = u(id,is,ni)
                     utilde = n(id,is,ie)                                  &
     &                        *(flall(1,id,is,ie)*u3(1)+flall(2,id,is,     &
     &                        ie)*u3(2)+flall(3,id,is,ie)*u3(3))
                     st(id,is,ni(1)) = st(id,is,ni(1))                     &
     &                                 + kelem(1,id,is,ie)                 &
     &                                 *(u3(1)-utilde)
                     st(id,is,ni(2)) = st(id,is,ni(2))                     &
     &                                 + kelem(2,id,is,ie)                 &
     &                                 *(u3(2)-utilde)
                     st(id,is,ni(3)) = st(id,is,ni(3))                     &
     &                                 + kelem(3,id,is,ie)                 &
     &                                 *(u3(3)-utilde)
                  ENDDO !IE
                  u(id,is,:) = MAX(0.0_JWRU,u(id,is,:)-dt4ai/SI*st(id,is,:)    &
     &                         *iobwb)*iobpd(id,:)
                  wild = u(id,is,:)
                  CALL exchange(wild)
                  u(id,is,:) = wild
               ENDDO ! IT----> End Iteration
            ENDDO !IS
         ENDDO  !ID
      ELSEIF ( ivector.EQ.2 ) THEN
         dt4ai = dt4a/iter_max
         DO it = 1 , iter_max
            DO id = 1 , mdc
               DO is = 1 , msc
                  st(id,is,:) = 0.0_JWRU
                  DO ie = 1 , mne
                     ni = INE(:,ie)
                     u3(:) = u(id,is,ni)
                     utilde = n(id,is,ie)                                  &
     &                        *(flall(1,id,is,ie)*u3(1)+flall(2,id,is,     &
     &                        ie)*u3(2)+flall(3,id,is,ie)*u3(3))
                     st(id,is,ni(1)) = st(id,is,ni(1))                     &
     &                                 + kelem(1,id,is,ie)                 &
     &                                 *(u3(1)-utilde)
                     st(id,is,ni(2)) = st(id,is,ni(2))                     &
     &                                 + kelem(2,id,is,ie)                 &
     &                                 *(u3(2)-utilde)
                     st(id,is,ni(3)) = st(id,is,ni(3))                     &
     &                                 + kelem(3,id,is,ie)                 &
     &                                 *(u3(3)-utilde)
                  ENDDO!IE
                  u(id,is,:) = MAX(0.0_JWRU,u(id,is,:)-dt4ai/SI*st(id,is,:)    &
     &                         *iobwb)*iobpd(id,:)
               ENDDO !IS
            ENDDO !ID
            CALL exchange(u)
         ENDDO  !IT
      ELSEIF ( ivector.EQ.3 ) THEN
         DO is = 1 , msc
            iter_max = ITER_EXPD(is)
            dt4ai = dt4a/iter_max
            DO it = 1 , iter_max
               DO id = 1 , mdc
                  st(id,is,:) = 0.0_JWRU
                  DO ie = 1 , mne
                     ni = INE(:,ie)
                     u3(:) = u(id,is,ni)
                     utilde = n(id,is,ie)                                  &
     &                        *(flall(1,id,is,ie)*u3(1)+flall(2,id,is,     &
     &                        ie)*u3(2)+flall(3,id,is,ie)*u3(3))
                     st(id,is,ni(1)) = st(id,is,ni(1))                     &
     &                                 + kelem(1,id,is,ie)                 &
     &                                 *(u3(1)-utilde)
                     st(id,is,ni(2)) = st(id,is,ni(2))                     &
     &                                 + kelem(2,id,is,ie)                 &
     &                                 *(u3(2)-utilde)
                     st(id,is,ni(3)) = st(id,is,ni(3))                     &
     &                                 + kelem(3,id,is,ie)                 &
     &                                 *(u3(3)-utilde)
                  ENDDO
                  u(id,is,:) = MAX(0.0_JWRU,u(id,is,:)-dt4ai/SI*st(id,is,:)    &
     &                         *iobwb)*iobpd(id,:)
               ENDDO
               wild2d = u(:,is,:)
               CALL exchange(wild2d)
               u(:,is,:) = wild2d
            ENDDO
         ENDDO
      ELSEIF ( ivector.EQ.4 ) THEN
         DO is = 1 , msc
            iter_max = ITER_EXPD(is)
            dt4ai = dt4a/iter_max
            DO it = 1 , iter_max
               DO id = 1 , mdc
                  DO ie = 1 , mne
                     ni = INE(:,ie)
                     u3(:) = u(id,is,ni)
                     utilde3(ie) = n(id,is,ie)                             &
     &                             *(flall(1,id,is,ie)*u3(1)+flall(2,id,   &
     &                             is,ie)*u3(2)+flall(3,id,is,ie)*u3(3))
                  ENDDO
                  st(id,is,:) = 0.0_JWRU
                  DO ip = 1 , mnp
                     DO i = 1 , CCON(ip)
                        ie = IE_CELL2(ip,i)
                        ipos = POS_CELL2(ip,i)
                        st(id,is,ip) = st(id,is,ip)                        &
     &                                 + kelem(ipos,id,is,ie)              &
     &                                 *(u(id,is,ip)-utilde3(ie))
                     ENDDO
                     u(id,is,ip) = MAX(0.0_JWRU,u(id,is,ip)-dt4ai/SI(ip)*st(   &
     &                             id,is,ip)*iobwb(ip))*iobpd(id,ip)
                  ENDDO
               ENDDO
               wild2d = u(:,is,:)
               CALL exchange(wild2d)
               u(:,is,:) = wild2d
            ENDDO !IT
         ENDDO  !IS
      ELSEIF ( ivector.EQ.5 ) THEN
         dt4ai = dt4a/iter_max
         DO it = 1 , iter_max
            DO is = 1 , msc
               DO id = 1 , mdc
                  DO ie = 1 , mne
                     ni = INE(:,ie)
                     u3(:) = u(id,is,ni)
                     utilde3(ie) = n(id,is,ie)                             &
     &                             *(flall(1,id,is,ie)*u3(1)+flall(2,id,   &
     &                             is,ie)*u3(2)+flall(3,id,is,ie)*u3(3))
                  ENDDO
                  st(id,is,:) = 0.0_JWRU
                  DO ip = 1 , mnp
                     DO i = 1 , CCON(ip)
                        ie = IE_CELL2(ip,i)
                        ipos = POS_CELL2(ip,i)
                        st(id,is,ip) = st(id,is,ip)                        &
     &                                 + kelem(ipos,id,is,ie)              &
     &                                 *(u(id,is,ip)-utilde3(ie))
                     ENDDO
                     u(id,is,ip) = MAX(0.0_JWRU,u(id,is,ip)-dt4ai/SI(ip)*st(   &
     &                             id,is,ip)*iobwb(ip))*iobpd(id,ip)
                  ENDDO
               ENDDO
            ENDDO
            CALL exchange(u)
         ENDDO
      ENDIF

        DO ip = 1,mnp
#ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'ip=', ip, ' mnp=', mnp
          FLUSH(740+MyRankGlobal)
#endif
          CALL REFRACTION_FREQSHIFT_EXPLICIT_SINGLE(U(:,:,ip), IP)
        END DO

      
      DO ip = 1 , mnp
         FLNEW(ip,:,:) = REAL(u(:,:,ip),JWRB)
      ENDDO

      END SUBROUTINE EXPLICIT_N_SCHEME_VECTOR
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EIMPS_ASPAR_BLOCK(ASPAR)
      USE YOWUNPOOL, ONLY : LCUR, LSPHE, CG, ONETHIRD, ONESIXTH, IEN, ONE, TWO, THR, TRIA, MDC, SI
      IMPLICIT NONE
      INTEGER(KIND=JWIM) :: POS_TRICK(3,2)
      INTEGER(KIND=JWIM) :: I1, I2, I3
      INTEGER(KIND=JWIM) :: IP, ID, IS, IE
      INTEGER(KIND=JWIM) :: I, idx

      REAL(KIND=JWRU), INTENT(INOUT) :: ASPAR(NANG, NFRE, NNZ)
      REAL(KIND=JWRU) :: FL11(NANG,NFRE), FL12(NANG,NFRE), FL21(NANG,NFRE), FL22(NANG,NFRE), FL31(NANG,NFRE), FL32(NANG,NFRE)
      REAL(KIND=JWRU) :: CRFS(NANG,NFRE,3), K1(NANG,NFRE), KM(NANG,NFRE,3), K(NANG,NFRE,3), TRIA03
      REAL(KIND=JWRU) :: CXY(2,NANG,NFRE,3)
      REAL(KIND=JWRU) :: DELTAL(NANG,NFRE,3)
      REAL(KIND=JWRU) :: KP(NANG,NFRE,3), NM(NANG,NFRE)
      REAL(KIND=JWRU) :: DTK(NANG,NFRE), TMP3(NANG,NFRE)
      REAL(KIND=JWRU) :: LAMBDA(2,NANG,NFRE)
      REAL(KIND=JWRU) sumCXY, sumTMP3, sumCG, sumASPARdiag
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2
!
!     Calculate countour integral quantities ...
!
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'sum(IOBWB)=', sum(IOBWB)
      WRITE(740+MyRankGlobal,*) 'sum(IOBPD)=', sum(IOBPD)
      WRITE(740+MyRankGlobal,*) 'LSPHE=', LSPHE
      sumTMP3 = 0
      sumCXY  = 0
      sumCG  = 0
#endif
      ASPAR = 0.0_JWRU
      DO IE = 1, MNE
        DO I=1,3
          IP = INE(I,IE)
          DO IS=1,NFRE
            DO ID=1,NANG
              IF (LCUR) THEN
                CXY(1,ID,IS,I) = CG(IS,IP) * SINTH(ID) + CURTXY(1,IP)
                CXY(2,ID,IS,I) = CG(IS,IP) * COSTH(ID) + CURTXY(2,IP)
              ELSE
                CXY(1,ID,IS,I) = CG(IS,IP) * SINTH(ID)
                CXY(2,ID,IS,I) = CG(IS,IP) * COSTH(ID)
              END IF
              IF (LSPHE) THEN
                CXY(1,ID,IS,I) = CXY(1,ID,IS,I) * INVTRANS1(IP)
                CXY(2,ID,IS,I) = CXY(2,ID,IS,I) * DGRTHM1
              END IF
            END DO
          END DO
        END DO
#ifdef DEBUG
        sumCXY  = sumCXY + sum(abs(CXY))
        sumCG   = sumCG  + sum(abs(CG ))
#endif

        LAMBDA(:,:,:) = ONESIXTH * (CXY(:,:,:,1) + CXY(:,:,:,2) + CXY(:,:,:,3))
        K(:,:,1)  = LAMBDA(1,:,:) * IEN(1,IE) + LAMBDA(2,:,:) * IEN(2,IE)
        K(:,:,2)  = LAMBDA(1,:,:) * IEN(3,IE) + LAMBDA(2,:,:) * IEN(4,IE)
        K(:,:,3)  = LAMBDA(1,:,:) * IEN(5,IE) + LAMBDA(2,:,:) * IEN(6,IE)
        FL11(:,:) = CXY(1,:,:,2)*IEN(1,IE)+CXY(2,:,:,2)*IEN(2,IE)
        FL12(:,:) = CXY(1,:,:,3)*IEN(1,IE)+CXY(2,:,:,3)*IEN(2,IE)
        FL21(:,:) = CXY(1,:,:,3)*IEN(3,IE)+CXY(2,:,:,3)*IEN(4,IE)
        FL22(:,:) = CXY(1,:,:,1)*IEN(3,IE)+CXY(2,:,:,1)*IEN(4,IE)
        FL31(:,:) = CXY(1,:,:,1)*IEN(5,IE)+CXY(2,:,:,1)*IEN(6,IE)
        FL32(:,:) = CXY(1,:,:,2)*IEN(5,IE)+CXY(2,:,:,2)*IEN(6,IE)
        CRFS(:,:,1) = - ONESIXTH *  (TWO *FL31(:,:) + FL32(:,:) + FL21(:,:) + TWO * FL22(:,:) )
        CRFS(:,:,2) = - ONESIXTH *  (TWO *FL32(:,:) + TWO * FL11(:,:) + FL12(:,:) + FL31(:,:) )
        CRFS(:,:,3) = - ONESIXTH *  (TWO *FL12(:,:) + TWO * FL21(:,:) + FL22(:,:) + FL11(:,:) )
        KM = MIN(0.0_JWRU,K)
        KP(:,:,:) = MAX(0.0_JWRU,K)
        DELTAL(:,:,:) = CRFS(:,:,:)- KP(:,:,:)
        NM(:,:)=ONE/MIN(-THR,KM(:,:,1) + KM(:,:,2) + KM(:,:,3))
        TRIA03 = ONETHIRD * TRIA(IE)
        DO I=1,3
          IP=INE(I,IE)
          I1=JA_IE(I,1,IE)
          I2=JA_IE(I,2,IE)
          I3=JA_IE(I,3,IE)
          K1(:,:) =  KP(:,:,I)
          DO ID=1,NANG
            DTK(ID,:) =  K1(ID,:) * DT4A * IOBPD(ID,IP) * IOBWB(IP)
          END DO
          TMP3(:,:)  =  DTK(:,:) * NM(:,:)
#ifdef DEBUG
          sumTMP3 = sumTMP3 + sum(abs(TMP3))
#endif
          ASPAR(:,:,I1) =  TRIA03 + DTK(:,:) - TMP3(:,:) * DELTAL(:,:,I             ) + ASPAR(:,:,I1)
          ASPAR(:,:,I2) =                    - TMP3(:,:) * DELTAL(:,:,POS_TRICK(I,1)) + ASPAR(:,:,I2)
          ASPAR(:,:,I3) =                    - TMP3(:,:) * DELTAL(:,:,POS_TRICK(I,2)) + ASPAR(:,:,I3)
        END DO
      END DO
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'sumTMP3=', sumTMP3
      WRITE(740+MyRankGlobal,*) 'sumCXY=', sumCXY
      WRITE(740+MyRankGlobal,*) 'sumCG =', sumCG
#endif
      IF (LBCWA) THEN
        DO idx = 1, IWBMNP
          IP = IWBNDLC(idx)
          ASPAR(:,:,I_DIAG(IP)) = SI(IP) ! Set boundary on the diagonal
        END DO
      END IF
#ifdef DEBUG
      sumASPARdiag=0
      DO IP=1,MNP
        sumASPARdiag = sumASPARdiag + sum(abs(ASPAR(:,:,I_DIAG(IP))))
      END DO
      WRITE(740+MyRankGlobal,*) 'sumASPARdiag=', sumASPARdiag
#endif
      END SUBROUTINE EIMPS_ASPAR_BLOCK
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ADD_FREQ_DIR_TO_ASPAR_COMP_CADS(ASPAR_JAC)
      USE YOWUNPOOL, ONLY : LSPHE
      USE YOWFRED,  ONLY : DELTH, FRATIO, FR
      USE YOWGRID , ONLY : NPROMA_WAM, NCHNK, ICHNKFROMIJ, IPRMFROMIJ
      USE YOWSTAT,  ONLY : IDELPRO, IREFRA
      USE YOWGRID,  ONLY : DELPHI, IJS, IJL
      USE YOWPCONS, ONLY : PI, ZPI
      USE YOWREFD,  ONLY : THDD, THDC, SDOT
      USE YOWSHAl,  ONLY : WVPRPT 
      USE yowpd, only : np_global
      USE YOWUNPOOL, ONLY : SI
      IMPLICIT NONE
      REAL(KIND=JWRU), intent(inout) :: ASPAR_JAC(NANG,NFRE,NNZ)
      REAL(KIND=JWRU) SS, SC, CC, SD, CD
      REAL(KIND=JWRU) DELPH0, DELTH0, DELFR0, DELPRO
      REAL(KIND=JWRU) DRCP(NANG), DRCM(NANG), DRDP(NANG), DRDM(NANG)
      INTEGER(KIND=JWIM) :: IP, IJ, K, M, KM1, KP1, MP1, MM1, IPRM, IK
      REAL(KIND=JWRU) DFP, DFM
      REAL(KIND=JWRU) DTHP, DTHM
      REAL(KIND=JWRU) CASS(0:NFRE+1)
      REAL(KIND=JWRU) CAS(NFRE,NANG)
      REAL(KIND=JWRU) CAD(NANG,NFRE)
      REAL(KIND=JWRU) CP_THE(NANG,NFRE), CM_THE(NANG,NFRE)
      REAL(KIND=JWRU) CP_SIG(0:NFRE+1), CM_SIG(0:NFRE+1)
      REAL(KIND=JWRU) eFact
      REAL(KIND=JWRU) B_SIG(NFRE)
      INTEGER(KIND=JWIM) :: ID, IS, TheVal
      DELPRO = REAL(IDELPRO,JWRU)
      DELPH0 = 0.25_JWRU*DELPRO/DELPHI
      DELTH0 = 0.25_JWRU*DELPRO/DELTH
      DELFR0 = 0.25_JWRU*DELPRO/((FRATIO-1)*ZPI)
      IF (REFRA_METHOD .eq. 1) THEN
        IF (IREFRA.NE.0) THEN
          DO IP=1,NP_RES
            IJ = IP
            IK = ICHNKFROMIJ(IJ)
            IPRM = IPRMFROMIJ(IJ)

            DO K=1,NANG
              KP1 = K+1
              IF (KP1.GT.NANG) KP1 = 1
              KM1 = K-1
              IF (KM1.LT.1) KM1 = NANG
              DRDP(K) = (THDD(IJ, K) + THDD(IJ, KP1))*DELTH0
              DRDM(K) = (THDD(IJ, K) + THDD(IJ, KM1))*DELTH0
            END DO

            IF ((IREFRA.EQ.2).OR.(IREFRA.EQ.3)) THEN
              DO K=1,NANG
                KP1 = K+1
                IF (KP1.GT.NANG) KP1 = 1
                KM1 = K-1
                IF (KM1.LT.1) KM1 = NANG
                DRCP(K) = (THDC(IJ,K) + THDC(IJ,KP1))*DELTH0
                DRCM(K) = (THDC(IJ,K) + THDC(IJ,KM1))*DELTH0
              END DO
            END IF
            DO M=1,NFRE
              DO K=1,NANG
                MP1 = MIN(NFRE,M+1)
                MM1 = MAX(1,M-1)
                DFP = DELFR0/FR(M)
                DFM = DELFR0/FR(MM1)
                !

                DTHP = REAL(WVPRPT%OMOSNH2KD(IPRM,M,IK),JWRU)*DRDP(K) + DRCP(K)
                DTHM = REAL(WVPRPT%OMOSNH2KD(IPRM,M,IK),JWRU)*DRDM(K) + DRCM(K)


                ASPAR_JAC(K,M,I_DIAG(IP)) = ASPAR_JAC(K,M,I_DIAG(IP)) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
                DTP_I(K,M,IJ) = -DTHP+ABS(DTHP)
                DTM_I(K,M,IJ) =  DTHM+ABS(DTHM)
                !
                DTHP = (SDOT(K,M,IJ) + SDOT(K,MP1,IJ))*DFP
                DTHM = (SDOT(K,M,IJ) + SDOT(K,MM1,IJ))*DFM
                ASPAR_JAC(K,M,I_DIAG(IP)) = ASPAR_JAC(K,M,I_DIAG(IP)) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
                DOP_I(K,M,IJ) = (-DTHP+ABS(DTHP))/FRATIO
                DOM_I(K,M,IJ) = ( DTHM+ABS(DTHM))*FRATIO
              END DO
            END DO
          END DO
        END IF
      END IF
      IF (REFRA_METHOD .eq. 2) THEN
        IF ((IREFRA .eq. 1).or.(IREFRA.eq.3).or.(LSPHE.eqv..TRUE.)) THEN
          DO IP=1,NP_RES
            TheVal=1
            IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3)) TheVal=0
            IF (DEP(IP) .LT. DMIN) TheVal=0
            IF (IOBP(IP) .EQ. 2) TheVal=0
            IF (TheVal .eq. 1) THEN
              CALL PROPTHETA(IP,CAD)
            ELSE
              CAD=ZERO
            END IF
            CP_THE = MAX(ZERO,CAD)
            CM_THE = MIN(ZERO,CAD)
            eFact=(DT4D/DDIR)*SI(IP)
            CAD_THE(:,:,IP)=CAD
            ASPAR_JAC(:,:,I_DIAG(IP)) = ASPAR_JAC(:,:,I_DIAG(IP)) + eFact * (CP_THE(:,:) - CM_THE(:,:))
          END DO
        END IF
        IF ((IREFRA.eq.2).or.(IREFRA.eq.3)) THEN
          DO IP=1,NP_RES
            TheVal=1
            IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3)) TheVal=0
            IF (DEP(IP) .LT. DMIN) TheVal=0
            IF (IOBP(IP) .EQ. 2) TheVal=0
            IF (TheVal .eq. 1) THEN
              CALL PROPSIGMA(IP,CAS)
            ELSE
              CAS=ZERO
            END IF
            CAS_SIG(:,:,IP) = CAS
            eFact=DT4F*SI(IP)
            DO ID = 1, NANG
              CASS(1:NFRE) = CAS(:,ID)
              CASS(0)     = 0._JWRU
              CASS(NFRE+1) = CASS(NFRE)
              CP_SIG = MAX(ZERO,CASS)
              CM_SIG = MIN(ZERO,CASS)
              DO IS=1,NFRE
                B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
              END DO
              B_SIG(NFRE) = B_SIG(NFRE) + eFact*CM_SIG(NFRE+1)/DS_INCR(NFRE) * PTAIL5
              ASPAR_JAC(:,ID,I_DIAG(IP)) = ASPAR_JAC(:,ID,I_DIAG(IP)) + B_SIG
            END DO
          END DO
        END IF
      END IF
      END SUBROUTINE ADD_FREQ_DIR_TO_ASPAR_COMP_CADS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_IMATRA_IMATDA(IP, AC, IMATRA, IMATDA)
      USE YOWUNPOOL
      USE yowpd, only : np_global, NP_RES => np, MNP=>npa
      USE YOWMPP   , ONLY : NINF, NSUP
      IMPLICIT NONE
      INTEGER(KIND=JWIM), INTENT(IN) :: IP
      REAL(KIND=JWRU), INTENT(IN)  :: AC(NANG,NFRE,MNP)
      REAL(KIND=JWRU), intent(inout) :: IMATRA(NANG,NFRE)
      REAL(KIND=JWRU), intent(inout) :: IMATDA(NANG,NFRE)
      REAL(KIND=JWRU) eVal
      IF (LNONL) THEN
         !*    Advanced computation, likely call to IMPLSCH
      ELSE
         !*    Put here the fields for IMATRA
      END IF
      eVal=SI(IP)
      IMATRA = IMATRA * eVal
      IMATDA = IMATDA * eVal
      END SUBROUTINE GET_IMATRA_IMATDA
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ACTION_LIMITER_local(IP,eSum,acloc)
      IMPLICIT NONE
      INTEGER(KIND=JWIM), INTENT(IN) :: IP
      REAL(KIND=JWRU), intent(INOUT) :: eSum(NANG,NFRE)
      REAL(KIND=JWRU), intent(IN) :: acloc(NANG,NFRE)


      Print *, 'This needs to be written'
      stop
      END SUBROUTINE ACTION_LIMITER_local
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_BLOCAL(IP, AC, BLOC)
      USE YOWUNPOOL
      USE UNSTRUCT_BOUND, ONLY : IWBNDLC_REV
      USE yowpd, only : np_global
      IMPLICIT NONE
      INTEGER(KIND=JWIM), INTENT(IN) :: IP
      REAL(KIND=JWRU), INTENT(IN)  :: AC(NANG,NFRE,MNP)
      REAL(KIND=JWRU), intent(out) :: BLOC(NANG,NFRE)
      INTEGER(KIND=JWIM) :: idx, ID
      idx=IWBNDLC_REV(IP)
      IF (LBCWA .and. (idx.gt.0)) THEN
        BLOC = WBAC(:,:,idx) * SI(IP)
      ELSE
        BLOC = AC(:,:,IP) * SI(IP)
        DO ID=1,NANG
          BLOC(ID,:) = BLOC(ID,:) * IOBPD(ID,IP) * IOBWB(IP)
        END DO
      END IF
      END SUBROUTINE GET_BLOCAL
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE IMPLICIT_N_SCHEME_BLOCK(FL1, FLNEW)
      USE yowMpiModule
      USE yowpd, only : comm
      USE yowpd, only : np_global
      USE YOWMPP, ONLY : NINF, NSUP
      USE YOWSTAT, ONLY : IREFRA
      USE yowunpool
      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: itmp(1), isend(1)
      INTEGER(KIND=JWIM) :: IS, ID, ID1, ID2, IP, J, idx, nbITer, TheVal, is_converged
      INTEGER(KIND=JWIM) :: I, K, IP_ADJ, IADJ, JDX
      INTEGER(KIND=JWIM) :: KP1, KM1, MP1, MM1, M, IJ, IJS
      INTEGER(KIND=JWIM) :: ierr

      REAL(KIND=JWRB), INTENT(IN)  :: FL1(NINF-1:NSUP,NANG,NFRE)
      REAL(KIND=JWRB), INTENT(INOUT) :: FLNEW(1:MNP,NANG,NFRE)

      REAL(KIND=JWRU) MaxNorm, SumNorm, p_is_converged
      REAL(KIND=JWRU) eSum(NANG,NFRE)
      REAL(KIND=JWRU) IMATRA(NANG,NFRE), IMATDA(NANG,NFRE)
      REAL(KIND=JWRU) Norm_L2(NANG,NFRE), Norm_LINF(NANG,NFRE)
      REAL(KIND=JWRU) ACLOC(NANG,NFRE)
      REAL(KIND=JWRU) CAD(NANG,NFRE), CAS(NANG,NFRE)
      REAL(KIND=JWRU) BLOC(NANG,NFRE)
      REAL(KIND=JWRU) ASPAR_DIAG(NANG,NFRE)
      REAL(KIND=JWRU) Norm_L2_gl(NANG,NFRE), Norm_LINF_gl(NANG,NFRE)
      REAL(KIND=JWRU) B_SIG(NFRE), eFact, lambda
      REAL(KIND=JWRU) Sum_new, Sum_prev, eVal, DiffNew, DiffOld
      REAL(KIND=JWRU) AC1(NANG,NFRE,MNP)
      REAL(KIND=JWRU) AC2(NANG,NFRE,MNP)
      REAL(KIND=JWRU) FLsing(NANG,NFRE)
      REAL(KIND=JWRU) sumASPARoff
      REAL(KIND=JWRU) CP_THE(NANG,NFRE), CM_THE(NANG,NFRE)
      REAL(KIND=JWRU) CP_SIG(NFRE,NANG), CM_SIG(NFRE,NANG)
#ifdef DEBUG
      REAL(KIND=JWRU) sumESUMbefore, sumESUMafter, sumESUMpre
#endif
      IJS = 1
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'Begin IMPLICIT SOLVER'
      WRITE(740+MyRankGlobal,*) 'JWRU=', JWRU
      WRITE(740+MyRankGlobal,*) 'JWRB=', JWRB
      WRITE(740+MyRankGlobal,*) 'IREFRA=', IREFRA
      FLUSH(740+MyRankGlobal)
#endif

      DO IP=1,MNP
        AC1(:,:,IP)=REAL(FL1(IP,:,:),JWRU) ! no this is not an error
        AC2(:,:,IP)=REAL(FL1(IP,:,:),JWRU) ! no this is not an error
      END DO
      CALL EIMPS_ASPAR_BLOCK(ASPAR_JAC)
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'After EIMPS_ASPAR_BLOCK'
      WRITE(740+MyRankGlobal,*) 'sum(abs(ASPAR_JAC))=', sum(abs(ASPAR_JAC))
      WRITE(740+MyRankGlobal,*) 'sum(    ASPAR_JAC )=', sum(ASPAR_JAC)
      FLUSH(740+MyRankGlobal)
#endif
      !
      CALL ADD_FREQ_DIR_TO_ASPAR_COMP_CADS(ASPAR_JAC)
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'After ADD_FREQ_DIR_TO_ASPAR_COMP_CADS'
      WRITE(740+MyRankGlobal,*) 'sum(abs(ASPAR_JAC))=', sum(abs(ASPAR_JAC))
      WRITE(740+MyRankGlobal,*) 'sum(    ASPAR_JAC )=', sum(ASPAR_JAC)
      WRITE(740+MyRankGlobal,*) 'SOURCE_IMPL=', SOURCE_IMPL
      WRITE(740+MyRankGlobal,*) 'LNONL=', LNONL
      FLUSH(740+MyRankGlobal)
#endif
      IF ((.NOT. LNONL) .AND. SOURCE_IMPL) THEN
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'In the loop'
        FLUSH(740+MyRankGlobal)
#endif
        DO IP=1,NP_RES
#ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'IP=', IP, ' NP_RES=', NP_RES
          FLUSH(740+MyRankGlobal)
#endif
          CALL GET_BLOCAL(IP, AC1, BLOC)
#ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'step 1'
          FLUSH(740+MyRankGlobal)
#endif
          CALL GET_IMATRA_IMATDA(IP, AC1, IMATRA, IMATDA)
#ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'step 2'
          FLUSH(740+MyRankGlobal)
#endif
          ASPAR_JAC(:,:,I_DIAG(IP)) = ASPAR_JAC(:,:,I_DIAG(IP)) + IMATDA
#ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'step 3'
          FLUSH(740+MyRankGlobal)
#endif
          B_JAC(:,:,IP)             = BLOC + IMATRA
#ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'step 4'
          FLUSH(740+MyRankGlobal)
#endif
        END DO
      END IF
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'After B_JAC and ASPAR_JAC assignations'
      FLUSH(740+MyRankGlobal)
#endif
      !
      ! Now the Gauss Seidel iterations
      !
      !SOLVERTHR=10E-8*AVETL!*TLMIN**2
      !
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'BLOCK_GAUSS_SEIDEL=', BLOCK_GAUSS_SEIDEL
      WRITE(740+MyRankGlobal,*) 'SOURCE_IMPL=', SOURCE_IMPL
      WRITE(740+MyRankGlobal,*) 'LLIMT=', LLIMT
      WRITE(740+MyRankGlobal,*) 'LCHKCONV=', LCHKCONV
      WRITE(740+MyRankGlobal,*) 'NANG=', NANG
      WRITE(740+MyRankGlobal,*) 'NFRE=', NFRE
      WRITE(740+MyRankGlobal,*) 'L_SOLVER_NORM=', L_SOLVER_NORM
#endif
      nbIter=0
      DO
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'nbIter=', nbIter
        WRITE(740+MyRankGlobal,*) 'Before iteration, sum(WBAC)=', sum(WBAC)
        WRITE(740+MyRankGlobal,*) 'Before iteration, sum(AC2)=', sum(AC2)
        WRITE(740+MyRankGlobal,*) 'Before iteration, sum(FL1)=', sum(FL1(1:MNP,:,:))
        WRITE(740+MyRankGlobal,*) 'Before iteration, minval(SI)=', minval(SI)
        WRITE(740+MyRankGlobal,*) 'Before iteration, maxval(SI)=', maxval(SI)
        sumESUMbefore=0
        sumESUMafter=0
        FLUSH(740+MyRankGlobal)
#endif 
        is_converged = 0
        DO IP=1,NP_RES
#ifdef DEBUG
!          WRITE(740+MyRankGlobal,*) 'Main implicit loop IP=', IP
!          FLUSH(740+MyRankGlobal)
#endif
          IJ = IJS - 1 + IP
          ACLOC = AC2(:,:,IP)
          ASPAR_DIAG=ASPAR_JAC(:,:,I_DIAG(IP))
          IF (SOURCE_IMPL) THEN
            IF (LNONL) THEN
              CALL GET_BLOCAL(IP, AC1, BLOC)
              CALL GET_IMATRA_IMATDA(IP, AC1, IMATRA, IMATDA)
              ASPAR_DIAG = ASPAR_DIAG + IMATDA
              eSum = BLOC + IMATRA
            ELSE
              eSum = B_JAC(:,:,IP)
            END IF
          ELSE
            CALL GET_BLOCAL(IP, AC1, eSum)
          END IF
!          eSum(:,:)= AC1(:,:,IP) * SI(IP)
!#ifdef DEBUG
!          sumESUMpre = sumESUMpre + sum(abs(eSum))
!#endif 
!          DO ID=1,NANG
!            eSum(ID,:) = eSum(ID,:) * IOBPD(ID,IP) * IOBWB(IP)
!          END DO
#ifdef DEBUG
          sumESUMbefore = sumESUMbefore + sum(abs(eSum))
#endif 
#ifdef DEBUG
!          WRITE(740+MyRankGlobal,*) 'Before ASPAR_JAC sums'
!          FLUSH(740+MyRankGlobal)
#endif 
          DO J=IA(IP),IA(IP+1)-1 
            IF (J .ne. I_DIAG(IP)) eSum = eSum - ASPAR_JAC(:,:,J) * AC2(:,:,JA(J))
          END DO
#ifdef DEBUG
!          WRITE(740+MyRankGlobal,*) 'After ASPAR_JAC sums'
!          FLUSH(740+MyRankGlobal)
!          WRITE(740+MyRankGlobal,*) 'REFRA_METHOD=', REFRA_METHOD
!          FLUSH(740+MyRankGlobal)
!          WRITE(740+MyRankGlobal,*) 'IREFRA=', IREFRA
!          FLUSH(740+MyRankGlobal)
#endif 
          IF (REFRA_METHOD .eq. 1) THEN
            IF (IREFRA .NE. 0) THEN
              DO K=1,NANG
                KP1 = K+1
                IF (KP1.GT.NANG) KP1 = 1
                KM1 = K-1
                IF (KM1.LT.1) KM1 = NANG
                DO M=1,NFRE
                  eSum(K,M) = eSum(K,M) - DTP_I(K,M,IJ)*AC2(KP1,M,IP)
                  eSum(K,M) = eSum(K,M) - DTM_I(K,M,IJ)*AC2(KM1,M,IP)
                END DO
              END DO
              DO K=1,NANG
                DO M=1,NFRE
                  MP1 = MIN(NFRE,M+1)
                  MM1 = MAX(1,M-1)
                  eSum(K,M) = eSum(K,M) - DOP_I(K,M,IJ)*AC2(K,MP1,IJ)
                  eSum(K,M) = eSum(K,M) - DOM_I(K,M,IJ)*AC2(K,MM1,IJ)
                END DO
              END DO
            END IF
          END IF
          IF (REFRA_METHOD .eq. 2) THEN
            IF ((IREFRA .eq. 1).or.(IREFRA.eq.3)) THEN
              CAD=CAD_THE(:,:,IP)
              CP_THE = MAX(ZERO,CAD)
              CM_THE = MIN(ZERO,CAD)
              eFact=(DT4D/DDIR)*SI(IP)
              DO ID=1,NANG
                ID1 = ID_PREV(ID)
                ID2 = ID_NEXT(ID)
                eSum(ID,:) = eSum(ID,:) + eFact*CP_THE(ID1,:)*AC2(ID1,:,IP)
                eSum(ID,:) = eSum(ID,:) - eFact*CM_THE(ID2,:)*AC2(ID2,:,IP)
              END DO
            END IF
            IF ((IREFRA.eq.2).or.(IREFRA.eq.3)) THEN
              CAS=CAS_SIG(:,:,IP)
              CP_SIG = MAX(ZERO,CAS)
              CM_SIG = MIN(ZERO,CAS)
              eFact=DT4F*SI(IP)
              DO ID=1,MDC
                DO IS=2,MSC
                  eSum(ID,IS)=eSum(ID,IS) + eFact*(CP_SIG(IS-1,ID)/DS_INCR(IS-1))*AC2(ID,IS-1,IP)
                END DO
                DO IS=1,MSC-1
                  eSum(ID,IS)=eSum(ID,IS) - eFact*(CM_SIG(IS+1,ID)/DS_INCR(IS))*AC2(ID,IS+1,IP)
                END DO
              END DO
            END IF
          END IF
#ifdef DEBUG
!          WRITE(740+MyRankGlobal,*) 'After the refraction business'
!          FLUSH(740+MyRankGlobal)
#endif 
!          WRITE(740+MyRankGlobal,*) 'IP=', IP, ' sum(ASPAR_diag)=', sum(ASPAR_diag)
#ifdef DEBUG
          sumESUMafter = sumESUMafter + sum(abs(eSum))
#endif 
          eSum=eSum/ASPAR_DIAG
#ifdef DEBUG
!          WRITE(740+MyRankGlobal,*) 'After division by ASPAR_DIAG'
!          FLUSH(740+MyRankGlobal)
!          WRITE(740+MyRankGlobal,*) 'LLIMT=', LLIMT
!          FLUSH(740+MyRankGlobal)
#endif 
          IF (LLIMT) CALL ACTION_LIMITER_LOCAL(IP,eSum,acloc)
          IF (BLOCK_GAUSS_SEIDEL) THEN
            AC2(:,:,IP)=eSum
          ELSE
            U_JACOBI(:,:,IP)=eSum
          END IF
!          WRITE(740+MyRankGlobal,*) 'LCHKCONV=', LCHKCONV
          IF (LCHKCONV) THEN
            Sum_new = sum(eSum)
!            WRITE(740+MyRankGlobal,*) 'IP=', IP, 'Sum_new=', Sum_new
            IF (Sum_new .gt. thr) THEN
              DiffNew=sum(abs(ACLOC - eSum))
!              WRITE(740+MyRankGlobal,*) 'DiffNew=', DiffNew
              p_is_converged = DiffNew/Sum_new
            ELSE
              p_is_converged = zero
            END IF
!            WRITE(740+MyRankGlobal,*) 'p_is_converged=', p_is_converged
!            WRITE(740+MyRankGlobal,*) 'JGS_DIFF_solverthr=', JGS_DIFF_SOLVERTHR
            IF (p_is_converged .lt. JGS_DIFF_SOLVERTHR) is_converged=is_converged+1
          ENDIF
        END DO
        IF (LCHKCONV) THEN
          isend(1)=is_converged
          CALL MPI_ALLREDUCE(isend, itmp, 1, MPI_INTEGER, MPI_SUM, COMM, ierr)
          is_converged = itmp(1)
          p_is_converged = (real(np_global) - real(is_converged))/real(np_global) * 100.
        ENDIF
        IF (BLOCK_GAUSS_SEIDEL) THEN
          CALL exchange(AC2)
        ELSE
          CALL exchange(U_JACOBI)
        END IF
        IF (.NOT. BLOCK_GAUSS_SEIDEL) THEN
          AC2 = U_JACOBI
        ENDIF
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'sumESUMpre = ', sumESUMpre
        WRITE(740+MyRankGlobal,*) 'sumESUMbefore = ', sumESUMbefore
        WRITE(740+MyRankGlobal,*) 'sumESUMafter  = ', sumESUMafter
#endif 
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) ' After iteration, sum(AC2)=', sum(AC2)
        FLUSH(740+MyRankGlobal)
#endif 
!
! The termination criterions several can be chosen
!
#ifdef DEBUG
        WRITE(740+MyRankGlobal,'(A10,3I10,E30.20,F10.5)') 'solver', nbiter, is_converged, (np_global - is_converged), p_is_converged, pmin
#endif
        !
        ! Number of iterations. If too large the exit.
        !
        nbIter=nbIter+1
        IF (nbiter .eq. maxiter) THEN
          EXIT
        ENDIF
        !
        ! Check via number of converged points
        !
        IF (LCHKCONV) THEN
          IF (p_is_converged .le. pmin) EXIT
        ENDIF
        !
        ! Check via the norm
        !
        IF (L_SOLVER_NORM) THEN
          Norm_L2=0
          DO IP=1,NP_RES
            IJ = IP
            ASPAR_DIAG=ASPAR_JAC(:,:,I_DIAG(IP))
            IF (SOURCE_IMPL) THEN
              IF (LNONL) THEN
                CALL GET_BLOCAL(IP, AC1, BLOC)
                CALL GET_IMATRA_IMATDA(IP, AC1, IMATRA, IMATDA)
                ASPAR_DIAG = ASPAR_DIAG + IMATDA
                eSum = BLOC + IMATRA
              ELSE
                eSum = B_JAC(:,:,IP)
              END IF
            ELSE
              CALL GET_BLOCAL(IP, AC1, eSum)
            END IF
            DO J=IA(IP),IA(IP+1)-1
              idx=JA(J)
              IF (J .eq. I_DIAG(IP)) THEN
                eSum=eSum - ASPAR_DIAG*AC2(:,:,idx)
              ELSE
                eSum=eSum - ASPAR_JAC(:,:,J)*AC2(:,:,idx)
              END IF
            END DO
            IF (IREFRA .NE. 0) THEN
              DO K=1,NANG
                KP1 = K+1
                IF (KP1.GT.NANG) KP1 = 1
                KM1 = K-1
                IF (KM1.LT.1) KM1 = NANG
                DO M=1,NFRE
                  eSum(K,M) = eSum(K,M) - DTP_I(K,M,IJ)*AC2(KP1,M,IJ)
                  eSum(K,M) = eSum(K,M) - DTM_I(K,M,IJ)*AC2(KM1,M,IJ)
                END DO
              END DO
              DO K=1,NANG
                DO M=1,NFRE
                  MP1 = MIN(NFRE,M+1)
                  MM1 = MAX(1,M-1)
                  eSum(K,M) = eSum(K,M) - DOP_I(K,M,IJ)*AC2(K,MP1,IJ)
                  eSum(K,M) = eSum(K,M) - DOM_I(K,M,IJ)*AC2(K,MM1,IJ)
                END DO
              END DO
            END IF
            Norm_L2 = Norm_L2 + (eSum**2)
            Norm_LINF = max(Norm_LINF, abs(eSum))
          END DO
          CALL MPI_ALLREDUCE(Norm_LINF, Norm_LINF_gl, NFRE*NANG,mpi_real8,MPI_MAX,comm,ierr)
          CALL MPI_ALLREDUCE(Norm_L2, Norm_L2_gl, NFRE*NANG, mpi_real8,MPI_SUM,comm,ierr)
          MaxNorm=maxval(Norm_L2_gl)
          SumNorm=sum(Norm_L2_gl)
          IF (sqrt(SumNorm) .le. WAE_SOLVERTHR) THEN
            EXIT
          END IF
       END IF
      END DO
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'Leaving SOLVER sum(FLNEW)=', sum(FLNEW(1:MNP,:,:))
#endif
      DO ip = 1,mnp
        FLsing=MAX(ZERO,AC2(:,:,IP))
        CALL REFRACTION_FREQSHIFT_EXPLICIT_SINGLE(FLsing, ip)
        FLNEW(IP,:,:) = REAL(FLsing,JWRB)
      END DO
      END SUBROUTINE IMPLICIT_N_SCHEME_BLOCK
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPLICIT_N_SCHEME_VECTOR_HPCF(FL1,FLNEW)
      USE yowMpiModule
      USE YOWUNPOOL
      USE yowpd, ONLY : comm
      USE yowexchangemodule
      IMPLICIT NONE
!*--********************************************************************
! calls       AC2      CCON     CG       COSTH    CURTXY
!             EXCHANGE_P4D_WWM  IEN      IE_CELL2 INVSPHTRANS
!             IOBPD    IOBWB    MPI_ALLREDUCE     MY_WTIME POS_CELL2
!             SI       SINTH
! called by   ** NOTHING **
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  CFLXY    COMM     CX       CY       DT4A     DT4AI
!             DTMAX_EXP         DTMAX_GLOBAL_EXP  DTMAX_GLOBAL_EXP_LOC
!             FL11     FL111    FL112    FL12     FL21     FL211
!             FL212    FL22     FL31     FL311    FL312    FL32
!             FLALL    I        I1       I2       I3       ID       IE
!             IERR     INE      IP       IPOS     IS       IT
!             ITER_MAX KELEM    KKSUM    KMAX     KSUM     LAMBDA
!             LCALC    MDC      MNP      MPI_MIN  MSC      N        NI
!             NP_RES   ONE      ONEHALF  ONESIXTH REST     RKIND
!             RTYPE    ST       STAT     THR      TIME1    TIME2
!             TWO      U        U3       UIP      UTILDE3  VERYLARGE
!             ZERO
! uses PARAMs *** NONE ****
!*++********************************************************************
!
! local INTEGER
!
      INTEGER(KIND=JWIM) :: IP, IE, IT, IS, ID
      INTEGER(KIND=JWIM) :: I1, I2, I3
      INTEGER(KIND=JWIM) :: NI(3), I, IPOS, IERR
!
! local double
!
      REAL(KIND=JWRB), INTENT(IN) :: FL1(:,:,:)
      REAL(KIND=JWRB), INTENT(INOUT) :: FLNEW(:,:,:)
      REAL(KIND=JWRU)  :: DTMAX_GLOBAL_EXP, DTMAX_EXP
      REAL(KIND=JWRU)  :: DTMAX_GLOBAL_EXP_LOC
      REAL(KIND=JWRU)  :: REST, CFLXY, DT4AI
      REAL(KIND=JWRU)  :: LAMBDA(2)
      REAL(KIND=JWRU)  :: FL11,FL12,FL21,FL22,FL31,FL32
      REAL(KIND=JWRU)  :: U3(3), UIP(MNP)
      REAL(KIND=JWRU)  :: KKSUM, ST, N
      REAL(KIND=JWRU)  :: CX(3), CY(3)
      REAL(KIND=JWRU)  :: U(NANG,NFRE,MNP)
      REAL(KIND=JWRU)  :: FLALL(3)
      REAL(KIND=JWRU)  :: KELEM(3)
      REAL(KIND=JWRU)  :: FL111, FL112, FL211, FL212, FL311, FL312
      REAL(KIND=JWRU)  :: UTILDE3
      REAL(KIND=JWRU)  :: KSUM(MNP), KMAX(MNP)
      REAL(KIND=JWRU)  :: TIME1, TIME2
!
!        Calculate phase speeds for the certain spectral comp1.d0nt ...
!
      flall = 0.0_JWRU
      kelem = 0.0_JWRU
      kksum = 0.0_JWRU
      st = 0.0_JWRU
      n = 0.0_JWRU
 
      IF ( lcalc ) THEN
         DO id = 1 , mdc
            DO is = 1 , msc
               kmax = 0.0_JWRU
               ksum = 0.0_JWRU
               DO ip = 1 , mnp
                  DO i = 1 , CCON(ip)
                     ie = IE_CELL2(ip,i)
                     ipos = POS_CELL2(ip,i)
! get node indices from the element table ...
                     ni = ine(:,ie)
! estimate speed in WAE
                     cx = (CG(is,ni)*SINTH(id)+CURTXY(1,ni))              &
     &                    *1.0_JWRU/(DEGRAD*REARTH*COS(YP(IP)*DEGRAD))
                     cy = (CG(is,ni)*COSTH(id)+CURTXY(2,ni))              &
     &                    *1.0_JWRU/(DEGRAD*REARTH)

! upwind indicators
                     lambda(1) = 1.0_JWRU/6.0_JWRU*SUM(cx)
                     lambda(2) = 1.0_JWRU/6.0_JWRU*SUM(cy)
! flux jacobians
                     kelem(1) = MAX(0.0_JWRU,lambda(1)*IEN(1,ie)+lambda(2)    &
     &                          *IEN(2,ie))                                             ! K
                     kelem(2) = MAX(0.0_JWRU,lambda(1)*IEN(3,ie)+lambda(2)    &
     &                          *IEN(4,ie))
                     kelem(3) = MAX(0.0_JWRU,lambda(1)*IEN(5,ie)+lambda(2)    &
     &                          *IEN(6,ie))
 
                     ksum(ip) = ksum(ip) + kelem(ipos)
!2do check if also stable when abs removed
                     IF ( kelem(ipos).GT.kmax(ip) ) kmax(ip)              &
     &                    = kelem(ipos)
                  ENDDO
               ENDDO
               dtmax_global_exp = 10.E10
               dtmax_global_exp_loc = 10.E10
               DO ip = 1 , mnp 
                  dtmax_exp = SI(ip)/MAX(thr,ksum(ip))
                  dtmax_global_exp_loc = MIN(dtmax_global_exp_loc,        &
     &               dtmax_exp)
               ENDDO
               CALL MPI_ALLREDUCE(dtmax_global_exp_loc,dtmax_global_exp,  &
     &                            1,MPI_REAL8,mpi_min,comm,ierr)
               cflxy = dt4a/dtmax_global_exp
               rest = ABS(MOD(cflxy,1.0_JWRU))
               IF ( rest.LT.thr ) THEN
                  ITER_EXP(id,is) = ABS(NINT(cflxy))
               ELSEIF ( rest.GT.thr .AND. rest.LT.0.5_JWRU ) THEN
                  ITER_EXP(id,is) = ABS(NINT(cflxy)) + 1
               ELSE
                  ITER_EXP(id,is) = ABS(NINT(cflxy))
               ENDIF
            ENDDO   ! IS
         ENDDO    ! ID
      ENDIF
 
      DO ip = 1 , mnp
        u(:,:,ip) = fl1(ip,:,:)
      ENDDO
 
      CALL exchange(u)
 
      iter_max = MAXVAL(ITER_EXP)
      dt4ai = dt4a/iter_max
 
      DO it = 1 , iter_max
         DO is = 1 , msc
            DO id = 1 , mdc
               uip = u(id,is,:)
!!$OMP DO PRIVATE(IP,I,IE,IPOS)
               DO ip = 1 , mnp
                  st = 0.0_JWRU
                  DO i = 1 , CCON(ip)
! get element and the position of IP in the element index
                     ie = IE_CELL2(ip,i)
                     ipos = POS_CELL2(ip,i)
! get node indices from the element table ...
                     ni = ine(:,ie)
                     i1 = ni(1)
                     i2 = ni(2)
                     i3 = ni(3)
! estimate speed in WAE
                     cx = (CG(is,ni)*SINTH(id)+CURTXY(1,ni))             &
     &                    *1.0_JWRU/(DEGRAD*REARTH*COS(YP(IP)*DEGRAD))
                     cy = (CG(is,ni)*COSTH(id)+CURTXY(2,ni))             &
     &                    *1.0_JWRU/(DEGRAD*REARTH)

! flux integration using simpson rule ...
                     fl11 = cx(2)*IEN(1,ie) + cy(2)*IEN(2,ie)
                     fl12 = cx(3)*IEN(1,ie) + cy(3)*IEN(2,ie)
                     fl21 = cx(3)*IEN(3,ie) + cy(3)*IEN(4,ie)
                     fl22 = cx(1)*IEN(3,ie) + cy(1)*IEN(4,ie)
                     fl31 = cx(1)*IEN(5,ie) + cy(1)*IEN(6,ie)
                     fl32 = cx(2)*IEN(5,ie) + cy(2)*IEN(6,ie)
                     fl111 = 2.0_JWRU*fl11 + fl12
                     fl112 = 2.0_JWRU*fl12 + fl11
                     fl211 = 2.0_JWRU*fl21 + fl22
                     fl212 = 2.0_JWRU*fl22 + fl21
                     fl311 = 2.0_JWRU*fl31 + fl32
                     fl312 = 2.0_JWRU*fl32 + fl31
! upwind indicators
                     lambda(1) = 1.0_JWRU/6.0_JWRU*SUM(cx)
                     lambda(2) = 1.0_JWRU/6.0_JWRU*SUM(cy)
 
! flux jacobians
                     kelem(1) = lambda(1)*IEN(1,ie) + lambda(2)       &
     &                          *IEN(2,ie)                                   ! K
                     kelem(2) = lambda(1)*IEN(3,ie) + lambda(2)       &
     &                          *IEN(4,ie)
                     kelem(3) = lambda(1)*IEN(5,ie) + lambda(2)       &
     &                          *IEN(6,ie)
 
! inverse of the positive sum ...
                     n = -1.0_JWRU/MIN(-thr,SUM(MIN(0.0_JWRU,kelem)))       ! N
 
! positive flux jacobians
                     kelem(1) = MAX(0.0_JWRU,kelem(1))
                     kelem(2) = MAX(0.0_JWRU,kelem(2))
                     kelem(3) = MAX(0.0_JWRU,kelem(3))
 
! simpson integration last step ...
                     flall(1) = (fl311+fl212)*1.0_JWRU/6.0_JWRU + kelem(1)
                     flall(2) = (fl111+fl312)*1.0_JWRU/6.0_JWRU + kelem(2)
                     flall(3) = (fl211+fl112)*1.0_JWRU/6.0_JWRU + kelem(3)
 
! flux conserving upwind contribution
                     u3 = uip(ni)
 
                     utilde3 = n*(flall(1)*u3(1)+flall(2)*u3(2)+flall(3)    &
     &                         *u3(3))
 
! coefficient for the integration in time
                     st = st + kelem(ipos)*(uip(ip)-utilde3)
 
                  ENDDO
! time stepping ...
                  uip(ip) = MAX(0.0_JWRU,uip(ip)-dt4ai/SI(ip)*st*IOBWB(ip))    &
     &                      *IOBPD(id,ip)
               ENDDO
               u(id,is,:) = uip
            ENDDO
         ENDDO
         CALL exchange(u)
      ENDDO
 
      DO ip = 1 , mnp
         flnew(ip,:,:) = REAL(u(:,:,ip),JWRB)
      ENDDO
      END SUBROUTINE EXPLICIT_N_SCHEME_VECTOR_HPCF
!**********************************************************************
!*                                                                    *
!**********************************************************************
      END MODULE UNWAM 
