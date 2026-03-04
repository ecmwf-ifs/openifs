! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE VDFDIFM(KIDIA,KFDIA,KLON,KLEV,KTOP,&
 & PTMST,PUM1,PVM1,PAPHM1,PCFM,PTODC,&
 & PSOTEU,PSOTEV,PSOC,&
 & PVOM,PVOL,PUCURR,PVCURR,PMFLX,PUUH,PVUH,&
 & PUDIF,PVDIF,PDTUTOFD,PDTVTOFD,PDTUSO,PDTVSO,PDTUVDF,PDTVVDF)
  
!     ------------------------------------------------------------------

!**   *VDFDIFM* - DOES THE IMLPLICIT CALCULATION FOR MOMENTUM DIFFUSION

!     PURPOSE
!     -------

!     SOLVE TRIDIAGONAL MATRICES FOR MOMENTUM DIFFUSION. 
!     The time step can contain the over implicit factor for 
!     extrapolation in time, or the factor for the predictor corrector 
!     method. The output is the transformed velocity (e.g. the extrapolated 
!     in time or the hat parameters in the documentation). Output also includes 
!     tendencies from TOFD, SO and VDF+TOFD.     

!     INTERFACE
!     ---------

!     *VDFDIFM* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KTOP*         INDEX FOR BOUNDARY LAYER TOP

!     INPUT PARAMETERS (REAL):

!     *PTMST*        TIME STEP including possible extrapolation in time
!     *PUM1*         X-VELOCITY COMPONENT AT T-1
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1
!     *PAPHM1*       PRESSURE AT T-1
!     *PCFM*         Rho*K/dz (C,K-STAR IN DOC.)
!     *PTODC*        TURB. OROGR. FORM DRAG (TOFD) IMPLICIT COEFFICIENTS 
!     *PSOTEU*       Explicit part of U-tendency from subgrid orography scheme    
!     *PSOTEV*       Explicit part of V-tendency from subgrid orography scheme     
!     *PSOC*         Implicit part of subgrid orography (df/dt=PSOTE-PSOC*f)
!     *PVOM*         U-TENDENCY
!     *PVOL*         V-TENDENCY
!     *PUCURR*       OCEAN CURRENT X-COMPONENT
!     *PVCURR*       OCEAN CURRENT Y-COMPONENT
!     *PMFLX*        CONVECTIVE MASS FLUX KG/M2S
!     *PUUH*         UPDRAUGHT ZONAL WIND SPEED
!     *PVUH*         UPDRAUGHT MERDIONAL WIND SPEED


!     OUTPUT PARAMETERS (REAL):

!     *PUDIF*        U-hat (velocity at new time level)
!     *PVDIF*        V-hat (velocity at new time level)
!     *PDTUTOFD*     U-Diag.Tend. due to TOFD (turb. orogr. form drag)
!     *PDTVTOFD*     V-Diag.Tend. due to TOFD
!     *PDTUSO*       U-Diag.Tend. due to SO (subgrid orography)
!     *PDTVSO*       V-Diag.Tend. due to SO
!     *PDTUVDF*      U-Diag.Tend. due to VDF (vertical diffusion) + TOFD
!     *PDTVVDF*      V-Diag.Tend. due to VDF + TOFD

!     METHOD
!     ------

!     *LU*-DECOMPOSITION AND BACK SUBSTITUTION IN ONE DOWNWARD SCAN
!     AND ONE UPWARD SCAN.

!     AUTHOR.
!     -------
!      A.C.M. BELJAARS       E.C.M.W.F.    10-11-89

!     MODIFICATIONS.
!     --------------
!      A BELJAARS              17-11-02 Turb. Orogr. Drag  
!      A.Beljaars              12/11/02 Ocean current b.c.
!      A BELJAARS              30-09-05 Include Subgr. Oro. in solver  
!      A.Beljaars              Jan 2014 Clean-up for flexible numerics
!      P.Bechtold+A. Beljaars  March 2019 Complete Recoding to same standard as vdfdifh
!      P.Bechtold              March 2019 Add mass flux transport - optional
!     ---------------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!USE YOEVDF   , ONLY : RVDIFTS
USE YOMCST   , ONLY : RG 

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTOP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTODC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOTEU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOTEV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVOM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVOL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFLX(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUUH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVUH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTUTOFD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTVTOFD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTUSO(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTVSO(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTUVDF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTVVDF(KLON,KLEV) 
!     ------------------------------------------------------------------

REAL(KIND=JPRB) ::    ZAA(KLON,KLEV) ,ZBB(KLON,KLEV) ,ZCC(KLON,KLEV) ,&
                    & ZUYY(KLON,KLEV),ZVYY(KLON,KLEV),ZGAM(KLON,KLEV),&
                    & Z1DP(KLON,KLEV), Z1BET(KLON), ZUVFLX(KLON,0:KLEV,2),&
                    & ZMFLX(KLON,0:KLEV)

INTEGER(KIND=JPIM) :: JK, JL, JKP

REAL(KIND=JPRB) :: ZHUCURR, ZCONS1,ZTMST,&
                  &ZDUDT, ZDVDT, ZHU2, ZHU3
REAL(KIND=JPRB) :: ZMF=0.0_JPRB ! Use or not convective mass flux

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('VDFDIFM',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

ZCONS1=RG*PTMST
ZTMST=1.0_JPRB/PTMST

!*         1.     FULL MODEL PHYSICS WITH MASS FLUX PBL
!                 -------------------------------------

DO JK=0,KLEV
  JKP=MAX(1,JK)
  DO JL=KIDIA,KFDIA
    ZMFLX(JL,JK) = PMFLX(JL,JK)*ZMF
    ZUVFLX (JL,JK,1) = ZMFLX(JL,JK)*PUUH(JL,JK)
    ZUVFLX (JL,JK,2) = ZMFLX(JL,JK)*PVUH(JL,JK)
  ENDDO
ENDDO

!*         1.1    SETTING OF THE MATRIX A, B AND C.

DO JK=KTOP+1,KLEV-1
  DO JL=KIDIA,KFDIA
    Z1DP(JL,JK)=ZCONS1/(PAPHM1(JL,JK)-PAPHM1(JL,JK-1))
    ZAA(JL,JK) =(-PCFM(JL,JK-1)-ZMFLX(JL,JK-1))*Z1DP(JL,JK)
    ZCC(JL,JK) =(-PCFM(JL,JK)                 )*Z1DP(JL,JK)
    ZBB(JL,JK) =1.0_JPRB+( PCFM(JL,JK-1)+PCFM(JL,JK)+ZMFLX(JL,JK) )*Z1DP(JL,JK)&
     & +(PTODC(JL,JK)+PSOC(JL,JK))*PTMST
  ENDDO
ENDDO

!          1.1a   THE SURFACE BOUNDARY CONDITION

DO JL=KIDIA,KFDIA
  Z1DP(JL,KLEV)=ZCONS1/(PAPHM1(JL,KLEV)-PAPHM1(JL,KLEV-1))
  ZCC(JL,KLEV) =0.0_JPRB
  ZAA(JL,KLEV) =(-PCFM(JL,KLEV-1)-ZMFLX(JL,JK-1))*Z1DP(JL,KLEV)
  ZBB(JL,KLEV) =1.0_JPRB+(PCFM(JL,KLEV-1)+PCFM(JL,KLEV))*Z1DP(JL,KLEV)&
   &+(PTODC(JL,KLEV)+PSOC(JL,KLEV))*PTMST  
ENDDO

!          1.1b   THE TOP BOUNDARY CONDITION    

DO JL=KIDIA,KFDIA
  Z1DP(JL,KTOP)=ZCONS1/(PAPHM1(JL,KTOP)-PAPHM1(JL,KTOP-1))
  ZAA(JL,KTOP) =0.0_JPRB
  ZCC(JL,KTOP) =         (-PCFM(JL,KTOP))*Z1DP(JL,KTOP)
  ZBB(JL,KTOP) =1.0_JPRB+( PCFM(JL,KTOP)+ZMFLX(JL,KTOP) )*Z1DP(JL,KTOP)&
                      &+(PTODC(JL,KTOP)+PSOC(JL,KTOP))*PTMST 
ENDDO

!*         1.2    SETTING OF RIGHT HAND SIDES.

DO JK=KTOP+1,KLEV-1
  DO JL=KIDIA,KFDIA
    ZUYY(JL,JK) = PUM1(JL,JK) &
     & + PTMST * (PVOM(JL,JK)+PSOTEU(JL,JK))&
     & + (ZUVFLX(JL,JK,1)-ZUVFLX(JL,JK-1,1)) * Z1DP(JL,JK)
    ZVYY(JL,JK) = PVM1(JL,JK) &
     & + PTMST * (PVOL(JL,JK)+PSOTEV(JL,JK))&
     & + (ZUVFLX(JL,JK,2)-ZUVFLX(JL,JK-1,2)) * Z1DP(JL,JK)
  ENDDO
ENDDO

!          1.2a   SURFACE

JK=KLEV
DO JL=KIDIA,KFDIA
  ZUYY(JL,JK) = PUM1(JL,JK) &
   & + PTMST * (PVOM(JL,JK)+PSOTEU(JL,JK))& 
   & - ZUVFLX(JL,JK-1,1) * Z1DP(JL,JK)
  ZVYY(JL,JK) = PVM1(JL,JK) &
   & + PTMST * (PVOL(JL,JK)+PSOTEV(JL,JK))& 
   & - ZUVFLX(JL,JK-1,2) * Z1DP(JL,JK)
ENDDO

!          1.2b   TOP

JK=KTOP
DO JL=KIDIA,KFDIA
  ZUYY(JL,JK) = PUM1(JL,JK) &
   & + PTMST * (PVOM(JL,JK)+PSOTEU(JL,JK))& 
   & + ZUVFLX(JL,JK,1) * Z1DP(JL,JK)
  ZVYY(JL,JK) = PVM1(JL,JK) &
   & + PTMST * (PVOL(JL,JK)+PSOTEV(JL,JK))& 
   & + ZUVFLX(JL,JK,2) * Z1DP(JL,JK)
ENDDO


!*         1.4    TOP LAYER ELIMINATION.

DO JL=KIDIA,KFDIA
  Z1BET(JL)=1.0_JPRB/ZBB(JL,KTOP)
  PUDIF(JL,KTOP)=ZUYY(JL,KTOP)*Z1BET(JL)
  PVDIF(JL,KTOP)=ZVYY(JL,KTOP)*Z1BET(JL)
ENDDO

!*         1.5    ELIMINATION FOR MIDDLE LAYERS.

DO JK=KTOP+1,KLEV-1
  DO JL=KIDIA,KFDIA
    ZGAM(JL,JK)=ZCC(JL,JK-1)*Z1BET(JL)
    Z1BET(JL)=1.0_JPRB/(ZBB(JL,JK)-ZAA(JL,JK)*ZGAM(JL,JK))
    PUDIF(JL,JK)=(ZUYY(JL,JK)-ZAA(JL,JK)*PUDIF(JL,JK-1))*Z1BET(JL)
    PVDIF(JL,JK)=(ZVYY(JL,JK)-ZAA(JL,JK)*PVDIF(JL,JK-1))*Z1BET(JL)
  ENDDO
ENDDO

!*         1.6    BOTTOM LAYER ELIMINATION

DO JL=KIDIA,KFDIA
  ZGAM(JL,KLEV)=ZCC(JL,KLEV-1)*Z1BET(JL)
  Z1BET(JL)=1.0_JPRB/(ZBB(JL,KLEV)-ZAA(JL,KLEV)*ZGAM(JL,KLEV))
  ZHUCURR=PCFM(JL,KLEV)*Z1DP(JL,KLEV)
  PUDIF(JL,KLEV)=(ZUYY(JL,KLEV)-ZAA(JL,KLEV)*PUDIF(JL,KLEV-1)+ZHUCURR*PUCURR(JL))*Z1BET(JL)
  PVDIF(JL,KLEV)=(ZVYY(JL,KLEV)-ZAA(JL,KLEV)*PVDIF(JL,KLEV-1)+ZHUCURR*PVCURR(JL))*Z1BET(JL)
ENDDO


!*         1.13   BACK-SUBSTITUTION.

DO JK=KLEV-1,KTOP,-1
  DO JL=KIDIA,KFDIA
    PUDIF(JL,JK)=PUDIF(JL,JK)-ZGAM(JL,JK+1)*PUDIF(JL,JK+1)
    PVDIF(JL,JK)=PVDIF(JL,JK)-ZGAM(JL,JK+1)*PVDIF(JL,JK+1)
  ENDDO
ENDDO


!*         1.6     COMPUTATION OF TENDENCIES.

PDTUTOFD(KIDIA:KFDIA,1:KTOP-1)=0.0_JPRB
PDTVTOFD(KIDIA:KFDIA,1:KTOP-1)=0.0_JPRB
PDTUSO  (KIDIA:KFDIA,1:KTOP-1)=0.0_JPRB
PDTVSO  (KIDIA:KFDIA,1:KTOP-1)=0.0_JPRB
PDTUVDF (KIDIA:KFDIA,1:KTOP-1)=0.0_JPRB
PDTVVDF (KIDIA:KFDIA,1:KTOP-1)=0.0_JPRB

DO JK=KTOP,KLEV
  DO JL=KIDIA,KFDIA

!   Compute total tendencies (dynamics + vertical diffusion + SO) 
    ZDUDT = ( PUDIF(JL,JK) - PUM1(JL,JK) ) * ZTMST
    ZDVDT = ( PVDIF(JL,JK) - PVM1(JL,JK) ) * ZTMST

!   TOFD-tendencies 
    ZHU2=-PTODC(JL,JK)
    PDTUTOFD(JL,JK)=PUDIF(JL,JK)*ZHU2
    PDTVTOFD(JL,JK)=PVDIF(JL,JK)*ZHU2

!   SO-tendencies 
    ZHU3=-PSOC(JL,JK)
    PDTUSO(JL,JK)=PUDIF(JL,JK)*ZHU3+PSOTEU(JL,JK)
    PDTVSO(JL,JK)=PVDIF(JL,JK)*ZHU3+PSOTEV(JL,JK)

!   VDF + TOFD tendencies (= total-SO)
    PDTUVDF(JL,JK)=ZDUDT-PVOM(JL,JK)-PDTUSO(JL,JK)
    PDTVVDF(JL,JK)=ZDVDT-PVOL(JL,JK)-PDTVSO(JL,JK)
  ENDDO
ENDDO


!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('VDFDIFM',1,ZHOOK_HANDLE)
END SUBROUTINE VDFDIFM
