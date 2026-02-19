! (C) Copyright 1990- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
#ifdef RS6K
@PROCESS HOT(NOVECTOR) 
#endif
SUBROUTINE VDFINCR(KIDIA  , KFDIA  , KLON   , KLEV  , KTOP  ,PTMST  , &
 & PUM1   , PVM1   , PSLGM1 , PTM1   , PQTM1  , PAPHM1 , &
 & PDTUTOFD,PDTVTOFD, PDTUSO, PDTVSO , PDTUVDF , PDTVVDF ,&
 & PVOM   , PVOL   , PSLGE  , PQTE   ,& 
 & PSLGEVDF,PQTEVDF,PSLGEWODIS, &
 & PVDIS  , PVDISG, PSTRTU , PSTRTV , PSTRSOU , PSTRSOV ,& 
 & PSTRTOFDU,PSTRTOFDV,PTOFDU , PTOFDV, &
 & PDISGW3D)  

!**   *VDFINCR* - INCREMENTS U,V,T AND Q-TENDENCIES; COMPUTE MULTILEVEL
!                 FLUXES AND DISSIPATION.

!     PURPOSE
!     -------

!     INCREMENT U,V,T AND Q; COMPUTE MULTILEVEL FLUXES AND DISSIPATION

!     INTERFACE
!     ---------

!     *VDFINCR* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KTOP*         FIRST LEVEL INDEX WITHOUT ZERO-DIFFUSION


!     INPUT PARAMETERS (REAL):

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)
!     *PUM1*         X-VELOCITY COMPONENT AT T-1
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1
!     *PSLGM1*       GENERALIZED LIQUID WATER STATIC ENERGY (SLG) AT T-1
!     *PTM1*         TEMPERATURE AT T-1
!     *PQTM1*        TOTAL WATER AT T-1
!     *PAPHM1*       PRESSURE AT T-1
!     *PDTUTOFD*     U-Diag.Tend. due to TOFD (turb. orogr. form drag)
!     *PDTVTOFD*     V-Diag.Tend. due to TOFD
!     *PDTUSO*       U-Diag.Tend. due to SO (subgrid orography)
!     *PDTVSO*       V-Diag.Tend. due to SO
!     *PDTUVDF*      U-Diag.Tend. due to VDF (vertical diffusion) + TOFD
!     *PDTVVDF*      V-Diag.Tend. due to VDF + TOFD

!     UPDATED PARAMETERS (REAL):

!     *PVOM*         U-TENDENCY
!     *PVOL*         V-TENDENCY
!     *PSLGE*        SLG-TENDENCY before VDF
!     *PQTE*         QT-TENDENCY before VDF
!     *PSLGEVDF*     SLG-TENDENCY VDF only
!     *PQTEVDF*      QT-TENDENCY VDF only

!     OUTPUT PARAMETERS (REAL):

!     *PVDIS*        TURBULENT DISSIPATION
!     *PVDISG*       SUBGRID OROGRAPHY DISSIPATION
!     *PSTRTU*       FLUX OF U-MOMEMTUM Turbulent + TOFD           KG*/(M*S2)
!     *PSTRTV*       FLUX OF V-MOMEMTUM Turbulent + TOFD           KG*/(M*S2)
!     *PSTRSOU*      FLUX OF U-MOMEMTUM SO                         KG*/(M*S2)
!     *PSTRSOV*      FLUX OF V-MOMEMTUM SO                         KG*/(M*S2)
!     *PSLGEWODIS*   SLG-TENDENCY MINUS (TOTAL) DISSIPATION
!     *PTOFDU*       TOFD COMP. OF TURBULENT FLUX OF U-MOMEMTUM    KG*(M/S)/(M2*S)
!     *PTOFDV*       TOFD COMP. OF TURBULENT FLUX OF V-MOMEMTUM    KG*(M/S)/(M2*S)
!     *PDISGW3D*     3-D DISSIPATION RATE FOR STOCHASTIC PHYSIC    (M2/S2)

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     AUTHOR.
!     -------
!      A.C.M. BELJAARS  18/01/90   DERIVED FROM VDIFF (CY34)

!     MODIFICATIONS.
!     --------------
!      OCEAN CURRENT B.C.    ACMB          12/11/02.
!      M. Ko"hler    3/12/2004 Conserved variables (qt and slg)
!      A. Beljaars   4/04/2005 TURBULENT OROGR. DRAG ACMB
!      A  Beljaars  30/09/2005 Include Subgr. Oro. in solver  
!      P. Lopez     02/06/2005 Removed option for linearized physics (now called separately)
!      P.de Rosnay/G.Balsamo   07/03/2009 Offline Jacobians EKF      
!      P.de Rosnay/G.Balsamo   October 2009 commented EKF (offline in callpar)      
!      A. Beljaars    Jan 2014 Clean-up to prepare for flexible numerics
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK

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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLGM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTUTOFD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTVTOFD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTUSO(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTVSO(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTUVDF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTVVDF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSLGE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQTE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLGEVDF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQTEVDF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSLGEWODIS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDIS(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDISG(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDISGW3D(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRTU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRTV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRSOU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRSOV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRTOFDU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRTOFDV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOFDU(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOFDV(KLON) 


INTEGER(KIND=JPIM) :: JK, JL, ILEV

REAL(KIND=JPRB) ::    ZCONS1, ZCONS2, &
                    & ZVDFDIS, ZSODIS, &
                    & ZDP, ZRG, ZGDPH, &
                    & ZU1,ZU2,ZU3,ZV1,ZV2,ZV3
REAL(KIND=JPHOOK) ::    ZHOOK_HANDLE


!     ------------------------------------------------------------------

!*         1.     INITIALIZE CONSTANTS
!                 --------------------

IF (LHOOK) CALL DR_HOOK('VDFINCR',0,ZHOOK_HANDLE)

ZCONS1  = 1.0_JPRB/PTMST
ZCONS2  = 1.0_JPRB/(RG*PTMST)
ZRG     = 1.0_JPRB/RG

  ILEV=1 
!     ------------------------------------------------------------------

!*         2.    COMPUTE TENDENCIES AND BUDGETS
!                ------------------------------

DO JL=KIDIA,KFDIA
  PVDIS(JL)    =0.0_JPRB
  PVDISG(JL)   =0.0_JPRB
  PSTRTU(JL,0) =0.0_JPRB
  PSTRTV(JL,0) =0.0_JPRB
  PSTRSOU(JL,0)=0.0_JPRB
  PSTRSOV(JL,0)=0.0_JPRB
  PSTRTOFDU(JL,0)=0.0_JPRB
  PSTRTOFDV(JL,0)=0.0_JPRB
  PTOFDU(JL)   =0.0_JPRB
  PTOFDV(JL)   =0.0_JPRB
ENDDO

!*         2.1  VERTICAL LOOP

DO JK=ILEV,KLEV
  DO JL=KIDIA,KFDIA

    ZDP            = PAPHM1(JL,JK)-PAPHM1(JL,JK-1)
    ZGDPH          = -ZDP*ZRG

!   Velocity after dynamics, before VDF 
    ZU1=PUM1(JL,JK)+PVOM(JL,JK)*PTMST
    ZV1=PVM1(JL,JK)+PVOL(JL,JK)*PTMST
!   Velocity after VDF 
    ZU2=ZU1+PDTUVDF(JL,JK)*PTMST
    ZV2=ZV1+PDTVVDF(JL,JK)*PTMST
!   Velocity after SO 
    ZU3=ZU2+PDTUSO(JL,JK)*PTMST   
    ZV3=ZV2+PDTVSO(JL,JK)*PTMST   
    ZVDFDIS=0.5_JPRB*(ZU1-ZU2)*(ZU1+ZU2) + 0.5_JPRB*(ZV1-ZV2)*(ZV1+ZV2)
    ZSODIS =0.5_JPRB*(ZU2-ZU3)*(ZU2+ZU3) + 0.5_JPRB*(ZV2-ZV3)*(ZV2+ZV3)
    PDISGW3D(JL,JK) = ZSODIS

!   Integrate VDF-tendencies (including TOFD) to find VDF-stress profile
    PSTRTU(JL,JK)  = PDTUVDF(JL,JK)*ZGDPH+PSTRTU(JL,JK-1)
    PSTRTV(JL,JK)  = PDTVVDF(JL,JK)*ZGDPH+PSTRTV(JL,JK-1)

!   Integrate SO-tendencies to find SO-stress profile
    PSTRSOU(JL,JK)  = PDTUSO(JL,JK)*ZGDPH+PSTRSOU(JL,JK-1)
    PSTRSOV(JL,JK)  = PDTVSO(JL,JK)*ZGDPH+PSTRSOV(JL,JK-1)

!   Integrate TOFD-tendencies to find TOFD-stress profile
    PSTRTOFDU(JL,JK)=PDTUTOFD(JL,JK)*ZGDPH+PSTRTOFDU(JL,JK-1)
    PSTRTOFDV(JL,JK)=PDTVTOFD(JL,JK)*ZGDPH+PSTRTOFDV(JL,JK-1)

    
    PVOM(JL,JK)    = PVOM(JL,JK)+PDTUVDF(JL,JK)+PDTUSO(JL,JK)
    PVOL(JL,JK)    = PVOL(JL,JK)+PDTVVDF(JL,JK)+PDTVSO(JL,JK)
    PVDIS(JL)      = PVDIS(JL) +ZVDFDIS*ZDP
    PVDISG(JL)     = PVDISG(JL)+ZSODIS*ZDP
   
    PQTE(JL,JK)    = PQTE(JL,JK)+PQTEVDF(JL,JK)

    PSLGEWODIS(JL,JK)= PSLGE(JL,JK)+PSLGEVDF(JL,JK)
    PSLGE(JL,JK)     = PSLGE(JL,JK)+PSLGEVDF(JL,JK)+(ZVDFDIS+ZSODIS)*ZCONS1
  ENDDO
ENDDO


DO JL=KIDIA,KFDIA
  PTOFDU(JL)=PSTRTOFDU(JL,KLEV)
  PTOFDV(JL)=PSTRTOFDV(JL,KLEV)
  PVDIS(JL) =PVDIS(JL) *ZCONS2
  PVDISG(JL)=PVDISG(JL)*ZCONS2
ENDDO


IF (LHOOK) CALL DR_HOOK('VDFINCR',1,ZHOOK_HANDLE)
END SUBROUTINE VDFINCR
