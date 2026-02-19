! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE VDFINCRS (KIDIA,KFDIA,KLON,KLEV,KTOP,PTMST,PVDIFTS,&
 & PUM1    ,PVM1  ,PTM1   ,PQM1  , PAPHM1, &
 & PGEOM1  ,PCFM  ,PTOFDC ,PCPTGZ ,&
 & PUDIF   ,PVDIF ,PTDIF  ,PQDIF , &
 & PVOM    ,PVOL  ,PTE    ,PQE   ,&
 & PTEWODIS,PVDIS ,PSTRTU ,PSTRTV)
!     ------------------------------------------------------------------

!**   *VDFINCRS* - INCREMENTS U,V,T AND Q-TENDENCIES; COMPUTE MULTILEVEL
!                   FLUXES AND DISSIPATION.
!                   (Nonlinear version for trajectory in adjoint) 

!     J.F. MAHFOUF          E.C.M.W.F.    02/10/95

!     Adapted from

!     DERIVED FROM VDIFF (CY34) BY
!     A.C.M. BELJAARS       E.C.M.W.F.    18-1-90

!     OBUKHOV-L UPDATE      ACMB          26/03/90.
!     TOFD                  ACMB          18/12/2005

!     PURPOSE
!     -------

!     INCREMENT U,V,T AND Q; COMPUTE MULTILEVEL FLUXES AND DISSIPATION

!     INTERFACE
!     ---------

!     *VDFINCRS* IS CALLED BY *VDFMAINS*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     OUTPUT PARAMETER (INTEGER):

!     *KTOP*         FIRST LEVEL INDEX WITHOUT ZERO-DIFFUSION

!     INPUT PARAMETERS (REAL):

!     *PUM1*        X-VELOCITY COMPONENT AT T-1     (Trajectory)
!     *PVM1*        Y-VELOCITY COMPONENT AT T-1     (Trajectory)
!     *PTM1*        TEMPERATURE AT T-1              (Trajectory)
!     *PQM1*        SPECIFIC HUMUDITY AT T-1        (Trajectory)
!     *PAPHM1*      PRESSURE AT T-1                 (Trajectory)
!     *PGEOM1*      GEOPOTENTIAL AT T-1             (Trajectory)
!     *PCFM*        PROP. TO EXCH. COEFF. FOR MOMENTUM (C-STAR IN DOC.)
!                                                    (Trajectory)
!     *PTOFDC*      TURBULENT OROGRAPHIC DRAG COEFFICIENT
!                                                    (Trajectory)
!     *PCPTGZ*      DRY STATIC ENERGY               (Trajectory)
!     *PUDIF*       U-DOUBLE TILDE DEVIDED BY ALFA  (Trajectory)
!     *PVDIF*       V-DOUBLE TILDE DEVIDED BY ALFA  (Trajectory)

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)

!     UPDATED PARAMETERS (REAL):

!     *PTDIF*       S-DOUBLE TILDE DEVIDED BY ALFA (ON ENTRY) (Trajectory)
!                    S-SINGLE TILDE                 (ON EXIT)  (Trajectory)
!     *PQDIF*       Q-DOUBLE TILDE DEVIDED BY ALFA (ON ENTRY) (Trajectory)
!                    Q-SINGLE TILDE                 (ON EXIT)  (Trajectory)
!     *PVOM*        U-TENDENCY
!     *PVOL*        V-TENDENCY
!     *PTE*         T-TENDENCY
!     *PQE*         Q-TENDENCY

!     OUTPUT PARAMETERS (REAL):

!     *PVDIS*       DISSIPATION                      (Trajectory)
!     *PSTRTU*      TURBULENT FLUX OF U-MOMEMTUM     (Trajectory)
!     *PSTRTV*      TURBULENT FLUX OF V-MOMEMTUM     (Trajectory)
!     *PTEWODIS     T-TENDENCY MINUS DISSIPATION     (Trajectory)

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     Modified by
!     -----------

!     P. LOPEZ    E.C.M.W.F.  25/02/05 (nonlinear version for
!                                       trajectory in adjoint)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RG       ,RCPD
USE YOETHF   , ONLY : RVTMP2

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTOP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVDIFTS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTOFDC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTGZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTEWODIS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRTU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRTV(KLON,0:KLEV) 

!*    LOCAL STORAGE
!     ----- -------

REAL(KIND=JPRB) ::    ZVIDIS(KLON),ZTAUU(KLON),ZTAUV(KLON)

INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPRB) :: Z1S, ZCONS13, ZCONS1, ZDQDT,&
 & ZDTDT, ZDTWODIS, ZDUDT, ZDVDT, ZLODIS, ZTPFAC2, &
 & ZTPFAC3, ZTPFAC4, ZDP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.     INITIALIZE CONSTANTS
!               ---------- ----------

IF (LHOOK) CALL DR_HOOK('VDFINCRS',0,ZHOOK_HANDLE)
ZTPFAC2 = 1.0_JPRB/PVDIFTS
ZTPFAC3 = 1.0_JPRB-ZTPFAC2
ZTPFAC4 = 1.0_JPRB+ZTPFAC3

ZCONS13 = 1.0_JPRB/PTMST
ZCONS1 = 1.0_JPRB/(RG*PTMST)

!     ------------------------------------------------------------------

!*       2.    COMPUTE TENDENCIES AND BUDGETS
!              ------- ---------- --- -------

DO JL = KIDIA, KFDIA
  ZVIDIS(JL) = 0.0_JPRB
  ZTAUU(JL)=0.0_JPRB
  ZTAUV(JL)=0.0_JPRB
ENDDO

!***  
DO JK = KTOP, KLEV
!***      
  DO JL = KIDIA, KFDIA
    ZDUDT = ZCONS13*(PUDIF(JL,JK)-ZTPFAC2*PUM1(JL,JK))
    ZDVDT = ZCONS13*(PVDIF(JL,JK)-ZTPFAC2*PVM1(JL,JK))
    ZLODIS = 0.5_JPRB*((ZTPFAC2*PUM1(JL,JK)-PUDIF(JL,JK)+&
     & PTMST*PVOM(JL,JK))*(ZTPFAC4*PUM1(JL,JK)+&
     & PUDIF(JL,JK))+(ZTPFAC2*PVM1(JL,JK)-&
     & PVDIF(JL,JK)+PTMST*PVOL(JL,JK))*&
     & (ZTPFAC4*PVM1(JL,JK)+PVDIF(JL,JK)))  
    PVOM(JL,JK) = ZDUDT
    PVOL(JL,JK) = ZDVDT
    ZDP         = PAPHM1(JL,JK)-PAPHM1(JL,JK-1)
    ZVIDIS(JL)  = ZVIDIS(JL) + ZLODIS*ZDP

    PQDIF(JL,JK) = PQDIF(JL,JK)+ZTPFAC3*PQM1(JL,JK)
    ZDQDT = ZCONS13*(PQDIF(JL,JK)-PQM1(JL,JK))
    PQE(JL,JK) = ZDQDT

    PTDIF(JL,JK) = PTDIF(JL,JK)+ZTPFAC3*PCPTGZ(JL,JK)
    ZDTDT = ZCONS13*((PTDIF(JL,JK)+ZLODIS-&
     & PGEOM1(JL,JK))/(RCPD*(1.0_JPRB+RVTMP2*PQDIF(JL,JK)))-&
     & PTM1(JL,JK))  
    PTE(JL,JK) = ZDTDT
    ZDTWODIS = ZCONS13*((PTDIF(JL,JK)-PGEOM1(JL,JK))/&
     & (RCPD*(1.0_JPRB+RVTMP2*PQDIF(JL,JK)))-PTM1(JL,JK))  
    PTEWODIS(JL,JK) = ZDTWODIS
    
!   INTEGRATE TOD TENDENCIES TO GET SURFACE STRESS
    ZTAUU(JL)=PUDIF(JL,JK)*PTOFDC(JL,JK)*ZDP*ZCONS1+ZTAUU(JL)
    ZTAUV(JL)=PVDIF(JL,JK)*PTOFDC(JL,JK)*ZDP*ZCONS1+ZTAUV(JL)

  ENDDO

!***
ENDDO
!***  

!*         2.2     COMPUTE SURFACE STRESSES AND COPY DISSIPATION

DO JL = KIDIA, KFDIA
  Z1S = ZCONS1*PCFM(JL,KLEV)
  PSTRTU(JL,KLEV) = Z1S*PUDIF(JL,KLEV)+ZTAUU(JL)
  PSTRTV(JL,KLEV) = Z1S*PVDIF(JL,KLEV)+ZTAUV(JL)
  PVDIS(JL) = ZCONS1*ZVIDIS(JL)
ENDDO

IF (LHOOK) CALL DR_HOOK('VDFINCRS',1,ZHOOK_HANDLE)
END SUBROUTINE VDFINCRS

