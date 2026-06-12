SUBROUTINE VDFINCR(KIDIA  , KFDIA  , KLON   , KLEV   , KTOP   , PTMST  , &
 & PUM1   , PVM1   , PSLGM1 , PTM1   , PQTM1  , PAPHM1 , PGEOM1 , &
 & PCFM   , PUDIF  , PVDIF  , PSLGDIF, PQTDIF , &
 & PVOM   , PVOL   , PSLGE  , PQTE   , PSLGEWODIS, &
 & PVDIS  , PSTRTU , PSTRTV)  
! (C) Copyright 1990- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!     ------------------------------------------------------------------

!**   *VDFINCR* - INCREMENTS U,V,T AND Q-TENDENCIES; COMPUTE MULTILEVEL
!                 FLUXES AND DISSIPATION.

!     A.C.M. BELJAARS  18/01/90   DERIVED FROM VDIFF (CY34)
!     A.C.M. BELJAARS  26/03/90   OBUKHOV-L UPDATE
!     M. Ko"hler        3/12/2004 Conserved variables (qt and slg)
!     P. Lopez         02/06/2005 Removed option for linearized
!                                 physics (now called separately)
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

!     OUTPUT PARAMETER (INTEGER):

!     *KTOP*         FIRST LEVEL INDEX WITHOUT ZERO-DIFFUSION

!     INPUT PARAMETERS (REAL):

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)
!     *PUM1*         X-VELOCITY COMPONENT AT T-1
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1
!     *PSLGM1*       GENERALIZED LIQUID WATER STATIC ENERGY (SLG) AT T-1
!     *PTM1*         TEMPERATURE AT T-1
!     *PQTM1*        TOTAL WATER AT T-1
!     *PAPHM1*       PRESSURE AT T-1
!     *PGEOM1*       GEOPOTENTIAL AT T-1
!     *PCFM*         PROP. TO EXCH. COEFF. FOR MOMENTUM (C-STAR IN DOC.)
!     *PUDIF*        U-DOUBLE TILDE DEVIDED BY ALFA
!     *PVDIF*        V-DOUBLE TILDE DEVIDED BY ALFA

!     UPDATED PARAMETERS (REAL):

!     *PSLGDIF*      SLG-DOUBLE TILDE DEVIDED BY ALFA (ON ENTRY)
!                    SLG-SINGLE TILDE                 (ON EXIT)
!     *PQTDIF*       QT-DOUBLE TILDE DEVIDED BY ALFA  (ON ENTRY)
!                    QT-SINGLE TILDE                  (ON EXIT)
!     *PVOM*         U-TENDENCY
!     *PVOL*         V-TENDENCY
!     *PSLGE*        SLG-TENDENCY
!     *PQTE*         QT-TENDENCY

!     OUTPUT PARAMETERS (REAL):

!     *PVDIS*        DISSIPATION
!     *PSTRTU*       TURBULENT FLUX OF U-MOMEMTUM         KG*(M/S)/(M2*S)
!     *PSTRTV*       TURBULENT FLUX OF V-MOMEMTUM         KG*(M/S)/(M2*S)
!     *PSLGEWODIS*   SLG-TENDENCY MINUS DISSIPATION

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB       ,JPRD
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RG
USE YOEVDF   , ONLY : RVDIFTS

IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSLGDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQTDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSLGE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQTE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSLGEWODIS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRTU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRTV(KLON,0:KLEV) 

!*         0.2    LOCAL VARIABLES

REAL(KIND=JPRB) ::    ZVIDIS(KLON)

INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPRB) ::    ZCONS1, ZCONS2, &
                    & ZDUDT, ZDVDT, ZLODIS, ZTPFAC2, ZTPFAC3, ZTPFAC4  
REAL(KIND=JPHOOK) ::    ZHOOK_HANDLE


!     ------------------------------------------------------------------

!*         1.     INITIALIZE CONSTANTS
!                 --------------------

IF (LHOOK) CALL DR_HOOK('VDFINCR',0,ZHOOK_HANDLE)
ZTPFAC2 = 1.0_JPRB/RVDIFTS
ZTPFAC3 = 1.0_JPRB-ZTPFAC2
ZTPFAC4 = 1.0_JPRB+ZTPFAC3

ZCONS1  = 1.0_JPRB/PTMST
ZCONS2  = 1.0_JPRB/(RG*PTMST)


!     ------------------------------------------------------------------

!*         2.    COMPUTE TENDENCIES AND BUDGETS
!                ------------------------------

DO JL=KIDIA,KFDIA
  ZVIDIS(JL)=0.0_JPRB
ENDDO

!*         2.1  VERTICAL LOOP

DO JK=KTOP,KLEV
  DO JL=KIDIA,KFDIA

    ZDUDT          = ( PUDIF(JL,JK) - ZTPFAC2 * PUM1(JL,JK) ) * ZCONS1
    ZDVDT          = ( PVDIF(JL,JK) - ZTPFAC2 * PVM1(JL,JK) ) * ZCONS1
    ZLODIS         = 0.5_JPRB * &
     & (  (ZTPFAC2*PUM1(JL,JK) - PUDIF(JL,JK) + PTMST*PVOM(JL,JK))&
     &   *(ZTPFAC4*PUM1(JL,JK) + PUDIF(JL,JK))&
     & +  (ZTPFAC2*PVM1(JL,JK) - PVDIF(JL,JK) + PTMST*PVOL(JL,JK))&
     &   *(ZTPFAC4*PVM1(JL,JK) + PVDIF(JL,JK)) )
    PVOM(JL,JK)    = ZDUDT
    PVOL(JL,JK)    = ZDVDT
    ZVIDIS(JL)     = ZVIDIS(JL) + ZLODIS * (PAPHM1(JL,JK)-PAPHM1(JL,JK-1))

    PQTDIF(JL,JK)  =   PQTDIF(JL,JK) + ZTPFAC3 * PQTM1(JL,JK)
    PQTE(JL,JK)    = ( PQTDIF(JL,JK) - PQTM1(JL,JK) ) * ZCONS1

!----------------------------------------------------
!  MPBL - Liquid water static energy
!    Input:  PSLGDIF liquid static energy first guess
!            PQTDIF  total water first guess
!    Output: PSLGE   liquid static energy tendency
!----------------------------------------------------

      PSLGDIF(JL,JK)   =   PSLGDIF(JL,JK) + ZTPFAC3 * PSLGM1(JL,JK)
      PSLGE(JL,JK)     = ( PSLGDIF(JL,JK) + ZLODIS - PSLGM1(JL,JK) ) * ZCONS1
      PSLGEWODIS(JL,JK)= ( PSLGDIF(JL,JK)          - PSLGM1(JL,JK) ) * ZCONS1

  ENDDO
ENDDO

!*         2.2  COMPUTE SURFACE STRESSES AND COPY DISSIPATION

DO JL=KIDIA,KFDIA
  PSTRTU(JL,KLEV) = ZCONS2*PCFM(JL,KLEV)*PUDIF(JL,KLEV)
  PSTRTV(JL,KLEV) = ZCONS2*PCFM(JL,KLEV)*PVDIF(JL,KLEV)
  PVDIS(JL)       = ZCONS2*ZVIDIS(JL)
ENDDO


IF (LHOOK) CALL DR_HOOK('VDFINCR',1,ZHOOK_HANDLE)
END SUBROUTINE VDFINCR
