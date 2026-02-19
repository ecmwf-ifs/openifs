! (C) Copyright 2015- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
MODULE SRFSN_DRIVER_MOD
CONTAINS
SUBROUTINE SRFSN_DRIVER(KIDIA   ,KFDIA   ,KLON   ,KLEVSN, PTMST, LDLAND,&
  & PSDOR, &
  ! input at T-1 prognostics
  & PSSNM1M, PTSNM1M, PRSNM1M  ,PWSNM1M, PASNM1M, &
  ! input at T-1 fluxes or constants 
  & PFRTI, PTSAM1M , PUSRF, PVSRF, PTSRF ,&
  & PSSFC, PSSFL   , PTSFC, PTSFL, &
  & PSLRFLTI, PSSRFLTI , PAHFSTI, PEVAPTI, PEVAPSNW, &
  & PWSAM1M , KSOTY, PAPRS, &
  ! input derived types (constants)
  & YDSOIL , YDCST, &
  ! output prognostics at T 
  & PSSN,  PTSN, PRSN, PWSN, PASN, &
  ! output fluxes 
  & PGSN, PMSN, PEMSSN, & 
  ! Diagnostics 
  & PDHTSS, PDHSSS)



USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOS_SOIL , ONLY : TSOIL 
USE YOS_CST  , ONLY : TCST

USE SRFSN_WEBAL_MOD
USE SRFSN_RSN_MOD
USE SRFSN_ASN_MOD
USE SRFSN_VGRID_MOD 
USE SRFSN_REGRID_MOD
USE SRFSN_SSRABS_MOD

USE ABORT_SURF_MOD

!**** *SRFSN_DRIVER* - Snow scheme driver 
!     PURPOSE.
!     --------
!          THIS ROUTINE CONTROLS THE SNOW SCHEME 

!**   INTERFACE.
!     ----------
!          *SRFSN_DRIVER* IS CALLED FROM *SURFTSTP*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----

!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*       NUMBER OF GRID POINTS PER PACKET
!    *KLEVSN*     NUMBER OF SNOW LAYERS

!     INPUT PARAMETERS (REAL):
!    *PTMST*      TIME STEP                                      S

!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*     LAND/SEA MASK (TRUE/FALSE) 
!    *PSDOR*      OROGRAPHIC PARAMETER                           m

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):
!    *PSSNM1M*    SNOW MASS (per unit area)                    kg m-2
!    *PTSNM1M*    SNOW TEMPERATURE                               K
!    *PASNM1M*    SNOW ALBEDO                                    -
!    *PRSNM1M*    SNOW DENSITY                                 kg m-3
!    *PWSNM1M*    SNOW LIQUID WATER CONTENT                    kg m-2

!    *PRFTI*      TILE FRACTIONS
!    *PTSAM1M*    SOIL TEMPERATURE                               K
!    *PUSRF*      WIND U LOWEST MODEL LEVEL                     m s-1
!    *PVSRF*      WIND V LOWEST MODEL LEVEL                     m s-1
!    *PTSRF*      AIR TEMPERATURE LOWEST MODEL LEVEL             K
!    *PAPRS*      AIR PRESSURE    LOWEST MODEL LEVEL             Pa

!    *PSSFC*      CONVECTIVE SNOWFALL                        kg m-2 s-1
!    *PSSFL*      LARGE-SCALE SNOWFALL                       kg m-2 s-1   
!    *PTSFC*      CONVECTIVE Throughfall                     kg m-2 s-1
!    *PTSFL*      LARGE-SCALE Throughfall                    kg m-2 s-1

!    *PSLRFLTI*   NET LW RADIATION TILED                        W m-2
!    *PSSRFLTI*   NET SW RADIATION TILED                        W m-2  
!    *PAHFSTI*    SENSIBLE HEAT FLUX TILED                      W m-2  
!    *PEVAPTI*    EVAPORATION  TILED                            kg m-2 s-1
!    *PEVAPSNW*   EVAPORATION FROM SNOW UNDER FOREST            kg m-2 s-1

!    *YDSOIL*    SOIL DERIVED TYPE WITH CONSTATNS 
!    *YDCST*     CONSTANTS 

!     OUTPUT PARAMETERS AT T+1 (UNFILTERED,REAL):
!    *PSSN*    SNOW MASS (per unit area)                      kg m-2
!    *PTSN*    SNOW TEMPERATURE                                 K
!    *PASN*    SNOW ALBEDO                                      -
!    *PRSN*    SNOW DENSITY                                   kg m-3
!    *PWSN*    SNOW LIQUID WATER CONTENT                      kg m-2


!     OUTPUT PARAMETERS (DIAGNOSTIC):
!    *PDHTSS*     Diagnostic array for snow T (see module yomcdh)
!    *PDHSSS*     Diagnostic array for snow mass (see module yomcdh)

!     METHOD.
!     -------
!          

!     EXTERNALS.
!     ----------
!          SRFSN_WEBAL - WATER/ENERGY BALANCE
!          SRFSN_RSN   - SNOW DENSITY 
!          SRFSN_ASN   - SNOW ALBEDO 

!     REFERENCE.
!     ----------
!          

!     Modifications:
!     Original   E. Dutra      ECMWF     04/12/2015
!                G. Arduini    ECMWF     01/09/2021

!     ------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments 
INTEGER(KIND=JPIM), INTENT(IN)   :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KLON
INTEGER(KIND=JPIM), INTENT(IN)   :: KLEVSN
REAL(KIND=JPRB)   , INTENT(IN)   :: PTMST
LOGICAL           , INTENT(IN)   :: LDLAND(:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PSDOR(:)

REAL(KIND=JPRB),    INTENT(IN)   :: PSSNM1M(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PTSNM1M(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PRSNM1M(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PWSNM1M(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PASNM1M(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PFRTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PTSAM1M(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PUSRF(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PVSRF(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PTSRF(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PAPRS(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PSSFC(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PSSFL(:)
REAL(KIND=JPRB),    INTENT(INOUT)   :: PTSFC(:)
REAL(KIND=JPRB),    INTENT(INOUT)   :: PTSFL(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PSLRFLTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PSSRFLTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PAHFSTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PEVAPTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PEVAPSNW(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PWSAM1M(:,:)
INTEGER(KIND=JPIM), INTENT(IN)   :: KSOTY(:)

TYPE(TSOIL)       , INTENT(IN)   :: YDSOIL
TYPE(TCST)        , INTENT(IN)   :: YDCST

REAL(KIND=JPRB),    INTENT(OUT)  :: PSSN(:,:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PTSN(:,:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PRSN(:,:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PwSN(:,:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PASN(:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PGSN(:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PMSN(:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PEMSSN(:)

REAL(KIND=JPRB),    INTENT(OUT)  :: PDHTSS(:,:,:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PDHSSS(:,:,:)

! Local variables 
REAL(KIND=JPRB) :: ZFRSN(KLON)   ! snow cover fraction 
REAL(KIND=JPRB) :: ZSNOWF(KLON)  ! total snowfall 
REAL(KIND=JPRB) :: ZRAINF(KLON)  ! total rainfall 
REAL(KIND=JPRB) :: ZHFLUX(KLON)  ! Net heat flux to the snow 
REAL(KIND=JPRB) :: ZEVAPSN(KLON) ! Snow evaporation 
REAL(KIND=JPRB) :: ZSURFCOND(KLON) ! under lying surface thermal conductivity 
REAL(KIND=JPRB) :: ZPHASE(KLON)
REAL(KIND=JPRB) :: ZMELTSN(KLON,KLEVSN) ! Snow melting in each layer
REAL(KIND=JPRB) :: ZFREZSN(KLON,KLEVSN) ! Snow freezing in each layer
REAL(KIND=JPRB) :: ZSSNM1M(KLON,KLEVSN)
REAL(KIND=JPRB) :: ZTSNM1M(KLON,KLEVSN)
REAL(KIND=JPRB) :: ZRSNM1M(KLON,KLEVSN)
REAL(KIND=JPRB) :: ZWSNM1M(KLON,KLEVSN)
REAL(KIND=JPRB) :: ZDSNOUT(KLON,KLEVSN)

INTEGER(KIND=JPIM) :: KLEVSNA(KLON)

! Temporary m1 arrays
REAL(KIND=JPRB) :: ZZSSNM1M(KLON,KLEVSN)
REAL(KIND=JPRB) :: ZZTSNM1M(KLON,KLEVSN)
REAL(KIND=JPRB) :: ZZRSNM1M(KLON,KLEVSN)
REAL(KIND=JPRB) :: ZZWSNM1M(KLON,KLEVSN)

REAL(KIND=JPRB) :: ZSNOTRS(KLON,KLEVSN+1) ! Solar rad abs in each snow layer 
                                           ! (+1 amount in the soil)

REAL(KIND=JPRB) :: ZFF              ! Frozen soil fraction
LOGICAL          :: LLNOSNOW(KLON)  ! FALSE to compute snow 

INTEGER(KIND=JPRB) :: JL,JK

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "fcsurf.h"

!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRFSN_DRIVER_MOD:SRFSN_DRIVER',0,ZHOOK_HANDLE)

ASSOCIATE(RTF1=>YDSOIL%RTF1, RTF2=>YDSOIL%RTF2, RTF3=>YDSOIL%RTF3, RTF4=>YDSOIL%RTF4)
!     ------------------------------------------------------------------
!*         1.1 Global computations 
!*             Snow fraction, total heat and precip/snow to the snow scheme 
!             -----------------------------------------------------------

DO JL=KIDIA,KFDIA
! This safety check must be put for DA 
  DO JK=1,KLEVSN
    ZZSSNM1M(JL,JK)=MAX( 0._JPRB, PSSNM1M(JL,JK) )
    ZZWSNM1M(JL,JK)=MIN(ZZSSNM1M(JL,JK), MAX( 0._JPRB, PWSNM1M(JL,JK) ))
  ENDDO

! SNOW COVER FRACTION 
  ZFRSN(JL)=MAX(PFRTI(JL,5)+PFRTI(JL,7),YDSOIL%RFRTINY)

! POINTS TO COMPUTE WATER&ENERGY BALANCE 
  IF (ZFRSN(JL) < YDSOIL%RFRSMALL) THEN
    LLNOSNOW(JL)=.TRUE.
  ELSE
    LLNOSNOW(JL)=.FALSE.
  ENDIF
  
! NET RAINFALL & SNOWFALL 
! this needs to be checked carefully for points with subgrid-scale lakes: 
! we're not sure if this partition is correct 
  IF (LDLAND(JL)) THEN
    ! Snowfall is all redirected here ! 
    ZSNOWF(JL) = PSSFC(JL) + PSSFL(JL)
  ELSE 
    ZSNOWF(JL) = 0.0_JPRB 
  ENDIF 
  
  ! Evap. is only fractional 
  ZEVAPSN(JL) = PFRTI(JL,5)*PEVAPTI(JL,5) + PFRTI(JL,7)*PEVAPSNW(JL)
  IF (LLNOSNOW(JL)) THEN
    ZHFLUX(JL) = 0.0_JPRB 
    ZRAINF(JL) = 0.0_JPRB
    ZSURFCOND(JL)=0._JPRB
    
  ELSE
!           NET HEAT FLUX AT SNOW SURFACE; EQUATIONS APPLY TO TOTAL SNOW MASS IN THE GRID SQUARE. 
!           HOWEVER THIS FLUX IS **NOT** DIVIDED BY TOTAL SNOW FRACTION, 
!           as the fraction is added to the left hand side of the snow energy
!           budget equation in srfsn_webal (see heat conductivity formulation)

    ZHFLUX(JL)=( PFRTI(JL,5)*(PAHFSTI(JL,5)+YDCST%RLSTT*PEVAPTI(JL,5))&
      & +PFRTI(JL,7)*(PAHFSTI(JL,7)+YDCST%RLVTT*(PEVAPTI(JL,7)-PEVAPSNW(JL))&
      & +YDCST%RLSTT*PEVAPSNW(JL))&
      & +PFRTI(JL,5)*PSSRFLTI(JL,5)&
      & +PFRTI(JL,7)*PSSRFLTI(JL,7)&
      & +PFRTI(JL,5)*PSLRFLTI(JL,5)&
      & +PFRTI(JL,7)*PSLRFLTI(JL,7)&
      & ) !/ZFRSN(JL) Average grid-box flux 
   
    
    ! Rainfall is only fractional 
    ZRAINF(JL) = ZFRSN(JL)*(PTSFC(JL)+PTSFL(JL))
    PTSFC(JL) = (1._JPRB - ZFRSN(JL) )*PTSFC(JL)
    PTSFL(JL) = (1._JPRB - ZFRSN(JL) )*PTSFL(JL)
   
 
! added fix to be consistent with Peters-Lidard et al. 1998 
      IF(PTSAM1M(JL,1) < RTF1.AND.PTSAM1M(JL,1) > RTF2) THEN
        ZFF=0.5_JPRB*(1.0_JPRB-SIN(RTF4*(PTSAM1M(JL,1)-RTF3)))
      ELSEIF (PTSAM1M(JL,1) <= RTF2) THEN
        ZFF=1.0_JPRB
      ELSE
        ZFF=0.0_JPRB
      ENDIF
      !  ZFF=0.0_JPRB
    ZSURFCOND(JL) = MAX(0.19_JPRB,MIN(2._JPRB,FSOILTCOND(PWSAM1M(JL,1),ZFF,KSOTY(JL))))

  ENDIF
  
  

ENDDO


!     ------------------------------------------------------------------
!*         1.1 RESET VERTICAL GRID 
!*             
!             -----------------------------------------------------------
  
CALL SRFSN_VGRID(KIDIA,KFDIA,KLON,KLEVSN,LLNOSNOW,PSDOR,&
                 ZZSSNM1M,PRSNM1M,&
                 YDSOIL%RLEVSNMIN,YDSOIL%RLEVSNMAX,&
                 YDSOIL%RLEVSNMIN_GL,YDSOIL%RLEVSNMAX_GL,&
                 ZDSNOUT,KLEVSNA)

!     ------------------------------------------------------------------
!*         1.1b REGRID FIELDS IF VERTICAL DISCRETIZATION CHANGED
!*             
!             -----------------------------------------------------------

CALL SRFSN_REGRID(KIDIA,KFDIA,KLON,KLEVSN,LLNOSNOW,&
                  ZZSSNM1M,PRSNM1M,ZZWSNM1M,PTSNM1M,ZDSNOUT,&
                  ZSSNM1M,ZRSNM1M,ZWSNM1M,ZTSNM1M)


!     ------------------------------------------------------------------
!*         2. Absorption of solar radiation by snow
ZSNOTRS(KIDIA:KFDIA,1:KLEVSN+1)=0._JPRB
CALL SRFSN_SSRABS(KIDIA,KFDIA,KLON,KLEVSN,&
                  LLNOSNOW,PFRTI,PSSRFLTI,&
                  ZSSNM1M,ZRSNM1M,&
                  YDSOIL,YDCST,&
                  ZSNOTRS)

!     ------------------------------------------------------------------
!*         3. Solve energy / water balance 
!*             
!             -----------------------------------------------------------

CALL SRFSN_WEBAL(KIDIA,KFDIA,KLON,KLEVSN, &
 & PTMST,LLNOSNOW,ZFRSN,&
 & ZSSNM1M,ZWSNM1M,ZRSNM1M,ZTSNM1M,&
 & PTSAM1M(:,1),ZHFLUX,ZSNOTRS,ZSNOWF,ZRAINF,ZEVAPSN,ZSURFCOND,&
 & PAPRS,&
 & YDSOIL,YDCST,&
 & PSSN,PWSN,PTSN,&
 & PGSN,PMSN,ZMELTSN,ZFREZSN,&
 & PDHTSS,PDHSSS)
PEMSSN(KIDIA:KFDIA) = 0._JPRB 
!     ------------------------------------------------------------------
!*         4. Update snow density 
!*             
!             -----------------------------------------------------------

CALL SRFSN_RSN(KIDIA,KFDIA,KLON,KLEVSN,PTMST,LLNOSNOW,&
              &ZRSNM1M,ZSSNM1M,ZTSNM1M,ZWSNM1M,PWSN,&
              &ZSNOWF,PUSRF,PVSRF,PTSRF,&
              &YDSOIL,YDCST,PRSN,PDHTSS)
 
!     ------------------------------------------------------------------
!*         5. Update snow albedo 
!*             
!             -----------------------------------------------------------
ZPHASE(KIDIA:KFDIA)=ZMELTSN(KIDIA:KFDIA,1)-ZFREZSN(KIDIA:KFDIA,1)
CALL SRFSN_ASN(KIDIA,KFDIA,KLON,PTMST,LLNOSNOW,PASNM1M,&
 & ZPHASE,ZTSNM1M,ZSNOWF,YDSOIL,YDCST,PASN)


!    ------------------------------------------------------------------
!*   6. DDH
!*
!    ------------------------------------------------------------------
IF (SIZE(PDHTSS) > 0 .AND. SIZE(PDHSSS) > 0) THEN
  DO JL=KIDIA,KFDIA
    IF (.NOT. LLNOSNOW(JL)) THEN      
      ! ls snowfall
      PDHSSS(JL,1,2)=PSSFL(JL)
      ! conv snowfall  
      PDHSSS(JL,1,3)=PSSFC(JL)
!     ! total rainfall (throughfall)
!     PDHSSS(JL,1,7)=ZRAINF(JL)
! Snow sensible heat flux
      PDHTSS(JL,1,11)=PFRTI(JL,5)*PAHFSTI(JL,5)+PFRTI(JL,7)*PAHFSTI(JL,7)
! Snow latent heat flux
      PDHTSS(JL,1,12)=PFRTI(JL,5)*YDCST%RLSTT*PEVAPTI(JL,5)+&
                    & PFRTI(JL,7)*(YDCST%RLVTT*(PEVAPTI(JL,7)-PEVAPSNW(JL))+&
                    & YDCST%RLSTT*PEVAPSNW(JL))
! Snow density
      DO JK=1,KLEVSN 
        PDHTSS(JL,JK,5)=PRSNM1M(JL,JK)
      ENDDO
      IF (KLEVSN > 1) THEN
        PDHSSS(JL,2:KLEVSN,2)=0._JPRB
        PDHSSS(JL,2:KLEVSN,3)=0._JPRB
        PDHTSS(JL,2:KLEVSN,11)=0._JPRB
        PDHTSS(JL,2:KLEVSN,12)=0._JPRB
      ENDIF
    ENDIF
  ENDDO
ENDIF 

!    ------------------------------------------------------------------
!*   7. Budget computations : TESTING ONLY !! 
!*   Should be comment out when running coupled to IFS ! 
!    ------------------------------------------------------------------

IF (YDSOIL%LESNCHECK) THEN
  CALL BUDGET_MASS_DDH
  CALL BUDGET_ENERGY_DDH
ENDIF 
END ASSOCIATE
!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRFSN_DRIVER_MOD:SRFSN_DRIVER',1,ZHOOK_HANDLE)

CONTAINS 

SUBROUTINE BUDGET_ENERGY_DDH 
IMPLICIT NONE
INTEGER(KIND=JPRB) :: JL,JK
REAL(KIND=JPRB)    :: ZFLUX,ZSTORAGE,ZRES
REAL(KIND=JPRB)    :: ZWTRH,ZEPSILON,ZTMP

ZEPSILON=10._JPRB*EPSILON(ZEPSILON)

ZWTRH=MAX(ZEPSILON,1.e-8_JPRB)   ! snowmip 3 wm-2
ZWTRH=MAX(ZEPSILON,1.e-6_JPRB)  ! added for glaciers   
DO JL=KIDIA,KFDIA
  IF (.NOT. LLNOSNOW(JL)) THEN    
    ZSTORAGE=PDHTSS(JL,1,15)/PTMST
    ZFLUX=ZHFLUX(JL)-PDHTSS(JL,1,13)-PDHTSS(JL,1,14)
    ZRES=ZSTORAGE-ZFLUX
    IF ( ABS(ZRES) > ZWTRH ) THEN
       write(*,'("DDH SEBAL:storage,flux,residual (kg/m2/s):",1X,I4,3(1X,E14.6E3))') &
            JL,ZSTORAGE,ZFLUX,ZRES
       write(*,'("DDH SEBAL1: snow mass",1X,20(1X,E14.6E3))') (PSSN(JL,JK),JK=1,KLEVSN)
       write(*,'("DDH SEBAL1: snowT t-1",1X,20(1X,E14.6E3))') (ZTSNM1M(JL,JK),JK=1,KLEVSN)
       write(*,'("DDH SEBAL1: snowT t+1",1X,20(1X,E14.6E3))') (PTSNM1M(JL,JK),JK=1,KLEVSN)
       write(*,'("DDH SEBAL1: - basal flux",1X,1(1X,E14.6E3))') -PDHTSS(JL,1,13)
       write(*,'("DDH SEBAL1: - phase changes",1X,1(1X,E14.6E3))') -PDHTSS(JL,1,14)
       write(*,'("DDH SEBAL1: ZHFLUX",1X,1(1X,E14.6E3))') ZHFLUX(JL)
       ZTMP=PFRTI(JL,5)*PAHFSTI(JL,5)+PFRTI(JL,7)*PAHFSTI(JL,7)
       write(*,'("DDH SEBAL1: Qh",1X,1(1X,E14.6E3))') ZTMP
       ZTMP=PFRTI(JL,5)*YDCST%RLSTT*PEVAPTI(JL,5)+PFRTI(JL,7)*(YDCST%RLVTT*(PEVAPTI(JL,7)-PEVAPSNW(JL))+YDCST%RLSTT*PEVAPSNW(JL))
       write(*,'("DDH SEBAL1: Qle",1X,1(1X,E14.6E3))') ZTMP
       ZTMP=PFRTI(JL,5)*PSSRFLTI(JL,5)+PFRTI(JL,7)*PSSRFLTI(JL,7)
       write(*,'("DDH SEBAL1: SWd,SWu,SWN",1X,1(1X,E14.6E3))') ZTMP
       ZTMP=PFRTI(JL,5)*PSLRFLTI(JL,5)+PFRTI(JL,7)*PSLRFLTI(JL,7)
       write(*,'("DDH SEBAL1: LWd,LWu,LWN",1X,1(1X,E14.6E3))') ZTMP
       IF (YDSOIL%LESNCHECKAbort) THEN
         CALL ABORT_SURF('SNOW ENERGY BALANCE CLOSURE ERROR')
       ENDIF
    ENDIF  
  ENDIF
  
ENDDO
END SUBROUTINE BUDGET_ENERGY_DDH 


SUBROUTINE BUDGET_MASS_DDH
IMPLICIT NONE 
INTEGER(KIND=JPRB) :: JL
REAL(KIND=JPRB)    :: ZFLUX,ZSTORAGE,ZRES
REAL(KIND=JPRB)    :: ZWTRH,ZEPSILON

ZEPSILON=10._JPRB*EPSILON(ZEPSILON)

ZWTRH=MAX(ZEPSILON,1.e-16_JPRB)  ! snow mip 1e-7 ! 
DO JL=KIDIA,KFDIA
  IF (.NOT. LLNOSNOW(JL)) THEN    
    ZSTORAGE=SUM(PSSN(JL,:)-ZSSNM1M(JL,:))/PTMST
    ZFLUX=ZSNOWF(JL)+ZRAINF(JL)+ZEVAPSN(JL)-PMSN(JL)
    ZRES=ZSTORAGE-ZFLUX
    ! add sum of snow mass to account for huge snow mass point that faill to pass the test 
    ZWTRH=MAX(1._JPRB,SUM(PSSN(JL,:)))*MAX(ZEPSILON,1.e-16_JPRB)
    IF ( ABS(ZRES) > ZWTRH ) THEN
      write(*,'("DDH WBAL:",1X,I4,3(1X,E14.6E3))') &
            JL,ZSTORAGE,ZFLUX,ZRES
      write(*,'("DDH WBAL1:",1X,I4,6(1X,E14.6E3))') &
            JL,ZSSNM1M(JL,:),PSSN(JL,:),ZSNOWF(JL),ZRAINF(JL),ZEVAPSN(JL),-PMSN(JL)
      IF (YDSOIL%LESNCHECKAbort) THEN
        CALL ABORT_SURF('SNOW WATER BALANCE CLOSURE ERROR')
      ENDIF
    ENDIF
   
!     write(*,'("DH WBAL:",1X,I4,3(1X,E14.6E3))') &
!             JL,ZSTORAGE,ZFLUX,ZRES
!       write(*,'("DH WBAL1:",1X,I4,7(1X,E14.6E3))') &
!             JL,PSSNM1M(JL,1),PSSN(JL,1),ZFRSN(JL),ZSNOWF(JL),ZRAINF(JL),ZEVAPSN(JL),-PMSN(JL)
    
  ENDIF
ENDDO

END SUBROUTINE BUDGET_MASS_DDH


END SUBROUTINE SRFSN_DRIVER


END MODULE SRFSN_DRIVER_MOD 


 
