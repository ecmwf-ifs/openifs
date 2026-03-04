! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUBASENS &
 & (YDEPHLI, YDECLDP,  YDECUMF2,&
 & KIDIA,    KFDIA,    KLON,     KLEV,   KINDEX,  LDMIXS,&
 & PTENH,    PQENH,    PGEOH,    PAPH,&
 & PQHFL,    PAHFS,    PKMFL,&
 & PTEN,     PQEN,     PQSEN,    PGEO,&
 & PTU,      PQU,      PLU,      PWU2H,  PWUBASE,&
 & KLAB,     LDCUM,    LDSC,     KCBOT,  KBOTSC,&
 & KCTOP,    KDPL,     PCAPE )  

!          THIS ROUTINE CALCULATES CLOUD BASE FIELDS
!          CLOUD BASE HEIGHT AND CLOUD TOP HEIGHT

!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD BASE AND CLOUD TOP VALUES FOR CU-PARAMETRIZATION

!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
!          IT RETURNS CLOUD FIELDS VALUES AND FLAGS AS FOLLOWS;
!                 KLAB=0 FOR STABLE LAYERS
!                 KLAB=1 FOR SUBCLOUD LEVELS
!                 KLAB=2 FOR CLOUD LEVELS LEVEL

!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD TOP
!          (ENTRAINING PLUME, WITH ENTRAINMENT PROPORTIONAL TO (1/Z))

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KINDEX*       TOP INDEX FOR OUTER VERTICAL LOOP

!    INPUT PARAMETERS (LOGICAL):
!    *LDMIXS*        WEAK (FALSE) OR STRONG (TRUE) CLOUD MIXING FOR SURFACE PARCEL ONLY

!    INPUT PARAMETERS (REAL):

! not used at the moment because we want to use linear intepolation
! for fields on the half levels.

!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG

!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
!    *PKMFL*        SURFACE KINEMATIC MOMENTUM FLUX              M2/S2  

!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PQSEN*        PROVISIONAL ENVIRONMENT SATU. HUMIDITY (T+1)  KG/KG
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2

!    UPDATED PARAMETERS (REAL):

!    *PTU*          TEMPERATURE IN UPDRAFTS                         K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *PWU2H*        KINETIC ENERGY IN UPDRAFTS  SURFACE PARCEL    M2/S2

!    UPDATED PARAMETERS (INTEGER):

!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!                        KLAB=2 FOR CLOUD LEVELS

!    OUTPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDSC*         FLAG: .TRUE. IF BL-CLOUDS EXIST

!    OUTPUT PARAMETERS (INTEGER):

!    *KCBOT*       CLOUD BASE LEVEL !    
!    *KCTOP*       CLOUD TOP LEVEL = HEIGHEST HALF LEVEL 
!                  WITH A NON-ZERO CLOUD UPDRAFT.
!    *KBOTSC*      CLOUD BASE LEVEL OF BL-CLOUDS
!    *KDPL*        DEPARTURE LEVEL
!    *PCAPE*       PSEUDOADIABATIQUE max CAPE (J/KG)

!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT

!     AUTHOR.
!     -------
!      A. Pier Siebesma   KNMI ********      

!     MODIFICATIONS.
!     --------------
!      modified C Jakob (ECMWF) (01/2001) 
!      modified P Bechtold (ECMWF) (08/2002) 
!      02-11-02 : Use fixed last possible departure level and 
!                 last updraft computation level for bit-reproducibility
!                 D.Salmond &  J. Hague
!      03-07-03 : Tuning for p690     J. Hague
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      N.Semane+P.Bechtold    04-10-2012 Add RPLRG/RPLDARE factors for small planet
!      Simplified version to be used in vertical diffusion to estimate PBL height.
!----------------------------------------------------------------------

USE YOEPHLI   , ONLY : TEPHLI
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST    , ONLY : RCPD, RETV, RD, RG, RLVTT, RLSTT, RTT  
USE PARPHY    , ONLY : RKAP
USE YOECUMF2  , ONLY : TECUMF2
USE YOECLDP   , ONLY : TECLDP
USE YOETHF    , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 &                     R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RALFDCP, RTWAT, RTICE, RTICECU, &
 &                     RTWAT_RTICECU_R, RTWAT_RTICE_R  
USE YOMDYNCORE, ONLY : RPLRG, RPLDARE

IMPLICIT NONE

TYPE(TECLDP)      ,INTENT(IN)    :: YDECLDP
TYPE(TECUMF2)     ,INTENT(IN)    :: YDECUMF2
TYPE(TEPHLI)      ,INTENT(IN)    :: YDEPHLI
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KINDEX
LOGICAL           ,INTENT(IN)    :: LDMIXS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQENH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQHFL(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAHFS(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKMFL(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEO(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWU2H(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWUBASE(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KLAB(KLON,KLEV) 
LOGICAL           ,INTENT(INOUT) :: LDCUM(KLON) 
LOGICAL           ,INTENT(OUT)   :: LDSC(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCBOT(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KBOTSC(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KCTOP(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KDPL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAPE(KLON) 
INTEGER(KIND=JPIM) ::  ICTOP(KLON),            ICBOT(KLON),&
 & IBOTSC(KLON),           ILAB(KLON,KLEV),&
 & IDPL(KLON)  

LOGICAL ::         LL_LDBASE(KLON),&
 & LLGO_ON(KLON),&
 & LLDEEP(KLON),    LLDCUM(KLON), &
 & LLDSC(KLON),     LLFIRST(KLON)  
LOGICAL ::     LLRESET,        LLRESETJL(KLON)

INTEGER(KIND=JPIM) :: ICALL, IK, IS, JK, JL, JKK, JKT1, JKT2, JKT, JKB

REAL(KIND=JPRB)    :: ZS(KLON,KLEV), ZSENH(KLON,KLEV+1), ZQENH(KLON,KLEV+1), ZSUH (KLON,KLEV),&
 & ZBUOH(KLON,KLEV), ZWU2H(KLON,KLEV) 
REAL(KIND=JPRB) :: ZQOLD(KLON),ZPH(KLON)
REAL(KIND=JPRB) :: ZMIX(KLON)
REAL(KIND=JPRB) :: ZDZ(KLON), ZCBASE(KLON)

REAL(KIND=JPRB) :: ZLU(KLON,KLEV), ZQU(KLON,KLEV), ZTU(KLON,KLEV)

REAL(KIND=JPRB) :: ZCAPE(KLON,KLEV) ! local for CAPE at every departure level

REAL(KIND=JPRB) :: ZBUOF     ! BUOYANCY
REAL(KIND=JPRB) :: ZRHO      ! DENSITY AT SURFACE (KG/M^3) 
REAL(KIND=JPRB) :: ZKHVFL    ! SURFACE BUOYANCY FLUX (K M/S)
REAL(KIND=JPRB) :: ZWS       ! SIGMA_W AT LOWEST MODEL HALFLEVEL (M/S)
REAL(KIND=JPRB) :: ZUST      ! U* 
REAL(KIND=JPRB) :: ZQEX(KLON), ZQEXC ! HUMIDITY EXCESS AT LOWEST MODEL HALFLEVEL (KG/KG)
REAL(KIND=JPRB) :: ZTEX(KLON), ZTEXC ! TEMPERATURE EXCESS AT LOWEST MODEL HALFLEVEL (K)
REAL(KIND=JPRB) :: ZEPS      ! FRACTIONAL ENTRAINMENT RATE   [M^-1]
REAL(KIND=JPRB) :: ZTVENH    ! ENVIRONMENT VIRTUAL TEMPERATURE AT HALF LEVELS (K)  
REAL(KIND=JPRB) :: ZTVUH     ! UPDRAFT VIRTUAL TEMPERATURE AT HALF LEVELS     (K)
REAL(KIND=JPRB) :: ZLGLAC    ! UPDRAFT LIQUID WATER FROZEN IN ONE LAYER
REAL(KIND=JPRB) :: ZQSU, ZCOR, ZDQ, ZALFAW, ZFACW, ZFACI, ZFAC,&
 & ZESDP, ZDQSDT, ZDTDP, ZDP,ZPDIFFTOP,ZPDIFFBOT,ZSF,ZQF,ZAW,ZBW  
REAL(KIND=JPRB) :: ZWORK1, ZWORK2! work arrays for T and w perturbations
REAL(KIND=JPRB) :: ZRCPD, ZRG, ZTMP
REAL(KIND=JPRB) :: ZXENTRORG

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "cuadjtq.intfb.h"
#include "fcttre.func.h"

!----------------------------------------------------------------------
!     0.           INITIALIZE CONSTANTS AND FIELDS
!                  -------------------------------
!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUBASENS',0,ZHOOK_HANDLE)
ASSOCIATE(RLMIN=>YDECLDP%RLMIN, &
 & ENTRORG2=>YDECUMF2%ENTRORG2, ENTSTPC12=>YDECUMF2%ENTSTPC12, ENTSTPC22=>YDECUMF2%ENTSTPC22, &
 & ENTSTPC32=>YDECUMF2%ENTSTPC32, NJKT22=>YDECUMF2%NJKT22, RDEPTHS2=>YDECUMF2%RDEPTHS2)
ZAW    = 1.0_JPRB
ZBW    = 1.0_JPRB

DO JL=KIDIA,KFDIA
  PWUBASE(JL)=0.0_JPRB
  LLGO_ON(JL)=.TRUE.
  LLFIRST(JL)=.TRUE.
  KDPL(JL)=KLEV
ENDDO

JKT1=KINDEX
JKT2=NJKT22
ZRG=1.0_JPRB/RG
ZRCPD=1.0_JPRB/RCPD

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZTU(JL,JK) = PTU(JL,JK)
    ZQU(JL,JK) = PQU(JL,JK)
    ZLU(JL,JK) = PLU(JL,JK)
    ILAB(JL,JK)= KLAB(JL,JK)
    ZCAPE(JL,JK)= 0.0_JPRB
  ENDDO
ENDDO

!----------------------------------------------------------------------
!       -----------------------------------------------------------
!       1.1  PREPARE FIELDS ON HALF LEVELS BY LINEAR INTERPOLATION
!             OF SPECIFIC HUMIDITY AND STATIC ENERGY
!       -----------------------------------------------------------

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PWU2H(JL,JK)=0.0_JPRB
    ZWU2H(JL,JK)=0.0_JPRB
    ZS   (JL,JK) = RCPD*PTEN(JL,JK) + PGEO(JL,JK)
    ZQENH(JL,JK) = PQENH(JL,JK)
    ZSENH(JL,JK) = RCPD*PTENH(JL,JK)+PGEOH(JL,JK)
  ENDDO
ENDDO

DO JKK=KLEV,JKT1,-1 ! Big external loop for level testing:
                    ! find first departure level that produces deepest cloud top
                    ! or take surface level for shallow convection and Sc
   !
   !        ---------------------------------------------------------
   !        1.2    INITIALISE FIELDS AT DEPARTURE HALF MODEL LEVEL
   !        ---------------------------------------------------------
   !
  IS=0
  DO JL=KIDIA,KFDIA
    IF (LLGO_ON(JL)) THEN
      IS=IS+1
      IDPL(JL)    =JKK      ! departure level
      ICBOT  (JL) =JKK      ! cloud base level for convection, (-1 if not found)
      IBOTSC (JL) =KLEV-1   ! sc    base level for sc-clouds , (-1 if not found)
      ICTOP(JL)   =KLEV-1   ! cloud top for convection (-1 if not found)
      LLDCUM(JL)  =.FALSE.  ! on exit: true if cloudbase=found
      LLDSC (JL)  =.FALSE.  ! on exit: true if cloudbase=found
      LL_LDBASE(JL)   =.FALSE. ! on exit: true if cloudbase=found
    ENDIF 
  ENDDO

  IF(IS /= 0) THEN

    IF(JKK == KLEV) THEN

      ZTEXC=0.2_JPRB
      ZQEXC=1.E-4_JPRB
      DO JL=KIDIA,KFDIA
        IF (LLGO_ON(JL)) THEN
          ZRHO  = PAPH(JL,JKK+1)/(RD*(PTEN(JL,JKK)*(1.0_JPRB+RETV*PQEN(JL,JKK))))
          ZKHVFL= (PAHFS(JL,JKK+1)*ZRCPD+RETV*PTEN(JL,JKK)*PQHFL(JL,JKK+1))/(ZRHO*RPLRG*RPLDARE)
          ZUST  = MAX( SQRT(PKMFL(JL)), 0.1_JPRB )
          ZWS   = ZUST**3._JPRB - 1.5_JPRB*RKAP*ZKHVFL*(PGEOH(JL,KLEV)-PGEOH(JL,KLEV+1))/PTEN(JL,KLEV)
          ZTEX(JL)= 0.0_JPRB
          ZQEX(JL)= 0.0_JPRB
          IF( ZKHVFL < 0.0_JPRB ) THEN
            ZWS =1.2_JPRB*ZWS**.3333_JPRB
            ILAB(JL,JKK)= 1
            ZTEX(JL)   = MAX(-1.5_JPRB*PAHFS(JL,JKK+1)/(ZRHO*ZWS*RCPD*RPLRG*RPLDARE),ZTEXC)
            ZQEX(JL)   = MAX(-1.5_JPRB*PQHFL(JL,JKK+1)/(ZRHO*ZWS*RPLRG*RPLDARE),ZQEXC)
            ZQU (JL,JKK) = ZQENH(JL,JKK) + ZQEX(JL)
            ZSUH (JL,JKK)= ZSENH(JL,JKK) + RCPD*ZTEX(JL)
            ZTU (JL,JKK) = (ZSENH(JL,JKK)-PGEOH(JL,JKK))*ZRCPD + ZTEX(JL)
            ZLU (JL,JKK) = 0.0_JPRB
            ZWU2H(JL,JKK) = ZWS**2+0.1_JPRB
            PWU2H(JL,JKK) = ZWU2H(JL,JKK)
        !
        !  determine buoyancy at lowest half level
        !
            ZTVENH            = (1.0_JPRB+RETV*ZQENH(JL,JKK)) &
             & *(ZSENH(JL,JKK)-PGEOH(JL,JKK))*ZRCPD  
            ZTVUH             = (1.0_JPRB+RETV*ZQU(JL,JKK))*ZTU(JL,JKK)
            ZBUOH(JL,JKK) = (ZTVUH-ZTVENH)*RG/ZTVENH
          ELSE
            LLGO_ON(JL)=.FALSE.      ! non-convective point
          ENDIF
        ENDIF
      ENDDO
   
    ELSE

      DO JL=KIDIA,KFDIA
        IF (LLGO_ON(JL)) THEN
          ZRHO  = PAPH(JL,JKK+1)/(RD*(PTEN(JL,JKK)*(1.+RETV*PQEN(JL,JKK))))
          ILAB(JL,JKK)= 1
          ZTEXC=0.2_JPRB
          ZQEXC=1.E-4_JPRB
          IF(JKK==KLEV-1) THEN
            ZTEXC=MAX(ZTEXC,ZTEX(JL))
            ZTEXC=MIN(ZTEXC, 1.0_JPRB)
            ZQEXC=MAX(ZQEXC,ZQEX(JL))
            ZQEXC=MIN(ZQEXC, 5.E-4_JPRB)
          ENDIF
          ZQU (JL,JKK) = ZQENH(JL,JKK) + ZQEXC
          ZSUH (JL,JKK) = ZSENH(JL,JKK) + RCPD*ZTEXC
          ZTU (JL,JKK) = (ZSENH(JL,JKK)-PGEOH(JL,JKK))*ZRCPD + ZTEXC
          ZLU (JL,JKK) = 0.0_JPRB
         ! construct mixed layer for parcels emanating in lowest 60 hPa
          IF (PAPH(JL,KLEV+1)-PAPH(JL,JKK-1)<60.E2_JPRB) THEN
            ZQU(JL,JKK) =0.0_JPRB
            ZSUH(JL,JKK)=0.0_JPRB
            ZWORK1      =0.0_JPRB
            DO JK=JKK+1,JKK-1,-1
              IF( ZWORK1 < 50.E2_JPRB ) THEN
                ZWORK2=PAPH(JL,JK)-PAPH(JL,JK-1)
                ZWORK1      =ZWORK1+ZWORK2
                ZQU(JL,JKK) =ZQU(JL,JKK) +ZQENH(JL,JK)*ZWORK2
                ZSUH(JL,JKK)=ZSUH(JL,JKK)+ZSENH(JL,JK)*ZWORK2
              ENDIF
            ENDDO
            ZQU(JL,JKK) =ZQU(JL,JKK) /ZWORK1+ZQEXC
            ZSUH(JL,JKK)=ZSUH(JL,JKK)/ZWORK1+RCPD*ZTEXC
            ZTU(JL,JKK) =(ZSUH(JL,JKK)-PGEOH(JL,JKK))*ZRCPD+ZTEXC
          ENDIF
          ZWU2H(JL,JKK) = 1.0_JPRB
        ! PWU2H(JL,JKK) = ZWU2H(JL,JKK)

      !
      !  determine buoyancy at lowest half level
      !
          ZTVENH            = (1.0_JPRB+RETV*ZQENH(JL,JKK)) &
           & *(ZSENH(JL,JKK)-PGEOH(JL,JKK))*ZRCPD  
          ZTVUH             = (1.0_JPRB+RETV*ZQU(JL,JKK))*ZTU(JL,JKK)
          ZBUOH(JL,JKK) = (ZTVUH-ZTVENH)*RG/ZTVENH
        ENDIF
      ENDDO
   
    ENDIF

  ENDIF
   
   !----------------------------------------------------------------------
   !     2.0          DO ASCENT IN SUBCLOUD AND LAYER,
   !                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
   !                  ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
   !                  CHECK FOR BUOYANCY AND SET FLAGS
   !                  -------------------------------------
   !       ------------------------------------------------------------
   !        1.2  DO THE VERTICAL ASCENT UNTIL VELOCITY BECOMES NEGATIVE
   !       ------------------------------------------------------------
  DO JK=JKK-1,JKT2,-1
    IS=0

    IF(JKK==KLEV.OR.KINDEX==KLEV-1) THEN ! 1/z mixing for shallow

      DO JL=KIDIA,KFDIA
        IF (LLGO_ON(JL)) THEN
          IS         = IS+1
          ZDZ(JL)    = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG
          ZEPS       = ENTSTPC12/((PGEO(JL,JK)-PGEOH(JL,KLEV+1))*ZRG*RPLRG) + ENTSTPC22
          IF(LDMIXS.AND.KINDEX==KLEV.AND.ZLU(JL,JK+1)>0.0_JPRB) ZEPS=ZEPS*ENTSTPC32
          ZMIX(JL)   = 0.5_JPRB*ZDZ(JL)*RPLRG*ZEPS
          ZMIX(JL)   = MIN(1.0_JPRB, ZMIX(JL))
          ZQF = (PQENH(JL,JK+1) + PQENH(JL,JK))*0.5_JPRB
          ZSF = (ZSENH(JL,JK+1) + ZSENH(JL,JK))*0.5_JPRB
          ZTMP = 1.0_JPRB/(1.0_JPRB+ZMIX(JL))
          ZQU(JL,JK)= (ZQU(JL,JK+1)*(1.0_JPRB-ZMIX(JL))&
         & +2.0_JPRB*ZMIX(JL)*ZQF) * ZTMP  
          ZSUH (JL,JK)= (ZSUH(JL,JK+1)*(1.0_JPRB-ZMIX(JL))&
         & +2.0_JPRB*ZMIX(JL)*ZSF) * ZTMP  
          ZQOLD(JL)  = ZQU(JL,JK)
          ZTU (JL,JK) = (ZSUH(JL,JK)-PGEOH(JL,JK))*ZRCPD
          ZPH  (JL)    = PAPH(JL,JK)
        ENDIF
      ENDDO

    ELSE

      ZXENTRORG=ENTRORG2
      DO JL=KIDIA,KFDIA
        IF (LLGO_ON(JL)) THEN
          ZDZ(JL)    = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG
          ZMIX(JL)=0.4_JPRB*ZXENTRORG*ZDZ(JL)*MIN(1.0_JPRB,(PQSEN(JL,JK)/PQSEN(JL,KLEV))**3)
        ENDIF
      ENDDO

      DO JL=KIDIA,KFDIA
        IF (LLGO_ON(JL)) THEN
          IS         = IS+1
          ZMIX(JL)=MIN(1.0_JPRB,ZMIX(JL))
          ZQF = (PQENH(JL,JK+1) + PQENH(JL,JK))*0.5_JPRB
          ZSF = (ZSENH(JL,JK+1) + ZSENH(JL,JK))*0.5_JPRB
          ZQU(JL,JK)= ZQU(JL,JK+1)*(1.0_JPRB-ZMIX(JL))+ ZQF*ZMIX(JL)
          ZSUH(JL,JK)= ZSUH(JL,JK+1)*(1.0_JPRB-ZMIX(JL))+ ZSF*ZMIX(JL)
          ZQOLD(JL)  = ZQU(JL,JK)
          ZTU (JL,JK)= (ZSUH(JL,JK)-PGEOH(JL,JK))*ZRCPD
          ZPH  (JL)  = PAPH(JL,JK)
        ENDIF
      ENDDO

    ENDIF

    IF (IS == 0) EXIT
     
    IK=JK
    ICALL=1
     
    CALL CUADJTQ &
     & ( YDEPHLI,  KIDIA,    KFDIA,    KLON,    KLEV,      IK,&
     &   ZPH,      ZTU,      ZQU,     LLGO_ON,   ICALL)  
   
   !DIR$ IVDEP
   !OCL NOVREC
   
    DO JL=KIDIA,KFDIA
      IF(LLGO_ON(JL)) THEN
   
   ! add condensation to water
   
        ZDQ=MAX(ZQOLD(JL)-ZQU(JL,JK),0.0_JPRB)
        ZLU(JL,JK)=ZLU(JL,JK+1)+ZDQ

   ! freezing
   
        ZLGLAC=ZDQ*((1.0_JPRB-FOEALFCU(ZTU(JL,JK)))-&
         & (1.0_JPRB-FOEALFCU(ZTU(JL,JK+1))))  
              
   
   ! pseudo-microphysics
   
        ZLU(JL,JK)=0.5_JPRB*ZLU(JL,JK) 
   
   ! update dry static energy after condensation + freezing
   
        ZTU(JL,JK)     = ZTU(JL,JK)+RALFDCP*ZLGLAC
        ZSUH(JL,JK)    = RCPD*ZTU(JL,JK)+PGEOH(JL,JK)
         
   ! Buoyancy on half and full levels
            
        ZTVUH           = (1.0_JPRB+RETV*ZQU(JL,JK)-ZLU(JL,JK))*ZTU(JL,JK)&
         & +RALFDCP*ZLGLAC  
        ZTVENH          = (1.0_JPRB+RETV*ZQENH(JL,JK)) &
         & *(ZSENH(JL,JK)-PGEOH(JL,JK))*ZRCPD  
        ZBUOH(JL,JK)   = (ZTVUH-ZTVENH)*RG/ZTVENH

        ZBUOF          = (ZBUOH(JL,JK) + ZBUOH(JL,JK+1))*0.5_JPRB
   
   ! solve kinetic energy equation
   
        ZTMP=1.0_JPRB/(1.0_JPRB+2.0_JPRB*ZBW*ZMIX(JL))
        ZWU2H(JL,JK) = (ZWU2H(JL,JK+1)*(1.0_JPRB-2.0_JPRB*ZBW*ZMIX(JL))&
         & +2.0_JPRB*ZAW*ZBUOF*ZDZ(JL)) * ZTMP  

        IF(JKK==KLEV) THEN
           PWU2H(JL,JK) = ZWU2H(JL,JK) !save surface parcel KE
        ENDIF
   
   ! compute CAPE for diagnostics
   
        ZCAPE(JL,JKK)  = ZCAPE(JL,JKK) + MAX(0.0_JPRB,ZBUOF*ZDZ(JL))
   
   ! first layer with liquid water - find exact cloud base
   
        IF(ZLU(JL,JK) >0.0_JPRB.AND.ILAB(JL,JK+1)==1) THEN
           
          IK=JK+1
          ZQSU=FOEEWM(ZTU(JL,IK))/PAPH(JL,IK)
          ZESDP=ZQSU
          ZQSU=MIN(0.5_JPRB,ZQSU)
          ZCOR=1.0_JPRB/(1.0_JPRB-RETV  *ZQSU)
          ZQSU=ZQSU*ZCOR
          ZDQ=MIN(0._JPRB,ZQU(JL,IK)-ZQSU)
          ZALFAW=FOEALFA(ZTU(JL,IK))
          ZFACW=R5LES/((ZTU(JL,IK)-R4LES)**2)
          ZFACI=R5IES/((ZTU(JL,IK)-R4IES)**2)
          ZFAC=ZALFAW*ZFACW+(1.-ZALFAW)*ZFACI
          ZCOR=1.0_JPRB/(1.0_JPRB-RETV*ZESDP)
          ZDQSDT=ZFAC*ZCOR*ZQSU
          ZDTDP=RD*ZTU(JL,IK)/(RCPD*PAPH(JL,IK))
          ZDP=ZDQ/(ZDQSDT*ZDTDP)
          ZCBASE(JL)=PAPH(JL,IK)+ZDP
           
   ! chose nearest half level as cloud base
   
          ZPDIFFTOP=ZCBASE(JL)-PAPH(JL,JK)
          ZPDIFFBOT=PAPH(JL,JK+1)-ZCBASE(JL)
           
          IF(ZPDIFFTOP > ZPDIFFBOT.AND.ZWU2H(JL,JK+1)>0.0_JPRB) THEN
            JKB=MIN(KLEV-1,JK+1)
            ILAB(JL,JKB)=2 
            ILAB(JL,JK)=2
            LL_LDBASE(JL) =.TRUE.
            LLDSC(JL)   =.TRUE.
            IBOTSC(JL) =JKB
            ICBOT(JL)  =JKB
            ZLU(JL,JK+1) = RLMIN
          ELSEIF(ZPDIFFTOP <= ZPDIFFBOT.AND.ZWU2H(JL,JK)>0.0_JPRB) THEN
            ILAB(JL,JK)=2
            LL_LDBASE(JL) =.TRUE.
            LLDSC(JL)   =.TRUE.
            IBOTSC(JL) =JK
            ICBOT(JL)  =JK
          ENDIF
   
        ENDIF
   
   ! decide on presence of convection, cloud base and cloud top based on
   ! kinetic energy
   
        IF (ZWU2H(JL,JK) < 0.0_JPRB) THEN
          LLGO_ON(JL) = .FALSE.             
          IF (ZLU(JL,JK+1)>0.0_JPRB) THEN
            ICTOP(JL)   = JK
            LLDCUM(JL)   = .TRUE.
          ELSE
            LLDCUM(JL)   = .FALSE.
          ENDIF
        ELSE
          IF (ZLU(JL,JK)>0.0_JPRB) THEN
            ILAB(JL,JK) = 2
          ELSE
            ILAB(JL,JK) = 1
          ENDIF
        ENDIF
      ENDIF
    ENDDO
   
!     IF (IS == 0) EXIT
  ENDDO
   
  IF( JKK==KLEV.OR. (JKK==KLEV-1.AND.KINDEX==KLEV-1) ) THEN

      ! set values for departure level for PBL clouds = first model level
    DO JL=KIDIA,KFDIA
      LDSC(JL)  = LLDSC(JL)
      IF(LDSC(JL)) THEN
        KBOTSC(JL)= IBOTSC(JL)
      ELSE
        KBOTSC(JL)=-1
      ENDIF
    
      JKT=ICTOP(JL)
      JKB=ICBOT(JL)
      LLDEEP(JL)=PAPH(JL,JKB)-PAPH(JL,JKT)>RDEPTHS2
      IF(LLDEEP(JL).AND.KINDEX<KLEV-1) LLDCUM(JL)=.FALSE. ! no deep allowed for KLEV
      LLDEEP(JL)=.FALSE. ! for deep convection start only at level KLEV-1
                         ! and form mixed layer, so go on
      ! test further for deep convective columns as not yet found
    ! IF ( LLDEEP(JL) ) LLFIRST(JL)=.FALSE.
      LLGO_ON(JL) = .NOT.LLDEEP(JL)
      IF(KINDEX==KLEV-1) LLGO_ON(JL) = .NOT.LLDCUM(JL)
      IF(LLDCUM(JL)) THEN
        KCBOT(JL)= ICBOT(JL)
        KCTOP(JL)= ICTOP(JL)
        KDPL(JL)  = IDPL(JL)
        LDCUM(JL) = LLDCUM(JL)
        PWUBASE(JL)=SQRT(MAX(ZWU2H(JL,JKB),0.0_JPRB))
      ELSE
        KCTOP(JL)=-1
        KCBOT(JL)=-1
        KDPL(JL) =KLEV-1
        LDCUM(JL)=.FALSE.
        PWUBASE(JL)=0.0_JPRB
      ENDIF
    ENDDO
    DO JK=KLEV,1,-1
      DO JL=KIDIA,KFDIA
        JKT=ICTOP(JL)
        IF ( JK>=JKT.OR. (KINDEX>=KLEV-1.AND.JKK==KLEV) ) THEN
          KLAB(JL,JK)=ILAB(JL,JK)
          PTU(JL,JK)=ZTU(JL,JK)
          PQU(JL,JK)=ZQU(JL,JK)
          PLU(JL,JK)=ZLU(JL,JK)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
   
  IF( JKK < KLEV ) THEN
    LLRESET=.FALSE.
    DO JL=KIDIA,KFDIA
      IF ( .NOT.LLDEEP(JL) ) THEN
        JKT=ICTOP(JL)
        JKB=ICBOT(JL)
           ! test on cloud thickness and buoyancy
        LLDEEP(JL)=PAPH(JL,JKB)-PAPH(JL,JKT)>=RDEPTHS2 
      ENDIF
      LLRESETJL(JL)=LLDEEP(JL).AND.LLFIRST(JL)
      LLRESET=LLRESET.OR.LLRESETJL(JL)
    ENDDO


    IF(LLRESET) THEN
      DO JK=KLEV,1,-1
        DO JL=KIDIA,KFDIA
          IF ( LLRESETJL(JL) ) THEN 
            JKT=ICTOP(JL)
            JKB=IDPL(JL)
            IF ( JK<=JKB .AND. JK>=JKT ) THEN
              KLAB(JL,JK)=ILAB(JL,JK)
              PTU(JL,JK)=ZTU(JL,JK)
              PQU(JL,JK)=ZQU(JL,JK)
              PLU(JL,JK)=ZLU(JL,JK)
            ELSE 
              KLAB(JL,JK)=1
              PTU(JL,JK)=PTENH(JL,JK)
              PQU(JL,JK)=PQENH(JL,JK)
              PLU(JL,JK)=0.0_JPRB
            ENDIF
            IF ( JK<JKT ) KLAB(JL,JK)=0
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    DO JL=KIDIA,KFDIA
      IF ( LLDEEP(JL) .AND. LLFIRST(JL) ) THEN
        KDPL(JL)  = IDPL(JL)
        KCTOP(JL) = ICTOP(JL)
        KCBOT(JL) = ICBOT(JL)
        LDCUM(JL) = LLDCUM(JL)
        LDSC(JL)  = .FALSE.
        KBOTSC(JL)= -1
        JKB=KCBOT(JL)
        PWUBASE(JL)=SQRT(MAX(ZWU2H(JL,JKB),0.0_JPRB))
!  no initialization of wind for deep here, this is done in
!  CUINI and CUASCN
        LLFIRST(JL)=.FALSE.
      ENDIF
      LLGO_ON(JL) = .NOT.LLDEEP(JL)
    ENDDO
  ENDIF

ENDDO ! end of big loop for search of departure level     

      ! chose maximum CAPE value
DO JL=KIDIA,KFDIA
  PCAPE(JL) = MAXVAL(ZCAPE(JL,:))
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUBASENS',1,ZHOOK_HANDLE)
END SUBROUTINE CUBASENS
