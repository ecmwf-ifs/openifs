! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUBASEN2 &
 & (YDECUMF2,YDEPHLI,  YDECLDP,&
 & KIDIA,    KFDIA,    KLON,    KLEV,   KINDEX,&
 & LDRAIN1D, LDCVOPT,  LDMIXS,&
 & PTENH,    PQENH,    PGEOH,    PAP,      PAPH,&
 & PQHFL,    PAHFS,&
 & PTEN,     PQEN,     PQSEN,    PGEO,&
 & PTU,      PQU,      PLU,&
 & PBUOH,    PWU2H,    PWUBASE,  PLGLAC,   PQPRCV,  PCAPE,&
 & KTYPE,    KLAB,     LDCUM,    LDSC,&
 & KCBOT,    KBOTSC,   KCTOP,    KSTUP,&
 & PTU5,     PQU5,     PLU5,&
 & PBUOH5,   PWU2H5,   PSUH5,    PWU2H15,&
 & PTVENH5,  PTVUH5,   PSUH25,   PQU25,&
 & PFACT35,  PFACT45,  PFACT55,  PWS5,&   
 & PKHVFL5,  PUST5,    PWUBASE5, PWORK15,  PWORK15T,&  
 & KCBOT5,   KLAB5,    LDCUM5,   LDTEST15,&
 & LDGO_ON5, LDMIDTEST5, LDCAPETEST5 )  

!          THIS ROUTINE CALCULATES CLOUD BASE FIELDS
!          CLOUD BASE HEIGHT AND CLOUD TOP HEIGHT

!          P. LOPEZ   ECMWF  (01/2002)
!          Inspired from A. Pier Siebesma  (KNMI)
!          and C. Jakob (ECMWF) (01/2001) 

!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD BASE AND CLOUD TOP VALUES FOR 
!          SIMPLIFIED CONVECTIVE PARAMETRIZATION

!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTRN2*.
!          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
!          IT RETURNS CLOUD FIELDS VALUES AND FLAGS AS FOLLOWS;
!                 KLAB=0 FOR STABLE LAYERS
!                 KLAB=1 FOR SUBCLOUD LEVELS
!                 KLAB=2 FOR CLOUD LEVELS LEVEL

!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD TOP.

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KINDEX*       TOP INDEX FOR OUTER VERTICAL LOOP

!    INPUT PARAMETERS (LOGICAL):
!    *LDRAIN1D*      ROUTINE IS CALLED FROM 1D-VAR.
!    *LDCVOPT*       STORE TRAJECTORY ARRAY TO AVOID RECOMPUTATIONS IN ADJOINT.
!    *LDMIXS*        WEAK (FALSE) OR STRONG (TRUE) CLOUD MIXING FOR SURFACE PARCEL ONLY

!    INPUT PARAMETERS (REAL):

! not used at the moment because we want to use linear intepolation
! for fields on the half levels.

!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG
!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
!    *PSSTRU*       KINEMATIC surface U-MOMENTUM FLUX             (M/S)^2
!    *PSSTRV*       KINEMATIC surface V-MOMENTUM FLUX             (M/S)^2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS             PA
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PQSEN*        PROVISIONAL ENVIRONMENT SAT. SPEC. HUMIDITY (T+1)  KG/KG
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2

!    UPDATED PARAMETERS (REAL):

!    *PTU*          TEMPERATURE IN UPDRAFTS                         K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG

!    UPDATED PARAMETERS (INTEGER):

!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!                        KLAB=2 FOR CLOUD LEVELS

!    OUTPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDSC*         FLAG: .TRUE. IF BL-CLOUDS EXIST

!    OUTPUT PARAMETERS (INTEGER):

!    *KTYPE*       TYPE OF CONVECTION    
!    *KSTUP*       LEVEL AT WHICH THE UPDRAFT STARTS
!    *KCBOT*       CLOUD BASE LEVEL !    
!    *KCTOP*       CLOUD TOP LEVEL = HEIGHEST HALF LEVEL 
!                  WITH A NON-ZERO CLOUD UPDRAFT.
!    *KBOTSC*      CLOUD BASE LEVEL OF BL-CLOUDS

!    OUTPUT PARAMETERS (REAL):

!    *PBUOH*       BUOYANCY ON HALF LEVELS                       M/S2
!    *PWU2H*       UPDRAFT KINETIC ENERGY * 2                   M2/S2
!    *PWUBASE*     UPDRAFT VERTICAL VELOCITY AT CLOUD BASE       M/S
!    *PLGLAC*      FROZEN CLOUD WATER CONTENT                   KG/KG 
!    *PQPRCV*      CONVECTIVE PRECIPITATION CONTENT             KG/KG 
!    *PCAPE*       CONVECTIVE AVAILABLE POTENTIAL ENERGY        J/KG 

!          EXTERNALS
!          ---------
!          *CUPDRA*  FOR COMPUTING CHARACTERISTICS OF UPDRAFT

!          MODIFICATIONS
!          -------------
!          P.Lopez      01-Oct-2005 New version of linearized convection
!          P. Lopez     11-01-2007  Added LDRAIN switch       
!          P.Lopez      11-Apr-2007 Optimization (LDCVOPT)
!          P.Lopez      19-Oct-2007 Revised version to improve match to Tiedtke scheme
!          A.Geer       01-OCt-2008 LDRAIN1D name change to reflect usage
!          P.Lopez      15-Oct-2015 Added CAPE computation (excl. liquid water loading)
!          P.Lopez      01-Dec-2020 Added arguments KINDEX, LDMIXS and PWU2H

!----------------------------------------------------------------------

USE YOECLDP   , ONLY : TECLDP
USE YOEPHLI   , ONLY : TEPHLI
USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST    , ONLY : RCPD, RETV, RD, RG
USE PARPHY    , ONLY : RKAP
USE YOECUMF2  , ONLY : TECUMF2

IMPLICIT NONE

TYPE(TECLDP)       ,INTENT(INOUT) :: YDECLDP
TYPE(TECUMF2)      ,INTENT(INOUT) :: YDECUMF2
TYPE(TEPHLI)       ,INTENT(INOUT) :: YDEPHLI
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM), INTENT(IN)    :: KINDEX
LOGICAL            ,INTENT(IN)    :: LDRAIN1D 
LOGICAL            ,INTENT(IN)    :: LDCVOPT
LOGICAL            ,INTENT(IN)    :: LDMIXS
REAL(KIND=JPRB)    ,INTENT(IN)    :: PTENH(KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PQENH(KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PQHFL(KLON,KLEV+1) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PAHFS(KLON,KLEV+1) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PTEN(KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PQEN(KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PQSEN(KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PGEO(KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PTU(KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PQU(KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PLU(KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PBUOH(KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PWU2H(KLON,KLEV)
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PWUBASE(KLON) 
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PLGLAC(KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PQPRCV(KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PCAPE(KLON) 
INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KTYPE(KLON) 
INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KLAB(KLON,KLEV) 
LOGICAL            ,INTENT(OUT)   :: LDCUM(KLON) 
LOGICAL            ,INTENT(OUT)   :: LDSC(KLON) 
INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KCBOT(KLON) 
INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KBOTSC(KLON) 
INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KCTOP(KLON) 
INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KSTUP(KLON) 
! arrays for optimization (storage of trajectory for adjoint)
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PTU5   (KLON,KLEV)    
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PQU5   (KLON,KLEV)    
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PLU5   (KLON,KLEV)    
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PBUOH5 (KLON,KLEV)
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PWU2H5 (KLON,KLEV)
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PSUH5  (KLON,KLEV)
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PWU2H15(KLON)
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PTVENH5(KLON,KLEV)    
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PTVUH5 (KLON,KLEV)    
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PSUH25 (KLON,KLEV)    
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PQU25  (KLON,KLEV)    
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PFACT35(KLON)    
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PFACT45(KLON)    
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PFACT55(KLON)    
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PWS5   (KLON)    
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PKHVFL5(KLON)    
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PUST5  (KLON)    
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PWUBASE5(KLON)    
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PWORK15(KLON,KLEV)    
REAL(KIND=JPRB)    ,INTENT(OUT)   :: PWORK15T(KLON,3,KLEV)    
INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KCBOT5 (KLON)
INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KLAB5  (KLON,KLEV)
LOGICAL            ,INTENT(OUT)   :: LDCUM5 (KLON)
LOGICAL            ,INTENT(OUT)   :: LDTEST15(KLON,KLEV)
LOGICAL            ,INTENT(OUT)   :: LDGO_ON5(KLON,KLEV)
LOGICAL            ,INTENT(OUT)   :: LDMIDTEST5(KLON,KLEV)
LOGICAL            ,INTENT(OUT)   :: LDCAPETEST5(KLON,KLEV)

!             LOCAL STORAGE
!             ----- -------

LOGICAL ::   LLGO_ON(KLON),  LLMID(KLON)

INTEGER(KIND=JPIM) :: IKSTART, JK, JL, IKB, JKK, JKT1

REAL(KIND=JPRB)    :: ZCAPE(KLON)
REAL(KIND=JPRB)    :: ZSENH(KLON,KLEV+1),  ZSUH(KLON,KLEV)

REAL(KIND=JPRB)     :: ZTU(KLON,KLEV) 
REAL(KIND=JPRB)     :: ZQU(KLON,KLEV) 
REAL(KIND=JPRB)     :: ZLU(KLON,KLEV) 
REAL(KIND=JPRB)     :: ZBUOH(KLON,KLEV) 
REAL(KIND=JPRB)     :: ZLGLAC(KLON,KLEV) 
REAL(KIND=JPRB)     :: ZQPRCV(KLON,KLEV)
REAL(KIND=JPRB)     :: ZZWU2H(KLON,KLEV)
REAL(KIND=JPRB)     :: ZZSUH(KLON,KLEV) 
INTEGER(KIND=JPIM)  :: ILAB(KLON,KLEV) 
INTEGER(KIND=JPIM)  :: ICBOT(KLON) 
INTEGER(KIND=JPIM)  :: IBOTSC(KLON) 
INTEGER(KIND=JPIM)  :: ICTOP(KLON) 
INTEGER(KIND=JPIM)  :: ISTUP(KLON)
LOGICAL             :: LLCUM(KLON) 
LOGICAL             :: LLSC(KLON) 
LOGICAL             :: LLFOUND_CV(KLON) 

REAL(KIND=JPRB)     :: ZTU0(KLON,KLEV)
REAL(KIND=JPRB)     :: ZQU0(KLON,KLEV)

REAL(KIND=JPRB) :: ZDEPTH    ! Convective depth
REAL(KIND=JPRB) :: ZRHO      ! DENSITY AT SURFACE                             (KG/M3) 
REAL(KIND=JPRB) :: ZUST      ! U_STAR                                         (M/S)
REAL(KIND=JPRB) :: ZKHVFL    ! SURFACE BUOYANCY FLUX                          (K M/S)
REAL(KIND=JPRB) :: ZWS       ! SIGMA_W            AT LOWEST MODEL HALFLEVEL   (M/S)
REAL(KIND=JPRB) :: ZQEXC     ! HUMIDITY    EXCESS AT LOWEST MODEL HALFLEVEL   (KG/KG)
REAL(KIND=JPRB) :: ZTEXC     ! TEMPERATURE EXCESS AT LOWEST MODEL HALFLEVEL   (K)
REAL(KIND=JPRB) :: ZTVENH    ! ENVIRONMENT VIRTUAL TEMPERATURE AT HALF LEVELS (K)  
REAL(KIND=JPRB) :: ZTVUH     ! UPDRAFT VIRTUAL TEMPERATURE AT HALF LEVELS     (K)
REAL(KIND=JPRB) :: ZWORK1, ZWORK2, ZRCPD, ZFACT3, ZFACT4, ZFACT5
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
LOGICAL         :: LLTEST1, LLTEST2

#include "cupdra.intfb.h"

!----------------------------------------------------------------------
!     0.           INITIALIZE CONSTANTS AND FIELDS
!                  -------------------------------
!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUBASEN2',0,ZHOOK_HANDLE)
ASSOCIATE(LMFMID2=>YDECUMF2%LMFMID2, RDEPTHS2=>YDECUMF2%RDEPTHS2)

JKT1=KINDEX
ZRCPD=1.0_JPRB/RCPD

DO JL=KIDIA,KFDIA
  LLGO_ON(JL)=.FALSE.
  LLMID(JL)=.FALSE.
ENDDO
 
DO JL=KIDIA,KFDIA
  PWUBASE(JL)=0.0_JPRB
  PCAPE(JL)=0.0_JPRB
ENDDO

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PBUOH(JL,JK)=0.0_JPRB
    PQPRCV(JL,JK)=0.0_JPRB
    PLGLAC(JL,JK)=0.0_JPRB   
  ENDDO
ENDDO

!----------------------------------------------------------------------
!       -----------------------------------------------------------
!       1.1  PREPARE FIELDS ON HALF LEVELS BY LINEAR INTERPOLATION
!             OF SPECIFIC HUMIDITY AND STATIC ENERGY
!       -----------------------------------------------------------

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZSUH(JL,JK)=0.0_JPRB
    PWU2H(JL,JK)=0.0_JPRB 
    ZSENH(JL,JK)=RCPD*PTENH(JL,JK)+PGEOH(JL,JK)
    ZTU0(JL,JK)=PTU(JL,JK)
    ZQU0(JL,JK)=PQU(JL,JK)
  ENDDO
ENDDO

! Initialize trajectory arrays (to be used in adjoint)
IF (LDCVOPT) THEN
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      LDTEST15(JL,JK)=.FALSE.
      LDMIDTEST5(JL,JK)=.FALSE.
      LDGO_ON5(JL,JK)=.FALSE.
      LDCAPETEST5(JL,JK)=.FALSE.
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
    LDCUM5(JL)=.FALSE.
  ENDDO
ENDIF

KCBOT  (:)  =-1      ! cloud base level for convection, (-1 if not found)
KSTUP  (:)  =-1      ! level at which the updraft starts (-1 if not found)
KBOTSC (:)  =-1      ! sc    base level for sc-clouds , (-1 if not found)
KCTOP  (:)  =-1      ! cloud top for convection (-1 if not found)
KTYPE  (:)  =0       ! convection type (0 if no convection)
LDCUM  (:)  =.FALSE. ! on exit: true if cloudbase=found
LDSC   (:)  =.FALSE. ! on exit: true if cloudbase=found

LLFOUND_CV(:)=.FALSE.

DO JK=KLEV-1,JKT1,-1

! Initialize local arrays
  ICBOT(:)=-1
  ISTUP(:)=-1
  IBOTSC(:)=-1
  ICTOP(:)=-1
  LLCUM(:)=.FALSE.
  LLSC(:)=.FALSE.
  
  DO JKK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZBUOH (JL,JKK)=0.0_JPRB
      ZQPRCV(JL,JKK)=0.0_JPRB
      ZLGLAC(JL,JKK)=0.0_JPRB  
    ENDDO
  ENDDO
  DO JKK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZTU (JL,JKK)=ZTU0(JL,JKK)
      ZQU (JL,JKK)=ZQU0(JL,JKK)
    ENDDO
  ENDDO
  DO JKK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZLU   (JL,JKK)=0.0_JPRB
      ZZWU2H(JL,JKK)=0.0_JPRB
      ZZSUH (JL,JKK)=0.0_JPRB
    ENDDO
  ENDDO

  IF (JK == KLEV-1) THEN

!        ---------------------------------------------------------
!        1.2    INITIALISE FIELDS AT LOWEST HALF MODEL LEVEL
!        ---------------------------------------------------------

!LOP    LL_LDO_ASC=.TRUE.

    DO JL=KIDIA,KFDIA
  
      ZRHO  = PAPH(JL,KLEV+1)/(RD*(PTEN(JL,KLEV)*(1.+RETV*PQEN(JL,KLEV))))

!LOP      ZUST  = MAX(SQRT(PSSTRU(JL)**2 + PSSTRV(JL)**2),ZREPUST)
      ZUST  = 0.1_JPRB

      ZKHVFL  = (PAHFS(JL,KLEV+1)*ZRCPD+&
       & RETV*PTEN(JL,KLEV)*PQHFL(JL,KLEV+1)&
       & )/ZRHO  
      ZFACT3 = ZUST**3 - 1.5_JPRB*RKAP*ZKHVFL*PGEOH(JL,KLEV)/PTEN(JL,KLEV)  
      ZWS = 1.2_JPRB * MAX(ZFACT3,0.0_JPRB)**(1.0_JPRB/3._JPRB)
! Store trajectory for adjoint (to save time in CUBASEN2AD)
      IF (LDCVOPT) THEN 
        PUST5(JL) = ZUST
        PKHVFL5(JL) = ZKHVFL
        PFACT35(JL) = ZFACT3
        PWS5(JL) = ZWS
      ENDIF

      IF (ZKHVFL < 0.0_JPRB) THEN
        LLGO_ON(JL)  = .TRUE.
        ILAB(JL,KLEV)= 1
        ZFACT4 = -1.5_JPRB*PAHFS(JL,KLEV+1)/(ZRHO*ZWS*RCPD)
        ZTEXC = MAX(ZFACT4,0.0_JPRB)
        ZFACT5 = -1.5_JPRB*PQHFL(JL,KLEV+1)/(ZRHO*ZWS)
        ZQEXC = MAX(ZFACT5,0.0_JPRB)
        ZQU (JL,KLEV) = PQENH(JL,KLEV) + ZQEXC
        ZZSUH (JL,KLEV) = ZSENH(JL,KLEV) + RCPD*ZTEXC
        ZTU (JL,KLEV) = (ZSENH(JL,KLEV)-PGEOH(JL,KLEV))*ZRCPD + ZTEXC
        ZZWU2H(JL,KLEV) = ZWS*ZWS

!  determine buoyancy at lowest half level

        ZTVENH = (1.0_JPRB+RETV*PQENH(JL,KLEV)) &
         & *(ZSENH(JL,KLEV)-PGEOH(JL,KLEV))*ZRCPD  
        ZTVUH = (1.0_JPRB+RETV*ZQU(JL,KLEV))*ZTU(JL,KLEV)
        ZBUOH(JL,KLEV) = (ZTVUH-ZTVENH)*RG/ZTVENH

! Store trajectory for adjoint (to save time in CUBASEN2AD)
        IF (LDCVOPT) THEN
          PFACT45(JL) = ZFACT4
          PFACT55(JL) = ZFACT5
          PTVENH5(JL,JK) = ZTVENH
          PTVUH5(JL,JK) = ZTVUH
        ENDIF
      ELSE
        LLGO_ON(JL)=.FALSE.    ! non-convective point
      ENDIF
    ENDDO

  ELSE

!        -------------------------------------------------
!        1.3    CHECK FOR "MID-LEVEL" CONVECTION 
!               AND INITIALISE FIELDS BEFORE ASCENT 
!               (SIMPLIFIED SCHEME)
!        -------------------------------------------------
 
!LOP    LL_LDO_ASC=.TRUE.

    DO JL=KIDIA,KFDIA
      
      LLTEST1 = .NOT.LLFOUND_CV(JL) .AND. PGEOH(JL,JK+1) < 1.5E+05_JPRB 
! Store trajectory for adjoint (to save time in CUBASEN2AD)
      IF (LDCVOPT) LDMIDTEST5(JL,JK)=LLTEST1

      IF (LLTEST1) THEN
      
        LLGO_ON(JL)=.TRUE.

        ZTEXC=0.2_JPRB
        ZQEXC=1.E-04_JPRB
        ZQU (JL,JK+1)=PQENH(JL,JK+1) + ZQEXC
        ZZSUH (JL,JK+1)=ZSENH(JL,JK+1) + RCPD*ZTEXC
        ZTU (JL,JK+1)=(ZSENH(JL,JK+1)-PGEOH(JL,JK+1))*ZRCPD + ZTEXC
        ZLU (JL,JK+1)=0.0_JPRB
        ZZWU2H(JL,JK+1)=1.0_JPRB
        ILAB(JL,JK+1)=1

! construct mixed layer for parcels emanating in lowest 60 hPa
        IF (PAPH(JL,KLEV+1)-PAPH(JL,JK)<60.E2_JPRB) THEN
          ZQU(JL,JK+1)  =0.0_JPRB
          ZZSUH(JL,JK+1)=0.0_JPRB
          ZWORK1        =0.0_JPRB
          DO JKK=JK+2,JK,-1
! Store trajectory for adjoint (to save time in CUBASEN2AD)
            IF (LDCVOPT) PWORK15T(JL,JKK-JK+1,JK)= ZWORK1
            IF( ZWORK1 < 50.E2_JPRB ) THEN
              ZWORK2=PAPH(JL,JKK)-PAPH(JL,JKK-1)
              ZWORK1        =ZWORK1+ZWORK2
              ZQU(JL,JK+1)  =ZQU(JL,JK+1) +PQENH(JL,JKK)*ZWORK2
              ZZSUH(JL,JK+1)=ZZSUH(JL,JK+1)+ZSENH(JL,JKK)*ZWORK2
            ENDIF
          ENDDO
! Store trajectory for adjoint (to save time in CUBASEN2AD)
          IF (LDCVOPT) THEN
            PSUH25(JL,JK) = ZZSUH(JL,JK+1)
            PQU25(JL,JK) = ZQU(JL,JK+1)
            PWORK15(JL,JK) = ZWORK1
          ENDIF
          ZQU(JL,JK+1)  =ZQU(JL,JK+1)/ZWORK1+ZQEXC
          ZZSUH(JL,JK+1)=ZZSUH(JL,JK+1)/ZWORK1+RCPD*ZTEXC
          ZTU(JL,JK+1)  =(ZZSUH(JL,JK+1)-PGEOH(JL,JK+1))*ZRCPD+ZTEXC
        ENDIF  
       
!  determine buoyancy 
 
        ZTVENH = (1.0_JPRB+RETV*PQENH(JL,JK+1)) &
         & *(ZSENH(JL,JK+1)-PGEOH(JL,JK+1))*ZRCPD  
        ZTVUH = (1.0_JPRB+RETV*ZQU(JL,JK+1))*ZTU(JL,JK+1)
        ZBUOH(JL,JK+1) = (ZTVUH-ZTVENH)*RG/ZTVENH

! Store trajectory for adjoint (to save time in CUBASEN2AD)
        IF (LDCVOPT) THEN
          PTVENH5(JL,JK) = ZTVENH
          PTVUH5(JL,JK) = ZTVUH
        ENDIF
      ELSE
        LLGO_ON(JL)=.FALSE.
      ENDIF

    ENDDO

  ENDIF

!----------------------------------------------------------------------

!     2.0          DO ASCENT IN SUBCLOUD AND LAYER,
!                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
!                  ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!                  CHECK FOR BUOYANCY AND SET FLAGS
!                  -------------------------------------
!       ------------------------------------------------------------
!        2.1  DO THE VERTICAL ASCENT UNTIL VELOCITY BECOMES NEGATIVE
!       ------------------------------------------------------------

!LOP  IF (LL_LDO_ASC) THEN

    IKSTART=JK

! Store trajectory for adjoint (to save time in CUBASEN2AD)
    IF (LDCVOPT) THEN
      KLAB5(:,JK) = ILAB(:,IKSTART+1)
      PTU5(:,JK)  = ZTU(:,IKSTART+1) 
      PQU5(:,JK)  = ZQU(:,IKSTART+1) 
      PLU5(:,JK)  = ZLU(:,IKSTART+1) 
      PSUH5(:,JK) = ZZSUH(:,IKSTART+1) 
      PWU2H5(:,JK)= ZZWU2H(:,IKSTART+1) 
      PBUOH5(:,JK)= ZBUOH(:,IKSTART+1)
      LDGO_ON5(:,JK)=LLGO_ON(:)
    ENDIF

    CALL CUPDRA &
     & ( YDECUMF2,YDEPHLI,YDECLDP,  KIDIA,    KFDIA,    KLON,    KLEV,   KINDEX,&
     & IKSTART,  LLGO_ON,  LDRAIN1D, LDMIXS,&
     & PQENH,    PGEO,     PGEOH,    PAPH,    PQSEN,&
     & ZTU,      ZQU,      ZLU,&
     & ZSENH,    ZZWU2H,   ZZSUH,&
     & ZBUOH,    ZLGLAC,   ZQPRCV,   ZCAPE,&
     & ILAB,     LLCUM,    LLSC,     ICBOT,&
     & IBOTSC,   ICTOP,    ISTUP  )  

    DO JL=KIDIA,KFDIA
      IF (LLCUM(JL)) THEN
        ZDEPTH=(PAPH(JL,ICBOT(JL))-PAPH(JL,ICTOP(JL)))
        LLTEST1=(IKSTART == KLEV-1 .AND. ZDEPTH <  RDEPTHS2) .OR. &
              & (IKSTART  < KLEV-1 .AND. ZDEPTH >= RDEPTHS2)
! Store trajectory for adjoint (to save time in CUBASEN2AD)
        IF (LDCVOPT) LDTEST15(JL,JK)=LLTEST1
        IF (LLTEST1) THEN
          DO JKK=1,KLEV
            PTU(JL,JKK)=ZTU(JL,JKK)
            PQU(JL,JKK)=ZQU(JL,JKK)
            PLU(JL,JKK)=ZLU(JL,JKK)
            PWU2H(JL,JKK)=ZZWU2H(JL,JKK)
            ZSUH(JL,JKK)=ZZSUH(JL,JKK)
            PBUOH(JL,JKK)=ZBUOH(JL,JKK)
            PLGLAC(JL,JKK)=ZLGLAC(JL,JKK)
            PQPRCV(JL,JKK)=ZQPRCV(JL,JKK)
            KLAB(JL,JKK)=ILAB(JL,JKK)
          ENDDO
          LDCUM(JL)=LLCUM(JL)
          LDSC(JL)=LLSC(JL)
          KCBOT(JL)=ICBOT(JL)
          KBOTSC(JL)=IBOTSC(JL)
          KCTOP(JL)=ICTOP(JL)
          KSTUP(JL)=ISTUP(JL)
          IF (IKSTART < KLEV-1) LLFOUND_CV(JL)=.TRUE.         
        ENDIF
      ENDIF
    ENDDO

! Update CAPE as maximum over all departure levels.
    DO JL=KIDIA,KFDIA
      LLTEST2=(LLCUM(JL) .AND. ZCAPE(JL) > PCAPE(JL))
! Store trajectory for adjoint (to save time in CUBASEN2AD)
      IF (LDCVOPT) LDCAPETEST5(JL,JK)=LLTEST2
      IF (LLTEST2) THEN
        PCAPE(JL) = ZCAPE(JL)
      ENDIF
    ENDDO

!LOP    LL_LDO_ASC=.FALSE.

!LOP  ENDIF

ENDDO

!----------------------------------------------------------------------

! Set flag for mid-level convection in simplified scheme

IF (LMFMID2) THEN
  DO JL=KIDIA,KFDIA
    IF (LDCUM(JL).AND.LLMID(JL)) KTYPE(JL)=3
  ENDDO
ENDIF

! Set updraft vertical velocity at cloud base 

DO JL=KIDIA,KFDIA
  IF (LDCUM(JL)) THEN
    IKB=KCBOT(JL)
    PWUBASE(JL)=SQRT(MAX(PWU2H(JL,IKB),0.0_JPRB))
  ENDIF
ENDDO

! Store trajectory for adjoint (to save time in CUBASEN2AD)
IF (LDCVOPT) THEN
  DO JL=KIDIA,KFDIA
    IF (LDCUM(JL)) THEN
      IKB=KCBOT(JL)
      PWU2H15(JL)=PWU2H(JL,IKB)
      PWUBASE5(JL)=PWUBASE(JL)
    ENDIF
    KCBOT5(JL)=KCBOT(JL)
    LDCUM5(JL)=LDCUM(JL)
  ENDDO
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUBASEN2',1,ZHOOK_HANDLE)
END SUBROUTINE CUBASEN2
