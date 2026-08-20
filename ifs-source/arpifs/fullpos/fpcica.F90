! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPCICA(YDML_PHY_MF,KST,KEND,KPROMA,KLEV,KCAPETYPE,LDMUL,PENTRA,PMLDEP,PTCLS,PRPCLS,PRHCLS,PT,PRP,&
 & PQV,PCAPE,PTCVS,PCIN,KLCL,KFCL,KEL,PMLCAPE)

! --------------------------------------------------------------
! **** *FPCICA* COMPUTE CAPE AND CIN.
! --------------------------------------------------------------
! SUBJECT:
!    ROUTINE COMPUTING CAPE AND CIN FOR SELECTED "TYPE" (PARCEL)  

! INTERFACE:
!    *CALL* *FPCICA*

! --------------------------------------------------------------
! -   INPUT ARGUMENTS
!     ---------------

! - DIMENSIONING

! KST      : FIRST INDEX OF LOOPS
! KEND     : LAST INDEX OF LOOPS
! KPROMA   : DEPTH OF THE VECTORIZATION ARRAYS
! KLEV     : END OF VERTICAL LOOP AND VERTICAL DIMENSION

! - VARIABLES
! KCAPETYPE: TYPE OF CAPE COMPUTATION
! LDMUL   : BOOLEAN TO COMPUTE MOST UNSTABLE LAYER
! PENTRA   : ENTRAINMENT
! PMLDEP   : Mean Layer DEPth for MLCAPE computation.
! PTCLS    : CLS TEMPERATURE (K)
! PRPCLS   : CLS PRESSURE (PA)
! PRHCLS   : CLS RELATIVE HUMIDITY (NO DIM)
! PT       : TEMPERATURE (K)
! PRP      : PRESSURE (PA)
! PQV      : WATER VAPOUR SPECIFIC HUMIDITY (NO DIM)

! --------------------------------------------------------------
! -   OUTPUT ARGUMENTS
!     ---------------
! - VARIABLES
! PCAPE    : CAPE - CONVECTIVE AVAILABLE POTENTIAL ENERGY (J/KG)
!                   (POTENTIALLY AVAILABLE CONVECTIVE KINETIC ENERGY)
! PCIN     : CIN - CONVECTIVE INHIBITION (J/KG)
! PTCVS    : CONVECTIVE TEMPERATURE AT SCREEN LEVEL (K)
! PMLCAPE  : MLCAPE ! MEAN LAYER CAPE JUST ABOVE SURFACE (J/KG)

! --------------------------------------------------------------
! -   IMPLICITE ARGUMENTS
!     -------------------
! YOMCAPE
! --------------------------------------------------------------
! EXTERNALS:

! METHOD:
!   different types of CAPE and CIN are calculated:
!    - CAPE for the parcel at the lowest model level
!    - CAPE for the most unstable parcel
!            CAPEs are calculated for the parcels released from
!            all model levels where p > GCAPEPSD*ps.
!    - CAPE for the CLS parcel 

! AUTEUR/AUTHOR:   2001-03, N. PRISTOV

! MODIFICATIONS:
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P. Marquet    16-Jan-2009 (CY33t1) : add the option
!                     (KCAPETYPE=5) = the merge of the options
!                      2 and 3. Correction for the computations
!                      of KLCL,KFCL,KEL for the options 2 and 4,
!                      where the value for the "maximum CAPE"
!                      ought to be retained (level=IMAX(1)).
!                  !! (KCAPETYPE=4) is equivalent to KCAPETYPE=3
!                      : see "endpos" or "phymfpos" => use of 
!                      KCAPETYPE=5 !!
!   2010-03-09, J.M. Piriou: change definition of GCAPEPSD.    
!        P.Marguinaud  10-Aug-2010 Run without SURFEX
!   2018-09, R. Brozkova: Added convective temperature; fixed CIN/CAPE
!     calculation from most unstable layer (no dependence on NPROC).
!   2020-07-07, J.M. Piriou: MLCAPE (KCAPETYPE=6).
!      R. El Khatib 08-Jul-2022 Contribution to the encapsulation of YOMCST and YOETHF
! --------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : YDCST=>YRCST ! allows use of included functions. REK.
USE YOMCAPE  , ONLY : NCAPEPSD

IMPLICIT NONE

TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN):: YDML_PHY_MF
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCAPETYPE 
LOGICAL           ,INTENT(IN)    :: LDMUL 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PENTRA
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMLDEP
REAL(KIND=JPRB)   ,INTENT(IN), TARGET    :: PTCLS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPCLS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHCLS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRP(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQV(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAPE(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTCVS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT),   OPTIONAL :: PCIN(KPROMA) 
INTEGER(KIND=JPIM)  ,INTENT(OUT), OPTIONAL :: KLCL(KPROMA)
INTEGER(KIND=JPIM)  ,INTENT(OUT), OPTIONAL :: KFCL(KPROMA)
INTEGER(KIND=JPIM)  ,INTENT(OUT), OPTIONAL :: KEL(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT),   OPTIONAL :: PMLCAPE(KPROMA)

INTEGER(KIND=JPIM):: JLON,JLEV,ILEV,IMAX(1)

REAL(KIND=JPRB):: ZDELTA, ZEW, ZQS
REAL(KIND=JPRB):: ZT(KPROMA,KLEV+1)
REAL(KIND=JPRB):: ZP(KPROMA,KLEV+1)
REAL(KIND=JPRB):: ZQV(KPROMA,KLEV+1)
REAL(KIND=JPRB):: ZCAPE(KPROMA,KLEV),ZCIN(KPROMA,KLEV)
INTEGER(KIND=JPIM):: ILEVST
REAL(KIND=JPRB) :: ZENTRA

REAL   (KIND=JPRB) :: Z2CAPE(KPROMA),    Z2CIN(KPROMA), ZQVCLS(KPROMA)
INTEGER(KIND=JPIM) :: ILCL(KPROMA,KLEV), I2LCL(KPROMA)
INTEGER(KIND=JPIM) :: IFCL(KPROMA,KLEV), I2FCL(KPROMA)
INTEGER(KIND=JPIM) :: IEL (KPROMA,KLEV), I2EL (KPROMA)
INTEGER(KIND=JPIM) :: IMX (KPROMA)
REAL   (KIND=JPRB) :: ZIND2(KPROMA)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZCINOPT(KPROMA)
INTEGER(KIND=JPIM):: ILCLOPT(KPROMA),IFCLOPT(KPROMA),IELOPT(KPROMA)

REAL(KIND=JPRB), POINTER :: PP(:)

#include "fpcincape.intfb.h"
!---------

#include "fcttrm.func.h"

!-------------------------------------------------
! INITIALIZE TO ZERO.
!-------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPCICA',0,ZHOOK_HANDLE)
ASSOCIATE(NSURFEXCTL=>YDML_PHY_MF%YRMSE%NSURFEXCTL, NSURFEXCTLMAX=>YDML_PHY_MF%YRMSE%NSURFEXCTLMAX,&
 & RCPV=>YDCST%RCPV, RETV=>YDCST%RETV, RCW=>YDCST%RCW, RCS=>YDCST%RCS, RLVTT=>YDCST%RLVTT, &
 & RLSTT=>YDCST%RLSTT, RTT=>YDCST%RTT, RALPW=>YDCST%RALPW, RBETW=>YDCST%RBETW, &
 & RGAMW=>YDCST%RGAMW, RALPS=>YDCST%RALPS, RBETS=>YDCST%RBETS, RGAMS=>YDCST%RGAMS, &
 & RALPD=>YDCST%RALPD, RBETD=>YDCST%RBETD, RGAMD=>YDCST%RGAMD, RKAPPA=>YDCST%RKAPPA, RV=>YDCST%RV, &
 & LNEIGE=>YDML_PHY_MF%YRPHY%LNEIGE)

PCAPE(KST:KEND)=0.0_JPRB
PTCVS (KST:KEND)=0.0_JPRB
ZCAPE(KST:KEND,1:KLEV)=0.0_JPRB
ZCIN (KST:KEND,1:KLEV)=0.0_JPRB
Z2CAPE(KST:KEND)=0.0_JPRB
Z2CIN (KST:KEND)=0.0_JPRB

IF (NSURFEXCTL < NSURFEXCTLMAX) THEN

  PP => PTCLS
  WHERE (PP <= 0._JPRB)
    PP = 200._JPRB
  ENDWHERE

ENDIF

!-------------------------------------------------
! Compute convective temperature at screen level.
!-------------------------------------------------

DO JLON=KST,KEND
  ZDELTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-PTCLS(JLON)))
  ZEW= FOEW (PTCLS(JLON),ZDELTA)
  ZQVCLS(JLON)=ZEW*PRHCLS(JLON) &
   & /((RETV+1.0_JPRB)*PRPCLS(JLON)-RETV*ZEW*PRHCLS(JLON))
ENDDO

DO JLEV=KLEV,1,-1
  DO JLON=KST,KEND
    ZDELTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-PT(JLON,JLEV)))
    ZQS=FOQS(FOEW(PT(JLON,JLEV),ZDELTA)/PRP(JLON,JLEV))
    IF(ZQS <= ZQVCLS(JLON) .AND. PTCVS(JLON) == 0.0_JPRB) THEN
      PTCVS(JLON)=PT(JLON,JLEV)*(PRPCLS(JLON)/PRP(JLON,JLEV))**RKAPPA
    ENDIF
  ENDDO
ENDDO

!========================
 IF (KCAPETYPE == 1) THEN
!========================

!    - 1) CAPE for the parcel at the lowest model level

CALL FPCINCAPE(YDML_PHY_MF%YRTOPH,KST,KEND,KPROMA,KLEV,KLEV,PENTRA,PMLDEP,PT,PRP,PQV,PCAPE,ZCINOPT,ILCLOPT,IFCLOPT,IELOPT)

!============================
 ELSEIF (KCAPETYPE == 2 .OR. KCAPETYPE == 6) THEN
!============================

!    KCAPETYPE=2 : MUCAPE : CAPE for the most unstable parcel
!    KCAPETYPE=6 : MLCAPE : Mean Layer CAPE. It is also a most unstable one: MUMLCAPE.

  ! Avoid vertical loop if most unstable CAPEs are not required
  IF (LDMUL) THEN
    ILEV=NCAPEPSD
  ELSE
    ILEV=KLEV
  ENDIF
  DO JLEV=KLEV,ILEV,-1
    IF(KCAPETYPE == 6) THEN
      ! The starting point of ascent is a vertical mean centered on current model level.
      ! This results in a Mean Layer CAPE (MLCAPE).
      ILEVST=-JLEV
      ZENTRA=PENTRA
    ELSE
      ! The starting point of ascent is the current model level.
      ILEVST=JLEV
      ZENTRA=PENTRA
    ENDIF
    CALL FPCINCAPE(YDML_PHY_MF%YRTOPH,KST,KEND,KPROMA,KLEV,ILEVST,ZENTRA,PMLDEP,PT,PRP,PQV,&
     & ZCAPE(:,JLEV),ZCIN(:,JLEV),&
     & ILCL(:,JLEV),IFCL(:,JLEV),IEL(:,JLEV))
    IF(KCAPETYPE == 6 .AND. JLEV == KLEV .AND. PRESENT(PMLCAPE)) THEN
      ! The starting point of ascent is a vertical mean starting from lowest model level.
      ! This results in a Mean Layer CAPE (MLCAPE) just above surface.
      DO JLON=KST,KEND
        PMLCAPE(JLON)=ZCAPE(JLON,JLEV)
      ENDDO
    ENDIF
  ENDDO
  IF (LDMUL) THEN
    DO JLON=KST,KEND
      IMAX=MAXLOC(ZCAPE(JLON,ILEV:KLEV))+ILEV-1
      PCAPE(JLON)=ZCAPE(JLON,IMAX(1))
      ILCLOPT (JLON)=ILCL (JLON,IMAX(1))
      IFCLOPT (JLON)=IFCL (JLON,IMAX(1))
      IELOPT  (JLON)=IEL  (JLON,IMAX(1))
      ! CIN can be estimated in case of no vertical mean on starting parcel (KCAPETYPE = 2),
      ! CIN cannot be estimated in case of vertical mean on starting parcel (KCAPETYPE = 6).
      ZCINOPT(JLON)=ZCIN(JLON,IMAX(1))
    ENDDO
  ENDIF

!============================
 ELSEIF (KCAPETYPE == 3) THEN
!============================

!    - 3) CAPE for the CLS parcel 

  ILEV=KLEV+1
  DO JLON=KST,KEND
    ZP(JLON,ILEV)=PRPCLS(JLON)
    ZT(JLON,ILEV)=PTCLS(JLON)
    IF (LNEIGE) THEN
      ZDELTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-PTCLS(JLON)))
    ELSE
      ZDELTA=0.0_JPRB
    ENDIF
    ZEW= FOEW (PTCLS(JLON),ZDELTA)
    ZQV(JLON,ILEV)=ZEW*PRHCLS(JLON)&
     & /((RETV+1.0_JPRB)*PRPCLS(JLON)-RETV*ZEW*PRHCLS(JLON))  
  ENDDO
  DO JLEV=1,KLEV
    DO JLON=KST,KEND
      ZP(JLON,JLEV)=PRP(JLON,JLEV)
      ZT(JLON,JLEV)=PT(JLON,JLEV)
      ZQV(JLON,JLEV)=PQV(JLON,JLEV)
    ENDDO
  ENDDO
  CALL FPCINCAPE(YDML_PHY_MF%YRTOPH,KST,KEND,KPROMA,ILEV,ILEV,PENTRA,PMLDEP,ZT,ZP,ZQV,PCAPE,ZCINOPT,ILCLOPT,IFCLOPT,IELOPT)

!============================
 ELSEIF (KCAPETYPE == 5) THEN
!============================

!    - 5) CAPE for the most unstable parcel among 
!         a) the CLS parcel
!         and b) levels where p > GCAPEPSD*ps.

!    - 5a) CAPE for the CLS parcel.
!          -----------------------

  ILEV=KLEV+1
  DO JLON=KST,KEND
    ZP(JLON,ILEV)=PRPCLS(JLON)
    ZT(JLON,ILEV)=PTCLS(JLON)
    IF (LNEIGE) THEN
      ZDELTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-PTCLS(JLON)))
    ELSE
      ZDELTA=0.0_JPRB
    ENDIF
    ZEW= FOEW (PTCLS(JLON),ZDELTA)
    ZQV(JLON,ILEV)=ZEW*PRHCLS(JLON)&
     & /((RETV+1.0_JPRB)*PRPCLS(JLON)-RETV*ZEW*PRHCLS(JLON))  
  ENDDO
  DO JLEV=1,KLEV
    DO JLON=KST,KEND
      ZP(JLON,JLEV)=PRP(JLON,JLEV)
      ZT(JLON,JLEV)=PT(JLON,JLEV)
      ZQV(JLON,JLEV)=PQV(JLON,JLEV)
    ENDDO
  ENDDO
  CALL FPCINCAPE(YDML_PHY_MF%YRTOPH,KST,KEND,KPROMA,ILEV,ILEV,PENTRA,PMLDEP,ZT,ZP,ZQV,&
     &           Z2CAPE,Z2CIN,I2LCL,I2FCL,I2EL)

!    - 5b) CAPE for the most unstable parcel.
!          ---------------------------------

! for how many levels CAPE should be computed
  ILEV=NCAPEPSD
  DO JLEV=KLEV,ILEV,-1
    CALL FPCINCAPE(YDML_PHY_MF%YRTOPH,KST,KEND,KPROMA,KLEV,JLEV,PENTRA,PMLDEP,PT,PRP,PQV,&
     & ZCAPE(:,JLEV),ZCIN(:,JLEV),&
     & ILCL(:,JLEV),IFCL(:,JLEV),IEL(:,JLEV))
    DO JLON=KST,KEND
      ZCAPE(JLON,JLEV)=ZCAPE(JLON,JLEV)*MAX(0.0_JPRB,SIGN(1.0_JPRB,REAL(JLEV-NCAPEPSD,KIND=JPRB)))
    ENDDO
  ENDDO

!    - 5c) The true maximum value for the CAPE. 
!          -----------------------------------
!          i.e. either for the CLS value (Z2CAPE) 
!          or for the Max-Upper-Air value
!          (ZCAPE(IMX(JLON))) :

  DO JLON=KST,KEND
  ! IMAX = the index array for the Maximum
  !        value among the Upper-Air data :
    IMAX=MAXLOC(ZCAPE(JLON,ILEV:KLEV))+ILEV-1
    IMX(JLON)=IMAX(1)
  ENDDO
  DO JLON=KST,KEND
  ! ZIND2=1 if the CLS value is greater
  !         than the Max-Upper-Air value : 
    ZIND2(JLON)= MAX(0.0_JPRB, SIGN(1.0_JPRB,&
     &           Z2CAPE(JLON)-ZCAPE(JLON,IMX(JLON))&
     &               )              )
  ENDDO
  DO JLON=KST,KEND
    PCAPE(JLON) = Z2CAPE(JLON)                   *ZIND2(JLON)&
     &          + ZCAPE (JLON,IMX(JLON))*(1._JPRB-ZIND2(JLON))
    ZCINOPT(JLON)=Z2CIN (JLON)                   *ZIND2(JLON)&
     &          + ZCIN  (JLON,IMX(JLON))*(1._JPRB-ZIND2(JLON))
    ILCLOPT(JLON)=I2LCL (JLON)                   *ZIND2(JLON)&
     &          + ILCL  (JLON,IMX(JLON))*(1._JPRB-ZIND2(JLON))
    IFCLOPT(JLON)=I2FCL (JLON)                   *ZIND2(JLON)&
     &          + IFCL  (JLON,IMX(JLON))*(1._JPRB-ZIND2(JLON))
    IELOPT (JLON)=I2EL  (JLON)                   *ZIND2(JLON)&
     &          + IEL   (JLON,IMX(JLON))*(1._JPRB-ZIND2(JLON))
  ENDDO

 ELSE

   ZCINOPT(KST:KEND)=0.0_JPRB

! End of the test on KCAPETYPE
!=====
 ENDIF
!=====

IF (PRESENT(PCIN)) PCIN=ZCINOPT
IF (PRESENT(KLCL)) KLCL=ILCLOPT
IF (PRESENT(KFCL)) KFCL=IFCLOPT
IF (PRESENT(KEL))  KEL =IELOPT

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FPCICA',1,ZHOOK_HANDLE)
END SUBROUTINE FPCICA
