! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPCICA(YDML_PHY_MF,KST,KEND,KPROMA,KLEV,KCAPETYPE,PENTRA,PTCLS,PRPCLS,PRHCLS,PT,PRP,&
 & PQV,PKAP,PCAPE,PCIN,PTCVS,KLCL,KFCL,KEL)  

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
! PENTRA   : ENTRAINMENT
! PTCLS    : CLS TEMPERATURE (K)
! PRPCLS   : CLS PRESSURE (PA)
! PRHCLS   : CLS RELATIVE HUMIDITY (NO DIM)
! PT       : TEMPERATURE (K)
! PRP      : PRESSURE (PA)
! PQV      : WATER VAPOUR SPECIFIC HUMIDITY (NO DIM)
! PKAP     : KAPPA (used for mixing for MLCAPE case)

! --------------------------------------------------------------
! -   OUTPUT ARGUMENTS
!     ---------------
! - VARIABLES
! PCAPE    : CAPE - CONVECTIVE AVAILABLE POTENTIAL ENERGY (J/KG)
!                   (POTENTIALLY AVAILABLE CONVECTIVE KINETIC ENERGY)
! PCIN     : CIN - CONVECTIVE INHIBITION (J/KG)
! PTCVS    : CONVECTIVE TEMPERATURE AT SCREEN LEVEL (K)

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
!    - CAPE for the mixed layer (ML CAPE), where the mixing depth pressure
!            is defined as  GCAPEPSD*ps above ground

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
!   2019-09, C.Wittmann, J.Cedilnik: Added KCAPETYPE=6 (mixed layer CAPE)   
! --------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RV       ,RCPV     ,RETV     ,&
 & RCW      ,RCS      ,RLVTT    ,RLSTT    ,RTT      ,&
 & RALPW    ,RBETW    ,RGAMW    ,RALPS    ,RBETS    ,&
 & RGAMS    ,RALPD    ,RBETD    ,RGAMD    ,RKAPPA   ,RATM
USE YOMCAPE  , ONLY : NCAPEPSD  

IMPLICIT NONE

TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN):: YDML_PHY_MF
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCAPETYPE 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PENTRA
REAL(KIND=JPRB)   ,INTENT(IN), TARGET    :: PTCLS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPCLS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHCLS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRP(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQV(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAP(KPROMA,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAPE(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCIN(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTCVS(KPROMA)
INTEGER(KIND=JPIM)  ,INTENT(OUT)   :: KLCL(KPROMA)
INTEGER(KIND=JPIM)  ,INTENT(OUT)   :: KFCL(KPROMA)
INTEGER(KIND=JPIM)  ,INTENT(OUT)   :: KEL(KPROMA)

INTEGER(KIND=JPIM):: JLON,JLEV,ILEV,IMAX(1),IL(KPROMA)

REAL(KIND=JPRB):: ZDELTA, ZEW, ZQS
REAL(KIND=JPRB):: ZT(KPROMA,KLEV+1)
REAL(KIND=JPRB):: ZP(KPROMA,KLEV+1)
REAL(KIND=JPRB):: ZQV(KPROMA,KLEV+1)
REAL(KIND=JPRB):: ZCAPE(KPROMA,KLEV),ZCIN(KPROMA,KLEV)

REAL   (KIND=JPRB) :: Z2CAPE(KPROMA),    Z2CIN(KPROMA), ZQVCLS(KPROMA)
INTEGER(KIND=JPIM) :: ILCL(KPROMA,KLEV), I2LCL(KPROMA)
INTEGER(KIND=JPIM) :: IFCL(KPROMA,KLEV), I2FCL(KPROMA)
INTEGER(KIND=JPIM) :: IEL (KPROMA,KLEV), I2EL (KPROMA)
INTEGER(KIND=JPIM) :: IMX (KPROMA)
REAL   (KIND=JPRB) :: ZIND2(KPROMA)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB), POINTER :: PP(:)

INTEGER(KIND=JPIM) :: ILEVCOUNT(KPROMA)                                                              
REAL(KIND=JPRB) :: ZTHSUM(KPROMA),ZQVSUM(KPROMA),ZTEST(KPROMA),ZUSRATM   

#include "fpcincape.intfb.h"
!---------

#include "fcttrm.func.h"

!-------------------------------------------------
! INITIALIZE TO ZERO.
!-------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPCICA',0,ZHOOK_HANDLE)
ASSOCIATE(NSURFEXCTL=>YDML_PHY_MF%YRMSE%NSURFEXCTL, NSURFEXCTLMAX=>YDML_PHY_MF%YRMSE%NSURFEXCTLMAX,&
 & LNEIGE=>YDML_PHY_MF%YRPHY%LNEIGE)
PCAPE(KST:KEND)=0.0_JPRB
PCIN (KST:KEND)=0.0_JPRB
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

CALL FPCINCAPE(YDML_PHY_MF%YRTOPH,KST,KEND,KPROMA,KLEV,KLEV,PENTRA,PT,PRP,PQV,PCAPE,PCIN,KLCL,KFCL,KEL)

!============================
 ELSEIF (KCAPETYPE == 2) THEN
!============================

!    - 2) CAPE for the most unstable parcel
  ILEV=NCAPEPSD
  DO JLEV=KLEV,ILEV,-1
    CALL FPCINCAPE(YDML_PHY_MF%YRTOPH,KST,KEND,KPROMA,KLEV,JLEV,PENTRA,PT,PRP,PQV,&
     & ZCAPE(:,JLEV),ZCIN(:,JLEV),&
     & ILCL(:,JLEV),IFCL(:,JLEV),IEL(:,JLEV))
    DO JLON=KST,KEND
      ZCAPE(JLON,JLEV)=ZCAPE(JLON,JLEV)*MAX(0.0_JPRB,SIGN(1.0_JPRB,REAL(JLEV-NCAPEPSD,KIND=JPRB)))
    ENDDO
  ENDDO
  DO JLON=KST,KEND
    IMAX=MAXLOC(ZCAPE(JLON,ILEV:KLEV))+ILEV-1
    PCAPE(JLON)=ZCAPE(JLON,IMAX(1))
    PCIN (JLON)=ZCIN (JLON,IMAX(1))
    KLCL (JLON)=ILCL (JLON,IMAX(1))
    KFCL (JLON)=IFCL (JLON,IMAX(1))
    KEL  (JLON)=IEL  (JLON,IMAX(1))
  ENDDO

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
  CALL FPCINCAPE(YDML_PHY_MF%YRTOPH,KST,KEND,KPROMA,ILEV,ILEV,PENTRA,ZT,ZP,ZQV,PCAPE,PCIN,KLCL,KFCL,KEL)

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
  CALL FPCINCAPE(YDML_PHY_MF%YRTOPH,KST,KEND,KPROMA,ILEV,ILEV,PENTRA,ZT,ZP,ZQV,&
     &           Z2CAPE,Z2CIN,I2LCL,I2FCL,I2EL)

!    - 5b) CAPE for the most unstable parcel.
!          ---------------------------------

! for how many levels CAPE should be computed
  ILEV=NCAPEPSD
  DO JLEV=KLEV,ILEV,-1
    CALL FPCINCAPE(YDML_PHY_MF%YRTOPH,KST,KEND,KPROMA,KLEV,JLEV,PENTRA,PT,PRP,PQV,&
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
    PCIN (JLON) = Z2CIN (JLON)                   *ZIND2(JLON)&
     &          + ZCIN  (JLON,IMX(JLON))*(1._JPRB-ZIND2(JLON))
    KLCL (JLON) = I2LCL (JLON)                   *ZIND2(JLON)&
     &          + ILCL  (JLON,IMX(JLON))*(1._JPRB-ZIND2(JLON))
    KFCL (JLON) = I2FCL (JLON)                   *ZIND2(JLON)&
     &          + IFCL  (JLON,IMX(JLON))*(1._JPRB-ZIND2(JLON))
    KEL  (JLON) = I2EL  (JLON)                   *ZIND2(JLON)&
     &          + IEL   (JLON,IMX(JLON))*(1._JPRB-ZIND2(JLON))
  ENDDO

  
!============================
 ELSEIF (KCAPETYPE == 6) THEN
!============================


  ! Mixed layer cape (ML CAPE); parcel entering CAPE computation is created by
  !  mixing of lowest GCAPEPSD*Ps Pascals (usually 50-100hPa)
  ! the temperature used for averaging is the potential temperature
  
  ! so let's mix the pbl (compute average for Theta and q):

  ILEVCOUNT(KST:KEND)=0
  ZTHSUM(KST:KEND)=0.0_JPRB
  ZQVSUM(KST:KEND)=0.0_JPRB
  ZTEST(KST:KEND)=1.0_JPRB
  
  ZUSRATM=1.0_JPRB/RATM
 
  DO JLEV=KLEV,1,-1
    
    DO JLON=KST,KEND
      ! conditional statements (if pressure above GCAPEPSD*Ps)
      ILEVCOUNT(JLON)=ILEVCOUNT(JLON)+NINT(ZTEST(JLON))

      ! compute average of theta
      ZTHSUM(JLON)=ZTHSUM(JLON)+ZTEST(JLON)*PT(JLON,JLEV)*(PRP(JLON,JLEV)/PRP(JLON,KLEV))**(-PKAP(JLON,JLEV)) 

      ! to compute average of qv
      ZQVSUM(JLON)=ZQVSUM(JLON)+ZTEST(JLON)*PQV(JLON,JLEV)

      ZTEST(JLON)=MAX(0.0_JPRB,SIGN(1.0_JPRB,REAL(JLEV-NCAPEPSD,KIND=JPRB)))

      ZP(JLON,JLEV)=PRP(JLON,JLEV)
      ZT(JLON,JLEV)=PT(JLON,JLEV)
      ZQV(JLON,JLEV)=PQV(JLON,JLEV) 
    ENDDO
  ENDDO

  ! build average: 
 
  DO JLON=KST,KEND
    
    ZP(JLON,KLEV+1)=PRPCLS(JLON)  
    ZT(JLON,KLEV+1)=ZTHSUM(JLON)/(ILEVCOUNT(JLON)*1.0_JPRB)
    ZQV(JLON,KLEV+1)=ZQVSUM(JLON)/(ILEVCOUNT(JLON)*1.0_JPRB)
       
  ENDDO
 
  ILEV=KLEV+1
  CALL FPCINCAPE(YDML_PHY_MF%YRTOPH,KST,KEND,KPROMA,ILEV,ILEV,PENTRA,ZT,ZP,ZQV,&
 &               PCAPE,PCIN,KLCL,KFCL,KEL)
  
! End of the test on KCAPETYPE
!=====
 ENDIF
!=====

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FPCICA',1,ZHOOK_HANDLE)
END SUBROUTINE FPCICA
