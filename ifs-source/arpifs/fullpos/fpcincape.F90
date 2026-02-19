! (C) Copyright 1989- Meteo-France.

!OPTIONS XOPT(NOEVAL)
SUBROUTINE FPCINCAPE(YDTOPH,KST,KEND,KPROMA,KLEV,KLEVST,PENTRA,PT,PRP,PQV,PCAPE,PCIN,KLCL,KFCL,KLNB)

! --------------------------------------------------------------
! **** *FPCINCAPE* COMPUTE CAPE AND CIN.
! --------------------------------------------------------------
! SUBJECT:
!    ROUTINE COMPUTING CAPE AND CIN  

! INTERFACE:
!    *CALL* *FPCINCAPE*

! --------------------------------------------------------------
! -   INPUT ARGUMENTS
!     ---------------

! - DIMENSIONING

! KST      : FIRST INDEX OF LOOPS
! KEND     : LAST INDEX OF LOOPS
! KPROMA   : DEPTH OF THE VECTORIZATION ARRAYS
! KLEV     : END OF VERTICAL LOOP AND VERTICAL DIMENSION

! - VARIABLES
! KLEVST   : LEVEL FROM WHICH PARCEL IS RAISED
! PENTRA   : ENTRAINMENT
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
! KLCL     : CONDENSATION LEVEL
! KFCL     : FREE CONVECTION LEVEL
! KLNB     : LEVEL OF NEUTRAL BUOYANCY

! --------------------------------------------------------------
! -   IMPLICITE ARGUMENTS
!     -------------------
! YOMCAPE 
! YOMCST
! FCTTRM
! FCTAST
! FCTTIM

! --------------------------------------------------------------
! EXTERNALS:

! METHOD:

!      THE PARCEL IS RAISED FROM LEVEL (LO which is equal KLEVST) 
!      TO THE LEVEL OF CONDENSATION (LC),
!      FURTHER TO ITS LEVEL OF FREE CONVECTION (LFC), WHERE THE 
!      PARCEL BECOMES BUOYANT,
!      THEN FURTHER TO THE LEVEL OF NEUTRAL BUOYANCY (LNB), WHERE THE
!      PARCEL BECOMES UNBUOYANT.

!      ALL COMPUTATIONS ARE DONE WITHOUT ENTRAINMENT OF ENVIRONMENTAL AIR.

!      CIN IS MASS SPECIFIC ENERGY TO RAISE THE PARCEL FROM 
!          from LO to LC further to LFC.
!          ONLY THE SUM OF NEGATIVE TERMS.
!      CAPE IS MASS SPECIFIC ENERGY  PROVIDED BY THE RAISE OF THE PARCEL
!          from LFC to LNB.
!          ONLY THE SUM OF POSITIVE TERMS.

! AUTEUR/AUTHOR:   2001-03, J.M. PIRIOU, N. PRISTOV.

! MODIFICATIONS:
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!    2010-02-04  J.M. Piriou. Bug correction: derivative of qs(T,p) from the Newton loops.
!    2010-02-11  J.M. Piriou. Bug correction: JLON dimension of ZFDERQS0 array.
!    2010-02-17  J.M. Piriou. Protections in case of cold temperature, low pressure, etc.
!    2010-02-17  J.M. Piriou. First Newton loop in a single DO loop.
!    2010-02-17  J.M. Piriou. Compute CAPE only below the convective reference level (NTCVIM).
! --------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCAPE  , ONLY :  NCAPEITER,NETAPES  ,GCAPERET
USE YOMCST   , ONLY :  RTT      ,RDAY     ,REPSM    ,RETV     ,&
 & RCW      ,REA      ,RD       ,RV       ,RCPD     ,RCPV     ,&
 & RCS      ,RLVTT    ,RLSTT    ,RBETS    ,RALPW    ,RBETW    ,&
 & RGAMW    ,RALPS    ,RGAMS    ,RALPD    ,RBETD    ,RGAMD    ,&
 & RG
USE YOMTOPH, ONLY : TTOPH

IMPLICIT NONE

TYPE(TTOPH)       ,INTENT(IN)    :: YDTOPH
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PENTRA
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRP(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQV(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAPE(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCIN(KPROMA) 
INTEGER(KIND=JPIM)   ,INTENT(OUT)   :: KLCL(KPROMA)
INTEGER(KIND=JPIM)   ,INTENT(OUT)   :: KFCL(KPROMA)
INTEGER(KIND=JPIM)   ,INTENT(OUT)   :: KLNB(KPROMA)
INTEGER(KIND=JPIM)      :: IPREVIOUS_NULL_ZRT(KPROMA)
INTEGER(KIND=JPIM) :: JLEV,JLON,JIT,JETAPES
REAL(KIND=JPRB) :: ZDLOG(KPROMA),ZBUOY,ZTV1,ZTV2,ZQV2,ZRT(KPROMA),ZT2,ZDT(KPROMA),ZST,ZDELARG,ZL,ZCP, &
 & ZFDERQS,ZADD,ZDERL,ZRDLOG,ZZQV
REAL(KIND=JPRB) :: ZBUOYPREC(KPROMA),ZQS(KPROMA),ZQV(KPROMA), ZQV1(KPROMA), &
 & ZT(KPROMA), ZT1(KPROMA),ZQL(KPROMA),ZQI(KPROMA),ZSQV(KPROMA), &
 & ZSQL(KPROMA),ZSQI(KPROMA),ZTDEPART(KPROMA),ZLOG(KPROMA), &
 & ZPARRIVEE(KPROMA),ZPDEPART(KPROMA),ZQVDEPART(KPROMA),ZZT(KPROMA)  
REAL(KIND=JPRB) :: ZFDERFOLH0(KPROMA),ZDELARG0(KPROMA),ZL0(KPROMA),ZCP0(KPROMA)
REAL(KIND=JPRB) :: ZFDERQS0, ZDZ, ZZFOEW(KPROMA)
REAL(KIND=JPRB) :: ZMAXT,ZMINT,ZMINDERI,ZMINQ,ZMAXQ
REAL(KIND=JPRB) :: ZTIN(KPROMA,KLEV) 
REAL(KIND=JPRB) :: ZQVIN(KPROMA,KLEV) 

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!---------

!  FUNCTIONS
#include "fctast.func.h"
#include "fcttrm.func.h"
#include "fcttim.func.h"

IF (LHOOK) CALL DR_HOOK('FPCINCAPE',0,ZHOOK_HANDLE)
ASSOCIATE(NTCVIM=>YDTOPH%NTCVIM)
!-------------------------------------------------
! INITIALIZE DEFAULT VALUES.
!-------------------------------------------------

PCIN(KST:KEND)=0.0_JPRB
PCAPE(KST:KEND)=0.0_JPRB
KLCL(KST:KEND)=-1
KFCL(KST:KEND)=-1
KLNB(KST:KEND)=-1

ZMINT=150._JPRB
ZMAXT=400._JPRB
ZMINDERI=1000._JPRB
ZMINQ=1.E-07_JPRB
ZMAXQ=1.0_JPRB-ZMINQ

DO JLEV=1,KLEV
  DO JLON=KST,KEND
    ZTIN(JLON,JLEV)=MAX(ZMINT,MIN(ZMAXT,PT(JLON,JLEV)))
    ZQVIN(JLON,JLEV)=MAX(ZMINQ,MIN(ZMAXQ,PQV(JLON,JLEV)))
  ENDDO
ENDDO

ZBUOYPREC(KST:KEND)=0.0_JPRB
ZT (KST:KEND)=ZTIN (KST:KEND,KLEVST)
ZQV(KST:KEND)=ZQVIN(KST:KEND,KLEVST)
ZQL(KST:KEND)=0.0_JPRB
ZQI(KST:KEND)=0.0_JPRB

IPREVIOUS_NULL_ZRT(KST:KEND)=999999

DO JLEV=KLEVST,NTCVIM+1,-1
  DO JLON=KST,KEND
    !
    !-------------------------------------------------
    ! SATURATION SPECIFIC HUMIDITY.
    !-------------------------------------------------
    !
    ZT(JLON)=MAX(ZMINT,MIN(ZMAXT,ZT(JLON)))
    ZQS(JLON)=FOQS(FOEW(ZT(JLON),MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZT(JLON))))/PRP(JLON,JLEV))  
    ZDLOG(JLON)=LOG(PRP(JLON,MIN(KLEVST,JLEV+1))/PRP(JLON,JLEV)) * MAX(0,-SIGN(1,JLEV-KLEVST))  
    !
  ENDDO
  DO JLON=KST,KEND
    !
    !-------------------------------------------------
    ! PRESSURE AND BUOYANCY.
    !-------------------------------------------------
    !
    ZQV(JLON)=MAX(ZMINQ,MIN(ZMAXQ,ZQV(JLON)))
    ZTV1=ZT(JLON)*(1.0_JPRB+ZQV(JLON)/(1.0_JPRB-ZQV(JLON))*RV/RD)/&
     & (1.0_JPRB+ZQV(JLON)/&
     & (1.0_JPRB-ZQV(JLON))+ZQL(JLON)/(1.0_JPRB-ZQL(JLON))+ZQI(JLON)/(1.0_JPRB-ZQI(JLON)))  
    ZTV2=ZTIN(JLON,JLEV)*(1.0_JPRB+ZQVIN(JLON,JLEV)/(1.0_JPRB-ZQVIN(JLON,JLEV))*RV/RD)/&
     & (1.0_JPRB+ZQVIN(JLON,JLEV)/(1.0_JPRB-ZQVIN(JLON,JLEV)) )  
    ZBUOY=(ZTV1/ZTV2-1.0_JPRB)*(RD+(RV-RD)*ZQVIN(JLON,JLEV))*ZTIN(JLON,JLEV)
    !
    !-------------------------------------------------
    ! CIN AND CAPE INTEGRALS.
    !-------------------------------------------------
    !
    ZRT(JLON)=0.5_JPRB*(ZBUOY+ZBUOYPREC(JLON))*ZDLOG(JLON)
    !------------------------------------------
    ! CUMULATE CAPE IF POSITIVE CONTRIBUTION.
    !------------------------------------------
    PCAPE(JLON)=PCAPE(JLON)+MAX(0.0_JPRB,ZRT(JLON))
    !------------------------------------------
    ! CUMULATE CIN IF NEGATIVE CONTRIBUTION AND BELOW LFC.
    !------------------------------------------
    IF(PCAPE(JLON) == 0.0_JPRB) PCIN(JLON)=PCIN(JLON)+MIN(0.0_JPRB,ZRT(JLON))
    ZBUOYPREC(JLON)=ZBUOY
    !
    !-------------------------------------------------
    ! IF THE PARCEL IS SUPERSATURATED, SUPERSATURATION IS REMOVED.
    ! THIS IS DONE THROUGH AN ISOBARIC TRANSFORMATION FROM (ZT,ZQV)
    ! TO (ZT1,ZQV1).
    !       Solution for T
    !            f(T)=cp*(T-T0)+L*(q-q0)=0
    !            with constraint q=qs(T,p0),
    !       it is solved by the method of Newton, iteration
    !       T --> T-f(T)/f'(T), with starting point ZT.
    !-------------------------------------------------
    !
    ZT1(JLON)=ZT(JLON)
  ENDDO
  DO JLON=KST,KEND
    ZQV1(JLON)=FOQS(FOEW(ZT1(JLON),MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZT1(JLON))))/PRP(JLON,JLEV))  
    IF ((ZQV1(JLON)<=ZQV(JLON)).AND.(KLCL(JLON)==-1)) THEN
      KLCL(JLON)=JLEV
    ENDIF
    IF ((PCAPE(JLON) > 0.0_JPRB).AND.(KFCL(JLON)==-1)) THEN
      KFCL(JLON)=JLEV
    ENDIF
    IF ((ZRT(JLON) <= 0.0_JPRB).AND.(JLEV<KFCL(JLON)).AND.(JLEV<IPREVIOUS_NULL_ZRT(JLON)-1)) THEN
      KLNB(JLON)=JLEV
    ENDIF
    IF (ZRT(JLON) <= 0.0_JPRB) THEN
      IPREVIOUS_NULL_ZRT(JLON)=JLEV
    ENDIF
  ENDDO

  DO JIT=1,NCAPEITER
    DO JLON=KST,KEND
      !
      !-------------------------------------------------
      ! LATENT HEAT
      !-------------------------------------------------
      !
      ZDELARG0(JLON)=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZT1(JLON)))
      ZL0(JLON)=FOLH(ZT1(JLON),ZDELARG0(JLON))
      ZFDERFOLH0(JLON)=-RV*(RGAMW+ZDELARG0(JLON)*RGAMD)
      ZCP0(JLON)=RCPD*(1.0_JPRB-ZQV1(JLON))+RCPV*ZQV1(JLON)
      !
      !-------------------------------------------------
      ! Newton's loop to solve the supersaturation.
      !-------------------------------------------------
      !
    ENDDO
    DO JLON=KST,KEND
      ZZFOEW(JLON)=FOEW(ZT1(JLON),ZDELARG0(JLON))
    ENDDO
    DO JLON=KST,KEND
      ZFDERQS0=FODQS(ZQV1(JLON),ZZFOEW(JLON)/PRP(JLON,JLEV),FODLEW(ZT1(JLON),ZDELARG0(JLON)))
      ZT1(JLON)=MAX(ZMINT,MIN(ZMAXT,ZT1(JLON)-(ZCP0(JLON)*(ZT1(JLON)-ZT(JLON))+ZL0(JLON)*(ZQV1(JLON)-ZQV(JLON)))&
       & /MAX(ZMINDERI,ZCP0(JLON)+((RCPV-RCPD)*(ZT1(JLON)-ZT(JLON))+ZL0(JLON))* ZFDERQS0&
       & +(ZQV1(JLON)-ZQV(JLON))*ZFDERFOLH0(JLON))))
      ZDELARG0(JLON)=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZT1(JLON)))
    ENDDO
    DO JLON=KST,KEND
      ZZFOEW(JLON)=FOEW(ZT1(JLON),ZDELARG0(JLON))
    ENDDO
    DO JLON=KST,KEND
      !
      !-------------------------------------------------
      ! SATURATION SPECIFIC HUMIDITY
      !-------------------------------------------------
      !
      ZQV1(JLON)=FOQS(ZZFOEW(JLON)/PRP(JLON,JLEV))  
    ENDDO
  ENDDO

  DO JLON=KST,KEND
    ZADD=GCAPERET*(ZQV(JLON)-ZQV1(JLON))
                        
    ZSQL(JLON)=ZQL(JLON)+ZADD*MAX(0.0_JPRB,SIGN(1.0_JPRB,ZT1(JLON)-RTT))
    ZSQI(JLON)=ZQI(JLON)+ZADD*MAX(0.0_JPRB,-SIGN(1.0_JPRB,ZT1(JLON)-RTT))
    !
    !-------------------------------------------------
    ! MOIST ADIABATIC ASCENT.
    ! TRANSFORMATION FROM (ZT1,ZQV1) TO (ZT2,ZQV2).
    !      Solution for T
    !             f(T)=cp*(T-T0)+L*(q-q0)+phi-phi0=0
    !      either f(T)=cp*(T-T0)+L*(q-q0)-R*T*log(p/p0)=0
    !            with constraint q=qs(T,p), knowing that q0=qs(T0,p0)
    !       it is solved by the method of Newton, iteration
    !       T --> T-f(T)/f'(T), with starting point ZT1.
    !-------------------------------------------------
    !
    ZZT(JLON)=ZT1(JLON)
  ENDDO
  !-------------------------------------------------
  ! CALCULATION IS DIVIDED IN MORE STEPS IN ORDER
  ! TO GET MORE PRECISE VALUES
  !-------------------------------------------------
  DO JETAPES=1,NETAPES
    DO JLON=KST,KEND
      ZPDEPART(JLON) =PRP(JLON,JLEV)+(PRP(JLON,JLEV-1)-PRP(JLON,JLEV))*REAL(JETAPES-1)/REAL(NETAPES)  
      ZPARRIVEE(JLON)=PRP(JLON,JLEV)+(PRP(JLON,JLEV-1)-PRP(JLON,JLEV))*REAL(JETAPES)/REAL(NETAPES)  
      ZTDEPART(JLON)=ZZT(JLON)
    ENDDO
    DO JLON=KST,KEND
      ZQVDEPART(JLON)=FOQS(FOEW(ZTDEPART(JLON),MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZTDEPART(JLON))))&
       & /ZPDEPART(JLON))
      ZLOG(JLON)=LOG(ZPARRIVEE(JLON)/ZPDEPART(JLON))
    ENDDO
    DO JIT=1,NCAPEITER
      DO JLON=KST,KEND
        ZZFOEW(JLON)=FOEW(ZZT(JLON),MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZZT(JLON))))
      ENDDO
      DO JLON=KST,KEND
        !
        !-------------------------------------------------
        ! SATURATION SPECIFIC HUMIDITY
        !-------------------------------------------------
        !
        ZZQV=FOQS(ZZFOEW(JLON)/ZPARRIVEE(JLON))  
        !
        !-------------------------------------------------
        ! LATENT HEAT
        !-------------------------------------------------
        !
        ZDELARG=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZTDEPART(JLON)))
        !
        ! CALCULATION WHERE IT IS CONSIDERED THAT LATENT HEAT RELEASE
        ! FROM CONDENSATION IS ABSORBTED ONLY BY THE GASEOUS PORTION 
        ! OF THE PARCEL, AND NOT BY THE CONDENSATE
        !
        ZL=FOLH(ZZT(JLON),ZDELARG)  
        ZDERL=-RV*(RGAMW+ZDELARG*RGAMD)
        !
        !-------------------------------------------------
        ! Newton's loop to solve the moist adiabatic ascent.
        !-------------------------------------------------
        !
        ZRDLOG=(RD+(RV-RD)*ZZQV)*ZLOG(JLON)
        ZCP=RCPD*(1.0_JPRB-ZZQV)+RCPV*ZZQV 
        ZFDERQS=FODQS(ZZQV,ZZFOEW(JLON)/ZPARRIVEE(JLON)&
          & ,FODLEW(ZZT(JLON),MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZZT(JLON)))))
        ZZT(JLON)=MAX(ZMINT,MIN(ZMAXT,ZZT(JLON)-(ZCP*(ZZT(JLON)-ZTDEPART(JLON))&
         & +ZL*(ZZQV-ZQVDEPART(JLON))-ZZT(JLON)&
         & *ZRDLOG)/MAX(ZMINDERI,ZCP+((RCPV-RCPD)*(ZZT(JLON)-ZTDEPART(JLON))&
         & +ZL-ZLOG(JLON)*ZZT(JLON)*(RV-RD))*ZFDERQS&
         & +(ZZQV-ZQVDEPART(JLON))*ZDERL-ZRDLOG)))
      ENDDO      ! JLON
    ENDDO     ! JIT
  ENDDO   !JETAPES
  
  DO JLON=KST,KEND
    ZZFOEW(JLON)=FOEW(ZZT(JLON),MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZZT(JLON))))
    !
    !-------------------------------------------------
    ! DRY ADIABATIC ASCENT.
    !-------------------------------------------------
    !
    ZDT(JLON)=ZT(JLON)*(PRP(JLON,JLEV-1)/PRP(JLON,JLEV))**((RD+(RV-RD)&
     & *ZQV(JLON))/(RCPD*(1.0_JPRB-ZQV(JLON))+RCPV*ZQV(JLON)))  
    !
  ENDDO
  DO JLON=KST,KEND
    ZT2=ZZT(JLON)

    ZQV2=FOQS(ZZFOEW(JLON)/PRP(JLON,JLEV-1))

    ZADD=GCAPERET*(ZQV1(JLON)-ZQV2)*MAX(0.0_JPRB,SIGN(1.0_JPRB,ZQV1(JLON)-ZQV2))
    ZSQL(JLON)=ZSQL(JLON)+ZADD*MAX(0.0_JPRB,SIGN(1.0_JPRB,ZT2-RTT))
    ZSQI(JLON)=ZSQI(JLON)+ZADD*MAX(0.0_JPRB,-SIGN(1.0_JPRB,ZT2-RTT))
    !
    !-------------------------------------------------
    ! UPDATE PARCEL STATE. The result of moist adiabatic ascent.
    !-------------------------------------------------
    !
    ZST=ZT2
    ZSQV(JLON)=ZQV2
    !-------------------------------------------------
    !  VALUES FROM DRY ADIABATIC OR MOIST ADIABATIC ASCENT ARE CHOSEN 
    !-------------------------------------------------
    !
    ZADD=MAX(0.0_JPRB,SIGN(1.0_JPRB,ZQV(JLON)-ZQS(JLON)))
    ZT(JLON)=ZST*ZADD + ZDT(JLON)*(1.0_JPRB-ZADD)
    ZQV(JLON)=ZSQV(JLON)*ZADD + ZQV(JLON)*(1.0_JPRB-ZADD)
    ZQL(JLON)=ZSQL(JLON)*ZADD
    ZQI(JLON)=ZSQI(JLON)*ZADD
    !-----------------------------------------------------
    ! Entrainement
    !-----------------------------------------------------
    IF (PENTRA > 0.0_JPRB) THEN
      ZDZ=(PRP(JLON,JLEV)-PRP(JLON,JLEV-1))/(0.5*(PRP(JLON,JLEV)+&
        & PRP(JLON,JLEV-1)))*(RD+(RV-RD)*(0.5*(ZQVIN(JLON,JLEV)+&
        & (ZQVIN(JLON,JLEV-1)))))*0.5*(ZTIN(JLON,JLEV)+(ZTIN(JLON,JLEV-1)))/RG
      IF(ZT(JLON) > ZTIN(JLON,JLEV-1)) THEN
        ZT(JLON)=MAX(ZTIN(JLON,JLEV-1),ZT(JLON)+PENTRA*ZDZ*(ZTIN(JLON,JLEV-1)-ZT(JLON)))
      ELSE
        ZT(JLON)=MIN(ZTIN(JLON,JLEV-1),ZT(JLON)+PENTRA*ZDZ*(ZTIN(JLON,JLEV-1)-ZT(JLON)))
      ENDIF
      IF(ZQV(JLON) > ZQVIN(JLON,JLEV-1)) THEN
        ZQV(JLON)=MAX(ZQVIN(JLON,JLEV-1),ZQV(JLON)+PENTRA*ZDZ*(ZQVIN(JLON,JLEV-1)-ZQV(JLON)))
      ELSE
        ZQV(JLON)=MIN(ZQVIN(JLON,JLEV-1),ZQV(JLON)+PENTRA*ZDZ*(ZQVIN(JLON,JLEV-1)-ZQV(JLON)))
      ENDIF
    ENDIF
  ENDDO                     ! JLON
ENDDO                     !JLEV

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FPCINCAPE',1,ZHOOK_HANDLE)
END SUBROUTINE FPCINCAPE
