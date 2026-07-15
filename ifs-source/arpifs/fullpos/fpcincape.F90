! (C) Copyright 1989- Meteo-France.

!OPTIONS XOPT(NOEVAL)
SUBROUTINE FPCINCAPE(YDTOPH,KST,KEND,KPROMA,KLEV,KLEVST,PENTRA,PMLDEP,PT,PRP,PQV,PCAPE,PCIN,KLCL,KFCL,KLNB)

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
! PMLDEP   : Mean Layer DEPth for MLCAPE computation (Pa).
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

!      CIN IS MASS SPECIFIC ENERGY TO RAISE THE PARCEL FROM 
!          from LO to LFC.
!          ONLY THE SUM OF NEGATIVE TERMS.
!      CAPE IS MASS SPECIFIC ENERGY  PROVIDED BY THE RAISE OF THE PARCEL
!          from LFC to LNB.
!          ONLY THE SUM OF POSITIVE TERMS.

! AUTHOR:   2001-03, J.M. PIRIOU, N. PRISTOV.

! MODIFICATIONS:
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!    2010-02-04  J.M. Piriou. Bug correction: derivative of qs(T,p) from the Newton loops.
!    2010-02-11  J.M. Piriou. Bug correction: JLON dimension of ZFDERQS0 array.
!    2010-02-17  J.M. Piriou. Protections in case of cold temperature, low pressure, etc.
!    2010-02-17  J.M. Piriou. First Newton loop in a single DO loop.
!    2010-02-17  J.M. Piriou. Compute CAPE only below the convective reference level (NTCVIM).
!    2018-10-11  J.M. Piriou. Mean level CAPE if KLEVST < 0. 
!    2020-07-06  J.M. Piriou. Compute entrainment and ascent in the same Newton loop.
!    2021-06-15  J.M. Piriou. Protect CIN from FA & GRIB compacting, in imposing CIN < GCINMAX.
!    R. El Khatib 08-Jul-2022 Contribution to the encapsulation of YOMCST and YOETHF
! --------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB  , JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCAPE  , ONLY :  NCAPEITER, GMISCINV, GCINMAX, LADAE, GCAPEMIN
USE YOMCST   , ONLY : YDCST=>YRCST ! allows use of included functions. REK.
USE YOMTOPH, ONLY : TTOPH
USE YOMLSFORC, ONLY : LMUSCLFA,NMUSCLFA

IMPLICIT NONE

TYPE(TTOPH)       ,INTENT(IN)    :: YDTOPH
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PENTRA
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMLDEP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRP(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQV(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAPE(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCIN(KPROMA) 
INTEGER(KIND=JPIM)   ,INTENT(OUT)   :: KLCL(KPROMA)
INTEGER(KIND=JPIM)   ,INTENT(OUT)   :: KFCL(KPROMA)
INTEGER(KIND=JPIM)   ,INTENT(OUT)   :: KLNB(KPROMA)
INTEGER(KIND=JPIM)      :: IPREVIOUS_NULL_ZRT(KPROMA)
INTEGER(KIND=JPIM) :: JLEV,JLON,JIT
REAL(KIND=JPRB) :: ZDLOG(KPROMA),ZBUOY,ZTV1,ZTV2,ZRT(KPROMA),ZDELARG,ZL,ZCP, &
 & ZFDERQS,ZDERL(KPROMA),ZRDLOG,ZZQV
REAL(KIND=JPRB) :: ZQS(KPROMA),ZQSENV(KPROMA,KLEV),ZQV(KPROMA), ZQV1(KPROMA), &
 & ZT(KPROMA), ZT1(KPROMA), &
 & ZLOG(KPROMA), &
 & ZZT(KPROMA),ZDELTA(KPROMA)
REAL(KIND=JPRB) :: ZDELARG0(KPROMA),ZL0(KPROMA)
REAL(KIND=JPRB) :: ZDZ(KPROMA,KLEV), ZZFOEW(KPROMA)
REAL(KIND=JPRB) :: ZMAXT,ZMINT,ZMINDERI,ZMINQ,ZMAXQ
REAL(KIND=JPRB) :: ZTIN(KPROMA,KLEV) 
REAL(KIND=JPRB) :: ZQVIN(KPROMA,KLEV) 
REAL(KIND=JPRB) :: ZMEAN_THETA(KPROMA),ZMEAN_QV(KPROMA),ZDENO(KPROMA)
REAL(KIND=JPRB) :: ZBIN(KPROMA,KLEV)
REAL(KIND=JPRB) :: ZBINCAPE
INTEGER(KIND=JPIM) :: ILEVST
REAL(KIND=JPRB) :: ZMIX(KPROMA),ZAUGM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZT_ASC(KPROMA,KLEV),ZQV_ASC(KPROMA,KLEV)
REAL(KIND=JPRB) :: ZTRS(KPROMA),ZQVRS(KPROMA)
REAL(KIND=JPRB) :: ZF,ZDERI,ZQVDRY,ZQVMOIST,ZENTRA,ZQSCL(KPROMA)
!---------

!  FUNCTIONS
#include "wrscmr.intfb.h"
#include "fctast.func.h"
#include "fcttrm.func.h"
#include "fcttim.func.h"

IF (LHOOK) CALL DR_HOOK('FPCINCAPE',0,ZHOOK_HANDLE)
ASSOCIATE(NTCVIM=>YDTOPH%NTCVIM, &
 & RTT=>YDCST%RTT, RDAY=>YDCST%RDAY, REPSM=>YDCST%REPSM, RETV=>YDCST%RETV, &
 & RCW=>YDCST%RCW, REA=>YDCST%REA, RCPD=>YDCST%RCPD, RCPV=>YDCST%RCPV, RCS=>YDCST%RCS, &
 & RLVTT=>YDCST%RLVTT, RLSTT=>YDCST%RLSTT, RBETS=>YDCST%RBETS, RALPW=>YDCST%RALPW, &
 & RBETW=>YDCST%RBETW, RGAMW=>YDCST%RGAMW, RALPS=>YDCST%RALPS, RGAMS=>YDCST%RGAMS, &
 & RALPD=>YDCST%RALPD, RBETD=>YDCST%RBETD, RGAMD=>YDCST%RGAMD, RG=>YDCST%RG, &
 & RATM=>YDCST%RATM, RV=>YDCST%RV, RD=>YDCST%RD)
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
!
!-------------------------------------------------
! Initialize T and qv.
!-------------------------------------------------
!
DO JLEV=1,KLEV
  DO JLON=KST,KEND
    ZTIN(JLON,JLEV)=MAX(ZMINT,MIN(ZMAXT,PT(JLON,JLEV)))
    ZQVIN(JLON,JLEV)=MAX(ZMINQ,MIN(ZMAXQ,PQV(JLON,JLEV)))
    ZQSENV(JLON,JLEV)=FOQS(FOEW(ZTIN(JLON,JLEV),MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZTIN(JLON,JLEV))))/PRP(JLON,JLEV))  
    ZT_ASC(JLON,JLEV)=ZTIN(JLON,JLEV)
    ZQV_ASC(JLON,JLEV)=ZQVIN(JLON,JLEV)
  ENDDO
ENDDO
DO JLON=KST,KEND
  ! ZQSCL is saturation qv at condensation level (CL). 
  ! Below this CL it is initialized to a low value
  ! so that the cubic entrainment function saturates to 1.
  ZQSCL(JLON)=FOQS(FOEW(ZTIN(JLON,NTCVIM),MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZTIN(JLON,NTCVIM))))/PRP(JLON,NTCVIM))  
ENDDO

!
!-------------------------------------------------
! Starting level and starting parcel.
!-------------------------------------------------
!
ILEVST=ABS(KLEVST)
ZT (KST:KEND)=ZTIN (KST:KEND,ILEVST)
ZQV(KST:KEND)=ZQVIN(KST:KEND,ILEVST)

IPREVIOUS_NULL_ZRT(KST:KEND)=999999

ZBIN(:,:)=0._JPRB ! 1. inside a PMLDEP pressure-depth layer, 0. outside.
IF(KLEVST < 0) THEN
  ! If KLEVST < 0, a mean layer CAPE is computed: 
  ! ascent starts from a (T,qv) value, valid at level -KLEVST, obtained as the
  ! mean value of theta and qv over a layer, centered on -KLEVST,
  ! whose depth is PMLDEP (in Pa).
  !
  ! Compute mean theta and qv over the layers.
  ZMEAN_THETA(:)=0._JPRB
  ZMEAN_QV(:)=0._JPRB
  ZDENO(:)=0._JPRB
  DO JLEV=KLEV,NTCVIM+1,-1
    DO JLON=KST,KEND
      ZBIN(JLON,JLEV)=MAX(0._JPRB,SIGN(1._JPRB,PMLDEP+PRP(JLON,JLEV)-PRP(JLON,ILEVST)))&
       & *MAX(0._JPRB,SIGN(1._JPRB,PRP(JLON,ILEVST)-PRP(JLON,JLEV)))
      ZMEAN_THETA(JLON)=ZMEAN_THETA(JLON)+ZBIN(JLON,JLEV)*ZTIN(JLON,JLEV)*(RATM/PRP(JLON,JLEV))**(RD/RCPD)
      ZMEAN_QV(JLON)=ZMEAN_QV(JLON)+ZBIN(JLON,JLEV)*ZQVIN(JLON,JLEV)
      ZDENO(JLON)=ZDENO(JLON)+ZBIN(JLON,JLEV)
    ENDDO
  ENDDO
  DO JLON=KST,KEND
    ZMEAN_THETA(JLON)=ZMEAN_THETA(JLON)/MAX(1._JPRB,ZDENO(JLON))
    ZMEAN_QV(JLON)=ZMEAN_QV(JLON)/MAX(1._JPRB,ZDENO(JLON))
    !
    ! Starting T and qv values.
    ZT(JLON)=ZMEAN_THETA(JLON)*(RATM/PRP(JLON,ILEVST))**(-RD/RCPD)
    ZQV(JLON)=ZMEAN_QV(JLON)
  ENDDO
ENDIF

DO JLEV=ILEVST,NTCVIM+1,-1
  DO JLON=KST,KEND
    !
    !-------------------------------------------------
    ! SATURATION SPECIFIC HUMIDITY ZQS TO DIAGNOSE WHETHER CONDENSATION IS REACHED OR NOT.
    !-------------------------------------------------
    !
    ZT(JLON)=MAX(ZMINT,MIN(ZMAXT,ZT(JLON)))
    ZQV(JLON)=MAX(ZMINQ,MIN(ZMAXQ,ZQV(JLON)))
    ZQS(JLON)=FOQS(FOEW(ZT(JLON),MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZT(JLON))))/PRP(JLON,JLEV))  
    ZDELTA(JLON)=MAX(0._JPRB,SIGN(1._JPRB,ZQV(JLON)-ZQS(JLON)))
    ZDLOG(JLON)=LOG(PRP(JLON,MIN(ILEVST,JLEV+1))/PRP(JLON,JLEV)) * MAX(0,-SIGN(1,JLEV-ILEVST))  
    !
  ENDDO
  DO JLON=KST,KEND
    !
    !-------------------------------------------------
    ! IF THE PARCEL IS SUPERSATURATED, SUPERSATURATION IS REMOVED.
    ! THIS IS DONE THROUGH AN ISOBARIC TRANSFORMATION FROM (ZT,ZQV)
    ! TO (ZTRS,ZQVRS): RS=Remove Sursaturation.
    ! Solution for T
    ! f(T)=cp*(T-T0)+L*(q-q0)=0
    ! with constraint q=qs(T,p0),
    ! solved by the method of Newton, iteration T --> T-f(T)/f'(T), with starting point ZT.
    !-------------------------------------------------
    !
    ZTRS(JLON)=ZT(JLON)
    !
    !-------------------------------------------------
    ! LATENT HEAT
    !-------------------------------------------------
    !
    ZDELARG0(JLON)=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZT(JLON)))
    ZL0(JLON)=FOLH(ZT(JLON),ZDELARG0(JLON))
    ZDERL(JLON)=-RV*(RGAMW+ZDELARG0(JLON)*RGAMD)
  ENDDO ! JLON
  !
  !-------------------------------------------------
  ! Newton's loop to solve the sursaturation.
  !-------------------------------------------------
  !
  DO JIT=1,NCAPEITER
    DO JLON=KST,KEND
      ZZFOEW(JLON)=FOEW(ZTRS(JLON),MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZTRS(JLON))))
    ENDDO
    DO JLON=KST,KEND
      !
      ! SATURATION SPECIFIC HUMIDITY ZQVRS.
      ZQVRS(JLON)=FOQS(ZZFOEW(JLON)/PRP(JLON,JLEV))  
      ZCP=RCPD*(1.0_JPRB-ZQVRS(JLON))+RCPV*ZQVRS(JLON) 
      ZFDERQS=FODQS(ZQVRS(JLON),ZZFOEW(JLON)/PRP(JLON,JLEV)&
        & ,FODLEW(ZTRS(JLON),MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZTRS(JLON)))))
      ZF=ZCP*(ZTRS(JLON)-ZT(JLON))+ZL0(JLON)*(ZQVRS(JLON)-ZQV(JLON))
      ZDERI=ZCP+((RCPV-RCPD)*(ZTRS(JLON)-ZT(JLON))+ZL0(JLON))*ZFDERQS+ZDERL(JLON)*(ZQVRS(JLON)-ZQV(JLON))
      ZTRS(JLON)=MAX(ZMINT,MIN(ZMAXT,ZTRS(JLON)-ZF/MAX(ZMINDERI,ZDERI)))
    ENDDO ! JLON
  ENDDO ! JIT
  DO JLON=KST,KEND
    !
    ! Choose RS (Remove Saturation) point or original point, depending on saturation ZDELTA.
    ZT(JLON)=ZDELTA(JLON)*ZTRS(JLON)+(1._JPRB-ZDELTA(JLON))*ZT(JLON)
    ZQV(JLON)=ZDELTA(JLON)*ZQVRS(JLON)+(1._JPRB-ZDELTA(JLON))*ZQV(JLON)
    !
    ! ZT_ASC and ZQV_ASC are stored as profile arrays for MUSC diagnostics only.
    ZT_ASC(JLON,JLEV)=ZT(JLON)
    ZQV_ASC(JLON,JLEV)=ZQV(JLON)
    !
    !-------------------------------------------------
    ! BUOYANCY.
    !-------------------------------------------------
    !
    ZTV1=ZT(JLON)*(1._JPRB+RETV*ZQV(JLON)) ! Tv ascent.
    ZTV2=ZTIN(JLON,JLEV)*(1._JPRB+RETV*ZQVIN(JLON,JLEV)) ! Tv environment.
    ZBUOY=RG*(ZTV1/ZTV2-1._JPRB)
    !
    !-------------------------------------------------
    ! CIN AND CAPE INTEGRALS.
    !-------------------------------------------------
    !
    ZRT(JLON)=ZBUOY/RG*(RD+(RV-RD)*ZQVIN(JLON,JLEV))*ZTIN(JLON,JLEV)*ZDLOG(JLON)
    !------------------------------------------
    ! CUMULATE CAPE IF POSITIVE CONTRIBUTION AND SATURATION.
    !------------------------------------------
    PCAPE(JLON)=PCAPE(JLON)+MAX(0.0_JPRB,ZRT(JLON))*MAX(0._JPRB,SIGN(1._JPRB,ZDELTA(JLON)-0.5_JPRB))
    !------------------------------------------
    ! CUMULATE CIN IF NEGATIVE CONTRIBUTION AND BELOW LFC.
    !------------------------------------------
    IF(PCAPE(JLON) < GCAPEMIN) PCIN(JLON)=PCIN(JLON)+MIN(0.0_JPRB,ZRT(JLON))
  ENDDO

  DO JLON=KST,KEND
    IF (ZDELTA(JLON) > 0.5_JPRB .AND. KLCL(JLON)==-1) THEN
      ! Parcel is saturated. LCL found.
      KLCL(JLON)=JLEV
      ZQSCL(JLON)=FOQS(FOEW(ZTIN(JLON,JLEV),MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZTIN(JLON,JLEV))))/PRP(JLON,JLEV))  
    ENDIF
    IF (PCAPE(JLON) > 0.0_JPRB .AND. KFCL(JLON)==-1) THEN
      ! Positive CAPE. FCL found.
      KFCL(JLON)=JLEV
    ENDIF
    IF (ZRT(JLON) <= 0.0_JPRB .AND. JLEV<KFCL(JLON) .AND. JLEV<IPREVIOUS_NULL_ZRT(JLON)-1) THEN
      ! Negative buoyancy. LNB found.
      KLNB(JLON)=JLEV
    ENDIF
    IF (ZRT(JLON) <= 0.0_JPRB) THEN
      IPREVIOUS_NULL_ZRT(JLON)=JLEV
    ENDIF
    !
    !-------------------------------------------------
    ! MOIST OR DRY ASCENT FROM JLEV TO JLEV-1.
    ! TRANSFORMATION FROM (ZT1,ZQV1) TO (ZZT,ZZQV).
    !      Solution for T
    !             f(T)=cp*(T-T0)+delta*L*(q-q0)+phi-phi0=0
    !      either f(T)=cp*(T-T0)+delta*L*(q-q0)-R*T*log(p/p0)=0
    !            with constraint q=qs(T,p), knowing that q0=qs(T0,p0)
    !       it is solved by the method of Newton, iteration
    !       T --> T-f(T)/f'(T), with starting point ZT1.
    !-------------------------------------------------
    !
    ZT1(JLON)=ZT(JLON)
    ZQV1(JLON)=ZQV(JLON)
    ZZT(JLON)=ZT1(JLON)
    ZLOG(JLON)=LOG(PRP(JLON,JLEV-1)/PRP(JLON,JLEV))
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
      ZZQV=FOQS(ZZFOEW(JLON)/PRP(JLON,JLEV-1))  
      !
      !-------------------------------------------------
      ! LATENT HEAT
      !-------------------------------------------------
      !
      ZDELARG=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZT1(JLON)))
      ZL=FOLH(ZZT(JLON),ZDELARG)  
      ZDERL(JLON)=-RV*(RGAMW+ZDELARG*RGAMD)
      IF(LADAE) THEN
        !-------------------------------------------------
        ! Entrainment is larger in low relative humidity and large saturation specific humidity.
        ! This entrainment formula follows the Tiedtke-Bechtold convection scheme.
        !-------------------------------------------------
        ZENTRA=(1._JPRB-ZBIN(JLON,JLEV))*PENTRA*(1.3_JPRB-MIN(1._JPRB,ZQVIN(JLON,JLEV)/ZQSENV(JLON,JLEV)))&
          & *MIN(1._JPRB,ZQSENV(JLON,JLEV)/ZQSCL(JLON))**3
      ELSE
        !-------------------------------------------------
        ! Uniform entrainment.
        !-------------------------------------------------
        ZENTRA=(1._JPRB-ZBIN(JLON,JLEV))*PENTRA
      ENDIF
      !
      !-------------------------------------------------
      ! Newton's loop to solve the moist or dry adiabatic ascent and entrainment.
      !-------------------------------------------------
      !
      ZRDLOG=(RD+(RV-RD)*ZZQV)*ZLOG(JLON)
      ZCP=RCPD*(1.0_JPRB-ZZQV)+RCPV*ZZQV 
      ZFDERQS=FODQS(ZZQV,ZZFOEW(JLON)/PRP(JLON,JLEV-1)&
        & ,FODLEW(ZZT(JLON),MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZZT(JLON)))))
      ! ZDZ: layer depth in meter.
      ZDZ(JLON,JLEV)=(PRP(JLON,JLEV)-PRP(JLON,JLEV-1))/(0.5*(PRP(JLON,JLEV)+ &
        & PRP(JLON,JLEV-1)))*(RD+(RV-RD)*(0.5*(ZQVIN(JLON,JLEV)+ &
        & (ZQVIN(JLON,JLEV-1)))))*0.5*(ZTIN(JLON,JLEV)+(ZTIN(JLON,JLEV-1)))/RG
      ZMIX(JLON)=MIN(1._JPRB,ZENTRA*ZDZ(JLON,JLEV))
      ZAUGM=1._JPRB+ZMIX(JLON)
      ZF=ZCP*(ZZT(JLON)*ZAUGM-ZT1(JLON)-ZMIX(JLON)*PT(JLON,JLEV))&
       & +ZDELTA(JLON)*ZL*(ZZQV*ZAUGM-ZQV1(JLON)-ZMIX(JLON)*PQV(JLON,JLEV))-ZZT(JLON)*ZRDLOG
      ZDERI=ZCP*ZAUGM+ZDELTA(JLON)*((RCPV-RCPD)*(ZZT(JLON)*ZAUGM-ZT1(JLON)-ZMIX(JLON)*PT(JLON,JLEV))&
       & +ZL*ZAUGM-(RV-RD)*ZZT(JLON)*ZLOG(JLON))*ZFDERQS&
       & +(ZZQV*ZAUGM-ZQV1(JLON)-ZMIX(JLON)*PQV(JLON,JLEV))*ZDELTA(JLON)*ZDERL(JLON)-ZRDLOG
      ZZT(JLON)=MAX(ZMINT,MIN(ZMAXT,ZZT(JLON)-ZF/MAX(ZMINDERI,ZDERI)))
    ENDDO ! JLON
  ENDDO ! JIT
  
  DO JLON=KST,KEND
    !
    !-------------------------------------------------
    ! UPDATE PARCEL STATE. T in moist or dry mode.
    !-------------------------------------------------
    !
    ZT(JLON)=ZZT(JLON)
    !
    !-------------------------------------------------
    ! UPDATE PARCEL STATE. qv in moist and dry modes.
    !-------------------------------------------------
    !
    ZQVMOIST=FOQS(FOEW(ZT(JLON),MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZT(JLON))))/PRP(JLON,JLEV-1))
    ZQVDRY=ZQV1(JLON)+ZMIX(JLON)*(ZQVIN(JLON,JLEV)-ZQV1(JLON))
    !
    !-------------------------------------------------
    !  VALUES FROM DRY OR MOIST ASCENT MODES ARE CHOSEN 
    !-------------------------------------------------
    !
    ZQV(JLON)=ZDELTA(JLON)*ZQVMOIST+(1._JPRB-ZDELTA(JLON))*ZQVDRY
  ENDDO ! JLON
ENDDO ! JLEV

DO JLON=KST,KEND
  ZBINCAPE=MAX(0._JPRB,SIGN(1._JPRB,PCAPE(JLON)-GCAPEMIN))
  ! 1. If ZBINCAPE=0., no LFC (Level of Free Convection) has been found, set CIN to a positive value GMISCINV, 
  ! to plot it as missing data.
  ! 2. If ZBINCAPE=1. : to avoid compacting problems that may create positive CIN values, one forces CIN <= GCINMAX.
  ! GCINMAX is < 0., its absolute value is choosen high enough to ensure 
  ! that CIN remains < 0 even after GRIB or FA file compacting.
  PCIN(JLON)=ZBINCAPE*MIN(GCINMAX,PCIN(JLON))+(1._JPRB-ZBINCAPE)*GMISCINV
ENDDO

IF(LMUSCLFA) THEN
  CALL WRSCMR(NMUSCLFA,'T_ASC',ZT_ASC,KPROMA,KLEV)
  CALL WRSCMR(NMUSCLFA,'QV_ASC',ZQV_ASC,KPROMA,KLEV)
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FPCINCAPE',1,ZHOOK_HANDLE)
END SUBROUTINE FPCINCAPE
