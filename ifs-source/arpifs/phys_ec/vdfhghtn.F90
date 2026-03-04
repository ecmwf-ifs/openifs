! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFHGHTN (YDEPHLI, YDECLDP  , YDECUMF , YDVDF, YDSPP_CONFIG, &
                   & KIDIA   , KFDIA   , KLON    , KLEV    , PTMST, LDLAND,&
                   & PTM1    , PQM1    , PSLGSM1 , PUM1    , PVM1   ,  PCPTGZ,&
                   & PAPHM1  , PAPM1   , PGEOM1  , PGEOH   , PBLH  ,&
                   & PKMFL   , PKHFL   , PKQFL   , PGP2DSPP,PKHVFL,  PMFLX,&
                   & PSLGUH  , PQTUH   , PUUH    , PVUH    ,&
                   & PZINV   , KINV    , KCLDBASE, KCLDTOP ,PZCLDBASE, PRI,  PEIS   , KPBLTYPE)  
!     ------------------------------------------------------------------

!**   *VDFHGHTN* - DETERMINES THE PBL-HEIGHT AND STRONG UPDRAFT FIELDS
!                  USING A ENTRAINING PARCEL ASCENT METHOD.

!     A.P. SIEBESMA    30/06/99   Original (dry)
!     M. Ko"hler        3/12/2004 Moist Version
!     P. Lopez         02/06/2005 Removed useless option LPHYLIN
!     P. Bechtold, I Sandu, N Semane 01/2019: complete revision=as shallow ascent
!     P. Bechtold      26/03/2019 Add updraught momentum

!     PURPOSE
!     -------

!     DETERMINE PBL HEIGHT AND UPDRAFT FIELDS

!     INTERFACE
!     ---------

!     *VDFHGHTN* IS CALLED BY *VDFMAIN*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)        S
!     *PTM1*         TEMPERATURE AT T-1                           K
!     *PQM1*         SPECIFIC HUMUDITY AT T-1                     KG/KG
!     *PSLGSM1*      MOIST STATIC ENERGY ENTROPY EQUIV. Pascal Marquet (J)
!     *PUM1*         WIND SPEED ZONAL                             M/S
!     *PVM1*         WIND SPEED MERIDIONAL                        M/S
!     *PAPHM1*       PRESSURE AT HALF LEVEL AT T-1                PA
!     *PAPM1*        PRESSURE AT FULL LEVEL AT T-1                PA
!     *PGEOM1*       GEOPOTENTIAL AT T-1                          M2/S2
!     *PGEOH*        GEOPOTENTIAL AT HALF LEVEL                   M2/S2
!     *PBLH*         GENERAL BOUNDARY-LAYER HEIGHT                M
!     *PKMFL*        SURFACE KINEMATIC MOMENTUM FLUX              M2/S2  
!     *PKHFL*        SURFACE KINEMATIC HEAT FLUX                  K*M/S
!     *PKQFL*        SURFACE KINEMATIC MOISTURE FLUX              M/S
!     *PCPTGZ*       DRY STATIC ENERGY

!     OUTPUT PARAMETERS (REAL):

!     *PSLGUH*       UPDRAFT GENERALIZED LIQUID STATIC ENERGY (SLG)
!                    AT HALF LEVEL                                M2/S2
!     *PQTUH*        UPDRAFT SPECIFIC TOTAL WATER AT HALF LEVEL   KG/KG
!     *PZINV*        PBL HEIGHT                                   M
!     *PZCLDBASE*    CLOUD HEIGHT                                 M
!     *PMFLX*        PBL MASS FLUX                                M/S
!     *PKHVFL*       w'Thv'                                       K M/S
!     *PEIS*         INVERSION STRENGTH PARAMETER                 K

!     OUTPUT PARAMETERS (INTEGER):

!     *KINV*         Level of Inversion height
!     *KCLDBASE*     Cloud base level
!     *KCLDTOP*      Cloud top level
!     *KPBLTYPE*     0: stable PBL
!                    1: dry convective PBL (no cloud below parcel top)
!                    2: cloudy PBL ("stratocumulus")
!                    3: dry PBL below convection ("cumulus")

!     METHOD
!     ------

!     ATTENTION HERE HALF-LEVELS GO GROM 1 TO KLEV+1, THIS EASES CONSISTENCY
!     WITH CUININ AND CUBASEN

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG, RD, RCPD, RETV, RLVTT, RLSTT, RTT      
USE PARPHY   , ONLY : RKAP, REPDU2
USE YOETHF   , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 & RVTMP2, R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RTWAT, RTICE, RTICECU, &
 & RTWAT_RTICE_R, RTWAT_RTICECU_R
USE YOEPHLI  , ONLY : TEPHLI
USE YOECLDP  , ONLY : TECLDP
USE YOECUMF  , ONLY : TECUMF
USE YOEVDF   , ONLY : TVDF
USE SPP_MOD  , ONLY : TSPP_CONFIG
IMPLICIT NONE

!*         0.1    GLOBAL VARIABLES

TYPE(TEPHLI)      ,INTENT(IN)    :: YDEPHLI
TYPE(TECLDP)      ,INTENT(IN)    :: YDECLDP
TYPE(TECUMF)      ,INTENT(IN)    :: YDECUMF
TYPE(TVDF)        ,INTENT(IN)    :: YDVDF
TYPE(TSPP_CONFIG) ,INTENT(IN)    :: YDSPP_CONFIG
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
LOGICAL           ,INTENT(IN)    :: LDLAND(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLGSM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTGZ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBLH(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKMFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKHFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKQFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP2DSPP(KLON,YDSPP_CONFIG%SM%NRFTOTAL)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFLX(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKHVFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSLGUH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQTUH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUUH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVUH(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZINV(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KPBLTYPE(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZCLDBASE(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEIS(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRI(KLON,KLEV)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KINV (KLON)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KCLDBASE(KLON)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KCLDTOP(KLON)

!*         0.2    LOCAL VARIABLES

LOGICAL :: LLCUM(KLON), LLSC(KLON), LLCONV(KLON)
LOGICAL :: LLEIS=.false. ! use Woods/Bretherton EIS or Pascal Marquet entropy

INTEGER(KIND=JPIM) :: JK, JL, JKM, IKB, IKDT, IKT
INTEGER(KIND=JPIM) :: ILAB(KLON,KLEV)
INTEGER(KIND=JPIM) :: IDPL(KLON), IBOTSC(KLON), IINV(KLON), ICLDTOP2(KLON)

REAL(KIND=JPRB) ::    ZQSAT(KLON,KLEV), ZKHFL(KLON,KLEV+1), ZKQFL(KLON,KLEV+1), &
                    & ZTENH(KLON,KLEV), ZQENH(KLON,KLEV), ZQSENH(KLON,KLEV), &
                    & ZTU(KLON,KLEV), ZQU(KLON,KLEV), ZLU(KLON,KLEV), &
                    & ZTD(KLON,KLEV), ZQD(KLON,KLEV), ZWU2H(KLON,KLEV), &
                    & ZUD(KLON,KLEV), ZVD(KLON,KLEV), &
                    & ZCAPE(KLON), ZWUBASE(KLON), ZRICRI(KLON),&
                    & ZINV(KLON), ZLCL(KLON), ZDPSUM(KLON)

REAL(KIND=JPRB) ::    ZDZ   , ZZ    ,ZCPM    , ZCONS10 , ZRG, ZRGOCPD,&
                    & ZRCPD, ZMFMAX, ZRHO  ,ZENTR   , ZCM     ,&
                    & ZDMFN, ZFAC, ZDU2, ZDRORO

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "cubasen.intfb.h"
#include "cuinin.intfb.h"
#include "satur.intfb.h"
#include "vdfeis.intfb.h"

!DIR$ VFUNCTION EXPHF
#include "fcttre.func.h"

!     ------------------------------------------------------------------

!*         1.     SPECIFY CONSTANTS AND INITIALIZE VARIABLES
!                 ------------------------------------------

IF (LHOOK) CALL DR_HOOK('VDFHGHTN',0,ZHOOK_HANDLE)

ASSOCIATE(  NJKT1=>YDECUMF%NJKT1, NJKT2=>YDECUMF%NJKT2,NJKT5=>YDECUMF%NJKT5, NJKT6=>YDECUMF%NJKT6,&
 &NJKT3=>YDECUMF%NJKT3,NJKT4=>YDECUMF%NJKT4,&
 & ENTRDD=>YDECUMF%ENTRDD, ENTSTPC1=>YDECUMF%ENTSTPC1, ENTSTPC2=>YDECUMF%ENTSTPC2, &
 & RDEPTHS=>YDECUMF%RDEPTHS, RMFCMIN=>YDECUMF%RMFCMIN, REISTHSC=>YDVDF%REISTHSC )

IKDT        = KLEV         ! top level for cubasen departure test
ZCM         = 0.1_JPRB     ! prefactor of the mass flux initialization
                           ! (consistency with vdfexcu!)

ZRG         = 1.0_JPRB/RG
ZRCPD       = 1.0_JPRB/RCPD 
ZRGOCPD     = RG/RCPD

!DIR$ LOOP_INFO EST_TRIPS(16)
DO JL=KIDIA,KFDIA
  PZCLDBASE(JL)  = -100._JPRB  ! default value: no PBL cloud
  PZINV (JL)     = 0.0_JPRB    ! PBL height: stable=0
  ZLCL (JL)      = PBLH(JL)    ! LCL height
  KPBLTYPE(JL)   = 0
  KINV(JL)       = KLEV
  IINV(JL)       = KLEV
  PEIS(JL)       = -999_JPRB   ! default for inversion strength
  ZDPSUM(JL)     = 0.0_JPRB
  ZRICRI(JL)     = 50.0_JPRB
ENDDO

!---------------------------------------------------------------------

DO JK=1,KLEV
!DIR$ LOOP_INFO EST_TRIPS(16)
  DO JL=KIDIA,KFDIA
    ZKHFL(JL,JK) = 0.0_JPRB
    ZKQFL(JL,JK) = 0.0_JPRB
    ZQSAT(JL,JK) = PQM1(JL,JK)
    PRI(JL,JK) = 0.0_JPRB
  ENDDO
ENDDO
DO JK=1,KLEV-1
  DO JL=KIDIA,KFDIA
    ZDZ=PGEOM1(JL,JK)-PGEOM1(JL,JK+1)
    ZDU2=MAX(REPDU2,(PUM1(JL,JK)-PUM1(JL,JK+1))**2&
                     & +(PVM1(JL,JK)-PVM1(JL,JK+1))**2)
    ZDRORO= 2.0_JPRB * (PCPTGZ(JL,JK)-PCPTGZ(JL,JK+1))&
     & / ( PCPTGZ(JL,JK)+PCPTGZ(JL,JK+1)&
     &   - PGEOM1(JL,JK)-PGEOM1(JL,JK+1))&
     & + RETV*(PQM1(JL,JK)-PQM1(JL,JK+1))
    PRI(JL,JK)=ZDZ*ZDRORO/ZDU2
  ENDDO
ENDDO

!*    2.           INITIALIZE VALUES AT DRAUGHTS=HALF LEVELS IN 'CUINI'
!                  ----------------------------------------------------

!Note that in convection routines updraught values go from [1,KLEV]

CALL SATUR (KIDIA , KFDIA , KLON  , NJKT2 , KLEV, YDEPHLI%LPHYLIN,&
 & PAPM1   , PTM1  , ZQSAT , 1  )

!---------------------------------------------------------------------
CALL CUININ &
 & ( YDEPHLI, YDECUMF, KIDIA,    KFDIA,    KLON,    KLEV,&
 & PTM1,     PQM1,     ZQSAT,    PUM1,     PVM1,&
 & PGEOM1,   PAPHM1,&
 & ILAB,&
 & ZTENH,    ZQENH,    ZQSENH,   PGEOH,&
 & ZTU,      ZQU,      ZTD,      ZQD,&
 & PUUH(:,1:KLEV),     PVUH(:,1:KLEV),      ZUD,      ZVD,&
 & ZLU     )
!---------------------------------------------------------------------
!*    3.0          CLOUD BASE CALCULATIONS
!                  -----------------------

JK=KLEV+1
JKM=JK-1
!DIR$ LOOP_INFO EST_TRIPS(16)
DO JL=KIDIA,KFDIA
  ZRHO=PAPHM1(JL,JK)/(RD*PTM1(JL,JKM)*(1.0_JPRB+RETV*PQM1(JL,JKM)))
  ZKHFL(JL,JK) = PKHFL(JL)*RCPD*ZRHO
  ZKQFL(JL,JK) = PKQFL(JL)*ZRHO
ENDDO
! use weak entrainment for true cloud top, only used for PBLTYPE 2, store kcldtop
! Nota: 
if(.false.) then
CALL CUBASEN &
 & ( YDEPHLI, YDECLDP, YDECUMF,  YDSPP_CONFIG, KIDIA,    KFDIA,    KLON,    KLEV,   IKDT,  .false.,&
 & ZTENH,    ZQENH,    PGEOH,    PAPHM1,&
 & ZKQFL,    ZKHFL,    PGP2DSPP, PKMFL,&
 & PTM1,     PQM1,     ZQSAT,    PGEOM1, &  
 & ZTU,      ZQU,      ZLU,      ZWU2H,  ZWUBASE,&
 & ILAB,     LLCUM,    LLSC,     KCLDBASE,    IBOTSC,&
 & KCLDTOP,  IDPL,     ZCAPE )
endif
! use stronger entrainment inside cloud only for quasi-dry PBL top and erase above results, ie ZWU2H 
! for all other PBL types. Nota: draft properties for dry mas flux not affected as limited to below cloud.
CALL CUBASEN &
 & ( YDEPHLI, YDECLDP, YDECUMF,  YDSPP_CONFIG, KIDIA,    KFDIA,    KLON,    KLEV,    IKDT,  .true.,&
 & ZTENH,    ZQENH,    PGEOH,    PAPHM1,&
 & ZKQFL,    ZKHFL,    PGP2DSPP, PKMFL,&
 & PTM1,     PQM1,     ZQSAT,    PGEOM1, &  
 & ZTU,      ZQU,      ZLU,      ZWU2H,  ZWUBASE,&
 & ILAB,     LLCUM,    LLSC,     KCLDBASE,    IBOTSC,&
!& ICLDTOP2,  IDPL,     ZCAPE )
 & KCLDTOP,  IDPL,     ZCAPE )

!     -----------------------------------------------------------------

!*    4.0          RETURN MOIST CONSERVED UPDRAUGHT QUANTITIES
!                  -------------------------------------------

DO JK=1,KLEV
!DIR$ LOOP_INFO EST_TRIPb(16)
  DO JL=KIDIA,KFDIA
    ZCPM         = RCPD * ( 1.0_JPRB + RVTMP2 * ZQU(JL,JK) )
    PSLGUH(JL,JK)= ZCPM * ZTU(JL,JK) + PGEOH(JL,JK) - FOELH(ZTU(JL,JK)) * ZLU(JL,JK)
    PQTUH(JL,JK) = ZQU(JL,JK) + ZLU(JL,JK)
    PMFLX(JL,JK) = 0.0_JPRB
  ENDDO
ENDDO
JK=KLEV+1
JKM=JK-1
DO JL=KIDIA,KFDIA
  PSLGUH(JL,JK)=PSLGUH(JL,JKM)
  PQTUH(JL,JK) =PQTUH(JL,JKM)
  PMFLX(JL,JK) = 0.0_JPRB
  PUUH(Jl,JK)  = PUM1(JL,JKM)
  PVUH(Jl,JK)  = PVM1(JL,JKM)
ENDDO

!     -----------------------------------------------------------------
 
!*         5.     DETERMINE INVERSION HEIGHT  AND PBL TYPE
!                 ----------------------------------------
!   Use ouput from CUBASEN: CLDBASE and TOP level and 
!   make sure to exclude deep convection (cloud depth>200 hPa)
!   use cloud base level for shallow clouds 
!   if no cloud base use last level where kinetic energy ZWU2H positive         

!DIR$ LOOP_INFO EST_TRIPS(16)
DO JL=KIDIA,KFDIA
  LLCONV(JL)=.FALSE. 
  PKHVFL(JL)  = (PKHFL(JL) + RETV * PTM1(JL,KLEV) * PKQFL(JL)) !w'theta,v'
  IF ( PKHVFL(JL) < 0.0_JPRB ) THEN
    LLCONV(JL) = .TRUE.
  ENDIF
ENDDO

! Use if necessary inversion strength to define either a KPBLTYPE=3
! or to avoid increasing cloud dissipation in CLOUDSC for SC (which are also convective type now)
IF(LLEIS) THEN
 !Woods/Bretherton
  CALL VDFEIS(YDECUMF, KIDIA   , KFDIA   , KLON    , KLEV   , .false.,&
            & PTM1   , PQM1    , PAPM1   , PGEOM1 , ZLCL,   PEIS)
ELSE
 !Pascal Marquet MSE (no moist gradient needed)
 !Nota: REISTHSC value should be ~2 lower than for Woods/Bretherton criterion/above
  DO JL=KIDIA,KFDIA
     PEIS(JL)=MAX(PSLGSM1(JL,NJKT6)-PSLGSM1(JL,NJKT3),PSLGSM1(JL,NJKT3)-PSLGSM1(JL,KLEV))*ZRCPD
  ENDDO
ENDIF

DO JK=KLEV,2,-1
!DO JK=KLEV-1,1,-1
!DIR$ LOOP_INFO EST_TRIPS(16)
  DO JL=KIDIA,KFDIA
    IF(LLCONV(JL).AND.ZWU2H(JL,JK)>0.0_JPRB) THEN
       ! overshoot limit
      IF(PRI(JL,JK)<ZRICRI(JL)) IINV(JL)=JK
    ENDIF
  ENDDO
ENDDO



!DIR$ LOOP_INFO EST_TRIPS(16)
DO JL=KIDIA,KFDIA
   ZINV(JL)=PZINV(JL)
   IF( LLCONV(JL)) THEN
    IKB = KCLDBASE(JL)
    IKT = KCLDTOP(JL)+1
    IF(KCLDBASE(JL)==KCLDTOP(JL)) IKT=IKB
    PZINV(JL)=PGEOH(JL,IINV(JL))*ZRG
    ZINV(JL)=PZINV(JL)
    KPBLTYPE(JL)=1
    IF(IKB>0) THEN
       KPBLTYPE(JL)=3
       PZCLDBASE(JL)=PGEOH(JL,IKB)*ZRG
       ZINV(JL)=PGEOH(JL,IKB)*ZRG
       IF (PEIS(JL)>REISTHSC.AND.(PAPHM1(JL,IKB)-PAPHM1(JL,IKT)<=RDEPTHS)) THEN
         KPBLTYPE(JL)  = 2                 !Sc PBL
         IINV(JL)=IKT
         PZINV(JL)=PGEOH(JL,IINV(JL))*ZRG
       ENDIF
    ENDIF
    PZINV(JL)=MIN(PGEOH(JL,NJKT5)*ZRG,PZINV(JL))
     ZINV(JL)=MIN(PGEOH(JL,NJKT5)*ZRG, ZINV(JL))
   ENDIF
!inversion level
   KINV(JL)=MAX(NJKT5,IINV(JL))
ENDDO
!     -----------------------------------------------------------------

 
 !*          6.     DRY PARCEL MASS FLUX
  !                 --------------------


!*         6.1  mass flux initialization
!               (updraft fraction * scale factor * sigma-w at L60 * rho)
!               (ignore u* term following Anton's suggestion
!               - consistency with similarity theory close to the surface)

ZCONS10 = 3.0_JPRB/(RG*PTMST)

!DIR$ LOOP_INFO EST_TRIPS(16)
  DO JL=KIDIA,KFDIA
    ZZ = PGEOH(JL,KLEV)*ZRG
    IF ( ZINV(JL) > ZZ .AND. KPBLTYPE(JL)<=1) THEN
      ZRHO=PAPM1(JL,KLEV)/(RD*PTM1(JL,KLEV)*(1.0_JPRB+RETV*PQM1(JL,KLEV)))
      PMFLX(JL,KLEV) = ZCM * 1.2_JPRB &
       & * (  1.5_JPRB * RKAP * MAX(0.0_JPRB,-PKHVFL(JL)) * RG * ZZ / PTM1(JL,KLEV) )**0.33333_JPRB &
       & * ( 1.0_JPRB - ZZ / ZINV(JL) ) ** 0.5_JPRB * ZRHO  
      ZMFMAX = (PAPM1(JL,KLEV)-PAPM1(JL,KLEV-1)) * ZCONS10
      PMFLX(JL,KLEV)=MIN(PMFLX(JL,KLEV),ZMFMAX)
    ENDIF
  ENDDO

!*         6.2  mass flux  profile

  DO JK=KLEV-1,1,-1
!DIR$ LOOP_INFO EST_TRIPS(16)
    DO JL=KIDIA,KFDIA
      ZZ = PGEOH(JL,KLEV)*ZRG
      IF ( ZINV(JL) > ZZ .AND. PGEOH(JL,JK)*ZRG < ZINV(JL) .AND. KPBLTYPE(JL)<=1) THEN 
          ZMFMAX = (PAPM1(JL,JK)-PAPM1(JL,JK-1)) * ZCONS10
          PMFLX(JL,JK)=PMFLX(JL,KLEV)*PGEOH(JL,JK)/PGEOH(JL,KLEV)*(1.0_JPRB-PGEOH(JL,JK)*ZRG/ZINV(JL))**2       !
          PMFLX(JL,JK)=MIN(PMFLX(JL,JK),ZMFMAX)
      ! compute updraught momentum implicit, arbitrary entr/detr split
          ZDMFN=PMFLX(JL,JK)-PMFLX(JL,JK+1)
          ZFAC=SIGN(1.0_JPRB,ZDMFN)
          ZENTR=MAX(0.0_JPRB,0.9_JPRB*ZFAC)+MAX(0.0_JPRB,-0.1_JPRB*ZFAC)*0.5_JPRB
          ZENTR=ZENTR*ABS(ZDMFN/(PMFLX(JL,JK)+1.E-10_JPRB))
          PUUH(JL,JK)= (PUUH(JL,JK+1)*(1.0_JPRB-ZENTR)&
         & +2.0_JPRB*ZENTR*PUM1(JL,JK))/(1.0_JPRB+ZENTR)
          PVUH(JL,JK)= (PVUH(JL,JK+1)*(1.0_JPRB-ZENTR)&
         & +2.0_JPRB*ZENTR*PVM1(JL,JK))/(1.0_JPRB+ZENTR)
      ENDIF
    ENDDO
  ENDDO

!----------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('VDFHGHTN',1,ZHOOK_HANDLE)
END SUBROUTINE VDFHGHTN
