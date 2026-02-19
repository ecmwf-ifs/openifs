! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFHGHTNS (YDEPHLI, YDECLDP , YDECUMF2, YDVDF   ,&
                    & KIDIA  , KFDIA   , KLON    , KLEV    , PTMST,&
                    & PTM1   , PQM1    , PSLGSM1 , PUM1    , PVM1 , PCPTGZ,&
                    & PAPHM1 , PAPM1   , PGEOM1  , PGEOH   ,&
                    & PKMFL  , PKHFL   , PKQFL   , PKHVFL  ,&
                    & PZINV  , KINV    , KCLDBASE, KCLDTOP ,&
                    & PZCLDBASE, PRI, PEIS, KPBLTYPE)  
!     ------------------------------------------------------------------

!**   *VDFHGHTNS* - DETERMINES THE PBL-HEIGHT AND STRONG UPDRAFT FIELDS
!                   USING A ENTRAINING PARCEL ASCENT METHOD.

!     A.P. SIEBESMA    30/06/99   Original (dry)
!     M. Ko"hler        3/12/2004 Moist Version
!     P. Lopez         02/06/2005 Removed useless option LPHYLIN
!     P. Bechtold, I Sandu, N Semane 01/2019: complete revision=as shallow ascent
!     P. Bechtold      26/03/2019 Add updraught momentum
!     P. Lopez         26/11/2020 Simplified version for linearized model

!     PURPOSE
!     -------

!     DETERMINE PBL HEIGHT

!     INTERFACE
!     ---------

!     *VDFHGHTNS* IS CALLED BY *VDFMAINS*, *VDFMAINSTL* AND *VDFMAINSAD*

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
!     *PKMFL*        SURFACE KINEMATIC MOMENTUM FLUX              M2/S2  
!     *PKHFL*        SURFACE KINEMATIC HEAT FLUX                  K*M/S
!     *PKQFL*        SURFACE KINEMATIC MOISTURE FLUX              M/S
!     *PCPTGZ*       DRY STATIC ENERGY

!     OUTPUT PARAMETERS (REAL):

!     *PZINV*        PBL HEIGHT                                   M
!     *PZCLDBASE*    CLOUD HEIGHT                                 M
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
USE PARPHY   , ONLY : REPDU2
USE YOETHF   , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 & R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RTWAT, RTICE, RTICECU, &
 & RTWAT_RTICE_R, RTWAT_RTICECU_R
USE YOEPHLI  , ONLY : TEPHLI
USE YOECLDP  , ONLY : TECLDP
USE YOECUMF2 , ONLY : TECUMF2
USE YOEVDF   , ONLY : TVDF

IMPLICIT NONE

!*         0.1    GLOBAL VARIABLES

TYPE(TEPHLI)      ,INTENT(INOUT) :: YDEPHLI
TYPE(TECLDP)      ,INTENT(IN)    :: YDECLDP
TYPE(TECUMF2)     ,INTENT(INOUT) :: YDECUMF2
TYPE(TVDF)        ,INTENT(IN)    :: YDVDF
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKMFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKHFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKQFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKHVFL(KLON) 
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

INTEGER(KIND=JPIM) :: JK, JL, JKM, IKB, IKT, IKDT
INTEGER(KIND=JPIM) :: ILAB(KLON,KLEV)
INTEGER(KIND=JPIM) :: IBOTSC(KLON), IINV(KLON), ISTUP(KLON)

REAL(KIND=JPRB) :: ZQSAT(KLON,KLEV), ZKHFL(KLON,KLEV+1), ZKQFL(KLON,KLEV+1), &
                 & ZTENH(KLON,KLEV), ZQENH(KLON,KLEV), ZQSENH(KLON,KLEV), &
                 & ZTU(KLON,KLEV), ZQU(KLON,KLEV), ZLU(KLON,KLEV), &
                 & ZTD(KLON,KLEV), ZQD(KLON,KLEV), ZWU2H(KLON,KLEV), &
                 & ZUU(KLON,KLEV), ZVU(KLON,KLEV), &
                 & ZUD(KLON,KLEV), ZVD(KLON,KLEV), &
                 & ZCAPE(KLON), ZWUBASE(KLON), ZRICRI(KLON)
!REAL(KIND=JPRB) :: ZINV(KLON)

REAL(KIND=JPRB)    :: ZDUM(KLON,KLEV)
INTEGER(KIND=JPIM) :: IDUM(KLON,KLEV)

REAL(KIND=JPRB) :: ZDZ, ZRG, ZRCPD, ZRHO, ZDU2, ZDRORO

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "cubasens.intfb.h"
#include "cuinin2.intfb.h"
#include "satur.intfb.h"

!DIR$ VFUNCTION EXPHF
#include "fcttre.func.h"

!     ------------------------------------------------------------------

!*         1.     SPECIFY CONSTANTS AND INITIALIZE VARIABLES
!                 ------------------------------------------

IF (LHOOK) CALL DR_HOOK('VDFHGHTNS',0,ZHOOK_HANDLE)

ASSOCIATE( NJKT12=>YDECUMF2%NJKT12, NJKT22=>YDECUMF2%NJKT22, NJKT32=>YDECUMF2%NJKT32, &
         & NJKT52=>YDECUMF2%NJKT52, NJKT62=>YDECUMF2%NJKT62,&
         & RDEPTHS2=>YDECUMF2%RDEPTHS2, REISTHSC=>YDVDF%REISTHSC )

IKDT = KLEV    ! top level for cubasen departure test.

ZRG   = 1.0_JPRB / RG
ZRCPD = 1.0_JPRB / RCPD 

!DIR$ LOOP_INFO EST_TRIPS(16)
DO JL=KIDIA,KFDIA
  PZCLDBASE(JL)  = -100._JPRB  ! default value: no PBL cloud
  PZINV (JL)     = 0.0_JPRB    ! PBL height: stable=0
  KPBLTYPE(JL)   = 0
  KINV(JL)       = KLEV
  IINV(JL)       = KLEV
  PEIS(JL)       = -999_JPRB   ! default for inversion strength
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
    ZDUM(JL,JK) = 0.0_JPRB
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

CALL SATUR (KIDIA , KFDIA , KLON  , NJKT22 , KLEV, YDEPHLI%LPHYLIN,&
          & PAPM1 , PTM1  , ZQSAT , 1  )

!---------------------------------------------------------------------

CALL CUININ2 &
 & (YDECUMF2, YDEPHLI, KIDIA, KFDIA, KLON, KLEV,&
 &  .FALSE.,&
 &  PTM1,  PQM1,   ZQSAT,  PUM1,  PVM1,&
 &  ZDUM,  PGEOM1, PAPHM1,&
 &  IDUM,  ILAB,&
 &  ZTENH, ZQENH, ZQSENH, PGEOH,&
 &  ZTU,   ZQU,   ZTD,    ZQD,&
 &  ZUU,   ZVU,   ZUD,    ZVD,&
 &  ZLU)  

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

CALL CUBASENS &
 & ( YDEPHLI,  YDECLDP,  YDECUMF2, &
 &   KIDIA,    KFDIA,    KLON,    KLEV,    IKDT,  .TRUE.,&
 &   ZTENH,    ZQENH,    PGEOH,   PAPHM1,&
 &   ZKQFL,    ZKHFL,    PKMFL,&
 &   PTM1,     PQM1,     ZQSAT,   PGEOM1, &  
 &   ZTU,      ZQU,      ZLU,     ZWU2H,     ZWUBASE,&
 &   ILAB,     LLCUM,    LLSC,    KCLDBASE,  IBOTSC,&
 &   KCLDTOP,  ISTUP,    ZCAPE )
 
!     -----------------------------------------------------------------

!*         5.     DETERMINE INVERSION HEIGHT  AND PBL TYPE
!                 ----------------------------------------
!   Use ouput from CUBASEN: CLDBASE and TOP level and 
!   make sure to exclude deep convection (cloud depth>200 hPa)
!   use cloud base level for shallow clouds 
!   if no cloud base use last level where kinetic energy ZWU2H positive         

!DIR$ LOOP_INFO EST_TRIPS(16)
DO JL=KIDIA,KFDIA
  LLCONV(JL) = .FALSE. 
  PKHVFL(JL) = (PKHFL(JL) + RETV * PTM1(JL,KLEV) * PKQFL(JL))     !w'theta,v'
  IF (PKHVFL(JL) < 0.0_JPRB) THEN
    LLCONV(JL) = .TRUE.
  ENDIF
ENDDO

! Use if necessary inversion strength to define either a KPBLTYPE=3
! or to avoid increasing cloud dissipation in CLOUDSC for SC (which are also convective type now)
! Pascal Marquet MSE (no moist gradient needed)
! Nota: REISTHSC value should be ~2 lower than for Woods/Bretherton criterion/above
DO JL=KIDIA,KFDIA
  PEIS(JL) = MAX(PSLGSM1(JL,NJKT62)-PSLGSM1(JL,NJKT32),PSLGSM1(JL,NJKT32)-PSLGSM1(JL,KLEV))*ZRCPD
ENDDO

DO JK=KLEV,2,-1
!DO JK=KLEV-1,1,-1
!DIR$ LOOP_INFO EST_TRIPS(16)
  DO JL=KIDIA,KFDIA
    IF (LLCONV(JL) .AND. ZWU2H(JL,JK) > 0.0_JPRB) THEN
       ! overshoot limit
      IF (PRI(JL,JK) < ZRICRI(JL)) IINV(JL)=JK
    ELSE
      CYCLE
    ENDIF
  ENDDO
ENDDO

!DIR$ LOOP_INFO EST_TRIPS(16)
DO JL=KIDIA,KFDIA
!   ZINV(JL)=PZINV(JL)
   IF (LLCONV(JL)) THEN
    IKB = KCLDBASE(JL)
    IKT = KCLDTOP(JL)+1
    IF (KCLDBASE(JL) == KCLDTOP(JL)) IKT=IKB
    PZINV(JL)=PGEOH(JL,IINV(JL))*ZRG
!    ZINV(JL)=PZINV(JL)
    KPBLTYPE(JL)=1
    IF (IKB > 0) THEN
       KPBLTYPE(JL)=3
       PZCLDBASE(JL)=PGEOH(JL,IKB)*ZRG
!       ZINV(JL)=PGEOH(JL,IKB)*ZRG
       IF (PEIS(JL) > REISTHSC .AND. PAPHM1(JL,IKB)-PAPHM1(JL,IKT) <= RDEPTHS2) THEN
         KPBLTYPE(JL) = 2                 !Sc PBL
         IINV(JL)=IKT
         PZINV(JL)=PGEOH(JL,IINV(JL))*ZRG
       ENDIF
    ENDIF
    PZINV(JL)=MIN(PGEOH(JL,NJKT52)*ZRG,PZINV(JL))
!     ZINV(JL)=MIN(PGEOH(JL,NJKT52)*ZRG, ZINV(JL))
   ENDIF
!inversion level
   KINV(JL) = MAX(NJKT52,IINV(JL))
ENDDO

!----------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('VDFHGHTNS',1,ZHOOK_HANDLE)
END SUBROUTINE VDFHGHTNS
