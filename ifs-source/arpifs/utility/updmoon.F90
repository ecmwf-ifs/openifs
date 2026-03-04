! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE UPDMOON(YDRIP)

!     ------------------------------------------------------------------
!**** UPDMOON* - Update moon position.
!     Calcul de la position lunaire (longitude et latitude)
!     et de l'eclairement lunaire au sommet de l'atmosphere.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *UPDMOON(...)*

!        Explicit arguments :
!        --------------------

!  En entree : RTIMTR: date courante
!                      (nombre de secondes ecoulees depuis le 1.1.2000 a 12 UTC)
!  En sortie : RTMOLT: pi - (longitude de la Lune) (radians)
!              RDECLU: latitude de la Lune (radians)
!              RIP0LU: eclairement lunaire hors atmosphere

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation.

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!      Beatrice Fradon *SCEM/PREVI/CDMA*
!      Original : 1999-04-19

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMCST   , ONLY : RPI, RA, RDAY
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMRIP   , ONLY : TRIP

!     ------------------------------------------------------------------

IMPLICIT NONE

!                 ZPHASE        reel      fraction eclairee de la Lune (0 a 1)
!                 ZDISMOON      reel      distance Terre - Lune (km)
!                 ZCORR         reel      facteur pour corriger l'eclairement
!                                         lunaire de la variation de distance
!                                         Terre - Lune
!                 ZRAS          reel      ascension droite du Soleil
!                 ZDECS         reel      declinaison du Soleil
!                 ZDISS         reel      distance Terre-Soleil
!                 ZRAM          reel      ascension droite de la Lune
!                 ZDECM         reel      declinaison de la Lune
!                 ZPHI          reel      elongation geocentrique de la Lune
!                                         par rapport au Soleil
!                 ZINC          reel      elongation de la Terre par rapport
!                                         au Soleil, vue de la Lune
!                 ZGST         reel       temps sideral local de Greenwich
!                 ZT,ZTT,ZT0,  reels      variables intermediaires pour le
!                 ZL,ZG,ZV,ZU,            calcul de la position de la Lune
!                 ZW,ZS,ZM,ZF,            et du Soleil
!                 ZD,ZN,ZC(4)
!                 ZECLMOON_MOY reel      eclairement moyen du a la Lune au
!                                         sommet de l'atmosphere (W/m2)
!                 ZDISMOON_MOY reel      distance moyenne Terre-Lune (km)
!                 ZDISSUN_MOY  reel      distance Terre-Soleil moyenne (km)

TYPE(TRIP),INTENT(INOUT):: YDRIP
REAL(KIND=JPRB) :: ZPHASE,ZDISMOON,ZCORR,ZRAS,ZDECS,ZDISS,ZRAM,ZDECM
REAL(KIND=JPRB) :: ZECLMOON_MOY,ZDISMOON_MOY
REAL(KIND=JPRB) :: Z2PI,ZDISSUN_MOY
REAL(KIND=JPRB) :: ZT,ZTT,ZT0,ZL,ZG,ZV,ZU,ZW,ZS,ZM,ZF,ZD,ZN,ZPHI,ZINC,ZC(4),ZGST,ZTMOLT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('UPDMOON',0,ZHOOK_HANDLE)
ASSOCIATE(RDECLU=>YDRIP%RDECLU, RIP0LU=>YDRIP%RIP0LU, RTIMTR=>YDRIP%RTIMTR, &
 & RTMOLT=>YDRIP%RTMOLT)
!     ------------------------------------------------------------------

!     INITIALISATIONS
!     ---------------

ZECLMOON_MOY=4.25158E-3_JPRB  ! UTILISATION DE L'ALBEDO DE LANE ET IRVINE(1973)
ZDISMOON_MOY=384402._JPRB
Z2PI = 2.0_JPRB*RPI
ZDISSUN_MOY = 149597870.61_JPRB

!     CALCULS DE DATE
!     ---------------

ZT = RTIMTR / RDAY ! nombre de jours ecoules depuis le 1.1.2000 a 12 UTC.
ZTT = ZT/36525._JPRB + 1.0_JPRB   ! NOMBRE DE SIECLES JULIENS DEPUIS 1900

!     CALCUL DE LA POSITION SOLAIRE (VAN FLANDERN ET PULKKINEN, 1979)
!     ---------------------------------------------------------------

ZL = 0.779072_JPRB + 0.00273790931_JPRB * ZT
ZL = ZL - INT(ZL)
ZL = ZL*Z2PI
ZG = 0.993126_JPRB + 0.0027377785_JPRB * ZT
ZG = ZG - INT(ZG)
ZG = ZG*Z2PI
ZV = 0.39785_JPRB * SIN(ZL)
ZV = ZV - 0.01_JPRB * SIN(ZL-ZG)
ZV = ZV + 0.00333_JPRB * SIN(ZL+ZG)
ZV = ZV - 0.00021_JPRB * ZTT * SIN(ZL)
ZU = 1.0_JPRB - 0.03349_JPRB * COS(ZG)
ZU = ZU - 0.00014_JPRB * COS(2.0_JPRB*ZL)
ZU = ZU + 0.00008_JPRB * COS(ZL)
ZW = -0.0001_JPRB - 0.04129_JPRB * SIN(2.0_JPRB*ZL)
ZW = ZW + 0.03211_JPRB * SIN(ZG)
ZW = ZW + 0.00104_JPRB * SIN(2.0_JPRB*ZL-ZG)
ZW = ZW - 0.00035_JPRB * SIN(2.0_JPRB*ZL+ZG)
ZW = ZW - 0.00008_JPRB * ZTT * SIN(ZG)

ZS = ZW / (ZU-ZV*ZV)**0.5_JPRB
ZRAS = ZL + ATAN(ZS/(1.0_JPRB-ZS*ZS)**0.5_JPRB)
ZS = ZV / ZU**0.5_JPRB
ZDECS = ATAN(ZS/(1.0_JPRB-ZS*ZS)**0.5_JPRB)
ZDISS = 1.00021_JPRB * ZU**0.5_JPRB * ZDISSUN_MOY

!     CALCUL DE LA POSITION LUNAIRE
!     -----------------------------

ZL = 0.606434_JPRB + 0.03660110129_JPRB * ZT
ZM = 0.374897_JPRB + 0.03629164709_JPRB * ZT
ZF = 0.259091_JPRB + 0.03674819520_JPRB * ZT
ZD = 0.827362_JPRB + 0.03386319198_JPRB * ZT
ZN = 0.347343_JPRB - 0.00014709391_JPRB * ZT
ZG = 0.993126_JPRB + 0.00273777850_JPRB * ZT
ZL = ZL - INT(ZL)
ZM = ZM - INT(ZM)
ZF = ZF - INT(ZF)
ZD = ZD - INT(ZD)
ZN = ZN - INT(ZN)
ZG = ZG - INT(ZG)
ZL = ZL * Z2PI
ZM = ZM * Z2PI
ZF = ZF * Z2PI
ZD = ZD * Z2PI
ZN = ZN * Z2PI
ZG = ZG * Z2PI
ZV = 0.39558_JPRB * SIN(ZF+ZN)
ZV = ZV + 0.08200_JPRB * SIN(ZF)
ZV = ZV + 0.03257_JPRB * SIN(ZM-ZF-ZN)
ZV = ZV + 0.01092_JPRB * SIN(ZM+ZF+ZN)
ZV = ZV + 0.00666_JPRB * SIN(ZM-ZF)
ZV = ZV - 0.00644_JPRB * SIN(ZM+ZF-2.0_JPRB*ZD+ZN)
ZV = ZV - 0.00331_JPRB * SIN(ZF-2.0_JPRB*ZD+ZN)
ZV = ZV - 0.00304_JPRB * SIN(ZF-2.0_JPRB*ZD)
ZV = ZV - 0.00240_JPRB * SIN(ZM-ZF-2.0_JPRB*ZD-ZN)
ZV = ZV + 0.00226_JPRB * SIN(ZM+ZF)
ZV = ZV - 0.00108_JPRB * SIN(ZM+ZF-2.0_JPRB*ZD)
ZV = ZV - 0.00079_JPRB * SIN(ZF-ZN)
ZV = ZV + 0.00078_JPRB * SIN(ZF+2.0_JPRB*ZD+ZN)
ZU = 1.0_JPRB - 0.10828_JPRB * COS(ZM)
ZU = ZU - 0.01880_JPRB * COS(ZM-2.0_JPRB*ZD)
ZU = ZU - 0.01479_JPRB * COS(2.0_JPRB*ZD)
ZU = ZU + 0.00181_JPRB * COS(2.0_JPRB*ZM-2.0_JPRB*ZD)
ZU = ZU - 0.00147_JPRB * COS(2.0_JPRB*ZM)
ZU = ZU - 0.00105_JPRB * COS(2.0_JPRB*ZD-ZG)
ZU = ZU - 0.00075_JPRB * COS(ZM-2.0_JPRB*ZD+ZG)
ZW = 0.10478_JPRB * SIN(ZM)
ZW = ZW - 0.04105_JPRB * SIN(2.0_JPRB*ZF+2.0_JPRB*ZN)
ZW = ZW - 0.02130_JPRB * SIN(ZM-2.0_JPRB*ZD)
ZW = ZW - 0.01779_JPRB * SIN(2.0_JPRB*ZF+ZN)
ZW = ZW + 0.01774_JPRB * SIN(ZN)
ZW = ZW + 0.00987_JPRB * SIN(2.0_JPRB*ZD)
ZW = ZW - 0.00338_JPRB * SIN(ZM-2.0_JPRB*ZF-2.0_JPRB*ZN)
ZW = ZW - 0.00309_JPRB * SIN(ZG)
ZW = ZW - 0.00190_JPRB * SIN(2.0_JPRB*ZF)
ZW = ZW - 0.00144_JPRB * SIN(ZM+ZN)
ZW = ZW - 0.00144_JPRB * SIN(ZM-2.0_JPRB*ZF-ZN)
ZW = ZW - 0.00113_JPRB * SIN(ZM+2.0_JPRB*ZF+2.0_JPRB*ZN)
ZW = ZW - 0.00094_JPRB * SIN(ZM-2.0_JPRB*ZD+ZG)
ZW = ZW - 0.00092_JPRB * SIN(2.0_JPRB*ZM-2.0_JPRB*ZD)

ZS = ZW / (ZU-ZV*ZV)**0.5_JPRB
ZRAM = ZL + ATAN(ZS/(1.0_JPRB-ZS*ZS)**0.5_JPRB)
ZS = ZV / ZU**0.5_JPRB
ZDECM = ATAN(ZS/(1.0_JPRB-ZS*ZS)**0.5_JPRB)
ZDISMOON = 60.40974_JPRB * ZU**0.5_JPRB
ZDISMOON = ZDISMOON * RA / 1000._JPRB

!     CALCUL DE LA PHASE DE LA LUNE
!     -----------------------------

ZPHI = ACOS(SIN(ZDECS)*SIN(ZDECM)+COS(ZDECS)*COS(ZDECM)*COS(ZRAS-ZRAM))
ZINC = ATAN2( ZDISS * SIN(ZPHI), ZDISMOON - ZDISS*COS(ZPHI) )
ZPHASE = (1.0_JPRB + COS(ZINC))/2.0_JPRB

!     PRISE EN COMPTE DES VARIATIONS DE LA DISTANCE TERRE-LUNE POUR
!     LE CALCUL DE L'ECLAIREMENT LUNAIRE HORS ATMOSPHERE
!     -------------------------------------------------------------

ZCORR=ZDISMOON_MOY/ZDISMOON
ZCORR=ZCORR*ZCORR
RIP0LU=ZECLMOON_MOY*ZPHASE*ZCORR

!     DETERMINATION DE LA POSITION DE LA LUNE EN LATITUDE ET LONGITUDE
!     ----------------------------------------------------------------

ZC(1) = 280.46061837_JPRB
ZC(2) = 360.98564736629_JPRB
ZC(3) = 0.000387933_JPRB
ZC(4) = 38710000._JPRB
ZT0 = ZTT - 1.0_JPRB

ZGST = ZC(1) + (ZC(2) * ZT) + ZT0**2*(ZC(3) - ZT0/ ZC(4) )

ZTMOLT=RPI-(ZRAM-ZGST*RPI/180._JPRB)
RTMOLT = MODULO(ZTMOLT,Z2PI)
RDECLU = ZDECM
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('UPDMOON',1,ZHOOK_HANDLE)

END SUBROUTINE UPDMOON
