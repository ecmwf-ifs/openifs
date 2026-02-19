! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUCLOPN(YDERAD,YDESWRT,KTSW,KSW,KLEV)

!**** *SUCLOP*  - INITIALIZE COMMON YOECLOP

!     PURPOSE.
!     --------
!           INITIALIZE YOMCLOP, WITH CLOUD OPTICAL PARAMETERS

!**   INTERFACE.
!     ----------
!        *CALL*  SUCLOPN
!        FROM *SUECRAD*

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOECLOP

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE
!     "INTEGRATED FORECASTING SYSTEM"

!     Fouquart, 1987: NATO ASI,  223-284
!     A. Slingo, 1989: J. Atmos. Sci., 46, 1419-1427
!     Ebert and Curry, 1992: J. Geophys. Res., 97D, 3831-3836
!     Sun and Shine, 1994: Quart. J. Roy. Meteor. Soc., 120, 111-138
!     Fu and Liou, 1993: J. Atmos. Sci., 50, 2008-2025
!     Fu, 1996: J. Climate, 9, 2058-2082
!     Fu et al., 1998: J. Climate, 11, 2223-2237
!     Lindner & Li, 2000: J. Climate, 13, 1797-1805

!     AUTHOR.
!     -------
!      JEAN-JACQUES MORCRETTE  *ECMWF*
!      ORIGINAL : 92-02-29

!     MODIFICATIONS.
!     --------------
!      01-01-16  JJ Morcrette  ice cloud from Fu 96, Fu et al 98
!      01-03-01  JJ Morcrette  water clouds from Lindner, Li 2000
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      2013-05-09 A Bozzo  bugfix in the parameter RFUETC from Fu et al 98 (bands 14,15,16 were left empty)
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOESW    , ONLY : TESWRT
USE YOERAD   , ONLY : TERAD
USE YOMLUN   , ONLY : NULOUT

IMPLICIT NONE

TYPE(TERAD)       ,INTENT(INOUT) :: YDERAD
TYPE(TESWRT)      ,INTENT(INOUT) :: YDESWRT
INTEGER(KIND=JPIM),INTENT(IN)    :: KTSW 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSW 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM) ::  JNU

!     -----------------------------------------------------------------
REAL(KIND=JPRB) :: ZEBCUA2(2)  ,ZEBCUB2(2)  ,ZEBCUC2(2)  ,ZEBCUD2(2)&
 & ,  ZEBCUE2(2)  ,ZEBCUF2(2)  ,ZYFWCA2(2)  ,ZYFWCB2(2)&
 & ,  ZYFWCC2(2)  ,ZYFWCD2(2)  ,ZYFWCE2(2)  ,ZYFWCF2(2)&
 & ,  ZASWCA2(2)  ,ZASWCB2(2)  ,ZASWCC2(2)  ,ZASWCD2(2)&
 & ,  ZASWCE2(2)  ,ZASWCF2(2) &
 & ,  ZKPNWCA2(2) ,ZKPNWCB2(2) ,ZKPNWCC2(2)&
 & ,  ZKPNWCD2(2) ,ZKPNWCE2(2) ,ZKPNWCF2(2) ,ZKPNWCG2(2) ,ZKPNWCH2(2)&
 & ,  ZSUSHE2(2)  ,ZSUSHF2(2)&
 & ,  ZSUSHH2(2)  ,ZSUSHK2(2)  ,ZSUSHA2(2)  ,ZSUSHG2(2)&
 & ,  ZFLAA02(2)  ,ZFLAA12(2)  ,ZFLBB02(2)  ,ZFLBB12(2)&
 & ,  ZFLBB22(2)  ,ZFLBB32(2)  ,ZFLCC02(2)  ,ZFLCC12(2)&
 & ,  ZFLCC22(2)  ,ZFLCC32(2)  ,ZFLDD02(2)  ,ZFLDD12(2)&
 & ,  ZFLDD22(2)  ,ZFLDD32(2)  

REAL(KIND=JPRB) :: ZEBCUA4(4)  ,ZEBCUB4(4)  ,ZEBCUC4(4)  ,ZEBCUD4(4)&
 & ,  ZEBCUE4(4)  ,ZEBCUF4(4)  ,ZYFWCA4(4)  ,ZYFWCB4(4)&
 & ,  ZYFWCC4(4)  ,ZYFWCD4(4)  ,ZYFWCE4(4)  ,ZYFWCF4(4)&
 & ,  ZASWCA4(4)  ,ZASWCB4(4)  ,ZASWCC4(4)  ,ZASWCD4(4)&
 & ,  ZASWCE4(4)  ,ZASWCF4(4)&
 & ,  ZKPNWCA4(4) ,ZKPNWCB4(4) ,ZKPNWCC4(4)&
 & ,  ZKPNWCD4(4) ,ZKPNWCE4(4) ,ZKPNWCF4(4) ,ZKPNWCG4(4) ,ZKPNWCH4(4)&
 & ,  ZSUSHE4(4)  ,ZSUSHF4(4)&
 & ,  ZSUSHH4(4)  ,ZSUSHK4(4)  ,ZSUSHA4(4)  ,ZSUSHG4(4)&
 & ,  ZFLAA04(4)  ,ZFLAA14(4)  ,ZFLBB04(4)  ,ZFLBB14(4)&
 & ,  ZFLBB24(4)  ,ZFLBB34(4)  ,ZFLCC04(4)  ,ZFLCC14(4)&
 & ,  ZFLCC24(4)  ,ZFLCC34(4)  ,ZFLDD04(4)  ,ZFLDD14(4)&
 & ,  ZFLDD24(4)  ,ZFLDD34(4)&
 & ,  ZFUAA04(4)  ,ZFUAA14(4)  ,ZFUBB04(4)  ,ZFUBB14(4)&
 & ,  ZFUBB24(4)  ,ZFUBB34(4)  ,ZFUCC04(4)  ,ZFUCC14(4)&
 & ,  ZFUCC24(4)  ,ZFUCC34(4)    
  
REAL(KIND=JPRB) :: ZAOVLP(3), ZBOVLP(3)  

REAL(KIND=JPRB) :: ZEBCUA6(6)  ,ZEBCUB6(6)  ,ZEBCUC6(6)  ,ZEBCUD6(6)&
 & ,  ZEBCUE6(6)  ,ZEBCUF6(6)  ,ZYFWCA6(6)  ,ZYFWCB6(6)&
 & ,  ZYFWCC6(6)  ,ZYFWCD6(6)  ,ZYFWCE6(6)  ,ZYFWCF6(6)&
 & ,  ZASWCA6(6)  ,ZASWCB6(6)  ,ZASWCC6(6)  ,ZASWCD6(6)&
 & ,  ZASWCE6(6)  ,ZASWCF6(6)&
 &  , ZKPNWCA6(6) ,ZKPNWCB6(6) ,ZKPNWCC6(6)&
 & ,  ZKPNWCD6(6) ,ZKPNWCE6(6) ,ZKPNWCF6(6) ,ZKPNWCG6(6) ,ZKPNWCH6(6)&
 & ,  ZSUSHE6(6)  ,ZSUSHF6(6)&
 & ,  ZSUSHH6(6)  ,ZSUSHK6(6)  ,ZSUSHA6(6)  ,ZSUSHG6(6)&
 & ,  ZFLAA06(6)  ,ZFLAA16(6)  ,ZFLBB06(6)  ,ZFLBB16(6)&
 & ,  ZFLBB26(6)  ,ZFLBB36(6)  ,ZFLCC06(6)  ,ZFLCC16(6)&
 & ,  ZFLCC26(6)  ,ZFLCC36(6)  ,ZFLDD06(6)  ,ZFLDD16(6)&
 & ,  ZFLDD26(6)  ,ZFLDD36(6)&
 & ,  ZFUAA06(6)  ,ZFUAA16(6)  ,ZFUBB06(6)  ,ZFUBB16(6)&
 & ,  ZFUBB26(6)  ,ZFUBB36(6)  ,ZFUCC06(6)  ,ZFUCC16(6)&
 & ,  ZFUCC26(6)  ,ZFUCC36(6)    

REAL(KIND=JPRB) :: ZEBCUA14(14)  ,ZEBCUB14(14)  ,ZEBCUC14(14)  ,ZEBCUD14(14)&
 & ,  ZEBCUE14(14)  ,ZEBCUF14(14)  ,ZYFWCA14(14)  ,ZYFWCB14(14)&
 & ,  ZYFWCC14(14)  ,ZYFWCD14(14)  ,ZYFWCE14(14)  ,ZYFWCF14(14)&
 & ,  ZASWCA14(14)  ,ZASWCB14(14)  ,ZASWCC14(14)  ,ZASWCD14(14)&
 & ,  ZASWCE14(14)  ,ZASWCF14(14)  ,ZSUSHE14(14)  ,ZSUSHF14(14)&
 & ,  ZSUSHH14(14)  ,ZSUSHK14(14)  ,ZSUSHA14(14)  ,ZSUSHG14(14)&
 & ,  ZFLAA014(14)  ,ZFLAA114(14)  ,ZFLBB014(14)  ,ZFLBB114(14)&
 & ,  ZFLBB214(14)  ,ZFLBB314(14)  ,ZFLCC014(14)  ,ZFLCC114(14)&
 & ,  ZFLCC214(14)  ,ZFLCC314(14)  ,ZFLDD014(14)  ,ZFLDD114(14)&
 & ,  ZFLDD214(14)  ,ZFLDD314(14)&
 & ,  ZFUAA014(14)  ,ZFUAA114(14)  ,ZFUBB014(14)  ,ZFUBB114(14)&
 & ,  ZFUBB214(14)  ,ZFUBB314(14)  ,ZFUCC014(14)  ,ZFUCC114(14)&
 & ,  ZFUCC214(14)  ,ZFUCC314(14)    
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     -----------------------------------------------------------------

!*          1.    SHORTWAVE CLOUD OPTICAL PROPERTIES
!                 ----------------------------------

!     ------------------------------------------------------------------

!*          1.1   TWO SPECTRAL INTERVALS
!                 ----------------------

! SW : 0.25 - 0.68 - 4.00 microns

!* Ice cloud properties - crystal: adapted from Ebert and Curry, 1992

!  optical properties
IF (LHOOK) CALL DR_HOOK('SUCLOPN',0,ZHOOK_HANDLE)
ASSOCIATE(LSRTM=>YDERAD%LSRTM, RAOVLP=>YDERAD%RAOVLP, RBOVLP=>YDERAD%RBOVLP, &
 & RYFWCA=>YDESWRT%RYFWCA, RYFWCB=>YDESWRT%RYFWCB, RYFWCC=>YDESWRT%RYFWCC, RYFWCD=>YDESWRT%RYFWCD, &
 & RYFWCE=>YDESWRT%RYFWCE, RYFWCF=>YDESWRT%RYFWCF, REBCUA=>YDESWRT%REBCUA, REBCUB=>YDESWRT%REBCUB, REBCUC=>YDESWRT%REBCUC, &
 & REBCUD=>YDESWRT%REBCUD, REBCUE=>YDESWRT%REBCUE, REBCUF=>YDESWRT%REBCUF, REBCUG=>YDESWRT%REBCUG, REBCUH=>YDESWRT%REBCUH, &
 & REBCUI =>YDESWRT%REBCUI, REBCUJ=>YDESWRT%REBCUJ, RASWCA=>YDESWRT%RASWCA, RASWCB=>YDESWRT%RASWCB, RASWCC=>YDESWRT%RASWCC, &
 & RASWCD=>YDESWRT%RASWCD, RASWCE=>YDESWRT%RASWCE, RASWCF=>YDESWRT%RASWCF ,&
 & RKPNWCA=>YDESWRT%RKPNWCA, RKPNWCB=>YDESWRT%RKPNWCB,RKPNWCC=>YDESWRT%RKPNWCC,&
 & RKPNWCD=>YDESWRT%RKPNWCD, RKPNWCE=>YDESWRT%RKPNWCE,RKPNWCF=>YDESWRT%RKPNWCF, RKPNWCG=>YDESWRT%RKPNWCG,&
 & RKPNWCH=>YDESWRT%RKPNWCH, REFFIA=>YDESWRT%REFFIA, RHSAVI=>YDESWRT%RHSAVI, &
 & RFUAA0=>YDESWRT%RFUAA0, RFUAA1=>YDESWRT%RFUAA1, RFUBB0=>YDESWRT%RFUBB0, RFUBB1=>YDESWRT%RFUBB1, RFUBB2=>YDESWRT%RFUBB2, &
 & RFUBB3=>YDESWRT%RFUBB3, RFUCC0=>YDESWRT%RFUCC0, RFUCC1=>YDESWRT%RFUCC1, RFUCC2=>YDESWRT%RFUCC2, RFUCC3=>YDESWRT%RFUCC3, &
 & RFUETA=>YDESWRT%RFUETA, RFUETB=>YDESWRT%RFUETB, RFUETC=>YDESWRT%RFUETC, RFULIO=>YDESWRT%RFULIO, RFLAA0=>YDESWRT%RFLAA0, &
 & RFLAA1=>YDESWRT%RFLAA1, RFLBB0=>YDESWRT%RFLBB0, &
 & RFLBB1=>YDESWRT%RFLBB1, RFLBB2=>YDESWRT%RFLBB2, RFLBB3=>YDESWRT%RFLBB3, RFLCC0=>YDESWRT%RFLCC0, RFLCC1=>YDESWRT%RFLCC1, &
 & RFLCC2=>YDESWRT%RFLCC2, RFLCC3=>YDESWRT%RFLCC3, RFLDD0=>YDESWRT%RFLDD0, RFLDD1=>YDESWRT%RFLDD1, RFLDD2=>YDESWRT%RFLDD2, &
 & RFLDD3=>YDESWRT%RFLDD3, RLILIA=>YDESWRT%RLILIA, RLILIB=>YDESWRT%RLILIB )
ZEBCUA2 = (/ 3.448E-04_JPRB , 3.448E-04_JPRB /)
ZEBCUB2 = (/ 2.431_JPRB     , 2.431_JPRB     /)
ZEBCUC2 = (/ 0.00001_JPRB   , 0.024366_JPRB  /)
ZEBCUD2 = (/ 0.0_JPRB         , 2.487E-04_JPRB /)
ZEBCUE2 = (/ 0.7661_JPRB    , 0.7866_JPRB    /)
ZEBCUF2 = (/ 5.851E-04_JPRB , 5.937E-04_JPRB /)

!  optical properties
!      ZEBCUA2 = (/ 3.448E-04 , 3.448E-04 /)
!      ZEBCUB2 = (/ 2.431     , 2.431     /)
!      ZEBCUC2 = (/ 0.00001   , 0.035589  /)
!      ZEBCUD2 = (/ 0.        , 2.757E-04 /)
!      ZEBCUE2 = (/ 0.7661    , 0.7921    /)
!      ZEBCUF2 = (/ 5.851E-04 , 5.893E-04 /)

!* Water cloud properties - from Fouquart (1987)

ZYFWCA2 = (/ 0.0_JPRB       , 0.0_JPRB       /)
ZYFWCB2 = (/ 1.5_JPRB       , 1.5_JPRB       /)
ZYFWCC2 = (/ 0.9999_JPRB    , 0.9988_JPRB    /)
ZYFWCD2 = (/ 5.000E-04_JPRB , 2.500E-03_JPRB /)
ZYFWCE2 = (/ 0.5_JPRB       , 0.05_JPRB      /)
ZYFWCF2 = (/ 0.865_JPRB     , 0.910_JPRB     /)

!* Water cloud properties - from Slingo (1989)

ZASWCA2 = (/ 2.817_JPRB     , 2.455_JPRB    /)
ZASWCB2 = (/ 1.305_JPRB     , 1.400_JPRB    /)
ZASWCC2 = (/-5.62E-08_JPRB  , 1.28E-02_JPRB /)
ZASWCD2 = (/ 1.63E-07_JPRB  , 8.97E-04_JPRB /)
ZASWCE2 = (/ 0.829_JPRB     , 0.783_JPRB    /)
ZASWCF2 = (/ 2.482_JPRB     , 4.988_JPRB    /)

!* Water cloud properties - from Nielsen (2013) - Approx. values.

ZKPNWCA2 = (/ 1.622_JPRB     , 1.73_JPRB     /)
ZKPNWCB2 = (/ 1.017_JPRB     , 1.03_JPRB     /)
ZKPNWCC2 = (/ 1.000_JPRB     , 0.999_JPRB    /)
ZKPNWCD2 = (/ 7.0E-08_JPRB   , 4.7E-04_JPRB  /)
ZKPNWCE2 = (/ 0.868_JPRB     , 0.866_JPRB    /)
ZKPNWCF2 = (/ 1.9E-04_JPRB   , 4.3E-04_JPRB  /)
ZKPNWCG2 = (/ 0.062_JPRB     , 0.106_JPRB    /)
ZKPNWCH2 = (/ 0.25_JPRB      , 0.194_JPRB    /)

!* Ice cloud properties - from Sun and Shine (1995)

ZSUSHE2 = (/ 0.0_JPRB       , 8.6822_JPRB    /)
ZSUSHF2 = (/ 0.0_JPRB       , 9.6277_JPRB    /)
ZSUSHH2 = (/ 0.8522_JPRB    , 0.8819_JPRB    /)
ZSUSHK2 = (/ 0.1620_JPRB    , 0.1630_JPRB    /)
ZSUSHA2 = (/ 0.0_JPRB       , 23.204_JPRB    /)
ZSUSHG2 = (/ 0.3270_JPRB    , 0.4180_JPRB    /)

!* Ice cloud properties - from Fu and Liou (1993)

ZFLAA02 = (/-6.656E-3_JPRB  ,-6.656E-3_JPRB  /)
ZFLAA12 = (/ 3.686_JPRB     , 3.686_JPRB     /)
ZFLBB02 = (/ .10998E-5_JPRB , .21136E-1_JPRB /)
ZFLBB12 = (/-.26101E-7_JPRB , .39150E-3_JPRB /)
ZFLBB22 = (/ .18096E-8_JPRB ,-.20740E-6_JPRB /)
ZFLBB32 = (/-.47387E-11_JPRB,-.28829E-8_JPRB /)
ZFLCC02 = (/ .22110E+1_JPRB , .22498E+1_JPRB /)
ZFLCC12 = (/-.10398E-2_JPRB , .23656E-3_JPRB /)
ZFLCC22 = (/ .65199E-4_JPRB , .51948E-4_JPRB /) 
ZFLCC32 = (/-.34498E-6_JPRB ,-.29768E-6_JPRB /)
ZFLDD02 = (/ .12495_JPRB    , .11716_JPRB    /)
ZFLDD12 = (/-.43582E-3_JPRB ,-.45208E-3_JPRB /)
ZFLDD22 = (/ .14092E-4_JPRB , .12772E-4_JPRB /)
ZFLDD32 = (/-.69565E-7_JPRB ,-.62779E-7_JPRB /) 

!* Ice cloud properties - from Fu (1996)

!ZFUAA02 = (/-.291721E-04_JPRB ,
!ZFUAA12 = (/ .251925E+01_JPRB ,
!ZFUBB02 = (/ .135403E-06_JPRB ,
!ZFUBB12 = (/ .992822E-07_JPRB ,
!ZFUBB22 = (/-.738432E-10_JPRB ,
!ZFUBB32 = (/ .331119E-12_JPRB ,
!ZFUCC02 = (/ .748127E+00_JPRB ,
!ZFUCC12 = (/ .956845E-03_JPRB ,
!ZFUCC22 = (/-.111517E-05_JPRB ,
!ZFUCC32 = (/-.815573E-08_JPRB ,
!ZFUDD02 = (/ .115730E+00_JPRB ,
!ZFUDD12 = (/ .256481E-03_JPRB ,
!ZFUDD22 = (/ .191313E-05_JPRB ,
!ZFUDD32 = (/-.124603E-07_JPRB ,
!     ------------------------------------------------------------------

!*          1.2    FOUR SPECTRAL INTERVALS
!                  -----------------------

! SW : 4 spectral intervals (0.25 - 0.69 - 1.19 - 2.38 - 4.00)

!* Ice cloud properties - crystal: adapted from Ebert and Curry, 1992

!      ZEBCUA4 = (/ 3.448E-03 , 3.448E-03 , 3.448E-03 , 3.448E-03 /)
!      ZEBCUB4 = (/ 2.431     , 2.431     , 2.431     , 2.431     /)
!      ZEBCUC4 = (/ 0.00001   , 0.00011   , 0.01861   , 0.46658   /)
!      ZEBCUD4 = (/ 0.        , 1.405E-05 , 8.328E-4  , 2.050E-05 /)
!      ZEBCUE4 = (/ 0.7661    , 0.7730    , 0.7940    , 0.9595    /)
!      ZEBCUF4 = (/ 5.851E-04 , 5.665E-04 , 7.267E-04 , 1.076E-04 /)

ZEBCUA4 = (/ 3.448E-03_JPRB , 3.448E-03_JPRB , 3.448E-03_JPRB , 3.448E-03_JPRB /)
ZEBCUB4 = (/ 2.431_JPRB     , 2.431_JPRB     , 2.431_JPRB     , 2.431_JPRB     /)
ZEBCUC4 = (/ 0.00001_JPRB   , 0.00011_JPRB   , 0.0197796_JPRB , 0.46658_JPRB   /)
ZEBCUD4 = (/ 0.0_JPRB       , 1.405E-05_JPRB , 7.95513E-4_JPRB, 2.050E-05_JPRB /)
ZEBCUE4 = (/ 0.7661_JPRB    , 0.7730_JPRB    , 0.795653_JPRB  , 0.9595_JPRB    /)
ZEBCUF4 = (/ 5.851E-04_JPRB , 5.665E-04_JPRB , 7.267E-04_JPRB , 1.076E-04_JPRB /)

!* Water cloud properties - from Fouquart (1987)

ZYFWCA4 = (/ 0.0_JPRB       , 0.0_JPRB       , 0.0_JPRB      , 0.0_JPRB        /)
ZYFWCB4 = (/ 1.5_JPRB       , 1.5_JPRB       , 1.5_JPRB      , 1.5_JPRB        /)
ZYFWCC4 = (/ 0.9999_JPRB    , 0.9988_JPRB    , 0.9988_JPRB   , 0.9988_JPRB     /)
ZYFWCD4 = (/ 5.000E-04_JPRB , 2.500E-03_JPRB , 2.500E-03_JPRB, 2.500E-03_JPRB  /)
ZYFWCE4 = (/ 0.5_JPRB       , 0.05_JPRB      , 0.05_JPRB     , 0.05_JPRB       /)
ZYFWCF4 = (/ 0.865_JPRB     , 0.910_JPRB     , 0.910_JPRB    , 0.910_JPRB      /)

!* Water cloud properties - from Slingo (1989)

ZASWCA4 = (/ 2.817_JPRB     , 2.682_JPRB     , 2.264_JPRB    , 1.281_JPRB      /)
ZASWCB4 = (/ 1.305_JPRB     , 1.346_JPRB     , 1.454_JPRB    , 1.641_JPRB      /)
ZASWCC4 = (/-5.62E-08_JPRB  ,-6.94E-06_JPRB  , 4.64E-04_JPRB , 2.01E-01_JPRB   /)
ZASWCD4 = (/ 1.63E-07_JPRB  , 2.35E-05_JPRB  , 1.24E-03_JPRB , 7.56E-03_JPRB   /)
ZASWCE4 = (/ 0.829_JPRB     , 0.794_JPRB     , 0.754_JPRB    , 0.826_JPRB      /)
ZASWCF4 = (/ 2.482_JPRB     , 4.226_JPRB     , 6.560_JPRB    , 4.353_JPRB      /)

!* Water cloud properties - from Nielsen (2013) - Approx. values.

ZKPNWCA4 = (/ 1.622_JPRB     , 1.685_JPRB     , 1.77_JPRB      , 1.87_JPRB     /)
ZKPNWCB4 = (/ 1.017_JPRB     , 1.024_JPRB     , 1.035_JPRB     , 1.046_JPRB    /)
ZKPNWCC4 = (/ 1.000_JPRB     , 0.99999_JPRB   , 0.9985_JPRB    , 0.823_JPRB    /)
ZKPNWCD4 = (/ 7.0E-08_JPRB     , 1.49E-05_JPRB  , 9.2E-04_JPRB   , 4.E-03_JPRB   /)
ZKPNWCE4 = (/ 0.868_JPRB     , 0.867_JPRB     , 0.864_JPRB     , 0.886_JPRB    /)
ZKPNWCF4 = (/ 1.9E-04_JPRB   , 3.1E-04_JPRB   , 5.4E-04_JPRB   , 1.1E-03_JPRB  /)
ZKPNWCG4 = (/ 0.062_JPRB     , 0.078_JPRB     , 0.133_JPRB     , 0.20_JPRB     /)
ZKPNWCH4 = (/ 0.25_JPRB      , 0.195_JPRB     , 0.194_JPRB     , 0.18_JPRB     /)

!* Ice cloud properties - from Sun and Shine (1995)

ZSUSHE4 = (/ 0.0_JPRB     , 7.2471E-02_JPRB , 17.5933_JPRB  , 48.7166_JPRB  /)
ZSUSHF4 = (/ 0.0_JPRB     , 4.01511_JPRB    , 21.1249_JPRB  , 1.20890_JPRB  /)
ZSUSHH4 = (/ 0.8522_JPRB  , 0.85841_JPRB    , 0.90778_JPRB  ,0.982046_JPRB  /)
ZSUSHK4 = (/ 0.1620_JPRB  , 0.160048_JPRB   , 0.188521_JPRB ,0.0411446_JPRB /)
ZSUSHA4 = (/ 0.0_JPRB     , 0.273455_JPRB   , 41.7675_JPRB  , 161.104_JPRB  /)
ZSUSHG4 = (/ 0.3270_JPRB  , 0.343668_JPRB   , 0.526192_JPRB , 0.574040_JPRB /)

!* Ice cloud properties - from Fu and Liou (1993)

ZFLAA04 = (/-6.656E-03_JPRB ,-6.656E-03_JPRB,-6.656E-03_JPRB,-6.656E-03_JPRB /)
ZFLAA14 = (/ 3.686_JPRB     , 3.686_JPRB    , 3.686_JPRB    , 3.686_JPRB     /)
ZFLBB04 = (/ .10998E-5_JPRB , .20208E-4_JPRB, .51557E-3_JPRB, .39517E+0_JPRB /)
ZFLBB14 = (/-.26101E-7_JPRB , .96483E-5_JPRB, .10731E-2_JPRB, .15787E-2_JPRB /)
ZFLBB24 = (/ .18096E-8_JPRB , .83009E-7_JPRB, .17753E-5_JPRB,-.14337E-4_JPRB /)
ZFLBB34 = (/-.47387E-11_JPRB,-.32217E-9_JPRB,-.18379E-7_JPRB, .46942E-7_JPRB /)
ZFLCC04 = (/ .22110E+1_JPRB , .22151E+1_JPRB, .22534E+1_JPRB, .26653E+1_JPRB /)
ZFLCC14 = (/-.10398E-2_JPRB ,-.77982E-3_JPRB, .16163E-2_JPRB, .56935E-2_JPRB /)
ZFLCC24 = (/ .65199E-4_JPRB , .63750E-4_JPRB, .44037E-4_JPRB,-.54316E-4_JPRB /)
ZFLCC34 = (/-.34498E-6_JPRB ,-.34466E-6_JPRB,-.27627E-6_JPRB, .17858E-6_JPRB /)
ZFLDD04 = (/ .12495_JPRB    , .12363_JPRB   , .11983_JPRB   , .21834E-1_JPRB /)
ZFLDD14 = (/-.43582E-3_JPRB ,-.44419E-3_JPRB,-.50108E-3_JPRB,-.29204E-3_JPRB /)
ZFLDD24 = (/ .14092E-4_JPRB , .14038E-4_JPRB, .11843E-4_JPRB, .18060E-5_JPRB /)
ZFLDD34 = (/-.69565E-7_JPRB ,-.68851E-7_JPRB,-.59367E-7_JPRB,-.46257E-8_JPRB /)

!* Ice cloud properties - from Fu (1996) as tabulated in Sun & Rikus (1999)

ZFUAA04 = (/-1.30817E-04_JPRB,-6.39479E-05_JPRB,-6.74730E-06_JPRB, 1.62674E-04_JPRB /)
ZFUAA14 = (/ 2.52883E+00_JPRB, 2.52393E+00_JPRB, 2.52056E+00_JPRB, 2.49823E+00_JPRB /)
ZFUBB04 = (/-1.55357E-08_JPRB,-7.90657E-07_JPRB, 9.64842E-04_JPRB, 2.25112E-01_JPRB /)
ZFUBB14 = (/ 1.95793E-07_JPRB, 7.79991E-06_JPRB, 9.09809E-04_JPRB, 3.05017E-03_JPRB /)
ZFUBB24 = (/-2.31234E-10_JPRB, 2.90894E-10_JPRB,-3.57557E-06_JPRB,-2.54236E-05_JPRB /)
ZFUBB34 = (/ 1.12247E-12_JPRB,-2.02818E-12_JPRB, 1.00197E-08_JPRB, 8.49116E-08_JPRB /)
ZFUCC04 = (/ 7.39781E-01_JPRB, 7.52335E-01_JPRB, 7.56307E-01_JPRB, 8.30812E-01_JPRB /)
ZFUCC14 = (/ 9.10564E-04_JPRB, 1.06211E-03_JPRB, 1.73364E-03_JPRB, 2.62788E-03_JPRB /)
ZFUCC24 = (/-4.62479E-07_JPRB,-2.45770E-06_JPRB,-8.92191E-06_JPRB,-2.43196E-05_JPRB /)
ZFUCC34 = (/-1.05910E-08_JPRB,-3.03712E-09_JPRB, 1.97757E-08_JPRB, 8.23543E-08_JPRB /)

!     ------------------------------------------------------------------

!*          1.3    SIX SPECTRAL INTERVALS
!                  ----------------------

! SW : 6 spectral intervals (0.185-0.25-0.44-0.69-1.19-2.38-4.00)

!* Ice cloud properties - crystal: adapted from Ebert and Curry, 1992

ZEBCUA6 = (/ 3.448E-03_JPRB , 3.448E-03_JPRB , 3.448E-03_JPRB , 3.448E-03_JPRB , 3.448E-03_JPRB , 3.448E-03_JPRB /)
ZEBCUB6 = (/ 2.431_JPRB     , 2.431_JPRB     , 2.431_JPRB     , 2.431_JPRB     , 2.431_JPRB     , 2.431_JPRB     /)
ZEBCUC6 = (/ 0.00001_JPRB   , 0.00001_JPRB   , 0.00001_JPRB   , 0.00011_JPRB   , 0.0197796_JPRB , 0.46658_JPRB   /)
ZEBCUD6 = (/ 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 1.405E-05_JPRB , 7.95513E-4_JPRB, 2.050E-05_JPRB /)
ZEBCUE6 = (/ 0.7661_JPRB    , 0.7661_JPRB    , 0.7661_JPRB    , 0.7730_JPRB    , 0.795653_JPRB  , 0.9595_JPRB    /)
ZEBCUF6 = (/ 5.851E-04_JPRB , 5.851E-04_JPRB , 5.851E-04_JPRB , 5.665E-04_JPRB , 7.267E-04_JPRB , 1.076E-04_JPRB /)

!* Water cloud properties - from Fouquart (1987)

ZYFWCA6 = (/ 0.0_JPRB       , 0.0_JPRB       , 0.0_JPRB       , 0.0_JPRB       , 0.0_JPRB      , 0.0_JPRB        /)
ZYFWCB6 = (/ 1.5_JPRB       , 1.5_JPRB       , 1.5_JPRB       , 1.5_JPRB       , 1.5_JPRB      , 1.5_JPRB        /)
ZYFWCC6 = (/ 0.9999_JPRB    , 0.9999_JPRB    , 0.9999_JPRB    , 0.9988_JPRB    , 0.9988_JPRB   , 0.9988_JPRB     /)
ZYFWCD6 = (/ 5.000E-04_JPRB , 5.000E-04_JPRB , 5.000E-04_JPRB , 2.500E-03_JPRB , 2.500E-03_JPRB, 2.500E-03_JPRB  /)
ZYFWCE6 = (/ 0.5_JPRB       , 0.5_JPRB       , 0.5_JPRB       , 0.05_JPRB      , 0.05_JPRB     , 0.05_JPRB       /)
ZYFWCF6 = (/ 0.865_JPRB     , 0.865_JPRB     , 0.865_JPRB     , 0.910_JPRB     , 0.910_JPRB    , 0.910_JPRB      /)

!* Water cloud properties - from Slingo (1989)

ZASWCA6 = (/ 3.094_JPRB     , 2.881_JPRB     , 2.790_JPRB     , 2.682_JPRB     , 2.254_JPRB    , 1.281_JPRB      /)
ZASWCB6 = (/ 1.252_JPRB     , 1.284_JPRB     , 1.313_JPRB     , 1.346_JPRB     , 1.456_JPRB    , 1.641_JPRB      /)
ZASWCC6 = (/ 7.90E-07_JPRB  , 1.77E-07_JPRB  ,-1.55E-07_JPRB  ,-6.94E-06_JPRB  , 4.75E-04_JPRB , 2.01E-01_JPRB   /)
ZASWCD6 = (/ 3.69E-07_JPRB  , 1.23E-07_JPRB  , 1.80E-07_JPRB  , 2.35E-05_JPRB  , 1.26E-03_JPRB , 7.56E-03_JPRB   /)
ZASWCE6 = (/ 0.844_JPRB     , 0.839_JPRB     , 0.824_JPRB     , 0.794_JPRB     , 0.754_JPRB    , 0.826_JPRB      /)
ZASWCF6 = (/ 1.558_JPRB     , 1.914_JPRB     , 2.723_JPRB     , 4.226_JPRB     , 6.561_JPRB    , 4.353_JPRB      /)

!* Water cloud properties - from Nielsen (2013)

ZKPNWCA6 = (/ 1.500_JPRB     , 1.606_JPRB     , 1.638_JPRB     , 1.685_JPRB     , 1.77_JPRB     , 1.87_JPRB     /)
ZKPNWCB6 = (/ 1.000_JPRB     , 1.015_JPRB     , 1.019_JPRB     , 1.024_JPRB     , 1.035_JPRB    , 1.046_JPRB    /)
ZKPNWCC6 = (/ 1.000_JPRB     , 1.000_JPRB     , 1.000_JPRB     , 0.99999_JPRB   , 0.9985_JPRB   , 0.823_JPRB    /)
ZKPNWCD6 = (/ 1.E-07_JPRB    , 3.3E-08_JPRB   , 1.E-07_JPRB    , 1.49E-05_JPRB  , 9.2E-04_JPRB  , 4.E-03_JPRB   /)
ZKPNWCE6 = (/ 0.868_JPRB     , 0.868_JPRB     , 0.868_JPRB     , 0.867_JPRB     , 0.864_JPRB    , 0.886_JPRB    /)
ZKPNWCF6 = (/ 1.4E-04_JPRB   , 1.4E-04_JPRB   , 2.5E-04_JPRB   , 3.1E-04_JPRB   , 5.4E-04_JPRB  , 1.1E-03_JPRB  /)
ZKPNWCG6 = (/ 0.061_JPRB     , 0.061_JPRB     , 0.063_JPRB     , 0.078_JPRB     , 0.133_JPRB    , 0.20_JPRB     /)
ZKPNWCH6 = (/ 0.25_JPRB      , 0.25_JPRB      , 0.25_JPRB      , 0.195_JPRB     , 0.194_JPRB    , 0.18_JPRB     /)


!* Ice cloud properties - from Sun and Shine (1995)

ZSUSHE6 = (/ 0.0_JPRB     , 0.0_JPRB     , 0.0_JPRB     , 7.2471E-02_JPRB , 17.5933_JPRB  , 48.7166_JPRB  /)
ZSUSHF6 = (/ 0.0_JPRB     , 0.0_JPRB     , 0.0_JPRB     , 4.01511_JPRB    , 21.1249_JPRB  , 1.20890_JPRB  /)
ZSUSHH6 = (/ 0.8522_JPRB  , 0.8522_JPRB  , 0.8522_JPRB  , 0.85841_JPRB    , 0.90778_JPRB  ,0.982046_JPRB  /)
ZSUSHK6 = (/ 0.1620_JPRB  , 0.1620_JPRB  , 0.1620_JPRB  , 0.160048_JPRB   , 0.188521_JPRB ,0.0411446_JPRB /)
ZSUSHA6 = (/ 0.0_JPRB     , 0.0_JPRB     , 0.0_JPRB     , 0.273455_JPRB   , 41.7675_JPRB  , 161.104_JPRB  /)
ZSUSHG6 = (/ 0.3270_JPRB  , 0.3270_JPRB  , 0.3270_JPRB  , 0.343668_JPRB   , 0.526192_JPRB , 0.574040_JPRB /)

!* Ice cloud properties - from Fu and Liou (1993)

ZFLAA06 = (/-6.656E-03_JPRB ,-6.656E-03_JPRB ,-6.656E-03_JPRB ,-6.656E-03_JPRB,-6.656E-03_JPRB,-6.656E-03_JPRB /)
ZFLAA16 = (/ 3.686_JPRB     , 3.686_JPRB     , 3.686_JPRB     , 3.686_JPRB    , 3.686_JPRB    , 3.686_JPRB     /)
ZFLBB06 = (/ .10998E-5_JPRB , .10998E-5_JPRB , .10998E-5_JPRB , .20208E-4_JPRB, .51557E-3_JPRB, .39517E+0_JPRB /)
ZFLBB16 = (/-.26101E-7_JPRB ,-.26101E-7_JPRB ,-.26101E-7_JPRB , .96483E-5_JPRB, .10731E-2_JPRB, .15787E-2_JPRB /)
ZFLBB26 = (/ .18096E-8_JPRB , .18096E-8_JPRB , .18096E-8_JPRB , .83009E-7_JPRB, .17753E-5_JPRB,-.14337E-4_JPRB /)
ZFLBB36 = (/-.47387E-11_JPRB,-.47387E-11_JPRB,-.47387E-11_JPRB,-.32217E-9_JPRB,-.18379E-7_JPRB, .46942E-7_JPRB /)
ZFLCC06 = (/ .22110E+1_JPRB , .22110E+1_JPRB , .22110E+1_JPRB , .22151E+1_JPRB, .22534E+1_JPRB, .26653E+1_JPRB /)
ZFLCC16 = (/-.10398E-2_JPRB ,-.10398E-2_JPRB ,-.10398E-2_JPRB ,-.77982E-3_JPRB, .16163E-2_JPRB, .56935E-2_JPRB /)
ZFLCC26 = (/ .65199E-4_JPRB , .65199E-4_JPRB , .65199E-4_JPRB , .63750E-4_JPRB, .44037E-4_JPRB,-.54316E-4_JPRB /)
ZFLCC36 = (/-.34498E-6_JPRB ,-.34498E-6_JPRB ,-.34498E-6_JPRB ,-.34466E-6_JPRB,-.27627E-6_JPRB, .17858E-6_JPRB /)
ZFLDD06 = (/ .12495_JPRB    , .12495_JPRB    , .12495_JPRB    , .12363_JPRB   , .11983_JPRB   , .21834E-1_JPRB /)
ZFLDD16 = (/-.43582E-3_JPRB ,-.43582E-3_JPRB ,-.43582E-3_JPRB ,-.44419E-3_JPRB,-.50108E-3_JPRB,-.29204E-3_JPRB /)
ZFLDD26 = (/ .14092E-4_JPRB , .14092E-4_JPRB , .14092E-4_JPRB , .14038E-4_JPRB, .11843E-4_JPRB, .18060E-5_JPRB /)
ZFLDD36 = (/-.69565E-7_JPRB ,-.69565E-7_JPRB ,-.69565E-7_JPRB ,-.68851E-7_JPRB,-.59367E-7_JPRB,-.46257E-8_JPRB /)

!* Ice cloud properties - from Fu (1996) as tabulated in Sun & Rikus (1999)

ZFUAA06 = (/-1.30817E-04_JPRB,-1.30817E-04_JPRB,-1.30817E-04_JPRB,-6.39479E-05_JPRB,-6.74730E-06_JPRB, 1.62674E-04_JPRB /)
ZFUAA16 = (/ 2.52883E+00_JPRB, 2.52883E+00_JPRB, 2.52883E+00_JPRB, 2.52393E+00_JPRB, 2.52056E+00_JPRB, 2.49823E+00_JPRB /)
ZFUBB06 = (/-1.55357E-08_JPRB,-1.55357E-08_JPRB,-1.55357E-08_JPRB,-7.90657E-07_JPRB, 9.64842E-04_JPRB, 2.25112E-01_JPRB /)
ZFUBB16 = (/ 1.95793E-07_JPRB, 1.95793E-07_JPRB, 1.95793E-07_JPRB, 7.79991E-06_JPRB, 9.09809E-04_JPRB, 3.05017E-03_JPRB /)
ZFUBB26 = (/-2.31234E-10_JPRB,-2.31234E-10_JPRB,-2.31234E-10_JPRB, 2.90894E-10_JPRB,-3.57557E-06_JPRB,-2.54236E-05_JPRB /)
ZFUBB36 = (/ 1.12247E-12_JPRB, 1.12247E-12_JPRB, 1.12247E-12_JPRB,-2.02818E-12_JPRB, 1.00197E-08_JPRB, 8.49116E-08_JPRB /)
ZFUCC06 = (/ 7.39781E-01_JPRB, 7.39781E-01_JPRB, 7.39781E-01_JPRB, 7.52335E-01_JPRB, 7.56307E-01_JPRB, 8.30812E-01_JPRB /)
ZFUCC16 = (/ 9.10564E-04_JPRB, 9.10564E-04_JPRB, 9.10564E-04_JPRB, 1.06211E-03_JPRB, 1.73364E-03_JPRB, 2.62788E-03_JPRB /)
ZFUCC26 = (/-4.62479E-07_JPRB,-4.62479E-07_JPRB,-4.62479E-07_JPRB,-2.45770E-06_JPRB,-8.92191E-06_JPRB,-2.43196E-05_JPRB /)
ZFUCC36 = (/-1.05910E-08_JPRB,-1.05910E-08_JPRB,-1.05910E-08_JPRB,-3.03712E-09_JPRB, 1.97757E-08_JPRB, 8.23543E-08_JPRB /)

!     ------------------------------------------------------------------

!*          1.4    FOURTEEN SPECTRAL INTERVALS FOR RRTM_SW
!                  ---------------------------------------

! SW : 14 spectral intervals 
!  3.846 -  3.077
!  3.077 -  2.500
!  2.500 -  2.150
!  2.150 -  1.942
!  1.942 -  1.626
!  1.626 -  1.299
!  1.299 -  1.242
!  1.242 -  0.7782
!  0.7782-  0.6250
!  0.6250-  0.4415
!  0.4415-  0.3448
!  0.3448-  0.2632
!  0.2632-  0.2000
! 12.195 -  3.846

!* Ice cloud properties - crystal: remapped from Ebert and Curry, 1992

ZEBCUA14 = (/ &
 & 3.448E-03_JPRB , 3.448E-03_JPRB , 3.448E-03_JPRB , 3.448E-03_JPRB , 3.448E-03_JPRB , 3.448E-03_JPRB &
 & , 3.448E-03_JPRB , 3.448E-03_JPRB , 3.448E-03_JPRB , 3.448E-03_JPRB , 3.448E-03_JPRB , 3.448E-03_JPRB &
 & , 3.448E-03_JPRB , 3.448E-03_JPRB /)  

ZEBCUB14 = (/ &
 & 2.431_JPRB     , 2.431_JPRB     , 2.431_JPRB     , 2.431_JPRB     , 2.431_JPRB     , 2.431_JPRB     &
 & , 2.431_JPRB     , 2.431_JPRB     , 2.431_JPRB     , 2.431_JPRB     , 2.431_JPRB     , 2.431_JPRB     &
 & , 2.431_JPRB     , 2.431_JPRB     /)  

ZEBCUC14 = (/ &
 & 0.46658_JPRB   , 0.46658_JPRB   , 0.03779_JPRB   , 0.03779_JPRB   , 0.01240_JPRB   , 0.01240_JPRB   &
 & , 0.00011_JPRB   , 0.00011_JPRB   , 0.00001_JPRB   , 0.00001_JPRB   , 0.00001_JPRB   , 0.00001_JPRB   &
 & , 0.00001_JPRB   , 0.46658_JPRB   /)  

ZEBCUD14 = (/ &
 & 2.050E-05_JPRB , 2.050E-05_JPRB , 1.284E-03_JPRB , 1.284E-03_JPRB , 6.867E-04_JPRB , 6.867E-04_JPRB &
 & , 1.405E-05_JPRB , 1.405E-05_JPRB , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB          &
 & , 0.0_JPRB         , 2.050E-05_JPRB /)  

ZEBCUE14 = (/ &
 & 0.9595_JPRB    , 0.9595_JPRB    , 0.8172_JPRB    , 0.8172_JPRB    , 0.7865_JPRB    , 0.7865_JPRB    &
 & , 0.7730_JPRB    , 0.7730_JPRB    , 0.7661_JPRB    , 0.7661_JPRB    , 0.7661_JPRB    , 0.7661_JPRB    &
 & , 0.7661_JPRB    , 0.9595_JPRB    /)  

ZEBCUF14 = (/ &
 & 1.076E-04_JPRB , 1.076E-04_JPRB , 7.463E-04_JPRB , 7.463E-04_JPRB , 7.204E-04_JPRB , 7.204E-04_JPRB &
 & , 5.665E-04_JPRB , 5.665E-04_JPRB , 5.851E-04_JPRB , 5.851E-04_JPRB , 5.851E-04_JPRB , 5.851E-04_JPRB &
 & , 5.851E-04_JPRB , 1.076E-04_JPRB /)  

!* Water cloud properties - from Fouquart (1987)

ZYFWCA14 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
 
ZYFWCB14 = (/ &
 & 1.5_JPRB       , 1.5_JPRB       , 1.5_JPRB       , 1.5_JPRB       , 1.5_JPRB      , 1.5_JPRB        &
 & , 1.5_JPRB       , 1.5_JPRB       , 1.5_JPRB       , 1.5_JPRB       , 1.5_JPRB      , 1.5_JPRB        &
 & , 1.5_JPRB       , 1.5_JPRB       /)  
 
ZYFWCC14 = (/ &
 & 0.9988_JPRB    , 0.9988_JPRB    , 0.9988_JPRB    , 0.9988_JPRB    , 0.9988_JPRB    , 0.9988_JPRB    &
 & , 0.9988_JPRB    , 0.9988_JPRB    , 0.9999_JPRB    , 0.9999_JPRB    , 0.9999_JPRB    , 0.9999_JPRB    &
 & , 0.9999_JPRB    , 0.9988_JPRB    /)  

ZYFWCD14 = (/ &
 & 2.500E-03_JPRB , 2.500E-03_JPRB , 2.500E-03_JPRB , 2.500E-03_JPRB , 2.500E-03_JPRB , 2.500E-03_JPRB &
 & , 2.500E-03_JPRB , 2.500E-03_JPRB , 5.000E-04_JPRB , 5.000E-04_JPRB , 5.000E-04_JPRB , 5.000E-04_JPRB &
 & , 5.000E-04_JPRB , 2.500E-03_JPRB /)  

ZYFWCE14 = (/ &
 & 0.05_JPRB      , 0.05_JPRB      , 0.05_JPRB      , 0.05_JPRB      , 0.05_JPRB      , 0.05_JPRB      &
 & , 0.05_JPRB      , 0.05_JPRB      , 0.5_JPRB         , 0.5_JPRB         , 0.5_JPRB         , 0.5_JPRB               &
 & , 0.5_JPRB         , 0.05_JPRB      /)  

ZYFWCF14 = (/ &
 & 0.910_JPRB     , 0.910_JPRB     , 0.910_JPRB     , 0.910_JPRB     , 0.910_JPRB     , 0.910_JPRB     &
 & , 0.910_JPRB     , 0.910_JPRB     , 0.865_JPRB     , 0.865_JPRB     , 0.865_JPRB     , 0.865_JPRB      &
 & , 0.865_JPRB     , 0.910_JPRB     /)  

!* Water cloud properties - from Slingo (1989)

ZASWCA14 = (/ &
 & -1.023_JPRB     , 1.950_JPRB     , 1.850_JPRB     , 1.970_JPRB     , 1.970_JPRB     , 2.463_JPRB     &
 & , 2.551_JPRB     , 2.622_JPRB     , 2.895_JPRB     , 2.672_JPRB     , 2.801_JPRB    , 2.944_JPRB     &
 & , 3.094_JPRB     ,-1.023_JPRB     /)  

ZASWCB14 = (/ &
 & 1.933_JPRB     , 1.540_JPRB     , 1.556_JPRB     , 1.501_JPRB     , 1.501_JPRB     , 1.420_JPRB     &
 & , 1.401_JPRB     , 1.362_JPRB     , 1.315_JPRB     , 1.320_JPRB     , 1.293_JPRB     , 1.270_JPRB     &
 & , 1.252_JPRB     , 1.933_JPRB     /)  
   
ZASWCC14 = (/ &
 & 2.50E-02_JPRB  , 4.49E-01_JPRB  , 1.90E-04_JPRB  , 1.20E-03_JPRB  , 1.20E-03_JPRB  , 2.40E-04_JPRB  &
 & , 6.20E-05_JPRB  , 3.30E-06_JPRB  ,-1.20E-07_JPRB  , 0.0_JPRB         , 1.00E-06_JPRB  ,-6.50E-07_JPRB  &
 & , 7.90E-07_JPRB  , 2.50E-02_JPRB  /)  
    
ZASWCD14 = (/ &
 & 1.22E-02_JPRB  , 1.54E-03_JPRB  , 2.54E-03_JPRB  , 2.16E-03_JPRB  , 2.16E-03_JPRB  , 8.56E-04_JPRB  &
 & , 2.60E-04_JPRB  , 2.80E-06_JPRB  , 4.40E-07_JPRB  , 0.0_JPRB         , 0.0_JPRB         , 4.33E-07_JPRB  &
 & , 3.69E-07_JPRB  , 1.22E-02_JPRB  /)  

ZASWCE14 = (/ &
 & 0.726_JPRB     , 0.831_JPRB     , 0.769_JPRB     , 0.740_JPRB     , 0.740_JPRB     , 0.754_JPRB     &
 & , 0.773_JPRB     , 0.806_JPRB     , 0.818_JPRB     , 0.828_JPRB     , 0.836_JPRB     , 0.841_JPRB     &
 & , 0.844_JPRB     , 0.726_JPRB     /)  

ZASWCF14 = (/ &
 & 6.652_JPRB     , 6.102_JPRB     , 5.171_JPRB     , 7.469_JPRB     , 7.469_JPRB     , 6.555_JPRB     &
 & , 5.405_JPRB     , 3.355_JPRB     , 2.989_JPRB     , 2.467_JPRB     , 2.153_JPRB     , 1.680_JPRB     &
 & , 1.558_JPRB     , 6.652_JPRB     /)  

!* Ice cloud properties - from Sun and Shine (1995)

ZSUSHE14 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZSUSHF14 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZSUSHH14 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZSUSHK14 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZSUSHA14 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZSUSHG14 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  

!* Ice cloud properties - from Fu and Liou (1993)

ZFLAA014 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZFLAA114 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZFLBB014 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZFLBB114 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZFLBB214 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZFLBB314 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZFLCC014 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZFLCC114 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZFLCC214 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZFLCC314 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZFLDD014 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZFLDD114 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZFLDD214 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  
ZFLDD314 = (/ &
 & 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB         , 0.0_JPRB        , 0.0_JPRB          &
 & , 0.0_JPRB         , 0.0_JPRB         /)  

!* Ice cloud properties - from Fu (1996) 

ZFUAA014 = (/ &
 & 1.87598E-04_JPRB, 2.97295E-04_JPRB, 4.89477E-04_JPRB,-8.37325E-06_JPRB,-8.37325E-06_JPRB,-8.05155E-04_JPRB  &
 & , 6.51659E-05_JPRB, 8.10443E-05_JPRB, 1.61983E-04_JPRB,-9.45458E-05_JPRB,-2.58858E-04_JPRB,-2.66955E-04_JPRB &
 & ,-2.36447E-04_JPRB,-2.54823E-04_JPRB /)  

ZFUAA114 = (/ &
 & 2.51396E+00_JPRB, 2.48895E+00_JPRB, 2.48776E+00_JPRB, 2.52504E+00_JPRB, 2.52504E+00_JPRB, 2.57600E+00_JPRB &
 & , 2.51660E+00_JPRB, 2.51619E+00_JPRB, 2.50746E+00_JPRB, 2.52061E+00_JPRB, 2.53815E+00_JPRB, 2.54179E+00_JPRB &
 & , 2.53817E+00_JPRB, 2.52909E+00_JPRB /)  

ZFUBB014 = (/ &
 & 1.96793E-01_JPRB, 4.64416E-01_JPRB, 5.83469E-04_JPRB, 2.53234E-03_JPRB, 2.53234E-03_JPRB,-2.85518E-05_JPRB &
 & ,-1.48012E-07_JPRB,-1.57963E-07_JPRB,-7.78001E-08_JPRB, 5.08447E-07_JPRB,-1.98529E-07_JPRB,-1.00570E-07_JPRB &
 & ,-2.69916E-07_JPRB, 2.60155E-01_JPRB /)   

ZFUBB114 = (/ &
 & 5.75235E-03_JPRB, 2.04716E-05_JPRB, 1.18127E-03_JPRB, 1.75078E-03_JPRB, 1.75078E-03_JPRB, 1.71993E-03_JPRB &
 & , 9.02355E-05_JPRB, 1.72475E-06_JPRB, 2.53360E-07_JPRB, 2.73206E-08_JPRB, 9.39480E-08_JPRB, 1.60441E-07_JPRB &
 & , 2.12909E-07_JPRB, 5.45547E-03_JPRB /)  

ZFUBB214 = (/ &
 & -5.29220E-05_JPRB,-4.60375E-07_JPRB,-3.40011E-06_JPRB,-8.00994E-06_JPRB,-8.00994E-06_JPRB,-7.43697E-06_JPRB &
 & ,-1.98190E-08_JPRB, 9.02156E-11_JPRB,-1.15489E-10_JPRB, 4.96553E-11_JPRB,-2.54540E-10_JPRB,-2.05663E-10_JPRB &
 & ,-2.65397E-10_JPRB,-5.58760E-05_JPRB /)  

ZFUBB314 = (/ &
 & 1.76618E-07_JPRB, 2.03198E-09_JPRB, 8.78549E-09_JPRB, 2.31309E-08_JPRB, 2.31309E-08_JPRB, 2.09647E-08_JPRB &
 & , 4.01914E-11_JPRB,-3.79423E-13_JPRB, 4.65084E-13_JPRB,-1.86001E-13_JPRB, 1.10876E-12_JPRB, 8.88595E-13_JPRB &
 & , 1.12983E-12_JPRB, 1.97086E-07_JPRB /)  

ZFUCC014 = (/ &
 & 0.759183_JPRB   , 0.919599_JPRB   , 0.775916_JPRB   , 0.758748_JPRB   , 0.758748_JPRB   , 0.752528_JPRB    &
 & , 0.751277_JPRB   , 0.752318_JPRB   , 0.749693_JPRB   , 0.749856_JPRB   , 0.743546_JPRB   , 0.737809_JPRB    &
 & , 0.733260_JPRB   , 0.799084_JPRB   /)  
   
ZFUCC114 = (/ &
 & 4.93765E-03_JPRB, 5.03025E-04_JPRB, 1.74517E-03_JPRB, 2.02709E-03_JPRB, 2.02709E-03_JPRB, 1.95748E-03_JPRB &
 & , 1.29824E-03_JPRB, 1.04224E-03_JPRB, 1.05446E-03_JPRB, 8.89161E-04_JPRB, 9.08674E-04_JPRB, 8.97515E-04_JPRB &
 & , 9.18317E-04_JPRB, 4.81706E-03_JPRB /)  

ZFUCC214 = (/ &
 & -4.84059E-05_JPRB,-5.74771E-06_JPRB,-9.21314E-06_JPRB,-1.17029E-05_JPRB,-1.17029E-05_JPRB,-1.02495E-05_JPRB &
 & ,-4.99075E-06_JPRB,-2.26618E-06_JPRB,-2.32576E-06_JPRB,-3.49578E-07_JPRB,-4.65326E-07_JPRB,-2.17099E-07_JPRB &
 & ,-4.22974E-07_JPRB,-5.13220E-05_JPRB /)  

ZFUCC314 = (/ &
 & 1.65801E-07_JPRB, 2.01731E-08_JPRB, 2.15003E-08_JPRB, 2.95195E-08_JPRB, 2.95195E-08_JPRB, 2.35479E-08_JPRB &
 & , 6.33757E-09_JPRB,-3.68283E-09_JPRB,-3.58307E-09_JPRB, 1.09913E-08_JPRB,-1.05786E-08_JPRB,-1.16090E-08_JPRB &
 & ,-1.07976E-08_JPRB, 1.84420E-07_JPRB /)  

!     ------------------------------------------------------------------

! LW : absorption coefficient as a function of effective radius in RRTM

! water clouds from Savijarvi

RHSAVI( 1, :) = (/  0.1651082_JPRB  , -0.003494839_JPRB  , 0.0_JPRB            /)  
RHSAVI( 2, :) = (/  0.327820597_JPRB, -0.0219634383_JPRB , 0.000506783898_JPRB /)   
RHSAVI( 3, :) = (/  0.504805453_JPRB, -0.0478602354_JPRB , 0.00141521102_JPRB  /)  
RHSAVI( 4, :) = (/  0.513169093_JPRB, -0.0508960145_JPRB , 0.00155498711_JPRB  /)  
RHSAVI( 5, :) = (/  0.448042082_JPRB, -0.0431857592_JPRB , 0.00130848978_JPRB  /)  
RHSAVI( 6, :) = (/  0.249547237_JPRB, -0.0185273835_JPRB , 0.00050361258_JPRB  /)  
RHSAVI( 7, :) = (/  0.135486796_JPRB, -0.00585852322_JPRB, 0.00011873119_JPRB  /)  
RHSAVI( 8, :) = (/  0.126668819_JPRB, -0.00497949082_JPRB, 0.0000927679172_JPRB/)  
RHSAVI( 9, :) = (/  0.130938752_JPRB, -0.00531814674_JPRB, 0.0000984953029_JPRB/)  
RHSAVI(10, :) = (/  0.147024519_JPRB, -0.00671655774_JPRB, 0.000130885091_JPRB /)  
RHSAVI(11, :) = (/  0.241463914_JPRB, -0.0174486461_JPRB , 0.000457756556_JPRB /)  
RHSAVI(12, :) = (/  0.07576579_JPRB , -0.001695588_JPRB  , 0.0_JPRB            /)  
RHSAVI(13, :) = (/  0.1032178_JPRB  , -0.00293412_JPRB   , 0.0_JPRB            /)  
RHSAVI(14, :) = (/  0.07342832_JPRB , -0.001775135_JPRB  , 0.0_JPRB            /)  
RHSAVI(15, :) = (/  0.04649514_JPRB , -0.0009165462_JPRB , 0.0_JPRB            /)  
RHSAVI(16, :) = (/  0.06893519_JPRB , -0.0001245402_JPRB , 0.0_JPRB            /)  
    
! water clouds from Lindner & Li (2000)

RLILIA( 1, :) = (/ 0.88116E-01_JPRB,-0.12857E-02_JPRB, 0.81658E+00_JPRB,-0.39428E+01_JPRB, 0.46652E+01_JPRB /)        
RLILIA( 2, :) = (/ 0.41307E-03_JPRB,-0.59631E-04_JPRB, 0.24275E+01_JPRB,-0.90838E+01_JPRB, 0.96069E+01_JPRB /)     
RLILIA( 3, :) = (/-0.57709E-01_JPRB, 0.99071E-03_JPRB, 0.31118E+01_JPRB,-0.95540E+01_JPRB, 0.90189E+01_JPRB /)     
RLILIA( 4, :) = (/-0.53069E-01_JPRB, 0.99992E-03_JPRB, 0.28045E+01_JPRB,-0.72836E+01_JPRB, 0.62573E+01_JPRB /)  
RLILIA( 5, :) = (/-0.23627E-01_JPRB, 0.55291E-03_JPRB, 0.21785E+01_JPRB,-0.54664E+01_JPRB, 0.47379E+01_JPRB /)       
RLILIA( 6, :) = (/ 0.29022E-01_JPRB,-0.39657E-03_JPRB, 0.14902E+01_JPRB,-0.50777E+01_JPRB, 0.52170E+01_JPRB /)  
RLILIA( 7, :) = (/-0.24901E-01_JPRB, 0.16195E-03_JPRB, 0.29375E+02_JPRB,-0.11437E+02_JPRB, 0.12273E+02_JPRB /)  
RLILIA( 8, :) = (/-0.14269E+00_JPRB, 0.22282E-02_JPRB, 0.46478E+01_JPRB,-0.16369E+02_JPRB, 0.16533E+02_JPRB /)  
RLILIA( 9, :) = (/-0.20398E+00_JPRB, 0.34708E-02_JPRB, 0.52858E+01_JPRB,-0.16603E+02_JPRB, 0.15392E+02_JPRB /)        
RLILIA(10, :) = (/-0.18318E+00_JPRB, 0.33080E-02_JPRB, 0.46120E+01_JPRB,-0.11550E+02_JPRB, 0.87086E+01_JPRB /)  
RLILIA(11, :) = (/-0.20420E+00_JPRB, 0.37167E-02_JPRB, 0.48566E+01_JPRB,-0.11972E+02_JPRB, 0.86344E+01_JPRB /)        
RLILIA(12, :) = (/-0.14037E+00_JPRB, 0.28058E-02_JPRB, 0.34969E+01_JPRB,-0.33770E+01_JPRB,-0.23541E+01_JPRB /)  
RLILIA(13, :) = (/-0.14037E+00_JPRB, 0.28058E-02_JPRB, 0.34969E+01_JPRB,-0.33770E+01_JPRB,-0.23541E+01_JPRB /)  
RLILIA(14, :) = (/-0.14037E+00_JPRB, 0.28058E-02_JPRB, 0.34969E+01_JPRB,-0.33770E+01_JPRB,-0.23541E+01_JPRB /)  
RLILIA(15, :) = (/-0.14037E+00_JPRB, 0.28058E-02_JPRB, 0.34969E+01_JPRB,-0.33770E+01_JPRB,-0.23541E+01_JPRB /)  
RLILIA(16, :) = (/-0.14037E+00_JPRB, 0.28058E-02_JPRB, 0.34969E+01_JPRB,-0.33770E+01_JPRB,-0.23541E+01_JPRB /)  

RLILIB( 1, :) = (/ 0.76110E+00_JPRB, 0.57048E+00_JPRB,-0.16721E-01_JPRB, 0.27886E-03_JPRB /)
RLILIB( 2, :) = (/ 0.51897E+00_JPRB, 0.86352E+00_JPRB,-0.29099E-02_JPRB, 0.63126E-04_JPRB /)
RLILIB( 3, :) = (/ 0.50964E+00_JPRB, 0.76057E+00_JPRB,-0.15630E-02_JPRB, 0.26245E-04_JPRB /)
RLILIB( 4, :) = (/ 0.55180E+00_JPRB, 0.65812E+00_JPRB,-0.41768E-02_JPRB, 0.64243E-04_JPRB /)
RLILIB( 5, :) = (/ 0.60890E+00_JPRB, 0.59525E+00_JPRB,-0.74557E-02_JPRB, 0.11757E-03_JPRB /)
RLILIB( 6, :) = (/ 0.49787E+00_JPRB, 0.74581E+00_JPRB,-0.37379E-02_JPRB, 0.90555E-04_JPRB /)
RLILIB( 7, :) = (/ 0.46323E-01_JPRB, 0.10131E+01_JPRB, 0.18692E-01_JPRB,-0.21898E-03_JPRB /)
RLILIB( 8, :) = (/-0.43354E-01_JPRB, 0.78454E+00_JPRB, 0.26449E-01_JPRB,-0.36385E-03_JPRB /) 
RLILIB( 9, :) = (/-0.27033E-01_JPRB, 0.54717E+00_JPRB, 0.27903E-01_JPRB,-0.41074E-03_JPRB /)
RLILIB(10, :) = (/ 0.11595E+00_JPRB, 0.30559E+00_JPRB, 0.22727E-01_JPRB,-0.36032E-03_JPRB /)
RLILIB(11, :) = (/-0.14764E-01_JPRB, 0.28331E+00_JPRB, 0.23515E-01_JPRB,-0.33662E-03_JPRB /)
RLILIB(12, :) = (/ 0.68056E-02_JPRB, 0.59837E-01_JPRB, 0.21125E-01_JPRB,-0.29740E-03_JPRB /)
RLILIB(13, :) = (/ 0.68056E-02_JPRB, 0.59837E-01_JPRB, 0.21125E-01_JPRB,-0.29740E-03_JPRB /)
RLILIB(14, :) = (/ 0.68056E-02_JPRB, 0.59837E-01_JPRB, 0.21125E-01_JPRB,-0.29740E-03_JPRB /)
RLILIB(15, :) = (/ 0.68056E-02_JPRB, 0.59837E-01_JPRB, 0.21125E-01_JPRB,-0.29740E-03_JPRB /)
RLILIB(16, :) = (/ 0.68056E-02_JPRB, 0.59837E-01_JPRB, 0.21125E-01_JPRB,-0.29740E-03_JPRB /)

! ice clouds from Fu & Liou (1993)

RFULIO( 1, :) = (/ -7.752E-03_JPRB, 4.624_JPRB, -42.01_JPRB /)  
RFULIO( 2, :) = (/ -1.741E-02_JPRB, 5.541_JPRB, -58.42_JPRB /)  
RFULIO( 3, :) = (/ -1.704E-02_JPRB, 4.830_JPRB,  16.27_JPRB /)  
RFULIO( 4, :) = (/ -1.151E-02_JPRB, 4.182_JPRB,  31.13_JPRB /)  
RFULIO( 5, :) = (/ -1.026E-02_JPRB, 4.105_JPRB,  16.36_JPRB /)  
RFULIO( 6, :) = (/ -8.294E-03_JPRB, 3.925_JPRB,  1.315_JPRB /)  
RFULIO( 7, :) = (/ -1.153E-02_JPRB, 4.109_JPRB,  17.32_JPRB /)  
RFULIO( 8, :) = (/ -9.609E-03_JPRB, 3.768_JPRB,  34.11_JPRB /)  
RFULIO( 9, :) = (/ -9.061E-03_JPRB, 3.741_JPRB,  26.48_JPRB /)  
RFULIO(10, :) = (/ -8.441E-03_JPRB, 3.715_JPRB,  19.48_JPRB /)  
RFULIO(11, :) = (/ -8.088E-03_JPRB, 3.717_JPRB,  17.17_JPRB /)  
RFULIO(12, :) = (/ -8.088E-03_JPRB, 3.717_JPRB,  17.17_JPRB /)  
RFULIO(13, :) = (/ -7.770E-03_JPRB, 3.734_JPRB,  11.85_JPRB /)  
RFULIO(14, :) = (/ -6.656E-03_JPRB, 3.686_JPRB,  0.0_JPRB   /)  
RFULIO(15, :) = (/ -6.656E-03_JPRB, 3.686_JPRB,  0.0_JPRB   /)  
RFULIO(16, :) = (/ -6.656E-03_JPRB, 3.686_JPRB,  0.0_JPRB   /)  
 
! ice clouds from Fu et al. (1998) 

!-- extinction coefficient
RFUETA( 1, :) = (/  4.919685E-03_JPRB, 2.327741E+00_JPRB,-1.390858E+01_JPRB /)  
RFUETA( 2, :) = (/  3.325756E-03_JPRB, 2.601360E+00_JPRB,-1.909602E+01_JPRB /)  
RFUETA( 3, :) = (/ -1.334860E-02_JPRB, 4.043808E+00_JPRB,-2.171029E+01_JPRB /)  
RFUETA( 4, :) = (/ -9.524174E-03_JPRB, 3.587742E+00_JPRB,-1.068895E+01_JPRB /)  
RFUETA( 5, :) = (/ -4.159424E-03_JPRB, 3.047325E+00_JPRB,-5.061568E+00_JPRB /)  
RFUETA( 6, :) = (/ -1.691632E-03_JPRB, 2.765756E+00_JPRB,-8.331033E+00_JPRB /)  
RFUETA( 7, :) = (/ -8.372696E-03_JPRB, 3.455018E+00_JPRB,-1.516692E+01_JPRB /)  
RFUETA( 8, :) = (/ -8.178608E-03_JPRB, 3.401245E+00_JPRB,-8.812820E+00_JPRB /)  
RFUETA( 9, :) = (/ -4.936610E-03_JPRB, 3.087764E+00_JPRB,-3.884262E+00_JPRB /)  
RFUETA(10, :) = (/ -3.034573E-03_JPRB, 2.900043E+00_JPRB,-1.849911E+00_JPRB /)  
RFUETA(11, :) = (/ -3.034573E-03_JPRB, 2.900043E+00_JPRB,-1.849911E+00_JPRB /)  
RFUETA(12, :) = (/ -2.465236E-03_JPRB, 2.833187E+00_JPRB,-4.227573E-01_JPRB /)  
RFUETA(13, :) = (/ -2.308881E-03_JPRB, 2.814002E+00_JPRB, 1.072211E+00_JPRB /)  
RFUETA(14, :) = (/ -2.308881E-03_JPRB, 2.814002E+00_JPRB, 1.072211E+00_JPRB /)  
RFUETA(15, :) = (/ -2.308881E-03_JPRB, 2.814002E+00_JPRB, 1.072211E+00_JPRB /)  
RFUETA(16, :) = (/ -2.308881E-03_JPRB, 2.814002E+00_JPRB, 1.072211E+00_JPRB /)  

!-- absorption coefficient
RFUETB( 1, :) = (/  8.869787E-01_JPRB, 2.118409E-02_JPRB,-2.781429E-04_JPRB, 1.094562E-06_JPRB /)
RFUETB( 2, :) = (/  2.005578E-01_JPRB, 2.132614E-02_JPRB,-1.751052E-04_JPRB, 5.355885E-07_JPRB /)
RFUETB( 3, :) = (/  3.003701E-01_JPRB, 2.051529E-02_JPRB,-1.931684E-04_JPRB, 6.583031E-07_JPRB /)
RFUETB( 4, :) = (/  9.551440E-01_JPRB, 1.309792E-02_JPRB,-1.793694E-04_JPRB, 7.313392E-07_JPRB /)
RFUETB( 5, :) = (/  1.466481E+00_JPRB,-2.129226E-03_JPRB,-1.361630E-05_JPRB, 1.193649E-07_JPRB /)
RFUETB( 6, :) = (/  1.195515E+00_JPRB, 3.350616E-03_JPRB,-5.266996E-05_JPRB, 2.233377E-07_JPRB /)
RFUETB( 7, :) = (/  5.409536E-01_JPRB, 1.949649E-02_JPRB,-2.050908E-04_JPRB, 7.364680E-07_JPRB /)
RFUETB( 8, :) = (/  5.874323E-01_JPRB, 1.876628E-02_JPRB,-2.045834E-04_JPRB, 7.510080E-07_JPRB /)
RFUETB( 9, :) = (/  7.152274E-01_JPRB, 1.621734E-02_JPRB,-1.868544E-04_JPRB, 7.078738E-07_JPRB /)
RFUETB(10, :) = (/  8.862434E-01_JPRB, 1.226538E-02_JPRB,-1.523076E-04_JPRB, 6.000892E-07_JPRB /)
RFUETB(11, :) = (/  8.862434E-01_JPRB, 1.226538E-02_JPRB,-1.523076E-04_JPRB, 6.000892E-07_JPRB /)
RFUETB(12, :) = (/  7.428957E-01_JPRB, 1.279601E-02_JPRB,-1.391803E-04_JPRB, 5.180104E-07_JPRB /)
RFUETB(13, :) = (/  4.346482E-01_JPRB, 1.721457E-02_JPRB,-1.623227E-04_JPRB, 5.561523E-07_JPRB /)
RFUETB(14, :) = (/  4.346482E-01_JPRB, 1.721457E-02_JPRB,-1.623227E-04_JPRB, 5.561523E-07_JPRB /)
RFUETB(15, :) = (/  4.346482E-01_JPRB, 1.721457E-02_JPRB,-1.623227E-04_JPRB, 5.561523E-07_JPRB /)
RFUETB(16, :) = (/  4.346482E-01_JPRB, 1.721457E-02_JPRB,-1.623227E-04_JPRB, 5.561523E-07_JPRB /)

!-- asymmetry factor
RFUETC( 1, :) = (/  4.949276E-01_JPRB, 1.186174E-02_JPRB,-1.267629E-04_JPRB, 4.603574E-07_JPRB /) 
RFUETC( 2, :) = (/  6.891414E-01_JPRB, 6.192281E-03_JPRB,-6.459514E-05_JPRB, 2.436963E-07_JPRB /)
RFUETC( 3, :) = (/  7.260484E-01_JPRB, 2.664334E-03_JPRB,-1.251136E-05_JPRB, 2.243377E-08_JPRB /)
RFUETC( 4, :) = (/  7.363466E-01_JPRB, 4.798266E-03_JPRB,-4.413293E-05_JPRB, 1.525774E-07_JPRB /)
RFUETC( 5, :) = (/  7.984021E-01_JPRB, 3.977117E-03_JPRB,-4.471984E-05_JPRB, 1.694919E-07_JPRB /)
RFUETC( 6, :) = (/  8.663385E-01_JPRB, 2.797934E-03_JPRB,-3.187011E-05_JPRB, 1.217209E-07_JPRB /)
RFUETC( 7, :) = (/  8.906280E-01_JPRB, 1.903269E-03_JPRB,-1.733552E-05_JPRB, 5.855071E-08_JPRB /) 
RFUETC( 8, :) = (/  8.609604E-01_JPRB, 2.200445E-03_JPRB,-1.748105E-05_JPRB, 5.176616E-08_JPRB /)
RFUETC( 9, :) = (/  8.522816E-01_JPRB, 2.523627E-03_JPRB,-2.149196E-05_JPRB, 6.685067E-08_JPRB /)
RFUETC(10, :) = (/  8.741665E-01_JPRB, 2.455409E-03_JPRB,-2.456935E-05_JPRB, 8.641223E-08_JPRB /)
RFUETC(11, :) = (/  8.741665E-01_JPRB, 2.455409E-03_JPRB,-2.456935E-05_JPRB, 8.641223E-08_JPRB /)
RFUETC(12, :) = (/  8.472918E-01_JPRB, 2.559953E-03_JPRB,-2.182660E-05_JPRB, 6.879977E-08_JPRB /)
RFUETC(13, :) = (/  7.962716E-01_JPRB, 3.003488E-03_JPRB,-2.082376E-05_JPRB, 5.366545E-08_JPRB /)
RFUETC(14, :) = (/  7.962716E-01_JPRB, 3.003488E-03_JPRB,-2.082376E-05_JPRB, 5.366545E-08_JPRB /)
RFUETC(15, :) = (/  7.962716E-01_JPRB, 3.003488E-03_JPRB,-2.082376E-05_JPRB, 5.366545E-08_JPRB /)
RFUETC(16, :) = (/  7.962716E-01_JPRB, 3.003488E-03_JPRB,-2.082376E-05_JPRB, 5.366545E-08_JPRB /)

!     ----------------------------------------------------------------
! Ebert-Curry

! LW : spectrally defined for EC-OPE

REBCUI = (/  1.136_JPRB,  1.338_JPRB,  1.166_JPRB,  1.166_JPRB,  1.118_JPRB,  &
 & 0.600_JPRB /)  
REBCUJ = (/ 0.0036_JPRB, 0.0003_JPRB, 0.0016_JPRB, 0.0016_JPRB, 0.0020_JPRB,  &
 & 0.0068_JPRB /)  

! LW : spectrally defined for RRTM
! mass-absorption coefficients for vertical path: no diffusivity factor

REBCUG = (/ 0.718_JPRB, 0.726_JPRB, 1.136_JPRB, 1.320_JPRB, 1.505_JPRB, &
 & 1.290_JPRB, 0.911_JPRB, 0.949_JPRB, 1.021_JPRB, 1.193_JPRB, &
 & 1.279_JPRB, 0.626_JPRB, 0.647_JPRB, 0.668_JPRB, 0.690_JPRB, &
 & 0.690_JPRB /)  

REBCUH = (/ 0.0069_JPRB, 0.0060_JPRB, 0.0024_JPRB, 0.0004_JPRB,-0.0016_JPRB, &
 & 0.0003_JPRB, 0.0043_JPRB, 0.0038_JPRB, 0.0030_JPRB, 0.0013_JPRB, &
 & 0.0005_JPRB, 0.0054_JPRB, 0.0052_JPRB, 0.0050_JPRB, 0.0048_JPRB, &
 & 0.0048_JPRB /)  


!*         2. 

! Ice particle Effective Radius as a function of LWC

REFFIA= 40._JPRB

!     ------------------------------------------------------------------

! SW : absorption coefficients

IF (KSW == 2) THEN
  DO JNU=1,KSW
    RASWCA(JNU)=ZASWCA2(JNU)*1.E-02_JPRB
    RASWCB(JNU)=ZASWCB2(JNU)
    RASWCC(JNU)=ZASWCC2(JNU)
    RASWCD(JNU)=ZASWCD2(JNU)
    RASWCE(JNU)=ZASWCE2(JNU)
    RASWCF(JNU)=ZASWCF2(JNU)*1.E-03_JPRB

    RKPNWCA(JNU)=ZKPNWCA2(JNU)
    RKPNWCB(JNU)=ZKPNWCB2(JNU)
    RKPNWCC(JNU)=ZKPNWCC2(JNU)
    RKPNWCD(JNU)=ZKPNWCD2(JNU)
    RKPNWCE(JNU)=ZKPNWCE2(JNU)
    RKPNWCF(JNU)=ZKPNWCF2(JNU)
    RKPNWCG(JNU)=ZKPNWCG2(JNU)
    RKPNWCH(JNU)=ZKPNWCH2(JNU)

    REBCUA(JNU)=ZEBCUA2(JNU)
    REBCUB(JNU)=ZEBCUB2(JNU)
    REBCUC(JNU)=ZEBCUC2(JNU)
    REBCUD(JNU)=ZEBCUD2(JNU)
    REBCUE(JNU)=ZEBCUE2(JNU)
    REBCUF(JNU)=ZEBCUF2(JNU)

    RYFWCA(JNU)=ZYFWCA2(JNU)
    RYFWCB(JNU)=ZYFWCB2(JNU)
    RYFWCC(JNU)=ZYFWCC2(JNU)
    RYFWCD(JNU)=ZYFWCD2(JNU)
    RYFWCE(JNU)=ZYFWCE2(JNU)
    RYFWCF(JNU)=ZYFWCF2(JNU)

    RFLAA0(JNU)=ZFLAA02(JNU)
    RFLAA1(JNU)=ZFLAA12(JNU)
    RFLBB0(JNU)=ZFLBB02(JNU)
    RFLBB1(JNU)=ZFLBB12(JNU)
    RFLBB2(JNU)=ZFLBB22(JNU)
    RFLBB3(JNU)=ZFLBB32(JNU)
    RFLCC0(JNU)=ZFLCC02(JNU)
    RFLCC1(JNU)=ZFLCC12(JNU)
    RFLCC2(JNU)=ZFLCC22(JNU)
    RFLCC3(JNU)=ZFLCC32(JNU)
    RFLDD0(JNU)=ZFLDD02(JNU)
    RFLDD1(JNU)=ZFLDD12(JNU)
    RFLDD2(JNU)=ZFLDD22(JNU)
    RFLDD3(JNU)=ZFLDD32(JNU)

! Fu (1996)'s optical properties are not available for 2 spectral intervals
!    RFUAA0(JNU)=ZFUAA02(JNU)
!    RFUAA1(JNU)=ZFUAA12(JNU)
!    RFUBB0(JNU)=ZFUBB02(JNU)
!    RFUBB1(JNU)=ZFUBB12(JNU)
!    RFUBB2(JNU)=ZFUBB22(JNU)
!    RFUBB3(JNU)=ZFUBB32(JNU)
!    RFUCC0(JNU)=ZFUCC02(JNU)
!    RFUCC1(JNU)=ZFUCC12(JNU)
!    RFUCC2(JNU)=ZFUCC22(JNU)
!    RFUCC3(JNU)=ZFUCC32(JNU)
    
  ENDDO
ELSEIF (KSW == 4) THEN
  DO JNU=1,KSW
    RASWCA(JNU)=ZASWCA4(JNU)*1.E-02_JPRB
    RASWCB(JNU)=ZASWCB4(JNU)
    RASWCC(JNU)=ZASWCC4(JNU)
    RASWCD(JNU)=ZASWCD4(JNU)
    RASWCE(JNU)=ZASWCE4(JNU)
    RASWCF(JNU)=ZASWCF4(JNU)*1.E-03_JPRB

    RKPNWCA(JNU)=ZKPNWCA4(JNU)
    RKPNWCB(JNU)=ZKPNWCB4(JNU)
    RKPNWCC(JNU)=ZKPNWCC4(JNU)
    RKPNWCD(JNU)=ZKPNWCD4(JNU)
    RKPNWCE(JNU)=ZKPNWCE4(JNU)
    RKPNWCF(JNU)=ZKPNWCF4(JNU)
    RKPNWCG(JNU)=ZKPNWCG4(JNU)
    RKPNWCH(JNU)=ZKPNWCH4(JNU)

    REBCUA(JNU)=ZEBCUA4(JNU)
    REBCUB(JNU)=ZEBCUB4(JNU)
    REBCUC(JNU)=ZEBCUC4(JNU)
    REBCUD(JNU)=ZEBCUD4(JNU)
    REBCUE(JNU)=ZEBCUE4(JNU)
    REBCUF(JNU)=ZEBCUF4(JNU)

    RYFWCA(JNU)=ZYFWCA4(JNU)
    RYFWCB(JNU)=ZYFWCB4(JNU)
    RYFWCC(JNU)=ZYFWCC4(JNU)
    RYFWCD(JNU)=ZYFWCD4(JNU)
    RYFWCE(JNU)=ZYFWCE4(JNU)
    RYFWCF(JNU)=ZYFWCF4(JNU)

    RFLAA0(JNU)=ZFLAA04(JNU)
    RFLAA1(JNU)=ZFLAA14(JNU)
    RFLBB0(JNU)=ZFLBB04(JNU)
    RFLBB1(JNU)=ZFLBB14(JNU)
    RFLBB2(JNU)=ZFLBB24(JNU)
    RFLBB3(JNU)=ZFLBB34(JNU)
    RFLCC0(JNU)=ZFLCC04(JNU)
    RFLCC1(JNU)=ZFLCC14(JNU)
    RFLCC2(JNU)=ZFLCC24(JNU)
    RFLCC3(JNU)=ZFLCC34(JNU)
    RFLDD0(JNU)=ZFLDD04(JNU)
    RFLDD1(JNU)=ZFLDD14(JNU)
    RFLDD2(JNU)=ZFLDD24(JNU)
    RFLDD3(JNU)=ZFLDD34(JNU)
    
    RFUAA0(JNU)=ZFUAA04(JNU)
    RFUAA1(JNU)=ZFUAA14(JNU)
    RFUBB0(JNU)=ZFUBB04(JNU)
    RFUBB1(JNU)=ZFUBB14(JNU)
    RFUBB2(JNU)=ZFUBB24(JNU)
    RFUBB3(JNU)=ZFUBB34(JNU)
    RFUCC0(JNU)=ZFUCC04(JNU)
    RFUCC1(JNU)=ZFUCC14(JNU)
    RFUCC2(JNU)=ZFUCC24(JNU)
    RFUCC3(JNU)=ZFUCC34(JNU)
    
  ENDDO
ELSEIF (KSW == 6) THEN
  DO JNU=1,KSW
    RASWCA(JNU)=ZASWCA6(JNU)*1.E-02_JPRB
    RASWCB(JNU)=ZASWCB6(JNU)
    RASWCC(JNU)=ZASWCC6(JNU)
    RASWCD(JNU)=ZASWCD6(JNU)
    RASWCE(JNU)=ZASWCE6(JNU)
    RASWCF(JNU)=ZASWCF6(JNU)*1.E-03_JPRB

    RKPNWCA(JNU)=ZKPNWCA6(JNU)
    RKPNWCB(JNU)=ZKPNWCB6(JNU)
    RKPNWCC(JNU)=ZKPNWCC6(JNU)
    RKPNWCD(JNU)=ZKPNWCD6(JNU)
    RKPNWCE(JNU)=ZKPNWCE6(JNU)
    RKPNWCF(JNU)=ZKPNWCF6(JNU)
    RKPNWCG(JNU)=ZKPNWCG6(JNU)
    RKPNWCH(JNU)=ZKPNWCH6(JNU)

    REBCUA(JNU)=ZEBCUA6(JNU)
    REBCUB(JNU)=ZEBCUB6(JNU)
    REBCUC(JNU)=ZEBCUC6(JNU)
    REBCUD(JNU)=ZEBCUD6(JNU)
    REBCUE(JNU)=ZEBCUE6(JNU)
    REBCUF(JNU)=ZEBCUF6(JNU)

    RYFWCA(JNU)=ZYFWCA6(JNU)
    RYFWCB(JNU)=ZYFWCB6(JNU)
    RYFWCC(JNU)=ZYFWCC6(JNU)
    RYFWCD(JNU)=ZYFWCD6(JNU)
    RYFWCE(JNU)=ZYFWCE6(JNU)
    RYFWCF(JNU)=ZYFWCF6(JNU)

    RFLAA0(JNU)=ZFLAA06(JNU)
    RFLAA1(JNU)=ZFLAA16(JNU)
    RFLBB0(JNU)=ZFLBB06(JNU)
    RFLBB1(JNU)=ZFLBB16(JNU)
    RFLBB2(JNU)=ZFLBB26(JNU)
    RFLBB3(JNU)=ZFLBB36(JNU)
    RFLCC0(JNU)=ZFLCC06(JNU)
    RFLCC1(JNU)=ZFLCC16(JNU)
    RFLCC2(JNU)=ZFLCC26(JNU)
    RFLCC3(JNU)=ZFLCC36(JNU)
    RFLDD0(JNU)=ZFLDD06(JNU)
    RFLDD1(JNU)=ZFLDD16(JNU)
    RFLDD2(JNU)=ZFLDD26(JNU)
    RFLDD3(JNU)=ZFLDD36(JNU)
    
    RFUAA0(JNU)=ZFUAA06(JNU)
    RFUAA1(JNU)=ZFUAA16(JNU)
    RFUBB0(JNU)=ZFUBB06(JNU)
    RFUBB1(JNU)=ZFUBB16(JNU)
    RFUBB2(JNU)=ZFUBB26(JNU)
    RFUBB3(JNU)=ZFUBB36(JNU)
    RFUCC0(JNU)=ZFUCC06(JNU)
    RFUCC1(JNU)=ZFUCC16(JNU)
    RFUCC2(JNU)=ZFUCC26(JNU)
    RFUCC3(JNU)=ZFUCC36(JNU)
    
  ENDDO
ELSEIF (KTSW == 14 .AND. LSRTM) THEN
    WRITE (NULOUT,*) 'SUCLOPN: 14-SPECTRAL INTERVALS --> RRTM_SW'
!  DO JNU=1,KSW
!    RASWCA(JNU)=ZASWCA14(JNU)*1.E-02_JPRB
!    RASWCB(JNU)=ZASWCB14(JNU)
!    RASWCC(JNU)=ZASWCC14(JNU)
!    RASWCD(JNU)=ZASWCD14(JNU)
!    RASWCE(JNU)=ZASWCE14(JNU)
!    RASWCF(JNU)=ZASWCF14(JNU)*1.E-03_JPRB

!    REBCUA(JNU)=ZEBCUA14(JNU)
!    REBCUB(JNU)=ZEBCUB14(JNU)
!    REBCUC(JNU)=ZEBCUC14(JNU)
!    REBCUD(JNU)=ZEBCUD14(JNU)
!    REBCUE(JNU)=ZEBCUE14(JNU)
!    REBCUF(JNU)=ZEBCUF14(JNU)

!    RYFWCA(JNU)=ZYFWCA14(JNU)
!    RYFWCB(JNU)=ZYFWCB14(JNU)
!    RYFWCC(JNU)=ZYFWCC14(JNU)
!    RYFWCD(JNU)=ZYFWCD14(JNU)
!    RYFWCE(JNU)=ZYFWCE14(JNU)
!    RYFWCF(JNU)=ZYFWCF14(JNU)

!    RFLAA0(JNU)=ZFLAA014(JNU)
!    RFLAA1(JNU)=ZFLAA114(JNU)
!    RFLBB0(JNU)=ZFLBB014(JNU)
!    RFLBB1(JNU)=ZFLBB114(JNU)
!    RFLBB2(JNU)=ZFLBB214(JNU)
!    RFLBB3(JNU)=ZFLBB314(JNU)
!    RFLCC0(JNU)=ZFLCC014(JNU)
!    RFLCC1(JNU)=ZFLCC114(JNU)
!    RFLCC2(JNU)=ZFLCC214(JNU)
!    RFLCC3(JNU)=ZFLCC314(JNU)
!    RFLDD0(JNU)=ZFLDD014(JNU)
!    RFLDD1(JNU)=ZFLDD114(JNU)
!    RFLDD2(JNU)=ZFLDD214(JNU)
!    RFLDD3(JNU)=ZFLDD314(JNU)

!    RFUAA0(JNU)=ZFUAA014(JNU)
!    RFUAA1(JNU)=ZFUAA114(JNU)
!    RFUBB0(JNU)=ZFUBB014(JNU)
!    RFUBB1(JNU)=ZFUBB114(JNU)
!    RFUBB2(JNU)=ZFUBB214(JNU)
!    RFUBB3(JNU)=ZFUBB314(JNU)
!    RFUCC0(JNU)=ZFUCC014(JNU)
!    RFUCC1(JNU)=ZFUCC114(JNU)
!    RFUCC2(JNU)=ZFUCC214(JNU)
!    RFUCC3(JNU)=ZFUCC314(JNU)
!!    RFUDD0(JNU)=ZFUDD014(JNU)
!!    RFUDD1(JNU)=ZFUDD114(JNU)
!!    RFUDD2(JNU)=ZFUDD214(JNU)
!!    RFUDD3(JNU)=ZFUDD314(JNU)

!    WRITE (NULOUT,*) 'SUCLOPN: 14-SPECTRAL INTERVALS --> RRTM_SW'     
!  ENDDO
ELSE
  CALL ABOR1('SUCLOPN: WRONG SW SPECTRAL RESOLUTION')
ENDIF
!     ------------------------------------------------------------------

!*          2.    CLOUD OVERLAP PARAMETERS
!                 ------------------------

ZAOVLP = (/ -2.250E-05_JPRB,-7.316E-06_JPRB,-1.966E-05_JPRB /) 
ZBOVLP = (/  0.7865_JPRB   , 0.8186_JPRB   , 0.8900_JPRB    /)

IF (KLEV == 19) THEN
  RAOVLP=ZAOVLP(1)
  RBOVLP=ZBOVLP(1)
ELSEIF (KLEV == 31) THEN
  RAOVLP=ZAOVLP(2)
  RBOVLP=ZBOVLP(2)
ELSEIF (KLEV == 60) THEN
  RAOVLP=ZAOVLP(3)
  RBOVLP=ZBOVLP(3)
ELSE  
  RAOVLP=ZAOVLP(3)
  RBOVLP=ZBOVLP(3)
ENDIF  

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUCLOPN',1,ZHOOK_HANDLE)
END SUBROUTINE SUCLOPN
