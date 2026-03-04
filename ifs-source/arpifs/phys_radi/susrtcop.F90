! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUSRTCOP

!**** *SUSRTCOP*  - INITIALIZE COMMON YOESRTCOP

!     PURPOSE.
!     --------
!           INITIALIZE YOESRTCOP, WITH CLOUD OPTICAL PARAMETERS

!**   INTERFACE.
!     ----------
!        *CALL*  SUSRTCOP
!        FROM *SUECRAD*

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOESRTCOP

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

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 03-03-06
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOESRTCOP , ONLY : &
 & RSYFWA   ,RSYFWB   ,RSYFWC   ,RSYFWD   ,RSYFWE   ,RSYFWF  &
 & , RSECIA   ,RSECIB   ,RSECIC   ,RSECID   ,RSECIE   ,RSECIF &
 & , RSASWA   ,RSASWB   ,RSASWC   ,RSASWD   ,RSASWE   ,RSASWF &
 & , RSFUA0   ,RSFUA1   ,RSFUB0   ,RSFUB1   ,RSFUB2   ,RSFUB3 &
 & , RSFUC0   ,RSFUC1   ,RSFUC2   ,RSFUC3   &
 & , RSFLA0   ,RSFLA1   ,RSFLB0   ,RSFLB1   ,RSFLB2   ,RSFLB3 &
 & , RSFLC0   ,RSFLC1   ,RSFLC2   ,RSFLC3   ,RSFLD0   ,RSFLD1 &
 & , RSFLD2   ,RSFLD3  

IMPLICIT NONE

INTEGER(KIND=JPIM) ::  JNU

!     -----------------------------------------------------------------

REAL(KIND=JPRB) :: ZEBCUA14(14)  ,ZEBCUB14(14)  ,ZEBCUC14(14)  ,ZEBCUD14(14)&
 & ,  ZEBCUE14(14)  ,ZEBCUF14(14)  ,ZYFWCA14(14)  ,ZYFWCB14(14)&
 & ,  ZYFWCC14(14)  ,ZYFWCD14(14)  ,ZYFWCE14(14)  ,ZYFWCF14(14)&
 & ,  ZASWCA14(14)  ,ZASWCB14(14)  ,ZASWCC14(14)  ,ZASWCD14(14)&
 & ,  ZASWCE14(14)  ,ZASWCF14(14)&
 & ,  ZFLAA014(14)  ,ZFLAA114(14)  ,ZFLBB014(14)  ,ZFLBB114(14)&
 & ,  ZFLBB214(14)  ,ZFLBB314(14)  ,ZFLCC014(14)  ,ZFLCC114(14)&
 & ,  ZFLCC214(14)  ,ZFLCC314(14)  ,ZFLDD014(14)  ,ZFLDD114(14)&
 & ,  ZFLDD214(14)  ,ZFLDD314(14)&
 & ,  ZFUAA014(14)  ,ZFUAA114(14)  ,ZFUBB014(14)  ,ZFUBB114(14)&
 & ,  ZFUBB214(14)  ,ZFUBB314(14)  ,ZFUCC014(14)  ,ZFUCC114(14)&
 & ,  ZFUCC214(14)  ,ZFUCC314(14)  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     -----------------------------------------------------------------

!*          1.    SHORTWAVE CLOUD OPTICAL PROPERTIES
!                 ----------------------------------

!     ------------------------------------------------------------------

!*          1.1    FOURTEEN SPECTRAL INTERVALS FOR RRTM_SW
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

IF (LHOOK) CALL DR_HOOK('SUSRTCOP',0,ZHOOK_HANDLE)
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
!     ------------------------------------------------------------------

! SW : absorption coefficients

DO JNU=1,14
  RSASWA(JNU)=ZASWCA14(JNU)*1.E-02_JPRB
  RSASWB(JNU)=ZASWCB14(JNU)
  RSASWC(JNU)=ZASWCC14(JNU)
  RSASWD(JNU)=ZASWCD14(JNU)
  RSASWE(JNU)=ZASWCE14(JNU)
  RSASWF(JNU)=ZASWCF14(JNU)*1.E-03_JPRB

  RSECIA(JNU)=ZEBCUA14(JNU)
  RSECIB(JNU)=ZEBCUB14(JNU)
  RSECIC(JNU)=ZEBCUC14(JNU)
  RSECID(JNU)=ZEBCUD14(JNU)
  RSECIE(JNU)=ZEBCUE14(JNU)
  RSECIF(JNU)=ZEBCUF14(JNU)

  RSYFWA(JNU)=ZYFWCA14(JNU)
  RSYFWB(JNU)=ZYFWCB14(JNU)
  RSYFWC(JNU)=ZYFWCC14(JNU)
  RSYFWD(JNU)=ZYFWCD14(JNU)
  RSYFWE(JNU)=ZYFWCE14(JNU)
  RSYFWF(JNU)=ZYFWCF14(JNU)

  RSFLA0(JNU)=ZFLAA014(JNU)
  RSFLA1(JNU)=ZFLAA114(JNU)
  RSFLB0(JNU)=ZFLBB014(JNU)
  RSFLB1(JNU)=ZFLBB114(JNU)
  RSFLB2(JNU)=ZFLBB214(JNU)
  RSFLB3(JNU)=ZFLBB314(JNU)
  RSFLC0(JNU)=ZFLCC014(JNU)
  RSFLC1(JNU)=ZFLCC114(JNU)
  RSFLC2(JNU)=ZFLCC214(JNU)
  RSFLC3(JNU)=ZFLCC314(JNU)
  RSFLD0(JNU)=ZFLDD014(JNU)
  RSFLD1(JNU)=ZFLDD114(JNU)
  RSFLD2(JNU)=ZFLDD214(JNU)
  RSFLD3(JNU)=ZFLDD314(JNU)
    
  RSFUA0(JNU)=ZFUAA014(JNU)
  RSFUA1(JNU)=ZFUAA114(JNU)
  RSFUB0(JNU)=ZFUBB014(JNU)
  RSFUB1(JNU)=ZFUBB114(JNU)
  RSFUB2(JNU)=ZFUBB214(JNU)
  RSFUB3(JNU)=ZFUBB314(JNU)
  RSFUC0(JNU)=ZFUCC014(JNU)
  RSFUC1(JNU)=ZFUCC114(JNU)
  RSFUC2(JNU)=ZFUCC214(JNU)
  RSFUC3(JNU)=ZFUCC314(JNU)
    
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUSRTCOP',1,ZHOOK_HANDLE)
END SUBROUTINE SUSRTCOP

