! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE CIGETDEAC

! ----------------------------------------------------------------------

!**** *CIGETDEAC* - DEFINE THE DIMENSIONLESS ENERGY ATTENUATION COEFFICIENT 

!*    PURPOSE.
!     --------

!       CIGETDEAC DEFINES THE DIMENSIONLESS ENERGY ATTENUATION COEFFICIENT

!**   INTERFACE.
!     ----------

!       *CALL* *CIGETDEAC

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCES.                                                       
!     -----------  

!     KOHOUT AND MEYLAN, 2008: JGR, 113, doi:10.1029/2007JC004424

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWICE   , ONLY : NICT    ,NICH    ,TICMIN  ,HICMIN   ,       &
     &              DTIC   ,DHIC    ,CIDEAC
      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: IT, IH
      REAL(KIND=JWRB) :: DHI, DCI
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CIGETDEAC',0,ZHOOK_HANDLE)

! I have kindly received the values from Figure 6 of Kohout and Meylan
! from Alison Kohout
! I have also extrapolated the values for thickness=0.2 and
! for periods between 1 and 5 sec by assuming that the values asymptotes to -3 to -1
! for 1 sec. linear dependence on ice thickness

!     ice thickness discretisation
      NICH=36
      HICMIN=0.2_JWRB
      DHIC=0.1_JWRB
      
!     wave period discretisation
      NICT=16
      TICMIN=1.0_JWRB
      DTIC=1.0_JWRB

      IF(ALLOCATED(CIDEAC)) DEALLOCATE(CIDEAC)
      ALLOCATE(CIDEAC(NICT,NICH))

!     assumption:
      CIDEAC( 1, 1)=-2.00_JWRB
      CIDEAC( 1,NICH)=-1.00_JWRB
      DHI=CIDEAC( 1,NICH)-CIDEAC( 1, 1)
      DO IH=2,NICH-1
        CIDEAC( 1,IH)=CIDEAC( 1, 1)+(IH-1)*DHI/(NICH-1)
      ENDDO

      !!!! all values for h=0.20 and 0.30 were extrapolated
! h =    0.2000000    
      CIDEAC( 6, 1)=  -5.50_JWRB
      CIDEAC( 7, 1)=  -7.20_JWRB
      CIDEAC( 8, 1)=  -8.90_JWRB
      CIDEAC( 9, 1)=  -10.10_JWRB
      CIDEAC(10, 1)=  -10.75_JWRB
      CIDEAC(11, 1)=  -11.10_JWRB
      CIDEAC(12, 1)=  -11.30_JWRB
      CIDEAC(13, 1)=  -11.50_JWRB
      CIDEAC(14, 1)= -11.70_JWRB
      CIDEAC(15, 1)= -12.35_JWRB
      CIDEAC(16, 1)= -13.00_JWRB

! h =    0.3000000    
      CIDEAC( 6, 2)=  -4.84_JWRB
      CIDEAC( 7, 2)=  -6.35_JWRB
      CIDEAC( 8, 2)=  -7.99_JWRB
      CIDEAC( 9, 2)=  -9.46_JWRB
      CIDEAC(10, 2)=  -10.43_JWRB
      CIDEAC(11, 2)=  -10.93_JWRB
      CIDEAC(12, 2)=  -11.15_JWRB
      CIDEAC(13, 2)=  -11.37_JWRB
      CIDEAC(14, 2)= -11.59_JWRB
      CIDEAC(15, 2)= -12.00_JWRB
      CIDEAC(16, 2)= -13.44_JWRB

! h =    0.4000000    
      CIDEAC( 6, 3)=   -4.179189_JWRB    
      CIDEAC( 7, 3)=   -5.494833_JWRB
      CIDEAC( 8, 3)=   -7.073360_JWRB
      CIDEAC( 9, 3)=   -8.809628_JWRB
      CIDEAC(10, 3)=   -10.11117_JWRB
      CIDEAC(11, 3)=   -10.75490_JWRB
      CIDEAC(12, 3)=   -11.01079_JWRB
      CIDEAC(13, 3)=   -11.23960_JWRB
      CIDEAC(14, 3)=   -11.46686_JWRB
      CIDEAC(15, 3)=   -11.64766_JWRB
      CIDEAC(16, 3)=   -11.88710_JWRB
 
! h =    0.5000000    
      CIDEAC( 6, 4)=   -3.671164_JWRB
      CIDEAC( 7, 4)=   -4.688528_JWRB
      CIDEAC( 8, 4)=   -6.063766_JWRB
      CIDEAC( 9, 4)=   -7.556777_JWRB
      CIDEAC(10, 4)=   -9.092730_JWRB
      CIDEAC(11, 4)=   -10.21499_JWRB
      CIDEAC(12, 4)=   -10.62448_JWRB
      CIDEAC(13, 4)=   -10.84617_JWRB
      CIDEAC(14, 4)=   -11.03255_JWRB
      CIDEAC(15, 4)=   -11.24719_JWRB
      CIDEAC(16, 4)=   -11.49634_JWRB
 
! h =    0.6000000    
      CIDEAC( 6, 5)=   -3.279724_JWRB
      CIDEAC( 7, 5)=   -4.138198_JWRB
      CIDEAC( 8, 5)=   -5.303435_JWRB
      CIDEAC( 9, 5)=   -6.634662_JWRB
      CIDEAC(10, 5)=   -8.086823_JWRB
      CIDEAC(11, 5)=   -9.434999_JWRB
      CIDEAC(12, 5)=   -10.24225_JWRB
      CIDEAC(13, 5)=   -10.53463_JWRB
      CIDEAC(14, 5)=   -10.77608_JWRB
      CIDEAC(15, 5)=   -10.94032_JWRB
      CIDEAC(16, 5)=   -11.14464_JWRB
 
! h =    0.7000000    
      CIDEAC( 6, 6)=   -2.986051_JWRB
      CIDEAC( 7, 6)=   -3.809261_JWRB
      CIDEAC( 8, 6)=   -4.791394_JWRB
      CIDEAC( 9, 6)=   -5.896924_JWRB
      CIDEAC(10, 6)=   -7.261892_JWRB
      CIDEAC(11, 6)=   -8.680384_JWRB
      CIDEAC(12, 6)=   -9.695383_JWRB
      CIDEAC(13, 6)=   -10.25629_JWRB
      CIDEAC(14, 6)=   -10.57562_JWRB
      CIDEAC(15, 6)=   -10.75792_JWRB
      CIDEAC(16, 6)=   -10.88391_JWRB
 
! h =    0.8000000    
      CIDEAC( 6, 7)=   -2.802249_JWRB
      CIDEAC( 7, 7)=   -3.473820_JWRB
      CIDEAC( 8, 7)=   -4.312738_JWRB
      CIDEAC( 9, 7)=   -5.404909_JWRB
      CIDEAC(10, 7)=   -6.585416_JWRB
      CIDEAC(11, 7)=   -7.936078_JWRB
      CIDEAC(12, 7)=   -9.055226_JWRB
      CIDEAC(13, 7)=   -9.891881_JWRB
      CIDEAC(14, 7)=   -10.28646_JWRB
      CIDEAC(15, 7)=   -10.50197_JWRB
      CIDEAC(16, 7)=   -10.69553_JWRB

! h =    0.9000000    
      CIDEAC( 6, 8)=   -2.649911_JWRB
      CIDEAC( 7, 8)=   -3.208950_JWRB
      CIDEAC( 8, 8)=   -3.984321_JWRB
      CIDEAC( 9, 8)=   -4.919455_JWRB
      CIDEAC(10, 8)=   -6.056835_JWRB
      CIDEAC(11, 8)=   -7.323155_JWRB
      CIDEAC(12, 8)=   -8.486164_JWRB
      CIDEAC(13, 8)=   -9.498883_JWRB
      CIDEAC(14, 8)=   -10.05186_JWRB
      CIDEAC(15, 8)=   -10.26166_JWRB
      CIDEAC(16, 8)=   -10.44009_JWRB
 
! h =     1.000000    
      CIDEAC( 6, 9)=   -2.494633_JWRB
      CIDEAC( 7, 9)=   -3.058290_JWRB
      CIDEAC( 8, 9)=   -3.697340_JWRB
      CIDEAC( 9, 9)=   -4.548424_JWRB
      CIDEAC(10, 9)=   -5.556311_JWRB
      CIDEAC(11, 9)=   -6.724046_JWRB
      CIDEAC(12, 9)=   -7.960689_JWRB
      CIDEAC(13, 9)=   -8.995591_JWRB
      CIDEAC(14, 9)=   -9.717956_JWRB
      CIDEAC(15, 9)=   -10.12654_JWRB
      CIDEAC(16, 9)=   -10.34523_JWRB
 
! h =     1.100000    
      CIDEAC( 6,10)=   -2.399249_JWRB
      CIDEAC( 7,10)=   -2.882891_JWRB
      CIDEAC( 8,10)=   -3.501052_JWRB
      CIDEAC( 9,10)=   -4.253071_JWRB
      CIDEAC(10,10)=   -5.189890_JWRB
      CIDEAC(11,10)=   -6.225914_JWRB
      CIDEAC(12,10)=   -7.413078_JWRB
      CIDEAC(13,10)=   -8.483686_JWRB
      CIDEAC(14,10)=   -9.360973_JWRB
      CIDEAC(15,10)=   -9.923790_JWRB
      CIDEAC(16,10)=   -10.21649_JWRB
 
! h =     1.200000    
      CIDEAC( 6,11)=   -2.314447_JWRB
      CIDEAC( 7,11)=   -2.792313_JWRB
      CIDEAC( 8,11)=   -3.336744_JWRB
      CIDEAC( 9,11)=   -4.020543_JWRB
      CIDEAC(10,11)=   -4.842176_JWRB
      CIDEAC(11,11)=   -5.899529_JWRB
      CIDEAC(12,11)=   -6.954353_JWRB
      CIDEAC(13,11)=   -8.057040_JWRB
      CIDEAC(14,11)=   -8.976917_JWRB
      CIDEAC(15,11)=   -9.661180_JWRB
      CIDEAC(16,11)=   -10.01354_JWRB
 
! h =     1.300000    
      CIDEAC( 6,12)=   -2.243207_JWRB
      CIDEAC( 7,12)=   -2.647049_JWRB
      CIDEAC( 8,12)=   -3.217140_JWRB
      CIDEAC( 9,12)=   -3.814464_JWRB
      CIDEAC(10,12)=   -4.590047_JWRB
      CIDEAC(11,12)=   -5.466197_JWRB
      CIDEAC(12,12)=   -6.520925_JWRB
      CIDEAC(13,12)=   -7.642193_JWRB
      CIDEAC(14,12)=   -8.617725_JWRB
      CIDEAC(15,12)=   -9.381130_JWRB
      CIDEAC(16,12)=   -9.848615_JWRB
 
! h =     1.400000    
      CIDEAC( 6,13)=   -2.169400_JWRB
      CIDEAC( 7,13)=   -2.576199_JWRB
      CIDEAC( 8,13)=   -3.010994_JWRB
      CIDEAC( 9,13)=   -3.624323_JWRB
      CIDEAC(10,13)=   -4.324925_JWRB
      CIDEAC(11,13)=   -5.187501_JWRB
      CIDEAC(12,13)=   -6.182012_JWRB
      CIDEAC(13,13)=   -7.190685_JWRB
      CIDEAC(14,13)=   -8.241459_JWRB
      CIDEAC(15,13)=   -9.042115_JWRB
      CIDEAC(16,13)=   -9.622642_JWRB

! h =     1.500000    
      CIDEAC( 6,14)=   -2.062091_JWRB
      CIDEAC( 7,14)=   -2.452876_JWRB
      CIDEAC( 8,14)=   -2.952611_JWRB
      CIDEAC( 9,14)=   -3.506958_JWRB
      CIDEAC(10,14)=   -4.177629_JWRB
      CIDEAC(11,14)=   -4.920455_JWRB
      CIDEAC(12,14)=   -5.839159_JWRB
      CIDEAC(13,14)=   -6.873662_JWRB
      CIDEAC(14,14)=   -7.827965_JWRB
      CIDEAC(15,14)=   -8.873965_JWRB
      CIDEAC(16,14)=   -9.413249_JWRB

! h =     1.600000    
      CIDEAC( 6,15)=   -1.996322_JWRB
      CIDEAC( 7,15)=   -2.392122_JWRB
      CIDEAC( 8,15)=   -2.823834_JWRB
      CIDEAC( 9,15)=   -3.359782_JWRB
      CIDEAC(10,15)=   -3.958120_JWRB
      CIDEAC(11,15)=   -4.667717_JWRB
      CIDEAC(12,15)=   -5.541727_JWRB
      CIDEAC(13,15)=   -6.549209_JWRB
      CIDEAC(14,15)=   -7.520937_JWRB
      CIDEAC(15,15)=   -8.429287_JWRB
      CIDEAC(16,15)=   -9.221933_JWRB
 
! h =     1.700000    
      CIDEAC( 6,16)=   -1.955446_JWRB
      CIDEAC( 7,16)=   -2.320779_JWRB
      CIDEAC( 8,16)=   -2.715107_JWRB
      CIDEAC( 9,16)=   -3.211328_JWRB
      CIDEAC(10,16)=   -3.808961_JWRB
      CIDEAC(11,16)=   -4.510885_JWRB
      CIDEAC(12,16)=   -5.320290_JWRB
      CIDEAC(13,16)=   -6.226068_JWRB
      CIDEAC(14,16)=   -7.185133_JWRB
      CIDEAC(15,16)=   -8.094091_JWRB
      CIDEAC(16,16)=   -8.907555_JWRB

! h =     1.800000    
      CIDEAC( 6,17)=   -1.919889_JWRB
      CIDEAC( 7,17)=   -2.260443_JWRB
      CIDEAC( 8,17)=   -2.674223_JWRB
      CIDEAC( 9,17)=   -3.108331_JWRB
      CIDEAC(10,17)=   -3.647704_JWRB
      CIDEAC(11,17)=   -4.258914_JWRB
      CIDEAC(12,17)=   -5.110578_JWRB
      CIDEAC(13,17)=   -5.947557_JWRB
      CIDEAC(14,17)=   -6.857039_JWRB
      CIDEAC(15,17)=   -7.801445_JWRB
      CIDEAC(16,17)=   -8.697256_JWRB
 
! h =     1.900000    
      CIDEAC( 6,18)=   -1.866967_JWRB
      CIDEAC( 7,18)=   -2.277320_JWRB
      CIDEAC( 8,18)=   -2.620277_JWRB
      CIDEAC( 9,18)=   -3.053688_JWRB
      CIDEAC(10,18)=   -3.530046_JWRB
      CIDEAC(11,18)=   -4.170623_JWRB
      CIDEAC(12,18)=   -4.870700_JWRB
      CIDEAC(13,18)=   -5.723269_JWRB
      CIDEAC(14,18)=   -6.651424_JWRB
      CIDEAC(15,18)=   -7.524860_JWRB
      CIDEAC(16,18)=   -8.409969_JWRB
 
! h =     2.000000    
      CIDEAC( 6,19)=   -1.822846_JWRB
      CIDEAC( 7,19)=   -2.088338_JWRB
      CIDEAC( 8,19)=   -2.531257_JWRB
      CIDEAC( 9,19)=   -2.942987_JWRB
      CIDEAC(10,19)=   -3.422226_JWRB
      CIDEAC(11,19)=   -4.036695_JWRB
      CIDEAC(12,19)=   -4.701579_JWRB
      CIDEAC(13,19)=   -5.505970_JWRB
      CIDEAC(14,19)=   -6.318645_JWRB
      CIDEAC(15,19)=   -7.290503_JWRB
      CIDEAC(16,19)=   -8.141773_JWRB

! h =     2.100000    
      CIDEAC( 6,20)=   -1.779950_JWRB
      CIDEAC( 7,20)=   -2.108976_JWRB
      CIDEAC( 8,20)=   -2.475516_JWRB
      CIDEAC( 9,20)=   -2.862222_JWRB
      CIDEAC(10,20)=   -3.349712_JWRB
      CIDEAC(11,20)=   -3.842627_JWRB
      CIDEAC(12,20)=   -4.558055_JWRB
      CIDEAC(13,20)=   -5.257342_JWRB
      CIDEAC(14,20)=   -6.072487_JWRB
      CIDEAC(15,20)=   -7.056084_JWRB
      CIDEAC(16,20)=   -7.865724_JWRB
 
! h =     2.200000    
      CIDEAC( 6,21)=   -1.758808_JWRB
      CIDEAC( 7,21)=   -2.097306_JWRB
      CIDEAC( 8,21)=   -2.397066_JWRB
      CIDEAC( 9,21)=   -2.773985_JWRB
      CIDEAC(10,21)=   -3.249429_JWRB
      CIDEAC(11,21)=   -3.731577_JWRB
      CIDEAC(12,21)=   -4.400121_JWRB
      CIDEAC(13,21)=   -5.056503_JWRB
      CIDEAC(14,21)=   -5.947149_JWRB
      CIDEAC(15,21)=   -6.758241_JWRB
      CIDEAC(16,21)=   -7.665355_JWRB
 
! h =     2.300000    
      CIDEAC( 6,22)=   -1.742336_JWRB
      CIDEAC( 7,22)=   -2.030318_JWRB
      CIDEAC( 8,22)=   -2.380843_JWRB
      CIDEAC( 9,22)=   -2.738259_JWRB
      CIDEAC(10,22)=   -3.130276_JWRB
      CIDEAC(11,22)=   -3.666918_JWRB
      CIDEAC(12,22)=   -4.254450_JWRB
      CIDEAC(13,22)=   -4.941644_JWRB
      CIDEAC(14,22)=   -5.660053_JWRB
      CIDEAC(15,22)=   -6.561189_JWRB
      CIDEAC(16,22)=   -7.462294_JWRB
 
! h =     2.400000    
      CIDEAC( 6,23)=   -1.682265_JWRB
      CIDEAC( 7,23)=   -1.974037_JWRB
      CIDEAC( 8,23)=   -2.300742_JWRB
      CIDEAC( 9,23)=   -2.670091_JWRB
      CIDEAC(10,23)=   -3.057821_JWRB
      CIDEAC(11,23)=   -3.594969_JWRB
      CIDEAC(12,23)=   -4.102533_JWRB
      CIDEAC(13,23)=   -4.751281_JWRB
      CIDEAC(14,23)=   -5.549400_JWRB
      CIDEAC(15,23)=   -6.337473_JWRB
      CIDEAC(16,23)=   -7.156771_JWRB
 
! h =     2.500000    
      CIDEAC( 6,24)=   -1.648687_JWRB
      CIDEAC( 7,24)=   -1.903361_JWRB
      CIDEAC( 8,24)=   -2.280064_JWRB
      CIDEAC( 9,24)=   -2.629137_JWRB
      CIDEAC(10,24)=   -3.024701_JWRB
      CIDEAC(11,24)=   -3.445110_JWRB
      CIDEAC(12,24)=   -3.987420_JWRB
      CIDEAC(13,24)=   -4.637731_JWRB
      CIDEAC(14,24)=   -5.365729_JWRB
      CIDEAC(15,24)=   -6.134120_JWRB
      CIDEAC(16,24)=   -6.950096_JWRB

! h =     2.600000    
      CIDEAC( 6,25)=   -1.641218_JWRB
      CIDEAC( 7,25)=   -1.956000_JWRB
      CIDEAC( 8,25)=   -2.195282_JWRB
      CIDEAC( 9,25)=   -2.547474_JWRB
      CIDEAC(10,25)=   -2.952879_JWRB
      CIDEAC(11,25)=   -3.382897_JWRB
      CIDEAC(12,25)=   -3.878702_JWRB
      CIDEAC(13,25)=   -4.551814_JWRB
      CIDEAC(14,25)=   -5.228397_JWRB
      CIDEAC(15,25)=   -5.965115_JWRB
      CIDEAC(16,25)=   -6.785750_JWRB
 
! h =     2.700000    
      CIDEAC( 6,26)=   -1.600150_JWRB
      CIDEAC( 7,26)=   -1.844258_JWRB
      CIDEAC( 8,26)=   -2.257075_JWRB
      CIDEAC( 9,26)=   -2.519011_JWRB
      CIDEAC(10,26)=   -2.883503_JWRB
      CIDEAC(11,26)=   -3.288872_JWRB
      CIDEAC(12,26)=   -3.823005_JWRB
      CIDEAC(13,26)=   -4.376071_JWRB
      CIDEAC(14,26)=   -5.016023_JWRB
      CIDEAC(15,26)=   -5.787023_JWRB
      CIDEAC(16,26)=   -6.621093_JWRB
 
! h =     2.800000    
      CIDEAC( 6,27)=   -1.586459_JWRB
      CIDEAC( 7,27)=   -1.823506_JWRB
      CIDEAC( 8,27)=   -2.191002_JWRB
      CIDEAC( 9,27)=   -2.443355_JWRB
      CIDEAC(10,27)=   -2.839071_JWRB
      CIDEAC(11,27)=   -3.269969_JWRB
      CIDEAC(12,27)=   -3.717166_JWRB
      CIDEAC(13,27)=   -4.261279_JWRB
      CIDEAC(14,27)=   -4.884097_JWRB
      CIDEAC(15,27)=   -5.582966_JWRB
      CIDEAC(16,27)=   -6.413867_JWRB
 
! h =     2.900000    
      CIDEAC( 6,28)=   -1.530613_JWRB
      CIDEAC( 7,28)=   -1.881014_JWRB
      CIDEAC( 8,28)=   -2.084314_JWRB
      CIDEAC( 9,28)=   -2.427499_JWRB
      CIDEAC(10,28)=   -2.796773_JWRB
      CIDEAC(11,28)=   -3.155616_JWRB
      CIDEAC(12,28)=   -3.628541_JWRB
      CIDEAC(13,28)=   -4.160691_JWRB
      CIDEAC(14,28)=   -4.748525_JWRB
      CIDEAC(15,28)=   -5.441663_JWRB
      CIDEAC(16,28)=   -6.239977_JWRB
 
! h =     3.000000    
      CIDEAC( 6,29)=   -1.520659_JWRB
      CIDEAC( 7,29)=   -1.827379_JWRB
      CIDEAC( 8,29)=   -2.099156_JWRB
      CIDEAC( 9,29)=   -2.403214_JWRB
      CIDEAC(10,29)=   -2.764410_JWRB
      CIDEAC(11,29)=   -3.155527_JWRB
      CIDEAC(12,29)=   -3.555711_JWRB
      CIDEAC(13,29)=   -4.072697_JWRB
      CIDEAC(14,29)=   -4.656217_JWRB
      CIDEAC(15,29)=   -5.335746_JWRB
      CIDEAC(16,29)=   -6.086100_JWRB
 
! h =     3.100000    
      CIDEAC( 6,30)=   -1.468517_JWRB
      CIDEAC( 7,30)=   -1.802037_JWRB
      CIDEAC( 8,30)=   -2.039742_JWRB
      CIDEAC( 9,30)=   -2.369589_JWRB
      CIDEAC(10,30)=   -2.647677_JWRB
      CIDEAC(11,30)=   -3.087644_JWRB
      CIDEAC(12,30)=   -3.489749_JWRB
      CIDEAC(13,30)=   -3.964992_JWRB
      CIDEAC(14,30)=   -4.532036_JWRB
      CIDEAC(15,30)=   -5.117778_JWRB
      CIDEAC(16,30)=   -5.921745_JWRB
 
! h =     3.200000    
      CIDEAC( 6,31)=   -1.462488_JWRB
      CIDEAC( 7,31)=   -1.748595_JWRB
      CIDEAC( 8,31)=   -1.976536_JWRB
      CIDEAC( 9,31)=   -2.330469_JWRB
      CIDEAC(10,31)=   -2.633314_JWRB
      CIDEAC(11,31)=   -2.967257_JWRB
      CIDEAC(12,31)=   -3.411971_JWRB
      CIDEAC(13,31)=   -3.881270_JWRB
      CIDEAC(14,31)=   -4.440731_JWRB
      CIDEAC(15,31)=   -5.080121_JWRB
      CIDEAC(16,31)=   -5.782035_JWRB
 
! h =     3.300000    
      CIDEAC( 6,32)=   -1.404058_JWRB
      CIDEAC( 7,32)=   -1.682052_JWRB
      CIDEAC( 8,32)=   -1.989879_JWRB
      CIDEAC( 9,32)=   -2.303284_JWRB
      CIDEAC(10,32)=   -2.609248_JWRB
      CIDEAC(11,32)=   -2.970038_JWRB
      CIDEAC(12,32)=   -3.335601_JWRB
      CIDEAC(13,32)=   -3.827919_JWRB
      CIDEAC(14,32)=   -4.264277_JWRB
      CIDEAC(15,32)=   -4.969507_JWRB
      CIDEAC(16,32)=   -5.624526_JWRB
 
! h =     3.400000    
      CIDEAC( 6,33)=   -1.421658_JWRB
      CIDEAC( 7,33)=   -1.694868_JWRB
      CIDEAC( 8,33)=   -1.939824_JWRB
      CIDEAC( 9,33)=   -2.256738_JWRB
      CIDEAC(10,33)=   -2.563078_JWRB
      CIDEAC(11,33)=   -2.854945_JWRB
      CIDEAC(12,33)=   -3.300406_JWRB
      CIDEAC(13,33)=   -3.748025_JWRB
      CIDEAC(14,33)=   -4.209355_JWRB
      CIDEAC(15,33)=   -4.798466_JWRB
      CIDEAC(16,33)=   -5.484742_JWRB
 
! h =     3.500000    
      CIDEAC( 6,34)=   -1.433080_JWRB
      CIDEAC( 7,34)=   -1.658592_JWRB
      CIDEAC( 8,34)=   -1.953619_JWRB
      CIDEAC( 9,34)=   -2.189728_JWRB
      CIDEAC(10,34)=   -2.483670_JWRB
      CIDEAC(11,34)=   -2.876348_JWRB
      CIDEAC(12,34)=   -3.227396_JWRB
      CIDEAC(13,34)=   -3.643594_JWRB
      CIDEAC(14,34)=   -4.174257_JWRB
      CIDEAC(15,34)=   -4.713649_JWRB
      CIDEAC(16,34)=   -5.365061_JWRB

! h =     3.600000    
      CIDEAC( 6,35)=   -1.404438_JWRB
      CIDEAC( 7,35)=   -1.603982_JWRB
      CIDEAC( 8,35)=   -1.976813_JWRB
      CIDEAC( 9,35)=   -2.186624_JWRB
      CIDEAC(10,35)=   -2.514024_JWRB
      CIDEAC(11,35)=   -2.848864_JWRB
      CIDEAC(12,35)=   -3.172881_JWRB
      CIDEAC(13,35)=   -3.602279_JWRB
      CIDEAC(14,35)=   -4.101719_JWRB
      CIDEAC(15,35)=   -4.649209_JWRB
      CIDEAC(16,35)=   -5.305425_JWRB
 
! h =     3.700000    
      CIDEAC( 6,36)=   -1.392771_JWRB
      CIDEAC( 7,36)=   -1.646687_JWRB
      CIDEAC( 8,36)=   -1.921211_JWRB
      CIDEAC( 9,36)=   -2.149506_JWRB
      CIDEAC(10,36)=   -2.463910_JWRB
      CIDEAC(11,36)=   -2.780142_JWRB
      CIDEAC(12,36)=   -3.117836_JWRB
      CIDEAC(13,36)=   -3.524015_JWRB
      CIDEAC(14,36)=   -4.022181_JWRB
      CIDEAC(15,36)=   -4.536423_JWRB
      CIDEAC(16,36)=   -5.194000_JWRB
 
!     extrapolation:
      DO IH=1,NICH
        DCI=CIDEAC(6,IH)-CIDEAC(1,IH)
        DO IT=2,5
          CIDEAC(IT,IH)=CIDEAC(1,IH)+DCI*(IT-1)*DTIC/(5*DTIC)
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('CIGETDEAC',1,ZHOOK_HANDLE)

      END SUBROUTINE CIGETDEAC
