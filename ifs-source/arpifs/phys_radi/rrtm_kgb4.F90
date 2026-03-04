SUBROUTINE RRTM_KGB4

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     R. Elkhatib 12-10-2005 Split for faster and more robust compilation.
!     G.Mozdzynski March 2011 read constants from files
!     ABozzo 201306 updated to rrtmg v4.85
!     T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN    ,ONLY : NULRAD
USE MPL_MODULE,ONLY : MPL_BROADCAST
USE YOMTAG    ,ONLY : MTAGRAD
USE YOMMP0    , ONLY : NPROC, MYPROC

USE YOERRTO4 , ONLY :  KAO     ,KBO     ,SELFREFO   ,FORREFO, FRACREFAO  ,FRACREFBO, &
      &  KAO_D, KBO_D


!     ------------------------------------------------------------------

IMPLICIT NONE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('RRTM_KGB4',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  READ(NULRAD,ERR=1001)  KAO_D,KBO_D
  KAO = REAL(KAO_D,JPRB)
  KBO = REAL(KBO_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KAO,MTAGRAD,1,CDSTRING='RRTM_KGB4:')
  CALL MPL_BROADCAST (KBO,MTAGRAD,1,CDSTRING='RRTM_KGB4:')
ENDIF


! Planck fraction mapping level : P = 142.5940 mbar, T = 215.70 K
      FRACREFAO(:, 1) = (/ &
     &   1.5572E-01_JPRB,1.4925E-01_JPRB,1.4107E-01_JPRB,1.3126E-01_JPRB,1.1791E-01_JPRB,1.0173E-01_JPRB, &
     &   8.2949E-02_JPRB,6.2393E-02_JPRB,4.2146E-02_JPRB,4.5907E-03_JPRB,3.7965E-03_JPRB,2.9744E-03_JPRB, &
     &   2.2074E-03_JPRB,1.4063E-03_JPRB,5.3012E-04_JPRB,7.4595E-05_JPRB/)
      FRACREFAO(:, 2) = (/ &
     &   1.5572E-01_JPRB,1.4925E-01_JPRB,1.4107E-01_JPRB,1.3126E-01_JPRB,1.1791E-01_JPRB,1.0173E-01_JPRB, &
     &   8.2949E-02_JPRB,6.2392E-02_JPRB,4.2146E-02_JPRB,4.5906E-03_JPRB,3.7965E-03_JPRB,2.9745E-03_JPRB, &
     &   2.2074E-03_JPRB,1.4063E-03_JPRB,5.3012E-04_JPRB,7.4595E-05_JPRB/)
      FRACREFAO(:, 3) = (/ &
     &   1.5572E-01_JPRB,1.4925E-01_JPRB,1.4107E-01_JPRB,1.3126E-01_JPRB,1.1791E-01_JPRB,1.0173E-01_JPRB, &
     &   8.2949E-02_JPRB,6.2393E-02_JPRB,4.2146E-02_JPRB,4.5907E-03_JPRB,3.7965E-03_JPRB,2.9745E-03_JPRB, &
     &   2.2074E-03_JPRB,1.4063E-03_JPRB,5.3012E-04_JPRB,7.4595E-05_JPRB/)
      FRACREFAO(:, 4) = (/ &
     &   1.5572E-01_JPRB,1.4925E-01_JPRB,1.4107E-01_JPRB,1.3126E-01_JPRB,1.1791E-01_JPRB,1.0173E-01_JPRB, &
     &   8.2949E-02_JPRB,6.2393E-02_JPRB,4.2146E-02_JPRB,4.5907E-03_JPRB,3.7964E-03_JPRB,2.9744E-03_JPRB, &
     &   2.2074E-03_JPRB,1.4063E-03_JPRB,5.3012E-04_JPRB,7.4595E-05_JPRB/)
      FRACREFAO(:, 5) = (/ &
     &   1.5572E-01_JPRB,1.4925E-01_JPRB,1.4107E-01_JPRB,1.3126E-01_JPRB,1.1791E-01_JPRB,1.0173E-01_JPRB, &
     &   8.2949E-02_JPRB,6.2393E-02_JPRB,4.2146E-02_JPRB,4.5907E-03_JPRB,3.7965E-03_JPRB,2.9744E-03_JPRB, &
     &   2.2074E-03_JPRB,1.4063E-03_JPRB,5.3012E-04_JPRB,7.4595E-05_JPRB/)
      FRACREFAO(:, 6) = (/ &
     &   1.5572E-01_JPRB,1.4925E-01_JPRB,1.4107E-01_JPRB,1.3126E-01_JPRB,1.1791E-01_JPRB,1.0173E-01_JPRB, &
     &   8.2949E-02_JPRB,6.2393E-02_JPRB,4.2146E-02_JPRB,4.5907E-03_JPRB,3.7965E-03_JPRB,2.9744E-03_JPRB, &
     &   2.2074E-03_JPRB,1.4063E-03_JPRB,5.3012E-04_JPRB,7.4595E-05_JPRB/)
      FRACREFAO(:, 7) = (/ &
     &   1.5572E-01_JPRB,1.4926E-01_JPRB,1.4107E-01_JPRB,1.3126E-01_JPRB,1.1791E-01_JPRB,1.0173E-01_JPRB, &
     &   8.2949E-02_JPRB,6.2393E-02_JPRB,4.2146E-02_JPRB,4.5908E-03_JPRB,3.7964E-03_JPRB,2.9745E-03_JPRB, &
     &   2.2074E-03_JPRB,1.4063E-03_JPRB,5.3012E-04_JPRB,7.4595E-05_JPRB/)
      FRACREFAO(:, 8) = (/ &
     &   1.5571E-01_JPRB,1.4926E-01_JPRB,1.4107E-01_JPRB,1.3125E-01_JPRB,1.1791E-01_JPRB,1.0173E-01_JPRB, &
     &   8.2949E-02_JPRB,6.2393E-02_JPRB,4.2146E-02_JPRB,4.5907E-03_JPRB,3.7964E-03_JPRB,2.9744E-03_JPRB, &
     &   2.2074E-03_JPRB,1.4063E-03_JPRB,5.3012E-04_JPRB,7.4595E-05_JPRB/)
      FRACREFAO(:, 9) = (/ &
     &   1.5952E-01_JPRB,1.5155E-01_JPRB,1.4217E-01_JPRB,1.3077E-01_JPRB,1.1667E-01_JPRB,1.0048E-01_JPRB, &
     &   8.1511E-02_JPRB,6.1076E-02_JPRB,4.1111E-02_JPRB,4.4432E-03_JPRB,3.6910E-03_JPRB,2.9076E-03_JPRB, &
     &   2.1329E-03_JPRB,1.3566E-03_JPRB,5.2235E-04_JPRB,7.9935E-05_JPRB/)

! Planck fraction mapping level : P = 95.58350 mb, T = 215.70 K
      FRACREFBO(:, 1) = (/ &
     &   1.5558E-01_JPRB,1.4931E-01_JPRB,1.4104E-01_JPRB,1.3124E-01_JPRB,1.1793E-01_JPRB,1.0160E-01_JPRB, &
     &   8.3142E-02_JPRB,6.2403E-02_JPRB,4.2170E-02_JPRB,4.5935E-03_JPRB,3.7976E-03_JPRB,2.9986E-03_JPRB, &
     &   2.1890E-03_JPRB,1.4061E-03_JPRB,5.3005E-04_JPRB,7.4587E-05_JPRB/)
      FRACREFBO(:, 2) = (/ &
     &   1.5558E-01_JPRB,1.4932E-01_JPRB,1.4104E-01_JPRB,1.3124E-01_JPRB,1.1792E-01_JPRB,1.0159E-01_JPRB, &
     &   8.3142E-02_JPRB,6.2403E-02_JPRB,4.2170E-02_JPRB,4.5935E-03_JPRB,3.7976E-03_JPRB,2.9986E-03_JPRB, &
     &   2.1890E-03_JPRB,1.4061E-03_JPRB,5.3005E-04_JPRB,7.4587E-05_JPRB/)
      FRACREFBO(:, 3) = (/ &
     &   1.5558E-01_JPRB,1.4933E-01_JPRB,1.4103E-01_JPRB,1.3124E-01_JPRB,1.1792E-01_JPRB,1.0159E-01_JPRB, &
     &   8.3142E-02_JPRB,6.2403E-02_JPRB,4.2170E-02_JPRB,4.5935E-03_JPRB,3.7976E-03_JPRB,2.9986E-03_JPRB, &
     &   2.1890E-03_JPRB,1.4061E-03_JPRB,5.3005E-04_JPRB,7.4587E-05_JPRB/)
      FRACREFBO(:, 4) = (/ &
     &   1.5569E-01_JPRB,1.4926E-01_JPRB,1.4102E-01_JPRB,1.3122E-01_JPRB,1.1791E-01_JPRB,1.0159E-01_JPRB, &
     &   8.3141E-02_JPRB,6.2403E-02_JPRB,4.2170E-02_JPRB,4.5935E-03_JPRB,3.7976E-03_JPRB,2.9986E-03_JPRB, &
     &   2.1890E-03_JPRB,1.4061E-03_JPRB,5.3005E-04_JPRB,7.4587E-05_JPRB/)
      FRACREFBO(:, 5) = (/ &
     &   1.5947E-01_JPRB,1.5132E-01_JPRB,1.4195E-01_JPRB,1.3061E-01_JPRB,1.1680E-01_JPRB,1.0054E-01_JPRB, &
     &   8.1785E-02_JPRB,6.1212E-02_JPRB,4.1276E-02_JPRB,4.4424E-03_JPRB,3.6628E-03_JPRB,2.8943E-03_JPRB, &
     &   2.1134E-03_JPRB,1.3457E-03_JPRB,5.1024E-04_JPRB,7.3998E-05_JPRB/)

!     ------------------------------------------------------------------

!     The array KAO contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels > ~100mb, temperatures, and ratios
!     of water vapor to CO2.  The first index in the array, JS, runs
!     from 1 to 9 and corresponds to different water vapor to CO2 ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = h2o/(h20 + (rat) * co2), where rat is the ratio of the integrated
!     line strength in the band of co2 to that of h2o.  For instance,
!     JS=1 refers to dry air (eta = 0), JS = 9 corresponds to eta = 1.0.
!     The 2nd index in the array, JT, which runs from 1 to 5, corresponds 
!     to different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this pressure 
!     level, JT = 2 refers to the temperature TREF-15,
!     JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the reference pressure level (e.g. JP = 1 is for a
!     pressure of 1053.63 mb).  The fourth index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.



!     The array KBO contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels  < ~100mb, temperatures, and ratios
!     of O3 to CO2.  The first index in the array, JS, runs from 1 to 6, 
!     and corresponds to different O3 to CO2 ratios, as expressed through 
!     the binary species parameter eta, defined as eta = O3/(O3+RAT*H2O), 
!     where RAT is the ratio of the integrated line strength in the band 
!     of CO2 to that of O3.  For instance, JS=1 refers to no O3 (eta = 0) 
!     and JS = 5 corresponds to eta = 1.0.  The second index, JT, which
!     runs from 1 to 5, corresponds to different temperatures.  More 
!     specifically, JT = 3 means that the data are for the corresponding 
!     reference temperature TREF for this  pressure level, JT = 2 refers 
!     to the TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and
!     JT = 5 is for TREF+30.  The third index, JP, runs from 13 to 59 and
!     refers to the corresponding pressure level in PREF (e.g. JP = 13 is
!     for a pressure of 95.5835 mb).  The fourth index, IG, goes from 1 to
!     16, and tells us which g-interval the absorption coefficients are for.



!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      FORREFO(1,:) = (/ &
     &3.3839E-04_JPRB,2.4739E-04_JPRB,2.2846E-04_JPRB,2.3376E-04_JPRB,2.2622E-04_JPRB,2.3188E-04_JPRB, &
     &2.2990E-04_JPRB,2.2532E-04_JPRB,2.1233E-04_JPRB,2.0593E-04_JPRB,2.0716E-04_JPRB,2.0809E-04_JPRB, &
     &2.0889E-04_JPRB,2.0932E-04_JPRB,2.0944E-04_JPRB,2.0945E-04_JPRB/)
      FORREFO(2,:) = (/ &
     &3.4391E-04_JPRB,2.6022E-04_JPRB,2.3449E-04_JPRB,2.4544E-04_JPRB,2.3831E-04_JPRB,2.3014E-04_JPRB, &
     &2.3729E-04_JPRB,2.2726E-04_JPRB,2.1892E-04_JPRB,1.9223E-04_JPRB,2.1291E-04_JPRB,2.1406E-04_JPRB, &
     &2.1491E-04_JPRB,2.1548E-04_JPRB,2.1562E-04_JPRB,2.1567E-04_JPRB/)
      FORREFO(3,:) = (/ &
     &3.4219E-04_JPRB,2.7334E-04_JPRB,2.3727E-04_JPRB,2.4515E-04_JPRB,2.5272E-04_JPRB,2.4212E-04_JPRB, &
     &2.3824E-04_JPRB,2.3615E-04_JPRB,2.2724E-04_JPRB,2.2381E-04_JPRB,1.9634E-04_JPRB,2.1625E-04_JPRB, &
     &2.1963E-04_JPRB,2.2032E-04_JPRB,2.2057E-04_JPRB,2.2058E-04_JPRB/)
      FORREFO(4,:) = (/ &
     &3.1684E-04_JPRB,2.4823E-04_JPRB,2.4890E-04_JPRB,2.4577E-04_JPRB,2.4106E-04_JPRB,2.4353E-04_JPRB, &
     &2.4038E-04_JPRB,2.3932E-04_JPRB,2.3604E-04_JPRB,2.3773E-04_JPRB,2.4243E-04_JPRB,2.2597E-04_JPRB, &
     &2.2879E-04_JPRB,2.2440E-04_JPRB,2.1104E-04_JPRB,2.1460E-04_JPRB/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      SELFREFO(:, 1) = (/ &
     & 2.62922E-01_JPRB, 2.29106E-01_JPRB, 1.99640E-01_JPRB, 1.73964E-01_JPRB, 1.51589E-01_JPRB, &
     & 1.32093E-01_JPRB, 1.15104E-01_JPRB, 1.00300E-01_JPRB, 8.74000E-02_JPRB, 7.61592E-02_JPRB/)
      SELFREFO(:, 2) = (/ &
     & 2.45448E-01_JPRB, 2.13212E-01_JPRB, 1.85210E-01_JPRB, 1.60886E-01_JPRB, 1.39756E-01_JPRB, &
     & 1.21401E-01_JPRB, 1.05457E-01_JPRB, 9.16070E-02_JPRB, 7.95759E-02_JPRB, 6.91249E-02_JPRB/)
      SELFREFO(:, 3) = (/ &
     & 2.41595E-01_JPRB, 2.09697E-01_JPRB, 1.82010E-01_JPRB, 1.57979E-01_JPRB, 1.37121E-01_JPRB, &
     & 1.19016E-01_JPRB, 1.03302E-01_JPRB, 8.96630E-02_JPRB, 7.78246E-02_JPRB, 6.75492E-02_JPRB/)
      SELFREFO(:, 4) = (/ &
     & 2.44818E-01_JPRB, 2.12172E-01_JPRB, 1.83880E-01_JPRB, 1.59360E-01_JPRB, 1.38110E-01_JPRB, &
     & 1.19694E-01_JPRB, 1.03733E-01_JPRB, 8.99010E-02_JPRB, 7.79131E-02_JPRB, 6.75238E-02_JPRB/)
      SELFREFO(:, 5) = (/ &
     & 2.43458E-01_JPRB, 2.10983E-01_JPRB, 1.82840E-01_JPRB, 1.58451E-01_JPRB, 1.37315E-01_JPRB, &
     & 1.18998E-01_JPRB, 1.03125E-01_JPRB, 8.93690E-02_JPRB, 7.74480E-02_JPRB, 6.71171E-02_JPRB/)
      SELFREFO(:, 6) = (/ &
     & 2.40186E-01_JPRB, 2.08745E-01_JPRB, 1.81420E-01_JPRB, 1.57672E-01_JPRB, 1.37032E-01_JPRB, &
     & 1.19095E-01_JPRB, 1.03505E-01_JPRB, 8.99560E-02_JPRB, 7.81806E-02_JPRB, 6.79467E-02_JPRB/)
      SELFREFO(:, 7) = (/ &
     & 2.42752E-01_JPRB, 2.10579E-01_JPRB, 1.82670E-01_JPRB, 1.58460E-01_JPRB, 1.37459E-01_JPRB, &
     & 1.19240E-01_JPRB, 1.03437E-01_JPRB, 8.97280E-02_JPRB, 7.78359E-02_JPRB, 6.75200E-02_JPRB/)
      SELFREFO(:, 8) = (/ &
     & 2.39620E-01_JPRB, 2.08166E-01_JPRB, 1.80840E-01_JPRB, 1.57101E-01_JPRB, 1.36479E-01_JPRB, &
     & 1.18563E-01_JPRB, 1.03000E-01_JPRB, 8.94790E-02_JPRB, 7.77332E-02_JPRB, 6.75292E-02_JPRB/)
      SELFREFO(:, 9) = (/ &
     & 2.38856E-01_JPRB, 2.07166E-01_JPRB, 1.79680E-01_JPRB, 1.55841E-01_JPRB, 1.35165E-01_JPRB, &
     & 1.17232E-01_JPRB, 1.01678E-01_JPRB, 8.81880E-02_JPRB, 7.64877E-02_JPRB, 6.63397E-02_JPRB/)
      SELFREFO(:,10) = (/ &
     & 2.29821E-01_JPRB, 2.00586E-01_JPRB, 1.75070E-01_JPRB, 1.52800E-01_JPRB, 1.33363E-01_JPRB, &
     & 1.16398E-01_JPRB, 1.01591E-01_JPRB, 8.86680E-02_JPRB, 7.73887E-02_JPRB, 6.75443E-02_JPRB/)
      SELFREFO(:,11) = (/ &
     & 2.39945E-01_JPRB, 2.08186E-01_JPRB, 1.80630E-01_JPRB, 1.56722E-01_JPRB, 1.35978E-01_JPRB, &
     & 1.17980E-01_JPRB, 1.02364E-01_JPRB, 8.88150E-02_JPRB, 7.70594E-02_JPRB, 6.68598E-02_JPRB/)
      SELFREFO(:,12) = (/ &
     & 2.40271E-01_JPRB, 2.08465E-01_JPRB, 1.80870E-01_JPRB, 1.56927E-01_JPRB, 1.36154E-01_JPRB, &
     & 1.18131E-01_JPRB, 1.02494E-01_JPRB, 8.89260E-02_JPRB, 7.71545E-02_JPRB, 6.69412E-02_JPRB/)
      SELFREFO(:,13) = (/ &
     & 2.40503E-01_JPRB, 2.08670E-01_JPRB, 1.81050E-01_JPRB, 1.57086E-01_JPRB, 1.36294E-01_JPRB, &
     & 1.18254E-01_JPRB, 1.02602E-01_JPRB, 8.90210E-02_JPRB, 7.72380E-02_JPRB, 6.70147E-02_JPRB/)
      SELFREFO(:,14) = (/ &
     & 2.40670E-01_JPRB, 2.08811E-01_JPRB, 1.81170E-01_JPRB, 1.57188E-01_JPRB, 1.36380E-01_JPRB, &
     & 1.18327E-01_JPRB, 1.02663E-01_JPRB, 8.90730E-02_JPRB, 7.72819E-02_JPRB, 6.70517E-02_JPRB/)
      SELFREFO(:,15) = (/ &
     & 2.40711E-01_JPRB, 2.08846E-01_JPRB, 1.81200E-01_JPRB, 1.57213E-01_JPRB, 1.36402E-01_JPRB, &
     & 1.18346E-01_JPRB, 1.02679E-01_JPRB, 8.90870E-02_JPRB, 7.72939E-02_JPRB, 6.70621E-02_JPRB/)
      SELFREFO(:,16) = (/ &
     & 2.40727E-01_JPRB, 2.08859E-01_JPRB, 1.81210E-01_JPRB, 1.57221E-01_JPRB, 1.36408E-01_JPRB, &
     & 1.18350E-01_JPRB, 1.02682E-01_JPRB, 8.90890E-02_JPRB, 7.72952E-02_JPRB, 6.70627E-02_JPRB/)

IF (LHOOK) CALL DR_HOOK('RRTM_KGB4',1,ZHOOK_HANDLE)
RETURN

1001 CONTINUE
CALL ABOR1("RRTM_KGB4:ERROR READING FILE RADRRTM")

END SUBROUTINE RRTM_KGB4
