SUBROUTINE RRTM_KGB14

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)
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
USE YOMMP0    , ONLY : NPROC, MYPROC
USE MPL_MODULE,ONLY : MPL_BROADCAST
USE YOMTAG    ,ONLY : MTAGRAD

USE YOERRTO14, ONLY : KAO     ,KBO     ,SELFREFO, FORREFO ,FRACREFAO  ,FRACREFBO, &
  &  KAO_D, KBO_D

!     ------------------------------------------------------------------

IMPLICIT NONE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('RRTM_KGB14',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  READ(NULRAD,ERR=1001) KAO_D,KBO_D
  KAO = REAL(KAO_D,JPRB)
  KBO = REAL(KBO_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KAO,MTAGRAD,1,CDSTRING='RRTM_KGB14:')
  CALL MPL_BROADCAST (KBO,MTAGRAD,1,CDSTRING='RRTM_KGB14:')
ENDIF

! Planck fraction mapping level : P = 142.5940 mb, T = 215.70 K
      FRACREFAO(:) = (/ &
     &  1.9360E-01_JPRB, 1.7276E-01_JPRB, 1.4811E-01_JPRB, 1.2238E-01_JPRB, &
     &  1.0242E-01_JPRB, 8.6830E-02_JPRB, 7.1890E-02_JPRB, 5.4030E-02_JPRB, &
     &  3.5075E-02_JPRB, 3.8052E-03_JPRB, 3.1458E-03_JPRB, 2.4873E-03_JPRB, &
     &  1.8182E-03_JPRB, 1.1563E-03_JPRB, 4.3251E-04_JPRB, 5.7744E-05_JPRB/)

! Planck fraction mapping level : P = 4.758820mb, T = 250.85 K
      FRACREFBO(:) = (/ &
     &  1.8599E-01_JPRB, 1.6646E-01_JPRB, 1.4264E-01_JPRB, 1.2231E-01_JPRB, &
     &  1.0603E-01_JPRB, 9.2014E-02_JPRB, 7.5287E-02_JPRB, 5.6758E-02_JPRB, &
     &  3.8386E-02_JPRB, 4.2139E-03_JPRB, 3.5399E-03_JPRB, 2.7381E-03_JPRB, &
     &  1.9202E-03_JPRB, 1.2083E-03_JPRB, 4.5395E-04_JPRB, 6.2699E-05_JPRB/)

!     ------------------------------------------------------------------

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels > ~100mb and temperatures.  The first
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the corresponding TREF for this  pressure level, 
!     JT = 2 refers to the temperatureTREF-15, JT = 1 is for TREF-30, 
!     JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  The second 
!     index, JP, runs from 1 to 13 and refers to the corresponding 
!     pressure level in PREF (e.g. JP = 1 is for a pressure of 1053.63 mb).  
!     The third index, IG, goes from 1 to 16, and tells us which 
!     g-interval the absorption coefficients are for.



!     The array KBO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels < ~100mb and temperatures. The first 
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the reference temperature TREF for this pressure 
!     level, JT = 2 refers to the temperature TREF-15, JT = 1 is for
!     TREF-30, JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  
!     The second index, JP, runs from 13 to 59 and refers to the JPth
!     reference pressure level (see taumol.f for the value of these
!     pressure levels in mb).  The third index, IG, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.


!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      FORREFO(1,:) = (/ &
     &2.7075E-06_JPRB,2.2609E-06_JPRB,1.5633E-06_JPRB,8.7484E-07_JPRB,5.5470E-07_JPRB,4.8456E-07_JPRB, &
     &4.7463E-07_JPRB,4.6154E-07_JPRB,4.4425E-07_JPRB,4.2960E-07_JPRB,4.2626E-07_JPRB,4.1715E-07_JPRB, &
     &4.2607E-07_JPRB,3.6616E-07_JPRB,2.6366E-07_JPRB,2.6029E-07_JPRB/)
      FORREFO(2,:) = (/ &
     &2.6759E-06_JPRB,2.2237E-06_JPRB,1.4466E-06_JPRB,9.3032E-07_JPRB,6.4927E-07_JPRB,5.4809E-07_JPRB, &
     &4.9504E-07_JPRB,4.6305E-07_JPRB,4.4873E-07_JPRB,4.2146E-07_JPRB,4.2176E-07_JPRB,4.2812E-07_JPRB, &
     &4.0529E-07_JPRB,4.0969E-07_JPRB,2.9442E-07_JPRB,2.6821E-07_JPRB/)
      FORREFO(3,:) = (/ &
     &2.6608E-06_JPRB,2.1140E-06_JPRB,1.4838E-06_JPRB,9.2083E-07_JPRB,6.3350E-07_JPRB,5.7195E-07_JPRB, &
     &6.2253E-07_JPRB,5.1783E-07_JPRB,4.4749E-07_JPRB,4.3261E-07_JPRB,4.2553E-07_JPRB,4.2175E-07_JPRB, &
     &4.1085E-07_JPRB,4.0358E-07_JPRB,3.5340E-07_JPRB,2.7191E-07_JPRB/)
      FORREFO(4,:) = (/ &
     &2.6412E-06_JPRB,1.9814E-06_JPRB,1.2672E-06_JPRB,8.1129E-07_JPRB,7.1447E-07_JPRB,7.5026E-07_JPRB, &
     &7.4386E-07_JPRB,7.2759E-07_JPRB,7.3583E-07_JPRB,7.6493E-07_JPRB,8.8959E-07_JPRB,7.5534E-07_JPRB, &
     &5.3734E-07_JPRB,4.5572E-07_JPRB,4.1676E-07_JPRB,3.6198E-07_JPRB/)


!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      SELFREFO(:, 1) = (/ &
     & 4.67262E-03_JPRB, 3.95211E-03_JPRB, 3.34270E-03_JPRB, 2.82726E-03_JPRB, 2.39130E-03_JPRB, &
     & 2.02256E-03_JPRB, 1.71069E-03_JPRB, 1.44690E-03_JPRB, 1.22379E-03_JPRB, 1.03508E-03_JPRB/)
      SELFREFO(:, 2) = (/ &
     & 4.42593E-03_JPRB, 3.73338E-03_JPRB, 3.14920E-03_JPRB, 2.65643E-03_JPRB, 2.24076E-03_JPRB, &
     & 1.89014E-03_JPRB, 1.59438E-03_JPRB, 1.34490E-03_JPRB, 1.13446E-03_JPRB, 9.56943E-04_JPRB/)
      SELFREFO(:, 3) = (/ &
     & 3.96072E-03_JPRB, 3.33789E-03_JPRB, 2.81300E-03_JPRB, 2.37065E-03_JPRB, 1.99786E-03_JPRB, &
     & 1.68369E-03_JPRB, 1.41893E-03_JPRB, 1.19580E-03_JPRB, 1.00776E-03_JPRB, 8.49286E-04_JPRB/)
      SELFREFO(:, 4) = (/ &
     & 3.71833E-03_JPRB, 3.10030E-03_JPRB, 2.58500E-03_JPRB, 2.15535E-03_JPRB, 1.79711E-03_JPRB, &
     & 1.49841E-03_JPRB, 1.24936E-03_JPRB, 1.04170E-03_JPRB, 8.68558E-04_JPRB, 7.24195E-04_JPRB/)
      SELFREFO(:, 5) = (/ &
     & 3.55755E-03_JPRB, 2.95355E-03_JPRB, 2.45210E-03_JPRB, 2.03578E-03_JPRB, 1.69015E-03_JPRB, &
     & 1.40320E-03_JPRB, 1.16497E-03_JPRB, 9.67180E-04_JPRB, 8.02973E-04_JPRB, 6.66646E-04_JPRB/)
      SELFREFO(:, 6) = (/ &
     & 3.47601E-03_JPRB, 2.88628E-03_JPRB, 2.39660E-03_JPRB, 1.99000E-03_JPRB, 1.65238E-03_JPRB, &
     & 1.37204E-03_JPRB, 1.13927E-03_JPRB, 9.45980E-04_JPRB, 7.85487E-04_JPRB, 6.52224E-04_JPRB/)
      SELFREFO(:, 7) = (/ &
     & 3.44479E-03_JPRB, 2.86224E-03_JPRB, 2.37820E-03_JPRB, 1.97602E-03_JPRB, 1.64185E-03_JPRB, &
     & 1.36420E-03_JPRB, 1.13350E-03_JPRB, 9.41810E-04_JPRB, 7.82539E-04_JPRB, 6.50204E-04_JPRB/)
      SELFREFO(:, 8) = (/ &
     & 3.40154E-03_JPRB, 2.82953E-03_JPRB, 2.35370E-03_JPRB, 1.95789E-03_JPRB, 1.62864E-03_JPRB, &
     & 1.35476E-03_JPRB, 1.12694E-03_JPRB, 9.37430E-04_JPRB, 7.79788E-04_JPRB, 6.48655E-04_JPRB/)
      SELFREFO(:, 9) = (/ &
     & 3.39380E-03_JPRB, 2.82288E-03_JPRB, 2.34800E-03_JPRB, 1.95301E-03_JPRB, 1.62446E-03_JPRB, &
     & 1.35119E-03_JPRB, 1.12389E-03_JPRB, 9.34820E-04_JPRB, 7.77560E-04_JPRB, 6.46755E-04_JPRB/)
      SELFREFO(:,10) = (/ &
     & 3.37185E-03_JPRB, 2.80654E-03_JPRB, 2.33600E-03_JPRB, 1.94435E-03_JPRB, 1.61837E-03_JPRB, &
     & 1.34704E-03_JPRB, 1.12120E-03_JPRB, 9.33220E-04_JPRB, 7.76759E-04_JPRB, 6.46530E-04_JPRB/)
      SELFREFO(:,11) = (/ &
     & 3.37924E-03_JPRB, 2.81172E-03_JPRB, 2.33950E-03_JPRB, 1.94659E-03_JPRB, 1.61967E-03_JPRB, &
     & 1.34765E-03_JPRB, 1.12132E-03_JPRB, 9.33000E-04_JPRB, 7.76306E-04_JPRB, 6.45930E-04_JPRB/)
      SELFREFO(:,12) = (/ &
     & 3.39658E-03_JPRB, 2.82289E-03_JPRB, 2.34610E-03_JPRB, 1.94984E-03_JPRB, 1.62051E-03_JPRB, &
     & 1.34680E-03_JPRB, 1.11933E-03_JPRB, 9.30270E-04_JPRB, 7.73146E-04_JPRB, 6.42561E-04_JPRB/)
      SELFREFO(:,13) = (/ &
     & 3.36070E-03_JPRB, 2.79913E-03_JPRB, 2.33140E-03_JPRB, 1.94183E-03_JPRB, 1.61735E-03_JPRB, &
     & 1.34709E-03_JPRB, 1.12199E-03_JPRB, 9.34510E-04_JPRB, 7.78354E-04_JPRB, 6.48292E-04_JPRB/)
      SELFREFO(:,14) = (/ &
     & 3.40428E-03_JPRB, 2.81994E-03_JPRB, 2.33590E-03_JPRB, 1.93495E-03_JPRB, 1.60282E-03_JPRB, &
     & 1.32770E-03_JPRB, 1.09980E-03_JPRB, 9.11020E-04_JPRB, 7.54645E-04_JPRB, 6.25111E-04_JPRB/)
      SELFREFO(:,15) = (/ &
     & 3.27075E-03_JPRB, 2.70783E-03_JPRB, 2.24180E-03_JPRB, 1.85597E-03_JPRB, 1.53655E-03_JPRB, &
     & 1.27210E-03_JPRB, 1.05317E-03_JPRB, 8.71910E-04_JPRB, 7.21849E-04_JPRB, 5.97615E-04_JPRB/)
      SELFREFO(:,16) = (/ &
     & 3.23123E-03_JPRB, 2.67891E-03_JPRB, 2.22100E-03_JPRB, 1.84136E-03_JPRB, 1.52661E-03_JPRB, &
     & 1.26567E-03_JPRB, 1.04932E-03_JPRB, 8.69960E-04_JPRB, 7.21256E-04_JPRB, 5.97970E-04_JPRB/)

IF (LHOOK) CALL DR_HOOK('RRTM_KGB14',1,ZHOOK_HANDLE)
RETURN

1001 CONTINUE
CALL ABOR1("RRTM_KGB14:ERROR READING FILE RADRRTM")

END SUBROUTINE RRTM_KGB14
