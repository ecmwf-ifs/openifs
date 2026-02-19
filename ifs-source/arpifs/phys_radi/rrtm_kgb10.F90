SUBROUTINE RRTM_KGB10

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 10:  1390-1480 cm-1 (low - H2O; high - H2O)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     R. Elkhatib 12-10-2005 Split for faster and more robust compilation.
!     G.Mozdzynski March 2011 read constants from files
!     ABozzo 101306 updated to rrtmg v4.85
!     T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN    ,ONLY : NULRAD
USE YOMMP0    , ONLY : NPROC, MYPROC
USE MPL_MODULE,ONLY : MPL_BROADCAST
USE YOMTAG    ,ONLY : MTAGRAD

USE YOERRTO10, ONLY : KAO, KBO, KAO_D, KBO_D, FRACREFAO, FRACREFBO, SELFREFO, FORREFO

!     ------------------------------------------------------------------

IMPLICIT NONE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('RRTM_KGB10',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  READ(NULRAD,ERR=1001) KAO_D,KBO_D
  KAO = REAL(KAO_D,JPRB)
  KBO = REAL(KBO_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KAO,MTAGRAD,1,CDSTRING='RRTM_KGB10:')
  CALL MPL_BROADCAST (KBO,MTAGRAD,1,CDSTRING='RRTM_KGB10:')
ENDIF

! Planck fraction mapping level : P = 212.7250, T = 223.06 K
      FRACREFAO(:) = (/ &
     &  1.6909E-01_JPRB, 1.5419E-01_JPRB, 1.3999E-01_JPRB, 1.2637E-01_JPRB, &
     &  1.1429E-01_JPRB, 9.9676E-02_JPRB, 8.0093E-02_JPRB, 6.0283E-02_JPRB, &
     &  4.1077E-02_JPRB, 4.4857E-03_JPRB, 3.6545E-03_JPRB, 2.9243E-03_JPRB, &
     &  2.0407E-03_JPRB, 1.2891E-03_JPRB, 4.8767E-04_JPRB, 6.7748E-05_JPRB/)

! Planck fraction mapping level : P = 95.58350 mb, T = 215.70 K
      FRACREFBO(:) = (/ &
     &  1.7391E-01_JPRB, 1.5680E-01_JPRB, 1.4419E-01_JPRB, 1.2672E-01_JPRB, &
     &  1.0708E-01_JPRB, 9.7034E-02_JPRB, 7.8545E-02_JPRB, 5.9784E-02_JPRB, &
     &  4.0879E-02_JPRB, 4.4704E-03_JPRB, 3.7150E-03_JPRB, 2.9038E-03_JPRB, &
     &  2.1454E-03_JPRB, 1.2802E-03_JPRB, 4.8328E-04_JPRB, 6.7378E-05_JPRB/)


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
     &1.0515E-02_JPRB,1.4860E-02_JPRB,1.7181E-02_JPRB,1.6642E-02_JPRB,1.6644E-02_JPRB,1.5649E-02_JPRB, &
     &1.7734E-02_JPRB,1.7521E-02_JPRB,1.7868E-02_JPRB,1.8400E-02_JPRB,1.9361E-02_JPRB,2.1487E-02_JPRB, &
     &2.0192E-02_JPRB,1.6545E-02_JPRB,2.0922E-02_JPRB,2.0922E-02_JPRB/)
      FORREFO(2,:) = (/ &
     &1.0423E-02_JPRB,1.4593E-02_JPRB,1.6329E-02_JPRB,1.7071E-02_JPRB,1.7252E-02_JPRB,1.6188E-02_JPRB, &
     &1.7752E-02_JPRB,1.7913E-02_JPRB,1.7551E-02_JPRB,1.8203E-02_JPRB,1.7946E-02_JPRB,1.9828E-02_JPRB, &
     &2.1566E-02_JPRB,1.9707E-02_JPRB,2.0944E-02_JPRB,2.0944E-02_JPRB/)
      FORREFO(3,:) = (/ &
     &9.2770E-03_JPRB,1.2818E-02_JPRB,1.7181E-02_JPRB,1.7858E-02_JPRB,1.7888E-02_JPRB,1.7121E-02_JPRB, &
     &1.8116E-02_JPRB,1.8230E-02_JPRB,1.7719E-02_JPRB,1.7833E-02_JPRB,1.8438E-02_JPRB,1.7995E-02_JPRB, &
     &2.0895E-02_JPRB,2.1525E-02_JPRB,2.0517E-02_JPRB,2.0954E-02_JPRB/)
      FORREFO(4,:) = (/ &
     &8.3290E-03_JPRB,1.3483E-02_JPRB,1.5432E-02_JPRB,2.0793E-02_JPRB,1.8404E-02_JPRB,1.7470E-02_JPRB, &
     &1.7253E-02_JPRB,1.7132E-02_JPRB,1.7119E-02_JPRB,1.7376E-02_JPRB,1.7030E-02_JPRB,1.6847E-02_JPRB, &
     &1.5562E-02_JPRB,1.6836E-02_JPRB,1.8746E-02_JPRB,2.1233E-02_JPRB/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      SELFREFO(:, 1) = (/ &
     & 2.41120E-01_JPRB, 2.27071E-01_JPRB, 2.13840E-01_JPRB, 2.01380E-01_JPRB, 1.89646E-01_JPRB, &
     & 1.78596E-01_JPRB, 1.68190E-01_JPRB, 1.58390E-01_JPRB, 1.49161E-01_JPRB, 1.40470E-01_JPRB/)
      SELFREFO(:, 2) = (/ &
     & 3.11156E-01_JPRB, 2.92249E-01_JPRB, 2.74490E-01_JPRB, 2.57810E-01_JPRB, 2.42144E-01_JPRB, &
     & 2.27430E-01_JPRB, 2.13610E-01_JPRB, 2.00630E-01_JPRB, 1.88439E-01_JPRB, 1.76988E-01_JPRB/)
      SELFREFO(:, 3) = (/ &
     & 3.37148E-01_JPRB, 3.17767E-01_JPRB, 2.99500E-01_JPRB, 2.82283E-01_JPRB, 2.66056E-01_JPRB, &
     & 2.50762E-01_JPRB, 2.36347E-01_JPRB, 2.22760E-01_JPRB, 2.09955E-01_JPRB, 1.97885E-01_JPRB/)
      SELFREFO(:, 4) = (/ &
     & 3.57139E-01_JPRB, 3.32763E-01_JPRB, 3.10050E-01_JPRB, 2.88888E-01_JPRB, 2.69170E-01_JPRB, &
     & 2.50798E-01_JPRB, 2.33680E-01_JPRB, 2.17730E-01_JPRB, 2.02869E-01_JPRB, 1.89022E-01_JPRB/)
      SELFREFO(:, 5) = (/ &
     & 3.60626E-01_JPRB, 3.35433E-01_JPRB, 3.12000E-01_JPRB, 2.90204E-01_JPRB, 2.69931E-01_JPRB, &
     & 2.51074E-01_JPRB, 2.33534E-01_JPRB, 2.17220E-01_JPRB, 2.02045E-01_JPRB, 1.87931E-01_JPRB/)
      SELFREFO(:, 6) = (/ &
     & 3.42420E-01_JPRB, 3.18795E-01_JPRB, 2.96800E-01_JPRB, 2.76323E-01_JPRB, 2.57258E-01_JPRB, &
     & 2.39509E-01_JPRB, 2.22985E-01_JPRB, 2.07600E-01_JPRB, 1.93277E-01_JPRB, 1.79942E-01_JPRB/)
      SELFREFO(:, 7) = (/ &
     & 3.65491E-01_JPRB, 3.41599E-01_JPRB, 3.19270E-01_JPRB, 2.98400E-01_JPRB, 2.78895E-01_JPRB, &
     & 2.60664E-01_JPRB, 2.43625E-01_JPRB, 2.27700E-01_JPRB, 2.12816E-01_JPRB, 1.98905E-01_JPRB/)
      SELFREFO(:, 8) = (/ &
     & 3.70354E-01_JPRB, 3.45005E-01_JPRB, 3.21390E-01_JPRB, 2.99392E-01_JPRB, 2.78899E-01_JPRB, &
     & 2.59809E-01_JPRB, 2.42026E-01_JPRB, 2.25460E-01_JPRB, 2.10028E-01_JPRB, 1.95652E-01_JPRB/)
      SELFREFO(:, 9) = (/ &
     & 3.60483E-01_JPRB, 3.37846E-01_JPRB, 3.16630E-01_JPRB, 2.96747E-01_JPRB, 2.78112E-01_JPRB, &
     & 2.60648E-01_JPRB, 2.44280E-01_JPRB, 2.28940E-01_JPRB, 2.14563E-01_JPRB, 2.01090E-01_JPRB/)
      SELFREFO(:,10) = (/ &
     & 3.71845E-01_JPRB, 3.48164E-01_JPRB, 3.25990E-01_JPRB, 3.05229E-01_JPRB, 2.85790E-01_JPRB, &
     & 2.67588E-01_JPRB, 2.50547E-01_JPRB, 2.34590E-01_JPRB, 2.19650E-01_JPRB, 2.05661E-01_JPRB/)
      SELFREFO(:,11) = (/ &
     & 3.60606E-01_JPRB, 3.40789E-01_JPRB, 3.22060E-01_JPRB, 3.04361E-01_JPRB, 2.87634E-01_JPRB, &
     & 2.71826E-01_JPRB, 2.56888E-01_JPRB, 2.42770E-01_JPRB, 2.29428E-01_JPRB, 2.16819E-01_JPRB/)
      SELFREFO(:,12) = (/ &
     & 3.90046E-01_JPRB, 3.68879E-01_JPRB, 3.48860E-01_JPRB, 3.29928E-01_JPRB, 3.12023E-01_JPRB, &
     & 2.95089E-01_JPRB, 2.79075E-01_JPRB, 2.63930E-01_JPRB, 2.49607E-01_JPRB, 2.36061E-01_JPRB/)
      SELFREFO(:,13) = (/ &
     & 4.38542E-01_JPRB, 4.05139E-01_JPRB, 3.74280E-01_JPRB, 3.45771E-01_JPRB, 3.19434E-01_JPRB, &
     & 2.95103E-01_JPRB, 2.72626E-01_JPRB, 2.51860E-01_JPRB, 2.32676E-01_JPRB, 2.14953E-01_JPRB/)
      SELFREFO(:,14) = (/ &
     & 4.19448E-01_JPRB, 3.81920E-01_JPRB, 3.47750E-01_JPRB, 3.16637E-01_JPRB, 2.88307E-01_JPRB, &
     & 2.62513E-01_JPRB, 2.39026E-01_JPRB, 2.17640E-01_JPRB, 1.98168E-01_JPRB, 1.80438E-01_JPRB/)
      SELFREFO(:,15) = (/ &
     & 4.20276E-01_JPRB, 3.92281E-01_JPRB, 3.66150E-01_JPRB, 3.41760E-01_JPRB, 3.18995E-01_JPRB, &
     & 2.97746E-01_JPRB, 2.77912E-01_JPRB, 2.59400E-01_JPRB, 2.42121E-01_JPRB, 2.25993E-01_JPRB/)
      SELFREFO(:,16) = (/ &
     & 4.20276E-01_JPRB, 3.92281E-01_JPRB, 3.66150E-01_JPRB, 3.41760E-01_JPRB, 3.18995E-01_JPRB, &
     & 2.97746E-01_JPRB, 2.77912E-01_JPRB, 2.59400E-01_JPRB, 2.42121E-01_JPRB, 2.25993E-01_JPRB/)

IF (LHOOK) CALL DR_HOOK('RRTM_KGB10',1,ZHOOK_HANDLE)
RETURN

1001 CONTINUE
CALL ABOR1("RRTM_KGB10:ERROR READING FILE RADRRTM")

END SUBROUTINE RRTM_KGB10
