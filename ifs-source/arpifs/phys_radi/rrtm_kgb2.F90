SUBROUTINE RRTM_KGB2

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     R. Elkhatib 12-10-2005 Split for faster and more robust compilation.
!     G.Mozdzynski March 2011 read constants from files
!     ABozzo May 2013 update to RRTMG v4.85
!     band 2:  350-500 cm-1 
!     T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN    ,ONLY : NULRAD
USE MPL_MODULE,ONLY : MPL_BROADCAST
USE YOMTAG    ,ONLY : MTAGRAD

USE YOERRTO2 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,&
 & FRACREFBO  ,FORREFO  ,KAO_D, KBO_D
USE YOMMP0    , ONLY : NPROC, MYPROC

!     ------------------------------------------------------------------

IMPLICIT NONE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('RRTM_KGB2',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  READ(NULRAD,ERR=1001) KAO_D,KBO_D
  KAO = REAL(KAO_D,JPRB)
  KBO = REAL(KBO_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KAO,MTAGRAD,1,CDSTRING='RRTM_KGB2:')
  CALL MPL_BROADCAST (KBO,MTAGRAD,1,CDSTRING='RRTM_KGB2:')
ENDIF


! Planck fraction mapping level: P = 1053.630 mbar, T = 294.2 K
      FRACREFAO(:) = (/ &
      &  1.6388E-01_JPRB, 1.5241E-01_JPRB, 1.4290E-01_JPRB, 1.2864E-01_JPRB, &
      &  1.1615E-01_JPRB, 1.0047E-01_JPRB, 8.0013E-02_JPRB, 6.0445E-02_JPRB, &
      &  4.0530E-02_JPRB, 4.3879E-03_JPRB, 3.5726E-03_JPRB, 2.7669E-03_JPRB, &
      &  2.0078E-03_JPRB, 1.2864E-03_JPRB, 4.7630E-04_JPRB, 6.9109E-05_JPRB/)

! Planck fraction mapping level: P = 3.206e-2 mb, T = 197.92 K
      FRACREFBO(:) = (/ &
      &  1.4697E-01_JPRB, 1.4826E-01_JPRB, 1.4278E-01_JPRB, 1.3320E-01_JPRB, &
      &  1.1965E-01_JPRB, 1.0297E-01_JPRB, 8.4170E-02_JPRB, 6.3282E-02_JPRB, &
      &  4.2868E-02_JPRB, 4.6644E-03_JPRB, 3.8619E-03_JPRB, 3.0533E-03_JPRB, &
      &  2.2359E-03_JPRB, 1.4226E-03_JPRB, 5.3642E-04_JPRB, 7.6316E-05_JPRB/)

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      FORREFO(1,:) = (/ &
     & 2.8549E-03_JPRB,4.8281E-03_JPRB,6.2570E-03_JPRB,8.2731E-03_JPRB,7.9056E-03_JPRB,7.7840E-03_JPRB, &
     & 1.0115E-02_JPRB,9.6599E-03_JPRB,1.0153E-02_JPRB,1.0921E-02_JPRB,1.2408E-02_JPRB,1.3496E-02_JPRB, &
     & 1.5059E-02_JPRB,1.4636E-02_JPRB,1.6483E-02_JPRB,1.2394E-02_JPRB/)
      FORREFO(2,:) = (/ &
     & 3.0036E-03_JPRB,5.1093E-03_JPRB,5.7317E-03_JPRB,9.2246E-03_JPRB,8.9829E-03_JPRB,8.6477E-03_JPRB, &
     & 1.1448E-02_JPRB,1.0391E-02_JPRB,1.0211E-02_JPRB,1.2921E-02_JPRB,1.2726E-02_JPRB,1.2426E-02_JPRB, &
     & 1.4609E-02_JPRB,1.5783E-02_JPRB,1.6617E-02_JPRB,1.6858E-02_JPRB/)
      FORREFO(3,:) = (/ &
     & 3.0771E-03_JPRB,5.1206E-03_JPRB,5.8426E-03_JPRB,9.5727E-03_JPRB,1.0338E-02_JPRB,9.3737E-03_JPRB, &
     & 1.2805E-02_JPRB,1.1272E-02_JPRB,1.1353E-02_JPRB,1.1837E-02_JPRB,1.1550E-02_JPRB,1.3020E-02_JPRB, &
     & 1.3536E-02_JPRB,1.6226E-02_JPRB,1.6039E-02_JPRB,2.2578E-02_JPRB/)
      FORREFO(4,:) = (/ &
     & 3.3072E-03_JPRB,5.0240E-03_JPRB,6.8474E-03_JPRB,8.2736E-03_JPRB,8.6151E-03_JPRB,8.6762E-03_JPRB, &
     & 1.1476E-02_JPRB,1.0246E-02_JPRB,1.0819E-02_JPRB,1.0640E-02_JPRB,1.0545E-02_JPRB,1.0533E-02_JPRB, &
     & 1.0496E-02_JPRB,1.0142E-02_JPRB,9.7979E-03_JPRB,1.5255E-02_JPRB/)


!     The following are parameters related to the reference water vapor
!     mixing ratios by REFPARAM(I) = REFH2O(I) / (.002+REFH2O(I)).
!     These parameters are used for the Planck function interpolation.
!REFPARAM( :) = (/&
! & 0.903661_JPRB   , 0.859386_JPRB   , 0.746542_JPRB   , 0.580496_JPRB   , 0.412889_JPRB   ,&
! & 0.275283_JPRB   , 0.162745_JPRB   , 7.63929E-02_JPRB, 1.82553E-02_JPRB, 3.72432E-03_JPRB, &
! & 2.14946E-03_JPRB, 1.66320E-03_JPRB, 1.59940E-03_JPRB/)  

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



!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

    SELFREFO(:, 1) = (/ &
     & 7.25695E-01_JPRB, 6.53591E-01_JPRB, 5.88650E-01_JPRB, 5.30162E-01_JPRB, 4.77485E-01_JPRB, &
     & 4.30042E-01_JPRB, 3.87313E-01_JPRB, 3.48830E-01_JPRB, 3.14170E-01_JPRB, 2.82954E-01_JPRB/)
      SELFREFO(:, 2) = (/ &
     & 9.61996E-01_JPRB, 8.77853E-01_JPRB, 8.01070E-01_JPRB, 7.31003E-01_JPRB, 6.67064E-01_JPRB, &
     & 6.08718E-01_JPRB, 5.55476E-01_JPRB, 5.06890E-01_JPRB, 4.62554E-01_JPRB, 4.22096E-01_JPRB/)
      SELFREFO(:, 3) = (/ &
     & 9.72584E-01_JPRB, 9.02658E-01_JPRB, 8.37760E-01_JPRB, 7.77527E-01_JPRB, 7.21626E-01_JPRB, &
     & 6.69743E-01_JPRB, 6.21591E-01_JPRB, 5.76900E-01_JPRB, 5.35423E-01_JPRB, 4.96927E-01_JPRB/)
      SELFREFO(:, 4) = (/ &
     & 1.24790E+00_JPRB, 1.14353E+00_JPRB, 1.04790E+00_JPRB, 9.60263E-01_JPRB, 8.79956E-01_JPRB, &
     & 8.06364E-01_JPRB, 7.38927E-01_JPRB, 6.77130E-01_JPRB, 6.20501E-01_JPRB, 5.68608E-01_JPRB/)
      SELFREFO(:, 5) = (/ &
     & 1.23574E+00_JPRB, 1.12928E+00_JPRB, 1.03200E+00_JPRB, 9.43096E-01_JPRB, 8.61851E-01_JPRB, &
     & 7.87605E-01_JPRB, 7.19755E-01_JPRB, 6.57750E-01_JPRB, 6.01087E-01_JPRB, 5.49305E-01_JPRB/)
      SELFREFO(:, 6) = (/ &
     & 1.20921E+00_JPRB, 1.10660E+00_JPRB, 1.01270E+00_JPRB, 9.26766E-01_JPRB, 8.48124E-01_JPRB, &
     & 7.76155E-01_JPRB, 7.10293E-01_JPRB, 6.50020E-01_JPRB, 5.94861E-01_JPRB, 5.44384E-01_JPRB/)
      SELFREFO(:, 7) = (/ &
     & 1.38112E+00_JPRB, 1.26727E+00_JPRB, 1.16280E+00_JPRB, 1.06694E+00_JPRB, 9.78990E-01_JPRB, &
     & 8.98287E-01_JPRB, 8.24236E-01_JPRB, 7.56290E-01_JPRB, 6.93945E-01_JPRB, 6.36739E-01_JPRB/)
      SELFREFO(:, 8) = (/ &
     & 1.30321E+00_JPRB, 1.20127E+00_JPRB, 1.10730E+00_JPRB, 1.02068E+00_JPRB, 9.40840E-01_JPRB, &
     & 8.67243E-01_JPRB, 7.99403E-01_JPRB, 7.36870E-01_JPRB, 6.79229E-01_JPRB, 6.26096E-01_JPRB/)
      SELFREFO(:, 9) = (/ &
     & 1.26713E+00_JPRB, 1.17927E+00_JPRB, 1.09750E+00_JPRB, 1.02140E+00_JPRB, 9.50575E-01_JPRB, &
     & 8.84662E-01_JPRB, 8.23319E-01_JPRB, 7.66230E-01_JPRB, 7.13099E-01_JPRB, 6.63653E-01_JPRB/)
      SELFREFO(:,10) = (/ &
     & 1.49824E+00_JPRB, 1.37053E+00_JPRB, 1.25370E+00_JPRB, 1.14683E+00_JPRB, 1.04908E+00_JPRB, &
     & 9.59651E-01_JPRB, 8.77849E-01_JPRB, 8.03020E-01_JPRB, 7.34569E-01_JPRB, 6.71954E-01_JPRB/)
      SELFREFO(:,11) = (/ &
     & 1.44786E+00_JPRB, 1.34594E+00_JPRB, 1.25120E+00_JPRB, 1.16313E+00_JPRB, 1.08125E+00_JPRB, &
     & 1.00514E+00_JPRB, 9.34392E-01_JPRB, 8.68620E-01_JPRB, 8.07477E-01_JPRB, 7.50639E-01_JPRB/)
      SELFREFO(:,12) = (/ &
     & 1.38460E+00_JPRB, 1.30437E+00_JPRB, 1.22880E+00_JPRB, 1.15760E+00_JPRB, 1.09053E+00_JPRB, &
     & 1.02735E+00_JPRB, 9.67825E-01_JPRB, 9.11750E-01_JPRB, 8.58924E-01_JPRB, 8.09159E-01_JPRB/)
      SELFREFO(:,13) = (/ &
     & 1.51953E+00_JPRB, 1.42822E+00_JPRB, 1.34240E+00_JPRB, 1.26173E+00_JPRB, 1.18592E+00_JPRB, &
     & 1.11465E+00_JPRB, 1.04768E+00_JPRB, 9.84720E-01_JPRB, 9.25548E-01_JPRB, 8.69932E-01_JPRB/)
      SELFREFO(:,14) = (/ &
     & 1.62608E+00_JPRB, 1.51021E+00_JPRB, 1.40260E+00_JPRB, 1.30266E+00_JPRB, 1.20983E+00_JPRB, &
     & 1.12363E+00_JPRB, 1.04356E+00_JPRB, 9.69200E-01_JPRB, 9.00138E-01_JPRB, 8.35998E-01_JPRB/)
      SELFREFO(:,15) = (/ &
     & 1.65383E+00_JPRB, 1.54808E+00_JPRB, 1.44910E+00_JPRB, 1.35644E+00_JPRB, 1.26971E+00_JPRB, &
     & 1.18853E+00_JPRB, 1.11254E+00_JPRB, 1.04140E+00_JPRB, 9.74813E-01_JPRB, 9.12484E-01_JPRB/)
      SELFREFO(:,16) = (/ &
     & 1.78105E+00_JPRB, 1.61421E+00_JPRB, 1.46300E+00_JPRB, 1.32595E+00_JPRB, 1.20174E+00_JPRB, &
     & 1.08917E+00_JPRB, 9.87141E-01_JPRB, 8.94670E-01_JPRB, 8.10861E-01_JPRB, 7.34904E-01_JPRB/)

IF (LHOOK) CALL DR_HOOK('RRTM_KGB2',1,ZHOOK_HANDLE)
RETURN

1001 CONTINUE
CALL ABOR1("RRTM_KGB2:ERROR READING FILE RADRRTM")

END SUBROUTINE RRTM_KGB2
