SUBROUTINE RRTM_KGB16

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)
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

USE YOERRTO16, ONLY : KAO,KBO ,SELFREFO,FORREFO ,FRACREFAO,FRACREFBO,KAO_D,KBO_D
USE YOMMP0    , ONLY : NPROC, MYPROC


!     ------------------------------------------------------------------

IMPLICIT NONE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('RRTM_KGB16',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  READ(NULRAD,ERR=1001) KAO_D,KBO_D
  CLOSE(NULRAD,ERR=1000)
  KAO = REAL(KAO_D,JPRB)
  KBO = REAL(KBO_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KAO,MTAGRAD,1,CDSTRING='RRTM_KGB16:')
  CALL MPL_BROADCAST (KBO,MTAGRAD,1,CDSTRING='RRTM_KGB16:')
ENDIF

! Planck fraction mapping level: P = 387.6100 mbar, T = 250.17 K
      FRACREFAO(:, 1) = (/ &
     &  1.1593E-01_JPRB,2.3390E-01_JPRB,1.9120E-01_JPRB,1.3121E-01_JPRB,1.0590E-01_JPRB,8.4852E-02_JPRB, &
     &  6.4168E-02_JPRB,4.2537E-02_JPRB,2.3220E-02_JPRB,2.1767E-03_JPRB,1.8203E-03_JPRB,1.3724E-03_JPRB, &
     &  9.5452E-04_JPRB,5.5015E-04_JPRB,1.9348E-04_JPRB,2.7344E-05_JPRB/)
      FRACREFAO(:, 2) = (/ &
     &  2.8101E-01_JPRB,1.9773E-01_JPRB,1.4749E-01_JPRB,1.1399E-01_JPRB,8.8190E-02_JPRB,7.0531E-02_JPRB, &
     &  4.6356E-02_JPRB,3.0774E-02_JPRB,1.7332E-02_JPRB,2.0054E-03_JPRB,1.5950E-03_JPRB,1.2760E-03_JPRB, &
     &  9.5034E-04_JPRB,5.4992E-04_JPRB,1.9349E-04_JPRB,2.7309E-05_JPRB/)
      FRACREFAO(:, 3) = (/ &
     &  2.9054E-01_JPRB,2.1263E-01_JPRB,1.4133E-01_JPRB,1.1083E-01_JPRB,8.5107E-02_JPRB,6.5247E-02_JPRB, &
     &  4.4542E-02_JPRB,2.7205E-02_JPRB,1.6495E-02_JPRB,1.8453E-03_JPRB,1.5222E-03_JPRB,1.1884E-03_JPRB, &
     &  8.1094E-04_JPRB,4.9173E-04_JPRB,1.9344E-04_JPRB,2.7286E-05_JPRB/)
      FRACREFAO(:, 4) = (/ &
     &  2.9641E-01_JPRB,2.1738E-01_JPRB,1.4228E-01_JPRB,1.0830E-01_JPRB,8.2837E-02_JPRB,6.1359E-02_JPRB, &
     &  4.4683E-02_JPRB,2.5027E-02_JPRB,1.6057E-02_JPRB,1.7558E-03_JPRB,1.4193E-03_JPRB,1.0970E-03_JPRB, &
     &  7.8281E-04_JPRB,4.3260E-04_JPRB,1.4837E-04_JPRB,2.2958E-05_JPRB/)
      FRACREFAO(:, 5) = (/ &
     &  2.9553E-01_JPRB,2.2139E-01_JPRB,1.4816E-01_JPRB,1.0601E-01_JPRB,8.0048E-02_JPRB,6.0082E-02_JPRB, &
     &  4.3952E-02_JPRB,2.3788E-02_JPRB,1.5734E-02_JPRB,1.6586E-03_JPRB,1.3434E-03_JPRB,1.0281E-03_JPRB, &
     &  7.0256E-04_JPRB,4.2577E-04_JPRB,1.2803E-04_JPRB,1.3315E-05_JPRB/)
      FRACREFAO(:, 6) = (/ &
     &  2.9313E-01_JPRB,2.2476E-01_JPRB,1.5470E-01_JPRB,1.0322E-01_JPRB,7.8904E-02_JPRB,5.8175E-02_JPRB, &
     &  4.3097E-02_JPRB,2.3618E-02_JPRB,1.5385E-02_JPRB,1.5942E-03_JPRB,1.2702E-03_JPRB,9.5566E-04_JPRB, &
     &  6.5421E-04_JPRB,4.0165E-04_JPRB,1.2805E-04_JPRB,1.3355E-05_JPRB/)
      FRACREFAO(:, 7) = (/ &
     &  2.9069E-01_JPRB,2.2823E-01_JPRB,1.5995E-01_JPRB,1.0170E-01_JPRB,7.7287E-02_JPRB,5.6780E-02_JPRB, &
     &  4.1752E-02_JPRB,2.3899E-02_JPRB,1.4937E-02_JPRB,1.4916E-03_JPRB,1.1909E-03_JPRB,9.1307E-04_JPRB, &
     &  6.3518E-04_JPRB,3.9866E-04_JPRB,1.2805E-04_JPRB,1.3298E-05_JPRB/)
      FRACREFAO(:, 8) = (/ &
     &  2.8446E-01_JPRB,2.2651E-01_JPRB,1.7133E-01_JPRB,1.0299E-01_JPRB,7.4231E-02_JPRB,5.6031E-02_JPRB, &
     &  4.1368E-02_JPRB,2.4318E-02_JPRB,1.4135E-02_JPRB,1.4216E-03_JPRB,1.1465E-03_JPRB,8.9800E-04_JPRB, &
     &  6.3553E-04_JPRB,3.9536E-04_JPRB,1.2749E-04_JPRB,1.3298E-05_JPRB/)
      FRACREFAO(:, 9) = (/ &
     &  2.0568E-01_JPRB,2.5049E-01_JPRB,2.0568E-01_JPRB,1.1781E-01_JPRB,7.5579E-02_JPRB,5.8136E-02_JPRB, &
     &  4.2397E-02_JPRB,2.6544E-02_JPRB,1.3067E-02_JPRB,1.4061E-03_JPRB,1.1455E-03_JPRB,8.9408E-04_JPRB, &
     &  6.3652E-04_JPRB,3.9450E-04_JPRB,1.2841E-04_JPRB,1.3315E-05_JPRB/)

! Planck fraction mapping level : P=95.58350 mb, T = 215.70 K
      FRACREFBO(:) = (/ &
     &  1.8111E-01_JPRB,2.2612E-01_JPRB,1.6226E-01_JPRB,1.1872E-01_JPRB,9.9048E-02_JPRB,8.0390E-02_JPRB, &
     &  6.1648E-02_JPRB,4.1704E-02_JPRB,2.2976E-02_JPRB,1.9263E-03_JPRB,1.4694E-03_JPRB,1.1498E-03_JPRB, &
     &  7.9906E-04_JPRB,4.8310E-04_JPRB,1.6188E-04_JPRB,2.2651E-05_JPRB/)


!     ------------------------------------------------------------------

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  The first 
!     index in the array, JS, runs from 1 to 9, and corresponds to 
!     different values of the binary species parameter.  For instance, 
!     JS=1 refers to dry air, JS = 2 corresponds to the paramter value 1/8, 
!     JS = 3 corresponds to the parameter value 2/8, etc.  The second index
!     in the array, JT, which runs from 1 to 5, corresponds to different
!     temperatures.  More specifically, JT = 3 means that the data are for
!     the reference temperature TREF for this  pressure level, JT = 2 refers
!     to TREF-15, JT = 1 is for TREF-30, JT = 4 is for TREF+15, and JT = 5
!     is for TREF+30.  The third index, JP, runs from 1 to 13 and refers
!     to the JPth reference pressure level (see taumol.f for these levels
!     in mb).  The fourth index, IG, goes from 1 to 16, and indicates
!     which g-interval the absorption coefficients are for.


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
     &5.1629E-06_JPRB,7.7578E-06_JPRB,1.9043E-05_JPRB,1.4802E-04_JPRB,2.2980E-04_JPRB,2.8057E-04_JPRB, &
     &3.2824E-04_JPRB,3.4913E-04_JPRB,3.6515E-04_JPRB,3.8271E-04_JPRB,3.7499E-04_JPRB,3.6966E-04_JPRB, &
     &3.7424E-04_JPRB,3.8884E-04_JPRB,3.7117E-04_JPRB,4.3710E-04_JPRB/)
      FORREFO(2,:) = (/ &
     &5.0804E-06_JPRB,1.3466E-05_JPRB,7.2606E-05_JPRB,1.6940E-04_JPRB,2.1022E-04_JPRB,2.5900E-04_JPRB, &
     &2.9106E-04_JPRB,3.2261E-04_JPRB,3.2066E-04_JPRB,3.5421E-04_JPRB,3.7128E-04_JPRB,3.8144E-04_JPRB, &
     &3.7854E-04_JPRB,3.8347E-04_JPRB,3.8921E-04_JPRB,3.7339E-04_JPRB/)
      FORREFO(3,:) = (/ &
     &5.4797E-05_JPRB,1.0026E-04_JPRB,1.2422E-04_JPRB,1.6386E-04_JPRB,1.8378E-04_JPRB,1.9616E-04_JPRB, &
     &2.0711E-04_JPRB,2.2492E-04_JPRB,2.5240E-04_JPRB,2.6187E-04_JPRB,2.6058E-04_JPRB,2.4892E-04_JPRB, &
     &2.6526E-04_JPRB,3.2105E-04_JPRB,3.6903E-04_JPRB,3.7213E-04_JPRB/)
      FORREFO(4,:) = (/ &
     &4.2782E-05_JPRB,1.4775E-04_JPRB,1.4588E-04_JPRB,1.6964E-04_JPRB,1.6667E-04_JPRB,1.7192E-04_JPRB, &
     &1.9057E-04_JPRB,2.0180E-04_JPRB,2.1177E-04_JPRB,2.2326E-04_JPRB,2.3801E-04_JPRB,2.9308E-04_JPRB, &
     &3.1130E-04_JPRB,3.1829E-04_JPRB,3.5035E-04_JPRB,3.7782E-04_JPRB/)


!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      SELFREFO(:, 1) = (/ &
     & 1.27793E-03_JPRB, 1.05944E-03_JPRB, 8.78300E-04_JPRB, 7.28133E-04_JPRB, 6.03641E-04_JPRB, &
     & 5.00434E-04_JPRB, 4.14873E-04_JPRB, 3.43940E-04_JPRB, 2.85135E-04_JPRB, 2.36384E-04_JPRB/)
      SELFREFO(:, 2) = (/ &
     & 1.42785E-03_JPRB, 1.17602E-03_JPRB, 9.68600E-04_JPRB, 7.97765E-04_JPRB, 6.57060E-04_JPRB, &
     & 5.41172E-04_JPRB, 4.45724E-04_JPRB, 3.67110E-04_JPRB, 3.02361E-04_JPRB, 2.49033E-04_JPRB/)
      SELFREFO(:, 3) = (/ &
     & 2.94095E-03_JPRB, 2.27102E-03_JPRB, 1.75370E-03_JPRB, 1.35422E-03_JPRB, 1.04574E-03_JPRB, &
     & 8.07525E-04_JPRB, 6.23577E-04_JPRB, 4.81530E-04_JPRB, 3.71841E-04_JPRB, 2.87138E-04_JPRB/)
      SELFREFO(:, 4) = (/ &
     & 3.94894E-03_JPRB, 3.48184E-03_JPRB, 3.07000E-03_JPRB, 2.70687E-03_JPRB, 2.38669E-03_JPRB, &
     & 2.10439E-03_JPRB, 1.85547E-03_JPRB, 1.63600E-03_JPRB, 1.44249E-03_JPRB, 1.27187E-03_JPRB/)
      SELFREFO(:, 5) = (/ &
     & 4.19971E-03_JPRB, 3.86333E-03_JPRB, 3.55390E-03_JPRB, 3.26925E-03_JPRB, 3.00740E-03_JPRB, &
     & 2.76652E-03_JPRB, 2.54494E-03_JPRB, 2.34110E-03_JPRB, 2.15359E-03_JPRB, 1.98110E-03_JPRB/)
      SELFREFO(:, 6) = (/ &
     & 4.95922E-03_JPRB, 4.57134E-03_JPRB, 4.21380E-03_JPRB, 3.88422E-03_JPRB, 3.58042E-03_JPRB, &
     & 3.30038E-03_JPRB, 3.04225E-03_JPRB, 2.80430E-03_JPRB, 2.58496E-03_JPRB, 2.38278E-03_JPRB/)
      SELFREFO(:, 7) = (/ &
     & 5.27379E-03_JPRB, 4.91005E-03_JPRB, 4.57140E-03_JPRB, 4.25611E-03_JPRB, 3.96256E-03_JPRB, &
     & 3.68925E-03_JPRB, 3.43480E-03_JPRB, 3.19790E-03_JPRB, 2.97734E-03_JPRB, 2.77199E-03_JPRB/)
      SELFREFO(:, 8) = (/ &
     & 5.75341E-03_JPRB, 5.31533E-03_JPRB, 4.91060E-03_JPRB, 4.53669E-03_JPRB, 4.19126E-03_JPRB, &
     & 3.87212E-03_JPRB, 3.57729E-03_JPRB, 3.30490E-03_JPRB, 3.05325E-03_JPRB, 2.82077E-03_JPRB/)
      SELFREFO(:, 9) = (/ &
     & 5.49849E-03_JPRB, 5.14295E-03_JPRB, 4.81040E-03_JPRB, 4.49935E-03_JPRB, 4.20842E-03_JPRB, &
     & 3.93629E-03_JPRB, 3.68177E-03_JPRB, 3.44370E-03_JPRB, 3.22102E-03_JPRB, 3.01275E-03_JPRB/)
      SELFREFO(:,10) = (/ &
     & 6.04962E-03_JPRB, 5.60945E-03_JPRB, 5.20130E-03_JPRB, 4.82285E-03_JPRB, 4.47194E-03_JPRB, &
     & 4.14656E-03_JPRB, 3.84485E-03_JPRB, 3.56510E-03_JPRB, 3.30570E-03_JPRB, 3.06518E-03_JPRB/)
      SELFREFO(:,11) = (/ &
     & 6.40108E-03_JPRB, 5.87551E-03_JPRB, 5.39310E-03_JPRB, 4.95029E-03_JPRB, 4.54385E-03_JPRB, &
     & 4.17077E-03_JPRB, 3.82833E-03_JPRB, 3.51400E-03_JPRB, 3.22548E-03_JPRB, 2.96065E-03_JPRB/)
      SELFREFO(:,12) = (/ &
     & 6.77938E-03_JPRB, 6.15713E-03_JPRB, 5.59200E-03_JPRB, 5.07874E-03_JPRB, 4.61259E-03_JPRB, &
     & 4.18922E-03_JPRB, 3.80472E-03_JPRB, 3.45550E-03_JPRB, 3.13834E-03_JPRB, 2.85029E-03_JPRB/)
      SELFREFO(:,13) = (/ &
     & 6.90020E-03_JPRB, 6.26766E-03_JPRB, 5.69310E-03_JPRB, 5.17121E-03_JPRB, 4.69717E-03_JPRB, &
     & 4.26658E-03_JPRB, 3.87546E-03_JPRB, 3.52020E-03_JPRB, 3.19750E-03_JPRB, 2.90439E-03_JPRB/)
      SELFREFO(:,14) = (/ &
     & 6.92759E-03_JPRB, 6.32882E-03_JPRB, 5.78180E-03_JPRB, 5.28206E-03_JPRB, 4.82552E-03_JPRB, &
     & 4.40843E-03_JPRB, 4.02740E-03_JPRB, 3.67930E-03_JPRB, 3.36129E-03_JPRB, 3.07076E-03_JPRB/)
      SELFREFO(:,15) = (/ &
     & 7.54539E-03_JPRB, 6.81161E-03_JPRB, 6.14920E-03_JPRB, 5.55120E-03_JPRB, 5.01136E-03_JPRB, &
     & 4.52402E-03_JPRB, 4.08407E-03_JPRB, 3.68690E-03_JPRB, 3.32836E-03_JPRB, 3.00468E-03_JPRB/)
      SELFREFO(:,16) = (/ &
     & 7.62039E-03_JPRB, 7.10834E-03_JPRB, 6.63070E-03_JPRB, 6.18515E-03_JPRB, 5.76955E-03_JPRB, &
     & 5.38186E-03_JPRB, 5.02023E-03_JPRB, 4.68290E-03_JPRB, 4.36823E-03_JPRB, 4.07471E-03_JPRB/)

IF (LHOOK) CALL DR_HOOK('RRTM_KGB16',1,ZHOOK_HANDLE)
RETURN

1000 CONTINUE
CALL ABOR1("RRTM_KGB16:ERROR CLOSING FILE RADRRTM")
1001 CONTINUE
CALL ABOR1("RRTM_KGB16:ERROR READING FILE RADRRTM")

END SUBROUTINE RRTM_KGB16
