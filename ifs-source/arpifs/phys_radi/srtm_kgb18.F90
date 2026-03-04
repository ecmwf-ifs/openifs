SUBROUTINE SRTM_KGB18

!     Originally by J.Delamere, Atmospheric & Environmental Research.
!     Revision: 2.4
!     BAND 18:  4000-4650 cm-1 (low - H2O,CH4; high - CH4)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     R. Elkhatib 12-10-2005 Split for faster and more robust compilation.
!     G.Mozdzynski March 2011 read constants from files
!     T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN    , ONLY : NULRAD
USE YOMMP0    , ONLY : NPROC, MYPROC
USE MPL_MODULE, ONLY : MPL_BROADCAST
USE YOMTAG    , ONLY : MTAGRAD
USE YOESRTA18 , ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, RAYL, STRRAT, LAYREFFR, &
   &  KA_D, KB_D 

!     ------------------------------------------------------------------

IMPLICIT NONE

! KURUCZ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SRTM_KGB18',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  READ(NULRAD,ERR=1001) KA_D,KB_D
  KA = REAL(KA_D,JPRB)
  KB = REAL(KB_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KA,MTAGRAD,1,CDSTRING='SRTM_KGB18:')
  CALL MPL_BROADCAST (KB,MTAGRAD,1,CDSTRING='SRTM_KGB18:')
ENDIF

SFLUXREF(:,1) = (/ &
 & 3.65840_JPRB    , 3.54375_JPRB    , 3.34481_JPRB    , 3.10534_JPRB    , &
 & 2.79879_JPRB    , 2.42841_JPRB    , 1.98748_JPRB    , 1.49377_JPRB    , &
 & 1.00196_JPRB    , 0.108342_JPRB   , 8.95099E-02_JPRB, 7.05199E-02_JPRB, &
 & 5.16432E-02_JPRB, 3.27635E-02_JPRB, 1.25133E-02_JPRB, 1.73001E-03_JPRB /)    
SFLUXREF(:,2) = (/ &
 & 3.86372_JPRB    , 3.48521_JPRB    , 3.30790_JPRB    , 3.08103_JPRB    , &
 & 2.77552_JPRB    , 2.40722_JPRB    , 1.97307_JPRB    , 1.48023_JPRB    , &
 & 0.993055_JPRB   , 0.107691_JPRB   , 8.84430E-02_JPRB, 6.99354E-02_JPRB, &
 & 5.07881E-02_JPRB, 3.24121E-02_JPRB, 1.19442E-02_JPRB, 1.57612E-03_JPRB /)  
SFLUXREF(:,3) = (/ &
 & 3.90370_JPRB    , 3.50657_JPRB    , 3.30629_JPRB    , 3.06046_JPRB    , &
 & 2.76982_JPRB    , 2.39907_JPRB    , 1.96358_JPRB    , 1.47458_JPRB    , &
 & 0.988475_JPRB   , 0.106698_JPRB   , 8.75242E-02_JPRB, 6.85898E-02_JPRB, &
 & 5.04798E-02_JPRB, 3.13718E-02_JPRB, 1.09533E-02_JPRB, 1.57612E-03_JPRB /)  
SFLUXREF(:,4) = (/ &
 & 3.93165_JPRB    , 3.52058_JPRB    , 3.31346_JPRB    , 3.04944_JPRB    , &
 & 2.76074_JPRB    , 2.39433_JPRB    , 1.95556_JPRB    , 1.46712_JPRB    , &
 & 0.984056_JPRB   , 0.105885_JPRB   , 8.73062E-02_JPRB, 6.84054E-02_JPRB, &
 & 4.87443E-02_JPRB, 2.99295E-02_JPRB, 1.09533E-02_JPRB, 1.57612E-03_JPRB /)  
SFLUXREF(:,5) = (/ &
 & 3.94082_JPRB    , 3.55221_JPRB    , 3.31863_JPRB    , 3.04730_JPRB    , &
 & 2.74918_JPRB    , 2.38328_JPRB    , 1.95212_JPRB    , 1.45889_JPRB    , &
 & 0.978888_JPRB   , 0.105102_JPRB   , 8.65732E-02_JPRB, 6.74563E-02_JPRB, &
 & 4.76592E-02_JPRB, 2.91017E-02_JPRB, 1.09533E-02_JPRB, 1.57612E-03_JPRB /)  
SFLUXREF(:,6) = (/ &
 & 3.94198_JPRB    , 3.58743_JPRB    , 3.32106_JPRB    , 3.05866_JPRB    , &
 & 2.74115_JPRB    , 2.36939_JPRB    , 1.94305_JPRB    , 1.45180_JPRB    , &
 & 0.971784_JPRB   , 1.04045E-01_JPRB, 8.53731E-02_JPRB, 6.60654E-02_JPRB, &
 & 4.63228E-02_JPRB, 2.91016E-02_JPRB, 1.09552E-02_JPRB, 1.57612E-03_JPRB /)  
SFLUXREF(:,7) = (/ &
 & 3.93596_JPRB    , 3.63366_JPRB    , 3.33144_JPRB    , 3.06252_JPRB    , &
 & 2.74054_JPRB    , 2.35492_JPRB    , 1.92769_JPRB    , 1.44300_JPRB    , &
 & 0.961809_JPRB   , 1.02867E-01_JPRB, 8.34164E-02_JPRB, 6.41005E-02_JPRB, &
 & 4.61826E-02_JPRB, 2.91006E-02_JPRB, 1.09553E-02_JPRB, 1.57612E-03_JPRB /)  
SFLUXREF(:,8) = (/ &
 & 3.92520_JPRB    , 3.69078_JPRB    , 3.35656_JPRB    , 3.07055_JPRB    , &
 & 2.73862_JPRB    , 2.34430_JPRB    , 1.90187_JPRB    , 1.42242_JPRB    , &
 & 0.946676_JPRB   , 9.96302E-02_JPRB, 8.14421E-02_JPRB, 6.38622E-02_JPRB, &
 & 4.61794E-02_JPRB, 2.91017E-02_JPRB, 1.09553E-02_JPRB, 1.57612E-03_JPRB /)  
SFLUXREF(:,9) = (/ &
 & 3.80721_JPRB    , 3.74437_JPRB    , 3.50205_JPRB    , 3.18009_JPRB    , &
 & 2.75757_JPRB    , 2.29188_JPRB    , 1.84382_JPRB    , 1.35694_JPRB    , &
 & 0.914040_JPRB   , 9.86811E-02_JPRB, 8.14321E-02_JPRB, 6.38541E-02_JPRB, &
 & 4.61795E-02_JPRB, 2.90960E-02_JPRB, 1.09613E-02_JPRB, 1.57612E-03_JPRB /)  

!     Rayleigh extinction coefficient at v = 4325 cm-1.
RAYL = 1.39E-09_JPRB

STRRAT = 38.9589_JPRB

LAYREFFR = 6

!     ------------------------------------------------------------------

!     The array KA contains absorption coefs at the 16 chosen g-values 
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
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
!     The array KB contains absorption coefs at the 16 chosen g-values 
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
!     -----------------------------------------------------------------

  
FORREF(:, 1) = (/ 0.860560E-06_JPRB, 0.130439E-05_JPRB, 0.382378E-05_JPRB /)
FORREF(:, 2) = (/ 0.817926E-06_JPRB, 0.158599E-05_JPRB, 0.658771E-04_JPRB /)
FORREF(:, 3) = (/ 0.129369E-05_JPRB, 0.824406E-05_JPRB, 0.952778E-04_JPRB /)
FORREF(:, 4) = (/ 0.438918E-05_JPRB, 0.375356E-04_JPRB, 0.119111E-03_JPRB /)
FORREF(:, 5) = (/ 0.306057E-04_JPRB, 0.622798E-04_JPRB, 0.100740E-03_JPRB /)
FORREF(:, 6) = (/ 0.891934E-04_JPRB, 0.856393E-04_JPRB, 0.635583E-04_JPRB /)
FORREF(:, 7) = (/ 0.171959E-03_JPRB, 0.173431E-03_JPRB, 0.611721E-04_JPRB /)
FORREF(:, 8) = (/ 0.357795E-03_JPRB, 0.247261E-03_JPRB, 0.488864E-04_JPRB /)
FORREF(:, 9) = (/ 0.326623E-03_JPRB, 0.289471E-03_JPRB, 0.548834E-04_JPRB /)
FORREF(:,10) = (/ 0.345103E-03_JPRB, 0.320898E-03_JPRB, 0.633214E-04_JPRB /)
FORREF(:,11) = (/ 0.392567E-03_JPRB, 0.325153E-03_JPRB, 0.744479E-04_JPRB /)
FORREF(:,12) = (/ 0.349277E-03_JPRB, 0.345610E-03_JPRB, 0.916479E-04_JPRB /)
FORREF(:,13) = (/ 0.425161E-03_JPRB, 0.348452E-03_JPRB, 0.125788E-03_JPRB /)
FORREF(:,14) = (/ 0.407594E-03_JPRB, 0.435836E-03_JPRB, 0.287583E-03_JPRB /)
FORREF(:,15) = (/ 0.521605E-03_JPRB, 0.486596E-03_JPRB, 0.483511E-03_JPRB /)
FORREF(:,16) = (/ 0.773790E-03_JPRB, 0.737247E-03_JPRB, 0.665939E-03_JPRB /)

!     -----------------------------------------------------------------
!     The array SELFREF contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).
     
SELFREF(:, 1) = (/ &
 & 0.750370E-03_JPRB, 0.644938E-03_JPRB, 0.554321E-03_JPRB, 0.476436E-03_JPRB, 0.409494E-03_JPRB, &
 & 0.351957E-03_JPRB, 0.302505E-03_JPRB, 0.260002E-03_JPRB, 0.223470E-03_JPRB, 0.192071E-03_JPRB /)  
SELFREF(:, 2) = (/ &
 & 0.136135E-02_JPRB, 0.113187E-02_JPRB, 0.941076E-03_JPRB, 0.782440E-03_JPRB, 0.650546E-03_JPRB, &
 & 0.540885E-03_JPRB, 0.449709E-03_JPRB, 0.373902E-03_JPRB, 0.310874E-03_JPRB, 0.258471E-03_JPRB /)  
SELFREF(:, 3) = (/ &
 & 0.333950E-02_JPRB, 0.256391E-02_JPRB, 0.196845E-02_JPRB, 0.151129E-02_JPRB, 0.116030E-02_JPRB, &
 & 0.890824E-03_JPRB, 0.683934E-03_JPRB, 0.525093E-03_JPRB, 0.403143E-03_JPRB, 0.309515E-03_JPRB /)  
SELFREF(:, 4) = (/ &
 & 0.793392E-02_JPRB, 0.589865E-02_JPRB, 0.438548E-02_JPRB, 0.326048E-02_JPRB, 0.242408E-02_JPRB, &
 & 0.180223E-02_JPRB, 0.133991E-02_JPRB, 0.996186E-03_JPRB, 0.740636E-03_JPRB, 0.550642E-03_JPRB /)  
SELFREF(:, 5) = (/ &
 & 0.828169E-02_JPRB, 0.703139E-02_JPRB, 0.596984E-02_JPRB, 0.506856E-02_JPRB, 0.430335E-02_JPRB, &
 & 0.365366E-02_JPRB, 0.310206E-02_JPRB, 0.263374E-02_JPRB, 0.223612E-02_JPRB, 0.189852E-02_JPRB /)  
SELFREF(:, 6) = (/ &
 & 0.834190E-02_JPRB, 0.780225E-02_JPRB, 0.729750E-02_JPRB, 0.682541E-02_JPRB, 0.638386E-02_JPRB, &
 & 0.597087E-02_JPRB, 0.558460E-02_JPRB, 0.522332E-02_JPRB, 0.488541E-02_JPRB, 0.456936E-02_JPRB /)  
SELFREF(:, 7) = (/ &
 & 0.119082E-01_JPRB, 0.112566E-01_JPRB, 0.106406E-01_JPRB, 0.100583E-01_JPRB, 0.950785E-02_JPRB, &
 & 0.898755E-02_JPRB, 0.849571E-02_JPRB, 0.803080E-02_JPRB, 0.759132E-02_JPRB, 0.717590E-02_JPRB /)  
SELFREF(:, 8) = (/ &
 & 0.144004E-01_JPRB, 0.141762E-01_JPRB, 0.139554E-01_JPRB, 0.137381E-01_JPRB, 0.135241E-01_JPRB, &
 & 0.133135E-01_JPRB, 0.131062E-01_JPRB, 0.129021E-01_JPRB, 0.127011E-01_JPRB, 0.125033E-01_JPRB /)  
SELFREF(:, 9) = (/ &
 & 0.186171E-01_JPRB, 0.175281E-01_JPRB, 0.165027E-01_JPRB, 0.155373E-01_JPRB, 0.146284E-01_JPRB, &
 & 0.137726E-01_JPRB, 0.129670E-01_JPRB, 0.122084E-01_JPRB, 0.114942E-01_JPRB, 0.108218E-01_JPRB /)  
SELFREF(:,10) = (/ &
 & 0.209396E-01_JPRB, 0.195077E-01_JPRB, 0.181737E-01_JPRB, 0.169309E-01_JPRB, 0.157731E-01_JPRB, &
 & 0.146945E-01_JPRB, 0.136897E-01_JPRB, 0.127535E-01_JPRB, 0.118814E-01_JPRB, 0.110689E-01_JPRB /)  
SELFREF(:,11) = (/ &
 & 0.203661E-01_JPRB, 0.193311E-01_JPRB, 0.183487E-01_JPRB, 0.174163E-01_JPRB, 0.165312E-01_JPRB, &
 & 0.156911E-01_JPRB, 0.148937E-01_JPRB, 0.141368E-01_JPRB, 0.134184E-01_JPRB, 0.127365E-01_JPRB /)  
SELFREF(:,12) = (/ &
 & 0.226784E-01_JPRB, 0.210210E-01_JPRB, 0.194848E-01_JPRB, 0.180608E-01_JPRB, 0.167409E-01_JPRB, &
 & 0.155174E-01_JPRB, 0.143834E-01_JPRB, 0.133322E-01_JPRB, 0.123579E-01_JPRB, 0.114547E-01_JPRB /)  
SELFREF(:,13) = (/ &
 & 0.221773E-01_JPRB, 0.210306E-01_JPRB, 0.199433E-01_JPRB, 0.189122E-01_JPRB, 0.179344E-01_JPRB, &
 & 0.170071E-01_JPRB, 0.161278E-01_JPRB, 0.152939E-01_JPRB, 0.145032E-01_JPRB, 0.137533E-01_JPRB /)  
SELFREF(:,14) = (/ &
 & 0.275920E-01_JPRB, 0.252595E-01_JPRB, 0.231241E-01_JPRB, 0.211693E-01_JPRB, 0.193797E-01_JPRB, &
 & 0.177415E-01_JPRB, 0.162417E-01_JPRB, 0.148687E-01_JPRB, 0.136117E-01_JPRB, 0.124610E-01_JPRB /)  
SELFREF(:,15) = (/ &
 & 0.288687E-01_JPRB, 0.269968E-01_JPRB, 0.252462E-01_JPRB, 0.236092E-01_JPRB, 0.220783E-01_JPRB, &
 & 0.206466E-01_JPRB, 0.193078E-01_JPRB, 0.180559E-01_JPRB, 0.168851E-01_JPRB, 0.157902E-01_JPRB /)  
SELFREF(:,16) = (/ &
 & 0.371842E-01_JPRB, 0.347595E-01_JPRB, 0.324929E-01_JPRB, 0.303741E-01_JPRB, 0.283934E-01_JPRB, &
 & 0.265419E-01_JPRB, 0.248112E-01_JPRB, 0.231933E-01_JPRB, 0.216809E-01_JPRB, 0.202671E-01_JPRB /)  

IF (LHOOK) CALL DR_HOOK('SRTM_KGB18',1,ZHOOK_HANDLE)
RETURN

1001 CONTINUE
CALL ABOR1("SRTM_KGB18:ERROR READING FILE RADSRTM")

END SUBROUTINE SRTM_KGB18
