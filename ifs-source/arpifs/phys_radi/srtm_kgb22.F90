SUBROUTINE SRTM_KGB22

!     Originally by J.Delamere, Atmospheric & Environmental Research.
!     Revision: 2.4
!     BAND 16:  7700-8050 cm-1 (low - H2O,O2; high - O2)
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
USE YOESRTA22 , ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, RAYL, STRRAT, LAYREFFR  ,&
  &  KA_D, KB_D

!     ------------------------------------------------------------------

IMPLICIT NONE

! KURUCZ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SRTM_KGB22',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  READ(NULRAD,ERR=1001) KA_D,KB_D
  KA = REAL(KA_D,JPRB)
  KB = REAL(KB_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KA,MTAGRAD,1,CDSTRING='SRTM_KGB22:')
  CALL MPL_BROADCAST (KB,MTAGRAD,1,CDSTRING='SRTM_KGB22:')
ENDIF

SFLUXREF(:, 1) = (/ &
 & 3.71641_JPRB    ,3.63190_JPRB    ,3.44795_JPRB    ,3.17936_JPRB    , &
 & 2.86071_JPRB    ,2.48490_JPRB    ,2.02471_JPRB    ,1.52475_JPRB    , &
 & 1.03811_JPRB    ,0.113272_JPRB   ,9.37115E-02_JPRB,7.38969E-02_JPRB, &
 & 5.44713E-02_JPRB,3.45905E-02_JPRB,1.30293E-02_JPRB,1.84198E-03_JPRB /)  
SFLUXREF(:, 2) = (/ &
 & 3.73933_JPRB    ,3.60360_JPRB    ,3.43370_JPRB    ,3.19749_JPRB    ,  &
 & 2.87747_JPRB    ,2.47926_JPRB    ,2.02175_JPRB    ,1.52010_JPRB    , &
 & 1.03612_JPRB    ,0.113265_JPRB   ,9.37145E-02_JPRB,7.38951E-02_JPRB, &
 & 5.44714E-02_JPRB,3.45906E-02_JPRB,1.30293E-02_JPRB,1.84198E-03_JPRB /)  
SFLUXREF(:, 3) = (/ &
 & 3.73889_JPRB    ,3.60279_JPRB    ,3.43404_JPRB    ,3.20560_JPRB    , &
 & 2.87367_JPRB    ,2.47515_JPRB    ,2.02412_JPRB    ,1.52315_JPRB    , &
 & 1.03146_JPRB    ,0.113272_JPRB   ,9.36707E-02_JPRB,7.39080E-02_JPRB, &
 & 5.44598E-02_JPRB,3.45906E-02_JPRB,1.30293E-02_JPRB,1.84198E-03_JPRB /)  
SFLUXREF(:, 4) = (/ &
 & 3.73801_JPRB    ,3.60530_JPRB    ,3.43659_JPRB    ,3.20640_JPRB    , &
 & 2.87039_JPRB    ,2.47330_JPRB    ,2.02428_JPRB    ,1.52509_JPRB    , &
 & 1.03037_JPRB    ,0.112553_JPRB   ,9.35352E-02_JPRB,7.39675E-02_JPRB, &
 & 5.43951E-02_JPRB,3.45669E-02_JPRB,1.30292E-02_JPRB,1.84198E-03_JPRB /)  
SFLUXREF(:, 5) = (/ &
 & 3.73809_JPRB    ,3.60996_JPRB    ,3.43602_JPRB    ,3.20364_JPRB    , &
 & 2.87005_JPRB    ,2.47343_JPRB    ,2.02353_JPRB    ,1.52617_JPRB    , &
 & 1.03138_JPRB    ,0.111172_JPRB   ,9.29885E-02_JPRB,7.35034E-02_JPRB, &
 & 5.42427E-02_JPRB,3.45732E-02_JPRB,1.30169E-02_JPRB,1.84550E-03_JPRB /)  
SFLUXREF(:, 6) = (/ &
 & 3.73872_JPRB    ,3.62054_JPRB    ,3.42934_JPRB    ,3.20110_JPRB    , &
 & 2.86886_JPRB    ,2.47379_JPRB    ,2.02237_JPRB    ,1.52754_JPRB    ,  &
 & 1.03228_JPRB    ,0.111597_JPRB   ,9.12252E-02_JPRB,7.33115E-02_JPRB, &
 & 5.35600E-02_JPRB,3.45187E-02_JPRB,1.30184E-02_JPRB,1.84551E-03_JPRB /)  
SFLUXREF(:, 7) = (/ &
 & 3.73969_JPRB    ,3.65461_JPRB    ,3.40646_JPRB    ,3.19082_JPRB    , &
 & 2.86919_JPRB    ,2.47289_JPRB    ,2.02312_JPRB    ,1.52629_JPRB    , &
 & 1.03329_JPRB    ,0.111611_JPRB   ,9.16275E-02_JPRB,7.14731E-02_JPRB, &
 & 5.31771E-02_JPRB,3.44980E-02_JPRB,1.30190E-02_JPRB,1.84551E-03_JPRB /)  
SFLUXREF(:, 8) = (/ &
 & 3.73995_JPRB    ,3.65348_JPRB    ,3.43707_JPRB    ,3.16351_JPRB    , &
 & 2.87003_JPRB    ,2.47392_JPRB    ,2.02114_JPRB    ,1.52548_JPRB    ,  &
 & 1.03306_JPRB    ,0.111088_JPRB   ,9.12422E-02_JPRB,7.11146E-02_JPRB, &
 & 5.31333E-02_JPRB,3.45302E-02_JPRB,1.30209E-02_JPRB,1.84554E-03_JPRB /)  
SFLUXREF(:, 9) = (/ &
 & 3.73788_JPRB    ,3.65004_JPRB    ,3.46938_JPRB    ,3.15236_JPRB    , &
 & 2.86381_JPRB    ,2.47393_JPRB    ,2.01715_JPRB    ,1.52134_JPRB    , &
 & 1.03163_JPRB    ,0.111259_JPRB   ,9.12948E-02_JPRB,7.09999E-02_JPRB, &
 & 5.31792E-02_JPRB,3.44955E-02_JPRB,1.30189E-02_JPRB,1.84551E-03_JPRB /)  

!     Rayleigh extinction coefficient at v = 8000 cm-1.
RAYL = 1.54E-08_JPRB

STRRAT = 0.022708_JPRB

LAYREFFR = 2

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


FORREF(:, 1) = (/ 0.351362E-07_JPRB, 0.341136E-07_JPRB, 0.181317E-06_JPRB /)
FORREF(:, 2) = (/ 0.109648E-06_JPRB, 0.344240E-06_JPRB, 0.139709E-05_JPRB /)
FORREF(:, 3) = (/ 0.374823E-06_JPRB, 0.103424E-05_JPRB, 0.188717E-05_JPRB /)
FORREF(:, 4) = (/ 0.580041E-06_JPRB, 0.116876E-05_JPRB, 0.121183E-05_JPRB /)
FORREF(:, 5) = (/ 0.115608E-05_JPRB, 0.148110E-05_JPRB, 0.836083E-06_JPRB /)
FORREF(:, 6) = (/ 0.181460E-05_JPRB, 0.133313E-05_JPRB, 0.500167E-06_JPRB /)
FORREF(:, 7) = (/ 0.199096E-05_JPRB, 0.115276E-05_JPRB, 0.432994E-06_JPRB /)
FORREF(:, 8) = (/ 0.183730E-05_JPRB, 0.122260E-05_JPRB, 0.433248E-06_JPRB /)
FORREF(:, 9) = (/ 0.198386E-05_JPRB, 0.100130E-05_JPRB, 0.269712E-06_JPRB /)
FORREF(:,10) = (/ 0.276382E-05_JPRB, 0.749215E-06_JPRB, 0.236919E-06_JPRB /)
FORREF(:,11) = (/ 0.298202E-05_JPRB, 0.629688E-06_JPRB, 0.228388E-06_JPRB /)
FORREF(:,12) = (/ 0.364604E-05_JPRB, 0.455336E-06_JPRB, 0.206130E-06_JPRB /)
FORREF(:,13) = (/ 0.373339E-05_JPRB, 0.245210E-06_JPRB, 0.201987E-06_JPRB /)
FORREF(:,14) = (/ 0.480378E-05_JPRB, 0.177591E-06_JPRB, 0.171458E-06_JPRB /)
FORREF(:,15) = (/ 0.521700E-05_JPRB, 0.203358E-06_JPRB, 0.189559E-06_JPRB /)
FORREF(:,16) = (/ 0.542717E-05_JPRB, 0.219022E-06_JPRB, 0.218271E-06_JPRB /)

!     -----------------------------------------------------------------
!     The array SELFREF contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

SELFREF(:, 1) = (/ &
 & 0.538526E-04_JPRB, 0.464603E-04_JPRB, 0.400828E-04_JPRB, 0.345807E-04_JPRB, 0.298339E-04_JPRB, &
 & 0.257386E-04_JPRB, 0.222055E-04_JPRB, 0.191574E-04_JPRB, 0.165277E-04_JPRB, 0.142590E-04_JPRB /)  
SELFREF(:, 2) = (/ &
 & 0.162409E-03_JPRB, 0.128347E-03_JPRB, 0.101430E-03_JPRB, 0.801571E-04_JPRB, 0.633460E-04_JPRB, &
 & 0.500607E-04_JPRB, 0.395616E-04_JPRB, 0.312645E-04_JPRB, 0.247075E-04_JPRB, 0.195257E-04_JPRB /)  
SELFREF(:, 3) = (/ &
 & 0.262882E-03_JPRB, 0.212793E-03_JPRB, 0.172247E-03_JPRB, 0.139427E-03_JPRB, 0.112860E-03_JPRB, &
 & 0.913557E-04_JPRB, 0.739487E-04_JPRB, 0.598584E-04_JPRB, 0.484529E-04_JPRB, 0.392206E-04_JPRB /)  
SELFREF(:, 4) = (/ &
 & 0.242873E-03_JPRB, 0.204225E-03_JPRB, 0.171726E-03_JPRB, 0.144399E-03_JPRB, 0.121421E-03_JPRB, &
 & 0.102099E-03_JPRB, 0.858516E-04_JPRB, 0.721899E-04_JPRB, 0.607022E-04_JPRB, 0.510426E-04_JPRB /)  
SELFREF(:, 5) = (/ &
 & 0.235614E-03_JPRB, 0.207814E-03_JPRB, 0.183293E-03_JPRB, 0.161666E-03_JPRB, 0.142591E-03_JPRB, &
 & 0.125766E-03_JPRB, 0.110927E-03_JPRB, 0.978381E-04_JPRB, 0.862939E-04_JPRB, 0.761119E-04_JPRB /)  
SELFREF(:, 6) = (/ &
 & 0.205508E-03_JPRB, 0.190174E-03_JPRB, 0.175985E-03_JPRB, 0.162854E-03_JPRB, 0.150702E-03_JPRB, &
 & 0.139458E-03_JPRB, 0.129052E-03_JPRB, 0.119423E-03_JPRB, 0.110513E-03_JPRB, 0.102267E-03_JPRB /)  
SELFREF(:, 7) = (/ &
 & 0.185027E-03_JPRB, 0.175148E-03_JPRB, 0.165796E-03_JPRB, 0.156944E-03_JPRB, 0.148565E-03_JPRB, &
 & 0.140633E-03_JPRB, 0.133124E-03_JPRB, 0.126016E-03_JPRB, 0.119288E-03_JPRB, 0.112919E-03_JPRB /)  
SELFREF(:, 8) = (/ &
 & 0.192634E-03_JPRB, 0.180192E-03_JPRB, 0.168554E-03_JPRB, 0.157668E-03_JPRB, 0.147484E-03_JPRB, &
 & 0.137959E-03_JPRB, 0.129048E-03_JPRB, 0.120713E-03_JPRB, 0.112917E-03_JPRB, 0.105624E-03_JPRB /)  
SELFREF(:, 9) = (/ &
 & 0.161632E-03_JPRB, 0.155919E-03_JPRB, 0.150408E-03_JPRB, 0.145092E-03_JPRB, 0.139963E-03_JPRB, &
 & 0.135016E-03_JPRB, 0.130244E-03_JPRB, 0.125640E-03_JPRB, 0.121199E-03_JPRB, 0.116915E-03_JPRB /)  
SELFREF(:,10) = (/ &
 & 0.120880E-03_JPRB, 0.125265E-03_JPRB, 0.129810E-03_JPRB, 0.134520E-03_JPRB, 0.139400E-03_JPRB, &
 & 0.144458E-03_JPRB, 0.149699E-03_JPRB, 0.155130E-03_JPRB, 0.160758E-03_JPRB, 0.166591E-03_JPRB /)  
SELFREF(:,11) = (/ &
 & 0.104705E-03_JPRB, 0.111761E-03_JPRB, 0.119291E-03_JPRB, 0.127330E-03_JPRB, 0.135910E-03_JPRB, &
 & 0.145068E-03_JPRB, 0.154843E-03_JPRB, 0.165277E-03_JPRB, 0.176414E-03_JPRB, 0.188302E-03_JPRB /)  
SELFREF(:,12) = (/ &
 & 0.846335E-04_JPRB, 0.951236E-04_JPRB, 0.106914E-03_JPRB, 0.120166E-03_JPRB, 0.135060E-03_JPRB, &
 & 0.151800E-03_JPRB, 0.170616E-03_JPRB, 0.191763E-03_JPRB, 0.215532E-03_JPRB, 0.242246E-03_JPRB /)  
SELFREF(:,13) = (/ &
 & 0.669754E-04_JPRB, 0.781902E-04_JPRB, 0.912829E-04_JPRB, 0.106568E-03_JPRB, 0.124413E-03_JPRB, &
 & 0.145245E-03_JPRB, 0.169566E-03_JPRB, 0.197959E-03_JPRB, 0.231107E-03_JPRB, 0.269805E-03_JPRB /)  
SELFREF(:,14) = (/ &
 & 0.597091E-04_JPRB, 0.722265E-04_JPRB, 0.873679E-04_JPRB, 0.105684E-03_JPRB, 0.127839E-03_JPRB, &
 & 0.154639E-03_JPRB, 0.187057E-03_JPRB, 0.226272E-03_JPRB, 0.273707E-03_JPRB, 0.331087E-03_JPRB /)  
SELFREF(:,15) = (/ &
 & 0.640410E-04_JPRB, 0.771879E-04_JPRB, 0.930338E-04_JPRB, 0.112133E-03_JPRB, 0.135152E-03_JPRB, &
 & 0.162897E-03_JPRB, 0.196338E-03_JPRB, 0.236644E-03_JPRB, 0.285225E-03_JPRB, 0.343778E-03_JPRB /)  
SELFREF(:,16) = (/ &
 & 0.666420E-04_JPRB, 0.801056E-04_JPRB, 0.962892E-04_JPRB, 0.115742E-03_JPRB, 0.139126E-03_JPRB, &
 & 0.167233E-03_JPRB, 0.201019E-03_JPRB, 0.241630E-03_JPRB, 0.290446E-03_JPRB, 0.349125E-03_JPRB /)  

IF (LHOOK) CALL DR_HOOK('SRTM_KGB22',1,ZHOOK_HANDLE)
RETURN

1001 CONTINUE
CALL ABOR1("SRTM_KGB22:ERROR READING FILE RADSRTM")

END SUBROUTINE SRTM_KGB22
