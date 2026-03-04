SUBROUTINE SRTM_KGB29

!     Originally by J.Delamere, Atmospheric & Environmental Research.
!     Revision: 2.4
!     BAND 29:   820-2600 cm-1 (low - H2O; high - CO2)
!     Reformatted for F90 by JJMorcrette, ECMWF
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
USE YOESRTA29 , ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, RAYL, &
 & ABSH2O, ABSCO2, LAYREFFR  , KA_D, KB_D

!     ------------------------------------------------------------------

IMPLICIT NONE

! KURUCZ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SRTM_KGB29',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  READ(NULRAD,ERR=1001) KA_D,KB_D
  CLOSE(NULRAD,ERR=1000)
  KA = REAL(KA_D,JPRB)
  KB = REAL(KB_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KA,MTAGRAD,1,CDSTRING='SRTM_KGB29:')
  CALL MPL_BROADCAST (KB,MTAGRAD,1,CDSTRING='SRTM_KGB29:')
ENDIF

SFLUXREF = (/ &
 & 1.32880_JPRB    , 2.14018_JPRB    , 1.97612_JPRB    , 1.79000_JPRB    , &
 & 1.51242_JPRB    , 1.22977_JPRB    , 1.06052_JPRB    , 0.800996_JPRB   , &
 & 0.748053_JPRB   , 8.64369E-02_JPRB, 7.10675E-02_JPRB, 5.62425E-02_JPRB, &
 & 4.46988E-02_JPRB, 3.07441E-02_JPRB, 1.16728E-02_JPRB, 1.65573E-03_JPRB /)  

ABSCO2 = (/ &
 & 2.90073E-06_JPRB, 2.12382E-05_JPRB, 1.03032E-04_JPRB, 1.86481E-04_JPRB, &
 & 4.31997E-04_JPRB, 6.08238E-04_JPRB, 2.17603E-03_JPRB, 4.64479E-02_JPRB, &
 & 2.96956_JPRB    , 14.9569_JPRB    , 28.4831_JPRB    , 61.3998_JPRB    , &
 & 164.129_JPRB    , 832.282_JPRB    , 4995.02_JPRB    , 12678.1_JPRB     /)  
     
ABSH2O = (/ &
 & 2.99508E-04_JPRB, 3.95012E-03_JPRB, 1.49316E-02_JPRB, 3.24384E-02_JPRB, &
 & 6.92879E-02_JPRB, 0.123523_JPRB   , 0.360985_JPRB   , 1.86434_JPRB    , &
 & 10.38157_JPRB   , 0.214129_JPRB   , 0.213914_JPRB   , 0.212781_JPRB   , &
 & 0.215562_JPRB   , 0.218087_JPRB   , 0.220918_JPRB   , 0.218546_JPRB    /)  
     
!     Rayleigh extinction coefficient at v = 2200 cm-1.
RAYL = 9.30E-11_JPRB

LAYREFFR = 49

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


FORREF(:, 1) = (/ 0.299818E-05_JPRB, 0.209282E-05_JPRB, 0.988353E-04_JPRB, 0.632178E-03_JPRB /)
FORREF(:, 2) = (/ 0.633648E-05_JPRB, 0.509214E-04_JPRB, 0.650535E-03_JPRB, 0.264019E-02_JPRB /)
FORREF(:, 3) = (/ 0.636782E-04_JPRB, 0.136577E-03_JPRB, 0.166500E-02_JPRB, 0.750821E-02_JPRB /)
FORREF(:, 4) = (/ 0.472314E-03_JPRB, 0.988296E-03_JPRB, 0.585751E-02_JPRB, 0.187352E-01_JPRB /)
FORREF(:, 5) = (/ 0.558635E-02_JPRB, 0.856489E-02_JPRB, 0.157438E-01_JPRB, 0.181471E-01_JPRB /)
FORREF(:, 6) = (/ 0.217395E-01_JPRB, 0.229156E-01_JPRB, 0.230125E-01_JPRB, 0.143821E-01_JPRB /)
FORREF(:, 7) = (/ 0.277222E-01_JPRB, 0.299252E-01_JPRB, 0.208929E-01_JPRB, 0.826748E-02_JPRB /)
FORREF(:, 8) = (/ 0.252119E-01_JPRB, 0.262911E-01_JPRB, 0.187663E-01_JPRB, 0.417110E-02_JPRB /)
FORREF(:, 9) = (/ 0.304941E-01_JPRB, 0.175545E-01_JPRB, 0.971224E-02_JPRB, 0.142023E-02_JPRB /)
FORREF(:,10) = (/ 0.327200E-01_JPRB, 0.215788E-01_JPRB, 0.346831E-02_JPRB, 0.157989E-02_JPRB /)
FORREF(:,11) = (/ 0.324955E-01_JPRB, 0.228571E-01_JPRB, 0.171749E-02_JPRB, 0.226853E-02_JPRB /)
FORREF(:,12) = (/ 0.326588E-01_JPRB, 0.198544E-01_JPRB, 0.532339E-06_JPRB, 0.279086E-02_JPRB /)
FORREF(:,13) = (/ 0.345157E-01_JPRB, 0.168679E-01_JPRB, 0.505361E-06_JPRB, 0.276647E-02_JPRB /)
FORREF(:,14) = (/ 0.448765E-01_JPRB, 0.123791E-02_JPRB, 0.488367E-06_JPRB, 0.122245E-02_JPRB /)
FORREF(:,15) = (/ 0.486925E-01_JPRB, 0.464371E-06_JPRB, 0.464241E-06_JPRB, 0.753846E-06_JPRB /)
FORREF(:,16) = (/ 0.530511E-01_JPRB, 0.376234E-06_JPRB, 0.409824E-06_JPRB, 0.470650E-06_JPRB /)

!     -----------------------------------------------------------------
!     The array SELFREF contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

SELFREF(:, 1) = (/ &
 & 0.118069E+00_JPRB, 0.713523E-01_JPRB, 0.431199E-01_JPRB, 0.260584E-01_JPRB, 0.157477E-01_JPRB, &
 & 0.951675E-02_JPRB, 0.575121E-02_JPRB, 0.347560E-02_JPRB, 0.210039E-02_JPRB, 0.126932E-02_JPRB /)  
SELFREF(:, 2) = (/ &
 & 0.137081E-01_JPRB, 0.139046E-01_JPRB, 0.141040E-01_JPRB, 0.143061E-01_JPRB, 0.145112E-01_JPRB, &
 & 0.147193E-01_JPRB, 0.149303E-01_JPRB, 0.151443E-01_JPRB, 0.153614E-01_JPRB, 0.155816E-01_JPRB /)  
SELFREF(:, 3) = (/ &
 & 0.166575E-01_JPRB, 0.164916E-01_JPRB, 0.163273E-01_JPRB, 0.161647E-01_JPRB, 0.160037E-01_JPRB, &
 & 0.158443E-01_JPRB, 0.156864E-01_JPRB, 0.155302E-01_JPRB, 0.153755E-01_JPRB, 0.152224E-01_JPRB /)  
SELFREF(:, 4) = (/ &
 & 0.597379E-01_JPRB, 0.509517E-01_JPRB, 0.434579E-01_JPRB, 0.370662E-01_JPRB, 0.316145E-01_JPRB, &
 & 0.269647E-01_JPRB, 0.229988E-01_JPRB, 0.196162E-01_JPRB, 0.167311E-01_JPRB, 0.142703E-01_JPRB /)  
SELFREF(:, 5) = (/ &
 & 0.227517E+00_JPRB, 0.198401E+00_JPRB, 0.173011E+00_JPRB, 0.150870E+00_JPRB, 0.131563E+00_JPRB, &
 & 0.114726E+00_JPRB, 0.100044E+00_JPRB, 0.872415E-01_JPRB, 0.760769E-01_JPRB, 0.663411E-01_JPRB /)  
SELFREF(:, 6) = (/ &
 & 0.453235E+00_JPRB, 0.414848E+00_JPRB, 0.379712E+00_JPRB, 0.347552E+00_JPRB, 0.318116E+00_JPRB, &
 & 0.291173E+00_JPRB, 0.266512E+00_JPRB, 0.243940E+00_JPRB, 0.223279E+00_JPRB, 0.204368E+00_JPRB /)  
SELFREF(:, 7) = (/ &
 & 0.569263E+00_JPRB, 0.516415E+00_JPRB, 0.468473E+00_JPRB, 0.424982E+00_JPRB, 0.385528E+00_JPRB, &
 & 0.349737E+00_JPRB, 0.317269E+00_JPRB, 0.287815E+00_JPRB, 0.261095E+00_JPRB, 0.236856E+00_JPRB /)  
SELFREF(:, 8) = (/ &
 & 0.490314E+00_JPRB, 0.448042E+00_JPRB, 0.409413E+00_JPRB, 0.374116E+00_JPRB, 0.341861E+00_JPRB, &
 & 0.312387E+00_JPRB, 0.285455E+00_JPRB, 0.260844E+00_JPRB, 0.238355E+00_JPRB, 0.217805E+00_JPRB /)  
SELFREF(:, 9) = (/ &
 & 0.258162E+00_JPRB, 0.265085E+00_JPRB, 0.272193E+00_JPRB, 0.279493E+00_JPRB, 0.286988E+00_JPRB, &
 & 0.294684E+00_JPRB, 0.302586E+00_JPRB, 0.310701E+00_JPRB, 0.319033E+00_JPRB, 0.327588E+00_JPRB /)  
SELFREF(:,10) = (/ &
 & 0.332019E+00_JPRB, 0.331902E+00_JPRB, 0.331784E+00_JPRB, 0.331666E+00_JPRB, 0.331549E+00_JPRB, &
 & 0.331431E+00_JPRB, 0.331314E+00_JPRB, 0.331197E+00_JPRB, 0.331079E+00_JPRB, 0.330962E+00_JPRB /)  
SELFREF(:,11) = (/ &
 & 0.357523E+00_JPRB, 0.353154E+00_JPRB, 0.348839E+00_JPRB, 0.344576E+00_JPRB, 0.340366E+00_JPRB, &
 & 0.336207E+00_JPRB, 0.332099E+00_JPRB, 0.328041E+00_JPRB, 0.324032E+00_JPRB, 0.320073E+00_JPRB /)  
SELFREF(:,12) = (/ &
 & 0.294662E+00_JPRB, 0.299043E+00_JPRB, 0.303488E+00_JPRB, 0.308000E+00_JPRB, 0.312579E+00_JPRB, &
 & 0.317226E+00_JPRB, 0.321941E+00_JPRB, 0.326727E+00_JPRB, 0.331585E+00_JPRB, 0.336514E+00_JPRB /)  
SELFREF(:,13) = (/ &
 & 0.227445E+00_JPRB, 0.241545E+00_JPRB, 0.256519E+00_JPRB, 0.272422E+00_JPRB, 0.289311E+00_JPRB, &
 & 0.307247E+00_JPRB, 0.326294E+00_JPRB, 0.346523E+00_JPRB, 0.368005E+00_JPRB, 0.390820E+00_JPRB /)  
SELFREF(:,14) = (/ &
 & 0.616203E-02_JPRB, 0.113523E-01_JPRB, 0.209144E-01_JPRB, 0.385307E-01_JPRB, 0.709852E-01_JPRB, &
 & 0.130776E+00_JPRB, 0.240929E+00_JPRB, 0.443865E+00_JPRB, 0.817733E+00_JPRB, 0.150651E+01_JPRB /)  
SELFREF(:,15) = (/ &
 & 0.279552E-03_JPRB, 0.808472E-03_JPRB, 0.233812E-02_JPRB, 0.676192E-02_JPRB, 0.195557E-01_JPRB, &
 & 0.565555E-01_JPRB, 0.163560E+00_JPRB, 0.473020E+00_JPRB, 0.136799E+01_JPRB, 0.395626E+01_JPRB /)  
SELFREF(:,16) = (/ &
 & 0.261006E-03_JPRB, 0.771043E-03_JPRB, 0.227776E-02_JPRB, 0.672879E-02_JPRB, 0.198777E-01_JPRB, &
 & 0.587212E-01_JPRB, 0.173470E+00_JPRB, 0.512452E+00_JPRB, 0.151385E+01_JPRB, 0.447209E+01_JPRB /)  
     
!     -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_KGB29',1,ZHOOK_HANDLE)
RETURN

1000 CONTINUE
CALL ABOR1("SRTM_KGB29:ERROR CLOSING FILE RADSRTM")
1001 CONTINUE
CALL ABOR1("SRTM_KGB29:ERROR READING FILE RADSRTM")

END SUBROUTINE SRTM_KGB29
