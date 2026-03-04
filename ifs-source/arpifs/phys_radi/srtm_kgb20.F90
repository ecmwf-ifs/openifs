SUBROUTINE SRTM_KGB20

!     Originally by J.Delamere, Atmospheric & Environmental Research.
!     Revision: 2.4
!     BAND 20:  5150-6150 cm-1 (low - H2O; high - H2O)
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
USE YOESRTA20 , ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, RAYL, ABSCH4, LAYREFFR  ,&
  & KA_D, KB_D

!     ------------------------------------------------------------------

IMPLICIT NONE

! KURUCZ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SRTM_KGB20',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  READ(NULRAD,ERR=1001) KA_D,KB_D
  KA = REAL(KA_D,JPRB)
  KB = REAL(KB_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KA,MTAGRAD,1,CDSTRING='SRTM_KGB20:')
  CALL MPL_BROADCAST (KB,MTAGRAD,1,CDSTRING='SRTM_KGB20:')
ENDIF

SFLUXREF = (/ &
 & 9.34081_JPRB , 8.93720_JPRB    , 8.19346_JPRB    , 7.39196_JPRB    , &
 & 6.12127_JPRB , 5.23956_JPRB    , 4.24941_JPRB    , 3.20013_JPRB    , &
 & 2.16047_JPRB , 0.234509_JPRB   , 0.194593_JPRB   , 0.151512_JPRB   , &
 & 0.110315_JPRB, 7.09959E-02_JPRB, 2.70573E-02_JPRB, 3.36042E-03_JPRB /)  
  
ABSCH4 = (/   &
 & 1.01381E-03_JPRB,6.33692E-03_JPRB,1.94185E-02_JPRB,4.83210E-02_JPRB, &
 & 2.36574E-03_JPRB,6.61973E-04_JPRB,5.64552E-04_JPRB,2.83183E-04_JPRB, &
 & 7.43623E-05_JPRB,8.90159E-07_JPRB,6.98728E-07_JPRB,6.51832E-08_JPRB, &
 & 2.96619E-08_JPRB,         0._JPRB,         0._JPRB,         0._JPRB /)  

!     Rayleigh extinction coefficient at v = 5670 cm-1.
RAYL = 4.12E-09_JPRB

LAYREFFR = 3

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


FORREF(:, 1) = (/ 0.214504E-06_JPRB, 0.460418E-06_JPRB, 0.357608E-05_JPRB, 0.192037E-05_JPRB /)
FORREF(:, 2) = (/ 0.142576E-05_JPRB, 0.364463E-05_JPRB, 0.117033E-04_JPRB, 0.112085E-04_JPRB /)
FORREF(:, 3) = (/ 0.101536E-04_JPRB, 0.124096E-04_JPRB, 0.509190E-04_JPRB, 0.565282E-04_JPRB /)
FORREF(:, 4) = (/ 0.143394E-03_JPRB, 0.154700E-03_JPRB, 0.466498E-03_JPRB, 0.918829E-03_JPRB /)
FORREF(:, 5) = (/ 0.251631E-02_JPRB, 0.241729E-02_JPRB, 0.240057E-02_JPRB, 0.350408E-02_JPRB /)
FORREF(:, 6) = (/ 0.410309E-02_JPRB, 0.416851E-02_JPRB, 0.390925E-02_JPRB, 0.383694E-02_JPRB /)
FORREF(:, 7) = (/ 0.445387E-02_JPRB, 0.448657E-02_JPRB, 0.432310E-02_JPRB, 0.370739E-02_JPRB /)
FORREF(:, 8) = (/ 0.458150E-02_JPRB, 0.460014E-02_JPRB, 0.450245E-02_JPRB, 0.336718E-02_JPRB /)
FORREF(:, 9) = (/ 0.465423E-02_JPRB, 0.465595E-02_JPRB, 0.467006E-02_JPRB, 0.368061E-02_JPRB /)
FORREF(:,10) = (/ 0.493955E-02_JPRB, 0.490181E-02_JPRB, 0.481941E-02_JPRB, 0.367577E-02_JPRB /)
FORREF(:,11) = (/ 0.511876E-02_JPRB, 0.490981E-02_JPRB, 0.493303E-02_JPRB, 0.357423E-02_JPRB /)
FORREF(:,12) = (/ 0.509845E-02_JPRB, 0.511556E-02_JPRB, 0.504031E-02_JPRB, 0.355915E-02_JPRB /)
FORREF(:,13) = (/ 0.523822E-02_JPRB, 0.530473E-02_JPRB, 0.523811E-02_JPRB, 0.414259E-02_JPRB /)
FORREF(:,14) = (/ 0.551133E-02_JPRB, 0.535831E-02_JPRB, 0.546702E-02_JPRB, 0.473875E-02_JPRB /)
FORREF(:,15) = (/ 0.609781E-02_JPRB, 0.589859E-02_JPRB, 0.561187E-02_JPRB, 0.528981E-02_JPRB /)
FORREF(:,16) = (/ 0.644958E-02_JPRB, 0.631718E-02_JPRB, 0.625201E-02_JPRB, 0.600448E-02_JPRB /)

!     -----------------------------------------------------------------
!     The array SELFREF contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

SELFREF(:, 1) = (/ &
 & 0.217058E-03_JPRB, 0.176391E-03_JPRB, 0.143342E-03_JPRB, 0.116486E-03_JPRB, 0.946614E-04_JPRB, &
 & 0.769257E-04_JPRB, 0.625131E-04_JPRB, 0.508007E-04_JPRB, 0.412828E-04_JPRB, 0.335481E-04_JPRB /)  
SELFREF(:, 2) = (/ &
 & 0.598055E-03_JPRB, 0.484805E-03_JPRB, 0.393000E-03_JPRB, 0.318580E-03_JPRB, 0.258252E-03_JPRB, &
 & 0.209348E-03_JPRB, 0.169705E-03_JPRB, 0.137569E-03_JPRB, 0.111518E-03_JPRB, 0.904008E-04_JPRB /)  
SELFREF(:, 3) = (/ &
 & 0.102691E-02_JPRB, 0.930281E-03_JPRB, 0.842740E-03_JPRB, 0.763437E-03_JPRB, 0.691596E-03_JPRB, &
 & 0.626516E-03_JPRB, 0.567560E-03_JPRB, 0.514152E-03_JPRB, 0.465769E-03_JPRB, 0.421940E-03_JPRB /)  
SELFREF(:, 4) = (/ &
 & 0.388569E-02_JPRB, 0.365098E-02_JPRB, 0.343045E-02_JPRB, 0.322324E-02_JPRB, 0.302854E-02_JPRB, &
 & 0.284561E-02_JPRB, 0.267372E-02_JPRB, 0.251222E-02_JPRB, 0.236047E-02_JPRB, 0.221789E-02_JPRB /)  
SELFREF(:, 5) = (/ &
 & 0.349845E-01_JPRB, 0.326678E-01_JPRB, 0.305045E-01_JPRB, 0.284845E-01_JPRB, 0.265982E-01_JPRB, &
 & 0.248369E-01_JPRB, 0.231921E-01_JPRB, 0.216563E-01_JPRB, 0.202222E-01_JPRB, 0.188831E-01_JPRB /)  
SELFREF(:, 6) = (/ &
 & 0.613705E-01_JPRB, 0.562676E-01_JPRB, 0.515890E-01_JPRB, 0.472994E-01_JPRB, 0.433665E-01_JPRB, &
 & 0.397606E-01_JPRB, 0.364545E-01_JPRB, 0.334233E-01_JPRB, 0.306442E-01_JPRB, 0.280961E-01_JPRB /)  
SELFREF(:, 7) = (/ &
 & 0.656981E-01_JPRB, 0.602660E-01_JPRB, 0.552830E-01_JPRB, 0.507120E-01_JPRB, 0.465190E-01_JPRB, &
 & 0.426726E-01_JPRB, 0.391443E-01_JPRB, 0.359077E-01_JPRB, 0.329387E-01_JPRB, 0.302153E-01_JPRB /)  
SELFREF(:, 8) = (/ &
 & 0.671782E-01_JPRB, 0.616461E-01_JPRB, 0.565695E-01_JPRB, 0.519110E-01_JPRB, 0.476361E-01_JPRB, &
 & 0.437132E-01_JPRB, 0.401134E-01_JPRB, 0.368100E-01_JPRB, 0.337787E-01_JPRB, 0.309970E-01_JPRB /)  
SELFREF(:, 9) = (/ &
 & 0.675902E-01_JPRB, 0.620888E-01_JPRB, 0.570351E-01_JPRB, 0.523928E-01_JPRB, 0.481284E-01_JPRB, &
 & 0.442110E-01_JPRB, 0.406125E-01_JPRB, 0.373069E-01_JPRB, 0.342703E-01_JPRB, 0.314809E-01_JPRB /)  
SELFREF(:,10) = (/ &
 & 0.708308E-01_JPRB, 0.651419E-01_JPRB, 0.599099E-01_JPRB, 0.550981E-01_JPRB, 0.506728E-01_JPRB, &
 & 0.466030E-01_JPRB, 0.428600E-01_JPRB, 0.394176E-01_JPRB, 0.362517E-01_JPRB, 0.333401E-01_JPRB /)  
SELFREF(:,11) = (/ &
 & 0.698445E-01_JPRB, 0.646584E-01_JPRB, 0.598573E-01_JPRB, 0.554128E-01_JPRB, 0.512982E-01_JPRB, &
 & 0.474892E-01_JPRB, 0.439630E-01_JPRB, 0.406986E-01_JPRB, 0.376766E-01_JPRB, 0.348791E-01_JPRB /)  
SELFREF(:,12) = (/ &
 & 0.743921E-01_JPRB, 0.682057E-01_JPRB, 0.625337E-01_JPRB, 0.573334E-01_JPRB, 0.525655E-01_JPRB, &
 & 0.481942E-01_JPRB, 0.441863E-01_JPRB, 0.405118E-01_JPRB, 0.371428E-01_JPRB, 0.340540E-01_JPRB /)  
SELFREF(:,13) = (/ &
 & 0.775758E-01_JPRB, 0.709818E-01_JPRB, 0.649484E-01_JPRB, 0.594277E-01_JPRB, 0.543764E-01_JPRB, &
 & 0.497544E-01_JPRB, 0.455253E-01_JPRB, 0.416556E-01_JPRB, 0.381149E-01_JPRB, 0.348751E-01_JPRB /)  
SELFREF(:,14) = (/ &
 & 0.776545E-01_JPRB, 0.714761E-01_JPRB, 0.657894E-01_JPRB, 0.605550E-01_JPRB, 0.557372E-01_JPRB, &
 & 0.513026E-01_JPRB, 0.472209E-01_JPRB, 0.434639E-01_JPRB, 0.400058E-01_JPRB, 0.368229E-01_JPRB /)  
SELFREF(:,15) = (/ &
 & 0.855675E-01_JPRB, 0.787337E-01_JPRB, 0.724456E-01_JPRB, 0.666598E-01_JPRB, 0.613360E-01_JPRB, &
 & 0.564374E-01_JPRB, 0.519301E-01_JPRB, 0.477827E-01_JPRB, 0.439666E-01_JPRB, 0.404552E-01_JPRB /)  
SELFREF(:,16) = (/ &
 & 0.934781E-01_JPRB, 0.855190E-01_JPRB, 0.782376E-01_JPRB, 0.715761E-01_JPRB, 0.654819E-01_JPRB, &
 & 0.599065E-01_JPRB, 0.548058E-01_JPRB, 0.501394E-01_JPRB, 0.458704E-01_JPRB, 0.419648E-01_JPRB /)  
     
!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SRTM_KGB20',1,ZHOOK_HANDLE)
RETURN

1001 CONTINUE
CALL ABOR1("SRTM_KGB20:ERROR READING FILE RADSRTM")

END SUBROUTINE SRTM_KGB20
