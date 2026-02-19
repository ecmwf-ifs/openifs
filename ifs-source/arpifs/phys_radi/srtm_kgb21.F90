SUBROUTINE SRTM_KGB21

!     Originally by J.Delamere, Atmospheric & Environmental Research.
!     Revision: 2.4
!     BAND 21:  6150-7700 cm-1 (low - H2O,CO2; high - H2O,CO2)
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
USE YOESRTA21 , ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, RAYL, STRRAT, LAYREFFR  ,&
  &  KA_D, KB_D

!     ------------------------------------------------------------------

IMPLICIT NONE

! KURUCZ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SRTM_KGB21',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  READ(NULRAD,ERR=1001) KA_D,KB_D
  KA = REAL(KA_D,JPRB)
  KB = REAL(KB_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KA,MTAGRAD,1,CDSTRING='SRTM_KGB21:')
  CALL MPL_BROADCAST (KB,MTAGRAD,1,CDSTRING='SRTM_KGB21:')
ENDIF

SFLUXREF(:, 1) = (/ &
 & 16.1643_JPRB , 15.5806_JPRB, 14.7254_JPRB    , 13.5541_JPRB    , &
 & 11.9519_JPRB ,10.44410_JPRB, 8.37884_JPRB    , 6.26384_JPRB    , &
 & 4.28435_JPRB ,0.465228_JPRB, 0.385095_JPRB   ,0.304226_JPRB    , &
 & 0.222479_JPRB,0.143286_JPRB, 5.58046E-02_JPRB, 7.84856E-03_JPRB /)  
SFLUXREF(:, 2) = (/ &
 & 15.6451_JPRB , 15.3170_JPRB, 14.6987_JPRB    , 13.7350_JPRB    , &
 & 12.2267_JPRB ,10.51646_JPRB, 8.47150_JPRB    , 6.38873_JPRB    , &
 & 4.33536_JPRB ,0.470610_JPRB,0.389426_JPRB    ,0.306461_JPRB    , &
 & 0.223537_JPRB,0.143273_JPRB, 5.58179E-02_JPRB, 7.84856E-03_JPRB /)  
SFLUXREF(:, 3) = (/ &
 & 15.6092_JPRB , 15.3293_JPRB, 14.6881_JPRB    , 13.6693_JPRB    , &
 & 12.2342_JPRB ,10.52010_JPRB, 8.49442_JPRB    , 6.42138_JPRB    , &
 & 4.35865_JPRB ,0.473349_JPRB,0.391349_JPRB    ,0.308861_JPRB    , &
 & 0.224666_JPRB,0.144799_JPRB, 5.58176E-02_JPRB, 7.84881E-03_JPRB /)  
SFLUXREF(:, 4) = (/ &
 & 15.5786_JPRB , 15.3422_JPRB, 14.6894_JPRB    , 13.6040_JPRB    , &
 & 12.2567_JPRB ,10.49400_JPRB, 8.53521_JPRB    , 6.44427_JPRB    , &
 & 4.37208_JPRB ,0.475709_JPRB,0.392956_JPRB    ,0.309737_JPRB    , &
 & 0.226274_JPRB,0.146483_JPRB, 5.59325E-02_JPRB, 7.84881E-03_JPRB /)  
SFLUXREF(:, 5) = (/ &
 & 15.5380_JPRB , 15.3826_JPRB, 14.6575_JPRB    , 13.5722_JPRB    , &
 & 12.2646_JPRB ,10.47672_JPRB, 8.57158_JPRB    , 6.46343_JPRB    , &
 & 4.38259_JPRB ,0.477647_JPRB,0.393982_JPRB    ,0.310686_JPRB    , &
 & 0.227620_JPRB,0.148376_JPRB, 5.60398E-02_JPRB, 7.83925E-03_JPRB /)  
SFLUXREF(:, 6) = (/ &
 & 15.5124_JPRB , 15.3986_JPRB, 14.6240_JPRB    , 13.5535_JPRB    , &
 & 12.2468_JPRB ,10.48891_JPRB, 8.60434_JPRB    , 6.47985_JPRB    , &
 & 4.39448_JPRB ,0.478267_JPRB,0.395618_JPRB    ,0.311043_JPRB    , &
 & 0.230927_JPRB,0.148774_JPRB, 5.61189E-02_JPRB, 7.83925E-03_JPRB /)  
SFLUXREF(:, 7) = (/ &
 & 15.4910_JPRB , 15.4028_JPRB, 14.5772_JPRB    , 13.5507_JPRB    , &
 & 12.2122_JPRB ,10.52735_JPRB, 8.62650_JPRB    , 6.49644_JPRB    , &
 & 4.41173_JPRB ,0.478627_JPRB,0.396433_JPRB    ,0.314199_JPRB    ,  &
 & 0.233125_JPRB,0.149052_JPRB, 5.62309E-02_JPRB, 7.83925E-03_JPRB /)  
SFLUXREF(:, 8) = (/ &
 & 15.4562_JPRB , 15.3928_JPRB, 14.5510_JPRB    , 13.5122_JPRB    , &
 & 12.1890_JPRB , 10.5826_JPRB, 8.65842_JPRB    , 6.51558_JPRB    , &
 & 4.42747_JPRB ,0.480669_JPRB,0.400143_JPRB    ,0.318144_JPRB    , &
 & 0.233937_JPRB,0.149119_JPRB, 5.62309E-02_JPRB, 7.83925E-03_JPRB /)  
SFLUXREF(:, 9) = (/ &
 & 15.0069_JPRB , 15.1479_JPRB, 14.7802_JPRB    , 13.6085_JPRB    , &
 & 12.2793_JPRB , 10.6929_JPRB, 8.72723_JPRB    , 6.57114_JPRB    , &
 & 4.46330_JPRB ,0.486724_JPRB,0.401446_JPRB    ,0.318879_JPRB    , &
 & 0.233959_JPRB,0.149119_JPRB, 5.62309E-02_JPRB, 7.83925E-03_JPRB /)  

!     Rayleigh extinction coefficient at v = 6925 cm-1.
RAYL = 9.41E-09_JPRB

STRRAT = 0.0045321_JPRB

LAYREFFR = 8

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


FORREF(:, 1) = (/ 0.110008E-06_JPRB, 0.630912E-06_JPRB, 0.363159E-05_JPRB, 0.616892E-05_JPRB /)
FORREF(:, 2) = (/ 0.429709E-05_JPRB, 0.789174E-05_JPRB, 0.217416E-04_JPRB, 0.639393E-04_JPRB /)
FORREF(:, 3) = (/ 0.436283E-04_JPRB, 0.526247E-04_JPRB, 0.116341E-03_JPRB, 0.205616E-03_JPRB /)
FORREF(:, 4) = (/ 0.215627E-03_JPRB, 0.234522E-03_JPRB, 0.280497E-03_JPRB, 0.838668E-03_JPRB /)
FORREF(:, 5) = (/ 0.529283E-03_JPRB, 0.620848E-03_JPRB, 0.935561E-03_JPRB, 0.171252E-02_JPRB /)
FORREF(:, 6) = (/ 0.212267E-02_JPRB, 0.218564E-02_JPRB, 0.222227E-02_JPRB, 0.199650E-02_JPRB /)
FORREF(:, 7) = (/ 0.291120E-02_JPRB, 0.281168E-02_JPRB, 0.259543E-02_JPRB, 0.210159E-02_JPRB /)
FORREF(:, 8) = (/ 0.316249E-02_JPRB, 0.310695E-02_JPRB, 0.279501E-02_JPRB, 0.208076E-02_JPRB /)
FORREF(:, 9) = (/ 0.354993E-02_JPRB, 0.336989E-02_JPRB, 0.298930E-02_JPRB, 0.180424E-02_JPRB /)
FORREF(:,10) = (/ 0.397729E-02_JPRB, 0.367409E-02_JPRB, 0.328982E-02_JPRB, 0.177807E-02_JPRB /)
FORREF(:,11) = (/ 0.408831E-02_JPRB, 0.398792E-02_JPRB, 0.352727E-02_JPRB, 0.192470E-02_JPRB /)
FORREF(:,12) = (/ 0.433926E-02_JPRB, 0.420667E-02_JPRB, 0.383894E-02_JPRB, 0.220836E-02_JPRB /)
FORREF(:,13) = (/ 0.436397E-02_JPRB, 0.433769E-02_JPRB, 0.425752E-02_JPRB, 0.237343E-02_JPRB /)
FORREF(:,14) = (/ 0.440525E-02_JPRB, 0.449018E-02_JPRB, 0.451881E-02_JPRB, 0.269169E-02_JPRB /)
FORREF(:,15) = (/ 0.491350E-02_JPRB, 0.481760E-02_JPRB, 0.475799E-02_JPRB, 0.362666E-02_JPRB /)
FORREF(:,16) = (/ 0.561641E-02_JPRB, 0.524553E-02_JPRB, 0.512473E-02_JPRB, 0.493802E-02_JPRB /)

!     -----------------------------------------------------------------
!     The array SELFREF contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

SELFREF(:, 1) = (/ &
 & 0.115887E-03_JPRB, 0.926537E-04_JPRB, 0.740783E-04_JPRB, 0.592270E-04_JPRB, 0.473530E-04_JPRB, &
 & 0.378596E-04_JPRB, 0.302694E-04_JPRB, 0.242010E-04_JPRB, 0.193491E-04_JPRB, 0.154700E-04_JPRB /)  
SELFREF(:, 2) = (/ &
 & 0.459557E-03_JPRB, 0.381962E-03_JPRB, 0.317469E-03_JPRB, 0.263866E-03_JPRB, 0.219313E-03_JPRB, &
 & 0.182283E-03_JPRB, 0.151505E-03_JPRB, 0.125924E-03_JPRB, 0.104662E-03_JPRB, 0.869904E-04_JPRB /)  
SELFREF(:, 3) = (/ &
 & 0.166821E-02_JPRB, 0.151103E-02_JPRB, 0.136866E-02_JPRB, 0.123970E-02_JPRB, 0.112290E-02_JPRB, &
 & 0.101710E-02_JPRB, 0.921266E-03_JPRB, 0.834463E-03_JPRB, 0.755839E-03_JPRB, 0.684623E-03_JPRB /)  
SELFREF(:, 4) = (/ &
 & 0.460175E-02_JPRB, 0.421372E-02_JPRB, 0.385842E-02_JPRB, 0.353307E-02_JPRB, 0.323516E-02_JPRB, &
 & 0.296236E-02_JPRB, 0.271257E-02_JPRB, 0.248385E-02_JPRB, 0.227440E-02_JPRB, 0.208262E-02_JPRB /)  
SELFREF(:, 5) = (/ &
 & 0.101589E-01_JPRB, 0.924742E-02_JPRB, 0.841772E-02_JPRB, 0.766247E-02_JPRB, 0.697497E-02_JPRB, &
 & 0.634917E-02_JPRB, 0.577951E-02_JPRB, 0.526096E-02_JPRB, 0.478893E-02_JPRB, 0.435926E-02_JPRB /)  
SELFREF(:, 6) = (/ &
 & 0.328043E-01_JPRB, 0.300853E-01_JPRB, 0.275917E-01_JPRB, 0.253048E-01_JPRB, 0.232075E-01_JPRB, &
 & 0.212839E-01_JPRB, 0.195198E-01_JPRB, 0.179020E-01_JPRB, 0.164182E-01_JPRB, 0.150574E-01_JPRB /)  
SELFREF(:, 7) = (/ &
 & 0.405936E-01_JPRB, 0.376032E-01_JPRB, 0.348331E-01_JPRB, 0.322671E-01_JPRB, 0.298901E-01_JPRB, &
 & 0.276883E-01_JPRB, 0.256486E-01_JPRB, 0.237591E-01_JPRB, 0.220089E-01_JPRB, 0.203876E-01_JPRB /)  
SELFREF(:, 8) = (/ &
 & 0.448362E-01_JPRB, 0.413811E-01_JPRB, 0.381923E-01_JPRB, 0.352492E-01_JPRB, 0.325329E-01_JPRB, &
 & 0.300259E-01_JPRB, 0.277121E-01_JPRB, 0.255766E-01_JPRB, 0.236056E-01_JPRB, 0.217866E-01_JPRB /)  
SELFREF(:, 9) = (/ &
 & 0.479741E-01_JPRB, 0.445389E-01_JPRB, 0.413497E-01_JPRB, 0.383889E-01_JPRB, 0.356400E-01_JPRB, &
 & 0.330880E-01_JPRB, 0.307188E-01_JPRB, 0.285191E-01_JPRB, 0.264770E-01_JPRB, 0.245812E-01_JPRB /)  
SELFREF(:,10) = (/ &
 & 0.519308E-01_JPRB, 0.484130E-01_JPRB, 0.451335E-01_JPRB, 0.420761E-01_JPRB, 0.392259E-01_JPRB, &
 & 0.365687E-01_JPRB, 0.340916E-01_JPRB, 0.317822E-01_JPRB, 0.296293E-01_JPRB, 0.276222E-01_JPRB /)  
SELFREF(:,11) = (/ &
 & 0.572039E-01_JPRB, 0.527780E-01_JPRB, 0.486945E-01_JPRB, 0.449270E-01_JPRB, 0.414510E-01_JPRB, &
 & 0.382439E-01_JPRB, 0.352849E-01_JPRB, 0.325549E-01_JPRB, 0.300361E-01_JPRB, 0.277122E-01_JPRB /)  
SELFREF(:,12) = (/ &
 & 0.601046E-01_JPRB, 0.554411E-01_JPRB, 0.511395E-01_JPRB, 0.471716E-01_JPRB, 0.435116E-01_JPRB, &
 & 0.401356E-01_JPRB, 0.370215E-01_JPRB, 0.341490E-01_JPRB, 0.314994E-01_JPRB, 0.290554E-01_JPRB /)  
SELFREF(:,13) = (/ &
 & 0.616595E-01_JPRB, 0.567145E-01_JPRB, 0.521662E-01_JPRB, 0.479826E-01_JPRB, 0.441346E-01_JPRB, &
 & 0.405951E-01_JPRB, 0.373395E-01_JPRB, 0.343450E-01_JPRB, 0.315906E-01_JPRB, 0.290571E-01_JPRB /)  
SELFREF(:,14) = (/ &
 & 0.647916E-01_JPRB, 0.592493E-01_JPRB, 0.541811E-01_JPRB, 0.495465E-01_JPRB, 0.453083E-01_JPRB, &
 & 0.414326E-01_JPRB, 0.378885E-01_JPRB, 0.346475E-01_JPRB, 0.316837E-01_JPRB, 0.289735E-01_JPRB /)  
SELFREF(:,15) = (/ &
 & 0.694231E-01_JPRB, 0.637703E-01_JPRB, 0.585777E-01_JPRB, 0.538079E-01_JPRB, 0.494265E-01_JPRB, &
 & 0.454019E-01_JPRB, 0.417050E-01_JPRB, 0.383091E-01_JPRB, 0.351897E-01_JPRB, 0.323244E-01_JPRB /)  
SELFREF(:,16) = (/ &
 & 0.761764E-01_JPRB, 0.701815E-01_JPRB, 0.646584E-01_JPRB, 0.595700E-01_JPRB, 0.548820E-01_JPRB, &
 & 0.505629E-01_JPRB, 0.465838E-01_JPRB, 0.429178E-01_JPRB, 0.395403E-01_JPRB, 0.364286E-01_JPRB /)  

IF (LHOOK) CALL DR_HOOK('SRTM_KGB21',1,ZHOOK_HANDLE)
RETURN

1001 CONTINUE
CALL ABOR1("SRTM_KGB21:ERROR READING FILE RADSRTM")

END SUBROUTINE SRTM_KGB21
