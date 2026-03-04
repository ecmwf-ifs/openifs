SUBROUTINE SRTM_KGB19

!     Originally by J.Delamere, Atmospheric & Environmental Research.
!     Revision: 2.4
!     BAND 16:  4650-5150 cm-1 (low - H2O,CO2; high - CO2)
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
USE YOESRTA19 , ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, RAYL, STRRAT, LAYREFFR, &
  &   KA_D, KB_D

!     ------------------------------------------------------------------

IMPLICIT NONE

! KURUCZ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SRTM_KGB19',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  READ(NULRAD,ERR=1001) KA_D,KB_D
  KA = REAL(KA_D,JPRB)
  KB = REAL(KB_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KA,MTAGRAD,1,CDSTRING='SRTM_KGB19:')
  CALL MPL_BROADCAST (KB,MTAGRAD,1,CDSTRING='SRTM_KGB19:')
ENDIF

SFLUXREF(:,1) = (/ &
 & 3.25791_JPRB    , 3.29697_JPRB    , 3.16031_JPRB    , 2.96115_JPRB    , &
 & 2.69238_JPRB    , 2.33819_JPRB    , 1.92760_JPRB    , 1.44918_JPRB    , &
 & 0.979764_JPRB   , 0.107336_JPRB   , 8.94523E-02_JPRB, 6.98325E-02_JPRB, &
 & 5.12051E-02_JPRB, 3.23645E-02_JPRB, 1.23401E-02_JPRB, 1.71339E-03_JPRB /)  
SFLUXREF(:,2) = (/ &
 & 3.22769_JPRB    , 3.28817_JPRB    , 3.16687_JPRB    , 2.97662_JPRB    , &
 & 2.69495_JPRB    , 2.34392_JPRB    , 1.92900_JPRB    , 1.45391_JPRB    , &
 & 0.982522_JPRB   , 0.107638_JPRB   , 8.92458E-02_JPRB, 6.99885E-02_JPRB, &
 & 5.09679E-02_JPRB, 3.23789E-02_JPRB, 1.22673E-02_JPRB, 1.56040E-03_JPRB /)  
SFLUXREF(:,3) = (/ &
 & 3.22294_JPRB    , 3.27780_JPRB    , 3.17424_JPRB    , 2.97143_JPRB    , &
 & 2.69785_JPRB    , 2.34993_JPRB    , 1.93155_JPRB    , 1.45196_JPRB    , &
 & 0.985329_JPRB   , 0.108027_JPRB   , 8.93552E-02_JPRB, 6.99937E-02_JPRB, &
 & 5.11678E-02_JPRB, 3.24846E-02_JPRB, 1.20636E-02_JPRB, 1.56040E-03_JPRB /)  
SFLUXREF(:,4) = (/ &
 & 3.22445_JPRB    , 3.26113_JPRB    , 3.18438_JPRB    , 2.96921_JPRB    , &
 & 2.69579_JPRB    , 2.35586_JPRB    , 1.93454_JPRB    , 1.44949_JPRB    , &
 & 0.987347_JPRB   , 0.108611_JPRB   , 8.91643E-02_JPRB, 7.02236E-02_JPRB, &
 & 5.12980E-02_JPRB, 3.25282E-02_JPRB, 1.21189E-02_JPRB, 1.56040E-03_JPRB /)  
SFLUXREF(:,5) = (/ &
 & 3.22497_JPRB    , 3.25109_JPRB    , 3.18741_JPRB    , 2.96970_JPRB    , &
 & 2.69460_JPRB    , 2.36020_JPRB    , 1.93301_JPRB    , 1.45224_JPRB    , &
 & 0.988564_JPRB   , 0.108255_JPRB   , 8.93830E-02_JPRB, 7.03655E-02_JPRB, &
 & 5.13017E-02_JPRB, 3.29414E-02_JPRB, 1.21189E-02_JPRB, 1.56040E-03_JPRB /)  
SFLUXREF(:,6) = (/ &
 & 3.22632_JPRB    , 3.24174_JPRB    , 3.18524_JPRB    , 2.97402_JPRB    , &
 & 2.69807_JPRB    , 2.35742_JPRB    , 1.93377_JPRB    , 1.45621_JPRB    , &
 & 0.988132_JPRB   , 0.108344_JPRB   , 8.93188E-02_JPRB, 7.04907E-02_JPRB, &
 & 5.17938E-02_JPRB, 3.31465E-02_JPRB, 1.21155E-02_JPRB, 1.56040E-03_JPRB /)  
SFLUXREF(:,7) = (/ &
 & 3.22793_JPRB    , 3.23589_JPRB    , 3.17720_JPRB    , 2.97869_JPRB    , &
 & 2.70293_JPRB    , 2.35436_JPRB    , 1.93557_JPRB    , 1.45868_JPRB    , &
 & 0.988654_JPRB   , 0.108198_JPRB   , 8.93375E-02_JPRB, 7.09790E-02_JPRB, &
 & 5.24733E-02_JPRB, 3.31298E-02_JPRB, 1.21126E-02_JPRB, 1.56040E-03_JPRB /)  
SFLUXREF(:,8) = (/ &
 & 3.22966_JPRB    , 3.24087_JPRB    , 3.15676_JPRB    , 2.98171_JPRB    , &
 & 2.70894_JPRB    , 2.34975_JPRB    , 1.93855_JPRB    , 1.46354_JPRB    , &
 & 0.988544_JPRB   , 0.108574_JPRB   , 9.02522E-02_JPRB, 7.12908E-02_JPRB, &
 & 5.24844E-02_JPRB, 3.31084E-02_JPRB, 1.21060E-02_JPRB, 1.56040E-03_JPRB /)  
SFLUXREF(:,9) = (/ &
 & 3.27240_JPRB    , 3.24666_JPRB    , 3.13886_JPRB    , 2.95238_JPRB    , &
 & 2.70190_JPRB    , 2.34460_JPRB    , 1.93948_JPRB    , 1.47111_JPRB    , &
 & 0.990821_JPRB   , 0.108730_JPRB   , 9.01625E-02_JPRB, 7.13261E-02_JPRB, &
 & 5.24813E-02_JPRB, 3.31083E-02_JPRB, 1.21126E-02_JPRB, 1.56040E-03_JPRB /)  

!     Rayleigh extinction coefficient at v = 4900 cm-1.
RAYL = 2.29E-09_JPRB

STRRAT = 5.49281_JPRB

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

  
FORREF(:, 1) = (/ 0.106275E-05_JPRB, 0.104185E-05_JPRB, 0.420154E-05_JPRB /)
FORREF(:, 2) = (/ 0.154343E-05_JPRB, 0.653193E-05_JPRB, 0.174596E-04_JPRB /)
FORREF(:, 3) = (/ 0.348917E-05_JPRB, 0.108420E-04_JPRB, 0.540849E-04_JPRB /)
FORREF(:, 4) = (/ 0.145822E-04_JPRB, 0.156027E-04_JPRB, 0.881263E-04_JPRB /)
FORREF(:, 5) = (/ 0.220204E-04_JPRB, 0.819892E-04_JPRB, 0.817937E-04_JPRB /)
FORREF(:, 6) = (/ 0.447840E-04_JPRB, 0.121116E-03_JPRB, 0.932635E-04_JPRB /)
FORREF(:, 7) = (/ 0.166516E-03_JPRB, 0.147640E-03_JPRB, 0.754029E-04_JPRB /)
FORREF(:, 8) = (/ 0.234756E-03_JPRB, 0.145934E-03_JPRB, 0.771734E-04_JPRB /)
FORREF(:, 9) = (/ 0.289207E-03_JPRB, 0.146768E-03_JPRB, 0.677806E-04_JPRB /)
FORREF(:,10) = (/ 0.334959E-03_JPRB, 0.125513E-03_JPRB, 0.636648E-04_JPRB /)
FORREF(:,11) = (/ 0.333755E-03_JPRB, 0.136575E-03_JPRB, 0.593651E-04_JPRB /)
FORREF(:,12) = (/ 0.340042E-03_JPRB, 0.116259E-03_JPRB, 0.595192E-04_JPRB /)
FORREF(:,13) = (/ 0.422470E-03_JPRB, 0.148691E-03_JPRB, 0.630266E-04_JPRB /)
FORREF(:,14) = (/ 0.440655E-03_JPRB, 0.461917E-04_JPRB, 0.108222E-04_JPRB /)
FORREF(:,15) = (/ 0.486207E-03_JPRB, 0.428458E-03_JPRB, 0.108086E-04_JPRB /)
FORREF(:,16) = (/ 0.657463E-03_JPRB, 0.657446E-03_JPRB, 0.126190E-04_JPRB /)

!     -----------------------------------------------------------------
!     The array SELFREF contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).
     
SELFREF(:, 1) = (/ &
 & 0.331728E-03_JPRB, 0.287480E-03_JPRB, 0.249135E-03_JPRB, 0.215904E-03_JPRB, 0.187106E-03_JPRB, &
 & 0.162149E-03_JPRB, 0.140520E-03_JPRB, 0.121777E-03_JPRB, 0.105534E-03_JPRB, 0.914573E-04_JPRB /)  
SELFREF(:, 2) = (/ &
 & 0.882628E-03_JPRB, 0.698914E-03_JPRB, 0.553439E-03_JPRB, 0.438244E-03_JPRB, 0.347026E-03_JPRB, &
 & 0.274795E-03_JPRB, 0.217598E-03_JPRB, 0.172306E-03_JPRB, 0.136442E-03_JPRB, 0.108042E-03_JPRB /)  
SELFREF(:, 3) = (/ &
 & 0.115461E-02_JPRB, 0.937203E-03_JPRB, 0.760730E-03_JPRB, 0.617486E-03_JPRB, 0.501215E-03_JPRB, &
 & 0.406837E-03_JPRB, 0.330231E-03_JPRB, 0.268049E-03_JPRB, 0.217576E-03_JPRB, 0.176607E-03_JPRB /)  
SELFREF(:, 4) = (/ &
 & 0.103450E-02_JPRB, 0.960268E-03_JPRB, 0.891360E-03_JPRB, 0.827397E-03_JPRB, 0.768024E-03_JPRB, &
 & 0.712911E-03_JPRB, 0.661754E-03_JPRB, 0.614267E-03_JPRB, 0.570188E-03_JPRB, 0.529272E-03_JPRB /)  
SELFREF(:, 5) = (/ &
 & 0.289040E-02_JPRB, 0.240129E-02_JPRB, 0.199495E-02_JPRB, 0.165737E-02_JPRB, 0.137692E-02_JPRB, &
 & 0.114392E-02_JPRB, 0.950351E-03_JPRB, 0.789535E-03_JPRB, 0.655933E-03_JPRB, 0.544938E-03_JPRB /)  
SELFREF(:, 6) = (/ &
 & 0.361772E-02_JPRB, 0.306611E-02_JPRB, 0.259861E-02_JPRB, 0.220239E-02_JPRB, 0.186659E-02_JPRB, &
 & 0.158198E-02_JPRB, 0.134077E-02_JPRB, 0.113634E-02_JPRB, 0.963078E-03_JPRB, 0.816234E-03_JPRB /)  
SELFREF(:, 7) = (/ &
 & 0.329878E-02_JPRB, 0.318245E-02_JPRB, 0.307021E-02_JPRB, 0.296194E-02_JPRB, 0.285749E-02_JPRB, &
 & 0.275671E-02_JPRB, 0.265950E-02_JPRB, 0.256571E-02_JPRB, 0.247522E-02_JPRB, 0.238793E-02_JPRB /)  
SELFREF(:, 8) = (/ &
 & 0.293562E-02_JPRB, 0.300077E-02_JPRB, 0.306737E-02_JPRB, 0.313544E-02_JPRB, 0.320503E-02_JPRB, &
 & 0.327615E-02_JPRB, 0.334886E-02_JPRB, 0.342318E-02_JPRB, 0.349915E-02_JPRB, 0.357680E-02_JPRB /)  
SELFREF(:, 9) = (/ &
 & 0.281453E-02_JPRB, 0.295894E-02_JPRB, 0.311076E-02_JPRB, 0.327038E-02_JPRB, 0.343818E-02_JPRB, &
 & 0.361459E-02_JPRB, 0.380006E-02_JPRB, 0.399504E-02_JPRB, 0.420002E-02_JPRB, 0.441553E-02_JPRB /)  
SELFREF(:,10) = (/ &
 & 0.239488E-02_JPRB, 0.262487E-02_JPRB, 0.287696E-02_JPRB, 0.315325E-02_JPRB, 0.345607E-02_JPRB, &
 & 0.378798E-02_JPRB, 0.415176E-02_JPRB, 0.455048E-02_JPRB, 0.498749E-02_JPRB, 0.546647E-02_JPRB /)  
SELFREF(:,11) = (/ &
 & 0.271001E-02_JPRB, 0.292235E-02_JPRB, 0.315134E-02_JPRB, 0.339826E-02_JPRB, 0.366453E-02_JPRB, &
 & 0.395167E-02_JPRB, 0.426131E-02_JPRB, 0.459521E-02_JPRB, 0.495527E-02_JPRB, 0.534354E-02_JPRB /)  
SELFREF(:,12) = (/ &
 & 0.206702E-02_JPRB, 0.232254E-02_JPRB, 0.260966E-02_JPRB, 0.293226E-02_JPRB, 0.329475E-02_JPRB, &
 & 0.370204E-02_JPRB, 0.415969E-02_JPRB, 0.467391E-02_JPRB, 0.525169E-02_JPRB, 0.590090E-02_JPRB /)  
SELFREF(:,13) = (/ &
 & 0.227023E-02_JPRB, 0.257331E-02_JPRB, 0.291685E-02_JPRB, 0.330626E-02_JPRB, 0.374766E-02_JPRB, &
 & 0.424799E-02_JPRB, 0.481511E-02_JPRB, 0.545794E-02_JPRB, 0.618660E-02_JPRB, 0.701253E-02_JPRB /)  
SELFREF(:,14) = (/ &
 & 0.851078E-03_JPRB, 0.111512E-02_JPRB, 0.146109E-02_JPRB, 0.191439E-02_JPRB, 0.250832E-02_JPRB, &
 & 0.328653E-02_JPRB, 0.430617E-02_JPRB, 0.564215E-02_JPRB, 0.739261E-02_JPRB, 0.968616E-02_JPRB /)  
SELFREF(:,15) = (/ &
 & 0.742711E-02_JPRB, 0.721347E-02_JPRB, 0.700598E-02_JPRB, 0.680446E-02_JPRB, 0.660873E-02_JPRB, &
 & 0.641863E-02_JPRB, 0.623400E-02_JPRB, 0.605468E-02_JPRB, 0.588052E-02_JPRB, 0.571137E-02_JPRB /)  
SELFREF(:,16) = (/ &
 & 0.107170E-01_JPRB, 0.101913E-01_JPRB, 0.969138E-02_JPRB, 0.921599E-02_JPRB, 0.876392E-02_JPRB, &
 & 0.833402E-02_JPRB, 0.792521E-02_JPRB, 0.753646E-02_JPRB, 0.716677E-02_JPRB, 0.681522E-02_JPRB /)  

IF (LHOOK) CALL DR_HOOK('SRTM_KGB19',1,ZHOOK_HANDLE)
RETURN

1001 CONTINUE
CALL ABOR1("SRTM_KGB19:ERROR READING FILE RADSRTM")

END SUBROUTINE SRTM_KGB19
