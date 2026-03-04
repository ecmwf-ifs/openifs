SUBROUTINE SRTM_KGB17

!     Originally by J.Delamere, Atmospheric & Environmental Research.
!     Revision: 2.4
!     BAND 17:  3250-4000 cm-1 (low - H2O,CO2; high - H2O, CO2)
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
USE YOESRTA17 , ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, RAYL, STRRAT, LAYREFFR, &
  &   KA_D, KB_D

!     ------------------------------------------------------------------

IMPLICIT NONE

! KURUCZ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SRTM_KGB17',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  READ(NULRAD,ERR=1001) KA_D,KB_D
  KA = REAL(KA_D,JPRB)
  KB = REAL(KB_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KA,MTAGRAD,1,CDSTRING='SRTM_KGB17:')
  CALL MPL_BROADCAST (KB,MTAGRAD,1,CDSTRING='SRTM_KGB17:')
ENDIF

SFLUXREF(:,1) = (/ &
 & 3.15613_JPRB  ,  3.03449_JPRB  ,  2.92069_JPRB  ,  2.63874_JPRB   , &
 & 2.34581_JPRB  ,  2.06999_JPRB  ,  1.70906_JPRB  ,  1.29085_JPRB   , &
 & 0.874851_JPRB ,  0.0955392_JPRB,  0.0787813_JPRB,  0.0621951_JPRB , &
 & 0.0459076_JPRB,  0.0294129_JPRB,  0.0110387_JPRB,  0.00159668_JPRB /)  
  
SFLUXREF(:,2) = (/ &
 & 2.83147_JPRB  ,  2.95919_JPRB  ,  2.96674_JPRB  ,  2.77677_JPRB   , &
 & 2.46826_JPRB  ,  2.11481_JPRB  ,  1.73243_JPRB  ,  1.30279_JPRB   , &
 & 0.882714_JPRB ,  0.0962350_JPRB,  0.0802122_JPRB,  0.0636194_JPRB , &
 & 0.0472620_JPRB,  0.0299051_JPRB,  0.0110785_JPRB,  0.00159668_JPRB /)  
  
SFLUXREF(:,3) = (/ &
 & 2.82300_JPRB  ,  2.94845_JPRB  ,  2.95887_JPRB  ,  2.77593_JPRB   , &
 & 2.47096_JPRB  ,  2.12596_JPRB  ,  1.73847_JPRB  ,  1.30796_JPRB   , &
 & 0.884395_JPRB ,  0.0966936_JPRB,  0.0801996_JPRB,  0.0640199_JPRB , &
 & 0.0472803_JPRB,  0.0300515_JPRB,  0.0112366_JPRB,  0.00160814_JPRB /)  
  
SFLUXREF(:,4) = (/ &
 & 2.81715_JPRB  ,  2.93789_JPRB  ,  2.95091_JPRB  ,  2.77046_JPRB   , &
 & 2.47716_JPRB  ,  2.13591_JPRB  ,  1.74365_JPRB  ,  1.31277_JPRB   , &
 & 0.887443_JPRB ,  0.0967016_JPRB,  0.0803391_JPRB,  0.0642442_JPRB , &
 & 0.0472909_JPRB,  0.0300720_JPRB,  0.0114817_JPRB,  0.00161875_JPRB /)  
  
SFLUXREF(:,5) = (/ &
 & 2.82335_JPRB  ,  2.93168_JPRB  ,  2.91455_JPRB  ,  2.75213_JPRB   , &
 & 2.49168_JPRB  ,  2.14408_JPRB  ,  1.75726_JPRB  ,  1.32401_JPRB   , &
 & 0.893644_JPRB ,  0.0969523_JPRB,  0.0805197_JPRB,  0.0639936_JPRB , &
 & 0.0475099_JPRB,  0.0305667_JPRB,  0.0115372_JPRB,  0.00161875_JPRB /)  
     
!     Rayleigh extinction coefficient at v = 3625 cm-1.
RAYL = 6.86E-10_JPRB

STRRAT = 0.364641_JPRB

LAYREFFR = 30

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


FORREF(:, 1) = (/ 0.553258E-03_JPRB, 0.555486E-03_JPRB, 0.601339E-03_JPRB, 0.708280E-03_JPRB /)
FORREF(:, 2) = (/ 0.158558E-02_JPRB, 0.162957E-02_JPRB, 0.204991E-02_JPRB, 0.475881E-02_JPRB /)
FORREF(:, 3) = (/ 0.772542E-02_JPRB, 0.784562E-02_JPRB, 0.111979E-01_JPRB, 0.229016E-01_JPRB /)
FORREF(:, 4) = (/ 0.255097E-01_JPRB, 0.256272E-01_JPRB, 0.270691E-01_JPRB, 0.259505E-01_JPRB /)
FORREF(:, 5) = (/ 0.323263E-01_JPRB, 0.324495E-01_JPRB, 0.305535E-01_JPRB, 0.263993E-01_JPRB /)
FORREF(:, 6) = (/ 0.346920E-01_JPRB, 0.348255E-01_JPRB, 0.323586E-01_JPRB, 0.276357E-01_JPRB /)
FORREF(:, 7) = (/ 0.366509E-01_JPRB, 0.366412E-01_JPRB, 0.344434E-01_JPRB, 0.319223E-01_JPRB /)
FORREF(:, 8) = (/ 0.378451E-01_JPRB, 0.375341E-01_JPRB, 0.374369E-01_JPRB, 0.320334E-01_JPRB /)
FORREF(:, 9) = (/ 0.407348E-01_JPRB, 0.396203E-01_JPRB, 0.393988E-01_JPRB, 0.318343E-01_JPRB /)
FORREF(:,10) = (/ 0.433035E-01_JPRB, 0.426488E-01_JPRB, 0.408085E-01_JPRB, 0.332749E-01_JPRB /)
FORREF(:,11) = (/ 0.428254E-01_JPRB, 0.441151E-01_JPRB, 0.408887E-01_JPRB, 0.327077E-01_JPRB /)
FORREF(:,12) = (/ 0.443226E-01_JPRB, 0.446690E-01_JPRB, 0.404676E-01_JPRB, 0.350492E-01_JPRB /)
FORREF(:,13) = (/ 0.466103E-01_JPRB, 0.460809E-01_JPRB, 0.401286E-01_JPRB, 0.370427E-01_JPRB /)
FORREF(:,14) = (/ 0.483928E-01_JPRB, 0.477284E-01_JPRB, 0.380684E-01_JPRB, 0.387940E-01_JPRB /)
FORREF(:,15) = (/ 0.506987E-01_JPRB, 0.490016E-01_JPRB, 0.467069E-01_JPRB, 0.368998E-01_JPRB /)
FORREF(:,16) = (/ 0.510836E-01_JPRB, 0.522771E-01_JPRB, 0.500130E-01_JPRB, 0.483406E-01_JPRB /)

!     -----------------------------------------------------------------
!     The array SELFREF contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

SELFREF(:, 1) = (/ &
 & 0.160537E-01_JPRB, 0.149038E-01_JPRB, 0.138363E-01_JPRB, 0.128452E-01_JPRB, 0.119251E-01_JPRB, &
 & 0.110709E-01_JPRB, 0.102779E-01_JPRB, 0.954175E-02_JPRB, 0.885829E-02_JPRB, 0.822379E-02_JPRB /)  
SELFREF(:, 2) = (/ &
 & 0.365753E-01_JPRB, 0.342267E-01_JPRB, 0.320288E-01_JPRB, 0.299720E-01_JPRB, 0.280474E-01_JPRB, &
 & 0.262463E-01_JPRB, 0.245609E-01_JPRB, 0.229837E-01_JPRB, 0.215078E-01_JPRB, 0.201267E-01_JPRB /)  
SELFREF(:, 3) = (/ &
 & 0.127419E+00_JPRB, 0.118553E+00_JPRB, 0.110304E+00_JPRB, 0.102629E+00_JPRB, 0.954883E-01_JPRB, &
 & 0.888442E-01_JPRB, 0.826624E-01_JPRB, 0.769107E-01_JPRB, 0.715593E-01_JPRB, 0.665802E-01_JPRB /)  
SELFREF(:, 4) = (/ &
 & 0.378687E+00_JPRB, 0.348961E+00_JPRB, 0.321568E+00_JPRB, 0.296325E+00_JPRB, 0.273064E+00_JPRB, &
 & 0.251629E+00_JPRB, 0.231876E+00_JPRB, 0.213674E+00_JPRB, 0.196901E+00_JPRB, 0.181444E+00_JPRB /)  
SELFREF(:, 5) = (/ &
 & 0.472822E+00_JPRB, 0.435018E+00_JPRB, 0.400236E+00_JPRB, 0.368236E+00_JPRB, 0.338794E+00_JPRB, &
 & 0.311706E+00_JPRB, 0.286783E+00_JPRB, 0.263854E+00_JPRB, 0.242757E+00_JPRB, 0.223348E+00_JPRB /)  
SELFREF(:, 6) = (/ &
 & 0.505620E+00_JPRB, 0.465050E+00_JPRB, 0.427736E+00_JPRB, 0.393416E+00_JPRB, 0.361849E+00_JPRB, &
 & 0.332815E+00_JPRB, 0.306111E+00_JPRB, 0.281550E+00_JPRB, 0.258959E+00_JPRB, 0.238181E+00_JPRB /)  
SELFREF(:, 7) = (/ &
 & 0.530488E+00_JPRB, 0.487993E+00_JPRB, 0.448902E+00_JPRB, 0.412943E+00_JPRB, 0.379864E+00_JPRB, &
 & 0.349434E+00_JPRB, 0.321443E+00_JPRB, 0.295694E+00_JPRB, 0.272007E+00_JPRB, 0.250218E+00_JPRB /)  
SELFREF(:, 8) = (/ &
 & 0.540222E+00_JPRB, 0.497746E+00_JPRB, 0.458610E+00_JPRB, 0.422551E+00_JPRB, 0.389327E+00_JPRB, &
 & 0.358716E+00_JPRB, 0.330511E+00_JPRB, 0.304524E+00_JPRB, 0.280580E+00_JPRB, 0.258519E+00_JPRB /)  
SELFREF(:, 9) = (/ &
 & 0.565727E+00_JPRB, 0.522899E+00_JPRB, 0.483313E+00_JPRB, 0.446724E+00_JPRB, 0.412905E+00_JPRB, &
 & 0.381646E+00_JPRB, 0.352753E+00_JPRB, 0.326048E+00_JPRB, 0.301365E+00_JPRB, 0.278550E+00_JPRB /)  
SELFREF(:,10) = (/ &
 & 0.610122E+00_JPRB, 0.562337E+00_JPRB, 0.518295E+00_JPRB, 0.477702E+00_JPRB, 0.440289E+00_JPRB, &
 & 0.405806E+00_JPRB, 0.374023E+00_JPRB, 0.344730E+00_JPRB, 0.317730E+00_JPRB, 0.292846E+00_JPRB /)  
SELFREF(:,11) = (/ &
 & 0.645176E+00_JPRB, 0.588957E+00_JPRB, 0.537636E+00_JPRB, 0.490788E+00_JPRB, 0.448022E+00_JPRB, &
 & 0.408982E+00_JPRB, 0.373344E+00_JPRB, 0.340812E+00_JPRB, 0.311114E+00_JPRB, 0.284004E+00_JPRB /)  
SELFREF(:,12) = (/ &
 & 0.651737E+00_JPRB, 0.596547E+00_JPRB, 0.546031E+00_JPRB, 0.499792E+00_JPRB, 0.457469E+00_JPRB, &
 & 0.418730E+00_JPRB, 0.383272E+00_JPRB, 0.350816E+00_JPRB, 0.321108E+00_JPRB, 0.293916E+00_JPRB /)  
SELFREF(:,13) = (/ &
 & 0.661086E+00_JPRB, 0.607954E+00_JPRB, 0.559093E+00_JPRB, 0.514159E+00_JPRB, 0.472836E+00_JPRB, &
 & 0.434834E+00_JPRB, 0.399886E+00_JPRB, 0.367747E+00_JPRB, 0.338191E+00_JPRB, 0.311011E+00_JPRB /)  
SELFREF(:,14) = (/ &
 & 0.692554E+00_JPRB, 0.635574E+00_JPRB, 0.583282E+00_JPRB, 0.535293E+00_JPRB, 0.491251E+00_JPRB, &
 & 0.450834E+00_JPRB, 0.413741E+00_JPRB, 0.379701E+00_JPRB, 0.348461E+00_JPRB, 0.319791E+00_JPRB /)  
SELFREF(:,15) = (/ &
 & 0.714646E+00_JPRB, 0.657179E+00_JPRB, 0.604334E+00_JPRB, 0.555737E+00_JPRB, 0.511049E+00_JPRB, &
 & 0.469954E+00_JPRB, 0.432164E+00_JPRB, 0.397412E+00_JPRB, 0.365455E+00_JPRB, 0.336068E+00_JPRB /)  
SELFREF(:,16) = (/ &
 & 0.782126E+00_JPRB, 0.710682E+00_JPRB, 0.645764E+00_JPRB, 0.586776E+00_JPRB, 0.533177E+00_JPRB, &
 & 0.484473E+00_JPRB, 0.440219E+00_JPRB, 0.400007E+00_JPRB, 0.363468E+00_JPRB, 0.330266E+00_JPRB /)  

IF (LHOOK) CALL DR_HOOK('SRTM_KGB17',1,ZHOOK_HANDLE)
RETURN

1001 CONTINUE
CALL ABOR1("SRTM_KGB17:ERROR READING FILE RADSRTM")

END SUBROUTINE SRTM_KGB17
