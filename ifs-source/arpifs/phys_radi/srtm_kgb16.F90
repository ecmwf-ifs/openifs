SUBROUTINE SRTM_KGB16

!     Originally by J.Delamere, Atmospheric & Environmental Research.
!     Revision: 2.4
!     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)
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
USE YOESRTA16 , ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, RAYL, STRRAT1, LAYREFFR, &
  &  KA_D, KB_D

!     ------------------------------------------------------------------

IMPLICIT NONE

! KURUCZ
CHARACTER(LEN = 512) :: CLZZZ
CHARACTER(LEN = 512) :: CLF1
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SRTM_KGB16',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  CALL GET_ENVIRONMENT_VARIABLE("DATA",CLZZZ)
  IF(CLZZZ /= " ") THEN
    CLF1=TRIM(CLZZZ)//"/ifsdata/RADSRTM"
    WRITE(0,'(1x,a)')'Reading '//TRIM(CLF1)
#ifdef LITTLE_ENDIAN
    OPEN(NULRAD,FILE=CLF1,FORM="UNFORMATTED",ACTION="READ",ERR=1000,CONVERT='BIG_ENDIAN')
#else
    OPEN(NULRAD,FILE=CLF1,FORM="UNFORMATTED",ACTION="READ",ERR=1000)
#endif
  ELSE
#ifdef LITTLE_ENDIAN
    OPEN(NULRAD,FILE='RADSRTM',FORM="UNFORMATTED",ACTION="READ",ERR=1000,CONVERT='BIG_ENDIAN')
#else
    OPEN(NULRAD,FILE='RADSRTM',FORM="UNFORMATTED",ACTION="READ",ERR=1000)
#endif
  ENDIF
  READ(NULRAD,ERR=1001) KA_D,KB_D
  KA = REAL(KA_D,JPRB)
  KB = REAL(KB_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KA,MTAGRAD,1,CDSTRING='SRTM_KGB16:')
  CALL MPL_BROADCAST (KB,MTAGRAD,1,CDSTRING='SRTM_KGB16:')
ENDIF

SFLUXREF = (/ &
 & 1.92269_JPRB    , 1.72844_JPRB    , 1.64326_JPRB    , 1.58451_JPRB     &
 & , 1.44031_JPRB    , 1.25108_JPRB    , 1.02724_JPRB    , 0.776759_JPRB    &
 & , 0.534444_JPRB   , 5.87755E-02_JPRB, 4.86706E-02_JPRB, 3.87989E-02_JPRB &
 & , 2.84532E-02_JPRB, 1.82431E-02_JPRB, 6.92320E-03_JPRB, 9.70770E-04_JPRB /)  

!     Rayleigh extinction coefficient at v = 2925 cm-1.
RAYL = 2.91E-10_JPRB

STRRAT1 = 252.131_JPRB

LAYREFFR = 18

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

FORREF(:, 1) = (/ 0.525585E-05_JPRB, 0.527618E-05_JPRB, 0.746929E-04_JPRB /)
FORREF(:, 2) = (/ 0.794660E-05_JPRB, 0.136902E-04_JPRB, 0.849878E-04_JPRB /)
FORREF(:, 3) = (/ 0.197099E-04_JPRB, 0.733094E-04_JPRB, 0.121687E-03_JPRB /)
FORREF(:, 4) = (/ 0.148274E-03_JPRB, 0.169776E-03_JPRB, 0.164848E-03_JPRB /)
FORREF(:, 5) = (/ 0.230296E-03_JPRB, 0.210384E-03_JPRB, 0.182028E-03_JPRB /)
FORREF(:, 6) = (/ 0.280575E-03_JPRB, 0.259217E-03_JPRB, 0.196080E-03_JPRB /)
FORREF(:, 7) = (/ 0.329034E-03_JPRB, 0.291575E-03_JPRB, 0.207044E-03_JPRB /)
FORREF(:, 8) = (/ 0.349989E-03_JPRB, 0.323471E-03_JPRB, 0.225712E-03_JPRB /)
FORREF(:, 9) = (/ 0.366097E-03_JPRB, 0.321519E-03_JPRB, 0.253150E-03_JPRB /)
FORREF(:,10) = (/ 0.383589E-03_JPRB, 0.355314E-03_JPRB, 0.262555E-03_JPRB /)
FORREF(:,11) = (/ 0.375933E-03_JPRB, 0.372443E-03_JPRB, 0.261313E-03_JPRB /)
FORREF(:,12) = (/ 0.370652E-03_JPRB, 0.382366E-03_JPRB, 0.250070E-03_JPRB /)
FORREF(:,13) = (/ 0.375092E-03_JPRB, 0.379542E-03_JPRB, 0.265794E-03_JPRB /)
FORREF(:,14) = (/ 0.389705E-03_JPRB, 0.384274E-03_JPRB, 0.322135E-03_JPRB /)
FORREF(:,15) = (/ 0.372084E-03_JPRB, 0.390422E-03_JPRB, 0.370035E-03_JPRB /)
FORREF(:,16) = (/ 0.437802E-03_JPRB, 0.373406E-03_JPRB, 0.373222E-03_JPRB /)

!     -----------------------------------------------------------------
!     The array SELFREF contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

SELFREF(:, 1) = (/ &
 & 0.126758E-02_JPRB, 0.105253E-02_JPRB, 0.873963E-03_JPRB, 0.725690E-03_JPRB, 0.602573E-03_JPRB, &
 & 0.500344E-03_JPRB, 0.415458E-03_JPRB, 0.344973E-03_JPRB, 0.286447E-03_JPRB, 0.237849E-03_JPRB /)  
SELFREF(:, 2) = (/ &
 & 0.144006E-02_JPRB, 0.118514E-02_JPRB, 0.975351E-03_JPRB, 0.802697E-03_JPRB, 0.660606E-03_JPRB, &
 & 0.543667E-03_JPRB, 0.447429E-03_JPRB, 0.368226E-03_JPRB, 0.303044E-03_JPRB, 0.249400E-03_JPRB /)  
SELFREF(:, 3) = (/ &
 & 0.294018E-02_JPRB, 0.227428E-02_JPRB, 0.175920E-02_JPRB, 0.136077E-02_JPRB, 0.105258E-02_JPRB, &
 & 0.814189E-03_JPRB, 0.629789E-03_JPRB, 0.487153E-03_JPRB, 0.376821E-03_JPRB, 0.291478E-03_JPRB /)  
SELFREF(:, 4) = (/ &
 & 0.395290E-02_JPRB, 0.348405E-02_JPRB, 0.307081E-02_JPRB, 0.270658E-02_JPRB, 0.238556E-02_JPRB, &
 & 0.210261E-02_JPRB, 0.185322E-02_JPRB, 0.163341E-02_JPRB, 0.143967E-02_JPRB, 0.126891E-02_JPRB /)  
SELFREF(:, 5) = (/ &
 & 0.419122E-02_JPRB, 0.385638E-02_JPRB, 0.354829E-02_JPRB, 0.326481E-02_JPRB, 0.300398E-02_JPRB, &
 & 0.276399E-02_JPRB, 0.254317E-02_JPRB, 0.234000E-02_JPRB, 0.215305E-02_JPRB, 0.198104E-02_JPRB /)  
SELFREF(:, 6) = (/ &
 & 0.495659E-02_JPRB, 0.456777E-02_JPRB, 0.420945E-02_JPRB, 0.387924E-02_JPRB, 0.357494E-02_JPRB, &
 & 0.329450E-02_JPRB, 0.303606E-02_JPRB, 0.279790E-02_JPRB, 0.257842E-02_JPRB, 0.237615E-02_JPRB /)  
SELFREF(:, 7) = (/ &
 & 0.526981E-02_JPRB, 0.490687E-02_JPRB, 0.456893E-02_JPRB, 0.425426E-02_JPRB, 0.396126E-02_JPRB, &
 & 0.368844E-02_JPRB, 0.343441E-02_JPRB, 0.319788E-02_JPRB, 0.297764E-02_JPRB, 0.277256E-02_JPRB /)  
SELFREF(:, 8) = (/ &
 & 0.575426E-02_JPRB, 0.531597E-02_JPRB, 0.491106E-02_JPRB, 0.453699E-02_JPRB, 0.419141E-02_JPRB, &
 & 0.387216E-02_JPRB, 0.357722E-02_JPRB, 0.330475E-02_JPRB, 0.305303E-02_JPRB, 0.282048E-02_JPRB /)  
SELFREF(:, 9) = (/ &
 & 0.549881E-02_JPRB, 0.514328E-02_JPRB, 0.481074E-02_JPRB, 0.449970E-02_JPRB, 0.420877E-02_JPRB, &
 & 0.393665E-02_JPRB, 0.368213E-02_JPRB, 0.344406E-02_JPRB, 0.322138E-02_JPRB, 0.301310E-02_JPRB /)  
SELFREF(:,10) = (/ &
 & 0.605357E-02_JPRB, 0.561246E-02_JPRB, 0.520349E-02_JPRB, 0.482432E-02_JPRB, 0.447278E-02_JPRB, &
 & 0.414686E-02_JPRB, 0.384469E-02_JPRB, 0.356453E-02_JPRB, 0.330479E-02_JPRB, 0.306398E-02_JPRB /)  
SELFREF(:,11) = (/ &
 & 0.640504E-02_JPRB, 0.587858E-02_JPRB, 0.539540E-02_JPRB, 0.495194E-02_JPRB, 0.454492E-02_JPRB, &
 & 0.417136E-02_JPRB, 0.382850E-02_JPRB, 0.351382E-02_JPRB, 0.322501E-02_JPRB, 0.295993E-02_JPRB /)  
SELFREF(:,12) = (/ &
 & 0.677803E-02_JPRB, 0.615625E-02_JPRB, 0.559152E-02_JPRB, 0.507859E-02_JPRB, 0.461271E-02_JPRB, &
 & 0.418957E-02_JPRB, 0.380524E-02_JPRB, 0.345617E-02_JPRB, 0.313913E-02_JPRB, 0.285116E-02_JPRB /)  
SELFREF(:,13) = (/ &
 & 0.690347E-02_JPRB, 0.627003E-02_JPRB, 0.569472E-02_JPRB, 0.517219E-02_JPRB, 0.469761E-02_JPRB, &
 & 0.426658E-02_JPRB, 0.387509E-02_JPRB, 0.351953E-02_JPRB, 0.319659E-02_JPRB, 0.290328E-02_JPRB /)  
SELFREF(:,14) = (/ &
 & 0.692680E-02_JPRB, 0.632795E-02_JPRB, 0.578087E-02_JPRB, 0.528109E-02_JPRB, 0.482452E-02_JPRB, &
 & 0.440742E-02_JPRB, 0.402638E-02_JPRB, 0.367828E-02_JPRB, 0.336028E-02_JPRB, 0.306977E-02_JPRB /)  
SELFREF(:,15) = (/ &
 & 0.754894E-02_JPRB, 0.681481E-02_JPRB, 0.615207E-02_JPRB, 0.555378E-02_JPRB, 0.501367E-02_JPRB, &
 & 0.452609E-02_JPRB, 0.408593E-02_JPRB, 0.368857E-02_JPRB, 0.332986E-02_JPRB, 0.300603E-02_JPRB /)  
SELFREF(:,16) = (/ &
 & 0.760689E-02_JPRB, 0.709755E-02_JPRB, 0.662232E-02_JPRB, 0.617891E-02_JPRB, 0.576519E-02_JPRB, &
 & 0.537917E-02_JPRB, 0.501899E-02_JPRB, 0.468293E-02_JPRB, 0.436938E-02_JPRB, 0.407682E-02_JPRB /)  

IF (LHOOK) CALL DR_HOOK('SRTM_KGB16',1,ZHOOK_HANDLE)
RETURN

1000 CONTINUE
CALL ABOR1("SRTM_KGB16:ERROR OPENING FILE RADSRTM")
1001 CONTINUE
CALL ABOR1("SRTM_KGB16:ERROR READING FILE RADSRTM")

END SUBROUTINE SRTM_KGB16
