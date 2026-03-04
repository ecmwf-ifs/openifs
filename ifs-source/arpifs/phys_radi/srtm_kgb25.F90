SUBROUTINE SRTM_KGB25

!     Originally by J.Delamere, Atmospheric & Environmental Research.
!     Revision: 2.4
!     BAND 25: 16000-22650 cm-1 (low - H2O; high - nothing)
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
USE YOESRTA25 , ONLY : KA, SFLUXREF, RAYL, ABSO3A, ABSO3B, LAYREFFR, KA_D

!     ------------------------------------------------------------------

IMPLICIT NONE

! KURUCZ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SRTM_KGB25',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  READ(NULRAD,ERR=1001) KA_D
  KA = REAL(KA_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KA,MTAGRAD,1,CDSTRING='SRTM_KGB25:')
ENDIF

SFLUXREF = (/ &
 & 42.6858_JPRB , 45.7720_JPRB, 44.9872_JPRB, 45.9662_JPRB    , &
 & 46.5458_JPRB , 41.6926_JPRB, 32.2893_JPRB, 24.0928_JPRB    , &
 & 16.7686_JPRB , 1.86048_JPRB, 1.54057_JPRB, 1.23503_JPRB    , &
 & 0.915085_JPRB,0.590099_JPRB,0.218622_JPRB, 3.21287E-02_JPRB /)  

!     Rayleigh extinction coefficient at v = 2925 cm-1.
RAYL = (/ &
 & 9.81132E-07_JPRB,8.25605E-07_JPRB,6.71302E-07_JPRB,5.53556E-07_JPRB,  &
 & 3.97383E-07_JPRB,3.68206E-07_JPRB,4.42379E-07_JPRB,4.57799E-07_JPRB, &
 & 4.22683E-07_JPRB,3.87113E-07_JPRB,3.79810E-07_JPRB,3.63192E-07_JPRB, &
 & 3.51921E-07_JPRB,3.34231E-07_JPRB,3.34294E-07_JPRB,3.32673E-07_JPRB /)  
     
ABSO3A = (/ &
 & 2.32664E-02_JPRB,5.76154E-02_JPRB,0.125389_JPRB,0.250158_JPRB, &
 & 0.378756_JPRB   ,0.402196_JPRB   ,0.352026_JPRB,0.352036_JPRB, &
 & 0.386253_JPRB   ,0.414598_JPRB   ,0.420079_JPRB,0.435471_JPRB, &
 & 0.445487_JPRB   ,0.459549_JPRB   ,0.452920_JPRB,0.456838_JPRB /)  

ABSO3B = (/      &
 & 1.76917E-02_JPRB,4.64185E-02_JPRB,1.03640E-01_JPRB,0.189469_JPRB, &
 & 0.303858_JPRB   ,0.400248_JPRB   ,0.447357_JPRB   ,0.470009_JPRB, &
 & 0.498673_JPRB   ,0.515696_JPRB   ,0.517053_JPRB   ,0.517930_JPRB, &
 & 0.518345_JPRB   ,0.524952_JPRB   ,0.508244_JPRB   ,0.468981_JPRB /)  

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
IF (LHOOK) CALL DR_HOOK('SRTM_KGB25',1,ZHOOK_HANDLE)
RETURN

1001 CONTINUE
CALL ABOR1("SRTM_KGB25:ERROR READING FILE RADSRTM")

END SUBROUTINE SRTM_KGB25
