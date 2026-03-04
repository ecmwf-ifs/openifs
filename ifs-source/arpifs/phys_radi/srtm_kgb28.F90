SUBROUTINE SRTM_KGB28

!     Originally by J.Delamere, Atmospheric & Environmental Research.
!     Revision: 2.4
!     BAND 28: 38000-50000 cm-1 (low - O3,O2; high - O3,O2)
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
USE YOESRTA28 , ONLY : KA, KB, SFLUXREF, RAYL, STRRAT, LAYREFFR, KA_D, KB_D

!     ------------------------------------------------------------------

IMPLICIT NONE

! KURUCZ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SRTM_KGB28',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  READ(NULRAD,ERR=1001) KA_D,KB_D
  KA = REAL(KA_D,JPRB)
  KB = REAL(KB_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KA,MTAGRAD,1,CDSTRING='SRTM_KGB28:')
  CALL MPL_BROADCAST (KB,MTAGRAD,1,CDSTRING='SRTM_KGB28:')
ENDIF

SFLUXREF(:,1) = (/ &
 & 1.06156_JPRB    , 0.599910_JPRB   , 0.422462_JPRB   , 0.400077_JPRB   , &
 & 0.282221_JPRB   , 0.187893_JPRB   , 6.77357E-02_JPRB, 3.04572E-02_JPRB, &
 & 2.00442E-02_JPRB, 2.30786E-03_JPRB, 2.08824E-03_JPRB, 1.42604E-03_JPRB, &
 & 9.67384E-04_JPRB, 6.35362E-04_JPRB, 1.47727E-04_JPRB, 6.87639E-06_JPRB /)  
SFLUXREF(:,2) = (/ &
 & 1.07598_JPRB    , 0.585099_JPRB   , 0.422852_JPRB   , 0.400077_JPRB   , &
 & 0.282221_JPRB   , 0.187893_JPRB   , 6.69686E-02_JPRB, 3.09070E-02_JPRB, &
 & 2.02400E-02_JPRB, 2.47760E-03_JPRB, 1.89411E-03_JPRB, 1.41122E-03_JPRB, &
 & 1.12449E-03_JPRB, 5.73505E-04_JPRB, 2.04160E-04_JPRB, 1.58371E-05_JPRB /)  
SFLUXREF(:,3) = (/ &
 & 0.461647_JPRB   , 0.406113_JPRB   , 0.332506_JPRB   , 0.307508_JPRB   , &
 & 0.211167_JPRB   , 0.235457_JPRB   , 0.495886_JPRB   , 0.363921_JPRB   , &
 & 0.192700_JPRB   , 2.04678E-02_JPRB, 1.55407E-02_JPRB, 1.03882E-02_JPRB, &
 & 1.10778E-02_JPRB, 1.00504E-02_JPRB, 4.93497E-03_JPRB, 5.73410E-04_JPRB /)  
SFLUXREF(:,4) = (/ &
 & 0.132669_JPRB   , 0.175058_JPRB   , 0.359263_JPRB   , 0.388142_JPRB   , &
 & 0.350359_JPRB   , 0.475892_JPRB   , 0.489593_JPRB   , 0.408437_JPRB   , &
 & 0.221049_JPRB   , 1.94514E-02_JPRB, 1.54848E-02_JPRB, 1.44999E-02_JPRB, &
 & 1.44568E-02_JPRB, 1.00527E-02_JPRB, 4.95897E-03_JPRB, 5.73327E-04_JPRB /)  
SFLUXREF(:,5) = (/ &
 & 7.54800E-02_JPRB, 0.232246_JPRB   , 0.359263_JPRB   , 0.388142_JPRB   , &
 & 0.350359_JPRB   , 0.426317_JPRB   , 0.493485_JPRB   , 0.432016_JPRB   , &
 & 0.239203_JPRB   , 1.74951E-02_JPRB, 1.74477E-02_JPRB, 1.83566E-02_JPRB, &
 & 1.44818E-02_JPRB, 1.01048E-02_JPRB, 4.97487E-03_JPRB, 5.66831E-04_JPRB /)  

!     Rayleigh extinction coefficient at v = ????? cm-1.
RAYL = 2.02E-05_JPRB

STRRAT = 6.67029E-07_JPRB

LAYREFFR = 58

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

IF (LHOOK) CALL DR_HOOK('SRTM_KGB28',1,ZHOOK_HANDLE)
RETURN

1001 CONTINUE
CALL ABOR1("SRTM_KGB28:ERROR READING FILE RADSRTM")

END SUBROUTINE SRTM_KGB28
