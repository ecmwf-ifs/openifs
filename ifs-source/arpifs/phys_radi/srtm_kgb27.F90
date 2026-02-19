SUBROUTINE SRTM_KGB27

!     Originally by J.Delamere, Atmospheric & Environmental Research.
!     Revision: 2.4
!     BAND 16: 29000-38000 cm-1 (low - O3; high - O3)
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
USE YOESRTA27 , ONLY : KA, KB, SFLUXREF, RAYL, SCALEKUR, LAYREFFR, &
  &  KA_D, KB_D

!     ------------------------------------------------------------------

IMPLICIT NONE

! KURUCZ
!     The following values were obtained using the "low resolution"
!     version of the Kurucz solar source function.  For unknown reasons,
!     the total irradiance in this band differs from the corresponding
!     total in the "high-resolution" version of the Kurucz function.
!     Therefore, below these values are scaled by the factor SCALEKUR.
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SRTM_KGB27',0,ZHOOK_HANDLE)

IF( MYPROC==1 )THEN
  READ(NULRAD,ERR=1001) KA_D,KB_D
  KA = REAL(KA_D,JPRB)
  KB = REAL(KB_D,JPRB)
ENDIF
IF( NPROC>1 )THEN
  CALL MPL_BROADCAST (KA,MTAGRAD,1,CDSTRING='SRTM_KGB27:')
  CALL MPL_BROADCAST (KB,MTAGRAD,1,CDSTRING='SRTM_KGB27:')
ENDIF

SFLUXREF = (/ &
 & 14.0526_JPRB    , 11.4794_JPRB    , 8.72590_JPRB    , 5.56966_JPRB    , &
 & 3.80927_JPRB    , 1.57690_JPRB    , 1.15099_JPRB    , 1.10012_JPRB    , &
 & 0.658212_JPRB   , 5.86859E-02_JPRB, 5.56186E-02_JPRB, 4.68040E-02_JPRB, &
 & 3.64897E-02_JPRB, 3.58053E-02_JPRB, 1.38130E-02_JPRB, 1.90193E-03_JPRB /)  

!     Rayleigh extinction coefficient at v = 2925 cm-1.
RAYL = (/ &
 & 3.44534E-06_JPRB,4.14480E-06_JPRB,4.95069E-06_JPRB,5.81204E-06_JPRB, &
 & 6.69748E-06_JPRB,7.56488E-06_JPRB,8.36344E-06_JPRB,9.04135E-06_JPRB, &
 & 9.58324E-06_JPRB,9.81542E-06_JPRB,9.75119E-06_JPRB,9.74533E-06_JPRB, &
 & 9.74139E-06_JPRB,9.73525E-06_JPRB,9.73577E-06_JPRB,9.73618E-06_JPRB /)  

SCALEKUR = 50.15_JPRB/48.37_JPRB

LAYREFFR = 32

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
  
     
!     -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_KGB27',1,ZHOOK_HANDLE)
RETURN

1001 CONTINUE
CALL ABOR1("SRTM_KGB27:ERROR READING FILE RADSRTM")

END SUBROUTINE SRTM_KGB27
