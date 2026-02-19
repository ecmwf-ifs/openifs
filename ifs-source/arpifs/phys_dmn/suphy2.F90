! (C) Copyright 1989- Meteo-France.

!OPTIONS XOPT(NOEVAL)
SUBROUTINE SUPHY2(YDVAB,YDDIMV,YDPHY2,KULOUT)

!**** *SUPHY2*   - Initialize common YOMPHY2 physics controlling
!                  constants

!     Purpose.
!     --------
!           Initialize YOMPHY2, the common that contains the parameters
!           for the control part of the physics of the model.

!**   Interface.
!     ----------
!        *CALL* *SUPHY2(...)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMPHY2

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Documentation ARPEGE

!     Author.
!     -------
!      J.-F. Geleyn .
!      Original : 90-9-1

!     Modifications.
!     --------------
!      J.M. Piriou  : 2002-01-10 set default values to operational ones.
!      Modified by R. EL Khatib : 02-03-29 Control XMULAF<0 ; add LMULAF
!      Modified by D. Banciu    : 02-12-09 Introduction of XDAMP
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      E. Bazile    : 07-05-07 Allowed  XMULAF>=1 for TKE (LECT)
!      F. Bouyssel  : 07-05-07 Remove control on XMULAF>0
!      S. Riette    : 2009-03-25 HTKERAF is added
!      E. Bazile    : 2009-07-20 LRAFTKE is added
!      K. Yessad (Jan 2010): remove useless variables.
!     ------------------------------------------------------------------

USE YOMVERT  , ONLY : TVAB, TVETA, TVFE
USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULNAM
USE YOMCT0   , ONLY : LECMWF
USE YOMPHY2  , ONLY : TPHY2

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVAB)         , INTENT(IN)   :: YDVAB
TYPE(TDIMV)        , INTENT(IN)   :: YDDIMV
TYPE(TPHY2), TARGET, INTENT(INOUT):: YDPHY2
INTEGER(KIND=JPIM) , INTENT(IN)   :: KULOUT 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV

REAL(KIND=JPRB) :: ZVETAF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"

REAL(KIND=JPRB) , POINTER ::  HCLP
REAL(KIND=JPRB) , POINTER ::  GZ0RAF
REAL(KIND=JPRB) , POINTER ::  HTSHM
REAL(KIND=JPRB) , POINTER ::  XMUCVPP
LOGICAL , POINTER ::  LRAFTUR
REAL(KIND=JPRB) , POINTER ::  XDAMP
LOGICAL , POINTER ::  LRAFTKE
LOGICAL , POINTER ::  LWMOCLOUD
REAL(KIND=JPRB) , POINTER ::  HTKERAF
REAL(KIND=JPRB) , POINTER ::  HTSML
REAL(KIND=JPRB) , POINTER ::  XMULAF
LOGICAL , POINTER ::  LMULAF
REAL(KIND=JPRB) , POINTER ::  FACRAF
REAL(KIND=JPRB) , POINTER ::  FACRAFCV
REAL(KIND=JPRB) , POINTER ::  FACRAFDCAPE
REAL(KIND=JPRB) , POINTER ::  GCAPERAF
#include "namphy2.nam.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUPHY2',0,ZHOOK_HANDLE)
ASSOCIATE(TSPHY=>YDPHY2%TSPHY, NTSML=>YDPHY2%NTSML, NTSHM=>YDPHY2%NTSHM, &
 & HVCLS=>YDPHY2%HVCLS, HWMOLOW=>YDPHY2%HWMOLOW, HTCLS=>YDPHY2%HTCLS, &
 & HWMOHIGH=>YDPHY2%HWMOHIGH, &
 & NFLEVG=>YDDIMV%NFLEVG)
LRAFTKE => YDPHY2%LRAFTKE
HCLP => YDPHY2%HCLP
LWMOCLOUD => YDPHY2%LWMOCLOUD
GZ0RAF => YDPHY2%GZ0RAF
LRAFTUR => YDPHY2%LRAFTUR
HTSML => YDPHY2%HTSML
XDAMP => YDPHY2%XDAMP
XMULAF => YDPHY2%XMULAF
LMULAF => YDPHY2%LMULAF
HTSHM => YDPHY2%HTSHM
HTKERAF => YDPHY2%HTKERAF
XMUCVPP => YDPHY2%XMUCVPP
FACRAF => YDPHY2%FACRAF
FACRAFCV => YDPHY2%FACRAFCV
FACRAFDCAPE => YDPHY2%FACRAFDCAPE
GCAPERAF => YDPHY2%GCAPERAF
!     ------------------------------------------------------------------

!*       1.    Set default values.
!              -------------------

!        1.1 Set implicit default values

XMULAF=-1.75_JPRB
XMUCVPP=0._JPRB
XDAMP=0._JPRB
HCLP=1500._JPRB
HTCLS=2._JPRB
HVCLS=10._JPRB
HTSHM=0.450_JPRB
HTSML=0.785_JPRB
TSPHY=1._JPRB
LRAFTUR=.FALSE.
LRAFTKE=.FALSE.
GZ0RAF=10.0_JPRB
FACRAF=15.0_JPRB
FACRAFCV=0._JPRB
FACRAFDCAPE=0._JPRB
GCAPERAF=100._JPRB
LMULAF=.FALSE.
LWMOCLOUD=.FALSE.
HWMOHIGH=5000._JPRB
HWMOLOW=2500._JPRB
HTKERAF=20.0_JPRB

!        1.2 Modify default values according to LECMWF

IF (LECMWF) THEN
ELSE
  LRAFTUR=.TRUE.
ENDIF

!     Remark : values for TSPHY, NTSHM/ML are calculated and not set up.

!     ------------------------------------------------------------------

!*       2.    Modify default values.
!              ----------------------

CALL POSNAM(NULNAM,'NAMPHY2')
READ(NULNAM,NAMPHY2)
!     ------------------------------------------------------------------

!*       3.    Compute cloud transition indexes.
!              ---------------------------------

NTSHM=0
NTSML=0
DO JLEV=1,NFLEVG
  ZVETAF=(YDVAB%VALH(JLEV)+YDVAB%VBH(JLEV)+YDVAB%VALH(JLEV-1)+YDVAB%VBH(JLEV-1))*0.5_JPRB
  IF (ZVETAF <= HTSHM) THEN
    NTSHM=JLEV
  ENDIF
  IF (ZVETAF <= HTSML) THEN
    NTSML=JLEV
  ENDIF
ENDDO

!     ------------------------------------------------------------------

!*       4.    Print final values.
!              -------------------

WRITE(UNIT=KULOUT,FMT='('' COMMON YOMPHY2 '')')
WRITE(UNIT=KULOUT,FMT='('' XMUCVPP = '',E10.4,'' XMULAF = '',E10.4&
 & ,'' XDAMP = '',E10.4&
 & ,'' LMULAF = '',L2,/,'' HTCLS = '',E10.4&
 & ,'' HVCLS = '',E10.4,'' HCLP = '',E10.4,/&
 & ,'' LWMOCLOUD = '',L2,'' HWMOHIGH = '',E10.4,'' HWMOLOW = '',E10.4&
 & ,'' LRAFTUR = '',L2,'' GZ0RAF = '',E10.4,'' FACRAF = '',E10.4&
 & ,'' HTSHM = '',F8.4,'' NTSHM = '',I3,'' HTSML = '',F8.4&
 & ,'' NTSML = '',I3,'' LRAFTKE = '',L2,'' HTKERAF = '',E10.4&
 & )')&
 & XMUCVPP,XMULAF,XDAMP,LMULAF,HTCLS,HVCLS,HCLP,&
 & LWMOCLOUD,HWMOHIGH,HWMOLOW,LRAFTUR,GZ0RAF,FACRAF,&
 & HTSHM,NTSHM,HTSML,NTSML,LRAFTKE,HTKERAF
WRITE(UNIT=KULOUT,FMT='(100(A,E10.4))') 'FACRAFDCAPE=',FACRAFDCAPE,'FACRAFCV=',FACRAFCV &
 & ,'GCAPERAF=',GCAPERAF

!*       5.    Control
!              -------

IF ((XDAMP /= 0.0_JPRB).AND.(XMUCVPP /= 0.0_JPRB)) THEN
  WRITE(UNIT=KULOUT,FMT='(A)') 'INCONSISTENCY BETWEEN XDAMP AND XMUCVPP !'
  CALL ABOR1('XDAMP/=0. IMPLIES XMUCVPP=0.!...')
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUPHY2',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHY2
