! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUCHET(YDDIMV)

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMMP0    ,ONLY : NPROC
USE YOMCHET   ,ONLY : MAXPHVAR ,GCHET ,GCHETN
USE YOMLUN    ,ONLY : NULNAM, NULOUT, NULUSR3, NULUSR4, NULUSR5

IMPLICIT NONE

!**** *SUCHET

!     Purpose.
!     --------
!         Initializes the calculations of logarithms of physical tendencies 
!         absolute values frequency distributions.

!**   Interface.
!     ----------
!        *CALL* *SUCHET*

!----------------------------------------------------------------
!        Explicit arguments :     
!        --------------------     

!        Implicit arguments :
!        --------------------
!     Namelist parameters (NAMCHET): LCTFREQ, LCTCOOR, LCTPROF
!     Global variable: GCHET!

!   References
!   ----------

!   Author
!   ------
!   2004-03-01: T.Kovacic, J.M.Piriou, F.Bouyssel

!   Modifications
!   -------------
!   2011-06: M. Jerczynski - some cleaning to meet norms
!   T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!----------------------------------------------------------------

!     Local variables

TYPE(TDIMV), INTENT(IN) :: YDDIMV

INTEGER(KIND=JPIM) :: I, J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
LOGICAL :: LLCHET

#include "abor1.intfb.h"
#include "posnam.intfb.h"

#include "namchet.nam.h"
!--------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUCHET',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG)

!--------------------------------------------------------
!    1.  Default values for the NAMCHET variables
!--------------------------------------------------------

GCHETN%LFREQD   = .FALSE.
GCHETN%LCOORD   = .FALSE.
GCHETN%LPROFV   = .FALSE.
GCHETN%CPHYST   = '1111000'
GCHETN%CPTSEP   = 'ALL '
GCHETN%NCOORMAX = 100
GCHETN%NPROFMAX = 10
GCHETN%NLEV1    = 1
GCHETN%NLEV2    = NFLEVG
GCHETN%APREC    = 1.1_JPRB
GCHETN%AWIDTH   = 14.0_JPRB

GCHET%CVARNA(1) = 'PTENDU              '
GCHET%CVARNA(2) = 'PTENDV              '
GCHET%CVARNA(3) = 'PTENDH              '
GCHET%CVARNA(4) = 'PTENDQ              '
GCHET%CVARNA(5) = 'PTENDQI             '
GCHET%CVARNA(6) = 'PTENDQL             '
GCHET%CVARNA(7) = 'PTENDQR             '
GCHET%CVARNA(8) = 'PTENDQS             '

GCHETN%STARTCLASS(1) = 1.0E-12_JPRB
GCHETN%STARTCLASS(2) = 1.0E-12_JPRB
GCHETN%STARTCLASS(3) = 1.0E-09_JPRB
GCHETN%STARTCLASS(4) = 1.0E-16_JPRB
GCHETN%STARTCLASS(5) = 1.0E-16_JPRB
GCHETN%STARTCLASS(6) = 1.0E-16_JPRB
GCHETN%STARTCLASS(7) = 1.0E-16_JPRB
GCHETN%STARTCLASS(8) = 1.0E-16_JPRB

GCHETN%TENDTHRESH(1) = 2.6E-3_JPRB
GCHETN%TENDTHRESH(2) = 2.6E-3_JPRB
GCHETN%TENDTHRESH(3) = 2.3E-0_JPRB
GCHETN%TENDTHRESH(4) = 8.2E-7_JPRB
GCHETN%TENDTHRESH(5) = 8.2E-7_JPRB
GCHETN%TENDTHRESH(6) = 8.2E-7_JPRB
GCHETN%TENDTHRESH(7) = 8.2E-7_JPRB
GCHETN%TENDTHRESH(8) = 8.2E-7_JPRB

GCHET%CFICFREQ  = 'CHET.FREQ.lfa'
GCHET%CFICCOOR  = 'CHET.COORD.lfa'
GCHET%MAXCLASS  = 1000
GCHET%NULFREQ   = NULUSR3
GCHET%NULCOOR   = NULUSR4
GCHET%NULPROF   = NULUSR5
GCHET%NVPHT     = MAXPHVAR
GCHET%NPHTEXTH  = 0

!--------------------------------------------------------  
!    2.  Reading the namelist
!--------------------------------------------------------

CALL POSNAM(NULNAM,'NAMCHET')
READ (NULNAM,NAMCHET)

!--------------------------------------------------------
!    3. Setup of GCHET's components
!--------------------------------------------------------

!    3.1. Setup of frequency distributions

IF (GCHETN%LFREQD) THEN
  GCHET%NWIDTH = INT(GCHETN%AWIDTH/LOG10(GCHETN%APREC))
  IF (GCHET%NWIDTH>GCHET%MAXCLASS) CALL ABOR1('suchet/ERROR: WIDTH>MAXCLASS!')

  ALLOCATE( GCHET%NQDIST(GCHET%NVPHT,GCHET%NWIDTH+1) )
  GCHET%NQDIST = 0

  ALLOCATE( GCHET%ACLASS(GCHET%NVPHT,GCHET%NWIDTH) )
  CALL LFAOUV(GCHET%NULFREQ, GCHET%CFICFREQ, 'W')
  DO I = 1, GCHET%NVPHT
   GCHET%ACLASS(I,:)=(/(GCHETN%STARTCLASS(I)*GCHETN%APREC**J,J=1,GCHET%NWIDTH)/)
   CALL LFAECRR(GCHET%NULFREQ,GCHET%CVARNA(I),GCHET%ACLASS(I,:),GCHET%NWIDTH)
  ENDDO
  CALL LFAFER(GCHET%NULFREQ)
ENDIF

!    3.2. Setup of geographical coordinates output

IF ((GCHETN%LCOORD).OR.(GCHETN%LPROFV)) THEN
  CALL LFAOUV(GCHET%NULCOOR, GCHET%CFICCOOR, 'W')
  CALL LFAFER(GCHET%NULCOOR)
ENDIF

!    3.3. Setup of SCM profiles output

IF (GCHETN%LPROFV) THEN
  ALLOCATE(GCHET%APOINT(GCHETN%NCOORMAX))
ENDIF
  
!--------------------------------------------------------
!    4. Print final values.
!--------------------------------------------------------

LLCHET = ( GCHETN%LFREQD .OR. GCHETN%LCOORD .OR. GCHETN%LPROFV )
IF ( NPROC /= 1 .AND. LLCHET ) THEN
  WRITE(NULOUT, *) 'Check on physical tend. is implemented only for 1 proc!'
  CALL ABOR1('SUCHET: ABOR1 CALLED')
ENDIF

WRITE(UNIT=NULOUT,FMT='('' COMMON YOMCHET '')')
WRITE(UNIT=NULOUT,FMT='('' LFREQD = '',L2 &
 & ,'' LCOORD = '',L2,'' LPROFV = '',L2 &
 & ,'' CPHYST = '',A10,/,'' CPTSEP = '',A4 &
 & ,'' NCOORMAX = '',I5,'' NPROFMAX = '',I3,/&
 & ,'' NLEV1 = '',I3,'' NLEV2 = '',I3,/&
 & ,'' APREC = '',E10.4,'' AWIDTH = '',E10.4,/&
 & ,'' TENDTHRESH = '',7(E9.4,1X),/&
 & ,'' STARTCLASS = '',7(E9.4,1X),/&
 & )')&
 & GCHETN%LFREQD, GCHETN%LCOORD, GCHETN%LPROFV, &
 & GCHETN%CPHYST, GCHETN%CPTSEP, &
 & GCHETN%NCOORMAX, GCHETN%NPROFMAX, GCHETN%NLEV1, GCHETN%NLEV2, &
 & GCHETN%APREC, GCHETN%AWIDTH, GCHETN%TENDTHRESH, GCHETN%STARTCLASS

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUCHET',1,ZHOOK_HANDLE)

END SUBROUTINE SUCHET
