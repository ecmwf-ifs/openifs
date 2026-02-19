! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SETRAN (KCONSEED,YD_RANDOM_STREAM,KOUTSEED,LDABSTIME,PTSTEP)

!*** *SETRAN* - Sets the seed for a random number stream

!   Purpose.
!   --------

!      To set the seed for a random number stream as a function of NINDAT, NSSSSS and KCONSEED
!      Also with option to use absolute time

!   Interface.
!   ----------

!      CALL SETRAN(...)

!      INPUT: KCONSEED                   - integer to control the seeding
!      OPTIONAL INPUT: LDABSTIME          - use absolute time
!      OPTIONAL INPUT: PTSTEP            - time step for absolute time calculation
!      OPTIONAL OUTPUT: YD_RANDOM_STREAM - initialized random number stream
!      OPTIONAL OUTPUT: KOUTSEED         - a seed constructed from KCONSEED and the model time

!   Externals.
!   ----------

!   Method.
!   -------

!      The seed is set to a function of the initial time of the run,
!      and of an input integer KCONSEED.
!      For dates chosen at random, the seeds are approximately
!      uniformly distributed between 1 and HUGE(0). A highly nonlinear
!      function is used to reduce the possibility of correlations
!      between random sequences generated for different initial dates.

!   Reference.
!   ----------
!      None yet!

!   Author.
!   -------
!    Mike Fisher *ECMWF*   01/02/94

!   Modifications.
!   --------------
!    Jan Haseler   10/04/02  -  fix for too large numbers on p690
!    M. Fisher     25/09/02  -  adapt for RANDOM_NUMBERS_MIX module
!    M.Hamrud      01-Oct-2003 CY28 Cleaning
!    M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!    M.Fisher      2004-05-06 - Added argument yd_random_stream
!    P.Bechtold    14/05/2012 replace 86400 by RDAY
!    K. Yessad (July 2014): Move some variables.
!    S. Lang   (March 2017): added option to use absolute time
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMRIP0  , ONLY : NINDAT   ,NSSSSS, RTIMST
USE YOMCST    ,ONLY : RDAY
USE YOMCT3   , ONLY : NSTEP
USE RANDOM_NUMBERS_MIX

USE YOMLUN   , ONLY : NULOUT

!-----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),       INTENT(IN)              :: KCONSEED 
LOGICAL,                  INTENT(IN),    OPTIONAL :: LDABSTIME
REAL(KIND=JPRB),          INTENT(IN),    OPTIONAL :: PTSTEP
TYPE(RANDOMNUMBERSTREAM), INTENT(INOUT), OPTIONAL :: YD_RANDOM_STREAM
INTEGER(KIND=JPIM),       INTENT(OUT),   OPTIONAL :: KOUTSEED
!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IDIGITS, IRADIX, IS, JDIGIT, ISCALE, ISEED

REAL(KIND=JPRD) :: ZIRR1, ZS, ZT, ZTIM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

LOGICAL         :: LLUSEABSTIME
!-----------------------------------------------------------------------

#include "fcttim.func.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SETRAN',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------------

IRADIX  = RADIX(ZTIM)
IDIGITS = DIGITS(ZTIM)
ZIRR1 = 0.5_JPRD*(SQRT(5._JPRD)-1.0_JPRD)

LLUSEABSTIME=.FALSE.

!--- generate a unique number from the date and the input KCONSEED

IF (PRESENT(LDABSTIME) .AND. PRESENT(PTSTEP)) THEN
  LLUSEABSTIME=LDABSTIME
ENDIF

IF (LLUSEABSTIME) THEN
  ZTIM = REAL(RTIMST,JPRD)+REAL(NSTEP,JPRD)*REAL(PTSTEP,JPRD)&
   & -1720994.5_JPRD &
   & -2581470.3_JPRD*KCONSEED  
ELSE
  ZTIM = RJUDAT(NCCAA(NINDAT),NMM(NINDAT),NDD(NINDAT))&
   & -1720994.5_JPRD + REAL(NSSSSS,JPRD)/RDAY &
   & -2581470.3_JPRD*KCONSEED  
ENDIF

!--- multiply by an irrational number to randomize the bits and scale
!--- to between 0 and 1.

ZTIM = FRACTION(ZIRR1*ABS(ZTIM))

!--- reverse the bits

ZS = 0.0_JPRD
ZT = ZTIM
DO JDIGIT=1,IDIGITS
  ZT = ZT*IRADIX
  IS = INT(ZT)
  ZT = ZT-IS
  ZS = (ZS+IS)/IRADIX
ENDDO

!--- Scale to an odd number between 0 and HUGE-100000000
!--- (Allow some headroom because some routines use setran to set an initial seed, 
!---  and then generate new seeds by incrementing.)

ISCALE=(HUGE(ISEED)-100000000)/2
ISEED = 1 + 2*INT( ISCALE*ZS )

!--- set/output the seed

IF (PRESENT(KOUTSEED)) THEN
  KOUTSEED=ISEED
ENDIF

IF (PRESENT(YD_RANDOM_STREAM)) THEN
  CALL INITIALIZE_RANDOM_NUMBERS (ISEED,YD_RANDOM_STREAM)
ENDIF

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SETRAN',1,ZHOOK_HANDLE)
END SUBROUTINE SETRAN
