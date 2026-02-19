! (C) Copyright 1989- Meteo-France.

SUBROUTINE GETFPLUN(KDOM,KUNIT)

!**** *GETFPLUN*  - GET FULLPOS LOGICAL UNIT NUMBERS

!     PURPOSE.
!     --------
!           To get the logical unit number of the file reserved to a Fullpos subdomain

!**   INTERFACE.
!     ----------
!       *CALL* *GETFPLUN(KDOM,KUNIT)*

!        EXPLICIT ARGUMENTS
!        --------------------
!            KDOM   : domain index
!            KUNIT  : logical unit number

!        IMPLICIT ARGUMENTS
!        --------------------
!         See modules above

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!        RYAD EL KHATIB *METEO-FRANCE*
!        ORIGINAL : 25-Feb-2016  from ini3wrfp

!     MODIFICATIONS.
!     --------------
!        R. El Khatib 17-May-2017 return only the logical unit number corresponding to the domain index given in argument
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULFP01  ,NULFP02  ,NULFP03 ,NULFP04  ,NULFP05  ,NULFP06  ,NULFP07  ,&
 & NULFP08 ,NULFP09  ,NULFP10  ,NULFP11  ,NULFP12 ,NULFP13  ,NULFP14  ,NULFP15  

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),  INTENT(IN)  :: KDOM
INTEGER(KIND=JPIM),  INTENT(OUT) :: KUNIT

!     ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GETFPLUN',0,ZHOOK_HANDLE)

SELECT CASE (KDOM)
  CASE (:0)
    CALL ABOR1('GETFPLUN : KDOM NEGATIVE OR NULL')
  CASE(01)
    KUNIT=NULFP01
  CASE(02)
    KUNIT=NULFP02
  CASE(03)
    KUNIT=NULFP03
  CASE(04)
    KUNIT=NULFP04
  CASE(05)
    KUNIT=NULFP05
  CASE(06)
    KUNIT=NULFP06
  CASE(07)
    KUNIT=NULFP07
  CASE(08)
    KUNIT=NULFP08
  CASE(09)
    KUNIT=NULFP09
  CASE(10)
    KUNIT=NULFP10
  CASE(11)
    KUNIT=NULFP11
  CASE(12)
    KUNIT=NULFP12
  CASE(13)
    KUNIT=NULFP13
  CASE(14)
    KUNIT=NULFP14
  CASE(15)
    KUNIT=NULFP15
  CASE DEFAULT
    CALL ABOR1('GETFPLUN : KDOM TOO BIG')
END SELECT

IF (LHOOK) CALL DR_HOOK('GETFPLUN',1,ZHOOK_HANDLE)

END SUBROUTINE GETFPLUN
