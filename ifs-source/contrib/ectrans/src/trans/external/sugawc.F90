! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE SUGAWC(KDGLG,PMU,PW)

!**** *SUGAWC* - Compute Gaussian latitudes and weights 

!     Purpose.
!     --------
!     Compute Gaussian latitudes and weights.

!**   Interface.
!     ----------
!     CALL SUGAWC(...)

!     Explicit arguments :
!     -------------------- 
!      INPUT:
!       KDGLG    - number of latitudes.

!      OUTPUT:
!       PMU      - sine of Gaussian latitudes.
!       PW       - Gaussian weights.

!     Method.
!     -------

!     Externals.  SUGAW
!     ----------  

!     Author.
!     -------
!        K. Yessad, from SUGAWA and SULEG (trans)
!        Original : May 2012

!     Modifications.
!     --------------
!      F. Vana  05-Mar-2015  Support for single precision

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRD, JPIM

!ifndef INTERFACE

USE SUGAW_MOD

!endif INTERFACE

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM) ,INTENT(IN)  :: KDGLG
REAL(KIND=JPRD)    ,INTENT(OUT) :: PMU(:)
REAL(KIND=JPRD)    ,INTENT(OUT) :: PW(:)

!ifndef INTERFACE

REAL(KIND=JPRD)    :: ZANM
INTEGER(KIND=JPIM) :: ISTART,IK,IODD,JN,JGL
REAL(KIND=JPRD) :: ZFN(0:KDGLG,0:KDGLG)
REAL(KIND=JPRD) :: ZFNN

!     ------------------------------------------------------------------

! * preliminary calculations to compute input quantities ZANM and ZFN
!   (k.y.: coded after what I found in tfl/module/suleg_mod.F90).
ISTART=1
! Belousov, Swarztrauber use ZFN(0,0)=SQRT(2._JPRD)
! IFS normalisation chosen to be 0.5*Integral(Pnm**2) = 1
ZFN(0,0)=2._JPRD
DO JN=ISTART,KDGLG
  ZFNN=ZFN(0,0)
  DO JGL=1,JN
    ZFNN=ZFNN*SQRT(1._JPRD-0.25_JPRD/REAL(JGL**2,JPRD))
  ENDDO
  IODD=MOD(JN,2)
  ZFN(JN,JN)=ZFNN
  DO JGL=2,JN-IODD,2
    ZFN(JN,JN-JGL)=ZFN(JN,JN-JGL+2)*REAL((JGL-1)*(2*JN-JGL+2),JPRD)/REAL(JGL*(2*JN-JGL+1),JPRD)
  ENDDO
ENDDO

ZANM=SQRT(REAL(2*KDGLG+1,JPRD)*REAL(KDGLG**2,JPRD)/REAL(2*KDGLG-1,JPRD))

! * call to SUGAW (output: PW, PMU):
CALL SUGAW(KDGLG,0,KDGLG,PMU,PW,ZANM,ZFN)

!     ------------------------------------------------------------------

!endif INTERFACE

END SUBROUTINE SUGAWC

