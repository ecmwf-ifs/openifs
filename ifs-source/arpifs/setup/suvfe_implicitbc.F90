! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE SUVFE_IMPLICITBC(&
& KTBC, KBBC, KORDER, KBASIS, PB )

!**** *SUVFE_IMPLICITBC* 

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SUVFE_IMPLICITBC

!     Explicit arguments :
!     --------------------
!      * INPUT:
!        KTBC            : top boundary condition type required
!        KBBC            : bottom boundary condition type required
!        Additional info about KTBC and KBBC:
!        K[TBC|BBC](1) = 1 - value of function at relevant point is zero
!        K[TBC|BBC](2) = N - all derivatives up to Nth order at relevant
!                            point are zero
!        KORDER          : order of basis
!        KBASIS          : number of basis functions needed to satisfy BC 
!                           and represent input vector

!      * INPUT/OUTPUT:
!        PB              : spline basis to be adjusted here to fullfill implicit BCs

!     Implicit arguments :
!     --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ALADIN-NH documentation.

!     Author.
!     -------
!     Jozef Vivoda, SHMU/LACE

!     Modifications.
!     --------------
!     Original : 2017-09
!     ------------------------------------------------------------------

USE PARKIND1, ONLY : JPRB, JPIM
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN  , ONLY : NULOUT

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KTBC(2)
INTEGER(KIND=JPIM), INTENT(IN) :: KBBC(2)
INTEGER(KIND=JPIM), INTENT(IN) :: KORDER
INTEGER(KIND=JPIM), INTENT(IN) :: KBASIS

REAL(KIND=JPRB), INTENT(INOUT) :: PB(KBASIS, KORDER, KORDER)

INTEGER(KIND=JPIM) :: IB, II, IS, ITBC, ISB, IBBC
INTEGER(KIND=JPIM) :: II_MAX, II_MIN
INTEGER(KIND=JPIM) :: IB_MAX, IB_MIN

REAL(KIND=JPRB)    :: Z_TMP(KORDER, KORDER)
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

INTEGER(KIND=JPIM),EXTERNAL ::  F_INTERVAL
INTEGER(KIND=JPIM),EXTERNAL ::  F_SEGMENT

!-----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVFE_IMPLICITBC',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------

!*  1. TBC
!   ------

ITBC = KTBC(2) + 1

IF( KTBC(2) > 0 )THEN
  II_MIN = MAX(F_INTERVAL(ITBC,      1), KORDER)   
  II_MAX = MIN(F_INTERVAL(ITBC, KORDER), KBASIS)
  Z_TMP = 0.0_JPRB

  DO II = II_MIN, II_MAX  ! intervals of ITBC 
    ISB = F_SEGMENT(ITBC,II)         ! segment of itbc function
    IB_MIN = MAX(1, II - KORDER + 1)
    IB_MAX = MIN(ITBC, II)
    DO IB = IB_MIN, IB_MAX           ! ib-th basis function
      IS = F_SEGMENT(IB,II)         ! is-th segment of piecewise spline basis function
      Z_TMP(ISB, : ) = Z_TMP(ISB, : ) + PB(IB, IS, : )
    ENDDO
  ENDDO
 
  PB(ITBC, : , : ) = Z_TMP( : , : )
  DO IB = 1, ITBC - 1
    PB(IB, :, :) = 0.0
  ENDDO
ENDIF

IF( KTBC(1) > 0 )THEN
  DO IB = 1, ITBC
    PB(IB, :, :) = 0.0
  ENDDO
ENDIF

!*  2. BBC
!   ------

IBBC = KBASIS - KBBC(2) 

IF( KBBC(2) > 0 )THEN
  II_MIN = MAX(F_INTERVAL(IBBC,      1), KORDER)   
  II_MAX = MIN(F_INTERVAL(IBBC, KORDER), KBASIS)
  Z_TMP = 0.0_JPRB
  DO II = II_MIN, II_MAX  ! intervals of IBBC 
    ISB = F_SEGMENT(IBBC, II)             ! segment of itbc function
    IB_MIN = MAX(IBBC, II - KORDER + 1)
    IB_MAX = MIN(KBASIS, II)
    DO IB = IB_MIN, IB_MAX                ! ib-th basis function
      IS = F_SEGMENT(IB, II)              ! is-th segment of piecewise spline basis function
        Z_TMP(ISB, : ) = Z_TMP(ISB, : ) + PB(IB, IS, : )
        WRITE(NULOUT,'("(IMPLICIT_BC) II IB IS ISB",10I5)') II, IB, IS, ISB
     ENDDO
   ENDDO

   PB(IBBC, :, :) = Z_TMP(:, :)
   DO IB = IBBC + 1, KBASIS
     PB(IB, :, :) = 0.0
   ENDDO
ENDIF

IF (KBBC(1) > 0) THEN
  DO IB = IBBC, KBASIS
    PB(IB, :, :) = 0.0
  ENDDO
ENDIF

!-----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVFE_IMPLICITBC',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------
END SUBROUTINE SUVFE_IMPLICITBC

!----------------------------------------------------------

INTEGER FUNCTION F_SEGMENT(KBASIS, KINTERVAL)

USE PARKIND1, ONLY : JPIM

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KBASIS     ! index of basis
INTEGER(KIND=JPIM), INTENT(IN) :: KINTERVAL  ! index of interval

F_SEGMENT = KINTERVAL + 1 - KBASIS

END FUNCTION F_SEGMENT


!----------------------------------------------------------

INTEGER FUNCTION F_INTERVAL(KBASIS, KSEGMENT)

USE PARKIND1, ONLY : JPIM

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KBASIS    ! index of basis
INTEGER(KIND=JPIM), INTENT(IN) :: KSEGMENT  ! index of segment

F_INTERVAL = KSEGMENT + KBASIS - 1

END FUNCTION F_INTERVAL
