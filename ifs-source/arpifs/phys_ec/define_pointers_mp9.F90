! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE DEFINE_POINTERS_MP9(YGFL,YDDYN, YDVARS)

!**** *DEFINE_POINTERS_MP9* - ECMWF PHYSICS

!     Purpose.
!     --------
!           Define local pointers KY[X]_MP9 for the input data used in
!           the EC lagged physics. These pointers are equal either to
!           Y[X]%MP9 or Y[X]%MP according to LTWOTL and NCURRENT_ITER.
!           For example, in a SL2TL scheme, the input data which is
!           not t+dt data is ALWAYS "t" data (predictor and corrector
!           steps if PC scheme). In this case:
!           - predictor step: "t" data match with the pointer Y[X]%MP;
!             so KY[X]_MP9=Y[X]%MP.
!           - corrector step: "t" data match with the pointer Y[X]%MP9;
!             so KY[X]_MP9=Y[X]%MP9.
!           In a leap-frog scheme, the input data which is
!           not t+dt data is ALWAYS "t-dt" data (LPC_FULL is not yet
!           implemented when a leap-frog scheme is used with lagged
!           EC physics). In this case KY[X]_MP9=Y[X]%MP9.

!**   Interface.
!     ----------
!        *CALL* *DEFINE_POINTERS_MP9*

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Called by GP_MODEL

!     Reference.
!     ----------

!     Author.
!     -------
!        K. Yessad (CNRM/GMAP) and Deborah Salmond *ECMWF* (FEB 2006).

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYN   , ONLY : TDYN
USE YOM_YGFL , ONLY : TYPE_GFLD
USE FIELD_VARIABLES_MOD, ONLY: FIELD_VARIABLES

IMPLICIT NONE

TYPE(TDYN),INTENT(INOUT):: YDDYN
TYPE(TYPE_GFLD),INTENT(INOUT):: YGFL
TYPE(FIELD_VARIABLES), INTENT(INOUT)  :: YDVARS
INTEGER(KIND=JPIM) :: JGFL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('DEFINE_POINTERS_MP9',0,ZHOOK_HANDLE)
ASSOCIATE(NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER)
!     ------------------------------------------------------------------

!*       1.    Interface to global arrays.
!              ---------------------------

! * For a SL2TL (resp. leap-frog scheme) scheme,
!   the following pointers must ALWAYS match with
!   quantities at time t (resp t-dt).

DO JGFL=1,NUMFLDS
  IF(YCOMP(JGFL)%LT9)THEN
    IF (NCURRENT_ITER == 0) THEN
      YCOMP(JGFL)%MP9_PH=YCOMP(JGFL)%MP
    ELSE
      YCOMP(JGFL)%MP9_PH=YCOMP(JGFL)%MP9
    ENDIF
  ENDIF
ENDDO

! Update the FIELD to reflect the above change
IF (NCURRENT_ITER == 0) THEN
  CALL YDVARS%GFL_PH9TOT0()
ELSE
  CALL YDVARS%GFL_PH9TOT9()
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('DEFINE_POINTERS_MP9',1,ZHOOK_HANDLE)
END SUBROUTINE DEFINE_POINTERS_MP9
