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

MODULE CONTROL_VECTORS_OPER_MOD

!   Purpose.
!   --------
!     This module contains the basic operators associated
!     to the CONTROL_VECTOR type.
!
!   Author.
!   -------
!     Y. Tremolet
!
!   Modifications.
!   --------------
!     Original   25-Jul-2004 Split from CONTROL_VECTORS
!     A. Deckmyn 01-Oct-2008 LAM wavelets
!     M. Fisher  15-Feb-2013 put globals into structures
!     Y. Michel, MF, June 2018 Extension of the control variable for sqrt EnVar
!     S. Massart 19-Feb-2019 Parameter optimisation
!     E. Holm    15-Jul-2019 Accommodate input control vector of higher as well as lower resolution
! ------------------------------------------------------------------

USE PARKIND1, ONLY: JPIM, JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE SPECTRAL_FIELDS_MOD, ONLY: SPECTRAL_FIELD, ASSIGNMENT(=), &
                         & ALLOCATE_SPEC, DEALLOCATE_SPEC
USE CONTROL_VECTORS_BASE_MIX, ONLY: CONTROL_VECTOR, CHECK, OPERATOR(.NEQV.)

IMPLICIT NONE
PRIVATE
PUBLIC ASSIGNMENT(=), INVERSE

! ------------------------------------------------------------------

INTERFACE ASSIGNMENT (=)
MODULE PROCEDURE ASSIGN_CV_CV, ASSIGN_SCALAR_CV, ASSIGN_AR_CV, ASSIGN_CV_AR
END INTERFACE

INTERFACE INVERSE
MODULE PROCEDURE INVERSE_CV, INVERSE_CV_CV
END INTERFACE

! ------------------------------------------------------------------
CONTAINS
! ------------------------------------------------------------------

SUBROUTINE ASSIGN_CV_CV(YDCV1,YDCV2)
TYPE(CONTROL_VECTOR), INTENT(INOUT) :: YDCV1
TYPE(CONTROL_VECTOR), INTENT(IN)    :: YDCV2
call abor1("assign_cv_cv should never be called")

END SUBROUTINE ASSIGN_CV_CV

! ------------------------------------------------------------------

SUBROUTINE ASSIGN_SCALAR_CV(YDCV,PSCALAR)
TYPE(CONTROL_VECTOR), INTENT(INOUT) :: YDCV
REAL(KIND=JPRB), INTENT(IN) :: PSCALAR
call abor1("assign_scalar_cv should never be called")

END SUBROUTINE ASSIGN_SCALAR_CV

! ------------------------------------------------------------------

SUBROUTINE ASSIGN_AR_CV(YDCV,PVEC)
TYPE(CONTROL_VECTOR), INTENT(INOUT) :: YDCV
REAL(KIND=JPRB), INTENT(IN) :: PVEC(:)
call abor1("assign_ar_cv should never be called")

END SUBROUTINE ASSIGN_AR_CV

! ------------------------------------------------------------------

SUBROUTINE ASSIGN_CV_AR(PVEC,YDCV)
REAL(KIND=JPRB), INTENT(OUT) :: PVEC(:)
TYPE(CONTROL_VECTOR), INTENT(IN) :: YDCV
call abor1("assign_cv_ar should never be called")

END SUBROUTINE ASSIGN_CV_AR

! ------------------------------------------------------------------

SUBROUTINE INVERSE_CV(YDCV)
TYPE(CONTROL_VECTOR), INTENT(INOUT) :: YDCV
call abor1("inverse_cv should never be called")

END SUBROUTINE INVERSE_CV

! ------------------------------------------------------------------

SUBROUTINE INVERSE_CV_CV(YDIN,YDOUT)
TYPE(CONTROL_VECTOR), INTENT(IN)    :: YDIN
TYPE(CONTROL_VECTOR), INTENT(INOUT) :: YDOUT
call abor1("inverse_cv_cv should never be called")

END SUBROUTINE INVERSE_CV_CV

! ------------------------------------------------------------------

END MODULE CONTROL_VECTORS_OPER_MOD
