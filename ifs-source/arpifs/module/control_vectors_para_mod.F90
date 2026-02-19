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

MODULE CONTROL_VECTORS_PARA_MOD

!   Purpose.
!   --------
!     This module contains the distributed routines associated
!     to the CONTROL_VECTOR type.
!
!   Author.
!   -------
!     Y. Tremolet
!
!   Modifications.
!   --------------
!     Original   25-Jul-2004 Split from CONTROL_VECTORS
!     M.Hamrud   25-Aug-2009 Remove use of DOT_PRODUCT_CTLVEC for most cases,
!                use simple dot-product or reproducible sum instead
!     M.Fisher   15-Feb-2013 Put globals into structures
!     Y. Michel, MF, June 2018 Extension of the control variable for sqrt EnVar
!     S. Massart 19-Feb-2019 Parameter optimisation
! ------------------------------------------------------------------

USE PARKIND1,   ONLY : JPRB, JPIM
USE YOMHOOK ,   ONLY : LHOOK, DR_HOOK, JPHOOK
USE MPL_MODULE, ONLY : MPL_ALLREDUCE
USE ORDER_INDEPENDENT_SUMMATION_MOD, ONLY  : ORDER_INDEP_DOT_PRODUCT
USE CONTROL_VECTORS_BASE_MIX, ONLY : CONTROL_VECTOR, CHECK
IMPLICIT NONE
PRIVATE
PUBLIC DOT_PRODUCT, MAXVAL, MINVAL, SUM, CTLVEC_SQNORM, CTLVEC_NORM

! ------------------------------------------------------------------

INTERFACE DOT_PRODUCT
MODULE PROCEDURE DOT_PRODUCT_CV_CV, DOT_PRODUCT_WEIGHT,DOT_PRODUCT_CVS_CV
END INTERFACE

INTERFACE MAXVAL
MODULE PROCEDURE MAXVAL_CV
END INTERFACE

INTERFACE MINVAL
MODULE PROCEDURE MINVAL_CV
END INTERFACE

INTERFACE SUM
MODULE PROCEDURE SUM_CV
END INTERFACE

INTERFACE CTLVEC_SQNORM
MODULE PROCEDURE CV_SQNORM1, CV_SQNORM4
END INTERFACE

! ------------------------------------------------------------------
CONTAINS
! ------------------------------------------------------------------

REAL(KIND=JPRB) FUNCTION DOT_PRODUCT_CV_CV(YDCV1,YDCV2)
TYPE (CONTROL_VECTOR), INTENT(IN) :: YDCV1,YDCV2
  call abor1("dot_product_cv_cv should never be called")

END FUNCTION DOT_PRODUCT_CV_CV

! ------------------------------------------------------------------
FUNCTION DOT_PRODUCT_CVS_CV(YDCV1,YDCV2)
TYPE (CONTROL_VECTOR), INTENT(IN) :: YDCV1(:),YDCV2
REAL(KIND=JPRB) :: DOT_PRODUCT_CVS_CV(SIZE(YDCV1))
  call abor1("dot_product_cvs_cv should never be called")

END FUNCTION DOT_PRODUCT_CVS_CV

! ------------------------------------------------------------------

REAL(KIND=JPRB) FUNCTION DOT_PRODUCT_WEIGHT(YDCV1,YDCV2,YDW)
TYPE (CONTROL_VECTOR), INTENT(IN) :: YDCV1,YDCV2,YDW
  call abor1("dot_product_weight should never be called")

END FUNCTION DOT_PRODUCT_WEIGHT

! ------------------------------------------------------------------

SUBROUTINE CV_SQNORM1(YDGEOMETRY,YDCV,PNORM)
USE GEOMETRY_MOD , ONLY : GEOMETRY
TYPE(GEOMETRY)       ,INTENT(IN)    :: YDGEOMETRY
TYPE (CONTROL_VECTOR),INTENT(IN)    :: YDCV
REAL(KIND=JPRB)      ,INTENT(OUT)   :: PNORM
  call abor1("cv_sqnorm1 should never be called")

END SUBROUTINE CV_SQNORM1

! ------------------------------------------------------------------

SUBROUTINE CV_SQNORM4(YDGEOMETRY,YDCV,PNORM1,PNORM2,PNORM3,PNORM4,PNORM5,PNORM6,PNORM7)
USE GEOMETRY_MOD , ONLY : GEOMETRY
TYPE(GEOMETRY)       ,INTENT(IN)    :: YDGEOMETRY
TYPE (CONTROL_VECTOR),INTENT(IN)    :: YDCV
REAL(KIND=JPRB)      ,INTENT(OUT)   :: PNORM1,PNORM2,PNORM3,PNORM4,PNORM5,PNORM6,PNORM7
  call abor1("cv_sqnorm4 should never be called")

END SUBROUTINE CV_SQNORM4

! ------------------------------------------------------------------

REAL(KIND=JPRB) FUNCTION CTLVEC_NORM(YDGEOMETRY,YDCV)
USE GEOMETRY_MOD , ONLY : GEOMETRY
TYPE(GEOMETRY)       , INTENT(IN) :: YDGEOMETRY
TYPE (CONTROL_VECTOR), INTENT(IN) :: YDCV
  call abor1("ctlvec_norm should never be called")

END FUNCTION CTLVEC_NORM

! ------------------------------------------------------------------

REAL(KIND=JPRB) FUNCTION MAXVAL_CV(YDCV)
TYPE (CONTROL_VECTOR), INTENT(IN) :: YDCV
  call abor1("maxval_cv should never be called")

END FUNCTION MAXVAL_CV

! ------------------------------------------------------------------

REAL(KIND=JPRB) FUNCTION MINVAL_CV(YDCV)
TYPE (CONTROL_VECTOR), INTENT(IN) :: YDCV
  call abor1("minval_cv should never be called")

END FUNCTION MINVAL_CV

! ------------------------------------------------------------------

REAL(KIND=JPRB) FUNCTION SUM_CV(YDCV)
TYPE (CONTROL_VECTOR), INTENT(IN) :: YDCV
  call abor1("sum_cv should never be called")

END FUNCTION SUM_CV

! ------------------------------------------------------------------

SUBROUTINE DPLENS(YDCV,KLEN,KLENG)

! Adjusts length of control variable used in dot-product to avoid
! double counting for the parts of the contol variable that is replicated
! on different mpi tasks (currently VarBC part and in a more complex way
! the mean winds for Aladin)

TYPE (CONTROL_VECTOR), INTENT(IN) :: YDCV ! Control variable
INTEGER(KIND=JPIM), INTENT(OUT)   :: KLEN   ! Adjusted local lenght of ctv
INTEGER(KIND=JPIM), INTENT(OUT)   :: KLENG  ! Adjusted global lenght of ctv
  call abor1("dplens should never be called")

END SUBROUTINE DPLENS

! ------------------------------------------------------------------

END MODULE CONTROL_VECTORS_PARA_MOD
