! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE PHYS_ARRAYS_FIN(DIMS, PAUX, ZDDHS, ZDIAG, PSURF, PRAD, FLUX)

!**** *Initialization of derived variables local to physics sub-space

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------

! Derived arrays           Reserved space
! -----------------------------------------------
! PAUX                     PAUX_ARR
! ZDDHS                    PDDHS_ARR
! ZDIAG                    PDIAG_ARR
! PSURF                    PSURF_ARR
!                          KSURF_ARR
! PRAD                     PRAD_ARR
! FLUX                     PFLUX_ARR

!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!      Original : 31-Jan-2013  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!     P. Lopez, ECMWF, July 2015 Added lightning fields in ZDIAG.
!     B. Ingleby,      2019-01-17 Remove PQCFL

!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER ,ONLY : AUX_RAD_TYPE, SURF_AND_MORE_TYPE, AUX_TYPE, &
  &  AUX_DIAG_TYPE, FLUX_TYPE, DDH_SURF_TYPE, DIMENSION_TYPE


!     -----------------------------------------------------------------------------

IMPLICIT NONE

TYPE (DIMENSION_TYPE)    ,INTENT(IN)    :: DIMS
TYPE (AUX_TYPE)          ,INTENT(INOUT) :: PAUX
TYPE (DDH_SURF_TYPE)     ,INTENT(INOUT) :: ZDDHS
TYPE (AUX_DIAG_TYPE)     ,INTENT(INOUT) :: ZDIAG
TYPE (SURF_AND_MORE_TYPE),INTENT(INOUT) :: PSURF
TYPE (AUX_RAD_TYPE)      ,INTENT(INOUT) :: PRAD
TYPE (FLUX_TYPE)         ,INTENT(INOUT) :: FLUX

!-----------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PHYS_ARRAYS_FIN',0,ZHOOK_HANDLE)

!-----------------------------------------------------------------------

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PHYS_ARRAYS_FIN',1,ZHOOK_HANDLE)
END SUBROUTINE PHYS_ARRAYS_FIN
