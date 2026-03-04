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

MODULE YOEOZOC

USE PARKIND1, ONLY : JPRB

IMPLICIT NONE

!     ------------------------------------------------------------------
!*     *YOEOZOC* SPECTRAL DISTRIBUTION OF OZONE
!     ------------------------------------------------------------------

TYPE :: TEOZOC

REAL(KIND=JPRB),ALLOCATABLE :: ZOZCL(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZOZCLAQUA(:,:)

REAL(KIND=JPRB) :: COZQC(21) = -HUGE(1._JPRB)
REAL(KIND=JPRB) :: COZQS(15) = -HUGE(1._JPRB)
REAL(KIND=JPRB) :: COZHC(21) = -HUGE(1._JPRB)
REAL(KIND=JPRB) :: COZHS(15) = -HUGE(1._JPRB)

REAL(KIND=JPRB) :: RSINC(64)     = -HUGE(1._JPRB)
REAL(KIND=JPRB) :: ROZT(64,0:60) = -HUGE(1._JPRB)
REAL(KIND=JPRB) :: RPROC(0:60)   = -HUGE(1._JPRB)

END TYPE TEOZOC

!*     *YOEOZOC* SPECTRAL DISTRIBUTION OF OZONE
!                     (TRIANGULAR *T5* TRUNCATION FOR OZONE).

!     R.G AND M.J        E.C.M.W.F.     29/11/82.

!     *Modifications*
!     O.Marsden Jan 2017 - define a type containing these (time-variable) arrays

!      NAME     TYPE      PURPOSE
!      ----     ----      -------

!     *COZ__*   REAL      *REFERS TO *OZONE.
!     *C___C*   REAL      *REFERS TO *COS COMPONENT.
!     *C___S*   REAL      *REFERS TO *SIN COMPONENT.
!       ----------------------------------------------------------------
END MODULE YOEOZOC
