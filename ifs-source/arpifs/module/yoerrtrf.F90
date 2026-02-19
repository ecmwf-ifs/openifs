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

MODULE YOERRTRF

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTRF* - RRTM REFERENCE ATMOSPHERE
!     -----------------------------------------------------------------

REAL(KIND=JPRB) , DIMENSION(59) :: PREF
REAL(KIND=JPRB) , DIMENSION(59) :: PREFLOG
REAL(KIND=JPRB) , DIMENSION(59) :: TREF
REAL(KIND=JPRB)  :: CHI_MLS(7,59)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! PREF   :  REAL    
! PREFLOG: REAL
! TREF   : REAL
!     -----------------------------------------------------------------
END MODULE YOERRTRF
