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

MODULE YOESRTAER

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOESRTAER* - AEROSOL OPTICAL PROPERTIES FOR SRTM 
!     ------------------------------------------------------------------

!-- AEROSOL OPTICAL PROPERTIES

!-- coefficients and parameters related to SRTM_SW 16 sp.int.

REAL(KIND=JPRB) :: RSRTAUA(14,6)
REAL(KIND=JPRB) :: RSRPIZA(14,6)
REAL(KIND=JPRB) :: RSRASYA(14,6)

!     -----------------------------------------------------------------

!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      03/03/07

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : -------
!RSRTAUA :  REAL     S.W. NORMALIZED OPTICAL THICKNESS AT 0.55 MICRON
!RSRPIZA :  REAL     S.W. SINGLE SCATTERING ALBEDO
!RSRASYA :  REAL     S.W. ASSYMETRY FACTOR

!     -----------------------------------------------------------------
END MODULE YOESRTAER

