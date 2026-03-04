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

MODULE YOERRTAB

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

!    -------------------------------------------------------------------

REAL(KIND=JPRB) , DIMENSION(0:5000) :: TRANS
REAL(KIND=JPRB) :: BPADE

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! TRANS  :  REAL   : TABULATED REFERENCE FOR LW TRANSMISSION    
! BPADE  :  REAL   : INVERSE OF PADE APPROXIMATION CONSTANT  (= 1./0.278)     
!     -----------------------------------------------------------------
END MODULE YOERRTAB

