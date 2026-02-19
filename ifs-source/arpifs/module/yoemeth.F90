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

MODULE YOEMETH

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!     ** YOEMETH - CONTROL PARAMETERS FOR METHANE OXIDATION
!     -----------------------------------------------------------------

!     * E.C.M.W.F. PHYSICS PACKAGE *

!     C. JAKOB        E.C.M.W.F.          98/04/07

!      NAME     TYPE      PURPOSE
!      ----     ----      -------

!     *RALPHA1* REAL      CONSTANT IN TIME SCALE 1 CALCULATIONS
!     *RALPHA2* REAL      CONSTANT IN TIME SCALE 2 CALCULATIONS
!     *RQLIM*   REAL      UPPERLIMIT OD SPECIFIC HUMIDITY FOR 
!                         METHANE OXIDATION TO BE ACTIVE
!     *RPBOTOX* REAL      PRESSURE BELOW WHICH METHANE OXIDATION
!                         IS CONSIDERED ACTIVE
!     *RPBOTPH* REAL      PRESSURE BELOW WHICH H2O PHOTOLYSIS
!                         IS CONSIDERED ACTIVE
!     *RPTOPOX* REAL      PRESSURE BELOW WHICH SHORTEST TIME SCALE
!                         IS USED IN METHANE OXIDATION
!     *RPTOPPH* REAL      PRESSURE BELOW WHICH SHORTEST TIME SCALE
!                         IS USED IN H2O PHOTOLYISIS
!     *RALPHA3* REAL      CONSTANT IN TIME SCALE 2 CALCULATIONS
!     *RLOGPPH* REAL      CONSTANT IN TIME SCALE 2 CALCULATIONS

REAL(KIND=JPRB) :: RALPHA1
REAL(KIND=JPRB) :: RALPHA2
REAL(KIND=JPRB) :: RQLIM
REAL(KIND=JPRB) :: RPBOTOX
REAL(KIND=JPRB) :: RPBOTPH
REAL(KIND=JPRB) :: RPTOPOX
REAL(KIND=JPRB) :: RPTOPPH
REAL(KIND=JPRB) :: RALPHA3
REAL(KIND=JPRB) :: RLOGPPH

END MODULE YOEMETH
