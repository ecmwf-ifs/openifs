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

!**** PRISM_DUMMY_MOD.F90
!
!     Purpose.
!     --------
!
!     
!     Contains   
!     ----------
!     
!      MODULE PRISM_CONSTANTS 
!
!     Reference:
!     ---------
!       S. Valcke, 2006: OASIS4 User Guide  
!       PRISM Support Initiative Report No 3,
!       CERFACS, Toulouse, France, 64 pp.
!
!     Author:
!     -------
!       Johannes Flemming 
!
!     Modifications.
!     --------------
!      F. Vana  05-Mar-2015  Support for single precision

!**************************************************************************

MODULE YOMPRISM

USE PARKIND1 , ONLY : JPRD, JPIM
IMPLICIT NONE
SAVE
  
   TYPE PRISM_TIME_STRUCT
      REAL(KIND=JPRD)   :: SECOND
      INTEGER(KIND=JPIM)            :: MINUTE
      INTEGER(KIND=JPIM)            :: HOUR
      INTEGER(KIND=JPIM)            :: DAY
      INTEGER(KIND=JPIM)            :: MONTH
      INTEGER(KIND=JPIM)            :: YEAR
   END TYPE PRISM_TIME_STRUCT

 
   TYPE(PRISM_TIME_STRUCT) :: PRISM_INITIAL_DATE
   TYPE(PRISM_TIME_STRUCT) :: PRISM_JOBSTART_DATE
   TYPE(PRISM_TIME_STRUCT) :: PRISM_JOBEND_DATE

   INTEGER(KIND=JPIM), PARAMETER :: PRISM_GRIDLESS                = 8
   INTEGER(KIND=JPIM), PARAMETER :: PRISM_GAUSSREDUCED_REGVRT   = 9
   INTEGER(KIND=JPIM), PARAMETER :: PRISM_DOUBLE_PRECISION =  5
   INTEGER(KIND=JPIM), PARAMETER :: PRISM_UNDEFINED = -65535

END MODULE YOMPRISM

