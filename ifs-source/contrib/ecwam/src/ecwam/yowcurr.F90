! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWCURR

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *CURRENT* - CURRENT FIELD.

      REAL(KIND=JWRB), PARAMETER :: CURRENT_MAX=1.5_JWRB
      REAL(KIND=JWRB), PARAMETER :: CURRENT_GRADIENT_MAX=0.00001_JWRB

      CHARACTER(LEN=14) :: CDTCUR
      CHARACTER(LEN=14) :: CDATECURA

      INTEGER(KIND=JWIM) :: IDELCUR

      LOGICAL :: LLCHKCFL
      LOGICAL :: LLCHKCFLA
      LOGICAL :: LLCFLCUROFF

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -----     --------
!      CURRENT_MAX REAL      MAXIMUM VALUE ALLOWED FOR EACH CURRENT COMPONENT (M/S)
!      CURRENT_GRADIENT_MAX  MAXIMUM VALUE FOR THE CURRENT GRADIENT (1/S)
!      *CDTCUR*    CHAR*14   REFERENCE DATE FOR U AN V.
!      *CDATECURA* CHAR*14   INITIAL DATE FOR INPUT OF U AND V.
!      *IDELCUR*   INTEGER   CURRENT INPUT STEP IN SECONDS.
!      *LLCHKCFL*  LOGICAL   TRUE IF THE CFL CRITERIA HAVE TO BE CHECKED 
!      *LLCHKCFLA* LOGICAL   TRUE IF THE CFL CRITERIA HAVE TO BE CHECKED
!                            AT ALL TIME STEPS (when currents are updated). 
!      *LLCFLCUROFF* LOGICAL IF TRUE THEN THE CURRENT REFRACTION TERMS WILL BE SET TO 0
!                            TO SATISFY THE CFL CRITERIA
!                            THIS SHOULD ONLY BE USED IN THE OPERATIONAL CONTEXT (to avoid crash)
!                            BUT NOT IN DEVELOPMENT MODE !!!!!!
!                            !!!!!!!!!!! ONLY CODED FOR IPROPAGS = 2 !!!!!!!!!!!!
! ----------------------------------------------------------------------
      END MODULE YOWCURR
