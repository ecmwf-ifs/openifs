! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWTEXT

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *TEXT* - FILE NAME INFORMATION.

      INTEGER(KIND=JWIM) :: ICPLEN 
      CHARACTER(LEN=  3) :: USERID
      CHARACTER(LEN=  3) :: RUNID
      CHARACTER(LEN= 60) :: PATH
      CHARACTER(LEN=256) :: CPATH
      LOGICAL            :: LRESTARTED 
      CHARACTER(LEN=263) :: CWI

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *ICPLEN*    INTEGER   LENGTH OF CPATH.
!      *USERID*    CHARACTER USERID FOR ECFS FILE NAMES.
!      *RUNID*     CHARACTER RUN IDENTIFIER FOR ECFS FILE NAMES.
!      *PATH*      CHARACTER ECFS PATH NAME FOR FILES.
!      *CPATH*     CHARACTER PATH FOR OUTPUT TO DISK.
!      *LRESTARTED*CHARACTER INDICATED WHETHER RUN IS IN RESTARTED MODE.
!      *CWI*       CHARACTER PATH AND NAME OF THE RESTART INFORMATION.
!                            FILE.
! ----------------------------------------------------------------------
      END MODULE YOWTEXT
