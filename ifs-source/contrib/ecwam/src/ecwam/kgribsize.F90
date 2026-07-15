! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE KGRIBSIZE(IU06, KLEN, ISIZE, CALLINGSUB) 

! ----------------------------------------------------------------------

!**** *KGRIBSIZE* -

!     J. BIDLOT     ECMWF  JULY 2000 


!*    PURPOSE.
!     --------
!     UTILITY USED TO DETERMINE APOSTERIORI THE MINIMUM SIZE OF GRIB
!     INPUT DATA AND TO RESET THE FILE POINTER TO THE BEGINNING OF
!     THE FILE.

!**   INTERFACE.
!     ----------
!     *CALL KGRIBSIZE(IU06, KLEN, ISIZE, CALLINGSUB)*

!     METHOD.
!     -------
!     EXTERNALS.
!     ----------
!     REFERENCE.
!     ----------
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWGRIB, ONLY : JPKSIZE_T


      IMPLICIT NONE

      INTEGER(KIND=JWIM) ::  IU06
      INTEGER(KIND=JWIM) ::  ISIZE

      INTEGER(KIND=JPKSIZE_T) ::  KLEN

      CHARACTER (LEN=*) :: CALLINGSUB

      WRITE(IU06,*) ' '
      WRITE(IU06,*) ' ***** WARNING IN ',CALLINGSUB
      WRITE(IU06,*) ' SIZE  OF THE GRIB BUFFER IS NOT BIG ENOUGH.'
      WRITE(IU06,*) ' IT WAS ', ISIZE

      ISIZE=(KLEN+KIND(ISIZE)-1)/KIND(ISIZE)+1

      WRITE(IU06,*) ' IT SHOULD AT LEAST BE ', ISIZE
      WRITE(IU06,*) ' THE SIZE WAS RESET AUTOMATICALLY'
      WRITE(IU06,*) ' AND THE FIELD READ WITH THE NEW SIZE'
      WRITE(IU06,*) ' IF THIS PROBLEM OCCURS TOO OFTEN'
      WRITE(IU06,*) ' MODIFY THE VALUE OF NBIT IN SOURCE'
      WRITE(IU06,*) ' ***** WARNING ****** WARNING *****'
      WRITE(IU06,*) ' '
      CALL FLUSH(IU06)

      END SUBROUTINE KGRIBSIZE
