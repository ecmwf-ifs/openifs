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

MODULE YEMMP

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!  ---------------------------------------------------------------------
!  C. Fischer : 99-06-11

!* store extra setup to perform distributed work along the efficient
!* wavenumber kstar in aladin 3d-var 

!  NEPROCN: number of the processor which possesses global kstar
!  NUEMP  : number of kstars treated in this processor
!  MYENS  : array containing those kstars (in 0:nsmax) treated in this proc
!  NUEMPP : array containing the number of kstars treated proc by proc
!  NEALLNS: array containing in a row the values of kstars proc by proc
!  NEPTRNS: pointer into neallns


TYPE :: TEMMP
  INTEGER(KIND=JPIM), POINTER :: NEPROCN(:) => NULL()
  INTEGER(KIND=JPIM) :: NUEMP
  INTEGER(KIND=JPIM), POINTER :: MYENS(:) => NULL()
  INTEGER(KIND=JPIM), POINTER :: NUEMPP(:) => NULL()
  INTEGER(KIND=JPIM), POINTER :: NEALLNS(:) => NULL()
  INTEGER(KIND=JPIM), POINTER :: NEPTRNS(:) => NULL()
END TYPE TEMMP

!  ---------------------------------------------------------------------
END MODULE YEMMP
