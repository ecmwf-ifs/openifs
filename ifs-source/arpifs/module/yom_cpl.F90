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

MODULE YOM_CPL

! --  1998-05-29
!  Modified : 00-02-18 JPh Piedelievre (MPI2)

!  Contents : variables for field symbolic names and restart files
!  --------

! -- cl_writ  : symbolic name for field to write

! -- cl_read  : symbolic name for field to READ

! -- cl_f_writ: restart file name for  field to WRITE

! -- cl_f_read: restart file name for  field to read

!     -------------------------------------------------------------------

USE PAR_COU   , ONLY : JPMAXFLD

IMPLICIT NONE
SAVE
CHARACTER(LEN=8), DIMENSION(JPMAXFLD) :: CL_WRIT, CL_READ
CHARACTER(LEN=8), DIMENSION(JPMAXFLD) :: CL_F_WRIT, CL_F_READ
!     -------------------------------------------------------------------
END MODULE YOM_CPL
