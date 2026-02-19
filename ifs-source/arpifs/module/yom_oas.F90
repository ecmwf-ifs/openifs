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

MODULE YOM_OAS

! -- 98-05-29       Author: S. Valcke

!  Contents : choice for the OASIS version: clim or pipe or sipc
!  --------
!  Modified : 00-02-18 JPh Piedelievre (MPI2)

IMPLICIT NONE
SAVE
CHARACTER(LEN=8), PARAMETER ::  CPHAN='SIPC' ! as $CHANNEL in namcouple
!     CHARACTER(len=3), PARAMETER ::  cpjobnam= 'IPC'

!     CHARACTER(len=8), PARAMETER ::  cphan='MPI2' ! as $CHANNEL in namcouple
CHARACTER(LEN=6), DIMENSION(2) :: CLBID

CHARACTER(LEN=6), PARAMETER ::  CPMODNAM= 'arpege'
! as $NBMODEL in namcouple
CHARACTER(LEN=3), PARAMETER ::  CPJOBNAM= 'CLI'
! as $JOBNAM in namcouple
CHARACTER(LEN=3), PARAMETER ::  CPMODINF= 'NOT'

END MODULE YOM_OAS
