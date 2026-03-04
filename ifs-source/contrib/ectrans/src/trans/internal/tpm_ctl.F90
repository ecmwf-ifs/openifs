! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_CTL

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE, INTRINSIC :: iso_c_binding, ONLY:  C_PTR, C_NULL_PTR
USE SHAREDMEM_MOD ,ONLY : SHAREDMEM
IMPLICIT NONE

SAVE


TYPE CTL_TYPE

LOGICAL :: LREAD_LEGPOL = .FALSE.
LOGICAL :: LWRITE_LEGPOL = .FALSE.
CHARACTER(LEN=256) :: CLEGPOLFNAME='legpol_file'
CHARACTER(LEN=4) :: CIO_TYPE='file'
TYPE(SHAREDMEM) :: STORAGE

END TYPE CTL_TYPE


TYPE(CTL_TYPE),ALLOCATABLE,TARGET :: CTL_RESOL(:)
TYPE(CTL_TYPE),POINTER     :: C


END MODULE TPM_CTL







