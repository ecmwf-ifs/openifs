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

MODULE CPLNG_TYPES_MOD

   USE PARKIND1, ONLY: JPIM
   USE PARKIND_OCEAN, ONLY : JPRO

   IMPLICIT NONE

   PRIVATE

   PUBLIC CPLNG_FLD_TYPE

   TYPE CPLNG_FLD_TYPE
      INTEGER(KIND=JPIM),ALLOCATABLE :: ID(:,:)
      CHARACTER(LEN=128)             :: NAME
      INTEGER(KIND=JPIM)             :: TYPE
      INTEGER(KIND=JPIM)             :: INOUT
      INTEGER(KIND=JPIM)             :: STAGE
      INTEGER(KIND=JPIM)             :: NUM_LVL
      INTEGER(KIND=JPIM)             :: NUM_CAT
      REAL(KIND=JPRO),   ALLOCATABLE :: D(:,:,:)
   END TYPE CPLNG_FLD_TYPE

END MODULE CPLNG_TYPES_MOD
