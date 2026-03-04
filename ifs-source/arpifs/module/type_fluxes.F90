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

MODULE TYPE_FLUXES

!     Purpose.
!     --------
!      To define the descriptors of families of model fluxes.
!      "FLUXES_DESCRIPTOR" describes :
!      - the "FA" prefix name
!      - the "FA" suffix name
!      - the index values of full/half levels if the field is upper air
!      - the total number of levels
!      - the vectorial parity

!     Author.
!     -------
!        Ryad El Khatib *METEO-FRANCE*

!     Modifications.
!     --------------
!        Original : 2001-01-09
!        M.Hamrud   01-Oct-2003 CY28 Cleaning
!       R El Khatib 05-03-11 Cleanups

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB
USE PARDIM    ,ONLY : JPMXLE

IMPLICIT NONE
SAVE

PRIVATE
PUBLIC FLUXES_DESCRIPTOR, COPY_FLUX

TYPE FLUXES_DESCRIPTOR
CHARACTER(LEN=4)   :: CLPREF(JPMXLE+1)
CHARACTER(LEN=16)  :: CLSUFF
INTEGER(KIND=JPIM) :: INUMLEV(JPMXLE+1)
INTEGER(KIND=JPIM) :: ILEV
INTEGER(KIND=JPIM) :: IPARITY
END TYPE FLUXES_DESCRIPTOR

CONTAINS

SUBROUTINE COPY_FLUX(PFLUXOUT,PFLUXIN)

USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

REAL(KIND=JPRB), ALLOCATABLE, INTENT(INOUT) :: PFLUXOUT(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE, INTENT(IN)    :: PFLUXIN(:,:,:)

INTEGER(KIND=JPIM) :: JBL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('TYPE_FLUXES:COPY_FLUX',0,ZHOOK_HANDLE)

IF (ALLOCATED(PFLUXOUT) .AND. ALLOCATED(PFLUXIN)) THEN
  IF (ALL(SHAPE(PFLUXOUT) == SHAPE(PFLUXIN))) THEN
!$OMP PARALLEL DO PRIVATE(JBL)
    DO JBL = 1, SIZE(PFLUXIN,DIM=3)
      PFLUXOUT(:,:,JBL) = PFLUXIN(:,:,JBL)
    ENDDO
!$OMP END PARALLEL DO
  ELSE
    CALL ABOR1('TYPE_FLUXES:COPY_FLUX different shapes')
  ENDIF
ELSE
  CALL ABOR1('TYPE_FLUXES:COPY_FLUX trying to copy into unallocated array')
ENDIF

IF (LHOOK) CALL DR_HOOK('TYPE_FLUXES:COPY_FLUX',1,ZHOOK_HANDLE)

END SUBROUTINE COPY_FLUX
!     ------------------------------------------------------------
END MODULE TYPE_FLUXES
