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

SUBROUTINE SPEC_SPLIT(YDSPIN,YDSP1,YDSP2)

!   Purpose.
!   --------
!     Split spectral fields into two spectral fields.

!   Author.
!   -------
!     Y. Tremolet

!   Modifications.
!   --------------
!     Original   02-Aug-2006
! ------------------------------------------------------------------

USE PARKIND1, ONLY: JPIM, JPRB
USE YOMHOOK , ONLY: LHOOK, DR_HOOK, JPHOOK
USE YOMCT0  , ONLY: LELAM
USE INDEXFIND_MOD, ONLY: INDXFIND
USE SPECTRAL_FIELDS_MOD, ONLY : ASSIGNMENT(=), SPECTRAL_FIELD

IMPLICIT NONE
TYPE(SPECTRAL_FIELD), INTENT(IN)    :: YDSPIN
TYPE(SPECTRAL_FIELD), INTENT(INOUT) :: YDSP1
TYPE(SPECTRAL_FIELD), INTENT(INOUT) :: YDSP2

LOGICAL, ALLOCATABLE :: LL2D(:), LL3D(:)
INTEGER(KIND=JPIM) :: JFO,JFI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

! ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPEC_SPLIT',0,ZHOOK_HANDLE)

IF (YDSP1%NS1D+YDSP2%NS1D/=YDSPIN%NS1D) CALL ABOR1('SPEC_SPLIT: Error 1D')
IF (YDSP1%NS2D+YDSP2%NS2D/=YDSPIN%NS2D) CALL ABOR1('SPEC_SPLIT: Error 2D')
IF (YDSP1%NS3D+YDSP2%NS3D/=YDSPIN%NS3D) CALL ABOR1('SPEC_SPLIT: Error 3D')
IF (YDSP1%NSPEC2/=YDSPIN%NSPEC2 .OR. YDSP2%NSPEC2/=YDSPIN%NSPEC2) &
 & CALL ABOR1('SPEC_SPLIT: Error resolution')

ALLOCATE(LL2D(YDSPIN%NS2D))
ALLOCATE(LL3D(YDSPIN%NS3D))
LL2D(:)=.FALSE.
LL3D(:)=.FALSE.

! 2D fields
DO JFO=1,YDSP1%NS2D
  JFI=INDXFIND(YDSPIN%NGRIB2,YDSP1%NGRIB2(JFO))
  IF (JFI==0) CALL ABOR1('SPEC_SPLIT: 2D field not found')
  IF (LL2D(JFI)) CALL ABOR1('SPEC_CONCAT: duplicate 2D field')
  YDSP1%SP2D(:,JFO)=YDSPIN%SP2D(:,JFI)
  LL2D(JFI)=.TRUE.
ENDDO

DO JFO=1,YDSP2%NS2D
  JFI=INDXFIND(YDSPIN%NGRIB2,YDSP2%NGRIB2(JFO))
  IF (JFI==0) CALL ABOR1('SPEC_SPLIT: 2D field not found')
  IF (LL2D(JFI)) CALL ABOR1('SPEC_CONCAT: duplicate 2D field')
  YDSP2%SP2D(:,JFO)=YDSPIN%SP2D(:,JFI)
  LL2D(JFI)=.TRUE.
ENDDO

! 3D fields
DO JFO=1,YDSP1%NS3D
  JFI=INDXFIND(YDSPIN%NGRIB3,YDSP1%NGRIB3(JFO))
  IF (JFI==0) CALL ABOR1('SPEC_SPLIT: 3D field not found')
  IF (LL3D(JFI)) CALL ABOR1('SPEC_CONCAT: duplicate 3D field')
  YDSP1%SP3D(:,:,JFO)=YDSPIN%SP3D(:,:,JFI)
  LL3D(JFI)=.TRUE.
ENDDO 

DO JFO=1,YDSP2%NS3D
  JFI=INDXFIND(YDSPIN%NGRIB3,YDSP2%NGRIB3(JFO))
  IF (JFI==0) CALL ABOR1('SPEC_SPLIT: 3D field not found')
  IF (LL3D(JFI)) CALL ABOR1('SPEC_CONCAT: duplicate 3D field')
  YDSP2%SP3D(:,:,JFO)=YDSPIN%SP3D(:,:,JFI)
  LL3D(JFI)=.TRUE.
ENDDO 

! 1D fields
IF (LELAM) THEN
  IF     (YDSP1%NS1D==YDSPIN%NS1D) THEN
    YDSP1%SP1D(:,:)=YDSPIN%SP1D(:,:)
  ELSEIF (YDSP2%NS1D==YDSPIN%NS1D) THEN
    YDSP2%SP1D(:,:)=YDSPIN%SP1D(:,:)
  ELSE
    CALL ABOR1 ('SPEC_SPLIT: Error LELAM')
  ENDIF
ENDIF

DEALLOCATE(LL2D,LL3D)

IF (LHOOK) CALL DR_HOOK('SPEC_SPLIT',1,ZHOOK_HANDLE)
END SUBROUTINE SPEC_SPLIT
