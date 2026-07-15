! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE MERGE_AER_SURFACE(YAER_RAD_DESC, KSTART, KEND, KCOL, &
     &                        KSOURCE, PNEWAERO_MMR, PAERO_MMR, LDCONTIG)

! MERGE_AER_SURFACE - merge prognostic (or parametric) aerosol species
!   at the surface into an array which may already contain
!   climatological species - this is needed for visibility calculations
!
! AUTHOR
! ------
!    Robin Hogan, ECMWF
!    Original: 2022-11-16
!
! MODIFICATIONS
! -------------
!

USE PARKIND1,      ONLY : JPIM, JPRB
USE YOMHOOK,       ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOEAERRADDESC, ONLY : TAER_RAD_DESC, IPARAMETRIC

IMPLICIT NONE

TYPE(TAER_RAD_DESC), INTENT(IN) :: YAER_RAD_DESC
! Indices of start and end columns
INTEGER(KIND=JPIM), INTENT(IN) :: KSTART, KEND
! Number of columns
INTEGER(KIND=JPIM), INTENT(IN) :: KCOL
! "Source" of new aerosols, either IPROGNOSTIC or ICLIMATOLOGICAL
INTEGER(KIND=JPIM), INTENT(IN) :: KSOURCE
! Mass mixing ratio of new aerosols to be inserted (column,newtype)
REAL(KIND=JPRB),    INTENT(IN) :: PNEWAERO_MMR(:,:)
! Existing mass mixing ratio with gaps for new species (column,type)
REAL(KIND=JPRB), INTENT(INOUT) :: PAERO_MMR(:,:)
! Are the aerosols in the new aerosol array contiguous, or arranged
! according to their expected location in a prognostic aerosol array
! (default)?
LOGICAL,  OPTIONAL, INTENT(IN) :: LDCONTIG

! Aerosol loop index, and index to the input aerosol species
INTEGER(KIND=JPIM) :: JAER, IAER

LOGICAL :: LLCONTIG

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('MERGE_AER_SURFACE',0,ZHOOK_HANDLE)

! Input data are assumed to be contiguous unless the source of new
! aerosols is prognostic in which case they're assumed to be
! non-contiguous (but in this case the LDCONTIG argument can state
! them to be contiguous).
IF (KSOURCE == IPARAMETRIC) THEN
  IF (PRESENT(LDCONTIG)) THEN
    LLCONTIG = LDCONTIG
  ELSE
    LLCONTIG = .FALSE.
  ENDIF
ELSE
  LLCONTIG = .TRUE.
ENDIF

IAER = 1
IF (LLCONTIG) THEN
  ! Contiguous case: loop through aerosol species in the output
  ! array...
  DO JAER = 1,YAER_RAD_DESC%NAEROSOL
    IF (YAER_RAD_DESC%DATA(JAER)%NSOURCE == KSOURCE) THEN
      ! ...and if it is to be filled then copy over the next aerosol
      ! in the input array
      PAERO_MMR(KSTART:KEND,JAER) = PNEWAERO_MMR(KSTART:KEND,IAER)
      IAER = IAER + 1
    ENDIF
  ENDDO
ELSE
  ! Non-contiguous case: loop through aerosol species in the output
  ! array...
  DO JAER = 1,YAER_RAD_DESC%NAEROSOL
    IF (YAER_RAD_DESC%DATA(JAER)%NSOURCE == KSOURCE) THEN
      ! ...and if it is to be filled copy over the prognostic aerosol
      ! at the appropriate location in the input array
      PAERO_MMR(KSTART:KEND,JAER) &
           &  = PNEWAERO_MMR(KSTART:KEND,YAER_RAD_DESC%PROG_INDEX(IAER))
      IAER = IAER + 1
    ENDIF
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('MERGE_AER_SURFACE',1,ZHOOK_HANDLE)
  
END SUBROUTINE MERGE_AER_SURFACE
