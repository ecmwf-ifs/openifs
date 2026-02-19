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

MODULE LFIMOD
! Jan-2011 P. Marguinaud Interface to thread-safe LFI
! Sep-2012 P. Marguinaud Initialize data + DrHook
USE PARKIND1, ONLY : JPIM, JPRB, JPIB, JPIA
USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK
USE LFI_PRECISION
IMPLICIT NONE

      TYPE LFICOM
      SEQUENCE
      CHARACTER (LEN=8)   :: CMAGIC = "LFI_FORT"
      INTEGER (KIND=JPIA) :: ILFICC = 0_JPIA
      END TYPE LFICOM


      TYPE (LFICOM), SAVE, TARGET :: LFICOM_DEFAULT
      LOGICAL, SAVE :: LFICOM_DEFAULT_INIT = .FALSE.

      CONTAINS

      SUBROUTINE NEW_LFI_DEFAULT ()
      INTEGER :: IERR
      REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK ('LFICOM:NEW_LFI_DEFAULT',0,ZHOOK_HANDLE)

      IF (.NOT. LFICOM_DEFAULT_INIT) THEN
        CALL NEW_LFI (LFICOM_DEFAULT, IERR)
        LFICOM_DEFAULT_INIT = .TRUE.
      ENDIF

      IF (LHOOK) CALL DR_HOOK ('LFICOM:NEW_LFI_DEFAULT',1,ZHOOK_HANDLE)

      END SUBROUTINE NEW_LFI_DEFAULT

      SUBROUTINE NEW_LFI (LFI, KERR, KPNXFI, KPFACX)
      TYPE (LFICOM) :: LFI
      INTEGER, INTENT(OUT) :: KERR
      INTEGER, OPTIONAL, INTENT(IN) :: KPNXFI
      INTEGER, OPTIONAL, INTENT(IN) :: KPFACX
      REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK ('LFICOM:NEW_LFI',0,ZHOOK_HANDLE)

      KERR = 0

      IF (LHOOK) CALL DR_HOOK ('LFICOM:NEW_LFI',1,ZHOOK_HANDLE)

      END SUBROUTINE NEW_LFI

      SUBROUTINE FREE_LFI (LFI, KERR)
      TYPE (LFICOM) :: LFI
      INTEGER, INTENT(OUT) :: KERR
      REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK ('LFICOM:FREE_LFI',0,ZHOOK_HANDLE)

      KERR = 0

      IF (LFI%ILFICC /= 0) CALL LFI_HNDL_FREE (LFI)

      IF (LHOOK) CALL DR_HOOK ('LFICOM:FREE_LFI',1,ZHOOK_HANDLE)

      END SUBROUTINE FREE_LFI

END MODULE LFIMOD


