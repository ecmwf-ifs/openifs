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

SUBROUTINE FPAVG(KASLB1,KSLWIDE,KFIELDS,KMASKS,KGPST,KGPEND,KFPROMA,KFLDBUF, &
               & KMASK,KS0,LDMASK,PBUF,PMASK,PUNDEF,PROW)

!**** *FPAVG*  - Fullpos interpolation based on average over the Fullpos halo

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE* 
!      Original : 01-10-2014

!     Description.
!     ------------
!     This routine makes interpolation by doing an average of valid (different 
!     from PUNDEF) values in the Fullpos halo. All arguments are very similar to
!     those of PINT4, except :
!     - KMASK : for each field, the rank of the mask to be used (the second
!     dimension of PMASKS)
!     - PMASKS : every possible mask
!     - KS0 : index in the halo, computed by SUHOX1, SUEHOX1, FPSCAX

!         KASLB1   : size of interpolation buffer (core + halo)
!         KSLWIDE  : width of halo part of interpolation buffer


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KASLB1
INTEGER(KIND=JPIM),INTENT(IN)    :: KSLWIDE
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS
INTEGER(KIND=JPIM),INTENT(IN)    :: KMASKS
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDBUF
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPST
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPEND
INTEGER(KIND=JPIM),INTENT(IN)    :: KMASK(KFIELDS)
INTEGER(KIND=JPIM),INTENT(IN)    :: KS0(KFPROMA,KSLWIDE*2)
LOGICAL           ,INTENT(IN)    :: LDMASK(KFIELDS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBUF(KASLB1*KFLDBUF)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMASK (KFPROMA,KMASKS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUNDEF
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PROW(KFPROMA,KFIELDS)

INTEGER(KIND=JPIM) :: IADD

INTEGER(KIND=JPIM) :: JF, JI, II, JLA, JLO
INTEGER(KIND=JPIM) :: JIND (KFPROMA,2*KSLWIDE*2*KSLWIDE)
REAL(KIND=JPRB)    :: ZVAL, ZWGT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FPAVG',0,ZHOOK_HANDLE)

! Compute all indices of the halo points, for each target point

DO JI = KGPST,KGPEND
  DO JLA = 1, 2*KSLWIDE
    DO JLO = 1, 2*KSLWIDE
      JIND (JI,(JLA-1)*2*KSLWIDE+JLO) = KS0(JI,JLA)+(JLO-KSLWIDE)
    ENDDO
  ENDDO
ENDDO

DO JF = 1, KFIELDS

  IF (.NOT.LDMASK(JF)) THEN

    IADD = KASLB1 * (JF-1)

    DO JI = KGPST, KGPEND

      PROW(JI,JF) = PUNDEF

      IF (PMASK (JI, KMASK (JF)) > 0._JPRB) THEN

        ZVAL = 0._JPRB
        ZWGT = 0._JPRB

! Compute the average over the halo, using only valid values

        DO II = 1, 2*KSLWIDE*2*KSLWIDE
          IF (PBUF (JIND (JI,II)+IADD) /= PUNDEF) THEN
            ZVAL = ZVAL + PBUF (JIND (JI,II)+IADD)
            ZWGT = ZWGT + 1._JPRB
          ENDIF
        ENDDO

        IF (ZWGT > 0) THEN
          PROW(JI,JF) = ZVAL / ZWGT
        ENDIF

      ENDIF

    ENDDO

  ENDIF

ENDDO

IF (LHOOK) CALL DR_HOOK('FPAVG',1,ZHOOK_HANDLE)

END SUBROUTINE FPAVG

