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

SUBROUTINE FPNEAR(KASLB1,KSLWIDE,KFIELDS,KMASKS,KGPST,KGPEND,KFPROMA,KFLDBUF, &
                & KMASK,KS0,LDMASK,PBUF,PMASK,PUNDEF,PROW)

!**** *FPNEAR*  - Interpolate field for Fullpos using the nearest point in the
!                 halo

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE* 
!      Original : 01-10-2014

!     Description.
!     ------------
!     This routine performs horizontal interpolations using the nearest point 
!     of the same kind in the Fullpos halo.
!     Arguments are the same as those of FPAVG.



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
INTEGER(KIND=JPIM) :: ISLW, JOFF
INTEGER(KIND=JPIM) :: JIND (KFPROMA,2*KSLWIDE*2*KSLWIDE)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FPNEAR',0,ZHOOK_HANDLE)

ISLW = KSLWIDE

! Compute point indexes, nearest first; we make larger 
! and larger squares around the target point (JOFF increases).

DO JI = KGPST, KGPEND

  II = 0

  DO JOFF = 1, ISLW

    JLA = ISLW-JOFF+1
    DO JLO = ISLW-JOFF+1, ISLW+JOFF-1
      II = II + 1
      JIND (JI,II) = KS0(JI,JLA)+(JLO-ISLW)
    ENDDO

    JLO = ISLW+JOFF
    DO JLA = ISLW-JOFF+1, ISLW+JOFF-1
      II = II + 1
      JIND (JI,II) = KS0(JI,JLA)+(JLO-ISLW)
    ENDDO

    JLA = ISLW+JOFF
    DO JLO = ISLW+JOFF, ISLW-JOFF+2, -1
      II = II + 1
      JIND (JI,II) = KS0(JI,JLA)+(JLO-ISLW)
    ENDDO
   
    JLO = ISLW-JOFF+1
    DO JLA = ISLW+JOFF, ISLW-JOFF+2, -1
      II = II + 1
      JIND (JI,II) = KS0(JI,JLA)+(JLO-ISLW)
    ENDDO

  ENDDO

  IF (II /= 2*ISLW*2*ISLW) &
  & CALL ABOR1 ('FPNEAR: COUNT MISMATCH')

ENDDO


DO JF = 1, KFIELDS

  IF (.NOT.LDMASK(JF)) THEN

    IADD = KASLB1 * (JF-1)

    DO JI = KGPST, KGPEND

      PROW(JI,JF) = PUNDEF

! Choose nearest point

      IF (PMASK (JI, KMASK (JF)) > 0._JPRB) THEN

        DO II = 1, 2*ISLW*2*ISLW
          IF (PBUF (JIND (JI,II)+IADD) /= PUNDEF) THEN
            PROW(JI,JF) = PBUF (JIND (JI,II)+IADD)
            EXIT
          ENDIF
        ENDDO

      ENDIF

    ENDDO

  ENDIF

ENDDO

IF (LHOOK) CALL DR_HOOK('FPNEAR',1,ZHOOK_HANDLE)

END SUBROUTINE FPNEAR

