! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPHALO(YDFPSTRUCT,KPROMA,KFIELDS,KSTGLO,KEND,KFLDCORE,PCORE,PHALO,KFLDPTP)

!**** *FPHALO*  - FULLPOS HALO FOR HORIZONTAL INTERPOLATIONS

!     PURPOSE.
!     --------
!        To fill the distributed halo from core data, prior to horizontal interpolations

!**   INTERFACE.
!     ----------
!       *CALL* *FPHALO*

!        EXPLICIT ARGUMENTS
!        --------------------

!        YDFPSTRUCT    - halo structure control
!        KPROMA  - horizontal segmentation
!        KFIELDS - number of fields
!        KSTGLO  - global first adress of points
!        KEND    - last active point in row
!        KFLDCORE- number of fields in core data
!        PCORE   - data core
!        PHALO   - data halo
!        KFLDPTP - field pointer for selection

!        IMPLICIT ARGUMENTS
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 21-Sep-2017 from fpmodprec

!     MODIFICATIONS.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE EINT_MOD, ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE    (SL_STRUCT), INTENT (IN)  :: YDFPSTRUCT
INTEGER (KIND=JPIM), INTENT (IN)  :: KPROMA
INTEGER (KIND=JPIM), INTENT (IN)  :: KFIELDS 
INTEGER (KIND=JPIM), INTENT (IN)  :: KSTGLO
INTEGER (KIND=JPIM), INTENT (IN)  :: KEND
INTEGER (KIND=JPIM), INTENT (IN)  :: KFLDCORE
REAL    (KIND=JPRB), INTENT (IN)  :: PCORE(KPROMA,KFLDCORE)
REAL    (KIND=JPRB), INTENT (INOUT) :: PHALO(YDFPSTRUCT%NASLB1,KFIELDS) ! inout because this subroutine fills only a part of it
INTEGER (KIND=JPIM), INTENT (IN), OPTIONAL  :: KFLDPTP(KFIELDS)

!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JF, JROF, IFLDPTP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FPHALO',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1. TRANSFER FROM CORE TO HALO
!           --------------------------

IF (PRESENT(KFLDPTP)) THEN
  IF (MINVAL(KFLDPTP) < 1 .OR. MAXVAL(KFLDPTP) > KFLDCORE) THEN
    CALL ABOR1('FPHALO : ERROR KFLDPTP VS KFLDCORE')
  ENDIF
  IF (YDFPSTRUCT%NSLWIDE > 0) THEN
    DO JF=1,KFIELDS
      IFLDPTP=KFLDPTP(JF)
!DIR$ NEXTSCALAR
!DIR$ IVDEP
!DIR$ NOINTERCHANGE
      DO JROF=1,KEND
        PHALO(YDFPSTRUCT%NSLCORE(JROF+KSTGLO-1),JF)=PCORE(JROF,IFLDPTP)
      ENDDO
    ENDDO
  ELSE
    ! Configuration of sampling ...
    DO JF=1,KFIELDS
      IFLDPTP=KFLDPTP(JF)
      PHALO(KSTGLO:KSTGLO-1+KEND,JF)=PCORE(1:KEND,IFLDPTP)
    ENDDO
  ENDIF
ELSE
  IF (YDFPSTRUCT%NSLWIDE > 0) THEN
    DO JF=1,KFIELDS
!DIR$ NEXTSCALAR
!DIR$ IVDEP
!DIR$ NOINTERCHANGE
      DO JROF=1,KEND
        PHALO(YDFPSTRUCT%NSLCORE(JROF+KSTGLO-1),JF)=PCORE(JROF,JF)
      ENDDO
    ENDDO
  ELSE
    ! Configuration of sampling ...
    DO JF=1,KFIELDS
      PHALO(KSTGLO:KSTGLO-1+KEND,JF)=PCORE(1:KEND,JF)
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FPHALO',1,ZHOOK_HANDLE)
END SUBROUTINE FPHALO
