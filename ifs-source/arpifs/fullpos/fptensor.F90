! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPTENSOR(KST,KEND,KPROMA,KFIELDS,KORDER,PGM,PBUF,PFPARITY)

!**** *FPTENSOR*  - FULLPOS TENSOR ASPECTS PRIOR TO HORIZONTAL POST-PROCESSING

!     PURPOSE.
!     --------
!        To compute the fields parity array and scale fields with map factor, depending
!        on the order of derivative of the fields.

!**   INTERFACE.
!     ----------
!       *CALL* *FPTENSOR*

!        EXPLICIT ARGUMENTS
!        --------------------
!         INPUT:
!          KST      : first address in horizontal dimension
!          KEND     : last  address in horizontal dimension
!          KPROMA   : horizontal dimension
!          KFIELDS  : number of fields
!          KORDER   : order of derivative of each field
!          PGM      : map factor

!         INPUT/OUTPUT:
!          PBUF     : interpolation buffer containing the fields to interpolate

!         OUTPUT:
!          PFPARITY : parity for interpolation buffer PBUF 1=scalar; -1=vector.

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
!        See documentation about FULL-POS.

!     AUTHOR.
!     -------
!        RYAD EL KHATIB *METEO-FRANCE*
!        ORIGINAL : 16-Sep-2016 from HPOS

!     MODIFICATIONS.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS
INTEGER(KIND=JPIM),INTENT(IN)    :: KORDER(KFIELDS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGM(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBUF(KPROMA,KFIELDS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPARITY(KFIELDS) 

REAL(KIND=JPRB) :: ZGM2(KPROMA)

INTEGER(KIND=JPIM) :: JFLD, JI

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPTENSOR',0,ZHOOK_HANDLE)

!*       1. PRELIMINARY INITIALISATIONS.
!           ----------------------------

IF (ANY(KORDER(:)==2)) THEN
  ZGM2(KST:KEND)=PGM(KST:KEND)**2
ENDIF

!*       2. TENSOR ASPECT
!           -------------

DO JFLD=1,KFIELDS

  IF (KORDER(JFLD) == 1) THEN
    PFPARITY(JFLD)=-1.0_JPRB
  ELSE
    PFPARITY(JFLD)=+1.0_JPRB
  ENDIF

  IF (KORDER(JFLD) > 0) THEN
    SELECT CASE (KORDER(JFLD))
    CASE (1) 
      ! Vectors => multiply by map factor if actual horizontal interpolations
      DO JI=KST,KEND
        PBUF(JI,JFLD)=PBUF(JI,JFLD)*PGM(JI)
      ENDDO
    CASE (2)
      ! 2nd order derivatives => multiply by square map factor 
      ! if actual horizontal interpolations
      DO JI=KST,KEND
        PBUF(JI,JFLD)=PBUF(JI,JFLD)*ZGM2(JI)
      ENDDO
    CASE (3:)
      DO JI=KST,KEND
        PBUF(JI,JFLD)=PBUF(JI,JFLD)*PGM(JI)**KORDER(JFLD)
      ENDDO
    END SELECT
  ENDIF 

ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPTENSOR',1,ZHOOK_HANDLE)
END SUBROUTINE FPTENSOR
