! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPCLIPHY(YDRQCLI,YDGFP,KFLD,KCOD,KST,KEND,KFPROMA,PCLI,PROW,LDCLI)

!**** *FPCLIPHY*  - FULL-POS writer of the output climatology

!     PURPOSE.
!     --------
!        To write the output climatology in the post-processing buffers
!        and initialize a logical climatology mask

!**   INTERFACE.
!     ----------
!       *CALL* *FPCLIPHY*

!        EXPLICIT ARGUMENTS
!        --------------------

!        INPUT:
!         YDRQCLI: fields descriptors in climatology row
!         KFLD   : number of fields in output row.
!         KST    : first output point in row.
!         KEND   : last output point in row.
!         KFPROMA: length of the output row.
!         PCLI   : output climatology row

!        OUTPUT:
!         PROW   : output fields row
!         LDCLI : logical mask : .TRUE. if the climatology
!                  field has been written

!        IMPLICIT ARGUMENTS
!        --------------------
!          See module above

!     METHOD.
!     -------
!          Scan the request ; if it is a climatology field
!          and the corresponding climatology is at disposal, 
!          then fill the row and set the logical mask to .TRUE.

!     EXTERNALS.
!     ----------
!      None

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        See documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 98-09-03 (partially from FPINTPHY)

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad: 31-01-2007 Cleanings + optional calcul. of PCLI,PROW,LDCLI
!      R. El Khatib 13-Dec-2012 Fullpos buffers reshaping
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMFP4L, ONLY : TRQFP, IFPSEARCH
USE YOMAFN   , ONLY : ALL_FPDSPHY_TYPES

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TRQFP),  INTENT(IN) :: YDRQCLI
TYPE(ALL_FPDSPHY_TYPES),  INTENT(IN) :: YDGFP
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLD 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCOD(KFLD) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPROMA 
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL :: PCLI(KFPROMA,0:YDRQCLI%NFIELDG) 
REAL(KIND=JPRB)   ,INTENT(OUT)  , OPTIONAL :: PROW(KFPROMA,KFLD) 
LOGICAL           ,INTENT(OUT)  , OPTIONAL :: LDCLI(KFLD) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JF, ISDOG, IACOT,IDPAT, ICOD, IPTR

LOGICAL :: LL, LLC, LLR

REAL(KIND=JPRB) :: ZU(KFPROMA), ZV(KFPROMA)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "cvlaniso.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FPCLIPHY',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1. Initializations
!           ---------------

! Compute LDCLI:
IF (PRESENT(LDCLI)) THEN
  LLC=.TRUE.
ELSE
  LLC=.FALSE.
ENDIF

! Compute PROW and use PCLI:
IF (PRESENT(PCLI).AND.(.NOT.PRESENT(PROW))) THEN
  CALL ABOR1(' FPCLIPHY: PCLI present => PROW should be present ')
ELSEIF (PRESENT(PROW).AND.(.NOT.PRESENT(PCLI))) THEN
  CALL ABOR1(' FPCLIPHY: PROW present => PCLI should be present ')
ELSEIF (PRESENT(PROW).AND.PRESENT(PCLI)) THEN
  LLR=.TRUE.
ELSE
  LLR=.FALSE.
ENDIF

ISDOG=IFPSEARCH(YDRQCLI,YDGFP%SDOG%ICOD)
IDPAT=IFPSEARCH(YDRQCLI,YDGFP%DPAT%ICOD)
IACOT=IFPSEARCH(YDRQCLI,YDGFP%ACOT%ICOD)

!     ------------------------------------------------------------------

!*       2. Loop on fields
!           --------------

DO JF=1,KFLD

  ICOD=KCOD(JF)
  IPTR=IFPSEARCH(YDRQCLI,ICOD)

  IF (IPTR > 0) THEN
    ! Climatology exists for this field :
    IF (ICOD == YDGFP%ST%ICOD .OR. ICOD == YDGFP%DST%ICOD .OR. &
     & ICOD == YDGFP%SD%ICOD) THEN  
      ! Actually Ts, Td, Snow depth are pronostic fields
      IF (LLC) LDCLI(JF)=.FALSE.
    ELSE
      ! Climatology exists, get it :
      IF (LLR) PROW(KST:KEND,JF)=PCLI(KST:KEND,IPTR)
      IF (LLC) LDCLI(JF)=.TRUE.
    ENDIF
  ELSE
    ! No actual climatology for this field :
    IF (IACOT*IDPAT*ISDOG > 0 .AND.  &
     & (ICOD == YDGFP%PADOU%ICOD .OR. ICOD == YDGFP%PADOV%ICOD)) THEN  
      ! Climatology exists indirectly for the vector of anisotropy,
      ! defined by its circular coordinates :
      ! module = (sigma**2)*(1-gamma)/(1+gamma) ; angle  = 2*alpha
      ! where sigma = standard deviation of orography
      ! gamma = Anisotropy coefficient of topography
      ! alpha = Direction of the principal axis of the topography
      IF (LLR) THEN
        LL=.TRUE.
        CALL CVLANISO(KST,KEND,KFPROMA,LL,LL,PCLI(1,IACOT),PCLI(1,IDPAT),ZU,ZV,PCLI(1,ISDOG))  
        IF (ICOD == YDGFP%PADOU%ICOD) THEN
          PROW(KST:KEND,JF)=ZU(KST:KEND)
        ELSEIF (ICOD == YDGFP%PADOV%ICOD) THEN
          PROW(KST:KEND,JF)=ZV(KST:KEND)
        ENDIF
      ENDIF
      IF (LLC) LDCLI(JF)=.TRUE.
    ELSE
      ! No climatology
      IF (LLC) LDCLI(JF)=.FALSE.
    ENDIF
  ENDIF
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPCLIPHY',1,ZHOOK_HANDLE)
END SUBROUTINE FPCLIPHY
