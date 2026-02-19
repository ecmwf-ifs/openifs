! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPUV2KP(YDQTYPE,YDTFP,KRESOL,PSPBFP)

!**** *FPUV2KP*  - CONVERT (U,V) INTO (PSI,KHI) - FULL POS

!     PURPOSE.
!     --------
!        To transforms wind into stream function and velocity potential 
!        (via vorticity and divergence).

!**   INTERFACE.
!     ----------
!       *CALL* *FPUV2KP*

!        EXPLICIT ARGUMENTS
!        --------------------
!          KRESOL : Resolution indicator
!          PSPBFP : spectral data arrray

!        IMPLICIT ARGUMENTS
!        --------------------
!          See #include below.

!     METHOD.
!     -------
!        See documentation about FULL-POS.

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 07-Aug-2012 from SPOS

!     MODIFICATIONS.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : MYSETV
USE YOMAFN   , ONLY : ALL_FULLPOS_TYPES
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(IN) :: YDQTYPE
TYPE(ALL_FULLPOS_TYPES), INTENT(IN) :: YDTFP
INTEGER(KIND=JPIM), INTENT(IN) :: KRESOL
REAL(KIND=JPRB), INTENT(INOUT) :: PSPBFP(:,:)

LOGICAL :: LLPSIKHI, LLETRANS
INTEGER(KIND=JPIM) :: JDOM, JSE, J, IJ, JL, ISPB, ISPEC2, ISMAX
REAL(KIND=JPRB),    ALLOCATABLE :: ZLAPIN(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: INVALUE(:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
#include "trans_inq.h"

#include "abor1.intfb.h"

#include "updtrans.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FPUV2KP',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    INQUIRIES
!              ----------

CALL UPDTRANS(KRESOL,LLETRANS)
ISPEC2=SIZE(PSPBFP,DIM=2)

DO J=1,YDQTYPE%NFPOSDYN
  IJ=YDQTYPE%NFPTRDYN(J)
  LLPSIKHI=(IJ == YDTFP%KHI%ICOD .OR. IJ == YDTFP%PSI%ICOD)
  IF (LLPSIKHI) EXIT
ENDDO

!*       2.    COMPUTATION OF PSI AND KHI
!              --------------------------

IF (LLPSIKHI) THEN

  IF (LLETRANS) THEN
    CALL ABOR1('FPUV2KP : COMPUTATION OF PSI/KHI NOT YET SUPPORTED IN LAM GEOMETRY')
  ELSE
    CALL TRANS_INQ(KRESOL=KRESOL,KSMAX=ISMAX)
    ALLOCATE(ZLAPIN(-1:ISMAX+2),INVALUE(ISPEC2))
    CALL TRANS_INQ(KRESOL=KRESOL,PLAPIN=ZLAPIN(-1:ISMAX+2),KNVALUE=INVALUE(:))
  ENDIF

  ISPB=1
  DO J=1,YDQTYPE%NFPOSDYN
    IJ=YDQTYPE%NFPTRDYN(J)
    IF (YDQTYPE%ISF(IJ) > 0) THEN
      DO JL=1,YDQTYPE%ILEV(IJ)
        IF (YDQTYPE%ISET(JL,IJ)==MYSETV) THEN
!CDIR NOVECTOR
          DO JDOM=1,YDQTYPE%ISPD(JL,IJ)
            IF (IJ == YDTFP%KHI%ICOD .OR. IJ == YDTFP%PSI%ICOD) THEN
              DO JSE=1,ISPEC2
                PSPBFP(ISPB,JSE) = PSPBFP(ISPB,JSE)*ZLAPIN(INVALUE(JSE))
              ENDDO
            ENDIF
            ISPB=ISPB+1
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  IF (ISPB-1 /= SIZE(PSPBFP,DIM=1)) CALL ABOR1('FPUV2KP:INTERNAL ERROR IN ISPB')

  DEALLOCATE(ZLAPIN,INVALUE)

ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FPUV2KP',1,ZHOOK_HANDLE)
END SUBROUTINE FPUV2KP
