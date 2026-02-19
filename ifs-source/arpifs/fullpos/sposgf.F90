! (C) Copyright 1989- Meteo-France.

SUBROUTINE SPOSGF(YDQTYPE,KRESOL,PK,PSPBFP,LDACTIVE)

!**** *SPOSGF*  - SPECTRAL POST-PROCESSING - GAUSSIAN FILTERS

!     PURPOSE.
!     --------
!        To modify the post-processed fields in the spectral space
!        by applying a gaussian filter on the transformed sphere

!**   INTERFACE.
!     ----------
!       *CALL* *SPOSGF*

!        EXPLICIT ARGUMENTS
!        --------------------
!          KRESOL : Resolution indicator
!          PK : tuning coefficient for the gaussian filter in the model spectral space,
!               for all fullpos fields
!          PSPBFP : spectral data arrray
!          LDACTIVE : .TRUE. if at least one field may have been filtered

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
!      ORIGINAL : 08-Aug-2012 from SPOS

!     MODIFICATIONS.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : MYSETV
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(IN) :: YDQTYPE
INTEGER(KIND=JPIM), INTENT(IN) :: KRESOL
REAL(KIND=JPRB),    INTENT(IN) :: PK(:)
REAL(KIND=JPRB), INTENT(INOUT) :: PSPBFP(:,:)
LOGICAL, INTENT(OUT) :: LDACTIVE

INTEGER(KIND=JPIM) :: JDOM, JSE, J, IJ, JL, ISPB, IFPSPB, ISPEC2, ISMAX, IMSMAX
INTEGER(KIND=JPIM) :: INVALUE(SIZE(PSPBFP,DIM=2))
INTEGER(KIND=JPIM) :: IMVALUE(SIZE(PSPBFP,DIM=2))
REAL(KIND=JPRB)    :: ZZ(SIZE(PSPBFP,DIM=2))
REAL(KIND=JPRB)    :: ZK, Z1, Z2, Z1SMX2, Z1SNX2
LOGICAL :: LLETRANS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "trans_inq.h"
#include "etrans_inq.h"

#include "abor1.intfb.h"
#include "updtrans.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPOSGF',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    INQUIRIES
!              ---------

ISPEC2=SIZE(PSPBFP,DIM=2)
IFPSPB=SIZE(PSPBFP,DIM=1)

DO J=1,YDQTYPE%NFPOSDYN
  IJ=YDQTYPE%NFPTRDYN(J)
  LDACTIVE=(YDQTYPE%ISF(IJ)==2)
  IF (LDACTIVE) EXIT
ENDDO

!*       2.    FILTER
!              ------

IF (LDACTIVE) THEN
  CALL UPDTRANS(KRESOL,LLETRANS)
  IF (LLETRANS) THEN
    CALL ETRANS_INQ(KRESOL=KRESOL,KSMAX=ISMAX,KMSMAX=IMSMAX,KNVALUE=INVALUE(:),KMVALUE=IMVALUE(:))
    IF (ISMAX > 0) THEN
      Z1SNX2=1.0_JPRB/REAL(ISMAX**2,JPRB)
    ELSE
      Z1SNX2=0.0_JPRB
    ENDIF
    IF (IMSMAX > 0) THEN
      Z1SMX2=1.0_JPRB/REAL(IMSMAX,JPRB)**2
    ELSE
      Z1SMX2=0.0_JPRB
    ENDIF
    DO JSE=1,ISPEC2
      Z1=REAL(INVALUE(JSE),JPRB)**2
      Z2=REAL(IMVALUE(JSE),JPRB)**2
      ZZ(JSE)=(Z1*Z1SNX2)+(Z2*Z1SMX2)
    ENDDO
  ELSE
    CALL TRANS_INQ(KRESOL=KRESOL,KSMAX=ISMAX,KNVALUE=INVALUE(:))
    Z1SNX2=1.0_JPRB/REAL(ISMAX**2,JPRB)
    DO JSE=1,ISPEC2
      Z2=REAL(INVALUE(JSE),JPRB)**2
      ZZ(JSE)=Z2*Z1SNX2
    ENDDO
  ENDIF

  ISPB=1
  DO J=1,YDQTYPE%NFPOSDYN
    IJ=YDQTYPE%NFPTRDYN(J)
    IF (YDQTYPE%ISF(IJ) > 0) THEN
      IF (YDQTYPE%ISF(IJ)==2) ZK=-0.5_JPRB*PK(IJ)
      DO JL=1,YDQTYPE%ILEV(IJ)
        IF (YDQTYPE%ISET(JL,IJ)==MYSETV) THEN
          DO JDOM=1,YDQTYPE%ISPD(JL,IJ) 
            IF (YDQTYPE%ISF(IJ)==2) THEN
              DO JSE=1,ISPEC2
                PSPBFP(ISPB,JSE)=PSPBFP(ISPB,JSE)*EXP(ZK*ZZ(JSE))
              ENDDO
            ENDIF
            ISPB=ISPB+1
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  IF (ISPB-1 /= IFPSPB) CALL ABOR1('SPOSGF:INTERNAL ERROR IN ISPB')

ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPOSGF',1,ZHOOK_HANDLE)
END SUBROUTINE SPOSGF
