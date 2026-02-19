! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPGEO(YDQTYPE,KST,KEND,KFPROMA,PGNXO,PGNYO,PGMO,KFIELDS,KORDER,PROW)

!**** *FPGEO*  - Full- POS geometry of output dynamical fields

!     PURPOSE.
!     --------
!        Transform the post-processed fields into the output geometry

!**   INTERFACE.
!     ----------
!       *CALL* *FPGEO*

!        EXPLICIT ARGUMENTS
!        --------------------
!       INPUT : 
!     KST    : first output point in row
!     KEND   : last output point in row
!     KFPROMA: length of the output subrow
!     PGNXO  : output X-component of rotation matrix devided by map factor
!     PGNYO  : output Y-component of rotation matrix devided by map factor
!     PGMO   : output map factor **2
!     KFIELDS: number of fields in row
!     KORDER : order of derivative of all fullpos fields
!       INPUT/OUTPUT : 
!     PROW   : modified output row

!        IMPLICIT ARGUMENTS
!        --------------------
!          None

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY   : JPIM     ,JPRB
USE YOMHOOK  , ONLY   : LHOOK,   DR_HOOK, JPHOOK

USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(IN) :: YDQTYPE
INTEGER(KIND=JPIM), INTENT(IN) :: KST 
INTEGER(KIND=JPIM), INTENT(IN) :: KEND 
INTEGER(KIND=JPIM), INTENT(IN) :: KFPROMA 
REAL(KIND=JPRB)   , INTENT(IN) :: PGNXO(KFPROMA) 
REAL(KIND=JPRB)   , INTENT(IN) :: PGNYO(KFPROMA) 
REAL(KIND=JPRB)   , INTENT(IN) :: PGMO(KFPROMA) 
INTEGER(KIND=JPIM), INTENT(IN) :: KFIELDS 
INTEGER(KIND=JPIM), INTENT(IN) :: KORDER(:)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PROW(KFPROMA,KFIELDS) 

! provisional copy of PROW. Needed because U/V in same row !! but could be optimized.
REAL(KIND=JPRB) :: ZROW(KFPROMA,KFIELDS)

INTEGER(KIND=JPIM) :: IJ, IPTR, INC, IU, IV, JF, JL, JF2, IJ2, IPTR2, JI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPGEO',0,ZHOOK_HANDLE)

ZROW(KST:KEND,:)=PROW(KST:KEND,:)

IPTR=1
DO JF=1,YDQTYPE%NFPOSDYN
  IJ=YDQTYPE%NFPTRDYN(JF)
  INC   =YDQTYPE%ILEV(IJ)
  SELECT CASE (KORDER(IJ))
  CASE (2)
    ! 2nd order derivative => multiply with square map factor
    DO JL=IPTR,IPTR+INC-1
      DO JI=KST,KEND
        PROW(JI,JL)=ZROW(JI,JL)*PGMO(JI)
      ENDDO
    ENDDO
  CASE (1)
    ! Vector/U => apply compass & map factor for both momenta
    ! Find second momentum :
    IPTR2=1   
    DO JF2=1,YDQTYPE%NFPOSDYN
      IJ2=YDQTYPE%NFPTRDYN(JF2)
      IF (YDQTYPE%ILED(IJ2) == -IJ) THEN
        EXIT
      ELSE
        IPTR2=IPTR2+YDQTYPE%ILEV(IJ2)
      ENDIF
    ENDDO
    IU=IPTR
    IV=IPTR2
    DO JL=0,INC-1
      DO JI=KST,KEND
        PROW(JI,IU+JL)=ZROW(JI,IU+JL)*PGNYO(JI)+ZROW(JI,IV+JL)*PGNXO(JI)
        PROW(JI,IV+JL)=ZROW(JI,IV+JL)*PGNYO(JI)-ZROW(JI,IU+JL)*PGNXO(JI)
      ENDDO
    ENDDO
  END SELECT
  IPTR=IPTR+INC
ENDDO
IF (IPTR-1 /= KFIELDS) CALL ABOR1('FPGEO : INTERNAL ERROR 1')

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPGEO',1,ZHOOK_HANDLE)
END SUBROUTINE FPGEO
