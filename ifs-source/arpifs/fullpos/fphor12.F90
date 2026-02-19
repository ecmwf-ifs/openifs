! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPHOR12(KFLD,KST,KEND,PROW,KFLDROW,KFPROMA,KORR,PSEA,&
 & PLSM,PMIN,PMAX,PSCAL,PTT)  

!**** *FPHOR12*  - FULL-POS corrector after horizontal interpolations,
!                  in order to prevent from overshoots due to the
!                  interpolations.

!     PURPOSE.
!     --------
!        Correct each point of a row, according to a rule defined by
!        KORR :
!        - option 1 is usually for ECMWF skin reservoir contents,snow depth :
!                   The output field is filtered in order to remain above
!                   a minimum value over land and set to a fixed value over sea
!        - option 2 is usually for snow cover :
!                   The output field is filtered in order to remain above
!                   zero.Then it is reset to zero over sea, and over land
!                   if the surface temperature at this point is below a
!                   threshold.
!        - option 3 is usually for wetness or constant fields :
!                   The output field is filtered in order to remain between
!                   a minimum and a maximum value ; then it is reset to a
!                   given value over sea
!        - option 4 is usually for positive fields :
!                   The output field is filtered in order to remain above
!                   a minimum value
!        - option 5 is usually for land-sea mask
!                   The output field is modified in order to be 0. or 1.
!        - option 6 is usually for normed fields :
!                   The output field is filtered in order to remain between
!                   a minimum value and a maximum value.
!        - option 7 is usually for snow cover :
!                   The output field is filtered in order to remain above
!                   zero.Then it is reset to zero over sea, but contradictory to 
!                   option 2 no threshold on surface temperature is used.
!
!**   INTERFACE.
!     ----------
!       *CALL* *FPHOR12*

!        EXPLICIT ARGUMENTS
!        --------------------

!        INPUT:
!         KFLD   : field index in output row.
!         KST    : first output point in row.
!         KEND   : last output point in row.
!         KFLDROW: total number of fields in row.
!         KFPROMA: length of the output row.
!         KORR   : kind of correction to perform after interpolations.
!         PSEA   : auxilary (scaled) field value on sea
!         PLSM   : land-sea mask.
!         PMIN   : minimum value of the (scaled) field in output.
!         PMAX   : maximum value of the (scaled) field in output.
!         PSCAL  : scaling factor.
!         PTT    : threshold.

!        OUTPUT:
!         PROW   : interpolated field.

!        IMPLICIT ARGUMENTS
!        --------------------
!          None

!     METHOD.
!     -------

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
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!      R. El Khatib 03-04-23 merge KCOR=10, 8 and 3 => new KCOR=3 + remove KORR=9
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!        T. Aspelien  11-04-11 Added option 7 for snow
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDROW 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLD 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PROW(KFPROMA,KFLDROW) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KORR 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSEA(KFPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KFPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMIN(KFPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMAX(KFPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSCAL(KFPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT(KFPROMA) 
INTEGER(KIND=JPIM) :: JFP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

!*       1. SPECIFIC CORRECTIONS
!           --------------------

IF (LHOOK) CALL DR_HOOK('FPHOR12',0,ZHOOK_HANDLE)
SELECT CASE (KORR)
CASE (1)
  DO JFP=KST,KEND
    PROW(JFP,KFLD) = PLSM(JFP)*MAX(PMIN(JFP),PROW(JFP,KFLD)) &
     & + (1.0_JPRB-PLSM(JFP))*PSEA(JFP)  
  ENDDO
CASE (2)
  DO JFP=KST,KEND
    PROW(JFP,KFLD) = PSCAL(JFP)*MAX(0.0_JPRB,PROW(JFP,KFLD))*PLSM(JFP) &
     & *MAX(0.0_JPRB,SIGN(1.0_JPRB,PTT(JFP)-PSEA(JFP)))  
  ENDDO
CASE (3)
  DO JFP=KST,KEND
    PROW(JFP,KFLD) = MIN(PMAX(JFP),MAX(PMIN(JFP),PROW(JFP,KFLD)))
    PROW(JFP,KFLD) = PSCAL(JFP)*(PLSM(JFP)*PROW(JFP,KFLD) &
     & + (1.0_JPRB-PLSM(JFP))*PSEA(JFP))  
  ENDDO
CASE (4)
  DO JFP=KST,KEND
    PROW(JFP,KFLD) = MAX(PMIN(JFP),PROW(JFP,KFLD))
  ENDDO
CASE (5)
  DO JFP=KST,KEND
    PROW(JFP,KFLD) = MAX(0.0_JPRB,SIGN(1.0_JPRB,PROW(JFP,KFLD)-PTT(JFP)))
  ENDDO
CASE (6)
  DO JFP=KST,KEND
    PROW(JFP,KFLD)=MIN(PMAX(JFP),MAX(PMIN(JFP),PROW(JFP,KFLD)))
  ENDDO
CASE (7)
  DO JFP=KST,KEND
    PROW(JFP,KFLD) = PSCAL(JFP)*MAX(0.0_JPRB,PROW(JFP,KFLD))*PLSM(JFP)
  ENDDO
CASE DEFAULT
  CALL ABOR1(' FPHOR12 : INTERNAL ERROR ON KORR')
END SELECT

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPHOR12',1,ZHOOK_HANDLE)
END SUBROUTINE FPHOR12
