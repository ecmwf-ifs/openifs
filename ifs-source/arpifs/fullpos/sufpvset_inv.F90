! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPVSET_INV(YDQTYPE,KFPTRV,KVSETSC,KVSETVD)

!**** *SUFPVSET_INV * - Set up V-set for Full-Pos inverse transforms.

!     Purpose
!     -------
!       To Set up V-set for Full-Pos inverse transforms

!**   Interface.  CALL SUFPVSET_INV(...)
!     ---------- 

!     Explicit arguments : 
!     --------------------
!        KFPTRV  : Number of MPI tasks for fields distribution
!        KVSETSC : V-set for each scalar field
!        KVSETVD : V-set for each pair of field (Vorticity,Divergence)
!                  (they MUST be on the same V-set)

!     Externals.
!     ----------
!        None.

!     Reference.
!     ----------
!        Fullpos technical & users guide.

!     Author.
!     -------
!        Ryad El Khatib *METEO-FRANCE*

!     Modifications.
!     --------------
!        Original : 03-02-11 From TRANSINV_FP
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        Ryad El Khatib 03-Aug-2012 Cleanings
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(IN) :: YDQTYPE
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPTRV
INTEGER(KIND=JPIM),INTENT(OUT)   :: KVSETSC(:) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KVSETVD(:) 
INTEGER(KIND=JPIM) :: J, INC, IJ, IPSCA, IPVEC, JLEV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFPVSET_INV',0,ZHOOK_HANDLE)
SELECT CASE (KFPTRV)
CASE (1)
  KVSETSC(:)=1
  KVSETVD(:)=1
CASE (2:)
  IPSCA=0
  IPVEC=0
  DO J=1,YDQTYPE%NFPOSDYN
    IJ=YDQTYPE%NFPTRDYN(J)
    SELECT CASE (YDQTYPE%ISF(IJ))
    CASE (1:)
      DO JLEV=1,YDQTYPE%ILEV(IJ)
        INC=YDQTYPE%ISPD(JLEV,IJ)
        SELECT CASE (YDQTYPE%ILED(IJ))
        CASE (:-1)
!               V (same V-set as U) :
          KVSETVD(IPVEC+1:IPVEC+INC)=YDQTYPE%ISET(JLEV,IJ)
          IPVEC=IPVEC+INC
        CASE (0,2,3)
!               "Scalar"
          KVSETSC(IPSCA+1:IPSCA+INC)=YDQTYPE%ISET(JLEV,IJ)
          IPSCA=IPSCA+INC
        END SELECT
      ENDDO
    END SELECT
  ENDDO
  IF (IPSCA /= SIZE(KVSETSC)) CALL ABOR1('SUFPVSET_INV : INTERNAL ERROR IPSCA')
  IF (IPVEC /= SIZE(KVSETVD)) CALL ABOR1('SUFPVSET_INV : INTERNAL ERROR IPVEC')
END SELECT

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFPVSET_INV',1,ZHOOK_HANDLE)

END SUBROUTINE SUFPVSET_INV
