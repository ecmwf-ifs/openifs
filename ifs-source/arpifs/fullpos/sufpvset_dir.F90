! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPVSET_DIR(YDQTYPE,KFPTRV,KVSETSC,KVSETVD)

!**** *SUFPVSET_DIR * - Setup V-set for Full-Pos direct transforms.

!     Purpose.
!     --------
!        To setup V-set for Full-Pos direct transforms.

!**   Interface.  CALL SUFPVSET_DIR(...)
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
!        Original : 03-02-11 From TRANSDIR_FP
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

INTEGER(KIND=JPIM) :: J, JF, JLEV, IJ, IPSCA, IPVEC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFPVSET_DIR',0,ZHOOK_HANDLE)
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
        SELECT CASE (YDQTYPE%ILED(IJ))
        CASE (:-1)
!               V (U is on the same V-set)
          DO JF=1,YDQTYPE%ISKP(IJ)
            IPVEC=IPVEC+1
            KVSETVD(IPVEC)=YDQTYPE%ISET(JLEV,IJ)
          ENDDO
        CASE (0)
!               Scalar
          DO JF=1,YDQTYPE%ISKP(IJ)
            IPSCA=IPSCA+1
            KVSETSC(IPSCA)=YDQTYPE%ISET(JLEV,IJ)
          ENDDO
        CASE (2,3)
!               Pair (U,V) for Vor or Div => treat 1 over 2
          DO JF=1,YDQTYPE%ISKP(IJ),2
            IPVEC=IPVEC+1
            KVSETVD(IPVEC)=YDQTYPE%ISET(JLEV,IJ)
          ENDDO
        END SELECT
      ENDDO
    END SELECT
  ENDDO 
  IF (IPSCA /= SIZE(KVSETSC)) CALL ABOR1('SUFPVSET_DIR : INTERNAL ERROR IPSCA')
  IF (IPVEC /= SIZE(KVSETVD)) CALL ABOR1('SUFPVSET_DIR : INTERNAL ERROR IPVEC')
CASE DEFAULT
  CALL ABOR1('SUFPVSET_DIR : INTERNAL ERROR KFPTRV')
END SELECT

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFPVSET_DIR',1,ZHOOK_HANDLE)
END SUBROUTINE SUFPVSET_DIR
