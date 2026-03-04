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

SUBROUTINE SPNORMAVE(KSMAX,PMFG,KDIM1,KLEVXG,LDKINER,PNORM,PAVE)

!**** *SPNORMAVE* - Compute averages of norms in spectral space

!     Purpose.
!     --------
!        Compute the average of norms in spectral space for each level
!        and for all levels

!**   Interface.
!     ----------
!        *CALL* *SPNORMAVE(...)

!        Explicit arguments :
!        --------------------
!        KSMAX   : truncation
!        PMFG    : spectrum for each n
!        KDIM1   : first dimension of input and output arrays PMFG and PNORM
!        KLEVXG  : number of levels to compute
!        LDKINER : true if kinetic energy computation
!        PNORM   : Average norm for each level
!        PAVE    : Average norm for all levels

!        Implicit arguments :  none.
!        --------------------

!     Method.
!     -------

!     Externals. None
!     ----------
!      Called by SPNORMBM.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Philippe Courtier  *ECMWF*
!      Original : 92-12-20

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KSMAX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIM1 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFG(KDIM1,0:KSMAX)
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVXG 
LOGICAL           ,INTENT(IN)    :: LDKINER 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNORM(KDIM1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAVE 

!     ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPNORMAVE',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    Computes norm.
!              --------------

!CDIR NOVECTOR
CALL GSTATS(1950,0)
PNORM(1:KLEVXG)=SUM(PMFG(1:KLEVXG,:),DIM=2)
IF (.NOT.LDKINER) THEN
  PNORM(1:KLEVXG)=SQRT(PNORM(1:KLEVXG))
ENDIF

PAVE=SUM(PNORM(1:KLEVXG))
PAVE=PAVE/REAL(KLEVXG,JPRB)
CALL GSTATS(1950,1)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPNORMAVE',1,ZHOOK_HANDLE)
END SUBROUTINE SPNORMAVE
