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

SUBROUTINE GPHLWI(YDDIMV,KPROMA,KSTART,KPROF,PLNPR,PALPH,PWW)

!**** *GPHLWI* - to half-levels interpolation weights

!     Purpose.
!     --------
!           Compute weights for interpolation of winds to
!           half levels (in non-hydrostatic dynamics)

!**   Interface.
!     ----------
!        *CALL* *GPHLWI()

!        Explicit arguments :
!        --------------------
!         INPUT:
!          KPROMA  - dimensioning.
!          KSTART  - start of work.
!          KPROF   - depth of work.
!          PLNPR   - "delta" (log(prehyd) depth) at full levels.
!          PALPH   - "alpha" at full levels.

!         OUTPUT:
!          PWW     - vertical weight.

!        Implicit arguments : none.
!        --------------------

!     Method.
!     -------
!        Interpolation from full-levels using logaritmic pressure profile
!        Then modify values on the top and bottom using boundary condition
!        (at the bottom free-slip condition)
!        Store the weights of vertical interpolation

!     Externals.
!     ----------
!     Reference.
!     ----------
!        ARPEGE/ALADIN documentation

!     Author.
!     -------
!        Radmila Bubnova,  CNRM/GMAP/EXT

!     Modifications.
!     --------------
!        Original : November 1997
!        Modified 02-07-02 by C. Fischer : remove pww5 computation + intents
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        K.Yessad      20-Mar-2006 Cleaning + optimisation
!        K. Yessad (Dec 2008): remove dummy CDLOCK
!     ----------------------------------------------------------------------

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIMV),       INTENT(IN)    :: YDDIMV
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLNPR(KPROMA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPH(KPROMA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWW(KPROMA,YDDIMV%NFLEVG) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPHLWI',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG)

!     ------------------------------------------------------------------

DO JLEV=1, NFLEVG - 1
  DO JROF=KSTART,KPROF
    PWW(JROF,JLEV)=(PLNPR(JROF,JLEV+1)-PALPH(JROF,JLEV+1))&
     & /(PLNPR(JROF,JLEV+1)-PALPH(JROF,JLEV+1)+PALPH(JROF,JLEV))  
  ENDDO
ENDDO

!     ------------------------------------------------------------------

DO JROF=KSTART,KPROF
  !* Bottom half level: 
  PWW(JROF,NFLEVG) = 0.0_JPRB
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPHLWI',1,ZHOOK_HANDLE)
END SUBROUTINE GPHLWI
