! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE GPINCOCH(YDCHEM,KSTART,KPROF,KLEV,KPROMA,KGPLAT,PKCO,PTOP)

!**** *GPINCOCH* - Linear CO chemistry.
 
!     Purpose.  
!     --------
!      Memory transfer: fills PKCO with the content of buffer TCO2DG.

!**   Interface.
!     ----------
!        *CALL* *GPINCOCH(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!         KSTART       : start of work
!         KPROF        : depth of work
!         KLEV         :  Number of Levels
!         KPROMA       : horizontal dimension
!         KGPLAT       : DM-global number of the latitude of point jrof=KSTART

!        OUTPUT:
!         PKCOO        : fields for photochemistery of carbon monoxide
!         PTOP         : fields for carbon monoxide top condition

!        Implicit arguments :
!        --------------------
  
!     Method.
!     -------
!        See documentation

!     Externals.  None. 
!     ----------
!      Called by EC_PHYS.
    
!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        S. Massart (CERFACS/ECMWF) from adiab/gpino3ch.F90

! Modifications.
! --------------
!   Original : 
! End Modifications
!------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCHEM  , ONLY : TCHEM

! -------------------------------------------------------------------------

IMPLICIT NONE

TYPE(TCHEM)       ,INTENT(INOUT) :: YDCHEM
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPLAT(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKCO(KPROMA,KLEV*YDCHEM%NCHEM_LCOCOEF) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOP(KPROMA)

! -------------------------------------------------------------------------

INTEGER(KIND=JPIM) ::  IGL, JJ, JVAR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


! -------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPINCOCH',0,ZHOOK_HANDLE)
ASSOCIATE(NCHEM_LCOCOEF=>YDCHEM%NCHEM_LCOCOEF, &
 & LCHEM_LCOCSTCLIM=>YDCHEM%LCHEM_LCOCSTCLIM)
! -------------------------------------------------------------------------

DO JVAR=1,KLEV*NCHEM_LCOCOEF
  DO JJ=KSTART,KPROF
    IGL=KGPLAT(JJ)
    PKCO(JJ,JVAR)=YDCHEM%TCO2DG(JVAR,IGL)
  ENDDO
ENDDO

IF (.NOT. LCHEM_LCOCSTCLIM) THEN
  DO JJ=KSTART,KPROF
    IGL=KGPLAT(JJ)
    PTOP(JJ) = YDCHEM%TCOTOP(IGL)
  ENDDO  
ENDIF



! -------------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPINCOCH',1,ZHOOK_HANDLE)
END SUBROUTINE GPINCOCH
