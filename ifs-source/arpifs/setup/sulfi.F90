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

SUBROUTINE SULFI(KULOUT,KNMAX,KNCMAX,PSTRET,KFACTD)

!**** *SULFI* - Routine to initialize characteristics for LFI files
!               (designed for contraction/dilatation matrixes)

!     Purpose / But
!     -------------
!           Read namelist NAMLFI

!**   Interface
!     ---------
!        *CALL* *SULFI(...)*

!        Explicit arguments / Arguments explicites :
!        -------------------------------------------
! input / entree :
!              KULOUT   :  Logical unit for output
!              KNMAX    :  Lower truncation
!              KNCMAX   :  Upper truncation
!              PSTRET   :  Stretching coefficient
! output / sortie :
!              KFACTD   :  Factor to enhance I/O on contr./dilat. matrixes

!        Implicit arguments / Arguments implicites :
!        -------------------------------------------

!     Reference
!     ---------
!        Note de travail ARPEGE n12 et additifs

!     Author / Auteur
!     ---------------
!      D. Giard  *METEO-FRANCE*
!      Original    : 92-12-08

!     Modifications
!     -------------
!      R. El Khatib : 01-08-07 Pruning options
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNMAX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNCMAX 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTRET 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KFACTD 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZFACTD,ZFACTDREF,ZSIZE,ZSIZEREF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SULFI',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

! * Granularity factor
IF (PSTRET > 1.0_JPRB) THEN
  ZSIZEREF=6232525._JPRB
  ZFACTDREF=10._JPRB
  ZSIZE=0.5_JPRB*REAL(KNMAX+1,JPRB)*REAL(KNMAX+2,JPRB)*REAL(KNCMAX+1,JPRB)&
   & -(1.0_JPRB/6._JPRB)*REAL(KNMAX,JPRB)*REAL(KNMAX+1,JPRB)*REAL(KNMAX+2,JPRB)
  ZFACTD=ZFACTDREF*ZSIZE/ZSIZEREF
  KFACTD=MAX(1,NINT(ZFACTD))
ELSE
  KFACTD=1
ENDIF

! * Upper bound to 120 to be consistent with the limits of the LFI software.
KFACTD=MIN(MAX(KFACTD,1),120)

! * Checkings.
WRITE(UNIT=KULOUT,FMT='(''   KFACTD ='',I3)') KFACTD

! ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SULFI',1,ZHOOK_HANDLE)
END SUBROUTINE SULFI
