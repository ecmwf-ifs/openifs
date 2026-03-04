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

MODULE YOPHLC

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOPHLC* - SWITCHES FOR SIMPLE PHYSICS
!     -----------------------------------------------------------------

!        * E.C.M.W.F. PHYSICS PACKAGE *

TYPE :: TPHLC
REAL(KIND=JPRB) :: ALPHA
REAL(KIND=JPRB) :: AH0
REAL(KIND=JPRB) :: USTARL
REAL(KIND=JPRB) :: USTARS
REAL(KIND=JPRB) :: ALANDZ0
REAL(KIND=JPRB) :: ASEAZ0
LOGICAL :: LSPHLC
LOGICAL :: LVDFLC
LOGICAL :: LSDRLC
LOGICAL :: LZMCON
LOGICAL :: LKEXP
LOGICAL :: LVDFDS
LOGICAL :: LSDRDS
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TPHLC
!============================================================================

!! TYPE(TPHLC), POINTER :: YRPHLC => NULL()

!     ------------------------------------------------------------------
!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! LSPHLC : LOGICAL : SWITCH TO TURN ON THE SIMPLE PHYSICS SCHEME
! LVDFLC : LOGICAL : TO TURN ON THE VERTICAL DIFFUSION ON U,V
! LVDFDS : LOGICAL : TO TURN ON THE VERT DIFF ALSO ON THE DRY STATIC ENERGY
! LSDRLC : LOGICAL : TO TURN ON THE SURFACE DRAG ON U,V
! LSDRDS : LOGICAL : TO TURN ON THE SURFACE DRAG ALSO ON THE DRY STATIC ENERGY
! LZMCON : LOGICAL : IF TRUE THE VERT DIFF COEFF K IS CONSTANT
! LKEXP  : LOGICAL : EXPONENTIAL APPROXIMATION OF THE CONSTANT K
! ALPHA  : REAL    : COEFFICIENT FOR THE IMPLICIT TIME STEP SCHEME
! AH0    : REAL    : REFERENCE HEIGHT FOR CONSTANT VERT DIFF (M)
! USTARL : REAL    : CONSTANT USTAR OVER LAND (M SEC-1)
! USTARS : REAL    : CONSTANT USTAR OVER SEA (M SEC-1)
! ALANDZ0: REAL    : CONSTANT Z0 OVER LAND (M)
! ASEAZ0 : REAL    : CONSTANT Z0 OVER SEA (M)
! NDIFLC : INTEGER : CONTROL DIFFERENT VERT DIFF SCHEMES
!                   - 1: CALL VDIFLCZ (DIRECT+TL+ADJ):
!                        X(T+1)-X(T-1)=K*X(T+1)
!                        K=CONSTANT
!                   - 2: CALL VDIFLCZ1 (DIRECT):
!                        X(T*)-X(T)=K(X(T))*(ALPHA*X(T*)+(1-ALPHA)*X(T))
!                        DX=X(T*)-X(T)
!                        X(T+1)=X(T-1)+2*DX
!     -----------------------------------------------------------------

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TPHLC), INTENT(IN) :: SELF
INTEGER     , INTENT(IN) :: KDEPTH
INTEGER     , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC

IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_slin%yrphlc : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'ALPHA = ', SELF%ALPHA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'AH0 = ', SELF%AH0
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'USTARL = ', SELF%USTARL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'USTARS = ', SELF%USTARS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'ALANDZ0 = ', SELF%ALANDZ0
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'ASEAZ0 = ', SELF%ASEAZ0
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSPHLC = ', SELF%LSPHLC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LVDFLC = ', SELF%LVDFLC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSDRLC = ', SELF%LSDRLC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LZMCON = ', SELF%LZMCON
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LKEXP = ', SELF%LKEXP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LVDFDS = ', SELF%LVDFDS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSDRDS = ', SELF%LSDRDS

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOPHLC
