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

MODULE PTRGPPC

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!*    Pointers for grid point array for PC schemes (obsolescent feature
!      which has not yet been moved into GMV, GMVS and GFL).
!     These pointers must later to be put in GMV+GMVS (new attributes
!      to create), and also GFL for MGPPCF_GFLP.
!     They are linked to the different PC schemes and they are no
!      longer purely NH pointers, so they have been renamed in Jan 2008.

!  NFGPPC       - number of (2D) fields in PC gridpoint buffer GPPCBUF

!  MGPPC        - whole GPPCBUF-pointer.
!  -- Pointers for LPC_FULL (hydrostatic or NH):
!  MGPPCF_U     - pointer for U.
!  MGPPCF_V     - pointer for V.
!  MGPPCF_T     - pointer for T.
!  MGPPCF_SPD   - pointer for pressure departure variable.
!  MGPPCF_SVD   - pointer for vertical divergence variable.
!  MGPPCF_SP    - pointer for log(prehyds).
!  -- Additional pointers for LPC_CHEAP (hydrostatic or NH):
!     (they are used to save interpolated quantities from the predictor
!     to the corrector step).
!  MGPPCF_CP    - pointer for continuity eqn (3D quantity).
!  MGPPCF_NHX   - pointer for NHX.
!  MGPPCF_UP    - pointer for split EC physics (U eqn).
!  MGPPCF_VP    - pointer for split EC physics (V eqn).
!  MGPPCF_TP    - pointer for split EC physics (T eqn).
!  MGPPCF_GFLP  - pointer for split EC physics (GFL eqn).
!  MGPPCF_BBC   - pointer for BBC quantity (LRDBBC=T).

!  The same names with a suffix "5": pointers for trajectory.

TYPE :: TPTRGPPC
INTEGER(KIND=JPIM) :: NFGPPC

INTEGER(KIND=JPIM) :: MGPPC

INTEGER(KIND=JPIM) :: MGPPCF_U
INTEGER(KIND=JPIM) :: MGPPCF_V
INTEGER(KIND=JPIM) :: MGPPCF_T
INTEGER(KIND=JPIM) :: MGPPCF_SPD
INTEGER(KIND=JPIM) :: MGPPCF_SVD
INTEGER(KIND=JPIM) :: MGPPCF_SP

INTEGER(KIND=JPIM) :: MGPPCF_CP
INTEGER(KIND=JPIM) :: MGPPCF_NHX
INTEGER(KIND=JPIM) :: MGPPCF_UP
INTEGER(KIND=JPIM) :: MGPPCF_VP
INTEGER(KIND=JPIM) :: MGPPCF_TP
INTEGER(KIND=JPIM) :: MGPPCF_GFLP
INTEGER(KIND=JPIM) :: MGPPCF_BBC

INTEGER(KIND=JPIM) :: MGPPC5

INTEGER(KIND=JPIM) :: MGPPCF_U5
INTEGER(KIND=JPIM) :: MGPPCF_V5
INTEGER(KIND=JPIM) :: MGPPCF_T5
INTEGER(KIND=JPIM) :: MGPPCF_SPD5
INTEGER(KIND=JPIM) :: MGPPCF_SVD5
INTEGER(KIND=JPIM) :: MGPPCF_SP5

CONTAINS
  
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
    
END TYPE TPTRGPPC

!!TYPE(TPTRGPPC), POINTER :: YRPTRGPPC => NULL()
  !---------------------------------------------------------------------

CONTAINS 
  
SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TPTRGPPC), INTENT(IN) :: SELF
INTEGER        , INTENT(IN) :: KDEPTH
INTEGER        , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC

IDEPTHLOC = KDEPTH + 2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH) // 'model%yrml_dyn%yrptrgppc : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFGPPC = ', SELF%NFGPPC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPC = ', SELF%MGPPC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_U = ', SELF%MGPPCF_U
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_V = ', SELF%MGPPCF_V
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_T = ', SELF%MGPPCF_T
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_SPD = ', SELF%MGPPCF_SPD
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_SVD = ', SELF%MGPPCF_SVD
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_SP = ', SELF%MGPPCF_SP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_CP = ', SELF%MGPPCF_CP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_NHX = ', SELF%MGPPCF_NHX
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_UP = ', SELF%MGPPCF_UP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_VP = ', SELF%MGPPCF_VP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_TP = ', SELF%MGPPCF_TP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_GFLP = ', SELF%MGPPCF_GFLP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_BBC = ', SELF%MGPPCF_BBC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPC5 = ', SELF%MGPPC5
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_U5 = ', SELF%MGPPCF_U5
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_V5 = ', SELF%MGPPCF_V5
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_T5 = ', SELF%MGPPCF_T5
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_SPD5 = ', SELF%MGPPCF_SPD5
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_SVD5 = ', SELF%MGPPCF_SVD5
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MGPPCF_SP5 = ', SELF%MGPPCF_SP5


END SUBROUTINE PRINT_CONFIGURATION


END MODULE PTRGPPC
