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

MODULE PTRSLB1

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     SEMI LAGRANGIAN BUFFER 1

!     Abbreviations:
!     "O" stands for "origin point", "M" stands for medium point.
!     "NH" stands for "non hydrostatic", "SI" stands for "semi-implicit". 

!     * NFLDSLB1  : Number of fields in semi-lagrangian buffer 1
!     * RPARSL1   : parities for extrapolar extensions for interpol. buffers.
!     * MSLBUF1   : start pointer for interpolation buffer 1.
!     * MSLB1U9   : pointer for U-wind eqn term to be interpolated at O.
!     * MSLB1V9   : pointer for V-wind eqn term to be interpolated at O.
!     * MSLB1Z9   : pointer for Z-wind eqn term to be interpolated at O.
!     * MSLB1T9   : pointer for T eqn term to be interpolated at O.
!     * MSLB1GFL9 : pointer for unified_treatment grid-point fields
!     * MSLB1GFLSP9 : pointer for spline rep. of gfl needed for spline int. at O.
!     * MSLB1PD9  : pointer for P-hat eqn term to be interpolated at O (NH).
!     * MSLB1VD9  : pointer for d-hat eqn term to be interpolated at O (NH).
!     * MSLB1NHX9 : pointer for "X" term eqn term to be interpolated at O (NH).
!     * MSLB1UR0  : pointer for U-wind to be interpolated for traj research.
!     * MSLB1VR0  : pointer for V-wind to be interpolated for traj research.
!     * MSLB1ZR0  : pointer for Z-wind to be interpolated for traj research.
!     * MSLB1WR0  : pointer for "etadot" to be interpolated for traj research.
!     * MSLB1UR00 : pointer for U-wind at t to be interpolated for traj research.
!     * MSLB1VR00 : pointer for V-wind at t to be interpolated for traj research.
!     * MSLB1ZR00 : pointer for Z-wind at t to be interpolated for traj research.
!     * MSLB1WR00 : pointer for "etadot" at t to be interpolated for traj research.
!     * MSLB1WRA  : pointer for "etadot" at the arrival point
!     * MSLB1UR9  : pointer for t-dt U-wind to be interpolated at O for 
!                   Coriolis term refined recalculation (L2TLFF=T).
!     * MSLB1VR9  : pointer for t-dt V-wind to be interpolated at O for
!                   Coriolis term refined recalculation (L2TLFF=T).
!     * MSLB1ZR9  : pointer for t-dt Z-wind to be interpolated at O for
!                   Coriolis term refined recalculation (L2TLFF=T).
!     * MSLB1DBBC9: pointer for intermediate quantity ((p/mRT)*Jacobian part
!                   of Laplacian in Lth level) to be interpolated at O
!                   for diagnostic BBC in NH.
!     * MSLB1DPHI9: pointer for (p/m Rd T) to be interpolated at O for
!                   diagnostic BBC in NH.
!     * MSLB1GWS9 : pointer for (gw_s) for high-order horizontal interp.
!                   at O(L) for diagnostic BBC in NH.
!     * MSLB1U0   : pointer for U-wind eq term to be interpolated linearly at O.
!     * MSLB1V0   : pointer for V-wind eq term to be interpolated linearly at O.
!     * MSLB1Z0   : pointer for Z-wind eq term to be interpolated linearly at O.
!     * MSLB1T0   : pointer for T eqn term to be interpolated linearly at O.
!     * MSLB1PD0  : pointer for P-hat eqn term to be interpolated linearly at O.
!     * MSLB1VD0  : pointer for d-hat eqn term to be interpolated linearly at O.
!     * MSLB1C9   : pointer for 3D continuity eqn term to be interpolated at O.
!     * MSLB1SP9  : pointer for 2D continuity eqn term to be interpolated at O.
!     * MSLB1SP0  : pointer for 2D continuity eqn term to be interpolated
!                   linearly at O.
!     * MSLB1C0   : pointer for 3D continuity eqn term to be interpolated
!                   linearly at O.
!     * MSLB1UP9  : pointer for U-wind eqn term to be interpolated at O:
!                   tendency in split physics (if LSLPHY=.T. only). 
!     * MSLB1VP9  : pointer for V-wind eqn term to be interpolated at O:
!                   tendency in split physics (if LSLPHY=.T. only). 
!     * MSLB1ZP9  : pointer for Z-wind eqn term to be interpolated at O:
!                   tendency in split physics (if LSLPHY=.T. only). 
!     * MSLB1TP9  : pointer for T eqn term to be interpolated at O:
!                   tendency in split physics (if LSLPHY=.T. only). 
!     * MSLB1GFLP9: pointer for unified_treatment grid-point fields
!                   tendency in split physics (if LSLPHY=.T. only). 
!     * MSLB1[X]9_SI: cf. MSLB1[X]9 ([X]=U,V,T,PD,VD,C) for cases where linear
!                   terms must be separately interpolated (contain separate
!                   linear terms).
!     * MSLB1[X]F9: pointers (applicable only when NSPLTHOI /= 0) for
!                   second part of the quantity to be interpolated at O
!                   by high order interpolation. ([X]=U,V,T,GFL,VD,GFLSP)
!     * MSLB1[X]9_NL: cf. MSLB1[X]9_NL ([X]=U,V,T,PD,VD,C) for cases where 
!                   nonlinear tendency must be separated if LSETTLS and LPC_CHEAP

TYPE :: TPTRSLB1
INTEGER(KIND=JPIM) :: NFLDSLB1

REAL(KIND=JPRB),ALLOCATABLE:: RPARSL1(:)

INTEGER(KIND=JPIM) :: MSLBUF1
INTEGER(KIND=JPIM) :: MSLB1U9
INTEGER(KIND=JPIM) :: MSLB1V9
INTEGER(KIND=JPIM) :: MSLB1Z9
INTEGER(KIND=JPIM) :: MSLB1T9
INTEGER(KIND=JPIM) :: MSLB1GFL9
INTEGER(KIND=JPIM) :: MSLB1GFLSP9
INTEGER(KIND=JPIM) :: MSLB1PD9
INTEGER(KIND=JPIM) :: MSLB1VD9
INTEGER(KIND=JPIM) :: MSLB1NHX9
INTEGER(KIND=JPIM) :: MSLB1UR0
INTEGER(KIND=JPIM) :: MSLB1VR0
INTEGER(KIND=JPIM) :: MSLB1ZR0
INTEGER(KIND=JPIM) :: MSLB1WR0
INTEGER(KIND=JPIM) :: MSLB1UR00
INTEGER(KIND=JPIM) :: MSLB1VR00
INTEGER(KIND=JPIM) :: MSLB1ZR00
INTEGER(KIND=JPIM) :: MSLB1WR00
INTEGER(KIND=JPIM) :: MSLB1WRA
INTEGER(KIND=JPIM) :: MSLB1UR9
INTEGER(KIND=JPIM) :: MSLB1VR9
INTEGER(KIND=JPIM) :: MSLB1ZR9
INTEGER(KIND=JPIM) :: MSLB1DBBC9
INTEGER(KIND=JPIM) :: MSLB1DPHI9
INTEGER(KIND=JPIM) :: MSLB1GWS9
INTEGER(KIND=JPIM) :: MSLB1U0
INTEGER(KIND=JPIM) :: MSLB1V0
INTEGER(KIND=JPIM) :: MSLB1Z0
INTEGER(KIND=JPIM) :: MSLB1T0
INTEGER(KIND=JPIM) :: MSLB1PD0
INTEGER(KIND=JPIM) :: MSLB1VD0
INTEGER(KIND=JPIM) :: MSLB1C9
INTEGER(KIND=JPIM) :: MSLB1SP9
INTEGER(KIND=JPIM) :: MSLB1SP0
INTEGER(KIND=JPIM) :: MSLB1C0
INTEGER(KIND=JPIM) :: MSLB1UP9
INTEGER(KIND=JPIM) :: MSLB1VP9
INTEGER(KIND=JPIM) :: MSLB1ZP9
INTEGER(KIND=JPIM) :: MSLB1TP9
INTEGER(KIND=JPIM) :: MSLB1GFLP9
INTEGER(KIND=JPIM) :: MSLB1U9_SI
INTEGER(KIND=JPIM) :: MSLB1V9_SI
INTEGER(KIND=JPIM) :: MSLB1T9_SI
INTEGER(KIND=JPIM) :: MSLB1PD9_SI
INTEGER(KIND=JPIM) :: MSLB1VD9_SI
INTEGER(KIND=JPIM) :: MSLB1C9_SI
INTEGER(KIND=JPIM) :: MSLB1UF9
INTEGER(KIND=JPIM) :: MSLB1VF9
INTEGER(KIND=JPIM) :: MSLB1TF9
INTEGER(KIND=JPIM) :: MSLB1GFLF9
INTEGER(KIND=JPIM) :: MSLB1VDF9
INTEGER(KIND=JPIM) :: MSLB1GFLSPF9
INTEGER(KIND=JPIM) :: MSLB1U9_NL
INTEGER(KIND=JPIM) :: MSLB1V9_NL
INTEGER(KIND=JPIM) :: MSLB1T9_NL
INTEGER(KIND=JPIM) :: MSLB1PD9_NL
INTEGER(KIND=JPIM) :: MSLB1VD9_NL
INTEGER(KIND=JPIM) :: MSLB1C9_NL
!--------------------------------------------------------------
CONTAINS
  
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 

END TYPE TPTRSLB1
!==============================================================

!!TYPE(TPTRSLB1), POINTER :: YRPTRSLB1 => NULL()

CONTAINS 
  
SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TPTRSLB1), INTENT(IN) :: SELF
  INTEGER        , INTENT(IN) :: KDEPTH
  INTEGER        , INTENT(IN) :: KOUTNO

  INTEGER :: IDEPTHLOC

  IDEPTHLOC = KDEPTH + 2
  
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_gconf%ptrslb1 : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFLDSLB1 = ', SELF%NFLDSLB1
  IF (ALLOCATED(SELF%RPARSL1)) THEN
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPARSL1 allocated of shape ', SHAPE(SELF%RPARSL1) ,' and sum ',SUM(SELF%RPARSL1)
  ELSE
    WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPARSL1 not allocated'
  ENDIF
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLBUF1 = ', SELF%MSLBUF1
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1U9 = ', SELF%MSLB1U9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1V9 = ', SELF%MSLB1V9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1T9 = ', SELF%MSLB1T9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1GFL9 = ', SELF%MSLB1GFL9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1GFLSP9 = ', SELF%MSLB1GFLSP9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1PD9 = ', SELF%MSLB1PD9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1VD9 = ', SELF%MSLB1VD9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1NHX9 = ', SELF%MSLB1NHX9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1UR0 = ', SELF%MSLB1UR0
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1VR0 = ', SELF%MSLB1VR0
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1WR0 = ', SELF%MSLB1WR0
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1WRA = ', SELF%MSLB1WRA
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1UR9 = ', SELF%MSLB1UR9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1VR9 = ', SELF%MSLB1VR9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1DBBC9 = ', SELF%MSLB1DBBC9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1DPHI9 = ', SELF%MSLB1DPHI9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1GWS9 = ', SELF%MSLB1GWS9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1U0 = ', SELF%MSLB1U0
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1V0 = ', SELF%MSLB1V0
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1T0 = ', SELF%MSLB1T0
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1PD0 = ', SELF%MSLB1PD0
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1VD0 = ', SELF%MSLB1VD0
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1C9 = ', SELF%MSLB1C9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1SP9 = ', SELF%MSLB1SP9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1SP0 = ', SELF%MSLB1SP0
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1C0 = ', SELF%MSLB1C0
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1UP9 = ', SELF%MSLB1UP9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1VP9 = ', SELF%MSLB1VP9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1ZP9 = ', SELF%MSLB1ZP9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1TP9 = ', SELF%MSLB1TP9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1GFLP9 = ', SELF%MSLB1GFLP9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1U9_SI = ', SELF%MSLB1U9_SI
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1V9_SI = ', SELF%MSLB1V9_SI
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1T9_SI = ', SELF%MSLB1T9_SI
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1PD9_SI = ', SELF%MSLB1PD9_SI
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1VD9_SI = ', SELF%MSLB1VD9_SI
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1C9_SI = ', SELF%MSLB1C9_SI
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1UF9 = ', SELF%MSLB1UF9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1VF9 = ', SELF%MSLB1VF9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1TF9 = ', SELF%MSLB1TF9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1GFLF9 = ', SELF%MSLB1GFLF9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1VDF9 = ', SELF%MSLB1VDF9
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MSLB1GFLSPF9 = ', SELF%MSLB1GFLSPF9
END SUBROUTINE PRINT_CONFIGURATION

END MODULE PTRSLB1
