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

MODULE YOE_MCICA

USE PARKIND1  ,ONLY : JPRD, JPIM, JPRB

IMPLICIT NONE

SAVE
!------------------------------------------------------------------------------

TYPE :: TEMCICA

REAL(KIND=JPRB) :: XCW(1000,140)
REAL(KIND=JPRD) :: XCW_D(1000,140)
INTEGER(KIND=JPIM) :: NMCI1, NMCI2
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TEMCICA
!============================================================================

!! TYPE(TEMCICA), POINTER :: YREMCICA => NULL()

!------------------------------------------------------------------------------
! Array (with associated dimension)s used in the Raisanen et al. (2004) cloud generator

! NMcI1 : INTEGER : 1st dimension of the XCW array 
!                   (= maximum number of increments in units of cumulative propability)
! NMcI2 : INTEGER : 2nd dimension of the XCW array (= number of g-points in RRTM_LW)
! XCW   : REAL    : resulting array, precomputed in su_mcica.F90
! XCW_D : REAL    : double prec array used for reading files
!------------------------------------------------------------------------------

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TEMCICA), INTENT(IN) :: SELF
INTEGER       , INTENT(IN) :: KDEPTH
INTEGER       , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC

IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_rad%yremcica : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'XCW SUM = ', SUM(SELF%XCW)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'XCW_D SUM = ', SUM(SELF%XCW_D)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NMCI1 = ', SELF%NMCI1
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NMCI2 = ', SELF%NMCI2

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOE_MCICA
