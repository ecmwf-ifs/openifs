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

MODULE YOELWRAD

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

! ------ CLOUD CHARACTERISTICS FOR SIMPLIFIED SCHEME

! NLEVLWC  : level up to which cloud effects on LW fluxes are considered
! NLWFR    : frequency of calling longwave radiation for TL/AD
! NLOOPLW  : required number of loops inside of nproma for TL/AD lw radiation 
! NLOOPLWO : optimized number of loops inside of nproma
! NPROMALW : optimized length of loops inside of nproma
! NLASTLW  : length of the last loop inside of nproma
! LOPTLWPR : .true. if shorter loops inside of nproma are required 
! LWLCLHR  : .true. if cloud effects on LW are only computed to the level
!                      defined using NLEVLWC and the current cloud top height 
!            .false. if computed for all levels          

TYPE :: TELWRAD
INTEGER(KIND=JPIM) :: NLEVLWC
INTEGER(KIND=JPIM) :: NLWFR
INTEGER(KIND=JPIM) :: NLOOPLW
INTEGER(KIND=JPIM) :: NLOOPLWO
INTEGER(KIND=JPIM) :: NPROMALW
INTEGER(KIND=JPIM) :: NLASTLW
 
LOGICAL :: LOPTLWPR
LOGICAL :: LWLCLHR
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TELWRAD
!============================================================================

!!TYPE(TELWRAD), POINTER :: YRELWRAD => NULL()

!     --------------------------------------------------------------------

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TELWRAD), INTENT(IN) :: SELF
INTEGER       , INTENT(IN) :: KDEPTH
INTEGER       , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC

IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_rad%yrelwrad : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLEVLWC = ', SELF%NLEVLWC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLWFR = ', SELF%NLWFR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLOOPLW = ', SELF%NLOOPLW
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLOOPLWO = ', SELF%NLOOPLWO
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NPROMALW = ', SELF%NPROMALW
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLASTLW = ', SELF%NLASTLW
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LOPTLWPR = ', SELF%LOPTLWPR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LWLCLHR = ', SELF%LWLCLHR

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOELWRAD
