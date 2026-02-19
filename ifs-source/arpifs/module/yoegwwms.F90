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

MODULE YOEGWWMS

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE  

! YOEGWWMS --  parameters/switches for the Warner McIntyre GW parameterization

!-----------------------------------------------------------------

TYPE :: TEGWWMS
LOGICAL            :: LOZPR     !If .TRUE. then enhancement of mom. flux over tropics
LOGICAL            :: LGACALC   !If .TRUE. recalculate vertical coordinate stretch for equal ratio
LOGICAL            :: LGSATL    !If .TRUE. then spectrum is saturated up to m* at launch
LOGICAL            :: LGINDL    !If .TRUE. then saturation and launch spectrum is independent of launch location
INTEGER(KIND=JPIM) :: NLAUNCHL(3)!launch levels of gravity wave spectrum (Pa)
INTEGER(KIND=JPIM) :: NLAUNCHLEV!Number of independent launch heights to be treated
INTEGER(KIND=JPIM) :: NLAUNCH   !single launch level as used for simplified
INTEGER(KIND=JPIM) :: NSLOPE    !slope at small-m end of spectrum
INTEGER(KIND=JPIM) :: NGAUSS    !if LOZPR=TRUE and GGAUSS=2 then gaussian distribution of GFLUXLAUN based on GGAUSSA and GGAUSSB
                                !if LOZPR=TRUE and GGAUSS=1 then FLUXLAUN=GCOEFF*PPRECIP
REAL(KIND=JPRB)    :: GFLUXLAUNL(3)!total launch momentum flux (Pa) in each azimuth for all launch evels
REAL(KIND=JPRB)    :: GFLUXLAUN !total launch momentum flux (Pa) in each azimuth at NLAUNCH
REAL(KIND=JPRB)    :: GMSTAR_L(3)!m* (expressed as length, see Scinocca 2003)
REAL(KIND=JPRB)    :: GCSTAR    !C* (see McLandress and Scinocca 2005)
REAL (KIND=JPRB)   :: GPTWO     !2*p where p is the exponent of omega for the expression of the launch energy density
REAL(KIND=JPRB)    :: GTPHYGWWMS!Time frequency (s) of call of scheme

REAL(KIND=JPRB)    :: GGAUSSA   !gaussian distribution half-width (used if GGAUSS=1)
REAL(KIND=JPRB)    :: GGAUSSB(3)!height of gaussian distribution (ie amplification factor)
REAL(KIND=JPRB)    :: GCOEFF    !if GGAUSS=1 then FLUXLAUN=GCOEFF*PPRECIP
!---------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TEGWWMS

!!TYPE(TEGWWMS), POINTER :: YREGWWMS => NULL()
!---------------------------------------------------------------------

CONTAINS 
  
SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TEGWWMS), INTENT(IN) :: SELF
INTEGER       , INTENT(IN) :: KDEPTH
INTEGER       , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC

IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH)    // 'model%yrml_phy_ec%yregwwms : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LOZPR = ', SELF%LOZPR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LGACALC = ', SELF%LGACALC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LGSATL = ', SELF%LGSATL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LGINDL = ', SELF%LGINDL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLAUNCHL sum ', SUM(SELF%NLAUNCHL)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLAUNCHLEV = ', SELF%NLAUNCHLEV
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLAUNCH = ', SELF%NLAUNCH
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSLOPE = ', SELF%NSLOPE
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NGAUSS = ', SELF%NGAUSS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GFLUXLAUNL sum ',SUM(SELF%GFLUXLAUNL)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GFLUXLAUN = ', SELF%GFLUXLAUN
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GMSTAR_L sum ',SUM(SELF%GMSTAR_L)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GCSTAR = ', SELF%GCSTAR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GPTWO = ', SELF%GPTWO
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GTPHYGWWMS = ', SELF%GTPHYGWWMS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GGAUSSA = ', SELF%GGAUSSA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GGAUSSB sum',SUM(SELF%GGAUSSB)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GCOEFF = ', SELF%GCOEFF

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOEGWWMS
