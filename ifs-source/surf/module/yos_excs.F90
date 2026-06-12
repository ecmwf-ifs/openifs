MODULE YOS_EXCS

USE PARKIND1  ,ONLY : JPIM     ,JPRB

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOS_EXCS* CONTAINS STABILITY FUNCTION TABLES FOR *V...*
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPRITBL=101

REAL(KIND=JPRB) :: RITBL(JPRITBL) ! Tabulated eta-values (z/l) as a function of
                                  ! Richardson number for stable cases.
REAL(KIND=JPRB) :: DRITBL   ! Increment of the Richardson number
REAL(KIND=JPRB) :: RIMAX    ! Maximum Richardson number tabulated
REAL(KIND=JPRB) :: RCHBA    ! Constant a in Holtslag and deBruin functions for stable situations
REAL(KIND=JPRB) :: RCHBB    ! Constant b in HB functions
REAL(KIND=JPRB) :: RCHBC    ! Constant c in HB functions
REAL(KIND=JPRB) :: RCHBD    ! Constant d in HB functions
REAL(KIND=JPRB) :: RCHB23A  ! 2./3.*a in HB functions 
REAL(KIND=JPRB) :: RCHBBCD  ! b*c/d in HB functions
REAL(KIND=JPRB) :: RCHBCD   ! c/d in HB functions
REAL(KIND=JPRB) :: RCHETA   ! Constant in the Hogstrom Ellison Turner functions for stably stratified turbulence 
REAL(KIND=JPRB) :: RCHETB   ! Constant in the HET functions 
REAL(KIND=JPRB) :: RCHBHDL  ! Maximum znlev/L for stable boundary layer
REAL(KIND=JPRB) :: RCDHALF  ! Constant in Dyer and Hicks formulae for unstable situations
REAL(KIND=JPRB) :: RCDHPI2  ! PI/2.
REAL(KIND=JPRB) :: RLPBB    ! CONSTANT FROM THE LOUIS ET AL. FORMULATION
                            !    (simplified physics)
REAL(KIND=JPRB) :: RLPCC    ! CONSTANT FROM THE LOUIS ET AL. FORMULATION
                            !    (simplified physics)
REAL(KIND=JPRB) :: RLPDD    ! CONSTANT FROM THE LOUIS ET AL. FORMULATION
                            !    (simplified physics)

!**   ** *YOS_EXCS* CONTAINS STABILITY FUNCTION TABLES FOR *v...*

!     A.C.M. BELJAARS   E.C.M.W.F.       26/03/90.
!     P. Viterbo        E.C.M.W.F.       09/06/2005

!      NAME      TYPE        PURPOSE
!      ----      ----        -------

!     ------------------------------------------------------------------
END MODULE YOS_EXCS

