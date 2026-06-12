MODULE YOS_MLM

USE PARKIND1  ,ONLY : JPIM     ,JPRB

! (C) Copyright 2011- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

IMPLICIT NONE

SAVE

!
!*    ** *YOS_MLM* - VARIABLES FOR APPLYING TKE BASED MIXED 
!                    LAYER MODEL


!     P. JANSSEN     E.C.M.W.F.       2/11/11  

!      NAME      TYPE       PURPOSE
!      ----      ----       -------
!      LOCMLTKE  LOGICAL    IF TRUE TKE BASED MIXED LAYER MODEL WILL BE SWITCHED ON
!      LLINLOG   LOGICAL    IF TRUE A LINEAR VERTICAL GRID IS USED, OTHERWISE
!                           A LOGARITHMIC GRID IS USED (THIS IS THE PREFERRED OPTION
!                           BECAUSE DETAILED FEATURES NEAR THE SURFACE ARE BETTER 
!                           TREATED).
!      SQ        REAL       MELLOR-YAMADA CONSTANT FOR TKE TRANSPORT
!      SM        REAL       MELLOR-YAMADA CONSTANT FOR MOMENTUM TRANSPORT
!      B         REAL       MELLOR-YAMADA CONSTANT FOR DISSIPATION OF TKE
!      XKAPPA    REAL       VON KARMAN CONSTANT
!      RHO       REAL       WATER DENSITY
!      CP        REAL       HEAT CAPACITY OF WATER
!      ZS        REAL       CONSTANT TO GENERATE VERTICAL GRID 
!      ZR1       REAL       PENETRATION DEPTH SHORT WAVE RADIATION
!      ZR2       REAL       PENETRATION DEPTH MEDIUM WAVE RADIATION
!      ZR3       REAL       PENETRATION DEPTH LONG WAVE RADIATION
!      R1        REAL       RELATIVE SHORT WAVE CONTRIBUTION
!      R2        REAL       RELATIVE MEDIUM WAVE CONTRIBUTION
!      R3        REAL       RELATIVE LONG WAVE CONTRIBUTION
!      R3        REAL       RELATIVE LONG WAVE CONTRIBUTION
!      DML       REAL       DEPTH OF MIXED LAYER
!
!      R         REAL ARRAY SOLAR ABSOPTION PROFILE
!      F         REAL ARRAY STRESS PROFILE
!      Z         REAL ARRAY DEPTH OF A LAYER
!      DELZ      REAL ARRAY THICKNESS OF THE LAYER
!
!-----------------------------------------------------------------------
!

TYPE :: TMLM
LOGICAL :: LOCMLTKE      
LOGICAL :: LLINLOG      
REAL(KIND=JPRB) :: SQ
REAL(KIND=JPRB) :: SM
REAL(KIND=JPRB) :: B
REAL(KIND=JPRB) :: XKAPPA
REAL(KIND=JPRB) :: RHO
REAL(KIND=JPRB) :: CP
REAL(KIND=JPRB) :: ZS
REAL(KIND=JPRB) :: ZR1
REAL(KIND=JPRB) :: ZR2
REAL(KIND=JPRB) :: ZR3
REAL(KIND=JPRB) :: R1
REAL(KIND=JPRB) :: R2
REAL(KIND=JPRB) :: R3
REAL(KIND=JPRB) :: DML

REAL(KIND=JPRB), ALLOCATABLE :: R(:)
REAL(KIND=JPRB), ALLOCATABLE :: F(:)
REAL(KIND=JPRB), ALLOCATABLE :: Z(:)
REAL(KIND=JPRB), ALLOCATABLE :: DELZ(:)
END TYPE TMLM

!     ------------------------------------------------------------------
END MODULE YOS_MLM
