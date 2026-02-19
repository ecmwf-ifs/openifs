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

MODULE YEMGSL

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

TYPE :: TEGSL
!*    Defining specific geometry variables for semi-Lagrangian advection in the LAM model.
!     These variables are set-up in SUEGEM_NAML and should not be modified elsewhere.
!     They are used in the ELARCHE.. routines.

!     IMPORTANT REMARK ABOUT UNITS OF ANGLES (RLON0R, RELONC, RELATC):
!      - Printings done in SUEGEM_NAML are in radians.
!      - Final values are computed in SUEGEM_NAML.
!      - When leaving SUEGEM_NAML, values SHOULD be in radians.
!      - Other routines using these values assume that they are in radians.

!     RIPORE : COPY OF XIPORE, X-location of the pole (Stereo/Lambert)
!     RJPORE : COPY OF XJPORE, Y-location of the pole (Stereo/Lambert)
!     REGPK  : COPY OF XGGPK , effective projection parameter
!     RNOBILE: QUANTITY TO SCALE THE DISTANCE TO THE POLE
!              r**2 * RNOBILE = [(1-sin(phi))/(1+sin(phi))]**REGPK  
!     RHSUD  : COPY OF HSUD                                      
!     RMERCAT: QUANTITY TO SCALE THE DISTANCE IN MERCATOR PROJECTION
!              x * RMERCAT = lon - lon0
!              y * RMERCAT = - ln(cos(phi)/(1.+sin(phi)))
!     RLON0R : ROTATED PROJECTION REFERENCE LONGITUDE (COPY OF XLON0R)                                   
!     RELONC and RELATC: longitude and latitude of the apparent pole,
!              which defines the Mercator rotated tilted projection.
!              These variables are used only for the Mercator rotated tilted projection.

REAL(KIND=JPRB) :: RIPORE
REAL(KIND=JPRB) :: RJPORE
REAL(KIND=JPRB) :: REGPK
REAL(KIND=JPRB) :: RNOBILE
REAL(KIND=JPRB) :: RHSUD
REAL(KIND=JPRB) :: RMERCAT
REAL(KIND=JPRB) :: RLON0R
REAL(KIND=JPRB) :: RELONC    
REAL(KIND=JPRB) :: RELATC 

END TYPE TEGSL

END MODULE YEMGSL
