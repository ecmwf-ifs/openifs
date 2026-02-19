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

MODULE YOMGEM

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

TYPE :: TGEM

!*    * Number of grid points

!     NGPTOT   : Total number of grid columns on a PE
!     NGPTOT_CAP  : Size of grid points arrays for ALADIN
!     NGPTOTMX : Maximum number of grid columns on any of the PEs
!     NGPTOTG  : Total number of grid columns on the Globe
!     NGPTOTL(NPRGPNS,NPRGPEW)  : Total number of grid columns on on each PE

INTEGER(KIND=JPIM) :: NGPTOT
INTEGER(KIND=JPIM) :: NGPTOT_CAP
INTEGER(KIND=JPIM) :: NGPTOTMX
INTEGER(KIND=JPIM) :: NGPTOTG
INTEGER(KIND=JPIM), ALLOCATABLE :: NGPTOTL(:,:)

!     ------------------------------------------------------------------
!*    Other horizontal geometry quantities:

!     ------------------------------------------------------------------

! RDELXN  :    Horizontal mesh size at a location where the mapping factor is 1.
! SLHDP   :    Factor of the actual model resolution, i.e.
!              - for plane geometry: = sqrt( (ELX/NMSMAX)^2 + (ELY/NSMAX)^2 )
!              - for spherical geometry: = 2 Pi sqrt( a**2 / (NSMAX*(NSMAX+1)) )

REAL(KIND=JPRB) :: RDELXN
REAL(KIND=JPRB) :: SLHDP

!     ------------------------------------------------------------------

!*    * Defining the transformed sphere

!     RMUCEN : MU OF THE POLE OF STRETCHING
!     RLOCEN : LONGITUDE OF THE POLE OF STRETCHING
!     RSTRET : STRETCHING FACTOR
!     NSTTYP : 1 = POLE OF STRETCHING, POLE OF THE COLLOCATION GRID
!                 AT THE NORTHERN POLE OF THE REAL EARTH.
!              2 = THE POLE OF STRETCHING IS ANYWHERE ON THE REAL EARTH
!             AND ON THE EQUATOR OF THE COLLOCATION GRID ON THE MERIDIAN PI.
!                  THE EQUATOR OF THE COLLOCATION GRID IS TANGENT
!             TO A PARALLEL OF THE EARTH.

!     NHTYP  : 0 = regular grid
!            : 2 = number of points read on namelist namrgri

!     RNLGINC: increment to get non-linear grid

!     R4JP    inverse de delta(teta) approche a l'ordre 1
!     RC2P1   RSTRET*RSTRET+1.
!     RC2M1   RSTRET*RSTRET-1.
!     RCOR0   COMPONENT (0,0) OF CORIOLIS
!     RCOR1   COMPONENT (0,1) OF CORIOLIS

!     NLOEN(NDGSAG:NDGENG)  : number of active points on a parallel
!     NLOENG(NDGSAG:NDGENG) : global version of NLOEN
!     NMEN(NDGSAG:NDGENG)   : associated cut-off wave number
!     NMENG(NDGSAG:NDGENG)   : global version of NMEN
!     NDGLU(0:MAX(NSMAX,NMSMAX)) : number of active points in an hemisphere
!                          for a given wave number m
!     NSTAGP(NGPTOT)     : start position of latitude data for boundary fields
!     NTSTAGP(NGPTOT)    : start position of latitude data for boundary fields
!                          in C+I+E zone (ALADIN).

REAL(KIND=JPRB)    :: RMUCEN
REAL(KIND=JPRB)    :: RLOCEN
REAL(KIND=JPRB)    :: RSTRET
INTEGER(KIND=JPIM) :: NSTTYP
INTEGER(KIND=JPIM) :: NHTYP
REAL(KIND=JPRB)    :: RNLGINC
REAL(KIND=JPRB)    :: R4JP
REAL(KIND=JPRB)    :: RC2P1
REAL(KIND=JPRB)    :: RC2M1
REAL(KIND=JPRB)    :: RCOR0
REAL(KIND=JPRB)    :: RCOR1

INTEGER(KIND=JPIM), ALLOCATABLE :: NLOEN(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: NLOENG(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: NMEN(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: NMENG(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: NDGLU(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: NSTAGP(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: NTSTAGP(:)

END TYPE TGEM


!     ------------------------------------------------------------------

END MODULE YOMGEM
