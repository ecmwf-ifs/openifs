! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

MODULE PARDIM1C

!     PURPOSE.
!     --------
!     Dimensions of SCM specific arrays.  PARDIM1C exists
!     additionally to PARDIM.

!     Modifications.
!     --------------
!     ORIGINAL SCM
!     Modified 2001-02 M.Koehler CY23R4
!     Modified 2006-01 M.Koehler CY30R1 
!             (includes large arrays for 2400 level runs)

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM

IMPLICIT NONE

SAVE

!**   for yomgt1c* (atmospheric variables)

INTEGER(KIND=JPIM), PARAMETER :: JPGT=10000
INTEGER(KIND=JPIM), PARAMETER :: JPFLEV=3000

!*    for yomgp1c* (surface variables)

INTEGER(KIND=JPIM), PARAMETER :: JPGPP=1000

!*    for yomgpd1c (diagnostic variables)

INTEGER(KIND=JPIM), PARAMETER :: JPVPD=2000
INTEGER(KIND=JPIM), PARAMETER :: JPCEXTR=3000
! Maximum number of tiles
INTEGER(KIND=JPIM), PARAMETER :: NTILESMX=9  ! 100

!*    for yomgf1c (large scale forcing)

INTEGER(KIND=JPIM), PARAMETER :: JPGF=10000

!     ------------------------------------------------------------------

END MODULE PARDIM1C
