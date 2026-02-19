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

MODULE YOERRTWN

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

INTEGER(KIND=JPIM) , DIMENSION(16) :: NG
INTEGER(KIND=JPIM) , DIMENSION(16) :: NSPA
INTEGER(KIND=JPIM) , DIMENSION(16) :: NSPB

REAL(KIND=JPRB) , DIMENSION(16) :: DELWAVE

REAL(KIND=JPRB) , DIMENSION(181,16) :: TOTPLNK

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15

!  NAME     TYPE     PURPOSE
!  ----   : ----    : -------
!  NG     : INTEGER : Number of k-coefficients in spectral intervals
!  NSPA   : INTEGER :
!  NSPB   : INTEGER :
! WAVENUM1: REAL    : Lower wavenumber spectral limit
! WAVENUM2: REAL    : Higher wavenumber spectral limit
! DELWAVE : REAL    : Spectral interval width
! TOTPLNK : REAL    :
!     -----------------------------------------------------------------
END MODULE YOERRTWN
