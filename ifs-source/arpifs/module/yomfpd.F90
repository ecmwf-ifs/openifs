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

MODULE YOMFPD

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE PARFPOS  ,ONLY : JPOSDOM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Horizontal dimensions and accuracy of the horizontal (sub)domain(s)
!     All these variables are DM-global.

!========== FOR ALL KINDS OF OUTPUT (SUB)DOMAINS =========

! NLAT : number of latitudes, poles are not included if Gaussian grid.
! NLON : number of longitudes.
! RLATC : domain centre  ! IN DEGREES !!!!!
! RLONC : domain centre  ! IN DEGREES !!!!!

! RDELX and RDELY: resolution in x and y:
!  given in metres if the arrival domain is LELAM
!  given in degrees otherwise

!========== FOR ALADIN SUBDOMAIN AS OUTPUT ONLY =============

! NFPLUX: actual last row of longitude.
! NFPGUX: actual last row of latitude.

! NFPNOEXTZL : alternative extension zone (E') zonal dimension
! NFPNOEXTZG : alternative extension zone (E') meridional dimension
! NFPBZONL : half-width of relaxation zone (I) zonal dimension
! NFPBZONG : half-width of relaxation zone (I) meridional dimension

!========== FOR FRANGP0025 without missing data =====

! NFPRLX : number of rows of I for interpolation and post-processing in x dir, first point
! NFPRLY : number of rows of I for interpolation and post-processing in y dir
! NFPRUX : number of rows of I for interpolation and post-processing in x dir, last point
! NFPRUY : number of rows of I for interpolation and post-processing in y dir

!========== DEFINITION OF BOYD DOMAIN =============

!     NFPBWX     : width of Boyd extension zone in x direction
!     NFPBWY     : width of Boyd extension zone in y direction


TYPE TNAMFPD

INTEGER(KIND=JPIM) :: NLAT(JPOSDOM)
INTEGER(KIND=JPIM) :: NLON(JPOSDOM)

INTEGER(KIND=JPIM) :: NFPLUX
INTEGER(KIND=JPIM) :: NFPGUX

REAL(KIND=JPRB) :: RLATC(JPOSDOM)
REAL(KIND=JPRB) :: RLONC(JPOSDOM)
REAL(KIND=JPRB) :: RDELX(JPOSDOM)
REAL(KIND=JPRB) :: RDELY(JPOSDOM)

INTEGER(KIND=JPIM) :: NFPBWX
INTEGER(KIND=JPIM) :: NFPBWY

INTEGER(KIND=JPIM) :: NFPNOEXTZL
INTEGER(KIND=JPIM) :: NFPNOEXTZG

INTEGER(KIND=JPIM) :: NFPBZONL
INTEGER(KIND=JPIM) :: NFPBZONG

INTEGER(KIND=JPIM) :: NFPRLX
INTEGER(KIND=JPIM) :: NFPRUX
INTEGER(KIND=JPIM) :: NFPRLY
INTEGER(KIND=JPIM) :: NFPRUY

END TYPE TNAMFPD
!     ------------------------------------------------------------------
END MODULE YOMFPD
