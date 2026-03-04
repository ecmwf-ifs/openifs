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

MODULE TYPE_FPUSERGEO

USE PARKIND1  ,ONLY : JPIM     ,JPRB   ,JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

SAVE
PUBLIC TFPUSERGEO, ALLOC_FPUSERGEO, DEALLOC_FPUSERGEO, PRINT_FPUSERGEO, CHECK_FPUSERGEO

!========== USER PARAMETERS TO DEFINE A GEOMETRY  ======

!========== FOR ALL KINDS OF DOMAINS =====================================

! CFPDOM : domain label

! CFPGRID : kind of geometry : 'GAUSS', 'LELAM' or 'LALON'

! NLAT  : number of latitudes or gridpoints along y axis
! NLON  : number of longitudes or gridpoints along x axis

!========== DEFINITION OF GAUSSIAN GRID ================================

! NFPHTYP : 0 = regular grid
!         : 1 = number of points proportional to sqrt(1-mu**2)
!         : 2 = number of points read on namelist namfpg
! NFPTTYP : 1 = POLE OF STRETCHING, POLE OF THE COLLOCATION GRID AT THE NORTHERN POLE OF THE REAL EARTH.
!           2 = POLE OF STRETCHING, POLE OF THE COLLOCATION GRID ANYWHERE ON THE REAL EARTH.
! FPMUCEN : MU OF THE POLE OF INTEREST
! FPLOCEN : LONGITUDE OF THE POLE OF INTEREST
! FPSTRET : STRETCHING FACTOR
! NFPRGRI : number of active points on a parallel
! FPMU    : Sine of the latitudes (**)
! NFPMEN  : Wavenumbers on each latitude (***)

! (**) : Though this variable is temporary and could be retrieved at any time
! from the spectral transforms data (it is used only to fill the output files headers
! and the interpolations weights), it is defined in the geometry because it could be used one day
! to define a more general grid than a gaussian one. REK

! (**) : Though this variable is temporary and could be retrieved at any time from
! the spectral transforms data (it is used only to fill the output files headers), it is defined
! in the geometry because it could be used one day to define a more general spectral geometry (?). REK

!========== DEFINITION OF LIMITED AREA DOMAIN ===================================

! RDELY : horizontal resolution along latitude or y axis (in metres if LAM ; in degrees otherwise)
! RDELX : horizontal resolution along longitude or x axis (in metres if LAM ; in degrees otherwise)

! NFPRLX : number of safety rows against missing data for interpolation along x dir, lower side
! NFPRLY : number of safety rows against missing data for interpolation along y dir, lower side
! NFPRUX : number of safety rows against missing data for interpolation along x dir, upper side
! NFPRUY : number of safety rows against missing data for interpolation along y dir, upper side

!========== DEFINITION OF LAM DOMAIN ======================================

! NFPGUX : Number of gridpoints along y axis, excluding extension zone
! NFPLUX : Number of gridpoints along x axis, excluding extension zone
! NFPBWX : width of Boyd window in x direction
! NFPBWY : width of Boyd window in y direction
! RLATC  : Center of domain along y axis (in degrees)
! RLONC  : Center of domain along x axis (in degrees)
! FPLON0 : geographical longitude of reference for the projection (in degrees)
! FPLAT0 : geographical latitude of reference for the projection (in degrees)
! LFPMRT : If .TRUE. in Mercator case => use of Rotated/Tilted option
! LFPMAP : .T./.F. if the domain is defined through its coordinates/wavelengths
! FPLX   : wavelength in x (if LFPMAP=.FALSE. only)
! FPLY   : wavelength in y (if LFPMAP=.FALSE. only
! NFPBZONL : half-width of relaxation zone (I) zonal dimension (*)
! NFPBZONG : half-width of relaxation zone (I) meridional dimension (*)
! NFPEDOM  : "Kind" of geometry (spectral, gridpoint C+I+E, or C+I only) ; actually for I/Os (weird !)

! (*) Though these variables are needed today only to fill the output files headers,
! they are defined in the geometry because they could be used one day to define the width
! of a gridpoint frame. REK.

!========== SPECTRAL SPACE DEFINITIONS (Gaussian grid & LAM) ================

! NFPMAX  : spectral truncation (along y axis for LAM)
! NMFPMAX : spectral truncation (along x axis for LAM)
! NFPNOEXTZL : alternative extension zone (E') in x direction - LAM only
! NFPNOEXTZG : alternative extension zone (E') in y direction - LAM only



!========== SELF-DETERMINED VARIABLES FOR A GIVEN GEOMETRY  ======

!========== SPECTRAL SPACE DEFINITIONS ================

! LFPMODELSPEC : .TRUE. if the spectral dimensions are equal to the model ones
! NSPEC2  : local  number of spectral coefficients (MPI distribution)
! NSPEC2G : global number of spectral coefficients

!========== GRIDPOINT SPACE DEFINITIONS ================

! LFPMODELGRID : .TRUE. if the gridpoint dimensions are equal to the model ones for the whole zone (C+I+E)
! LFPMODELCORE : .TRUE. if the gridpoint dimensions are equal to the model ones for core zone only (C+I)
! LFPCOORD     : .TRUE. if actual horizontal interpolation on this grid (new gridpoint coodinates)
! LFPBIPER     : .TRUE. if actual biperiodicization on this grid (new extension zone)
! LFPLAMCOREXT : .TRUE. if LAM core extraction only
! LFPMAPF      : .TRUE. if map factor is smoothed in spectral space (and geometry is plane) for spectral outputs
! LFPOSBUFSHAPE: .TRUE. if Fullpos buffer shape should be used (otherwise model buffer shape can be used)
! NGPTOT  : local  number of grid points (MPI distribution)
! NGPTOTX : maximum local number of grid points (MPI distribution)
! NFPSIZEG: global number of grid points of the output grid
! NFPSIZEG_DEP : global number of grid points of the interpolation grid

!========== FOR ALL KINDS OF DOMAINS ====================

! NFPRESOL : resolution tag of the spectral transforms
! NFPDIST  : Kind of additional transposition between departure geometry and arrival geometry :
!  0 : no additional transposition 
!  1 : transposition based on a straightforward load balancing (can't enable spectral transforms afterwards)
!  2 : transposition based on spectral transforms package distribution


TYPE TFPUSERGEO

CHARACTER(LEN=32) :: CFPDOM = ' '

CHARACTER(LEN=5) ::  CFPGRID

INTEGER(KIND=JPIM) :: NLAT = 0
INTEGER(KIND=JPIM) :: NLON = 0

INTEGER(KIND=JPIM) :: NFPHTYP = 0
INTEGER(KIND=JPIM) :: NFPTTYP = 1
REAL(KIND=JPRB)    :: FPMUCEN = 0._JPRB
REAL(KIND=JPRB)    :: FPLOCEN = 0._JPRB
REAL(KIND=JPRB)    :: FPSTRET = 0._JPRB
INTEGER(KIND=JPIM), ALLOCATABLE :: NFPRGRI(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: NFPMEN(:)
REAL(KIND=JPRD),    ALLOCATABLE :: FPMU(:)

REAL(KIND=JPRB) :: RDELY = 0._JPRB
REAL(KIND=JPRB) :: RDELX = 0._JPRB

INTEGER(KIND=JPIM) :: NFPRLX = 0
INTEGER(KIND=JPIM) :: NFPRUX = 0
INTEGER(KIND=JPIM) :: NFPRLY = 0
INTEGER(KIND=JPIM) :: NFPRUY = 0

INTEGER(KIND=JPIM) :: NFPGUX = 0
INTEGER(KIND=JPIM) :: NFPLUX = 0
INTEGER(KIND=JPIM) :: NFPBWX = 0
INTEGER(KIND=JPIM) :: NFPBWY = 0
REAL(KIND=JPRB)    :: RLATC  = 0._JPRB
REAL(KIND=JPRB)    :: RLONC  = 0._JPRB
REAL(KIND=JPRB)    :: FPLON0 = 0._JPRB
REAL(KIND=JPRB)    :: FPLAT0 = 0._JPRB
LOGICAL            :: LFPMRT = .FALSE.
LOGICAL            :: LFPMAP = .FALSE.
REAL(KIND=JPRB)    :: FPLX   = 0._JPRB
REAL(KIND=JPRB)    :: FPLY   = 0._JPRB
INTEGER(KIND=JPIM) :: NFPBZONL = 0
INTEGER(KIND=JPIM) :: NFPBZONG = 0
INTEGER(KIND=JPIM) :: NFPEDOM  = 1

INTEGER(KIND=JPIM) :: NFPMAX  = 0
INTEGER(KIND=JPIM) :: NMFPMAX = 0
INTEGER(KIND=JPIM) :: NFPNOEXTZL = 0
INTEGER(KIND=JPIM) :: NFPNOEXTZG = 0

LOGICAL :: LFPMODELSPEC
LOGICAL :: LFPMODELGRID
LOGICAL :: LFPMODELCORE
LOGICAL :: LFPCOORD
LOGICAL :: LFPBIPER
LOGICAL :: LFPLAMCOREXT
LOGICAL :: LFPMAPF
LOGICAL :: LFPOSBUFSHAPE

INTEGER(KIND=JPIM) :: NFPDIST = 0

INTEGER(KIND=JPIM) :: NFPRESOL = 0
INTEGER(KIND=JPIM) :: NSPEC2   = 0
INTEGER(KIND=JPIM) :: NSPEC2G  = 0 
INTEGER(KIND=JPIM) :: NGPTOT   = 0 
INTEGER(KIND=JPIM) :: NGPTOTX  = 0
INTEGER(KIND=JPIM) :: NFPSIZEG = 0 
INTEGER(KIND=JPIM) :: NFPSIZEG_DEP = 0 

END TYPE TFPUSERGEO

CONTAINS

!-----------------------------------------------------------------------------

SUBROUTINE ALLOC_FPUSERGEO(YD)

TYPE(TFPUSERGEO), INTENT(INOUT) :: YD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TYPE_FPUSERGEO:ALLOC_FPUSERGEO',0,ZHOOK_HANDLE)

ALLOCATE(YD%NFPRGRI(YD%NLAT))
ALLOCATE(YD%NFPMEN(YD%NLAT))
ALLOCATE(YD%FPMU(YD%NLAT))
YD%NFPRGRI(:)=HUGE(1)
YD%NFPMEN(:)=HUGE(1)
YD%FPMU(:)=HUGE(1._JPRD)

IF (LHOOK) CALL DR_HOOK('TYPE_FPUSERGEO:ALLOC_FPUSERGEO',1,ZHOOK_HANDLE)

END SUBROUTINE ALLOC_FPUSERGEO

SUBROUTINE DEALLOC_FPUSERGEO(YD)

TYPE(TFPUSERGEO), INTENT(INOUT) :: YD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TYPE_FPUSERGEO:DEALLOC_FPUSERGEO',0,ZHOOK_HANDLE)

DEALLOCATE(YD%NFPRGRI)
DEALLOCATE(YD%NFPMEN)
DEALLOCATE(YD%FPMU)

IF (LHOOK) CALL DR_HOOK('TYPE_FPUSERGEO:DEALLOC_FPUSERGEO',1,ZHOOK_HANDLE)

END SUBROUTINE DEALLOC_FPUSERGEO

SUBROUTINE PRINT_FPUSERGEO(YD,KULOUT)

TYPE(TFPUSERGEO), INTENT(IN) :: YD
INTEGER(KIND=JPIM), INTENT(IN) :: KULOUT

INTEGER(KIND=JPIM) :: J

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TYPE_FPUSERGEO:PRINT_FPUSERGEO',0,ZHOOK_HANDLE)

WRITE(UNIT=KULOUT,FMT='('' CFPDOM = '',A9)') YD%CFPDOM
WRITE(UNIT=KULOUT,FMT='('' CFPGRID = '',A5, 5X,'' NLAT = '',I5, '' NLON = '',I5)') YD%CFPGRID, YD%NLAT, YD%NLON
WRITE(UNIT=KULOUT,FMT='('' NFPMAX = '',I5, '' NMFPMAX = '',I5)') YD%NFPMAX, YD%NMFPMAX
IF (YD%CFPGRID == 'GAUSS') THEN
  WRITE(UNIT=KULOUT,FMT='('' NFPTTYP = '',I6,'' NFPHTYP = '',I6)') &
   & YD%NFPTTYP, YD%NFPHTYP  
  WRITE(UNIT=KULOUT,FMT='('' FPMUCEN = '',E15.7,'' FPLOCEN = '',E15.7,'' FPSTRET = '',F17.14)') &
   & YD%FPMUCEN, YD%FPLOCEN, YD%FPSTRET
  WRITE(UNIT=KULOUT,FMT='('' (J,NFPRGRI) '')')
  WRITE(UNIT=KULOUT,FMT='(4(1X,''('',I5,I5,I5,E13.6,'')''))') &
   & (J,YD%NFPRGRI(J),YD%NFPMEN(J),YD%FPMU(J),J = 1,(YD%NLAT+1)/2)  
ELSE
  WRITE(UNIT=KULOUT,FMT='(''In degrees : RLONC = '',F8.2,'' RLATC = '',F8.2, &
   & '' RDELX = '',F10.2,'' RDELY = '',F10.2)') YD%RLONC,YD%RLATC,YD%RDELX,YD%RDELY
  IF (YD%NFPRLX > 0 .OR. YD%NFPRUX > 0 .OR. YD%NFPRLY > 0 .OR. YD%NFPRUY > 0) THEN
    WRITE(UNIT=KULOUT,FMT='('' NFPRLX = '',I2,'' NFPRUX = '',I2, '' NFPRLY = '',I2,'' NFPRUY = '',I2)') &
     & YD%NFPRLX, YD%NFPRUX, YD%NFPRLY, YD%NFPRUY
  ENDIF
ENDIF
IF (YD%CFPGRID == 'LELAM') THEN
  WRITE(UNIT=KULOUT,FMT='('' NFPLUX = '',I4,'' NFPGUX = '',I4, '' NFPBWX = '',I5,'' NFPBWY = '',I5)') &
   & YD%NFPLUX, YD%NFPGUX, YD%NFPBWX, YD%NFPBWY
  WRITE(UNIT=KULOUT,FMT='('' FPLON0 = '',E15.7,'' FPLAT0 = '',E15.7)') YD%FPLON0, YD%FPLAT0
  WRITE(UNIT=KULOUT,FMT='('' LFPMRT = '',L2,'' LFPMAP = '',L2, '' NFPEDOM = '', I2)') YD%LFPMRT, YD%LFPMAP, YD%NFPEDOM
  IF (.NOT.YD%LFPMAP) THEN
    WRITE(UNIT=KULOUT,FMT='('' FPLX = '',E15.7,'' FPLY = '',E15.7)') YD%FPLX, YD%FPLY
  ENDIF
  WRITE(UNIT=KULOUT,FMT='('' NFPBZONL = '',I5,'' NFPBZONG = '',I5,'' NFPNOEXTZL = '',I5,'' NFPNOEXTZG = '',I5)') &
   & YD%NFPBZONL, YD%NFPBZONG, YD%NFPNOEXTZL, YD%NFPNOEXTZG
ENDIF

WRITE(KULOUT,'('' LFPMODELSPEC = '', L2, '' LFPMODELGRID = '',L2,'' LFPMODELCORE = '',L2,'' LFPCOORD = '',L2, &
 & '' LFPBIPER = '',L2, '' LFPLAMCOREXT = '',L2,'' LFPMAPF = '',L2,'' LFPOSBUFSHAPE = '',L2,'' NFPDIST = '', I2)') &
 & YD%LFPMODELSPEC, YD%LFPMODELGRID, YD%LFPMODELCORE, YD%LFPCOORD, YD%LFPBIPER, YD%LFPLAMCOREXT, YD%LFPMAPF, YD%LFPOSBUFSHAPE, &
 & YD%NFPDIST
WRITE(KULOUT,'('' NFPSIZEG = '',I9,'' NFPSIZEG_DEP = '',I9)') YD%NFPSIZEG, YD%NFPSIZEG_DEP

IF (LHOOK) CALL DR_HOOK('TYPE_FPUSERGEO:PRINT_FPUSERGEO',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_FPUSERGEO

SUBROUTINE CHECK_FPUSERGEO(YD,KULOUT)

TYPE(TFPUSERGEO), INTENT(INOUT) :: YD
INTEGER(KIND=JPIM), INTENT(IN) :: KULOUT

INTEGER(KIND=JPIM) :: IERR, J
REAL(KIND=JPRB) :: ZEPS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('TYPE_FPUSERGEO:CHECK_FPUSERGEO',0,ZHOOK_HANDLE)

! Finalize :
IF (YD%NFPTTYP == 1) THEN
  YD%FPMUCEN = 1._JPRB
  YD%FPLOCEN = 0._JPRB
ENDIF
IF (YD%CFPGRID == 'LELAM') THEN
  IF (YD%LFPMRT) THEN
    YD%FPLAT0 = 0._JPRB
  ENDIF
  YD%NFPBZONL=MAX(0,MIN(YD%NFPBZONL,(YD%NFPLUX-1)/2))
  YD%NFPBZONG=MAX(0,MIN(YD%NFPBZONG,(YD%NFPGUX-1)/2))
ELSE
  YD%NFPLUX=0
  YD%NFPGUX=0
  YD%FPLON0 = 0._JPRB
  YD%FPLAT0 = 0._JPRB
  YD%LFPMRT = .FALSE.
  YD%LFPMAP =.TRUE.
  YD%NFPBWX=0
  YD%NFPBWY=0
ENDIF
IF (YD%CFPGRID == 'GAUSS') THEN
  YD%RLONC = YD%FPLOCEN
  YD%RLATC = YD%FPMUCEN
  YD%NFPTTYP=MAX(1,MIN(2,YD%NFPTTYP))
  YD%NFPHTYP=MAX(0,MIN(2,YD%NFPHTYP))
ELSE
  YD%NFPHTYP = 0
  YD%NFPTTYP = 1
  YD%FPSTRET = 1._JPRB
  YD%FPLX=YD%RDELX*REAL(YD%NLON,JPRB)
  YD%FPLY=YD%RDELY*REAL(YD%NLAT,JPRB)
ENDIF

IF (YD%NFPHTYP == 0) THEN
  YD%NFPRGRI(1:YD%NLAT) = YD%NLON
ELSE
  ! Control the range of values
  DO J = 1, (YD%NLAT+1)/2
    YD%NFPRGRI(J) = MIN(YD%NLON,MAX(1,YD%NFPRGRI(J)))
  ENDDO
  ! Ensure symetry of the grid
  DO J = 1, YD%NLAT/2
    YD%NFPRGRI(YD%NLAT+1-J) = YD%NFPRGRI(J)
  ENDDO
ENDIF

YD%NFPMAX=MAX(0,YD%NFPMAX)
YD%NMFPMAX=MAX(0,YD%NMFPMAX)

! Control validity of variables :
ZEPS=EPSILON(1.0_JPRB)*100._JPRB
IERR=0
IF (YD%NFPHTYP==1) THEN
  IERR=IERR+1
  WRITE(UNIT=KULOUT,FMT='('' NFPHTYP = 1 IS NOT SUPPORTED ANYMORE !'')')
ENDIF
IF (YD%CFPGRID /= 'GAUSS') THEN
  IF ((YD%RLATC < (-90._JPRB-ZEPS)).OR.(YD%RLATC > (90._JPRB+ZEPS))) THEN
    IERR=IERR+1
    WRITE(UNIT=KULOUT,FMT='('' RLATC('',I2,'') = '',F6.1,'' OUT OF [-90.,90.]'')') YD%RLATC  
  ENDIF
  IF ((YD%RLONC < (-360._JPRB-ZEPS)).OR.(YD%RLONC > (360._JPRB+ZEPS)))THEN
    IERR=IERR+1
    WRITE(UNIT=KULOUT,FMT='('' RLONC('',I2,'') = '',F6.1,'' OUT OF [-360.,360.]'')') YD%RLONC  
  ENDIF
ELSE
  IF ((YD%FPSTRET-1.0_JPRB) < -ZEPS) THEN
    IERR = IERR+1
    WRITE(UNIT=KULOUT,FMT='('' FPSTRET = '',F17.14,'' < 1. FORBIDDEN, CHANGE THE POLE !'')') YD%FPSTRET
  ENDIF
ENDIF
IF (YD%NLAT < 1) THEN
  IERR=IERR+1
  WRITE(UNIT=KULOUT,FMT='('' NLAT = '',I5,'' LESS THAN 1 !'')') YD%NLAT  
ELSE
  IF (YD%CFPGRID == 'LELAM' .AND. YD%NFPGUX > YD%NLAT) THEN
    IERR=IERR+1
    WRITE(UNIT=KULOUT,FMT='(''ERROR : NFPGUX > NLAT '')') YD%NFPGUX, YD%NLAT
  ENDIF
ENDIF
IF (YD%NLON < 1) THEN
  IERR=IERR+1
  WRITE(UNIT=KULOUT,FMT='('' NLON = '',I5,'' LESS THAN 1 !'')') YD%NLON  
ELSE
  IF (YD%CFPGRID == 'LELAM' .AND. YD%NFPLUX > YD%NLON) THEN
    IERR=IERR+1
    WRITE(UNIT=KULOUT,FMT='(''ERROR : NFPLUX > NLON '')') YD%NFPLUX, YD%NLON
  ENDIF
ENDIF
IF (IERR > 0) CALL ABOR1('TYPE_FPUSERGEO:CHECK_FPUSERGEO : ABOR1 CALLED')

IF (LHOOK) CALL DR_HOOK('TYPE_FPUSERGEO:CHECK_FPUSERGEO',1,ZHOOK_HANDLE)

END SUBROUTINE CHECK_FPUSERGEO

!     ------------------------------------------------------------------
END MODULE TYPE_FPUSERGEO
