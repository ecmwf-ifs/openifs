! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOMDIM

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

TYPE :: TDIM

!*    Dimensions of model working arrays

! === COLLOCATION GRID OF THE DYNAMICS ========================================

! NDGLG  : number of rows of latitudes
! NDGLL  : number of rows of latitudes for which this process is
!          performing Fourier Space calculations
! NDGNH  : number of rows in the northern hemisphere
! NDGSUR : number of additional rows at each pole for horizontal
!          interpolations.
! NDGSAG = -NDGSUR+1
! NDGSAL = Local version of NDGSAG.
! NDGSAH = 1-YRSL%NSLWIDE in DM version.
! NDGSAFPH=1-YRFP%NSLWIDE in DM version.
! NDGENG = NDGLG+NDGSUR
! NDGENL = Number of latitude rows for which this process has grid
!          point calculations to perform.
! NDGENH = NDGENL+YRSL%NSLWIDE in DM version.
! NDGENFPH=NDGENL+YRFP%NSLWIDE in DM version.
! NDGUNG : first row of the area of interest in Aladin
!        = NDGSAG in the global model
! NDGUXG : last  row of the area of interest in Aladin
!        = NDGENG in the global model
! NDGUNL : local first row in C+I zone in distributed memory Aladin
! NDGUXL : local last row in C+I zone in distributed memory Aladin
! NDLON  : length of a row of latitude near equator
! NDSUR1 : over dimensioning of NDLON for technical reasons (at least 2)
! NSTENCILWIDE : max stencil width / 2, default = 2
! NDLSUR = NDLON+NDSUR1
! NDLSM  = NDLSUR-1
! NDLUNG : first meridian of the area of interest in Aladin
!        = 1 in the global model
! NDLUXG : last  meridian of the area of interest in Aladin
!        = NDLON in the global model
! NDLUNL : local first meridian in C+I zone in distributed memory Aladin
! NDLUXL : local last meridian in C+I zone in distributed memory Aladin
! NPROMA : working dimension for grid-point computations
! NPROMM : working dimension for Meteo-France physics computations
! NPROMNH: working dimension for arrays used only in the non hydrostatic model

! NGPBLKS: number of grid point NPROMA-blocks.
! LOPTPROMA : .TRUE. NPROMA will be optimised
!           : .FALSE. NPROMA will not be optimised (forced by
!           : negative NPROMA in namelist)

INTEGER(KIND=JPIM) :: NDGLG
INTEGER(KIND=JPIM) :: NDGLL
INTEGER(KIND=JPIM) :: NDGNH
INTEGER(KIND=JPIM) :: NDGSUR
INTEGER(KIND=JPIM) :: NDGSAG
INTEGER(KIND=JPIM) :: NDGSAL
INTEGER(KIND=JPIM) :: NDGSAH
INTEGER(KIND=JPIM) :: NDGSAFPH
INTEGER(KIND=JPIM) :: NDGENG
INTEGER(KIND=JPIM) :: NDGENL
INTEGER(KIND=JPIM) :: NDGENH
INTEGER(KIND=JPIM) :: NDGENFPH
INTEGER(KIND=JPIM) :: NDGUNG
INTEGER(KIND=JPIM) :: NDGUXG
INTEGER(KIND=JPIM) :: NDGUNL
INTEGER(KIND=JPIM) :: NDGUXL
INTEGER(KIND=JPIM) :: NDLON
INTEGER(KIND=JPIM) :: NDSUR1
INTEGER(KIND=JPIM) :: NSTENCILWIDE
INTEGER(KIND=JPIM) :: NDLSUR
INTEGER(KIND=JPIM) :: NDLSM
INTEGER(KIND=JPIM) :: NDLUNG
INTEGER(KIND=JPIM) :: NDLUXG
INTEGER(KIND=JPIM), ALLOCATABLE :: NDLUNL(:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: NDLUXL(:,:)
INTEGER(KIND=JPIM) :: NPROMA
INTEGER(KIND=JPIM) :: NPROMM
INTEGER(KIND=JPIM) :: NPROMM9
INTEGER(KIND=JPIM) :: NPROMNH
INTEGER(KIND=JPIM) :: NGPBLKS
LOGICAL :: LOPTPROMA

! === SPECTRAL SPACE ==========================================================

! NRESOL  : resolution identifier
! NSMAX   : truncation order
! NMSMAX  : truncation order in longitude
! NVARMAX: truncation order in 3d-var distributed direction
!          this is a priori longitude, so that nvarmax = nsmax in Arp/IFS
!          and nvarmax = nmsmax in Aladin
! NSEFRE : number of degrees of freedom in the spectral space
! NSPECG : number of complex spectral coefficients (global)
! NSPEC2G = 2*NSPECG
! NSPEC  : number of complex spectral coefficients (local, i.e. on this PE)
! NSPEC2 = 2*NSPEC
! NSPEC2MX : maximun NSPEC2 among all PEs
! NCMAX  : upper trunc. order for dilatation matrices (used in TRAGEO, fullpos)

INTEGER(KIND=JPIM) :: NRESOL
INTEGER(KIND=JPIM) :: NSMAX
INTEGER(KIND=JPIM) :: NMSMAX
INTEGER(KIND=JPIM) :: NVARMAX
INTEGER(KIND=JPIM) :: NSEFRE
INTEGER(KIND=JPIM) :: NSPECG
INTEGER(KIND=JPIM) :: NSPEC2G
INTEGER(KIND=JPIM) :: NSPEC
INTEGER(KIND=JPIM) :: NSPEC2
INTEGER(KIND=JPIM) :: NSPEC2MX
INTEGER(KIND=JPIM) :: NCMAX

! === DISTRIBUTED MEMORY DIMENSIONS ===========================================

! NUMP  :  Number of spectral waves handled by this processor
! NUMCP :  Same as NUMP, but related to NCMAX

INTEGER(KIND=JPIM) :: NUMP
INTEGER(KIND=JPIM) :: NUMCP

END TYPE TDIM

!     ------------------------------------------------------------------

END MODULE YOMDIM
