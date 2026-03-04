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

MODULE YOMWFPB

! Horizontal interpolations : weights management

!       Numbering of neighbouring points on the input grid :

!                     13       5       6

!          ILAT ->     7       1       2       8
!                                  x
!                      9       3       4      10

!                     15      11      12

! FOR ALL :
! =======

! Indexes and weights for horizontal interpolations:
! ML0   : index of the four western points of the 16 points interpolation grid (resp. 13,7,9,15) in a cache-blocked buffer
! for surfex
! MS1   : indexes of a column from north to south (MS1 is supposed to contains ML0, but extends over the whole halo height).

! LWSTD04 : .TRUE. if any bilinear interpolation (*) without any mask
! LWSTD12 : .TRUE. if any quadratic interpolation without any mask

! WSTD04 : weights for bilinear interpolations, without mask
! WSTD12 : weights for quadratic interpolations, without mask
! WSTDML : additional weights for multi-linear interpolations, without mask
! RINVGMS: inverse of input model geographic resolution


! FOR SURFACE PHYSICAL FIELDS ONLY :
! ================================

! LWLAN04 : .TRUE. if any bilinear interpolation (*) with land mask
! LWLAN12 : .TRUE. if any quadratic interpolation with land mask

! WLAN04 : weights for bilinear interpolations, with land-sea mask
! WLAN12 : weights for quadratic interpolations, with land-sea mask
! WLANML : additional weights for multi-linear interpolations, with land-sea mask

! LWSEA04 : .TRUE. if any bilinear interpolation (*) with sea mask
! LWSEA12 : .TRUE. if any quadratic interpolation with sea mask 

! WSEA04 : weights for bilinear interpolations, with sea mask
! WSEA12 : weights for quadratic interpolations, with sea mask
! WSEAML : additional weights for multi-linear interpolations, with sea mask

! MWIC   : indicator of a created isolated lake or island (target grid)
!         =0 : at least one surrounding point of the same nature
!         =1 : isolated lake point
!         =2 : isolated island point


! FOR SURFEX ONLY :
! ===============

! WSFXPS : Weights for PREP fields

! LPSLSFX/NPSLSFX indices for validity mask selection

! SFXMSK : Masks (target grid) for PREP fields

!------------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

TYPE TFPWSTD

! to be re-worked. REK

LOGICAL :: LWSTD04
LOGICAL :: LWSTD12

REAL(KIND=JPRB), ALLOCATABLE :: WSTD04(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: WSTD12(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: WSTDML(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: RINVGMS(:)

INTEGER(KIND=JPIM), ALLOCATABLE :: ML0(:,:,:)
! for surfex 
INTEGER(KIND=JPIM), ALLOCATABLE :: MS1(:,:,:)

END TYPE TFPWSTD


TYPE TFPSUW

! Weights for ISBA or TESSEL fields :

LOGICAL :: LWLAN04
LOGICAL :: LWLAN12

REAL(KIND=JPRB), ALLOCATABLE :: WLAN04(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: WLAN12(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: WLANML(:,:,:)

LOGICAL :: LWSEA04
LOGICAL :: LWSEA12

REAL(KIND=JPRB), ALLOCATABLE :: WSEA04(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: WSEA12(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: WSEAML(:,:,:)

INTEGER(KIND=JPIM), ALLOCATABLE :: MWIC(:,:)

! Weights for SURFEX fields :

LOGICAL, ALLOCATABLE :: LPSLSFX (:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: NPSLSFX (:,:,:,:)

REAL(KIND=JPRB), ALLOCATABLE :: WSFXPS(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SFXMSK(:,:,:)

END TYPE TFPSUW

!------------------------------------------------------------------------------

END MODULE YOMWFPB
