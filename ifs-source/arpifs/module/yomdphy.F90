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
MODULE YOMDPHY

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

TYPE :: TDPHY

!     ------------------------------------------------------------------

!     DIMENSION DES TABLEAUX POINT DE GRILLE PHYSIQUE

!     NCSS : number of deep soil layers

!     NVXP : number of variables in the generic EXTRP.
!     NVXP2: number of variables in the generic XTRP2.

!     NCXP : number of levels in EXTRP
!     NCSI : number of sea-ice levels
!     NCOM : number of ocean mix-layer model levels
!     NCSNEC: number of snow levels in EC physics (specific for ECMWF) 
!     NTILES: number of surface tiles

!     NVEXTR : number of variables in the generic VEXTRA (extra-fields)
!     NVEXTRRAD : number of variables in the radiation VEXTRA (extra-fields)
!     NVEXTRDYN : number of extra-fields comming from the dynamics
!     NVXTR2 : number of variables in the generic VEXTR2

!     NVECOUT : number of model variable put out via passive DDH pathway

!     NCEXTR : number of levels in the generic VEXTRA

!     NTSV : nombre de types de vegetations.
!     NTOZ1D: 1 si representation 1D des NVCLIS variables , 0 sinon
!     NTOZ2D: 1 si representation 2D des NVCLIS variables , 0 sinon
!     NTOZ3D: 1 si representation 3D des NVCLIS variables , 0 sinon
!     NTSSG : number of surface temperatures for subgrid diagnostics

!     LTPROF: .T. if more than 1 vertical layer in deep soil
!     LDIRCLSMOD: L to take 2m model equivalent directly from input file

INTEGER(KIND=JPIM) :: NCSS
INTEGER(KIND=JPIM) :: NVXP
INTEGER(KIND=JPIM) :: NVXP2
INTEGER(KIND=JPIM) :: NCXP
INTEGER(KIND=JPIM) :: NCSI
INTEGER(KIND=JPIM) :: NCOM
INTEGER(KIND=JPIM) :: NCSNEC
INTEGER(KIND=JPIM) :: NTILES
INTEGER(KIND=JPIM) :: NVEXTR
INTEGER(KIND=JPIM) :: NVEXTRDI
INTEGER(KIND=JPIM) :: NVEXTRRAD
INTEGER(KIND=JPIM) :: NVEXTRDYN
INTEGER(KIND=JPIM) :: NVXTR2
INTEGER(KIND=JPIM) :: NVECOUT
INTEGER(KIND=JPIM) :: NCEXTR
INTEGER(KIND=JPIM) :: NVCLIS
INTEGER(KIND=JPIM) :: NTOZ1D
INTEGER(KIND=JPIM) :: NTOZ2D
INTEGER(KIND=JPIM) :: NTOZ3D
INTEGER(KIND=JPIM) :: NTSSG
LOGICAL            :: LTPROF
LOGICAL            :: LDIRCLSMOD
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TDPHY
!     --------------------------------------------------------------------------------

TYPE(TDPHY), POINTER :: YRDPHY => NULL()

!     ------------------------------------------------------------------
CONTAINS 
  
SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TDPHY), INTENT(IN) :: SELF
  INTEGER     , INTENT(IN) :: KDEPTH
  INTEGER     , INTENT(IN) :: KOUTNO

  INTEGER :: IDEPTHLOC

  IDEPTHLOC = KDEPTH+2
  
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_g%yrdphy : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NCSS = ', SELF%NCSS
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVXP = ', SELF%NVXP
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVXP2 = ', SELF%NVXP2
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NCXP = ', SELF%NCXP
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NCSI = ', SELF%NCSI
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NCOM = ', SELF%NCOM
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NCSNEC = ', SELF%NCSNEC
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NTILES = ', SELF%NTILES
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVEXTR = ', SELF%NVEXTR
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVEXTRDI = ', SELF%NVEXTRDI
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVEXTRRAD = ', SELF%NVEXTRRAD
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVEXTRDYN = ', SELF%NVEXTRDYN
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVXTR2 = ', SELF%NVXTR2
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVECOUT = ', SELF%NVECOUT
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NCEXTR = ', SELF%NCEXTR
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVCLIS = ', SELF%NVCLIS
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NTOZ1D = ', SELF%NTOZ1D
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NTOZ2D = ', SELF%NTOZ2D
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NTOZ3D = ', SELF%NTOZ3D
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NTSSG = ', SELF%NTSSG
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LTPROF = ', SELF%LTPROF
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LDIRCLSMOD = ', SELF%LDIRCLSMOD

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOMDPHY
