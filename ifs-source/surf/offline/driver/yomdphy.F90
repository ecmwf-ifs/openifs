MODULE YOMDPHY
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE ISO_C_BINDING
IMPLICIT NONE
SAVE

!     ------------------------------------------------------------------

!     DIMENSION DES TABLEAUX POINT DE GRILLE PHYSIQUES

!     NGP  : total number of fields in the physics array
!     NGPP : number of fields in the prognostic physics array
!     NGPA : number of fields in accumulated prognostics array
!     NVOC : number of variables in the generic OCEAN.
!     NVSO : number of variables in the generic SOILB.
!     NVEG : number of variables in the generic VEGET.
!     NVSC : number of variables in the generic SNOWC.
!     NVSG : number of variables in the generic SNOWG.
!     NVRS : number of variables in the generic RESVR.
!     NVWV : number of variables in the generic WAVES.
!     NVXP : number of variables in the generic EXTRP.
!     NVXP2: number of variables in the generic XTRP2.
!     NVCAN: number of variables in the generic CANRI (specific for CANARI OI).
!     NVLA : number of variables in the generic FLAKEG  ! ENDUTRA

!     NSOTY: number of soil types
!     NCOC : number of levels in OCEAN
!     NCOM : number of levels in OCEAN MIXED LAYER   !KPP
!     NCSS : number of levels in SOILB
!     NCSV : number of levels in VEGET
!     NCWV : number of levels in WAVES
!     NCXP : number of levels in EXTRP
!     NCSI : number of sea-ice levels
!     NCSNEC: number of snow levels in EC physics (specific for ECMWF)
!     NTILES: number of surface tiles
!     NVHILO: number of vegetation classes
!     NVTYPES: number of land use types (including vegetation)
!     NTRAC: number of tracers

!     NVPD : number of variables in the generic GPD.
!     NVSF : number of variables in the generic VARSF.
!     NVCLIP : number of variables in the generic VCLIP
!     NVCLIV : number of variables in the generic VCLIV
!     NVCLIN : number of variables in the generic VCLIN
!     NVCLIH : number of variables in the generic VCLIH
!     NVCLIA : number of variables in the generic VCLIA
!     NVCLIG : number of variables in the generic VCLIG
!     NVRADF : number of variables in the generic VRADF
!     NVEXTR : number of variables in the generic VEXTRA (extra-fields)
!     NVEXTRDYN : number of extra-fields comming from the dynamics
!     NVXTR2 : number of variables in the generic VEXTR2
!     NVDIAG : number of variables in the generic VDIAG
!     NVCLIX : number of variables in the generic VCLIX
!     NVCLIU : number of variables in the generic VCLIU
!     NVO3ABC: number of variables in the generic VO3ABC
!     NVRMOON: number of variables in the generic VRMOON

!     NCRADF : number of levels in the generic VRADF
!     NCEXTR : number of levels in the generic VEXTRA

!     NTSL : nombre de types de sols nus.
!     NTSV : nombre de types de vegetations.
!     NTOZ1D: 1 si representation 1D des NVCLIS variables , 0 sinon
!     NTOZ2D: 1 si representation 2D des NVCLIS variables , 0 sinon
!     NTOZ3D: 1 si representation 3D des NVCLIS variables , 0 sinon
!     NLOA : nombre de longueurs d'ondes pour le spectre d'albedo.
!     NLOE : nombre de longueurs d'ondes pour le spectre d'emissivite.
!     NTSSG : number of surface temperatures for subgrid diagnostics

!     NCHAC: nombre de champs a accumuler.
!     NCHIN: nombre de champs instantanes.
!     NSIRA: nombre d'intervales spectraux pour les diagnostiques de
!            flux radiatif a p=0 et p=ps.

!     NPH0 : 1 ==> physic's fields at time t    allocated, 0 ==> not
!     NPH9 : 1 ==> physic's fields at time t-dt allocated  0 ==> not
!     NPHS : 1 ==> physic's climatalogical flds allocated  0 ==> not

!     NVTEND: number of tendencies used in 2.order scheme

!     NGCC : number of variables in the generic GPCC
!                  (seasonally varying fields)
!     NCLIMDAT : number of seasonally varying data per year ( 12 months, 36 decades)


!     NVELO  : dimension of velocity of ocean mixed layer model !KPP
!     NSCLRO : dimension of scalars of ocean mixed layer model  !KPP

!     NPOI : number of gridpoints not masked out
!     NLON : number of longitudes
!     NLAT : number of latitudes
!     NLALO: number of gridpoints (NLAT*NLON)
!     NLEV : number of vertical levels
!     NSLAB: max number of grid points per slab (when NDIMCDF = 1)    ! DEPRECATED
!            or max nr of latitude lines per slab (when NDIMCDF = 2)  ! DEPRECATED
!     ISLAB: Current slab number                                      ! DEPRECATED
!     NCOOR: Selected grid point index
!     NPDONE:Number of points done so far

!     DZLAT: array with latitudes
!     DZLON: array with longitudes

INTEGER(KIND=JPIM) :: NGP
INTEGER(KIND=JPIM) :: NGPP
INTEGER(KIND=JPIM) :: NGPA
INTEGER(KIND=JPIM) :: NVPD
INTEGER(KIND=JPIM) :: NGCC
INTEGER(KIND=JPIM) :: NVOC
INTEGER(KIND=JPIM) :: NVSO
INTEGER(KIND=JPIM) :: NVEG
INTEGER(KIND=JPIM) :: NVSC
INTEGER(KIND=JPIM) :: NVSG
INTEGER(KIND=JPIM) :: NVRS
INTEGER(KIND=JPIM) :: NVWV
INTEGER(KIND=JPIM) :: NVXP
INTEGER(KIND=JPIM) :: NVXP2
INTEGER(KIND=JPIM) :: NVCAN
INTEGER(KIND=JPIM) :: NCOC
INTEGER(KIND=JPIM) :: NCOM
INTEGER(KIND=JPIM) :: NCSS
INTEGER(KIND=JPIM) :: NSOTY
INTEGER(KIND=JPIM) :: NCSV
INTEGER(KIND=JPIM) :: NCWV
INTEGER(KIND=JPIM) :: NCXP
INTEGER(KIND=JPIM) :: NCSI
INTEGER(KIND=JPIM) :: NCSNEC
INTEGER(KIND=JPIM) :: NSCLRO !KPP
INTEGER(KIND=JPIM) :: NTILES
INTEGER(KIND=JPIM) :: NVHILO !CTESSEL
INTEGER(KIND=JPIM) :: NVTYPES !CTESSEL
INTEGER(KIND=JPIM) :: NTRAC
INTEGER(KIND=JPIM) :: NTSL
INTEGER(KIND=JPIM) :: NVSF
INTEGER(KIND=JPIM) :: NVCLIP
INTEGER(KIND=JPIM) :: NVCLIV
INTEGER(KIND=JPIM) :: NVCLIN
INTEGER(KIND=JPIM) :: NVCLIH
INTEGER(KIND=JPIM) :: NVCLIA
INTEGER(KIND=JPIM) :: NVCLIG
INTEGER(KIND=JPIM) :: NVCLIU
INTEGER(KIND=JPIM) :: NVELO  !KPP
INTEGER(KIND=JPIM) :: NVO3ABC
INTEGER(KIND=JPIM) :: NVRMOON
INTEGER(KIND=JPIM) :: NVRADF
INTEGER(KIND=JPIM) :: NVEXTR
INTEGER(KIND=JPIM) :: NVEXTRDYN
INTEGER(KIND=JPIM) :: NVXTR2
INTEGER(KIND=JPIM) :: NVDIAG
INTEGER(KIND=JPIM) :: NCRADF
INTEGER(KIND=JPIM) :: NCEXTR
INTEGER(KIND=JPIM) :: NVCLIX
INTEGER(KIND=JPIM) :: NVCLIS
INTEGER(KIND=JPIM) :: NTOZ1D
INTEGER(KIND=JPIM) :: NTOZ2D
INTEGER(KIND=JPIM) :: NTOZ3D
INTEGER(KIND=JPIM) :: NTSSG
INTEGER(KIND=JPIM) :: NTVG
INTEGER(KIND=JPIM) :: NLOA
INTEGER(KIND=JPIM) :: NLOE
INTEGER(KIND=JPIM) :: NCHAC
INTEGER(KIND=JPIM) :: NCHIN
INTEGER(KIND=JPIM) :: NSIRA
INTEGER(KIND=JPIM) :: NPH0
INTEGER(KIND=JPIM) :: NPH9
INTEGER(KIND=JPIM) :: NPHS
INTEGER(KIND=JPIM) :: NVTEND

INTEGER(KIND=JPIM) :: NPOI
INTEGER(KIND=JPIM) :: NLON
INTEGER(KIND=JPIM) :: NLAT
INTEGER(KIND=JPIM) :: NLEV
INTEGER(KIND=JPIM) :: NLALO
INTEGER(KIND=JPIM) :: NCOOR
INTEGER(KIND=JPIM) :: NPDONE
INTEGER(KIND=JPIM) :: NCLIMDAT 
     
REAL(KIND=JPRB),ALLOCATABLE:: DZLAT(:)
REAL(KIND=JPRB),ALLOCATABLE:: DZLON(:)

INTEGER(KIND=JPIM) :: NVLA ! ENDUTRA

TYPE(C_PTR)      :: YDSURF

INTEGER(KIND=JPIM) :: NPOIG
INTEGER(KIND=JPIM),ALLOCATABLE :: NPOIOFF(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: NPOIP(:)
!     ------------------------------------------------------------------
END MODULE YOMDPHY
