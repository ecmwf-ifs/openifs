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

MODULE GOM_PLUS
!
! Alan Geer   1-June-2015     ECMWF
!
! Purpose.
! --------
!
! Maintains and fills a data structure containing the model fields interpolated to
! observation locations, arranged in KSET order, and with additional re-generated
! model fields such as geopotential, full and half pressure. This makes
! extensive use of forecast model code like GPGEO - now we need to work out how
! to OOPSify this!
!
! Replaces code previously known as "the preint routines".
!
! AJGDB could be made more automatic in the future (using the "model fields" concept
! of a few years ago) but for now it is quite hand-coded.
!
! AJGDB2 we might aspire to change the variable names to something more meaningful
! and/or consistent with the names in the GOMS, but for now these names are, to
! prevent confusion, identical to the old variables used around hop and elsewhere
!
! (GOM**6) Copy GOM variables from/to hop*.F90 "packets"
! 1) First follow the instructions at the top of gom_mod.F90 to add the new variable in the GOMs.
! 2) Add a pointer of the appropriate name, shape and type in the definition of type_gom_plus
! 3) In gom_plus_alloc, point the new variable into a storage array of the appropriate shape, at the next available index
! 4) Also in gom_plus_alloc, increase the size of the storage area as appropriate (e.g. insfc, inua1 etc.)
! 5) Nullify the new pointer in gom_plus_destroy
! 6) Extract the variable from the GOMs in gom_plus_fill
! 7) Put back to the GOMs in gom_plus_fill_ad
!
! Possibly out of date but it might help: (see gom_mod for up-to-date accurate info)
!
!        WHERE KDLEN   = FIRST DIMENSION (OBSERVATIONS)
!              KSET    = NUMBER (NUMERO) OF THE OBS SET TO BE TREATED
!              POROG   = OROGRAPHY AT OBS. PONTS
!              PLS5    = Land/Sea Mask
!              PCI5    = Sea-ice fraction
!              PAL5    = Albedo
!              PTS5    = Surface skin temperature
!              PWS5    = Soil wetness
!              PSN5    = Snow depth
!              PU10N5  = U-COMP OF 10M NEUTRAL WIND
!              PV10N5  = V-COMP OF 10M NEUTRAL WIND
!              PSST5   = SURFACE TEMPERATURE OVER SEA (duplicated obs)
!              PU05    = U-COMP OF OCEAN CURRENT
!              PV05    = V-COMP OF OCEAN CURRENT
!              PCHN5   = WAVE CHARNOCK
!              PRF5    = R
!              PVEG5   = VEGETATION
!              QN      = Negative values of humidity set to zero (AJGDB - commonality with hop call to qneglim?)
!              PWL5    = SKIN RESERVOIR
!              PZ05    = ROUGHNESS LENGTH
!              PCPTGZ5 = STATIC ENERGY (LOWER LEVEL)
!              PWET5   = WETNESS
!              PQS5    = SATURATION HUMIDITY
!              PCPTS5  = STATIC ENERGY AT THE SURFACE
!              PRI5    = RICHARDSON NO.
!              PBN5    = STABILITY FUNCTION (NEUTRAL)
!              PBM5    = DRAG COEFFICIENT FOR MOMENTUM
!              PBH5    = DRAG COEFFICIENT FOR HEAT
!              PLNPR5  = LOG. OF THE RATIO OF PRESSURE
!              PALPH5  = COEFFICIENTS OF THE HYDROSTATICS
!              PTCLS5    = CLS TEMPERATURE
!              PRESH5  = HALF LEVEL PRESS. AT OBS. POINTS
!              PRESF5  = FULL LEVEL PRESS. AT OBS. POINTS
!              PGEOPH5 = HALF LEVEL GEOP. AT OBS. POINTS
!              PGEOPF5 = FULL LEVEL GEOP. AT OBS. POINTS
!              PAPHI   = SOME SORT OF MUTANT GEOPOTENTIAL USED BY SOME SURFACE OPERATORS
!              PXP5    = INTERMEDIATE ARRAY (FOR INTERPOLATION IN P)
!              PXPD5   = ...........................................
!              PXZ5    = INTERMEDIATE ARRAY (FOR INTERPOLATION IN Z)
!              PXZD5   = ...........................................
!              PCP5    = Cp
!              PRF5    = R
!              PKAP5   = Kappa whatever that is
!              PUF5    = U VALUES (MODEL) AT OB. POINS
!              PVF5    = V VALUES (MODEL) AT OB. POINS
!              PTF5    = T VALUES (MODEL) AT OB. POINS
!              PQF5    = Q VALUES (MODEL) AT OB. POINS
!              PO3F5   = O3 VALUES (MODEL) AT OB. POINS
!              PCLWF5  = CLW VALUES (MODEL) AT OB. POINS
!              PCIWF5  = CIW VALUES (MODEL) AT OB. POINS
!              PCCF5   = CC VALUES (MODEL) AT OB. POINS
!              PCSWF5    = S (snow) VALUES (MODEL) AT OB. POINS
!              PCRWF5    = R (rain) VALUES (MODEL) AT OB. POINS
!              PCGWF5    = G (graupel) VALUES (MODEL) AT OB. POINS
!              PCLWD5  = Diagnostic CLW from callpar
!              PCIWD5  =      "     CIW      "
!              PCCD5   =      "     CC       "
!              PRFL5   = Diagnostic stratiform rain flux from callpar
!              PSFL5   =      "     stratiform snow flux      "
!              PCRFL5  =      "     convective rain flux      "
!              PCSFL5  =      "     convective snow flux      "
!              PPFRC5  =      "     large-scale precipitation fraction
!              PTSTAR5 = T*
!              PT05    = T0
!              PXYB5   = Output from GPHPRE
!              PAL     = Albedo AT OBS. PONTS
!              VGEOM   = vertical model geometry profiles

! Modifications.
! --------------
!     A.Geer      17-Apr-2015   First version
!     A.Geer      31-Dec-2015   Finally removed the spurious vertical dimension 0:NFLEVG
!     N.Bormann    1-Feb-2016   Slant-path and add lat/lon/hgt
!     A.Geer      06-Apr-2016   Copy, dot-product for test harness/hop_driver
!     A.Geer      01-Jul-2016   Last of 3 supergom/slant-path/vgeom merge bugs
!     P.Lean      17-Nov-2017   Optimisation for slant-path ; change loop ordering
!     P.Lopez     24-May-2017   Moved accumulated precipitation from model to GOM.
!     A.Geer      10-Apr-2018   Rationalise treatment of -ve Q using QNEGLIM which preserve TL/AD sensitivity
!     F.Duruisseau 24-Sep-2018  Add hydrometeors sum (resol + conv) for BAYRAD
!     B.Ingleby   09-Jan-2019   Add 2m humidity and 10m wind to GOM.
!     Y.Hirahara  01-Apr-2019   Add GOM variables for CMEM (Soil, Snow, Lake)
!     A.Geer      21-May-2019   NetCDF format gom_plus_dump to create colocated model datasets
!     N.Bormann   06-Aug-2020   Modifications required for slant-path for all-sky
!----------------------------------------------------------------------------

use parkind1  , only : jpim, jpib, jprb
use yomhook   , only : lhook, dr_hook, jphook
use yomancs   , only : rmdi
use intdyn_mod, only : yytxyb
use yomsta    , only : nlextrap
use yomct0    , only : lecmwf, l_oops, lelam
!use yoephy    , only : yrephy
!use yomgem    , only : tgem
use yomcst    , only : rg, rpi
use yomcoctp  , only : nscatt
!use yomphy    , only : yrphy
use yomvert   , only : tvertical_geom,alloc_copy_vertical_geom,dealloc_vertical_geom
use pardim    , only : jpnppm
use yomvar    , only : lecv
use yomjbecv  , only : lsktecv

implicit none

private
public type_gom_plus, ih


type type_gom_plus
  integer(kind=jpim) :: ndlen, nflevg, nhoriz, nppm, nxyb, ngems, nlev_soil
  logical            :: lsrg, lhydro_sum, lphys, lcanari
  real(kind=jprb) :: missing_value
  real(kind=jprb), pointer, contiguous :: lat   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: lon   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: timestep(:,:)  => null()
  real(kind=jprb), pointer, contiguous :: orog  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: al    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: tstar (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: t0    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: ts    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: ts_ir (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: ts_mw (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: ls    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: ci    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: ws    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: sn    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: u10n  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: v10n  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: u0    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: v0    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: chn   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: tcls  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: hucls (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: ucls  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: vcls  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: nucls (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: nvcls (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: z0    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: wl    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: sst   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: veg   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: cptgz (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: wet   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: qs    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: cpts  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: ri    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: bn    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: bm    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: bh    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: upd   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: arg   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: hv    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: z0h   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: wsi   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: rrr   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: rsn   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: es    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: prec_accum (:,:) => null()
  real(kind=jprb), pointer, contiguous :: sdfor (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: t2m   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: d2m   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: q2m   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: u10m  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: v10m  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: soty  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: lail  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: laih  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: sndt  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: sndp  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: cvl   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: cvh   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: tvl   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: tvh   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: tlk   (:,:)    => null()
!canari
  real(kind=jprb), pointer, contiguous :: cori  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: eh    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: ez    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: est   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: esn   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: et2   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: eh2   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: ev1   (:,:)    => null()
!canari
  real(kind=jprb), pointer, contiguous :: presh (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: presf (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: geoph (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: geopf (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: aphi  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: cp    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: rf    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: kap   (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: uf    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: vf    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: tf    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: qf    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: qn    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: o3f   (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: clwf  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: ciwf  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: ccf   (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: cswf  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: crwf  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: cgwf  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: clwd  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: ciwd  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: clw_sum  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: ciw_sum  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: crw_sum  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: csw_sum  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: ccd   (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: phys_type(:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: rfl   (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: sfl   (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: crfl  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: csfl  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: pfrc  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: hgt   (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: tsoil (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: qsoil (:,:,:)  => null()
  real(kind=jprb), allocatable :: gems  (:,:,:,:)
  integer(kind=jpim), allocatable :: gems_igrib (:)
  integer(kind=jpim), allocatable :: gems_type (:)
  real(kind=jprb), pointer, contiguous :: pxp   (:,:,:,:) => null()
  real(kind=jprb), pointer, contiguous :: pxpd  (:,:,:,:) => null()
  real(kind=jprb), pointer, contiguous :: pxz   (:,:,:,:) => null()
  real(kind=jprb), pointer, contiguous :: pxzd  (:,:,:,:) => null()
  real(kind=jprb), allocatable :: pxyb  (:,:,:,:)
! Zonal + meridonal component of geographical north direction (angle):
  real(kind=jprb), pointer, contiguous :: gnordl  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: gnordm  (:,:)    => null()
  real(kind=jprb), private, allocatable :: store_sfc(:,:,:)
  real(kind=jprb), private, allocatable :: store_ua0(:,:,:,:)
  real(kind=jprb), private, allocatable :: store_ua1(:,:,:,:)
  real(kind=jprb), private, allocatable :: store_soil(:,:,:,:)
  real(kind=jprb), private, allocatable :: store_ppm(:,:,:,:,:)
  type(tvertical_geom)    , allocatable :: vgeom
end type
save

! In lieu of anything more complicated (a future enhancement...) this is
! how we reduce the 2D-permitting arrays to 1D for the majority of obsops:
integer(kind=jpim), parameter :: ih=1

integer(kind=jpim), parameter :: gems_type_ghg = 1
integer(kind=jpim), parameter :: gems_type_aero = 2
integer(kind=jpim), parameter :: gems_type_chem = 3

#include "abor1.intfb.h"
end module gom_plus
