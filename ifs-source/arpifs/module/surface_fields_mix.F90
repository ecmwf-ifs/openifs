! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE SURFACE_FIELDS_MIX

!     Purpose.
!     --------

!      SURFACE_FIELDS_MIX contains data structures and manipulation routines
!      for the surface (physics) fields in the IFS

!      This module is a mix of declarations, type definitions and
!       subroutines linked with surface fields. There are four parts:
!       1/ Declaration of dimensions (including some parameter variables).
!       2/ Definition of types.
!       3/ Declarations:
!          Declaration of variables SP_[group], YSP_[group]D, YSP_[group]
!           (prognostic surface fields).
!          Declaration of variables SD_[group], YSD_[group]D, YSD_[group]
!           (diagnostic surface fields).
!       4/ Some routines linked to the surface data flow:
!          * INI_SFLP3: Initialize 3-D surface field group
!          * SETUP_SFLP3: Setup 3-D surface field
!          * FINALISE_SFLP3: Finalise 3-D surface field
!          * INI_SFLP2: Initialize 2-D surface field group
!          * SETUP_SFLP2: Setup 2-D surface field
!          * FINALISE_SFLP2: Finalise 2-D surface field
!          * GPPOPER: Operations on prognostic surface fields
!          * GPOPER: Operations on ALL surface groups
!          * GPOPER_2: Operations on 2-D surface groups
!          * GPOPER_3: Operations on 3-D surface groups
!          * SURF_STORE: Store all surface fields
!          * SURF_RESTORE: Restore all surface fields
!          * ALLO_SURF: Allocate surface field arrays

!     Author.
!     -------
!     Mats Hamrud(ECMWF)

!     Modifications.
!     --------------
!        Original : 2006-07-01
!        Modifications:
!        K. Yessad (25 Oct 2006): rephase ALARO0 contribution.
!        K. Yessad (26 Oct 2006): add missing comments.
!        G. Balsamo (14 Mar 2007): add soil type allocatable (SOTY).
!        Jean Bidlot (June 2007):  named pointer for fields to wave model.
!        S. Serrar (17 Jul 2007) methane surface fields pointers
!        Y. Takaya (?? Nov 2008): add ocean mixed layer model fields
!        A. Alias  (07 Aug 2009) move field Nudging mask from group VCLIA to VARSF
!        JJMorcrette 20091201 Total and clear-sky direct SW radiation flux at surface
!        H. Hersbach (04-Dec-2009): 10m-neutral wind and friction velocity,
!                                   introduce YDUPD, SETPERTOVAL in GPOPEN
!        S. Boussetta/G.Balsamo (May 2009): add high/low vegetation LAI pointer (LAIH/LAIL)
!        Y. Takaya/P. de Rosnay May 2010: SSS for SMOS
!        R. Forbes (01-Mar-2010): Added TCRW,TCSW diagnostics
!        JJMorcrette 20100212 PP of CBASE, 0DEGL and VISIH
!        T. Wilhelmsson (25 Mar 2010) add 6 hourly min/max fields
!        A. Alias  (07 Aug 2009) move field Nudging mask from group VCLIA to VARSF
!        Y. Bouteloup (20 Oct 2010) Surface forcing for 1D model : SFORC
!        Y. Bouteloup (04 Jan 2011) Surface flux for EDKF and 1D model surface forcing : SFLUX
!        G.Balsamo/S.Boussetta (Apr 2011): add land carbon dioxide fluxes
!        H. Hersbach (01 April 2011): auxiliary diagnostic radiation fields
!        P. Marguinaud (07 November 2012): Fix uninitialized component
!        P. Bechtold (9 Aug 2011): add CIN, Convective Indices
!        P.Marguinaud (11 Sep 2012): Initialize TYPE_SURF_MTL_2D%CNAME (avoid valgrind warning)
!        M. Ahlgrimm 31 Oct 2011: clear-sky downward radiation at surface
!        L. Jones (25-Oct-2011): Created FINALISE_SFLP2/3, removing need for
!                                pre-counting number of SETUP_SFLP2/3 calls
!        M. Fisher   7-March-2012 Use DEALLOCATE_IF_ASSOCIATED
!        JJMorcrette 20130213 PP optical depths GEMS/MACC aerosols
!        G. Balsamo  14-Jun-2013 Introduce lake buffer group SL=LAKEB
!        R. Forbes 01-March-2014 Add precip rates/type,TCSLW,I10FG,PEV
!        M. Ahlgrimm Apr 2014: Add precip fraction for DDH output
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!        JJMorcrette 20130730 15-variable aerosol model + oceanic DMS
!        R. Forbes 10-Jan-2015 Add freezing rain FZRA
!        S. Rémy, 23/1/2015, injection height for biomass burning aerosol emissions
!        A. Agusti-Panareda (31 Oct 2013): add GPP/REC flux adjustment coefficient pointer (CGPP,CREC)
!        P. Bechtold 10-Nov-2015: add CBASEA, CTOPC, ZTWETB0,1
!        A. Bozzo Jan 2016 Add surface direct beam solar radiation
!        S. Rémy, 09/3/2016 add SO2 dry dep velocity pointer SO2DD
!        L. Gerard (Mar 2016) add VK group for deep convection prognostic fields
!        R. El Khatib 09-Sep-2016 Extend VCLIX
!        S. Rémy, 16/9/2016 add nitrate and ammonium optical depths
!        S. Rémy, 21/4/2017 add altitude of volcanoes
!        S. Rémy, 25/4/2017 add calcite fraction over dust 1 and 2
!        S. Rémy, 13/11/2017 SOA optical depth
!        S. Rémy, 16/11/2017 add dust source function
!        E.Dutra/G.Arduini, Jan 2018 change SNOWG to 3D and add new field for snow liquid water (YW)
!        B. Ingleby, 14-01-2019 add Y2SH - specific humidity at 2m
!        R. Hogan 14/01/2019 MODIS Albedo 2x3-components
!        H Petithomme, Oct-2019 add fields to vclih for water conservation
!          could we please avoid too long lines by not indenting? is there any common rule for this?
!        P. Bechtold 11-Feb-2020 add -10C and more CAPE/CIN convection parameters
!        M. Chrust 6/08/2021 Add a copy consructor (used in STEPX)
!        R. Forbes 01-06-2022 Increased JPMAXSFLDS from 622 to 694 for 6*12 prectype freq/sev fields
!--------------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDIM   , ONLY : TDIM
USE YOMLUN   , ONLY : NULOUT, NULERR
USE YOMCT0   , ONLY : NUNDEFLD
USE YOMDYNA  , ONLY : TDYNA

IMPLICIT NONE
SAVE

#include "abor1.intfb.h"
!     -------------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPMAXSFLDS=694 ! Max number of fields in individual group  !KPP
INTEGER(KIND=JPIM), PARAMETER :: JPMAXSTRAJ=100 ! Dimension of NSTRAJGRIB
INTEGER(KIND=JPIM) :: NOFFTRAJ                  ! Offset in surf trajectory
INTEGER(KIND=JPIM) :: NOFFTRAJ_CST              ! Offset in "constant" surf trajectory
INTEGER(KIND=JPIM) :: NSTRAJGRIB(JPMAXSTRAJ)    ! Used in trajectory setup
INTEGER(KIND=JPIM), PRIVATE :: NPTRSURF         ! Used by routine GPOPER

! General type defintions

TYPE TYPE_SURF_MTL
  INTEGER(KIND=JPIM) :: MP           ! Basic field pointer
  INTEGER(KIND=JPIM) :: MP0          ! Field pointer timelevel  0 (prognostic fields)
  INTEGER(KIND=JPIM) :: MP9          ! Field pointer timelevel -1 (prognostic fields)
  INTEGER(KIND=JPIM) :: MP1          ! Field pointer timelevel +1 (prognostic fields)
  INTEGER(KIND=JPIM) :: ITRAJ        !  0 not in trajectory (default)
                                     !  1 in trajectory
                                     !  2 in "constant" trajectory
  INTEGER(KIND=JPIM) :: IBUGFINDER   ! Debug use only. Can be removed if confident.
  LOGICAL            :: LSET         ! True if structure has been set up
END TYPE TYPE_SURF_MTL
! 2D surface field structure
TYPE,EXTENDS(TYPE_SURF_MTL) :: TYPE_SURF_MTL_2D
  INTEGER(KIND=JPIM) :: IGRBCODE     ! GRIB parameter code (default: -999)
  CHARACTER(LEN=16)  :: CNAME = '                '
                                     ! ARPEGE field name   (default: all spaces)

  REAL(KIND=JPRB)    :: REFVALI      ! Default value       (default: 0.0)
  INTEGER(KIND=JPIM) :: NREQIN       ! -1 - initial value from default (default)
                                     ! +1 - initial value from reading file
                                     !  0 - no initial value
END TYPE TYPE_SURF_MTL_2D

! 3D surface field structure
TYPE ,EXTENDS(TYPE_SURF_MTL) :: TYPE_SURF_MTL_3D
  INTEGER(KIND=JPIM), ALLOCATABLE :: IGRBCODE(:)  ! GRIB parameter code (default: -999)
  CHARACTER(LEN=16) , ALLOCATABLE :: CNAME(:)     ! ARPEGE field name   (default: all spaces)
  REAL(KIND=JPRB)   , ALLOCATABLE :: REFVALI(:)   ! Default value       (default: 0.0)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NREQIN(:)    ! -1 - initial value from default (default)
                                                  ! +1 - initial value from reading file
                                                  !  0 - no initial value
END TYPE TYPE_SURF_MTL_3D

! Descriptor pertaining to group
TYPE TYPE_SURF_GEN
  INTEGER(KIND=JPIM) :: NUMFLDS         ! Number of field in group
  INTEGER(KIND=JPIM) :: NDIM            ! Field dimension
  INTEGER(KIND=JPIM) :: NLEVS           ! Number of levels (for multi level groups)
  INTEGER(KIND=JPIM) :: IPTR            ! Internal use
  INTEGER(KIND=JPIM) :: IPTR5           ! Internal use
  INTEGER(KIND=JPIM) :: NDIM5           ! Dimension of trajectory array
  INTEGER(KIND=JPIM) :: NOFFTRAJ        ! Internal use
  INTEGER(KIND=JPIM) :: NOFFTRAJ_CST    ! Internal use
  CHARACTER(LEN=16)  :: CGRPNAME        ! Name of group (for prints)
  LOGICAL            :: L3D             ! TRUE if multi-level field (3-D)
  LOGICAL            :: LMTL            ! TRUE if prognostic field (multi time level)
  LOGICAL            :: FINALISED       ! TRUE if group finalised & no more fields allowed
END TYPE TYPE_SURF_GEN

! Type descriptor for derived type for communicating with GPOPER (see below)
TYPE TYPE_SFL_COMM
  INTEGER(KIND=JPIM) :: IGRBCODE
  LOGICAL            :: L_OK
  CHARACTER(LEN=16)  :: CNAME
  INTEGER(KIND=JPIM) :: IFLDNUM
  REAL(KIND=JPRB)    :: VALUE
  INTEGER(KIND=JPIM) :: IPTRSURF
  INTEGER(KIND=JPIM) :: ICODES(JPMAXSFLDS)
  INTEGER(KIND=JPIM) :: ICODESML(JPMAXSFLDS)
  INTEGER(KIND=JPIM) :: ICOUNT
  INTEGER(KIND=JPIM) :: ICOUNTML
END TYPE TYPE_SFL_COMM

! Group specific type definitions: pronostic groups (SB, SG, RR, CL, OM, EP, X2, CI)

! * Group SB=SOILB: soil prognostic quantities for the different reservoirs
!    (four reservoirs at ECMWF, deep reservoir at METEO-FRANCE):
TYPE TYPE_SFL_SOILB
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YT => NULL()   ! temperature
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YQ => NULL()   ! liquid water content
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YTL=> NULL()   ! ice water content (for MF)
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YSB(:) => NULL()
END TYPE TYPE_SFL_SOILB

! * Group SG=SNOWG: surface snow prognostic quantities:
TYPE TYPE_SFL_SNOWG
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YF =>NULL()   ! content of surface snow
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YA =>NULL()   ! snow albedo
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YR =>NULL()   ! snow density
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YT =>NULL()   ! total albedo (diagnostic for MF for LVGSN)
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YW =>NULL()   ! Liquid water content
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YSG(:) => NULL()
END TYPE TYPE_SFL_SNOWG

! * Group SL=LAKEB: Lake (FLAKE Model) prognostic quantities:
TYPE TYPE_SFL_LAKEB
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YLICT =>NULL() ! lake ice temperature
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YLMLT =>NULL() ! lake mixed-layer temperature
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YLTLT =>NULL() ! lake total layer temperature
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YLBLT =>NULL() ! lake bottom layer temperature
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YLSHF =>NULL() ! lake shape factor
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YLICD =>NULL() ! lake ice depth
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YLMLD =>NULL() ! lake mixed-layer depth
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YSL(:) => NULL()
END TYPE TYPE_SFL_LAKEB

! * Group RR=RESVR: surface prognostic quantities (ECMWF) or
!   surface + superficial reservoir prognostic quantities (MF):
!   Remark:
!    at ECMWF there are 4 soil reservoirs and there is a
!    clear distinction between the soil reservoirs (group SOILB)
!    and the surface (group RESVR);
!    at METEO-FRANCE there is a deep reservoir (group SOILB) and a
!    superficial reservoir (group RESVR):
!    - there is a skin surface temperature (Ts) which is the temperature at the
!      interface surface/superficial reservoir (and not two separate quantities
!      for superficial reservoir and surface)
!    - there is a skin surface water content (denoted by Wl) and a superficial
!      reservoir water content (denoted by Ws).
!    - there is a superficial reservoir ice content but no surface ice content.
!    (remark k.y.: it would have been more logical to use group name
!    RESVR for internal reservoirs and group name SOILB for surface!).
TYPE TYPE_SFL_RESVR
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YT  =>NULL()  ! skin temperature (Ts)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YW  =>NULL()  ! skin water content (Wskin) at ECMWF
                                                   ! superficial reservoir water content (Ws) at MF
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YFC =>NULL()  ! skin water content (Wl) at MF
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YIC =>NULL()  ! superficial reservoir ice
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YFP1=>NULL()  ! interpolated Ts for 2nd part of 927-FULLPOS
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YRR(:) => NULL()
END TYPE TYPE_SFL_RESVR

! * Group CL=CLS: surface boundary layer prognostic quantities:
TYPE TYPE_SFL_CLS
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCLS =>NULL()  ! 2m temperature
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YHUCLS=>NULL()  ! 2m humidity
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YUCLS =>NULL()  ! 10m U-wind
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVCLS =>NULL()  ! 10m V-wind
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YNUCLS=>NULL()  ! 10m neutral U-wind
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YNVCLS=>NULL()  ! 10m neutral V-wind
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCL(:) => NULL()
END TYPE TYPE_SFL_CLS

! * Group OM=OML: prognostic quantities for ocean mixed layer model (KPP/TKE):
TYPE TYPE_SFL_OML
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YTO=>NULL()     ! temperature
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YSO=>NULL()     ! salinity
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YUO=>NULL()     ! U velocity
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YVO=>NULL()     ! V velocity
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YOM(:) => NULL()
END TYPE TYPE_SFL_OML

! * Group EP=EXTRP: extra 3-d prognostic fields:
TYPE TYPE_SFL_EXTRP
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YEP(:) => NULL()
END TYPE TYPE_SFL_EXTRP

! * Group X2=XTRP2: extra 2-d prognostic fields:
!   (is used for precipitation fields in CANARI)
TYPE TYPE_SFL_XTRP2
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YX2(:) => NULL()
END TYPE TYPE_SFL_XTRP2

! * Group CI=CANRI: 2-d prognostic fields for CANARI:
TYPE TYPE_SFL_CANRI
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCI(:) => NULL()
END TYPE TYPE_SFL_CANRI

! Group specific type definitions: diagnostic groups

! * Group VF=VARSF: climatological/geographical diagnostic fields:
TYPE TYPE_SFL_VARSF
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YZ0F  =>NULL()  ! gravity * surface roughness length
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALBF =>NULL()  ! surface shortwave albedo
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YEMISF=>NULL()  ! surface longwave emissivity
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YGETRL=>NULL()  ! standard deviation of orography
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLSM  =>NULL()  ! land-sea mask
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVEG  =>NULL()  ! vegetation cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVRLAN=>NULL()  ! anisotropy of the sub-grid scale orography
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVRLDI=>NULL()  ! angle of the direction of orography with the x axis
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSIG  =>NULL()  ! characteristic orographic slope
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALBSF=>NULL()  ! soil shortwave albedo
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLAN  =>NULL()  ! fraction of land
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSST  =>NULL()  ! (open) sea surface temperature
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSSS  =>NULL()  ! sea surface salinity
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLZ0H =>NULL()  ! logarithm of roughness length for heat
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCVL  =>NULL()  ! low vegetation cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCO2TYP=>NULL() ! co2 photosynthesis type (c3/c4) for low vegetation cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCVH  =>NULL()  ! high vegetation cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCUR  =>NULL()  ! urban cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTVL  =>NULL()  ! low vegetation type
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTVH  =>NULL()  ! high vegetation type
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YFWET =>NULL()  ! wetland fraction
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLAIL =>NULL()  ! low vegetation LAI
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLAIH =>NULL()  ! high vegetation LAI
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSOTY =>NULL()  ! soil type
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCLK  =>NULL()  ! lake cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YDL   =>NULL()  ! lake depth
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCI   =>NULL()  ! sea ice fraction
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YUCUR =>NULL()  ! U-component of the ocean current
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVCUR =>NULL()  ! V-component of the ocean current
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YZ0RLF=>NULL()  ! gravity * vegetation roughness length
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCGPP =>NULL()  ! GPP bias correction factor
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCREC =>NULL()  ! REC bias correction factor
  ! 4-component MODIS land albedo climatology
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSDFOR=>NULL()  ! SD filtered orography
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALUVP=>NULL()  ! MODIS-derived parallel albedo for shortwave radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALUVD=>NULL()  ! MODIS-derived diffuse albedo for shortwave radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALNIP=>NULL()  ! MODIS-derived parallel albedo for longwave radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALNID=>NULL()  ! MODIS-derived diffuse albedo for longwave radiation
  ! 6-component MODIS land albedo climatology
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALUVI=>NULL()  ! MODIS-derived parallel albedo for UV-visible radiation (iso.)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALNII=>NULL()  ! MODIS-derived parallel albedo for near-infrared radiation (iso.)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALUVV=>NULL()  ! MODIS-derived parallel albedo for UV-visible radiation (volu.)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALNIV=>NULL()  ! MODIS-derived parallel albedo for near-infrared radiation (volu.)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALUVG=>NULL()  ! MODIS-derived parallel albedo for UV-visible radiation (geom.)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALNIG=>NULL()  ! MODIS-derived parallel albedo for near-infrared radiation (geom.)
  ! Indices to start/end albedo climatology coefficients currently in use
  INTEGER(KIND=JPIM)              :: IALSTART, IALEND
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YFP1  =>NULL()  ! surface orography in the 2nd part of FULLPOS-927
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSO2DD=>NULL()  ! sulphate dry dep velocity
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YDMSO =>NULL()  ! oceanic DMS
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YURBF =>NULL()  ! Urban fraction
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YFCA1  =>NULL() ! Fraction of calcite over dust 1st bin
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YFCA2  =>NULL() ! Fraction of calcite over dust 2nd bin
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YAERDEP=>NULL() ! dust emission potential
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YAERLTS=>NULL() ! dust lifting threshold speed
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YAERSCC=>NULL() ! dust soil clay content
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YDSF   =>NULL() ! dust source function
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YDSZ   =>NULL() ! dust size distribution modulation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCHEMFLXO(:)=>NULL() ! total chemistry flux (emissions + deposition)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCHEMWDFLX(:)=>NULL() ! chemistry wet deposition flux 
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCHEMDDFLX(:)=>NULL() ! chemistry dry deposition flux 
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCHEMDV(:)=>NULL() ! chemistry deposition velocity
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YEMIS2D(:)=>NULL()! 2D emission fluxes for composition
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YEMIS2DAUX(:)=>NULL()! 2D emission auxiliary fields for composition (e.g. injection heights)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YNUDM=>NULL()   ! nudging mask
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVF(:) => NULL()
END TYPE TYPE_SFL_VARSF

! * Group VP=VCLIP: deep soil diagnostic fields:
TYPE TYPE_SFL_VCLIP
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTPC=>NULL()    ! climatological deep layer temperature
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YWPC=>NULL()    ! climatological deep layer moisture
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVP(:) => NULL()
END TYPE TYPE_SFL_VCLIP

! * Group VV=VCLIV: vegetation diagnostic fields:
TYPE TYPE_SFL_VCLIV
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YARG  =>NULL()  ! silt percentage within soil
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSAB  =>NULL()  ! percentage of sand within the soil
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YD2   =>NULL()  ! soil depth
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YIVEG =>NULL()  ! type of vegetation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YRSMIN=>NULL()  ! stomatal minimum resistance
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLAI  =>NULL()  ! leaf area index
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YHV   =>NULL()  ! resistance to evapotranspiration
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YZ0H  =>NULL()  ! gravity * roughness length for heat
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALS  =>NULL()  ! albedo of bare ground
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALV  =>NULL()  ! albedo of vegetation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVV(:)=> NULL()
END TYPE TYPE_SFL_VCLIV

! * Group VN=VCLIN: cloudiness diagnostic predictors:
TYPE TYPE_SFL_VCLIN
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTOP   =>NULL() ! index of convective cloud top
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YBAS   =>NULL() ! index of convective cloud base
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YACPR  =>NULL() ! averaged convective precipitaion rate
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YACCPR =>NULL() ! accumulated total precipitaion for assimilation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YACCPR5=>NULL() ! accumulated total precipitaion for assimilation (trajectory)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVN(:) => NULL()
END TYPE TYPE_SFL_VCLIN

! * Group VH=VCLIH: convective cloud diagnostic fields:
TYPE TYPE_SFL_VCLIH
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCCH  => NULL() ! total convective cloudiness
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSCCH  => NULL() ! convective cloud summit
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YBCCH  => NULL() ! convective cloud base
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YPBLH  => NULL() ! PBL height
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSPSH  => NULL() ! variable for prognostic convection scheme (ALARO)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YQSH   => NULL() ! surface moisture historic variable (used by TOUCANS)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVH(:) => NULL()
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YPCL=>NULL()
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YPSL=>NULL()
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YPCN=>NULL()
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YPSN=>NULL()
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YEVA=>NULL()
END TYPE TYPE_SFL_VCLIH

! * Group VK=VCLIK: Convective cloud pseudo-historic fields:
TYPE TYPE_SFL_VCLIK
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YUDGRO => NULL() ! ud top position (accsu)
  TYPE(TYPE_SURF_MTL_2D),POINTER :: YVK(:) => NULL()
END TYPE TYPE_SFL_VCLIK

! * Group VA=VCLIA: aerosol diagnostic fields:
TYPE TYPE_SFL_VCLIA
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSEA=>NULL()    ! aerosol: sea
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLAN=>NULL()    ! aerosol: land
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSOO=>NULL()    ! aerosol: soot
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YDES=>NULL()    ! aerosol: desert
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSUL=>NULL()    ! aerosol: sulfate
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVOL=>NULL()    ! aerosol: volcano
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVA(:) => NULL()
END TYPE TYPE_SFL_VCLIA

! * Group VG=VCLIG: ice-coupler diagnostic fields:
!   ky: currently not used, missing setup in su_surf_flds.F90
TYPE TYPE_SFL_VCLIG
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YICFR=>NULL()   ! sea-ice fraction
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSOUP=>NULL()   ! upward solar flux over sea-ice
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YIRUP=>NULL()   ! upward IR flux over sea-ice
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCHSS=>NULL()   ! sensible heat over sea-ice
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YEVAP=>NULL()   ! evaporation over sea-ice
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTAUX=>NULL()   ! U-component of stress over sea-ice
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTAUY=>NULL()   ! V-component of stress over sea-ice
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVG(:) => NULL()
END TYPE TYPE_SFL_VCLIG

! * Group VC=VO3ABC: A,B and C (Climatological ozone profiles) diagnostic fields:
TYPE TYPE_SFL_VO3ABC
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YA=>NULL()      ! A climatological ozone profile
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YB=>NULL()      ! B climatological ozone profile
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YC=>NULL()      ! C climatological ozone profile
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVC(:) => NULL()
END TYPE TYPE_SFL_VO3ABC

! * Group V2=VDIAGO2: 2-D climatological/diagnostic fields for an ocean mixed layer model (KPP):
TYPE TYPE_SFL_VDIAGO2
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YOCDEP=>NULL()  ! bottom layer depth
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YUSTRC=>NULL()  ! taux clim.
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YVSTRC=>NULL()  ! tauy clim.
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YV2(:) => NULL()
END TYPE TYPE_SFL_VDIAGO2

! * Group V3=VDIAGO3: 3-D climatological/diagnostic fields for an ocean mixed layer model (KPP or TKE):
TYPE TYPE_SFL_VDIAGO3
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YDIFM  =>NULL() ! viscosity
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YDIFT  =>NULL() ! diff. coef. of temp
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YDIFS  =>NULL() ! diff. coef. of salinity
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YADVT  =>NULL() ! correction term for temp.
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YADVS  =>NULL() ! correction term for sal.
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YTRI0  =>NULL() ! coef. for solving matrix.
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YTRI1  =>NULL() ! coef. for solving matrix.
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YSWDK  =>NULL() ! radiation term
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YZO    =>NULL() ! depth of layer
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YHO    =>NULL() ! depth of interface layer
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YDO    =>NULL() ! layer thickness
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YHO_INV=>NULL() ! 1 / YHO
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YUOC   =>NULL() ! U velocity clim.
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YVOC   =>NULL() ! V velocity clim.
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YOTKE  =>NULL() ! Ocean Turbulent kinetic energy
  TYPE (TYPE_SURF_MTL_3D), POINTER :: YV3(:) => NULL()
END TYPE TYPE_SFL_VDIAGO3

! * Group VD=VDIAG: (ECMWF) diagnostic fields:
TYPE TYPE_SFL_VDIAG
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLSP      =>NULL()  ! Large scale precipitation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCP       =>NULL()  ! Convective precipitation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSF       =>NULL()  ! Snowfall
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YFZRA     =>NULL()  ! Freezing rain
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YBLD      =>NULL()  ! Boundary layer dissipation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSSHF     =>NULL()  ! Surface sensible heat flux
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSLHF     =>NULL()  ! Surface latent heat flux
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YNEE      =>NULL()  ! Surface net ecosystem exchange of CO2
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YGPP      =>NULL()  ! Surface gross primary production of CO2
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YREC      =>NULL()  ! Surface ecosystem respiration of CO2
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMSL      =>NULL()  ! Mean sea level pressure
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSP       =>NULL()  ! Surface pressure
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCC      =>NULL()  ! Total cloud cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y10U      =>NULL()  ! U-wind at 10 m
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y10V      =>NULL()  ! V-wind at 10 m
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y2T       =>NULL()  ! Temperature at 2 m
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y2D       =>NULL()  ! Dewpoint temperature at 2 m
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y2SH      =>NULL()  ! Specific humidity at 2 m
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSSR      =>NULL()  ! Surface solar radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSTR      =>NULL()  ! Surface thermal radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTSR      =>NULL()  ! Top solar radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTTR      =>NULL()  ! Top thermal radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YEWSS     =>NULL()  ! Instantaneous surface U-wind stress
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YNSSS     =>NULL()  ! Instantaneous surface V-wind stress
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YE        =>NULL()  ! Water evaporation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YPEV      =>NULL()  ! Potential evaporation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCCC      =>NULL()  ! Convective cloud cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLCC      =>NULL()  ! Low cloud cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMCC      =>NULL()  ! Medium cloud cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YHCC      =>NULL()  ! High cloud cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLGWS     =>NULL()  ! Zonal gravity wave stress
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMGWS     =>NULL()  ! Meridian gravity wave stress
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YGWD      =>NULL()  ! Gravity wave dissipation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMX2T     =>NULL()  ! Maximum temperature at 2 m
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMN2T     =>NULL()  ! Minimum temperature at 2 m
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMX2T6(:) =>NULL()  ! Bins for maximum temperature at 2 m since last 6 hours
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMN2T6(:) =>NULL()  ! Bins for minimum temperature at 2 m since last 6 hours
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YRO       =>NULL() ! Runoff (total)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSRO      =>NULL() ! Runoff surface
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSSRO     =>NULL() ! Runoff sub-surface
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YALB      =>NULL() ! (surface shortwave) albedo
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YIEWSS    =>NULL() ! Instantaneous surface zonal component of stress
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YINSSS    =>NULL() ! Instantaneous surface meridian component of stress
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YISSHF    =>NULL() ! Instantaneous surface heat flux
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YIE       =>NULL() ! Instantaneous surface moisture flux
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YINEE     =>NULL() ! Instantaneous net ecosystem exchange of CO2
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YIGPP     =>NULL() ! Instantaneous gross primary production of CO2
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YIREC     =>NULL() ! Instantaneous ecosystem respiration of CO2
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YICH4     =>NULL() ! Instantaneous wetland CH4 flux
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCSF      =>NULL() ! Convective snow fall
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLSSF     =>NULL() ! Large scale snowfall
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMXTPR    =>NULL() ! Max precip rate since last post-processing
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMNTPR    =>NULL() ! Min precip rate since last post-processing
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMXTPR6(:)=>NULL() ! Max precip rate in last 6 hours
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMNTPR6(:)=>NULL() ! Min precip rate in last 6 hours
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTPR      =>NULL() ! Total precipitation rate
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLSRR     =>NULL() ! Large scale rain rate
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCRR      =>NULL() ! Convective rain rate
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLSSFR    =>NULL() ! Large scale snowfall rate
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCSFR     =>NULL() ! Convective snowfall rate
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YPTYPE    =>NULL() ! Precipitation type
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YILSPF    =>NULL() ! Large-scale precipitation fraction (inst.)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YZ0F      =>NULL() ! Gravity * surface roughness length
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLZ0H     =>NULL() ! Logarithm of z0 times heat flux
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVIWVE    =>NULL() ! Vertical integral of eastward water vapour flux
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVIWVN    =>NULL() ! Vertical integral of northward water vapour flux
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCW      =>NULL() ! Total water content in a vertical column
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCWV     =>NULL() ! Total water vapor content in a vertical column
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCLW     =>NULL() ! Total liquid water content in a vertical column
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCIW     =>NULL() ! Total ice water content in a vertical column
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCRW     =>NULL() ! Total rain water content in a vertical column
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCSW     =>NULL() ! Total snow water content in a vertical column
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCSLW    =>NULL() ! Total supercooled liquid water content in a vertical column
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSSRD     =>NULL() ! Downward surface solar radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSTRD     =>NULL() ! Downward surface thermic radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSSRDC    =>NULL() ! Clear-sky downward surface solar radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSTRDC    =>NULL() ! Claer-sky downward surface thermal radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YBLH      =>NULL() ! Height of boundary layer
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSUND     =>NULL() ! Sunshine duration
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSPAR     =>NULL() ! Surface downward PARadiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSUVB     =>NULL() ! Surface downward UV-B radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSFDIR    =>NULL() ! Surface total sky direct downward SW radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSCDIR    =>NULL() ! Surface clear-sky direct downward SW radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSDSRP    =>NULL() ! Surface total-sky direct beam downward SW radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCAPE     =>NULL() ! Conv.avail.potential energy (CAPE)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMUCAPE   =>NULL() ! most unstable CAPE
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMLCAPE50 =>NULL() ! CAPE from 50 hPa mixed-layer
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMLCAPE100=>NULL() ! CAPE from 100 hPa mixed-layer
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YPDEPL    =>NULL() ! Pressure at parcel departure level for MUCAPE
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCAPES    =>NULL() ! CAPE-Shear
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMXCAP6(:)=>NULL() ! Bins for maximum CAPE in last 6 hours
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMXCAPS6(:)=>NULL()! Bins for maximum CAPE-Shear in last 6 hours
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTROPOTP  =>NULL() ! Pressure of thermal tropopause
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTSRC     =>NULL() ! Top solar radiation clear sky
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTTRC     =>NULL() ! Top thermal radiation clear sky
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSSRC     =>NULL() ! Surface solar radiation clear sky
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSTRC     =>NULL() ! Surface thermal radiation clear sky
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YES       =>NULL() ! Evaporation of snow
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSMLT     =>NULL() ! Snow melt
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y10FG     =>NULL() ! Wind gust at 10 m (max since previous pp)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y10FG6(:) =>NULL() ! Bins for wind gust at 10 m (max since last 6 hours)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y10FGCV   =>NULL() ! convective wind gust at 10m for current time level (m/s)
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YI10FG    =>NULL() ! Wind gust at 10 m ("instantaneous")
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLSPF     =>NULL() ! Large scale precipitation fraction
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCO3     =>NULL() ! Total ozone content in a vertical column
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVIMD     =>NULL() ! Vertically integrated mass divergence
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSPARC    =>NULL() ! Surface clear-sky parallel radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSTINC    =>NULL() ! Top of atmosphere incident solar radiation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCBASE    =>NULL() ! Cloud base level
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y0DEGL    =>NULL() ! Zero deg. level
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YM10DEGL  =>NULL() ! -10 deg. level
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVISIH    =>NULL() ! Horizontal visibility
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCIN      =>NULL() ! CIN
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMLCIN50  =>NULL() ! CIN from 50 hpa mixed-layer
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YMLCIN100 =>NULL() ! CIN from 100 hpa mixed-layer
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YKINDEX   =>NULL() ! Convective K-Index
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTTINDEX  =>NULL() ! Convective TT-Index
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCBASEA   =>NULL() ! Cloud base aviation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCTOPC    =>NULL() ! Cloud top convective
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YZTWETB0  =>NULL() ! Height of 0 deg wet bulb temperature
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YZTWETB1  =>NULL() ! Height of 1 deg wet bulb temperature
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCGHG(:) =>NULL() ! Total column greenhouse gases
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTCCHEM(:)=>NULL() ! Total column chemistry
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YAERODIAG(:,:)=>NULL()! Per-aerosol-type diagnostics
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YAERO_WVL_DIAG(:,:)=>NULL()! Per-wavelength aerosol optical diagnostics
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y100U     =>NULL() ! 100m zonal wind
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y100V     =>NULL() ! 100m meridional wind
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y200U     =>NULL() ! 200m zonal wind
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y200V     =>NULL() ! 200m meridional wind
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YZUST     =>NULL() ! Friction velocity
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y10NU     =>NULL() ! 10m zonal neutral wind
  TYPE(TYPE_SURF_MTL_2D), POINTER :: Y10NV     =>NULL() ! 10m meridional neutral wind
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YDNDZN    =>NULL() ! Minimum vertical refractivity gradient
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YDNDZA    =>NULL() ! Mean vertical refractivity gradient
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YDCTB     =>NULL() ! Duct base height
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTPLB     =>NULL() ! Trapping layer base height
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTPLT     =>NULL() ! Trapping layer top height
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODSS     =>NULL() ! optical depth sea salt aerosols
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODDU     =>NULL() ! optical depth dust aerosols
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODOM     =>NULL() ! optical depth organic m. aerosols
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODBC     =>NULL() ! optical depth black C aerosols
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODSU     =>NULL() ! optical depth sulphate aerosols
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODNI     =>NULL() ! optical depth nitrate aerosols
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODAM     =>NULL() ! optical depth ammonium aerosols
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODSOA    =>NULL() ! optical depth secondary organic aerosols
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODVFA    =>NULL() ! optical depth volcanic flying ash
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODVSU    =>NULL() ! optical depth volcanic sulphate aerosols
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YODTOACC  =>NULL() ! optical depth total aerosol accumulated
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YAEPM1    =>NULL() ! particulate matter le 1 um
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YAEPM25   =>NULL() ! particulate matter le 2.5um
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YAEPM10   =>NULL() ! particulate matter le 10 um
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YUVBED    =>NULL() ! UV biologically effective dose
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YUVBEDCS  =>NULL() ! UV biologically effective dose clear sky
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLITOTI   =>NULL() ! instantaneous total lightning flash density
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLICGI    =>NULL() ! instantaneous cloud-to-ground lightning flash density
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLITOTA6(:)=>NULL()! Bins for averaged total lightning over last 6 hours
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLICGA6(:)=>NULL() ! Bins for averaged cloud-to-ground lightning over last 6 hours
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSDSL     =>NULL() ! Snow depth total (aggr from ml) in m of water eq
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVD(:) => NULL()
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YPTYPEOCC6(:,:)=>NULL() ! Bins for precip type accum occurrence frequency over last 6 hours
END TYPE TYPE_SFL_VDIAG

! * Group SM=SATSIM: (ECMWF) simulated satellite images:
TYPE TYPE_SFL_SATSIM
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YCLBT=>NULL()      ! Cloudy brightness temperature
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YCSBT=>NULL()      ! Clear-sky brightness temperature
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YSM(:) => NULL()
END TYPE TYPE_SFL_SATSIM

! * Group WS=WAVES: surface prognostic quantities over sea (used by IFS):
TYPE TYPE_SFL_WAVES
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YWS(:) => NULL()
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCHAR   =>NULL()   ! Charnock parameter as modified by the wave model.
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCHARHQ =>NULL()   ! Charnock parameter for heat and moisture as modified by the wave model.
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YUSTOKES=>NULL()   ! U-component of the surface Stokes drift.
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVSTOKES=>NULL()   ! V-component of the surface Stokes drift.
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTAUOCX =>NULL()   ! U-component of the Momentum flux to ocean.
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTAUOCY =>NULL()   ! V-component of the Momentum flux to ocean.
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YPHIOC  =>NULL()   ! Energy flux to ocean.
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YWSEMEAN=>NULL()   ! Windsea variance.
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YWSFMEAN=>NULL()   ! Windsea mean frequency.
END TYPE TYPE_SFL_WAVES

! * Group WW=WAM: surface prognostic quantities over sea (used by WAM):
TYPE TYPE_SFL_WAM
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YWW(:) => NULL()
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YU10N =>NULL()     ! 10m neutral wind U-component passed to the wave model (WAM).
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YV10N =>NULL()     ! 10m neutral wind V-component passed to the wave model (WAM).
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YRHO  =>NULL()     ! surface density passed to the wave model (WAM).
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YZIL  =>NULL()     ! ZI/L passed to the wave model (used for gustiness in WAM).
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCIF  =>NULL()     ! Sea ice fraction passed to the wave model (WAM).
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YCLK  =>NULL()     ! Lake cover passed to the wave model (WAM).
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YUCURW=>NULL()     ! Ocean current    U-component passed to the wave model (WAM).
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVCURW=>NULL()     ! Ocean current    V-component passed to the wave model (WAM).
END TYPE TYPE_SFL_WAM

! * Group VX=VCLIX: auxilary climatological diagnostic fields:
TYPE TYPE_SFL_VCLIX
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YORO  =>NULL()  ! climatological surface geopotential
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTSC  =>NULL()  ! climatological surface temperature
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YPWS  =>NULL()  ! climatological surface max. prop. moisture
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YPWP  =>NULL()  ! climatological deep soil max. prop. moisture
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSNO  =>NULL()  ! climatological snow cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YTPC  =>NULL()  ! climatological deep soil temperature
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSAB  =>NULL()  ! climatologic percentage of sand within the soil
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YXD2  =>NULL()  ! climatologic soil depth
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLSM  =>NULL()  ! climatologic land sea mask
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YIVEG =>NULL()  ! climatologic type of vegetation
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YARG  =>NULL()  ! silt percentage within soil
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YRSMIN=>NULL()  ! climatologic stomatal minimum resistance
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YLAI  =>NULL()  ! leaf area index
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVEG  =>NULL()  ! vegetation cover
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YVX(:) => NULL()
END TYPE TYPE_SFL_VCLIX

! * Group XA=VEXTRA: extra 3-d diagnostic fields:
TYPE TYPE_SFL_VEXTRA
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YXA(:) => NULL()
END TYPE TYPE_SFL_VEXTRA

! * Group XA=VEXTRDI: targeted 3-d diagnostic fields:
TYPE TYPE_SFL_VEXTRDI
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YXA(:) => NULL()
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YXEDR=>NULL() ! Eddy diffusivity rate
END TYPE TYPE_SFL_VEXTRDI

! * Group XA=VPRECIP:  3-d diagnostic fields:
TYPE TYPE_SFL_VPRECIP
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YPRECIP ! Diagnostic of precipitations type
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YXP1(:) => NULL()
END TYPE TYPE_SFL_VPRECIP

TYPE TYPE_SFL_VPRECIP2
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YPRECIP2 ! Diagnostic of precipitations type
  TYPE(TYPE_SURF_MTL_3D), POINTER :: YXP2(:) => NULL()
END TYPE TYPE_SFL_VPRECIP2

! * Group X2=VEXTR2: extra 2-d diagnostic fields:
TYPE TYPE_SFL_VEXTR2
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YX2(:) => NULL()
END TYPE TYPE_SFL_VEXTR2

! * Group SFL:SFLUX Surface flux for EDKF
TYPE TYPE_SFL_SFLUX
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSFL(:) => NULL()
END TYPE TYPE_SFL_SFLUX

! * Group SFO:SFORC Surface forcing for 1D model (MUSC)
TYPE TYPE_SFL_SFORC
  TYPE(TYPE_SURF_MTL_2D), POINTER :: YSFO(:) => NULL()
END TYPE TYPE_SFL_SFORC

! * Group OC=OCE: ocean model fields :
TYPE TYPE_SFL_OCE
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YSSTM =>NULL()  ! sst masked
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YCIM  =>NULL()  ! ci masked
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YUCURM=>NULL()  ! ucur masked
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YVCURM=>NULL()  ! vcur masked
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YICTH =>NULL()  ! ice thickness
  TYPE (TYPE_SURF_MTL_2D), POINTER :: Y20D  =>NULL()  ! 20 degree isotherm
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YSSH  =>NULL()  ! sea level
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YMLD  =>NULL()  ! mixed layer depth
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YSSS  =>NULL()  ! sea surface salinity
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YTEM3 =>NULL()  ! pottemp avg 300m
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YSAL3 =>NULL()  ! salt avg 300m
  TYPE (TYPE_SURF_MTL_2D), POINTER :: YOC(:) => NULL()
END TYPE TYPE_SFL_OCE

! End of type definitions

TYPE :: TSURF

INTEGER(KIND=JPIM) :: NSURF=0                   ! Number of surf var.
INTEGER(KIND=JPIM) :: NSURFL=0                  ! Number of surf flds (fields*levels)
INTEGER(KIND=JPIM) :: NDIMSURF=0                ! Total of surf var (includes timelevels etc)
INTEGER(KIND=JPIM) :: NDIMSURFL=0               ! Total dimension of all surface variables
INTEGER(KIND=JPIM) :: NPROGSURF=0               ! Number of prognostic surf var.
INTEGER(KIND=JPIM) :: NPROGSURFL=0              ! Number of prognostic surf flds (fields*levels)

REAL(KIND=JPRB), ALLOCATABLE :: STORE_ARRAY(:,:,:) ! Backup array for surf (see routine SURF_STORE )

! Information on status of tl perturbations
TYPE(TYPE_SFL_COMM) :: YDUPD

! Data structures

! Prognostic (multi time level) fields (SB, SG, SL, RR, CL, OM, EP, X2, CI)

! Soilb
REAL(KIND=JPRB), ALLOCATABLE :: SP_SB (:,:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSP_SBD
TYPE(TYPE_SFL_SOILB)         :: YSP_SB

! Snowg
REAL(KIND=JPRB), ALLOCATABLE :: SP_SG (:,:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSP_SGD
TYPE(TYPE_SFL_SNOWG)         :: YSP_SG

! Lakeb
REAL(KIND=JPRB),ALLOCATABLE :: SP_SL (:,:,:)
TYPE(TYPE_SURF_GEN)     :: YSP_SLD
TYPE(TYPE_SFL_LAKEB)    :: YSP_SL

! Resvr
REAL(KIND=JPRB), ALLOCATABLE :: SP_RR (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSP_RRD
TYPE(TYPE_SFL_RESVR)         :: YSP_RR

! Cls
REAL(KIND=JPRB), ALLOCATABLE :: SP_CL (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSP_CLD
TYPE(TYPE_SFL_CLS)           :: YSP_CL

! Oml (used by ocean mixed layer model (KPP))
REAL(KIND=JPRB), ALLOCATABLE :: SP_OM (:,:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSP_OMD
TYPE(TYPE_SFL_OML)           :: YSP_OM

! Extrp
REAL(KIND=JPRB), ALLOCATABLE :: SP_EP (:,:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSP_EPD
TYPE(TYPE_SFL_EXTRP)         :: YSP_EP

! Xtrp2
REAL(KIND=JPRB), ALLOCATABLE :: SP_X2 (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSP_X2D
TYPE(TYPE_SFL_XTRP2)         :: YSP_X2

! Canri
REAL(KIND=JPRB), ALLOCATABLE :: SP_CI (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSP_CID
TYPE(TYPE_SFL_CANRI)         :: YSP_CI

! Diagnostic 3D ( 2D spatial and 1D time)

!Vprecip
REAL(KIND=JPRB), ALLOCATABLE :: SD_XP(:,:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_XPD
TYPE(TYPE_SFL_VPRECIP)       :: YSD_XP

REAL(KIND=JPRB), ALLOCATABLE :: SD_XP2(:,:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_XP2D
TYPE(TYPE_SFL_VPRECIP2)      :: YSD_XP2

! Diagnostic (one time level) fields

! Varsf
REAL(KIND=JPRB), ALLOCATABLE :: SD_VF (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VFD
TYPE(TYPE_SFL_VARSF)         :: YSD_VF

! Vclip
REAL(KIND=JPRB), ALLOCATABLE :: SD_VP (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VPD
TYPE(TYPE_SFL_VCLIP)         :: YSD_VP

! Vcliv
REAL(KIND=JPRB), ALLOCATABLE :: SD_VV (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VVD
TYPE(TYPE_SFL_VCLIV)         :: YSD_VV

! Vclin
REAL(KIND=JPRB), ALLOCATABLE :: SD_VN (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VND
TYPE(TYPE_SFL_VCLIN)         :: YSD_VN

! Vclih
REAL(KIND=JPRB), ALLOCATABLE :: SD_VH (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VHD
TYPE(TYPE_SFL_VCLIH)         :: YSD_VH

! Vclik
REAL(KIND=JPRB), ALLOCATABLE :: SD_VK (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VKD
TYPE(TYPE_SFL_VCLIK)         :: YSD_VK

! Vclia
REAL(KIND=JPRB), ALLOCATABLE :: SD_VA (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VAD
TYPE(TYPE_SFL_VCLIA)         :: YSD_VA

! Vclig
! currently nothing declared

! Vo3abc
REAL(KIND=JPRB), ALLOCATABLE :: SD_VC (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VCD
TYPE(TYPE_SFL_VO3ABC)        :: YSD_VC

! Vdiago2 (used by ocean mixed layer model (KPP))
REAL(KIND=JPRB), ALLOCATABLE :: SD_V2 (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_V2D
TYPE(TYPE_SFL_VDIAGO2)       :: YSD_V2

! Vdiago3 (used by ocean mixed layer model (KPP))
REAL(KIND=JPRB), ALLOCATABLE :: SD_V3 (:,:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_V3D
TYPE(TYPE_SFL_VDIAGO3)       :: YSD_V3

! Vdiag
REAL(KIND=JPRB), ALLOCATABLE :: SD_VD (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VDD
TYPE(TYPE_SFL_VDIAG)         :: YSD_VD

! Satsim
REAL(KIND=JPRB), ALLOCATABLE :: SD_SM (:,:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_SMD
TYPE(TYPE_SFL_SATSIM)        :: YSD_SM

! Waves (used by IFS)
REAL(KIND=JPRB), ALLOCATABLE :: SD_WS (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_WSD
TYPE(TYPE_SFL_WAVES)         :: YSD_WS

! Waves (used by WAM)
REAL(KIND=JPRB), ALLOCATABLE :: SD_WW (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_WWD
TYPE(TYPE_SFL_WAM)           :: YSD_WW

! Vclix
REAL(KIND=JPRB), ALLOCATABLE :: SD_VX (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_VXD
TYPE(TYPE_SFL_VCLIX)         :: YSD_VX

! Vextra
REAL(KIND=JPRB), ALLOCATABLE :: SD_XA (:,:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_XAD
TYPE(TYPE_SFL_VEXTRA)        :: YSD_XA

! Phys diagnostics
REAL(KIND=JPRB),ALLOCATABLE :: SD_DI (:,:,:,:)
TYPE(TYPE_SURF_GEN)     :: YSD_DID
TYPE(TYPE_SFL_VEXTRDI)   :: YSD_DI

! Vextra-radiation
REAL(KIND=JPRB), ALLOCATABLE :: SD_XR (:,:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_XRD
TYPE(TYPE_SFL_VEXTRA)        :: YSD_XR

! Vextr2
REAL(KIND=JPRB), ALLOCATABLE :: SD_X2 (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_X2D
TYPE(TYPE_SFL_VEXTR2)        :: YSD_X2

! SFLUX
REAL(KIND=JPRB), ALLOCATABLE :: SD_SFL (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_SFLD
TYPE(TYPE_SFL_SFLUX)         :: YSD_SFL

! SFORC
REAL(KIND=JPRB), ALLOCATABLE :: SD_SFO (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_SFOD
TYPE(TYPE_SFL_SFORC)         :: YSD_SFO

! Precip fraction
REAL(KIND=JPRB),ALLOCATABLE :: SD_PF (:,:,:,:)
TYPE(TYPE_SURF_GEN)     :: YSD_PFD
TYPE(TYPE_SFL_VEXTRA)   :: YSD_PF

! Ocean fields from the ocean model
REAL(KIND=JPRB), ALLOCATABLE :: SD_OC (:,:,:)
TYPE(TYPE_SURF_GEN)          :: YSD_OCD
TYPE(TYPE_SFL_OCE)           :: YSD_OC

! Number of fields to be put in TRAJEC%SRFC
INTEGER(KIND=JPIM) :: NGP5_oops
! Number of fields to be put in TRAJEC%CST
INTEGER(KIND=JPIM) :: NTRAJ_CST_oops

END TYPE TSURF

!-------------------------------------------------------------------------


CONTAINS

!=========================================================================

SUBROUTINE INI_SFLP3(YDSC,YD,KLEVS,LDMTL,CDGRPNAME)
! Initialize 3-D surface field group
TYPE(TYPE_SURF_GEN),INTENT(INOUT)    :: YDSC
TYPE(TYPE_SURF_MTL_3D),INTENT(INOUT) :: YD(:)
INTEGER(KIND=JPIM),INTENT(IN)        :: KLEVS
LOGICAL,INTENT(IN)                   :: LDMTL
CHARACTER(LEN=*),INTENT(IN)          :: CDGRPNAME

INTEGER(KIND=JPIM) :: JFLD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:INI_SFLP3',0,ZHOOK_HANDLE)

YDSC%FINALISED = .FALSE.
YDSC%NLEVS = KLEVS
YDSC%IPTR  = 1
YDSC%LMTL  = LDMTL
YDSC%CGRPNAME = CDGRPNAME
YDSC%NDIM5 = 0
YDSC%NOFFTRAJ = NOFFTRAJ
YDSC%NOFFTRAJ_CST = NOFFTRAJ_CST

DO JFLD=1,SIZE(YD)
  YD(JFLD)%MP  = NUNDEFLD
  YD(JFLD)%MP0 = NUNDEFLD
  YD(JFLD)%MP9 = NUNDEFLD
  YD(JFLD)%MP1 = NUNDEFLD
  YD(JFLD)%ITRAJ = 0
  YD(JFLD)%LSET = .FALSE.
  YD(JFLD)%IBUGFINDER = JFLD  ! For debug use only
ENDDO

WRITE(NULOUT,*) 'INITIALIZING 3-D SURFACE FIELD GROUP ', YDSC%CGRPNAME,YDSC%NOFFTRAJ

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:INI_SFLP3',1,ZHOOK_HANDLE)
END SUBROUTINE INI_SFLP3

!=========================================================================

SUBROUTINE SETUP_SFLP3(YDSC,YD,KGRIB,CDNAME,PDEFAULT,KTRAJ,KREQIN)
! Setup 3-D surface field
TYPE(TYPE_SURF_GEN),INTENT(INOUT)      :: YDSC
TYPE(TYPE_SURF_MTL_3D),INTENT(INOUT)   :: YD
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KGRIB(:)
CHARACTER(LEN=16) ,OPTIONAL,INTENT(IN) :: CDNAME(:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN) :: PDEFAULT(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KTRAJ
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KREQIN(:)

INTEGER(KIND=JPIM) :: IPTR,JLEV,J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SETUP_SFLP3',0,ZHOOK_HANDLE)

IPTR = YDSC%IPTR
IF(IPTR > JPMAXSFLDS-1) THEN
  ! We are about to set the last element which is generally used as a
  ! "missing field" placeholder.
  WRITE(NULERR,*) 'SURFACE FIELDS UNDER-DIMENSIONED - GROUP ',&
   & YDSC%CGRPNAME,KGRIB(1),CDNAME(1)
  CALL ABOR1('IPTR > JPMAXSFLDS-1')
ENDIF

IF (YDSC%FINALISED) THEN
  WRITE(NULERR,*) 'ERROR - SETUP_SFLP3 CALLED FOR A GROUP WHICH HAS BEEN FINALISED',&
 &  YDSC%CGRPNAME,KGRIB(1),CDNAME(1)
  CALL ABOR1('SETUP_SFLP3 CALLED FOR A GROUP WHICH HAS BEEN FINALISED')
ENDIF

! If the YD(:) elements are not being set in order it's a possible sign of a
! bug - the caller may have passed in a pointer to the wrong element. Or it
! might be entirely intentional and an abort is unnecessary.
! This if-block and all other references to IBUGFINDER can be commented if
! they cause a problem.
IF (YD%IBUGFINDER /= IPTR) THEN
  WRITE(NULERR,*) 'SETUP_SFLP3 NOT BEING CALLED FOR YD ELEMENTS IN ORDER. '//&
    & '-POSSIBLY- A BUG.',YDSC%CGRPNAME,KGRIB(1),CDNAME(1)
  CALL ABOR1('SETUP_SFLP3 NOT CALLED FOR YD ELEMENTS IN ORDER. POSSIBLE BUG.')
ENDIF

ALLOCATE(YD%IGRBCODE(YDSC%NLEVS))
ALLOCATE(YD%CNAME(YDSC%NLEVS))
ALLOCATE(YD%REFVALI(YDSC%NLEVS))
ALLOCATE(YD%NREQIN(YDSC%NLEVS))

IF(PRESENT(KGRIB)) THEN
  DO J=1,SIZE(YD%IGRBCODE)
    YD%IGRBCODE(J) = KGRIB(J)
  ENDDO
ELSE
  YD%IGRBCODE(:) = -999
ENDIF

IF(PRESENT(KREQIN)) THEN
  DO J=1,SIZE(YD%NREQIN)
    YD%NREQIN(J) = KREQIN(J)
  ENDDO
ELSE
  YD%NREQIN(:) = -1
ENDIF

IF(PRESENT(CDNAME)) THEN
  DO J=1,SIZE(YD%CNAME)
    YD%CNAME(J) = CDNAME(J)
  ENDDO
ELSE
  YD%CNAME(:) = ''
ENDIF

IF(PRESENT(PDEFAULT)) THEN
  DO J=1,SIZE(YD%REFVALI)
    YD%REFVALI(J) = PDEFAULT(J)
  ENDDO
ELSE
  YD%REFVALI(:) = 0.0_JPRB
ENDIF

IF(PRESENT(KTRAJ)) THEN
  IF(KTRAJ == 1) THEN
    IF (NOFFTRAJ+YDSC%NLEVS > JPMAXSTRAJ) CALL ABOR1('JPMAXSTRAJ MUST BE INCREASED!')
    DO JLEV=1,YDSC%NLEVS
      NSTRAJGRIB(NOFFTRAJ+JLEV) = YD%IGRBCODE(JLEV)
    ENDDO
    NOFFTRAJ = NOFFTRAJ+YDSC%NLEVS
  ELSEIF(KTRAJ == 2) THEN
    NOFFTRAJ_CST = NOFFTRAJ_CST+YDSC%NLEVS
  ELSEIF(KTRAJ /= 0) THEN
    CALL ABOR1('SURFACE_FIELDS_MIX:SETUP_SFLP3 - UNKNOWN KTRAJ')
  ENDIF
  ! MC to check with Mats
  !IF(KTRAJ /= 0 .AND. ANY(YD%NREQIN(:) == 0)) THEN
  !  WRITE(NULERR,*) 'ERROR SURFACE_FIELDS_MIX:SETUP_SFLP3 ',YD%IGRBCODE
  !  CALL ABOR1('SURFACE_FIELDS_MIX:SETUP_SFLP3 - FOR TRAJECTORY FIELDS MUST BE INITIALIZED')
  !ENDIF
  YD%ITRAJ = KTRAJ
  YDSC%NDIM5 = YDSC%NDIM5+1
ELSE
  YD%ITRAJ = 0
ENDIF

YD%LSET=.TRUE.
YDSC%IPTR = YDSC%IPTR+1

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SETUP_SFLP3',1,ZHOOK_HANDLE)
END SUBROUTINE SETUP_SFLP3

!=========================================================================

SUBROUTINE FINALISE_SFLP3(YDDYNA,YDSURF,YDSC,YD,LDVERBOSE)
! Finalise 3-D surface field group
TYPE(TDYNA), INTENT(IN)              :: YDDYNA
TYPE(TSURF), INTENT(INOUT)           :: YDSURF
TYPE(TYPE_SURF_GEN),INTENT(INOUT)    :: YDSC
TYPE(TYPE_SURF_MTL_3D),INTENT(INOUT) :: YD(:)
LOGICAL,OPTIONAL,INTENT(IN)          :: LDVERBOSE
INTEGER(KIND=JPIM) :: JFLD, JLEV
LOGICAL :: LLVERBOSE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:FINALISE_SFLP3',0,ZHOOK_HANDLE)

LLVERBOSE =  .TRUE.
IF(PRESENT(LDVERBOSE)) LLVERBOSE=LDVERBOSE
YDSC%FINALISED=.TRUE.
YDSC%NUMFLDS = YDSC%IPTR-1

WRITE(NULOUT,*) '3-D SURFACE FIELD GROUP INITIALIZED ', YDSC%CGRPNAME
WRITE(NULOUT,*) 'NUMFLDS=',YDSC%NUMFLDS,' NLEVS=',YDSC%NLEVS,' LMTL=',YDSC%LMTL

YDSURF%NSURF = YDSURF%NSURF+YDSC%NUMFLDS
YDSURF%NSURFL = YDSURF%NSURFL+YDSC%NUMFLDS*YDSC%NLEVS
IF(YDSC%LMTL) THEN
  YDSURF%NPROGSURF = YDSURF%NPROGSURF+YDSC%NUMFLDS
  YDSURF%NPROGSURFL = YDSURF%NPROGSURFL+YDSC%NUMFLDS*YDSC%NLEVS
ENDIF

IF(YDSC%LMTL) THEN
  IF (YDDYNA%LTWOTL) THEN
    YDSC%NDIM = 2*YDSC%NUMFLDS
  ELSE
    YDSC%NDIM = 3*YDSC%NUMFLDS
  ENDIF
ELSE
  YDSC%NDIM = YDSC%NUMFLDS
ENDIF
YDSURF%NDIMSURF = YDSURF%NDIMSURF + YDSC%NDIM
YDSURF%NDIMSURFL = YDSURF%NDIMSURFL + YDSC%NDIM*YDSC%NLEVS

DO JFLD=1,YDSC%NUMFLDS

  YD(JFLD)%MP = JFLD
  IF (YDSC%LMTL) THEN
    YD(JFLD)%MP0 = YD(JFLD)%MP
    IF(YDDYNA%LTWOTL) THEN
      YD(JFLD)%MP9 = YD(JFLD)%MP0
      YD(JFLD)%MP1 = YD(JFLD)%MP0+YDSC%NUMFLDS
    ELSE
      YD(JFLD)%MP9 = YD(JFLD)%MP0+YDSC%NUMFLDS
      YD(JFLD)%MP1 = YD(JFLD)%MP0+2*YDSC%NUMFLDS
    ENDIF
  ELSE
    YD(JFLD)%MP0 = NUNDEFLD
    YD(JFLD)%MP9 = NUNDEFLD
    YD(JFLD)%MP1 = NUNDEFLD
  ENDIF

  IF(LLVERBOSE) THEN
    DO JLEV=1,YDSC%NLEVS
      IF(YDSC%LMTL) THEN
        WRITE(NULOUT,'(1X,A,2I4,1X,A,6I6)')&
         & YDSC%CGRPNAME(1:6),JFLD,JLEV,YD(JFLD)%CNAME(JLEV),YD(JFLD)%IGRBCODE(JLEV),&
         & YD(JFLD)%MP0,YD(JFLD)%MP9,YD(JFLD)%MP1,YD(JFLD)%ITRAJ,YD(JFLD)%NREQIN(JLEV)
      ELSE
        WRITE(NULOUT,'(1X,A,2I4,1X,A,4I6)')&
         & YDSC%CGRPNAME(1:6),JFLD,JLEV,YD(JFLD)%CNAME(JLEV),YD(JFLD)%IGRBCODE(JLEV),&
         & YD(JFLD)%MP,YD(JFLD)%ITRAJ,YD(JFLD)%NREQIN(JLEV)
      ENDIF
    ENDDO
  ENDIF

ENDDO

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:FINALISE_SFLP3',1,ZHOOK_HANDLE)
END SUBROUTINE FINALISE_SFLP3

!=========================================================================

SUBROUTINE INI_SFLP2(YDSC,YD,LDMTL,CDGRPNAME)
! Initialize 2-D surface field group
TYPE(TYPE_SURF_GEN),INTENT(INOUT)    :: YDSC
TYPE(TYPE_SURF_MTL_2D),INTENT(INOUT) :: YD(:)
LOGICAL,INTENT(IN)                   :: LDMTL
CHARACTER(LEN=*),INTENT(IN)          :: CDGRPNAME

INTEGER(KIND=JPIM) :: JFLD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:INI_SFLP2',0,ZHOOK_HANDLE)

YDSC%FINALISED = .FALSE.
YDSC%NLEVS = -1
YDSC%IPTR  = 1
YDSC%LMTL  = LDMTL
YDSC%CGRPNAME = CDGRPNAME
YDSC%NDIM5 = 0
YDSC%NOFFTRAJ = NOFFTRAJ
YDSC%NOFFTRAJ_CST = NOFFTRAJ_CST

DO JFLD=1,SIZE(YD)
  YD(JFLD)%MP  = NUNDEFLD
  YD(JFLD)%MP0 = NUNDEFLD
  YD(JFLD)%MP9 = NUNDEFLD
  YD(JFLD)%MP1 = NUNDEFLD
  YD(JFLD)%ITRAJ = 0
  YD(JFLD)%LSET = .FALSE.
  YD(JFLD)%IBUGFINDER = JFLD  ! For debug use only.
ENDDO

WRITE(NULOUT,*) 'INITIALIZING 2-D SURFACE FIELD GROUP ', YDSC%CGRPNAME,YDSC%NOFFTRAJ

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:INI_SFLP2',1,ZHOOK_HANDLE)
END SUBROUTINE INI_SFLP2

!=========================================================================

SUBROUTINE SETUP_SFLP2(YDSC,YD,KGRIB,CDNAME,PDEFAULT,KTRAJ,KREQIN)
! Setup 2-D surface field
TYPE(TYPE_SURF_GEN),INTENT(INOUT)      :: YDSC
TYPE(TYPE_SURF_MTL_2D),INTENT(INOUT)   :: YD
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KGRIB
CHARACTER(LEN=16) ,OPTIONAL,INTENT(IN) :: CDNAME
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN) :: PDEFAULT
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KTRAJ
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KREQIN

INTEGER(KIND=JPIM) :: IPTR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SETUP_SFLP2',0,ZHOOK_HANDLE)

IPTR = YDSC%IPTR
IF(IPTR > JPMAXSFLDS-1) THEN
  ! We are about to set the last element which is generally used as a
  ! "missing field" placeholder.
  WRITE(NULERR,*) 'SURFACE FIELDS UNDER-DIMENSINED - GROUP ',YDSC%CGRPNAME,&
   &  KGRIB,CDNAME
  CALL ABOR1('IPTR > JPMAXSFLDS-1')
ENDIF

IF (YDSC%FINALISED) THEN
  WRITE(NULERR,*) 'ERROR - SETUP_SFLP2 CALLED FOR A GROUP WHICH HAS BEEN FINALISED',&
 &  YDSC%CGRPNAME,KGRIB,CDNAME
  CALL ABOR1('SETUP_SFLP2 CALLED FOR A GROUP WHICH HAS BEEN FINALISED')
ENDIF

! If the YD(:) elements are not being set in order it's a possible sign of a
! bug - the caller may have passed in a pointer to the wrong element. Or it
! might be entirely intentional and an abort is unnecessary.
! This if-block and all other references to IBUGFINDER can be commented if
! they cause a problem.
IF (YD%IBUGFINDER /= IPTR) THEN
  WRITE(NULERR,*) 'SETUP_SFLP2 NOT BEING CALLED FOR YD ELEMENTS IN ORDER. '//&
    & '-POSSIBLY- A BUG.',YDSC%CGRPNAME,KGRIB,CDNAME
  CALL ABOR1('SETUP_SFLP2 NOT CALLED FOR YD ELEMENTS IN ORDER. POSSIBLE BUG.')
ENDIF

IF(PRESENT(KGRIB)) THEN
  YD%IGRBCODE = KGRIB
ELSE
  YD%IGRBCODE = -999
ENDIF

IF(PRESENT(KREQIN)) THEN
  YD%NREQIN = KREQIN
ELSE
  YD%NREQIN = -1
ENDIF

IF(PRESENT(CDNAME)) THEN
  YD%CNAME = CDNAME
ELSE
  YD%CNAME = ''
ENDIF

IF(PRESENT(PDEFAULT)) THEN
  YD%REFVALI = PDEFAULT
ELSE
  YD%REFVALI = 0.0_JPRB
ENDIF

IF(PRESENT(KTRAJ)) THEN
  IF(KTRAJ == 1) THEN
    IF (NOFFTRAJ+1 > JPMAXSTRAJ) CALL ABOR1('JPMAXSTRAJ MUST BE INCREASED!')
    NSTRAJGRIB(NOFFTRAJ+1) = YD%IGRBCODE
    NOFFTRAJ = NOFFTRAJ+1
  ELSEIF(KTRAJ == 2) THEN
    NOFFTRAJ_CST = NOFFTRAJ_CST+1
  ELSEIF(KTRAJ /= 0) THEN
    CALL ABOR1('SURFACE_FIELDS_MIX:SETUP_SFLP2 - UNKNOWN KTRAJ')
  ENDIF
  YD%ITRAJ = KTRAJ
  YDSC%NDIM5 = YDSC%NDIM5+1
  ! MC to check with Mats
  !IF(KTRAJ /= 0 .AND. YD%NREQIN == 0) THEN
  !  WRITE(NULERR,*) 'ERROR SURFACE_FIELDS_MIX:SETUP_SFLP2 ',YD%IGRBCODE
  !  CALL ABOR1('SURFACE_FIELDS_MIX:SETUP_SFLP2 - FOR TRAJECTORY FIELDS MUST BE INITIALIZED')
  !ENDIF
ELSE
  YD%ITRAJ = 0
ENDIF

YD%LSET=.TRUE.
YDSC%IPTR = YDSC%IPTR+1

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SETUP_SFLP2',1,ZHOOK_HANDLE)
END SUBROUTINE SETUP_SFLP2

!=========================================================================

SUBROUTINE FINALISE_SFLP2(YDDYNA,YDSURF,YDSC,YD,LDVERBOSE)
! Finalise 2-D surface field group
TYPE(TDYNA), INTENT(IN)              :: YDDYNA
TYPE(TSURF), INTENT(INOUT)           :: YDSURF
TYPE(TYPE_SURF_GEN),INTENT(INOUT)    :: YDSC
TYPE(TYPE_SURF_MTL_2D),INTENT(INOUT) :: YD(:)
LOGICAL,OPTIONAL,INTENT(IN)          :: LDVERBOSE

INTEGER(KIND=JPIM) :: JFLD
LOGICAL :: LLVERBOSE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:FINALISE_SFLP2',0,ZHOOK_HANDLE)

LLVERBOSE =  .TRUE.
IF(PRESENT(LDVERBOSE)) LLVERBOSE=LDVERBOSE
YDSC%FINALISED=.TRUE.
YDSC%NUMFLDS = YDSC%IPTR-1

WRITE(NULOUT,*) '2-D SURFACE FIELD GROUP INITIALIZED ', YDSC%CGRPNAME
WRITE(NULOUT,*) 'NUMFLDS=',YDSC%NUMFLDS,' LMTL=',YDSC%LMTL

YDSURF%NSURF = YDSURF%NSURF+YDSC%NUMFLDS
YDSURF%NSURFL = YDSURF%NSURFL+YDSC%NUMFLDS
IF(YDSC%LMTL) THEN
  YDSURF%NPROGSURF = YDSURF%NPROGSURF+YDSC%NUMFLDS
  YDSURF%NPROGSURFL = YDSURF%NPROGSURFL+YDSC%NUMFLDS
ENDIF

IF(YDSC%LMTL) THEN
  IF (YDDYNA%LTWOTL) THEN
    YDSC%NDIM = 2*YDSC%NUMFLDS
  ELSE
    YDSC%NDIM = 3*YDSC%NUMFLDS
  ENDIF
ELSE
  YDSC%NDIM = YDSC%NUMFLDS
ENDIF
YDSURF%NDIMSURF = YDSURF%NDIMSURF + YDSC%NDIM
YDSURF%NDIMSURFL = YDSURF%NDIMSURFL + YDSC%NDIM

DO JFLD=1,YDSC%NUMFLDS

  YD(JFLD)%MP = JFLD
  IF (YDSC%LMTL) THEN
    YD(JFLD)%MP0 = YD(JFLD)%MP
    IF(YDDYNA%LTWOTL) THEN
      YD(JFLD)%MP9 = YD(JFLD)%MP0
      YD(JFLD)%MP1 = YD(JFLD)%MP0+YDSC%NUMFLDS
    ELSE
      YD(JFLD)%MP9 = YD(JFLD)%MP0+YDSC%NUMFLDS
      YD(JFLD)%MP1 = YD(JFLD)%MP0+2*YDSC%NUMFLDS
    ENDIF
  ELSE
    YD(JFLD)%MP0 = NUNDEFLD
    YD(JFLD)%MP9 = NUNDEFLD
    YD(JFLD)%MP1 = NUNDEFLD
  ENDIF

  IF(LLVERBOSE) THEN
    IF(YDSC%LMTL) THEN
      WRITE(NULOUT,'(1X,A,I4,1X,A,1x,I6,x,5I6)')&
       & YDSC%CGRPNAME(1:6),JFLD,YD(JFLD)%CNAME,YD(JFLD)%IGRBCODE,&
       & YD(JFLD)%MP0,YD(JFLD)%MP9,YD(JFLD)%MP1,YD(JFLD)%ITRAJ,YD(JFLD)%NREQIN
    ELSE
      WRITE(NULOUT,'(1X,A,I4,1X,A,1x,I6,x,3I6)')&
       & YDSC%CGRPNAME(1:6),JFLD,YD(JFLD)%CNAME,YD(JFLD)%IGRBCODE,&
       & YD(JFLD)%MP,YD(JFLD)%ITRAJ,YD(JFLD)%NREQIN
    ENDIF
  ENDIF

ENDDO

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:FINALISE_SFLP2',1,ZHOOK_HANDLE)
END SUBROUTINE FINALISE_SFLP2

!=========================================================================

SUBROUTINE GPPOPER(YDDYN,CDACT,YDSURF,KBL,PSP_SB,PSP_SG,PSP_SL,PSP_RR,PSP_CL,PSP_OM,PSP_EP,PSP_X2,PSP_OC,YDCOM)
! Operations on prognostic surface fields
USE YOMDYN , ONLY : TDYN
TYPE(TDYN)                  ,INTENT(INOUT) :: YDDYN
CHARACTER(LEN=*),            INTENT(IN)    :: CDACT
TYPE(TSURF),                 INTENT(INOUT) :: YDSURF
INTEGER(KIND=JPIM), OPTIONAL,INTENT(IN)    :: KBL
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_SB(:,:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_SG(:,:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_SL(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_RR(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_CL(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_OM(:,:,:)                          !KPP
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_EP(:,:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_OC(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_X2(:,:)
TYPE(TYPE_SFL_COMM),OPTIONAL,INTENT(INOUT) :: YDCOM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------


IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:GPPOPER',0,ZHOOK_HANDLE)
IF(PRESENT(KBL)) THEN
  CALL GPOPER_3(YDDYN,CDACT,YDSURF%SP_SB(:,:,:,KBL),YDSURF%YSP_SBD,YDSURF%YSP_SB%YSB,YDCOM)
  CALL GPOPER_3(YDDYN,CDACT,YDSURF%SP_SG(:,:,:,KBL),YDSURF%YSP_SGD,YDSURF%YSP_SG%YSG,YDCOM,K3D=1)
  CALL GPOPER_2(YDDYN,CDACT,YDSURF%SP_SL(:,:,KBL)  ,YDSURF%YSP_SLD,YDSURF%YSP_SL%YSL,YDCOM)
  CALL GPOPER_2(YDDYN,CDACT,YDSURF%SP_RR(:,:,KBL)  ,YDSURF%YSP_RRD,YDSURF%YSP_RR%YRR,YDCOM)
  CALL GPOPER_2(YDDYN,CDACT,YDSURF%SP_CL(:,:,KBL)  ,YDSURF%YSP_CLD,YDSURF%YSP_CL%YCL,YDCOM)
  CALL GPOPER_3(YDDYN,CDACT,YDSURF%SP_OM(:,:,:,KBL),YDSURF%YSP_OMD,YDSURF%YSP_OM%YOM,YDCOM)
  CALL GPOPER_3(YDDYN,CDACT,YDSURF%SP_EP(:,:,:,KBL),YDSURF%YSP_EPD,YDSURF%YSP_EP%YEP,YDCOM)
  CALL GPOPER_2(YDDYN,CDACT,YDSURF%SP_X2(:,:,KBL)  ,YDSURF%YSP_X2D,YDSURF%YSP_X2%YX2,YDCOM)
ELSE
  IF(PRESENT(PSP_SB)) CALL GPOPER_3(YDDYN,CDACT,PSP_SB(:,:,:),YDSURF%YSP_SBD,YDSURF%YSP_SB%YSB,YDCOM)
  IF(PRESENT(PSP_SG)) CALL GPOPER_3(YDDYN,CDACT,PSP_SG(:,:,:),YDSURF%YSP_SGD,YDSURF%YSP_SG%YSG,YDCOM,K3D=1)
  IF(PRESENT(PSP_SL)) CALL GPOPER_2(YDDYN,CDACT,PSP_SL(:,:)  ,YDSURF%YSP_SLD,YDSURF%YSP_SL%YSL,YDCOM)
  IF(PRESENT(PSP_RR)) CALL GPOPER_2(YDDYN,CDACT,PSP_RR(:,:)  ,YDSURF%YSP_RRD,YDSURF%YSP_RR%YRR,YDCOM)
  IF(PRESENT(PSP_CL)) CALL GPOPER_2(YDDYN,CDACT,PSP_CL(:,:)  ,YDSURF%YSP_CLD,YDSURF%YSP_CL%YCL,YDCOM)
  IF(PRESENT(PSP_OM)) CALL GPOPER_3(YDDYN,CDACT,PSP_OM(:,:,:),YDSURF%YSP_OMD,YDSURF%YSP_OM%YOM,YDCOM)
  IF(PRESENT(PSP_EP)) CALL GPOPER_3(YDDYN,CDACT,PSP_EP(:,:,:),YDSURF%YSP_EPD,YDSURF%YSP_EP%YEP,YDCOM)
  IF(PRESENT(PSP_X2)) CALL GPOPER_2(YDDYN,CDACT,PSP_X2(:,:)  ,YDSURF%YSP_X2D,YDSURF%YSP_X2%YX2,YDCOM)
ENDIF
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:GPPOPER',1,ZHOOK_HANDLE)
END SUBROUTINE GPPOPER

!=========================================================================

SUBROUTINE GPOPER(YDDIM,YDDYN,CDACT,YDSURF,KBL,K2D,PSP_SB,PSP_SG,PSP_SL,PSP_RR,PSP_CL,PSP_OM,&
 & PSD_VF,PSD_VD,PSD_SM,PSD_VV,PSD_WS,PSD_V2,PSD_V3,PSD_OC,YDCOM,PFIELD,PFIELD2)

!Operations on ALL surface groups

USE YOMDYN , ONLY : TDYN
TYPE(TDIM)        ,          INTENT(IN)    :: YDDIM
TYPE(TDYN)        ,          INTENT(IN)    :: YDDYN
CHARACTER(LEN=*)  ,          INTENT(IN)    :: CDACT
TYPE(TSURF)       ,          INTENT(INOUT) :: YDSURF
INTEGER(KIND=JPIM), OPTIONAL,INTENT(IN)    :: KBL
INTEGER(KIND=JPIM), OPTIONAL,INTENT(IN)    :: K2D
! pronostic groups:
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_SB(:,:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_SG(:,:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_SL(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_RR(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_CL(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSP_OM(:,:,:)
! diagnostic groups:
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSD_VF(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSD_VD(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSD_SM(:,:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSD_VV(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSD_WS(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSD_V2(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSD_V3(:,:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PSD_OC(:,:)
TYPE(TYPE_SFL_COMM),OPTIONAL,INTENT(INOUT) :: YDCOM
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PFIELD(:,:)
REAL(KIND=JPRB),    OPTIONAL,INTENT(INOUT) :: PFIELD2(:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:GPOPER',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDDIM%NPROMA)
IF(CDACT == 'PUTALLFLDS' .OR. CDACT == 'GETALLFLDS'   .OR.&
 & CDACT == 'TRAJSTORE'  .OR. CDACT == 'TRAJSTORECST' .OR.&
 & CDACT == 'SET0TOTRAJ' .OR. CDACT == 'GETTRAJ'          ) THEN
  IF(.NOT.PRESENT(PFIELD)) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER - PFIELD MISSING')
  IF(SIZE(PFIELD,1) < NPROMA)  CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER - SIZE(PFIELD,1) < NPROMA)')
ENDIF
IF(CDACT == 'PUTALLFLDS' .OR. CDACT == 'GETALLFLDS') THEN
  IF(SIZE(PFIELD,2) < YDSURF%NPROGSURFL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER - SIZE(PFIELD,2) < NPROGSURFL)')
ENDIF
IF(CDACT == 'GETTRAJ') THEN
  IF(.NOT.PRESENT(PFIELD2)) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER - PFIELD2 MISSING')
  IF(SIZE(PFIELD2,1) < NPROMA)  CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER - SIZE(PFIELD2,1) < NPROMA)')
ENDIF
IF(PRESENT(YDCOM)) THEN
  YDCOM%L_OK = .FALSE.
  YDCOM%IPTRSURF = 0
  YDCOM%ICOUNT   = 0
  YDCOM%ICOUNTML = 0
ENDIF

NPTRSURF = 0

IF(PRESENT(KBL)) THEN
  ! pronostic groups:
  IF(YDSURF%YSP_SBD%NDIM > 0) THEN
    CALL GPOPER_3(YDDYN,CDACT,YDSURF%SP_SB(:,:,:,KBL),YDSURF%YSP_SBD,YDSURF%YSP_SB%YSB,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_SGD%NDIM > 0) THEN
    IF (.NOT.PRESENT(K2D))THEN
      CALL GPOPER_3(YDDYN,CDACT,YDSURF%SP_SG(:,:,:,KBL)  ,YDSURF%YSP_SGD,YDSURF%YSP_SG%YSG,YDCOM,PFIELD,PFIELD2,K3D=1)
    ELSE
      CALL GPOPER_3(YDDYN,CDACT,YDSURF%SP_SG(:,:,:,KBL)  ,YDSURF%YSP_SGD,YDSURF%YSP_SG%YSG,YDCOM,PFIELD,PFIELD2)
    ENDIF
  ENDIF
  IF(YDSURF%YSP_SLD%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SP_SL(:,:,KBL)  ,YDSURF%YSP_SLD,YDSURF%YSP_SL%YSL,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_RRD%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SP_RR(:,:,KBL)  ,YDSURF%YSP_RRD,YDSURF%YSP_RR%YRR,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_CLD%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SP_CL(:,:,KBL)  ,YDSURF%YSP_CLD,YDSURF%YSP_CL%YCL,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_OMD%NDIM > 0) THEN !KPP
    CALL GPOPER_3(YDDYN,CDACT,YDSURF%SP_OM(:,:,:,KBL),YDSURF%YSP_OMD,YDSURF%YSP_OM%YOM,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_EPD%NDIM > 0) THEN
    CALL GPOPER_3(YDDYN,CDACT,YDSURF%SP_EP(:,:,:,KBL),YDSURF%YSP_EPD,YDSURF%YSP_EP%YEP,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_X2D%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SP_X2(:,:,KBL)  ,YDSURF%YSP_X2D,YDSURF%YSP_X2%YX2,YDCOM,PFIELD,PFIELD2)
  ENDIF
  ! diagnostic groups:
  IF(YDSURF%YSD_VFD%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_VF(:,:,KBL)  ,YDSURF%YSD_VFD,YDSURF%YSD_VF%YVF,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VPD%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_VP(:,:,KBL)  ,YDSURF%YSD_VPD,YDSURF%YSD_VP%YVP,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VVD%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_VV(:,:,KBL),YDSURF%YSD_VVD,YDSURF%YSD_VV%YVV,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VND%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_VN(:,:,KBL),YDSURF%YSD_VND,YDSURF%YSD_VN%YVN,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VHD%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_VH(:,:,KBL),YDSURF%YSD_VHD,YDSURF%YSD_VH%YVH,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VKD%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_VK(:,:,KBL),YDSURF%YSD_VKD,YDSURF%YSD_VK%YVK,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VAD%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_VA(:,:,KBL),YDSURF%YSD_VAD,YDSURF%YSD_VA%YVA,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VCD%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_VC(:,:,KBL),YDSURF%YSD_VCD,YDSURF%YSD_VC%YVC,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_V2D%NDIM > 0) THEN !KPP
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_V2(:,:,KBL),YDSURF%YSD_V2D,YDSURF%YSD_V2%YV2,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_V3D%NDIM > 0) THEN !KPP
    CALL GPOPER_3(YDDYN,CDACT,YDSURF%SD_V3(:,:,:,KBL),YDSURF%YSD_V3D,YDSURF%YSD_V3%YV3,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_OCD%NDIM > 0) THEN !KPP
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_OC(:,:,KBL)  ,YDSURF%YSD_OCD,YDSURF%YSD_OC%YOC,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VDD%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_VD(:,:,KBL),YDSURF%YSD_VDD,YDSURF%YSD_VD%YVD,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_SMD%NDIM > 0) THEN
    CALL GPOPER_3(YDDYN,CDACT,YDSURF%SD_SM(:,:,:,KBL),YDSURF%YSD_SMD,YDSURF%YSD_SM%YSM,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_WSD%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_WS(:,:,KBL),YDSURF%YSD_WSD,YDSURF%YSD_WS%YWS,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_WWD%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_WW(:,:,KBL),YDSURF%YSD_WWD,YDSURF%YSD_WW%YWW,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VXD%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_VX(:,:,KBL),YDSURF%YSD_VXD,YDSURF%YSD_VX%YVX,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_XAD%NDIM > 0) THEN
    CALL GPOPER_3(YDDYN,CDACT,YDSURF%SD_XA(:,:,:,KBL),YDSURF%YSD_XAD,YDSURF%YSD_XA%YXA,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_XRD%NDIM > 0) THEN
    CALL GPOPER_3(YDDYN,CDACT,YDSURF%SD_XR(:,:,:,KBL),YDSURF%YSD_XRD,YDSURF%YSD_XR%YXA,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_DID%NDIM > 0) THEN
    CALL GPOPER_3(YDDYN,CDACT,YDSURF%SD_DI(:,:,:,KBL),YDSURF%YSD_DID,YDSURF%YSD_DI%YXA,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_X2D%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_X2(:,:,KBL),YDSURF%YSD_X2D,YDSURF%YSD_X2%YX2,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_SFLD%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_SFL(:,:,KBL),YDSURF%YSD_SFLD,YDSURF%YSD_SFL%YSFL,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_SFOD%NDIM > 0) THEN
    CALL GPOPER_2(YDDYN,CDACT,YDSURF%SD_SFO(:,:,KBL),YDSURF%YSD_SFOD,YDSURF%YSD_SFO%YSFO,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_PFD%NDIM > 0) THEN
    CALL GPOPER_3(YDDYN,CDACT,YDSURF%SD_PF(:,:,:,KBL),YDSURF%YSD_PFD,YDSURF%YSD_PF%YXA,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_XPD%NDIM > 0) THEN
    CALL GPOPER_3(YDDYN,CDACT,YDSURF%SD_XP(:,:,:,KBL),YDSURF%YSD_XPD,YDSURF%YSD_XP%YXP1,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_XP2D%NDIM > 0) THEN
    CALL GPOPER_3(YDDYN,CDACT,YDSURF%SD_XP2(:,:,:,KBL),YDSURF%YSD_XP2D,YDSURF%YSD_XP2%YXP2,YDCOM,PFIELD,PFIELD2)
  ENDIF
ELSE
  ! pronostic groups:
  IF(YDSURF%YSP_SBD%NDIM > 0) THEN
    IF(PRESENT(PSP_SB))&
     & CALL GPOPER_3(YDDYN,CDACT,PSP_SB,YDSURF%YSP_SBD,YDSURF%YSP_SB%YSB,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_SGD%NDIM > 0) THEN
    IF(PRESENT(PSP_SG))&
     & CALL GPOPER_3(YDDYN,CDACT,PSP_SG,YDSURF%YSP_SGD,YDSURF%YSP_SG%YSG,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_SLD%NDIM > 0) THEN             !FLAKE
    IF(PRESENT(PSP_SL))&
     & CALL GPOPER_2(YDDYN,CDACT,PSP_SL,YDSURF%YSP_SLD,YDSURF%YSP_SL%YSL,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_RRD%NDIM > 0) THEN
    IF(PRESENT(PSP_RR))&
     & CALL GPOPER_2(YDDYN,CDACT,PSP_RR,YDSURF%YSP_RRD,YDSURF%YSP_RR%YRR,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_CLD%NDIM > 0) THEN
    IF(PRESENT(PSP_CL))&
     & CALL GPOPER_2(YDDYN,CDACT,PSP_CL,YDSURF%YSP_CLD,YDSURF%YSP_CL%YCL,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSP_OMD%NDIM > 0) THEN
    IF(PRESENT(PSP_OM))&
    & CALL GPOPER_3(YDDYN,CDACT,PSP_OM,YDSURF%YSP_OMD,YDSURF%YSP_OM%YOM,YDCOM,PFIELD,PFIELD2)
  ENDIF
  ! diagnostic groups:
  IF(YDSURF%YSD_VFD%NDIM > 0) THEN
    IF(PRESENT(PSD_VF))&
     & CALL GPOPER_2(YDDYN,CDACT,PSD_VF,YDSURF%YSD_VFD,YDSURF%YSD_VF%YVF,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VVD%NDIM > 0) THEN
    IF(PRESENT(PSD_VV))&
     & CALL GPOPER_2(YDDYN,CDACT,PSD_VV,YDSURF%YSD_VVD,YDSURF%YSD_VV%YVV,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_VDD%NDIM > 0) THEN
    IF(PRESENT(PSD_VD))&
     & CALL GPOPER_2(YDDYN,CDACT,PSD_VD,YDSURF%YSD_VDD,YDSURF%YSD_VD%YVD,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_SMD%NDIM > 0) THEN
    IF(PRESENT(PSD_SM))&
    & CALL GPOPER_3(YDDYN,CDACT,PSD_SM,YDSURF%YSD_SMD,YDSURF%YSD_SM%YSM,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_WSD%NDIM > 0) THEN
    IF(PRESENT(PSD_WS))&
     & CALL GPOPER_2(YDDYN,CDACT,PSD_WS,YDSURF%YSD_WSD,YDSURF%YSD_WS%YWS,YDCOM,PFIELD,PFIELD2)
  ENDIF
  IF(YDSURF%YSD_OCD%NDIM > 0) THEN
    IF(PRESENT(PSD_OC))&
     & CALL GPOPER_2(YDDYN,CDACT,PSD_OC,YDSURF%YSD_OCD,YDSURF%YSD_OC%YOC,YDCOM,PFIELD,PFIELD2)
  ENDIF
ENDIF
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:GPOPER',1,ZHOOK_HANDLE)
END SUBROUTINE GPOPER

!=========================================================================

SUBROUTINE GPOPER_2(YDDYN,CDACT,PFLD,YDSC,YD,YDCOM,PFIELD,PFIELD2)
! Operations on 2-D surface groups
USE YOMDYN , ONLY : TDYN
TYPE(TDYN)      ,INTENT(IN)                :: YDDYN
CHARACTER(LEN=*),INTENT(IN)                :: CDACT
REAL(KIND=JPRB) ,INTENT(INOUT)             :: PFLD(:,:)
TYPE(TYPE_SURF_GEN),INTENT(IN)             :: YDSC
TYPE(TYPE_SURF_MTL_2D),INTENT(IN)          :: YD(:)
TYPE(TYPE_SFL_COMM),OPTIONAL,INTENT(INOUT) :: YDCOM
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT)     :: PFIELD(:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT)     :: PFIELD2(:,:)

INTEGER(KIND=JPIM) :: J,IPTR,IPTR2
REAL(KIND=JPRB) :: ZZPHY
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:GPOPER_2',0,ZHOOK_HANDLE)
ASSOCIATE(REPSP1=>YDDYN%REPSP1)
IF(CDACT == 'SET9TO0') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP9) = PFLD(:,YD(J)%MP0)
  ENDDO
ELSEIF(CDACT == 'SET1TO0') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP1) = PFLD(:,YD(J)%MP0)
  ENDDO
ELSEIF(CDACT == 'SET1TO9') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP1) = PFLD(:,YD(J)%MP9)
  ENDDO
ELSEIF(CDACT == 'SET1TO9AD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP0) = PFLD(:,YD(J)%MP9)+PFLD(:,YD(J)%MP1)
    PFLD(:,YD(J)%MP1) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET0TO1') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP0) = PFLD(:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'SET0TO1AD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP1) = PFLD(:,YD(J)%MP1)+PFLD(:,YD(J)%MP0)
    PFLD(:,YD(J)%MP0) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET9TO1') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP9) = PFLD(:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'PHTFILT') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  ZZPHY=1.0_JPRB-REPSP1
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP9) = REPSP1*PFLD(:,YD(J)%MP1)+ZZPHY*PFLD(:,YD(J)%MP0)
    PFLD(:,YD(J)%MP0) = PFLD(:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'PHTFILTAD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  ZZPHY=1.0_JPRB-REPSP1
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP1) = PFLD(:,YD(J)%MP1)+PFLD(:,YD(J)%MP0)
    PFLD(:,YD(J)%MP0) = 0.0_JPRB
    PFLD(:,YD(J)%MP1) = PFLD(:,YD(J)%MP1)+REPSP1*PFLD(:,YD(J)%MP9)
    PFLD(:,YD(J)%MP0) = PFLD(:,YD(J)%MP0)+ZZPHY *PFLD(:,YD(J)%MP9)
    PFLD(:,YD(J)%MP9) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET0TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP0) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SET9TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP9) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SET1TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_2 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP1) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SETALLTOVAL') THEN
  DO J=1,YDSC%NDIM
    PFLD(:,J) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SETPERTOVAL') THEN
  DO J=1,YDSC%NUMFLDS
    PFLD(:,YD(J)%MP) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SETDEFAULT') THEN
  DO J=1,YDSC%NUMFLDS
    IF(YD(J)%NREQIN == -1) THEN
      PFLD(:,YD(J)%MP) = YD(J)%REFVALI
    ENDIF
  ENDDO
ELSEIF(CDACT == 'TRAJSTORE') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        IPTR = IPTR+1
        PFIELD(:,IPTR) = PFLD(:,YD(J)%MP)
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'TRAJSTORECST') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR2 = YDSC%NOFFTRAJ_CST
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 2) THEN
        IPTR2 = IPTR2+1
        PFIELD(:,IPTR2) = PFLD(:,YD(J)%MP)
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'SET0TOTRAJ') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        IPTR = IPTR+1
        PFLD(:,YD(J)%MP) = PFIELD(:,IPTR)
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'GETTRAJ') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    IPTR2 = YDSC%NOFFTRAJ_CST
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        IPTR = IPTR+1
        PFLD(:,YD(J)%MP) = PFIELD(:,IPTR)
      ELSEIF(YD(J)%ITRAJ == 2) THEN
        IPTR2 = IPTR2+1
        PFLD(:,YD(J)%MP) = PFIELD2(:,IPTR2)
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'GETALLFLDS') THEN
  DO J=1,YDSC%NDIM
    NPTRSURF = NPTRSURF+1
    PFIELD(:,NPTRSURF) = PFLD(:,J)
  ENDDO
ELSEIF(CDACT == 'PUTALLFLDS') THEN
  DO J=1,YDSC%NDIM
    NPTRSURF = NPTRSURF+1
    PFLD(:,J) = PFIELD(:,NPTRSURF)
  ENDDO
ELSEIF(CDACT == 'GETGRIBPOS') THEN
  DO J=1,YDSC%NUMFLDS
    YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
    IF(YD(J)%IGRBCODE == YDCOM%IGRBCODE) THEN
      YDCOM%IFLDNUM  = YDCOM%IPTRSURF
      YDCOM%L_OK = .TRUE.
    ENDIF
  ENDDO
ELSEIF(CDACT == 'GETFIELD') THEN
  DO J=1,YDSC%NUMFLDS
    YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
    IF(YDCOM%IPTRSURF == YDCOM%IFLDNUM) THEN
      PFIELD(:,1) = PFLD(:,J)
      YDCOM%L_OK = .TRUE.
    ENDIF
  ENDDO
ELSEIF(CDACT == 'GRIBIN') THEN
  DO J=1,YDSC%NUMFLDS
    YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
    IF(YD(J)%NREQIN == 1) THEN
      YDCOM%ICOUNT = YDCOM%ICOUNT+1
      YDCOM%ICODES(YDCOM%ICOUNT) = YD(J)%IGRBCODE
    ENDIF
  ENDDO
ELSE
  WRITE(NULOUT,*) 'SURFACE_FIELD:GPOPER_2 UNKNOWN ACTION - ',CDACT
  CALL ABOR1('SURFACE_FIELD:GPOPER_2 - UNKNOWN ACTION')
ENDIF
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:GPOPER_2',1,ZHOOK_HANDLE)
END SUBROUTINE GPOPER_2

!=========================================================================

SUBROUTINE GPOPER_3(YDDYN,CDACT,PFLD,YDSC,YD,YDCOM,PFIELD,PFIELD2, K3D)
! Operations on 3-D surface groups
USE YOMDYN , ONLY : TDYN
TYPE(TDYN)      ,INTENT(IN)                :: YDDYN
CHARACTER(LEN=*),INTENT(IN)                :: CDACT
REAL(KIND=JPRB) ,INTENT(INOUT)             :: PFLD(:,:,:)
TYPE(TYPE_SURF_GEN),INTENT(IN)             :: YDSC
TYPE(TYPE_SURF_MTL_3D),INTENT(IN)          :: YD(:)
TYPE(TYPE_SFL_COMM),OPTIONAL,INTENT(INOUT) :: YDCOM
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT)     :: PFIELD(:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT)     :: PFIELD2(:,:)
INTEGER,        OPTIONAL,INTENT(IN)        :: K3D

INTEGER(KIND=JPIM) :: J,JLEV,IPTR,IPTR2
REAL(KIND=JPRB) :: ZZPHY
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:GPOPER_3',0,ZHOOK_HANDLE)
ASSOCIATE(REPSP1=>YDDYN%REPSP1)
IF(CDACT == 'SET9TO0') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP9) = PFLD(:,:,YD(J)%MP0)
  ENDDO
ELSEIF(CDACT == 'SET1TO0') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP1) = PFLD(:,:,YD(J)%MP0)
  ENDDO
ELSEIF(CDACT == 'SET1TO9') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP1) = PFLD(:,:,YD(J)%MP9)
  ENDDO
ELSEIF(CDACT == 'SET1TO9AD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP9) = PFLD(:,:,YD(J)%MP9)+PFLD(:,:,YD(J)%MP1)
    PFLD(:,:,YD(J)%MP1) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET0TO1') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP0) = PFLD(:,:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'SET0TO1AD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP1) = PFLD(:,:,YD(J)%MP1)+PFLD(:,:,YD(J)%MP0)
    PFLD(:,:,YD(J)%MP0) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET9TO1') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP9) = PFLD(:,:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'PHTFILT') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  ZZPHY=1.0_JPRB-REPSP1
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP9) = REPSP1*PFLD(:,:,YD(J)%MP1)+ZZPHY*PFLD(:,:,YD(J)%MP0)
    PFLD(:,:,YD(J)%MP0) = PFLD(:,:,YD(J)%MP1)
  ENDDO
ELSEIF(CDACT == 'PHTFILTAD') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  ZZPHY=1.0_JPRB-REPSP1
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP1) = PFLD(:,:,YD(J)%MP1)+PFLD(:,:,YD(J)%MP0)
    PFLD(:,:,YD(J)%MP0) = 0.0_JPRB
    PFLD(:,:,YD(J)%MP1) = PFLD(:,:,YD(J)%MP1)+REPSP1*PFLD(:,:,YD(J)%MP9)
    PFLD(:,:,YD(J)%MP0) = PFLD(:,:,YD(J)%MP0)+ZZPHY *PFLD(:,:,YD(J)%MP9)
    PFLD(:,:,YD(J)%MP9) = 0.0_JPRB
  ENDDO
ELSEIF(CDACT == 'SET0TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP0) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SET9TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP9) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SET1TOVAL') THEN
  IF( .NOT. YDSC%LMTL) CALL ABOR1('SURFACE_FIELDS_MIX:GPOPER_3 : FIELD NOT MULTI-TIME LEVEL')
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP1) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SETALLTOVAL') THEN
  DO J=1,YDSC%NDIM
    PFLD(:,:,J) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SETPERTOVAL') THEN
  DO J=1,YDSC%NUMFLDS
    PFLD(:,:,YD(J)%MP) = YDCOM%VALUE
  ENDDO
ELSEIF(CDACT == 'SETDEFAULT') THEN
  DO J=1,YDSC%NUMFLDS
    DO JLEV=1,YDSC%NLEVS
      IF(YD(J)%NREQIN(JLEV) == -1) THEN
        PFLD(:,JLEV,YD(J)%MP) = YD(J)%REFVALI(JLEV)
      ENDIF
    ENDDO
  ENDDO
ELSEIF(CDACT == 'TRAJSTORE') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        DO JLEV=1,YDSC%NLEVS
          IPTR = IPTR+1
          PFIELD(:,IPTR) = PFLD(:,JLEV,YD(J)%MP)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'TRAJSTORECST') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR2 = YDSC%NOFFTRAJ_CST
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 2) THEN
        DO JLEV=1,YDSC%NLEVS
          IPTR2 = IPTR2+1
          PFIELD(:,IPTR2) = PFLD(:,JLEV,YD(J)%MP)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'SET0TOTRAJ') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        DO JLEV=1,YDSC%NLEVS
          IPTR = IPTR+1
          PFLD(:,JLEV,YD(J)%MP) = PFIELD(:,IPTR)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'GETTRAJ') THEN
  IF(YDSC%NDIM5 > 0 ) THEN
    IPTR = YDSC%NOFFTRAJ
    IPTR2 = YDSC%NOFFTRAJ_CST
    DO J=1,YDSC%NUMFLDS
      IF(YD(J)%ITRAJ == 1) THEN
        DO JLEV=1,YDSC%NLEVS
          IPTR = IPTR+1
          PFLD(:,JLEV,YD(J)%MP) = PFIELD(:,IPTR)
        ENDDO
      ELSEIF(YD(J)%ITRAJ == 2) THEN
        DO JLEV=1,YDSC%NLEVS
          IPTR2 = IPTR2+1
          PFLD(:,JLEV,YD(J)%MP) = PFIELD2(:,IPTR2)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
ELSEIF(CDACT == 'GETALLFLDS') THEN
  DO J=1,YDSC%NDIM
    DO JLEV=1,YDSC%NLEVS
      NPTRSURF = NPTRSURF+1
      PFIELD(:,NPTRSURF) = PFLD(:,JLEV,J)
    ENDDO
  ENDDO
ELSEIF(CDACT == 'PUTALLFLDS') THEN
  DO J=1,YDSC%NDIM
    DO JLEV=1,YDSC%NLEVS
      NPTRSURF = NPTRSURF+1
      PFLD(:,JLEV,J) = PFIELD(:,NPTRSURF)
    ENDDO
  ENDDO
ELSEIF(CDACT == 'GETGRIBPOS') THEN
  DO J=1,YDSC%NUMFLDS
    DO JLEV=1,YDSC%NLEVS
      YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
      IF(YD(J)%IGRBCODE(JLEV) == YDCOM%IGRBCODE) THEN
        YDCOM%IFLDNUM  = YDCOM%IPTRSURF
        YDCOM%L_OK = .TRUE.
      ENDIF
    ENDDO
  ENDDO
ELSEIF(CDACT == 'GETFIELD') THEN
  DO J=1,YDSC%NUMFLDS
    DO JLEV=1,YDSC%NLEVS
      YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
      IF(YDCOM%IPTRSURF == YDCOM%IFLDNUM) THEN
        PFIELD(:,1) = PFLD(:,JLEV,J)
        YDCOM%L_OK = .TRUE.
      ENDIF
    ENDDO
  ENDDO
ELSEIF(CDACT == 'GRIBIN') THEN
  IF (.NOT.PRESENT(K3D)) THEN 
    DO J=1,YDSC%NUMFLDS
      DO JLEV=1,YDSC%NLEVS
        YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
        IF(YD(J)%NREQIN(JLEV) == 1) THEN
          YDCOM%ICOUNT = YDCOM%ICOUNT+1
          YDCOM%ICODES(YDCOM%ICOUNT) = YD(J)%IGRBCODE(JLEV)
        ENDIF
      ENDDO
    ENDDO
  ELSE !take only first grib code, same for all NLEVS
    DO J=1,YDSC%NUMFLDS  
      DO JLEV=1,YDSC%NLEVS
        YDCOM%IPTRSURF = YDCOM%IPTRSURF+1
      ENDDO
      IF(YD(J)%NREQIN(1) == 1) THEN
        YDCOM%ICOUNTML = YDCOM%ICOUNTML+1
        YDCOM%ICODESML(YDCOM%ICOUNTML) = YD(J)%IGRBCODE(1)
      ENDIF
    ENDDO
  ENDIF

ELSE
  WRITE(NULOUT,*) 'SURFACE_FIELD:GPOPER_3 UNKNOWN ACTION - ',CDACT
  CALL ABOR1('SURFACE_FIELD:GPOPER_3 - UNKNOWN ACTION')
ENDIF
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:GPOPER_3',1,ZHOOK_HANDLE)
END SUBROUTINE GPOPER_3

!=========================================================================

SUBROUTINE SURF_STORE(YDDIM,YDSURF,YDDYN)
! Store all surface fields
USE YOMDYN , ONLY : TDYN
TYPE(TDIM) , INTENT(IN)    :: YDDIM
TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(TDYN) , INTENT(IN)    :: YDDYN

INTEGER(KIND=JPIM) :: JBL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SURF_STORE',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, NGPBLKS=>YDDIM%NGPBLKS)
ALLOCATE(YDSURF%STORE_ARRAY(NPROMA,YDSURF%NDIMSURFL,NGPBLKS))
DO JBL=1,NGPBLKS
  CALL GPOPER(YDDIM,YDDYN,'GETALLFLDS',YDSURF,KBL=JBL,PFIELD=YDSURF%STORE_ARRAY(:,:,JBL))
ENDDO
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SURF_STORE',1,ZHOOK_HANDLE)
END SUBROUTINE SURF_STORE

!=========================================================================

SUBROUTINE SURF_RESTORE(YDDIM,YDSURF,YDDYN)
! Restore all surface fields
USE YOMDYN , ONLY : TDYN
TYPE(TDIM) , INTENT(IN)    :: YDDIM
TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(TDYN) , INTENT(IN)    :: YDDYN

INTEGER(KIND=JPIM) :: JBL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SURF_RESTORE',0,ZHOOK_HANDLE)
ASSOCIATE(NGPBLKS=>YDDIM%NGPBLKS)
IF(.NOT. ALLOCATED(YDSURF%STORE_ARRAY))&
 & CALL ABOR1('SURFACE_FIELDS_MIX:SURF_RESTORE - SURF_STORE NOT ALLOCATED')
DO JBL=1,NGPBLKS
  CALL GPOPER(YDDIM,YDDYN,'PUTALLFLDS',YDSURF,KBL=JBL,PFIELD=YDSURF%STORE_ARRAY(:,:,JBL))
ENDDO
DEALLOCATE(YDSURF%STORE_ARRAY)
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SURF_RESTORE',1,ZHOOK_HANDLE)

END SUBROUTINE SURF_RESTORE

!=========================================================================

SUBROUTINE ALLO_SURF(YDDIM,YDSURF)
! Allocate surface field arrays
!!TYPE(TSURF), POINTER, INTENT(INOUT) :: YDSURF
TYPE(TDIM) , INTENT(IN)    :: YDDIM
TYPE(TSURF), INTENT(INOUT) :: YDSURF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:ALLO_SURF',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, NGPBLKS=>YDDIM%NGPBLKS)

! pronostic groups:
ALLOCATE(YDSURF%SP_SB(NPROMA,YDSURF%YSP_SBD%NLEVS,YDSURF%YSP_SBD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SP_SG(NPROMA,YDSURF%YSP_SGD%NLEVS,YDSURF%YSP_SGD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SP_SL(NPROMA,YDSURF%YSP_SLD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SP_RR(NPROMA,YDSURF%YSP_RRD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SP_CL(NPROMA,YDSURF%YSP_CLD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SP_OM(NPROMA,YDSURF%YSP_OMD%NLEVS,YDSURF%YSP_OMD%NDIM,NGPBLKS)) !KPP
ALLOCATE(YDSURF%SP_EP(NPROMA,YDSURF%YSP_EPD%NLEVS,YDSURF%YSP_EPD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SP_X2(NPROMA,YDSURF%YSP_X2D%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SP_CI(NPROMA,YDSURF%YSP_CID%NDIM,NGPBLKS))
! diagnostic groups:
ALLOCATE(YDSURF%SD_VF(NPROMA,YDSURF%YSD_VFD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_VP(NPROMA,YDSURF%YSD_VPD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_VV(NPROMA,YDSURF%YSD_VVD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_VN(NPROMA,YDSURF%YSD_VND%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_VH(NPROMA,YDSURF%YSD_VHD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_VK(NPROMA,YDSURF%YSD_VKD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_VA(NPROMA,YDSURF%YSD_VAD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_VC(NPROMA,YDSURF%YSD_VCD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_V2(NPROMA,YDSURF%YSD_V2D%NDIM,NGPBLKS)) !KPP
ALLOCATE(YDSURF%SD_V3(NPROMA,YDSURF%YSD_V3D%NLEVS,YDSURF%YSD_V3D%NDIM,NGPBLKS)) !KPP
ALLOCATE(YDSURF%SD_VD(NPROMA,YDSURF%YSD_VDD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_SM(NPROMA,YDSURF%YSD_SMD%NLEVS,YDSURF%YSD_SMD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_WS(NPROMA,YDSURF%YSD_WSD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_WW(NPROMA,YDSURF%YSD_WWD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_VX(NPROMA,YDSURF%YSD_VXD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_XA(NPROMA,YDSURF%YSD_XAD%NLEVS,YDSURF%YSD_XAD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_XR(NPROMA,YDSURF%YSD_XRD%NLEVS,YDSURF%YSD_XRD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_X2(NPROMA,YDSURF%YSD_X2D%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_SFL(NPROMA,YDSURF%YSD_SFLD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_SFO(NPROMA,YDSURF%YSD_SFOD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_PF(NPROMA,YDSURF%YSD_PFD%NLEVS,YDSURF%YSD_PFD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_DI(NPROMA,YDSURF%YSD_DID%NLEVS,YDSURF%YSD_DID%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_OC(NPROMA,YDSURF%YSD_OCD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_XP(NPROMA,YDSURF%YSD_XPD%NLEVS,YDSURF%YSD_XPD%NDIM,NGPBLKS))
ALLOCATE(YDSURF%SD_XP2(NPROMA,YDSURF%YSD_XP2D%NLEVS,YDSURF%YSD_XP2D%NDIM,NGPBLKS))

IF(ALLOCATED(YDSURF%SP_SL)) YDSURF%SP_SL = 0.0_JPRB
IF(ALLOCATED(YDSURF%SP_RR)) YDSURF%SP_RR = 0.0_JPRB
IF(ALLOCATED(YDSURF%SP_SG)) YDSURF%SP_SG = 0.0_JPRB
IF(ALLOCATED(YDSURF%SD_VF)) YDSURF%SD_VF = 0.0_JPRB

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:ALLO_SURF',1,ZHOOK_HANDLE)
END SUBROUTINE ALLO_SURF

!=========================================================================

SUBROUTINE DEALLO_SURF(YDSURF)
! Deallocate surface field arrays

IMPLICIT NONE
TYPE(TSURF), INTENT(INOUT) :: YDSURF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:DEALLO_SURF',0,ZHOOK_HANDLE)

! pronostic groups:
DEALLOCATE(YDSURF%SP_SB)
DEALLOCATE(YDSURF%SP_SG)
DEALLOCATE(YDSURF%SP_SL)
DEALLOCATE(YDSURF%SP_RR)
DEALLOCATE(YDSURF%SP_CL)
DEALLOCATE(YDSURF%SP_OM)
DEALLOCATE(YDSURF%SP_EP)
DEALLOCATE(YDSURF%SP_X2)
DEALLOCATE(YDSURF%SP_CI)
! diagnostic groups:
DEALLOCATE(YDSURF%SD_VF)
DEALLOCATE(YDSURF%SD_VP)
DEALLOCATE(YDSURF%SD_VV)
DEALLOCATE(YDSURF%SD_VN)
DEALLOCATE(YDSURF%SD_VH)
DEALLOCATE(YDSURF%SD_VK)
DEALLOCATE(YDSURF%SD_VA)
DEALLOCATE(YDSURF%SD_VC)
DEALLOCATE(YDSURF%SD_V2)
DEALLOCATE(YDSURF%SD_V3)
DEALLOCATE(YDSURF%SD_VD)
DEALLOCATE(YDSURF%SD_SM)
DEALLOCATE(YDSURF%SD_WS)
DEALLOCATE(YDSURF%SD_WW)
DEALLOCATE(YDSURF%SD_VX)
DEALLOCATE(YDSURF%SD_XA)
DEALLOCATE(YDSURF%SD_XR)
DEALLOCATE(YDSURF%SD_X2)
DEALLOCATE(YDSURF%SD_SFL)
DEALLOCATE(YDSURF%SD_SFO)
DEALLOCATE(YDSURF%SD_PF)
DEALLOCATE(YDSURF%SD_DI)
DEALLOCATE(YDSURF%SD_OC)
DEALLOCATE(YDSURF%SD_XP)
DEALLOCATE(YDSURF%SD_XP2)

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:DEALLO_SURF',1,ZHOOK_HANDLE)

END SUBROUTINE DEALLO_SURF

!=========================================================================

SUBROUTINE ZERO_SURF(SELF)

IMPLICIT NONE
TYPE(TSURF), INTENT(INOUT) :: SELF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:ZERO_SURF',0,ZHOOK_HANDLE)

! pronostic groups:
IF (ALLOCATED(SELF%SP_SB))  SELF%SP_SB  = 0.0_JPRB
IF (ALLOCATED(SELF%SP_SG))  SELF%SP_SG  = 0.0_JPRB
IF (ALLOCATED(SELF%SP_SL))  SELF%SP_SL  = 0.0_JPRB
IF (ALLOCATED(SELF%SP_RR))  SELF%SP_RR  = 0.0_JPRB
IF (ALLOCATED(SELF%SP_CL))  SELF%SP_CL  = 0.0_JPRB
IF (ALLOCATED(SELF%SP_OM))  SELF%SP_OM  = 0.0_JPRB
IF (ALLOCATED(SELF%SP_EP))  SELF%SP_EP  = 0.0_JPRB
IF (ALLOCATED(SELF%SP_X2))  SELF%SP_X2  = 0.0_JPRB
IF (ALLOCATED(SELF%SP_CI))  SELF%SP_CI  = 0.0_JPRB

! diagnostic groups:
IF (ALLOCATED(SELF%SD_VF))  SELF%SD_VF  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VP))  SELF%SD_VP  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VV))  SELF%SD_VV  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VN))  SELF%SD_VN  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VH))  SELF%SD_VH  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VK))  SELF%SD_VK  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VA))  SELF%SD_VA  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VC))  SELF%SD_VC  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_V2))  SELF%SD_V2  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_V3))  SELF%SD_V3  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_OC))  SELF%SD_OC  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VD))  SELF%SD_VD  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_SM))  SELF%SD_SM  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_WS))  SELF%SD_WS  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_WW))  SELF%SD_WW  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_VX))  SELF%SD_VX  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_XA))  SELF%SD_XA  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_DI))  SELF%SD_DI  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_XR))  SELF%SD_XR  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_X2))  SELF%SD_X2  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_SFL)) SELF%SD_SFL = 0.0_JPRB
IF (ALLOCATED(SELF%SD_SFO)) SELF%SD_SFO = 0.0_JPRB
IF (ALLOCATED(SELF%SD_PF))  SELF%SD_PF  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_XP))  SELF%SD_XP  = 0.0_JPRB
IF (ALLOCATED(SELF%SD_XP2)) SELF%SD_XP2 = 0.0_JPRB

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:ZERO_SURF',1,ZHOOK_HANDLE)

END SUBROUTINE ZERO_SURF

!=========================================================================

SUBROUTINE COPY_CTOR_SURF(SELF,RHS)

IMPLICIT NONE
TYPE(TSURF), INTENT(INOUT) :: SELF
TYPE(TSURF), INTENT(IN)    :: RHS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:COPY_CTOR_SURF',0,ZHOOK_HANDLE)

! A straight copy of the whole structure is not possible
! because of the existence of pointers inside. REK
! example : copy_surf(A,B) ; delete(B) : then the pointers of A are lost.

! pronostic groups:
IF (ALLOCATED(RHS%SP_SB)) THEN
  SELF%SP_SB = RHS%SP_SB
ENDIF

IF (ALLOCATED(RHS%SP_SG)) THEN
  SELF%SP_SG = RHS%SP_SG
ENDIF

IF (ALLOCATED(RHS%SP_SL)) THEN
  SELF%SP_SL = RHS%SP_SL
ENDIF

IF (ALLOCATED(RHS%SP_RR)) THEN
  SELF%SP_RR = RHS%SP_RR
ENDIF

IF (ALLOCATED(RHS%SP_CL)) THEN
  SELF%SP_CL = RHS%SP_CL
ENDIF

IF (ALLOCATED(RHS%SP_OM)) THEN
  SELF%SP_OM = RHS%SP_OM
ENDIF

IF (ALLOCATED(RHS%SP_EP)) THEN
  SELF%SP_EP = RHS%SP_EP
ENDIF

IF (ALLOCATED(RHS%SP_X2)) THEN
  SELF%SP_X2 = RHS%SP_X2
ENDIF

IF (ALLOCATED(RHS%SP_CI)) THEN
  SELF%SP_CI = RHS%SP_CI
ENDIF

! diagnostic groups:
IF (ALLOCATED(RHS%SD_VF)) THEN
  SELF%SD_VF= RHS%SD_VF
ENDIF

IF (ALLOCATED(RHS%SD_VP)) THEN
  SELF%SD_VP= RHS%SD_VP
ENDIF

IF (ALLOCATED(RHS%SD_VV)) THEN
  SELF%SD_VV = RHS%SD_VV
ENDIF

IF (ALLOCATED(RHS%SD_VN)) THEN
  SELF%SD_VN = RHS%SD_VN
ENDIF

IF (ALLOCATED(RHS%SD_VH)) THEN
  SELF%SD_VH = RHS%SD_VH
ENDIF

IF (ALLOCATED(RHS%SD_VK)) THEN
  SELF%SD_VK = RHS%SD_VK
ENDIF

IF (ALLOCATED(RHS%SD_VA)) THEN
  SELF%SD_VA = RHS%SD_VA
ENDIF

IF (ALLOCATED(RHS%SD_VC)) THEN
  SELF%SD_VC = RHS%SD_VC
ENDIF

IF (ALLOCATED(RHS%SD_V2)) THEN
  SELF%SD_V2 = RHS%SD_V2
ENDIF

IF (ALLOCATED(RHS%SD_V3)) THEN
  SELF%SD_V3 = RHS%SD_V3
ENDIF

IF (ALLOCATED(RHS%SD_OC)) THEN
  SELF%SD_OC = RHS%SD_OC
ENDIF

IF (ALLOCATED(RHS%SD_VD)) THEN
  SELF%SD_VD = RHS%SD_VD
ENDIF

IF (ALLOCATED(RHS%SD_WS)) THEN
  SELF%SD_WS = RHS%SD_WS
ENDIF

IF (ALLOCATED(RHS%SD_WW)) THEN
  SELF%SD_WW = RHS%SD_WW
ENDIF

IF (ALLOCATED(RHS%SD_VX)) THEN
  SELF%SD_VX = RHS%SD_VX
ENDIF

IF (ALLOCATED(RHS%SD_XA)) THEN
  SELF%SD_XA = RHS%SD_XA
ENDIF

IF (ALLOCATED(RHS%SD_DI)) THEN
  SELF%SD_DI = RHS%SD_DI
ENDIF

IF (ALLOCATED(RHS%SD_XR)) THEN
  SELF%SD_XR = RHS%SD_XR
ENDIF

IF (ALLOCATED(RHS%SD_X2)) THEN
  SELF%SD_X2 = RHS%SD_X2
ENDIF

IF (ALLOCATED(RHS%SD_SFL)) THEN
  SELF%SD_SFL = RHS%SD_SFL
ENDIF

IF (ALLOCATED(RHS%SD_SFO)) THEN
  SELF%SD_SFO = RHS%SD_SFO
ENDIF

IF (ALLOCATED(RHS%SD_PF)) THEN
  SELF%SD_PF = RHS%SD_PF
ENDIF


IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:COPY_CTOR_SURF',1,ZHOOK_HANDLE)

END SUBROUTINE COPY_CTOR_SURF

!=========================================================================

SUBROUTINE COPY_SURF(SELF,RHS)

IMPLICIT NONE
TYPE(TSURF), INTENT(INOUT) :: SELF
TYPE(TSURF), INTENT(IN)    :: RHS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:COPY_SURF',0,ZHOOK_HANDLE)

! A straight copy of the whole structure is not possible
! because of the existence of pointers inside. REK
! example : copy_surf(A,B) ; delete(B) : then the pointers of A are lost.

! pronostic groups:
IF (ALLOCATED(SELF%SP_SB) .AND. ALLOCATED(RHS%SP_SB)) THEN
  IF (ALL(SHAPE(SELF%SP_SB) == SHAPE(RHS%SP_SB))) THEN
    SELF%SP_SB = RHS%SP_SB
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_SB different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_SG) .AND. ALLOCATED(RHS%SP_SG)) THEN
  IF (ALL(SHAPE(SELF%SP_SG) == SHAPE(RHS%SP_SG))) THEN
    SELF%SP_SG = RHS%SP_SG
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_SG different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_SL) .AND. ALLOCATED(RHS%SP_SL)) THEN
  IF (ALL(SHAPE(SELF%SP_SL) == SHAPE(RHS%SP_SL))) THEN
    SELF%SP_SL = RHS%SP_SL
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_SL different shapes")
  ENDIF
ENDIF


IF (ALLOCATED(SELF%SP_RR) .AND. ALLOCATED(RHS%SP_RR)) THEN
  IF (ALL(SHAPE(SELF%SP_RR) == SHAPE(RHS%SP_RR))) THEN
    SELF%SP_RR = RHS%SP_RR
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_RR different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_CL) .AND. ALLOCATED(RHS%SP_CL)) THEN
  IF (ALL(SHAPE(SELF%SP_CL) == SHAPE(RHS%SP_CL))) THEN
    SELF%SP_CL = RHS%SP_CL
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_CL different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_OM) .AND. ALLOCATED(RHS%SP_OM)) THEN
  IF (ALL(SHAPE(SELF%SP_OM) == SHAPE(RHS%SP_OM))) THEN
    SELF%SP_OM = RHS%SP_OM
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_OM different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_EP) .AND. ALLOCATED(RHS%SP_EP)) THEN
  IF (ALL(SHAPE(SELF%SP_EP) == SHAPE(RHS%SP_EP))) THEN
    SELF%SP_EP = RHS%SP_EP
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_EP different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_X2) .AND. ALLOCATED(RHS%SP_X2)) THEN
  IF (ALL(SHAPE(SELF%SP_X2) == SHAPE(RHS%SP_X2))) THEN
    SELF%SP_X2 = RHS%SP_X2
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_X2 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_CI) .AND. ALLOCATED(RHS%SP_CI)) THEN
  IF (ALL(SHAPE(SELF%SP_CI) == SHAPE(RHS%SP_CI))) THEN
    SELF%SP_CI = RHS%SP_CI
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SP_CI different shapes")
  ENDIF
ENDIF

! diagnostic groups:
IF (ALLOCATED(SELF%SD_VF) .AND. ALLOCATED(RHS%SD_VF)) THEN
  IF (ALL(SHAPE(SELF%SD_VF) == SHAPE(RHS%SD_VF))) THEN
    SELF%SD_VF= RHS%SD_VF
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VF different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VP) .AND. ALLOCATED(RHS%SD_VP)) THEN
  IF (ALL(SHAPE(SELF%SD_VP) == SHAPE(RHS%SD_VP))) THEN
    SELF%SD_VP= RHS%SD_VP
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VP different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VV) .AND. ALLOCATED(RHS%SD_VV)) THEN
  IF (ALL(SHAPE(SELF%SD_VV) == SHAPE(RHS%SD_VV))) THEN
    SELF%SD_VV = RHS%SD_VV
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VV different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VN) .AND. ALLOCATED(RHS%SD_VN)) THEN
  IF (ALL(SHAPE(SELF%SD_VN) == SHAPE(RHS%SD_VN))) THEN
    SELF%SD_VN = RHS%SD_VN
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VN different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VH) .AND. ALLOCATED(RHS%SD_VH)) THEN
  IF (ALL(SHAPE(SELF%SD_VH) == SHAPE(RHS%SD_VH))) THEN
    SELF%SD_VH = RHS%SD_VH
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VH different shapes")
  ENDIF
ENDIF
IF (ALLOCATED(SELF%SD_VK) .AND. ALLOCATED(RHS%SD_VK)) THEN
  IF (ALL(SHAPE(SELF%SD_VK) == SHAPE(RHS%SD_VK))) THEN
    SELF%SD_VK = RHS%SD_VK
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VK different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VA) .AND. ALLOCATED(RHS%SD_VA)) THEN
  IF (ALL(SHAPE(SELF%SD_VA) == SHAPE(RHS%SD_VA))) THEN
    SELF%SD_VA = RHS%SD_VA
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VA different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VC) .AND. ALLOCATED(RHS%SD_VC)) THEN
  IF (ALL(SHAPE(SELF%SD_VC) == SHAPE(RHS%SD_VC))) THEN
    SELF%SD_VC = RHS%SD_VC
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VC different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_V2) .AND. ALLOCATED(RHS%SD_V2)) THEN
  IF (ALL(SHAPE(SELF%SD_V2) == SHAPE(RHS%SD_V2))) THEN
    SELF%SD_V2 = RHS%SD_V2
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_V2 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_V3) .AND. ALLOCATED(RHS%SD_V3)) THEN
  IF (ALL(SHAPE(SELF%SD_V3) == SHAPE(RHS%SD_V3))) THEN
    SELF%SD_V3 = RHS%SD_V3
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_V3 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_OC) .AND. ALLOCATED(RHS%SD_OC)) THEN
  IF (ALL(SHAPE(SELF%SD_OC) == SHAPE(RHS%SD_OC))) THEN
    SELF%SD_OC = RHS%SD_OC
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_OC different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VD) .AND. ALLOCATED(RHS%SD_VD)) THEN
  IF (ALL(SHAPE(SELF%SD_VD) == SHAPE(RHS%SD_VD))) THEN
    SELF%SD_VD = RHS%SD_VD
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VD different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_WS) .AND. ALLOCATED(RHS%SD_WS)) THEN
  IF (ALL(SHAPE(SELF%SD_WS) == SHAPE(RHS%SD_WS))) THEN
    SELF%SD_WS = RHS%SD_WS
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_WS different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_WW) .AND. ALLOCATED(RHS%SD_WW)) THEN
  IF (ALL(SHAPE(SELF%SD_WW) == SHAPE(RHS%SD_WW))) THEN
    SELF%SD_WW = RHS%SD_WW
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_WW different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VX) .AND. ALLOCATED(RHS%SD_VX)) THEN
  IF (ALL(SHAPE(SELF%SD_VX) == SHAPE(RHS%SD_VX))) THEN
    SELF%SD_VX = RHS%SD_VX
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_VX different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_XA) .AND. ALLOCATED(RHS%SD_XA)) THEN
  IF (ALL(SHAPE(SELF%SD_XA) == SHAPE(RHS%SD_XA))) THEN
    SELF%SD_XA = RHS%SD_XA
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_XA different shapes")
  ENDIF
ENDIF
IF (ALLOCATED(SELF%SD_DI) .AND. ALLOCATED(RHS%SD_DI)) THEN
  IF (ALL(SHAPE(SELF%SD_DI) == SHAPE(RHS%SD_DI))) THEN
    SELF%SD_DI = RHS%SD_DI
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_DI different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_XP) .AND. ALLOCATED(RHS%SD_XP)) THEN
  IF (ALL(SHAPE(SELF%SD_XP) == SHAPE(RHS%SD_XP))) THEN
    SELF%SD_XP = RHS%SD_XP
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_XP different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_XP2) .AND. ALLOCATED(RHS%SD_XP2)) THEN
  IF (ALL(SHAPE(SELF%SD_XP2) == SHAPE(RHS%SD_XP2))) THEN
    SELF%SD_XP2 = RHS%SD_XP2
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_XP2 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_XR) .AND. ALLOCATED(RHS%SD_XR)) THEN
  IF (ALL(SHAPE(SELF%SD_XR) == SHAPE(RHS%SD_XR))) THEN
    SELF%SD_XR = RHS%SD_XR
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_XR different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_X2) .AND. ALLOCATED(RHS%SD_X2)) THEN
  IF (ALL(SHAPE(SELF%SD_X2) == SHAPE(RHS%SD_X2))) THEN
    SELF%SD_X2 = RHS%SD_X2
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_X2 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_SFL) .AND. ALLOCATED(RHS%SD_SFL)) THEN
  IF (ALL(SHAPE(SELF%SD_SFL) == SHAPE(RHS%SD_SFL))) THEN
    SELF%SD_SFL = RHS%SD_SFL
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_SFL different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_SFO) .AND. ALLOCATED(RHS%SD_SFO)) THEN
  IF (ALL(SHAPE(SELF%SD_SFO) == SHAPE(RHS%SD_SFO))) THEN
    SELF%SD_SFO = RHS%SD_SFO
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_SFO different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_PF) .AND. ALLOCATED(RHS%SD_PF)) THEN
  IF (ALL(SHAPE(SELF%SD_PF) == SHAPE(RHS%SD_PF))) THEN
    SELF%SD_PF = RHS%SD_PF
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_PF different shapes")
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:COPY_SURF',1,ZHOOK_HANDLE)

END SUBROUTINE COPY_SURF

!=========================================================================

SUBROUTINE AXPBY_SURF(SELF,PA,RHS,PB)

IMPLICIT NONE
TYPE(TSURF)    , INTENT(INOUT) :: SELF
REAL(KIND=JPRB), INTENT(IN)    :: PA
TYPE(TSURF)    , INTENT(IN)    :: RHS
REAL(KIND=JPRB), INTENT(IN)    :: PB

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:AXPBY_SURF',0,ZHOOK_HANDLE)

! pronostic groups:
IF (ALLOCATED(SELF%SP_SB) .AND. ALLOCATED(RHS%SP_SB)) THEN
  IF (ALL(SHAPE(SELF%SP_SB) == SHAPE(RHS%SP_SB))) THEN
    SELF%SP_SB = PA*SELF%SP_SB + PB*RHS%SP_SB
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SP_SB different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_SG) .AND. ALLOCATED(RHS%SP_SG)) THEN
  IF (ALL(SHAPE(SELF%SP_SG) == SHAPE(RHS%SP_SG))) THEN
    SELF%SP_SG = PA*SELF%SP_SG + PB*RHS%SP_SG
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SP_SG different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_SL) .AND. ALLOCATED(RHS%SP_SL)) THEN
  IF (ALL(SHAPE(SELF%SP_SL) == SHAPE(RHS%SP_SL))) THEN
    SELF%SP_SL = PA*SELF%SP_SL + PB*RHS%SP_SL
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SP_SL different shapes")
  ENDIF
ENDIF


IF (ALLOCATED(SELF%SP_RR) .AND. ALLOCATED(RHS%SP_RR)) THEN
  IF (ALL(SHAPE(SELF%SP_RR) == SHAPE(RHS%SP_RR))) THEN
    SELF%SP_RR = PA*SELF%SP_RR + PB*RHS%SP_RR
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SP_RR different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_CL) .AND. ALLOCATED(RHS%SP_CL)) THEN
  IF (ALL(SHAPE(SELF%SP_CL) == SHAPE(RHS%SP_CL))) THEN
    SELF%SP_CL = PA*SELF%SP_CL + PB*RHS%SP_CL
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SP_CL different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_OM) .AND. ALLOCATED(RHS%SP_OM)) THEN
  IF (ALL(SHAPE(SELF%SP_OM) == SHAPE(RHS%SP_OM))) THEN
    SELF%SP_OM = PA*SELF%SP_OM + PB*RHS%SP_OM
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SP_OM different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_EP) .AND. ALLOCATED(RHS%SP_EP)) THEN
  IF (ALL(SHAPE(SELF%SP_EP) == SHAPE(RHS%SP_EP))) THEN
    SELF%SP_EP = PA*SELF%SP_EP + PB*RHS%SP_EP
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SP_EP different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_X2) .AND. ALLOCATED(RHS%SP_X2)) THEN
  IF (ALL(SHAPE(SELF%SP_X2) == SHAPE(RHS%SP_X2))) THEN
    SELF%SP_X2 = PA*SELF%SP_X2 + PB*RHS%SP_X2
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SP_X2 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SP_CI) .AND. ALLOCATED(RHS%SP_CI)) THEN
  IF (ALL(SHAPE(SELF%SP_CI) == SHAPE(RHS%SP_CI))) THEN
    SELF%SP_CI = PA*SELF%SP_CI + PB*RHS%SP_CI
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SP_CI different shapes")
  ENDIF
ENDIF

! diagnostic groups:
IF (ALLOCATED(SELF%SD_VF) .AND. ALLOCATED(RHS%SD_VF)) THEN
  IF (ALL(SHAPE(SELF%SD_VF) == SHAPE(RHS%SD_VF))) THEN
    SELF%SD_VF= PA*SELF%SD_VF + PB*RHS%SD_VF
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_VF different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VP) .AND. ALLOCATED(RHS%SD_VP)) THEN
  IF (ALL(SHAPE(SELF%SD_VP) == SHAPE(RHS%SD_VP))) THEN
    SELF%SD_VP= PA*SELF%SD_VP + PB*RHS%SD_VP
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_VP different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VV) .AND. ALLOCATED(RHS%SD_VV)) THEN
  IF (ALL(SHAPE(SELF%SD_VV) == SHAPE(RHS%SD_VV))) THEN
    SELF%SD_VV = PA*SELF%SD_VV + PB*RHS%SD_VV
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_VV different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VN) .AND. ALLOCATED(RHS%SD_VN)) THEN
  IF (ALL(SHAPE(SELF%SD_VN) == SHAPE(RHS%SD_VN))) THEN
    SELF%SD_VN = PA*SELF%SD_VN + PB*RHS%SD_VN
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_VN different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VH) .AND. ALLOCATED(RHS%SD_VH)) THEN
  IF (ALL(SHAPE(SELF%SD_VH) == SHAPE(RHS%SD_VH))) THEN
    SELF%SD_VH = PA*SELF%SD_VH + PB*RHS%SD_VH
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_VH different shapes")
  ENDIF
ENDIF
IF (ALLOCATED(SELF%SD_VK) .AND. ALLOCATED(RHS%SD_VK)) THEN
  IF (ALL(SHAPE(SELF%SD_VK) == SHAPE(RHS%SD_VK))) THEN
    SELF%SD_VK = PA*SELF%SD_VK + PB*RHS%SD_VK
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_VK different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VA) .AND. ALLOCATED(RHS%SD_VA)) THEN
  IF (ALL(SHAPE(SELF%SD_VA) == SHAPE(RHS%SD_VA))) THEN
    SELF%SD_VA = PA*SELF%SD_VA + PB*RHS%SD_VA
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_VA different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VC) .AND. ALLOCATED(RHS%SD_VC)) THEN
  IF (ALL(SHAPE(SELF%SD_VC) == SHAPE(RHS%SD_VC))) THEN
    SELF%SD_VC = PA*SELF%SD_VC + PB*RHS%SD_VC
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_VC different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_V2) .AND. ALLOCATED(RHS%SD_V2)) THEN
  IF (ALL(SHAPE(SELF%SD_V2) == SHAPE(RHS%SD_V2))) THEN
    SELF%SD_V2 = PA*SELF%SD_V2 + PB*RHS%SD_V2
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_V2 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_V3) .AND. ALLOCATED(RHS%SD_V3)) THEN
  IF (ALL(SHAPE(SELF%SD_V3) == SHAPE(RHS%SD_V3))) THEN
    SELF%SD_V3 = PA*SELF%SD_V3 + PB*RHS%SD_V3
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_V3 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_OC) .AND. ALLOCATED(RHS%SD_OC)) THEN
  IF (ALL(SHAPE(SELF%SD_OC) == SHAPE(RHS%SD_OC))) THEN
    SELF%SD_OC = PA*SELF%SD_OC + PB*RHS%SD_OC
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_OC different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VD) .AND. ALLOCATED(RHS%SD_VD)) THEN
  IF (ALL(SHAPE(SELF%SD_VD) == SHAPE(RHS%SD_VD))) THEN
    SELF%SD_VD = PA*SELF%SD_VD + PB*RHS%SD_VD
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_VD different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_WS) .AND. ALLOCATED(RHS%SD_WS)) THEN
  IF (ALL(SHAPE(SELF%SD_WS) == SHAPE(RHS%SD_WS))) THEN
    SELF%SD_WS = PA*SELF%SD_WS + PB*RHS%SD_WS
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_WS different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_WW) .AND. ALLOCATED(RHS%SD_WW)) THEN
  IF (ALL(SHAPE(SELF%SD_WW) == SHAPE(RHS%SD_WW))) THEN
    SELF%SD_WW = PA*SELF%SD_WW + PB*RHS%SD_WW
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_WW different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_VX) .AND. ALLOCATED(RHS%SD_VX)) THEN
  IF (ALL(SHAPE(SELF%SD_VX) == SHAPE(RHS%SD_VX))) THEN
    SELF%SD_VX = PA*SELF%SD_VX + PB*RHS%SD_VX
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_VX different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_XA) .AND. ALLOCATED(RHS%SD_XA)) THEN
  IF (ALL(SHAPE(SELF%SD_XA) == SHAPE(RHS%SD_XA))) THEN
    SELF%SD_XA = PA*SELF%SD_XA + PB*RHS%SD_XA
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_XA different shapes")
  ENDIF
ENDIF
IF (ALLOCATED(SELF%SD_DI) .AND. ALLOCATED(RHS%SD_DI)) THEN
  IF (ALL(SHAPE(SELF%SD_DI) == SHAPE(RHS%SD_DI))) THEN
    SELF%SD_DI = PA*SELF%SD_DI + PB*RHS%SD_DI
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_DI different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_XP) .AND. ALLOCATED(RHS%SD_XP)) THEN
  IF (ALL(SHAPE(SELF%SD_XP) == SHAPE(RHS%SD_XP))) THEN
    SELF%SD_XP = PA*SELF%SD_XP + PB*RHS%SD_XP
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_XP different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_XP2) .AND. ALLOCATED(RHS%SD_XP2)) THEN
  IF (ALL(SHAPE(SELF%SD_XP2) == SHAPE(RHS%SD_XP2))) THEN
    SELF%SD_XP2 = PA*SELF%SD_XP2 + PB*RHS%SD_XP2
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:COPY SD_XP2 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_XR) .AND. ALLOCATED(RHS%SD_XR)) THEN
  IF (ALL(SHAPE(SELF%SD_XR) == SHAPE(RHS%SD_XR))) THEN
    SELF%SD_XR = PA*SELF%SD_XR + PB*RHS%SD_XR
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_XR different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_X2) .AND. ALLOCATED(RHS%SD_X2)) THEN
  IF (ALL(SHAPE(SELF%SD_X2) == SHAPE(RHS%SD_X2))) THEN
    SELF%SD_X2 = PA*SELF%SD_X2 + PB*RHS%SD_X2
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_X2 different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_SFL) .AND. ALLOCATED(RHS%SD_SFL)) THEN
  IF (ALL(SHAPE(SELF%SD_SFL) == SHAPE(RHS%SD_SFL))) THEN
    SELF%SD_SFL = PA*SELF%SD_SFL + PB*RHS%SD_SFL
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_SFL different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_SFO) .AND. ALLOCATED(RHS%SD_SFO)) THEN
  IF (ALL(SHAPE(SELF%SD_SFO) == SHAPE(RHS%SD_SFO))) THEN
    SELF%SD_SFO = PA*SELF%SD_SFO + PB*RHS%SD_SFO
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_SFO different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%SD_PF) .AND. ALLOCATED(RHS%SD_PF)) THEN
  IF (ALL(SHAPE(SELF%SD_PF) == SHAPE(RHS%SD_PF))) THEN
    SELF%SD_PF = PA*SELF%SD_PF + PB*RHS%SD_PF
  ELSE
    CALL ABOR1("SURFACE_FIELDS_MIX:ADD SD_PF different shapes")
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:AXPBY_SURF',1,ZHOOK_HANDLE)

END SUBROUTINE AXPBY_SURF
!=========================================================================

SUBROUTINE ADD_SURF(SELF,RHS)

IMPLICIT NONE
TYPE(TSURF), INTENT(INOUT) :: SELF
TYPE(TSURF), INTENT(IN)    :: RHS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:ADD_SURF',0,ZHOOK_HANDLE)

CALL AXPBY_SURF(SELF,1.0_JPRB,RHS,1.0_JPRB)

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:ADD_SURF',1,ZHOOK_HANDLE)

END SUBROUTINE ADD_SURF

!=========================================================================

SUBROUTINE DIFF_SURF(SELF,RHS,DIFF)

IMPLICIT NONE
TYPE(TSURF), INTENT(IN) :: SELF
TYPE(TSURF), INTENT(IN)    :: RHS
TYPE(TSURF), INTENT(INOUT) :: DIFF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:DIFF_SURF',0,ZHOOK_HANDLE)

DIFF%SP_SB  = SELF%SP_SB  - RHS%SP_SB
DIFF%SP_SG  = SELF%SP_SG  - RHS%SP_SG
DIFF%SP_SL  = SELF%SP_SL  - RHS%SP_SL
DIFF%SP_RR  = SELF%SP_RR  - RHS%SP_RR
DIFF%SP_CL  = SELF%SP_CL  - RHS%SP_CL
DIFF%SP_OM  = SELF%SP_OM  - RHS%SP_OM
DIFF%SP_EP  = SELF%SP_EP  - RHS%SP_EP
DIFF%SP_X2  = SELF%SP_X2  - RHS%SP_X2
DIFF%SP_CI  = SELF%SP_CI  - RHS%SP_CI
DIFF%SD_VF  = SELF%SD_VF  - RHS%SD_VF
DIFF%SD_VP  = SELF%SD_VP  - RHS%SD_VP
DIFF%SD_VV  = SELF%SD_VV  - RHS%SD_VV
DIFF%SD_VN  = SELF%SD_VN  - RHS%SD_VN
DIFF%SD_VH  = SELF%SD_VH  - RHS%SD_VH
DIFF%SD_VK  = SELF%SD_VK  - RHS%SD_VK
DIFF%SD_VA  = SELF%SD_VA  - RHS%SD_VA
DIFF%SD_VC  = SELF%SD_VC  - RHS%SD_VC
DIFF%SD_V2  = SELF%SD_V2  - RHS%SD_V2
DIFF%SD_V3  = SELF%SD_V3  - RHS%SD_V3
DIFF%SD_OC  = SELF%SD_OC  - RHS%SD_OC
DIFF%SD_VD  = SELF%SD_VD  - RHS%SD_VD
DIFF%SD_WS  = SELF%SD_WS  - RHS%SD_WS
DIFF%SD_WW  = SELF%SD_WW  - RHS%SD_WW
DIFF%SD_VX  = SELF%SD_VX  - RHS%SD_VX
DIFF%SD_XA  = SELF%SD_XA  - RHS%SD_XA
DIFF%SD_DI  = SELF%SD_DI  - RHS%SD_DI
DIFF%SD_XP  = SELF%SD_XP  - RHS%SD_XP
DIFF%SD_XP2 = SELF%SD_XP2 - RHS%SD_XP2
DIFF%SD_XR  = SELF%SD_XR  - RHS%SD_XR
DIFF%SD_X2  = SELF%SD_X2  - RHS%SD_X2
DIFF%SD_SFL = SELF%SD_SFL - RHS%SD_SFL
DIFF%SD_SFO = SELF%SD_SFO - RHS%SD_SFO
DIFF%SD_PF  = SELF%SD_PF  - RHS%SD_PF

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:DIFF_SURF',1,ZHOOK_HANDLE)

END SUBROUTINE DIFF_SURF

!=========================================================================

SUBROUTINE CLONE_SURF(YDDIM,YDSURF,YDCLONE)
! Allocate and copy a surface field structure from a reference surface field structure
TYPE(TDIM) , INTENT(IN) :: YDDIM
TYPE(TSURF), INTENT(IN) :: YDSURF
TYPE(TSURF), INTENT(OUT) :: YDCLONE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:CLONE_SURF',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, NGPBLKS=>YDDIM%NGPBLKS)

! pronostic groups:
ALLOCATE(YDCLONE%SP_SB(NPROMA,YDSURF%YSP_SBD%NLEVS,YDSURF%YSP_SBD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SP_SG(NPROMA,YDSURF%YSP_SGD%NLEVS,YDSURF%YSP_SGD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SP_SL(NPROMA,YDSURF%YSP_SLD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SP_RR(NPROMA,YDSURF%YSP_RRD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SP_CL(NPROMA,YDSURF%YSP_CLD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SP_OM(NPROMA,YDSURF%YSP_OMD%NLEVS,YDSURF%YSP_OMD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SP_EP(NPROMA,YDSURF%YSP_EPD%NLEVS,YDSURF%YSP_EPD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SP_X2(NPROMA,YDSURF%YSP_X2D%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SP_CI(NPROMA,YDSURF%YSP_CID%NDIM,NGPBLKS))
! diagnostic groups:
ALLOCATE(YDCLONE%SD_VF(NPROMA,YDSURF%YSD_VFD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_VP(NPROMA,YDSURF%YSD_VPD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_VV(NPROMA,YDSURF%YSD_VVD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_VN(NPROMA,YDSURF%YSD_VND%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_VH(NPROMA,YDSURF%YSD_VHD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_VK(NPROMA,YDSURF%YSD_VKD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_VA(NPROMA,YDSURF%YSD_VAD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_VC(NPROMA,YDSURF%YSD_VCD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_V2(NPROMA,YDSURF%YSD_V2D%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_V3(NPROMA,YDSURF%YSD_V3D%NLEVS,YDSURF%YSD_V3D%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_VD(NPROMA,YDSURF%YSD_VDD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_SM(NPROMA,YDSURF%YSD_SMD%NLEVS,YDSURF%YSD_SMD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_WS(NPROMA,YDSURF%YSD_WSD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_WW(NPROMA,YDSURF%YSD_WWD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_VX(NPROMA,YDSURF%YSD_VXD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_XA(NPROMA,YDSURF%YSD_XAD%NLEVS,YDSURF%YSD_XAD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_XR(NPROMA,YDSURF%YSD_XRD%NLEVS,YDSURF%YSD_XRD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_X2(NPROMA,YDSURF%YSD_X2D%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_SFL(NPROMA,YDSURF%YSD_SFLD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_SFO(NPROMA,YDSURF%YSD_SFOD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_PF(NPROMA,YDSURF%YSD_PFD%NLEVS,YDSURF%YSD_PFD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_DI(NPROMA,YDSURF%YSD_DID%NLEVS,YDSURF%YSD_DID%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_OC(NPROMA,YDSURF%YSD_OCD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_XP(NPROMA,YDSURF%YSD_XPD%NLEVS,YDSURF%YSD_XPD%NDIM,NGPBLKS))
ALLOCATE(YDCLONE%SD_XP2(NPROMA,YDSURF%YSD_XP2D%NLEVS,YDSURF%YSD_XP2D%NDIM,NGPBLKS))

YDCLONE=YDSURF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:CLONE_SURF',1,ZHOOK_HANDLE)

END SUBROUTINE CLONE_SURF
!-------------------------------------------------------------------------

SUBROUTINE SURF_NORM(YDGEOMETRY,YDSURF)
USE GEOMETRY_MOD       , ONLY : GEOMETRY

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF), INTENT(IN) :: YDSURF
INTEGER(KIND=JPIM) :: JV, ICLI ,J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "gpnorm2.intfb.h"
#include "gpnorm3.intfb.h"

IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SURF_NORM',0,ZHOOK_HANDLE)
ASSOCIATE( &
 & SD_VA=>YDSURF%SD_VA, SD_VC=>YDSURF%SD_VC, SD_VF=>YDSURF%SD_VF, &
 & SD_VV=>YDSURF%SD_VV, SD_WS=>YDSURF%SD_WS, SD_WW=>YDSURF%SD_WW, &
 & SP_RR=>YDSURF%SP_RR, SP_SB=>YDSURF%SP_SB, SP_SG=>YDSURF%SP_SG, &
 & SP_SL=>YDSURF%SP_SL, SP_CI=>YDSURF%SP_CI, SP_CL=>YDSURF%SP_CL, &
 & SP_X2=>YDSURF%SP_X2, SD_VX=>YDSURF%SD_VX, &
 & YSD_VAD=>YDSURF%YSD_VAD, YSD_VCD=>YDSURF%YSD_VCD, &
 & YSD_VFD=>YDSURF%YSD_VFD, YSD_VVD=>YDSURF%YSD_VVD, YSD_VXD=>YDSURF%YSD_VXD, &
 & YSD_WSD=>YDSURF%YSD_WSD, YSD_WWD=>YDSURF%YSD_WWD, YSP_RRD=>YDSURF%YSP_RRD, &
 & YSP_SBD=>YDSURF%YSP_SBD, YSP_SGD=>YDSURF%YSP_SGD, YSP_SLD=>YDSURF%YSP_SLD, &
 & YSP_CID=>YDSURF%YSP_CID, YSP_CLD=>YDSURF%YSP_CLD, YSP_X2D=>YDSURF%YSP_X2D )

WRITE(NULOUT,'(A)') ' SURF_NORM: STATISTICS FOR ALL SURFACE FIELDS'
IF (YSP_SBD%NUMFLDS>0) THEN
  WRITE(NULOUT,'(A,I1,A)') ' SOILB   ',YSP_SBD%NUMFLDS,' FIELDS '
  DO JV=1,YSP_SBD%NUMFLDS
    CALL GPNORM3(YDGEOMETRY,SP_SB,YSP_SBD%NUMFLDS,JV,KLEVS=YSP_SBD%NLEVS,LDLEVELS=.TRUE.)
  ENDDO
ENDIF
IF (YSP_SGD%NUMFLDS>0) THEN
  DO JV=1,YSP_SGD%NUMFLDS
    WRITE(NULOUT,'(A,I4)') ' SNOWG, FIELD ',JV
    CALL GPNORM3(YDGEOMETRY,SP_SG,YSP_SGD%NUMFLDS,JV,KLEVS=YSP_SGD%NLEVS,LDLEVELS=.TRUE.)
  ENDDO
ENDIF
IF (YSP_CID%NUMFLDS>0) THEN
  WRITE(NULOUT,'(A)') ' CANRI '
  CALL GPNORM2(YDGEOMETRY,YSP_CID%NUMFLDS,1,SP_CI)
ENDIF
IF (YSP_CLD%NUMFLDS>0) THEN
  WRITE(NULOUT,'(A)') ' CLS '
  CALL GPNORM2(YDGEOMETRY,YSP_CLD%NUMFLDS,1,SP_CL)
ENDIF
IF (YSP_X2D%NUMFLDS>0) THEN
  WRITE(NULOUT,'(A)') ' XTRP2 '
  CALL GPNORM2(YDGEOMETRY,YSP_X2D%NUMFLDS,1,SP_X2)
ENDIF
IF (YSP_SLD%NUMFLDS>0) THEN
  WRITE(NULOUT,'(A)') ' LAKEB '
  CALL GPNORM2(YDGEOMETRY,YSP_SLD%NUMFLDS,1,SP_SL)
ENDIF
IF (YSP_RRD%NUMFLDS>0) THEN
  WRITE(NULOUT,'(A)') ' RESVR '
  CALL GPNORM2(YDGEOMETRY,YSP_RRD%NUMFLDS,1,SP_RR)
ENDIF
IF (YSD_WSD%NUMFLDS>0) THEN
  WRITE(NULOUT,'(A)') ' FROM WAM '
  CALL GPNORM2(YDGEOMETRY,YSD_WSD%NUMFLDS,1,SD_WS)
ENDIF
IF (YSD_WWD%NUMFLDS>0) THEN
  WRITE(NULOUT,'(A)') ' TO WAM '
  CALL GPNORM2(YDGEOMETRY,YSD_WWD%NUMFLDS,1,SD_WW)
ENDIF


IF (YSD_VFD%NUMFLDS > 0) THEN
  WRITE(NULOUT,'(A)') ' VARSF'
  CALL GPNORM2(YDGEOMETRY,YSD_VFD%NUMFLDS,1,SD_VF)
ENDIF
IF (YSD_VVD%NUMFLDS > 0) THEN
  WRITE(NULOUT,'(A)') ' VCLIV'
  CALL GPNORM2(YDGEOMETRY,YSD_VVD%NUMFLDS,1,SD_VV)
ENDIF
IF (YSD_VAD%NUMFLDS > 0) THEN
  WRITE(NULOUT,'(A)') ' VCLIA'
  CALL GPNORM2(YDGEOMETRY,YSD_VAD%NUMFLDS,1,SD_VA)
ENDIF
IF (YSD_VCD%NUMFLDS > 0) THEN
  WRITE(NULOUT,'(A)') ' VO3ABC'
  CALL GPNORM2(YDGEOMETRY,YSD_VCD%NUMFLDS,1,SD_VC)
ENDIF
IF (YSD_VXD%NUMFLDS > 0) THEN
  WRITE(NULOUT,'(A)') ' VCLIX'
  CALL GPNORM2(YDGEOMETRY,YSD_VXD%NUMFLDS,1,SD_VX)
ENDIF
CALL FLUSH(NULOUT)

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFACE_FIELDS_MIX:SURF_NORM',1,ZHOOK_HANDLE)
END SUBROUTINE SURF_NORM
END MODULE SURFACE_FIELDS_MIX
