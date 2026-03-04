! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      SUBROUTINE SUD1C_NC(YDDIMV,YDMODEL,YDSURF)

#ifdef DOC

!**** *SUD1C_NC * - Routine to initialize diagnostic NetCDF output

!     Purpose.
!     --------
!           Initialize NetCDF outputfile diagvar.nc

!**   Interface.
!     ----------
!        *CALL* *SUD1C_NC

!     Explicit arguments :
!     --------------------
!        none

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!        Bart vd HURK (KNMI)

!     Modifications.
!     --------------
!        Original     2000-7-14
!        M. Ko"hler      9-2000  adoption to single column model
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 
!        G. Carver     Aug/2012  Bugfixes for gfortran compiler
!        F. Vana       Sep/2014  New fields for convection
!     ------------------------------------------------------------------
#endif

USE YOMDIMV  , ONLY : TDIMV
USE TYPE_MODEL, ONLY : MODEL
USE PARKIND1 , ONLY : JPIM
USE YOMGPD1C , ONLY : VEXTR2   ,VEXTRA
USE YOMLOG1C , ONLY : CMODID   ,CSIMID   ,NPOSDIA  ,NPOSDIA2
USE YOMRIP0  , ONLY : NINDAT   ,NSSSSS
USE SURFACE_FIELDS_MIX, ONLY : TSURF

IMPLICIT NONE

TYPE(TDIMV), INTENT(INOUT) :: YDDIMV
TYPE(MODEL), INTENT(INOUT) :: YDMODEL
TYPE(TSURF), INTENT(INOUT) :: YDSURF

INTEGER(KIND=JPIM) :: ICSS
INTEGER(KIND=JPIM) :: INCID(2), N, I,&
 & NLEVDID, NLEVP1DID, NLEVSDID, NTILESDID, NORGDID, NTIMDID, NVARID,&
 & IDIMID2(2), IDIMID3(2), IDIMID4(2), IDIMID5(2), IDIMID6(2),&
 & ILEV(YDDIMV%NFLEVG), ILEVP1(YDDIMV%NFLEVG+1), ILEVS(YDSURF%YSP_SBD%NLEVS),&
 & ISTATUS, IACCUR,NCEXTRID,IDIMIDE(2)

CHARACTER (LEN=40) :: TITLE
CHARACTER (LEN=2)  :: NTXT

#include "netcdf.inc"

!     ------------------------------------------------------------------
#include "handle_err_nc.intfb.h"
#include "varsetup1c_nc.intfb.h"
!     ------------------------------------------------------------------
ASSOCIATE(YDPHYDS=>YDMODEL%YRML_PHY_MF%YRPHYDS,YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY)

ICSS = YDSURF%YSP_SBD%NLEVS

!        1.    open NetCDF files.
!              ------------------

ISTATUS = NF_CREATE ('diagvar.nc',  NF_CLOBBER, INCID(1))
CALL HANDLE_ERR_NC(ISTATUS)
WRITE(*,*) 'NETCDF-FILE diagvar.nc OPENED ON UNIT ', INCID(1)
NPOSDIA = INCID(1)       !store unit number for later use

ISTATUS = NF_CREATE ('diagvar2.nc', NF_CLOBBER, INCID(2))
CALL HANDLE_ERR_NC(ISTATUS)
WRITE(*,*) 'NETCDF-FILE diagvar2.nc OPENED ON UNIT ',INCID(2)
NPOSDIA2= INCID(2)       !store unit number for later use

!        output accuracy
 IACCUR=NF_FLOAT
!iaccur=NF_DOUBLE


!        2.    meta data set up.
!              -----------------

DO N=1,2

!        title
TITLE = 'SCM: ' // TRIM(CMODID) // '  Sim: ' // TRIM(CSIMID)
ISTATUS = NF_PUT_ATT_TEXT (INCID(N), NF_GLOBAL, 'title', LEN(TRIM(TITLE)), TRIM(TITLE))
CALL HANDLE_ERR_NC(ISTATUS)

!        model identification
  ISTATUS = NF_PUT_ATT_TEXT (INCID(N), NF_GLOBAL, 'modelID', LEN(TRIM(CMODID)), TRIM(CMODID))
  CALL HANDLE_ERR_NC(ISTATUS)

!        simulation identification
  ISTATUS = NF_PUT_ATT_TEXT (INCID(N), NF_GLOBAL, 'simulationID', LEN(TRIM(CSIMID)), TRIM(CSIMID))
  CALL HANDLE_ERR_NC(ISTATUS)

!        dataID for MetView
  ISTATUS = NF_PUT_ATT_TEXT (INCID(N), NF_GLOBAL, 'dataID', 10, 'SCM_OUTPUT' ) 
  CALL HANDLE_ERR_NC(ISTATUS)

!        model startup
  ISTATUS = NF_PUT_ATT_INT (INCID(N), NF_GLOBAL,  'start_day',  NF_INT, 1, NINDAT)
  CALL HANDLE_ERR_NC(ISTATUS)
  ISTATUS = NF_PUT_ATT_INT (INCID(N), NF_GLOBAL,  'start_hour', NF_INT, 1, NSSSSS)
  CALL HANDLE_ERR_NC(ISTATUS)

!        create dimensions
  ISTATUS = NF_DEF_DIM (INCID(N), 'nlev',   YDDIMV%NFLEVG,       NLEVDID)   !atmosphere
  CALL HANDLE_ERR_NC(ISTATUS)
  ISTATUS = NF_DEF_DIM (INCID(N), 'nlevp1', YDDIMV%NFLEVG+1,     NLEVP1DID) !atmosphere + 1
  CALL HANDLE_ERR_NC(ISTATUS)
  ISTATUS = NF_DEF_DIM (INCID(N), 'nlevs',  ICSS,         NLEVSDID)  !land/sea-ice
  CALL HANDLE_ERR_NC(ISTATUS)
  ISTATUS = NF_DEF_DIM (INCID(N), 'ntiles', YDDPHY%NTILES,       NTILESDID) !tiles
  CALL HANDLE_ERR_NC(ISTATUS)
  ISTATUS = NF_DEF_DIM (INCID(N), 'norg',   4,            NORGDID)   !orog. var.
  CALL HANDLE_ERR_NC(ISTATUS)
  ISTATUS = NF_DEF_DIM (INCID(N), 'time',   NF_UNLIMITED, NTIMDID)   !time
  CALL HANDLE_ERR_NC(ISTATUS)
  ISTATUS = NF_DEF_DIM (INCID(N), 'ncextr',   YDDPHY%NCEXTR,     NCEXTRID)  !extra variable levels !Maike
  CALL HANDLE_ERR_NC(ISTATUS)
! istatus = NF_DEF_DIM (incid(n), 'timestp',NF_UNLIMITED, ntimdid)   !time steps
! call handle_err_nc(istatus)


!        2.1   create dimensional variables.

!        atmospheric model levels
  CALL VARSETUP1C_NC (INCID(N), NF_INT, 1, (/ NLEVDID /), &
       'nlev',       'Atmospheric Model Levels',            'count')
!        atmospheric model levels plus 1
  CALL VARSETUP1C_NC (INCID(N), NF_INT, 1, (/ NLEVP1DID /), &
       'nlevp1',     'Atmospheric Model Levels',            'count')
!        soil/sea-ice model levels
  CALL VARSETUP1C_NC (INCID(N), NF_INT, 1, (/ NLEVSDID /), &
       'nlevs',      'Soil/Sea-Ice Model Levels',           'count')
!        tiles
  CALL VARSETUP1C_NC (INCID(N), NF_INT, 1, (/ NTILESDID /), &
       'ntiles',     'Tiles',                               'count')
!        orographic variances
  CALL VARSETUP1C_NC (INCID(N), NF_INT, 1, (/ NORGDID /), &
       'norg',       'Orographic Variances',                'count')
!        time
  CALL VARSETUP1C_NC (INCID(N), IACCUR, 1, (/ NTIMDID /), &
       'time',       'Time',                                'seconds')
!        extra variable levels
  CALL VARSETUP1C_NC (INCID(N), NF_INT, 1, (/NCEXTRID/), &
       'ncextr',     'Extra Variable Levels',               'count')
!        timestep
! call varsetup1c_nc (incid(n), NF_INT, 1, (/ ntimdid /), &
!      'timestp',    'Model Time Step',                     'count')

ENDDO

IDIMID2(1) = NLEVDID
IDIMID3(1) = NLEVP1DID
IDIMID5(1) = NTILESDID
IDIMID6(1) = NORGDID
IDIMID2(2) = NTIMDID
IDIMID3(2) = NTIMDID
IDIMID5(2) = NTIMDID
IDIMID6(2) = NTIMDID
IDIMIDE(1) = NCEXTRID 
IDIMIDE(2) = NTIMDID  

!        2.2   create coordinate.

CALL VARSETUP1C_NC (INCID(1), IACCUR, 2, IDIMID2, &
   'pressure_f',     'Pressure - full levels',              'Pa')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 2, IDIMID3, &
   'pressure_h',     'Pressure - half levels',              'Pa')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 2, IDIMID2, &
   'height_f',       'Height - full levels',                'm')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 2, IDIMID3, &
   'height_h',       'Height - half levels',                'm')

!        2.3   create scalar diagnostic variables.

CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'wat_vap_path',   'Water Vapor Path',                    'kg/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'liq_wat_path',   'Liquid Water Path',                   'kg/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'ice_wat_path',   'Ice Water Path',                      'kg/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'pbl_height',     'Boundary Layer Height',               'm')

CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'total_cloud',    'Total Cloud Cover',                   '0-1')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'low_cloud',      'Low Cloud Cover',                     '0-1')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'middle_cloud',   'Middle Cloud Cover',                  '0-1')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'high_cloud',     'High Cloud Cover',                    '0-1')

CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'conv_rain',      'Convective Rain',                     'kg/m^2/s')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_conv_rain',  'Accu. Convective Rain',               'kg/m^2/s')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'conv_snow',      'Convective Snow',                     'kg/m^2/s')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_conv_snow',  'Accu. Convective Snow',               'kg/m^2/s')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'stra_rain',      'Stratiform Rain',                     'kg/m^2/s')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_stra_rain',  'Accu. Stratiform Rain',               'kg/m^2/s')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'stra_snow',      'Stratiform Snow',                     'kg/m^2/s')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_stra_snow',  'Accu. Stratiform Snow',               'kg/m^2/s')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'top_swrad_inc',  'Top SW Radiation Incoming',           'W/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'conv_type',      'Deep/Shallow/Midlevel Convection',    '0/1/2/3')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'pbl_type',      'Dry stable/Dry conv/Sc/Cu',    '0/1/2/3')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'n_conv_base',    'Index of Convective Cloud Base',      '1')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'n_conv_top',     'Index of Convective Cloud Top',       '1')

CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'top_swrad',      'Top SW Radiation',                    'W/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'top_swrad_clr',  'Top SW Radiation Clear-Sky',          'W/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'top_lwrad',      'Top LW Radiation',                    'W/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'top_lwrad_clr',  'Top LW Radiation Clear-Sky',          'W/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'sfc_swrad',      'Surface SW Radiation',                'W/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'sfc_swrad_clr',  'Surface SW Radiation Clear-Sky',      'W/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'sfc_lwrad',      'Surface LW Radiation',                'W/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'sfc_lwrad_clr',  'Surface LW Radiation Clear-Sky',      'W/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'sfc_swrad_down', 'Surface Downward SW Radiation',       'W/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'sfc_albedo',     'Surface SW Albedo',                   '0-1')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'sfc_lwrad_down', 'Surface Downward LW Radiation',       'W/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'sfc_emiss',      'Surface Emissivity',                  '0-1')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_top_swrad',  'Accu. Top SW Radiation',              'Ws/m^2')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_sfc_swrad',  'Accu. Surface SW Radiation',          'Ws/m^2')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_top_lwrad',  'Accu. Top LW Radiation',              'Ws/m^2')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_sfc_lwrad',  'Accu. Surface LW Radiation',          'Ws/m^2')

CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'sfc_sen_flx',    'Surface Sensible Heat Flux',          'W/m^2')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_sfc_sen_flx','Accu. Surface Sensible Heat Flux',    'W/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'sfc_lat_flx',    'Surface Latent Heat Flux',            'W/m^2')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_sfc_lat_flx','Accu. Surface Latent Heat Flux',      'W/m^2')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'eva_flx',        'Evaporation Flux',                    'kg/m^2/s')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_eva_flx',    'Accu. Evaporation Flux',              'kg/m^2/s')

CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'u_sfc_strss',    'U-Surface Stress', '')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_u_sfc_strss','Accu. U-Surface Stress', '')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'v_sfc_strss',    'V-Surface Stress',                    '')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_v_sfc_strss','Accu. Surface Stress',                '')

CALL VARSETUP1C_NC (INCID(1), IACCUR, 2, IDIMID2, &
   'sfc_sen_flx_ti', 'Tiled Surface Sensible Heat Flux',    'W/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 2, IDIMID2, &
   'sfc_lat_flx_ti', 'Tiled Surface Latent Heat Flux',      'W/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 2, IDIMID2, &
   'u_sfc_strss_ti', 'Tiled U Surface Stress',              'N/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 2, IDIMID2, &
   'v_sfc_strss_ti', 'Tiled V Surface Stress',              'N/m^2')

CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'sst',            'Sea Surface Temperature',             'K')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 2, IDIMID2, &
   't_skin_ti',      'Tiled Skin Temperature',              'K')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'land_sea_mask',  'Land-Sea Mask',                       '0-1')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'sea_ice_cov',    'Sea Ice Fraction',                    '0-1')

CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'rough_len_mom',  'Surface Roughness length',            'm')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'rough_len_heat', 'Surface Roughness Length for Heat',   'm')

CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'low_veg_cov',    'Low Vegetation Cover',                '0-1')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'high_veg_cov',   'High Vegetation Cover',               '0-1')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'low_veg_type',   'Low Vegetation Type',                 '1-20')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'high_veg_type',  'High Vegetation Type',                '1-20')

CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'turb_diss',      'Turbulent Dissipation',               '')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_turb_diss',  'Accu. Turbulent Dissipation',         '')

CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'gwd_diss',       'Gravity Wave Drag Dissipation',       '')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_gwd_diss',   'Accu. Gravity Wave Drag Dissipation', '')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'gwd_u_strs',     'Gravity Wave Drag Surface U-Stress',  '')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_gwd_u_strs', 'Accu. Gravity Wave Drag Surface U-Stress', '')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'gwd_v_strs',     'Gravity Wave Drag Surface V-Stress',  '')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_gwd_v_strs', 'Accu. Gravity Wave Drag Surface V-Stress', '')

CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'u_wind_10m',     'U-Wind 10m',                         'm/s')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'v_wind_10m',     'V-Wind 10m',                         'm/s')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'temperature_2m', 'Temperature 2m',                     'K')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'dew_point_2m',   'Dew Point Temperature 2m',           'K')
!call varsetup1c_nc (incid(1), iaccur, 1, (/ ntimdid /), &
!   '10m_u_wind',     '10 m U-Wind',                         'm/s')
!call varsetup1c_nc (incid(1), iaccur, 1, (/ ntimdid /), &
!   '10m_v_wind',     '10 m V-Wind',                         'm/s')
!call varsetup1c_nc (incid(1), iaccur, 1, (/ ntimdid /), &
!   '2m_temperature', '2 m Temperature',                     'K')
!call varsetup1c_nc (incid(1), iaccur, 1, (/ ntimdid /), &
!   '2m_dew_point',   '2 m Dew Point Temperature',           'K')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'wind_gust',      'Wind Gust at 10m',                    'm/s')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'sun_duration',   'Sunshine Duration',                   '0-1')
!call varsetup1c_nc (incid(1), iaccur, 1, (/ ntimdid /), &
!   'ls_prec_fct',    'Large-scale Precipitation Fraction',  '0-1')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'runoff',         'Run-off',                             '')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'acc_runoff',     'Accu. Run-off',                       '')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'snow_melt',      'Snow Melt',                           '')
CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
   'snow_evap',      'Snow Evaporation',                    '')

CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'sfc_sen_flx_ext','External Sensible Heat Flux',         'W/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 1, (/ NTIMDID /), &
   'sfc_lat_flx_ext','External Latent Heat Flux',           'W/m^2')
CALL VARSETUP1C_NC (INCID(1), IACCUR, 2, IDIMID6, &
   'dir_orog_var',   'Directional Orographic Variances',    '')


!        2.4   create tendencies.

CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'ls_prec_fct',    'Large-scale Precipitation Fraction',  '0-1') !Maike
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_u_wind',    'Tendency of U-Wind',                  'm/s^2')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_u_wind_d',  'Tendency of U-Wind (dyn)',            'm/s^2')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_u_wind_p',  'Tendency of U-Wind (phy)',            'm/s^2')

CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_v_wind',    'Tendency of V-Wind',                  'm/s^2')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_v_wind_d',  'Tendency of V-Wind (dyn)',            'm/s^2')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_v_wind_p',  'Tendency of V-Wind (phy)',            'm/s^2')

CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_temp',      'Tendency of Temperature',             'K/s')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_temp_d',    'Tendency of Temperature (dyn)',       'K/s')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_temp_p',    'Tendency of Temperature (phy)',       'K/s')

CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_wat_vap',   'Tendency of Water Vapor',             'kg/kg/s')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_wat_vap_d', 'Tendency of Water Vapor (dyn)',       'kg/kg/s')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_wat_vap_p', 'Tendency of Water Vapor (phy)',       'kg/kg/s')

CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_cld_fract', 'Tendency of Cloud Fraction',          'kg/kg/s')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_cld_fract_d','Tendency of Cloud Fraction (dyn)',   'kg/kg/s')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_cld_fract_p', 'Tendency of Cloud Fraction (phy)',  'kg/kg/s')

CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_cld_liq',   'Tendency of Cloud Liq. Water',        'kg/kg/s')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_cld_liq_d', 'Tendency of Cloud Liq. Water (dyn)',  'kg/kg/s')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_cld_liq_p', 'Tendency of Cloud Liq. Water (phy)',  'kg/kg/s')

CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_cld_ice',   'Tendency of Cloud Ice Water',         'kg/kg/s')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_cld_ice_d', 'Tendency of Cloud Ice Water (dyn)',   'kg/kg/s')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'tend_cld_ice_p', 'Tendency of Cloud Ice Water (phy)',   'kg/kg/s')


!        2.5   create column diagnostic variables.

CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'stra_rain_3d',   'Stratiform Rain 3D',                  'kg/m-2/s')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'stra_snow_3d',   'Stratiform Snow 3D',                  'kg/m-2/s')
CALL VARSETUP1C_NC (INCID(2),  IACCUR, 2, IDIMID3, &
   's_flux_str_rain', 's-Flux of Stratiform Rain',          '')
CALL VARSETUP1C_NC (INCID(2),  IACCUR, 2, IDIMID3, &
   's_flux_str_snow', 's-Flux of Stratiform Snow',          '')
CALL VARSETUP1C_NC (INCID(2),  IACCUR, 2, IDIMID3, &
   'ice_flx_str_cond', 'Stratiform Condensation Flux of Ice Water', '')
CALL VARSETUP1C_NC (INCID(2),  IACCUR, 2, IDIMID3, &
   'liq_flx_str_cond', 'Stratiform Condensation Flux of Liq. Water', '')

CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'conv_rain_3d',   'Convective Rain 3D',                  'kg/m-2/s')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'conv_snow_3d',   'Convective Snow 3D',                  'kg/m-2/s')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'conv_flx_wv',    'Conv. Flux of WV',                    'kg/m-2/s')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'conv_flx_s',     'Conv. Flux of Dry Static Energy',     'W/m^2')
CALL VARSETUP1C_NC (INCID(2),  IACCUR, 2, IDIMID3, &
   's_flux_cnv_rain', 's-Flux of Convective Rain',          '')
CALL VARSETUP1C_NC (INCID(2),  IACCUR, 2, IDIMID3, &
   's_flux_cnv_snow', 's-Flux of Convective Snow',          '')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'conv_flx_u',     'Convective U-Momentum Flux',          '')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'conv_flx_v',     'Convective V-Momentum Flux',          '')
CALL VARSETUP1C_NC (INCID(2),  IACCUR, 2, IDIMID3, &
   'ice_flx_cnv_cond', 'Convective Condensation Flux of Ice Water', '')
CALL VARSETUP1C_NC (INCID(2),  IACCUR, 2, IDIMID3, &
   'liq_flx_cnv_cond', 'Convective Condensation Flux of Liq. Water', '')

! New variables for convection
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'cond_detrained_updr',   'Detrained liquid water',         'kg/kg')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID2, &
   'cond_updraft',   'Liquid water content in updrafts',      'kg/kg')

CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'turb_flx_wv',    'Turb. Flux (inc. neg q) of WV',       '')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'turb_flx_s',     'Turb. Flux of Dry Static Energy',     'W/m^2')
CALL VARSETUP1C_NC (INCID(2),  IACCUR, 2, IDIMID3, &
   'turb_flx_wv_qpos', 'Pseudo-Flux of WV to correct q<0',  '')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'turb_flx_u',     'Turbulent Flux of U-Momentum',        '')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'turb_flx_v',     'Turbulent Flux of V-Momentum',        '')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'turb_flx_ice',   'Turbulent Flux (inc. neg. q) of Ice Water', '')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'turb_flx_liq',   'Turbulent Flux (inc. neg. q) of Liq. Water', '')
CALL VARSETUP1C_NC (INCID(2),  IACCUR, 2, IDIMID3, &
   'liq_flx_turb_pos', 'Pseudo-Flux of Liq. Water to Correct fo ql<0', '')
CALL VARSETUP1C_NC (INCID(2),  IACCUR, 2, IDIMID3, &
   'ice_flx_turb_pos', 'Pseudo-Flux of Ice Water to Correct fo qi<0', '')

CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'gwd_flx_u',      'Gravity Wave Drag U-Flux',            '')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'gwd_flx_v',      'Gravity Wave Drag V-Flux',            '')

CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'sw_rad_flux',    'SW radiative Flux',                   'W/m^2')
CALL VARSETUP1C_NC (INCID(1),  IACCUR, 2, IDIMID3, &
   'lw_rad_flux',    'LW radiative Flux',                   'W/m^2')

DO N=1,YDDPHY%NVXTR2
  WRITE(NTXT,"(2i1)") INT(N/10), MOD(N,10)
  CALL VARSETUP1C_NC (INCID(2), IACCUR, 1, (/ NTIMDID /), &
     'extra_sfc_'//NTXT,  YDPHYDS%CVEXTR2(N),  '')
ENDDO

DO N=1,YDDPHY%NVEXTR
  WRITE(NTXT,"(2i1)") INT(N/10), MOD(N,10)
  CALL VARSETUP1C_NC (INCID(2), IACCUR, 2, IDIMIDE, & !Maike
     'extra_col_'//NTXT,  YDPHYDS%CVEXTRA(N),  '')
ENDDO


!        3.    write model level data
!              ----------------------

DO N=1,2

!        return to data mode
  ISTATUS = NF_ENDDEF(INCID(N))
  CALL HANDLE_ERR_NC(ISTATUS)

!        atmospheric model levels
  DO I=1,YDDIMV%NFLEVG
    ILEV(I)=I
  ENDDO
  ISTATUS = NF_INQ_VARID   (INCID(N), 'nlev', NVARID)
  CALL HANDLE_ERR_NC(ISTATUS)
  ISTATUS = NF_PUT_VAR_INT (INCID(N), NVARID, ILEV)
  CALL HANDLE_ERR_NC(ISTATUS)

!        atmospheric model levels + 1
  DO I=1,YDDIMV%NFLEVG+1
    ILEVP1(I)=I
  ENDDO
  ISTATUS = NF_INQ_VARID   (INCID(N), 'nlevp1', NVARID)
  CALL HANDLE_ERR_NC(ISTATUS)
  ISTATUS = NF_PUT_VAR_INT (INCID(N), NVARID,   ILEVP1)
  CALL HANDLE_ERR_NC(ISTATUS)

!        soil/sea-ice model levels
  DO I=1,ICSS
    ILEVS(I)=I
  ENDDO
  ISTATUS = NF_INQ_VARID   (INCID(N), 'nlevs', NVARID)
  CALL HANDLE_ERR_NC(ISTATUS)
  ISTATUS = NF_PUT_VAR_INT (INCID(N), NVARID, ILEVS)
  CALL HANDLE_ERR_NC(ISTATUS)

ENDDO

END ASSOCIATE
END SUBROUTINE SUD1C_NC
