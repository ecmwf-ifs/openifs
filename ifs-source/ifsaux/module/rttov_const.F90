! Description:
!> @file
!!   Defines various constants used by RTTOV
!
!> @brief
!!   Defines various constants used by RTTOV
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 7 December 2016, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, DWD and MeteoFrance.
!
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
MODULE rttov_const

  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE

  ! -------------------------------------------------------
  ! General
  ! -------------------------------------------------------
  ! Version number of the current code
  INTEGER(jpim), PARAMETER :: version = 13
  INTEGER(jpim), PARAMETER :: release = 0
  INTEGER(jpim), PARAMETER :: minor_version = 0

  ! Min/max version numbers of compatible coefficient files:
  !   rtcoef files with "id_comp_lvl" outside range will be rejected
  INTEGER(jpim), PARAMETER :: version_compatible_min = 12
  INTEGER(jpim), PARAMETER :: version_compatible_max = 13
  !   scaer/sccld files with "version" outside range will be rejected
  INTEGER(jpim), PARAMETER :: scaercld_version_compatible_min = 1
  INTEGER(jpim), PARAMETER :: scaercld_version_compatible_max = 1
  !   MFASIS LUT files with "version" outside range will be rejected
  INTEGER(jpim), PARAMETER :: mfasislut_version_compatible_min = 0
  INTEGER(jpim), PARAMETER :: mfasislut_version_compatible_max = 2
  !   PC-RTTOV files with "fmv_pc_comp_pc" outside range will be rejected
  INTEGER(jpim), PARAMETER :: pccoef_version_compatible_min = 3
  INTEGER(jpim), PARAMETER :: pccoef_version_compatible_max = 5
  !   RTTOV-SCATT hydrotable with "rttov_version" outside range will be rejected
  INTEGER(jpim), PARAMETER :: hydrotable_version_compatible_min = 13
  INTEGER(jpim), PARAMETER :: hydrotable_version_compatible_max = 13

  ! Values for identifying RTTOV coefficient files
  CHARACTER(LEN=16), PARAMETER :: rttov_magic_string = '%RTTOV_COEFF    '
  REAL(jprb),        PARAMETER :: rttov_magic_number = 1.2345E+12_jprb

  ! Standard error unit number (for HPUX standard error unit number is 7)
  INTEGER(jpim), PARAMETER :: default_err_unit = 0

  ! Error status values
  INTEGER(jpim), PARAMETER :: errorstatus_success = 0
  INTEGER(jpim), PARAMETER :: errorstatus_fatal   = 1

  ! -------------------------------------------------------
  ! Precision and numerical constants
  ! -------------------------------------------------------
  ! Try to ensure this is large enough to avoid overflows in reciprocals but small enough to not affect the calculations.
  ! these parameters are defined at bottom of module (because they use
  REAL(jprb), PARAMETER :: max_exp_exponent = 50._jprb   ! approx 1e22 which should be sufficiently big for most purposes
  REAL(jprb), PARAMETER :: min_exponent     = 1E-16_jprb ! approx log_10(1+2^-52) - anything raised to this power or smaller
                                                         ! should be approx equal to 1
  REAL(jprb), PARAMETER :: realtol          = 1E-32_jprb ! Tolerance for real (in)equality tests
  ! small_val is defined in rttov_transmit to avoid compiler incompatibility
  ! small_val is used in rttov_transmit to ensure small values do not result in underflows. In subsequent calculations
  ! these values are multiplied together hence the exponent of 1/3.
  ! REAL(jprb)            :: small_val = (tiny(min_exponent)) ** (0.333333_jprb) ! XLF doesn't like 1/3

  ! Non-zero coefficient bounds not set
  INTEGER(jpim), PARAMETER :: upper_layer_bound_not_set = 9999_jpim

  ! -------------------------------------------------------
  ! Physical constants
  ! -------------------------------------------------------
  ! Molecular weights  (g/mole) are calculated by adding NIST Standard Atomic Weights
  ! Molecular weight of dry air refers to US standard atmosphere 1976
  ! NIST  Standard Atomic Weight are:
  ! H    1.00794   (7)
  ! C   12.0107    (8)
  ! N   14.0067    (2)
  ! O   15.9994    (3)
  ! http://webbook.nist.gov/chemistry/form-ser.html
  REAL(jprb), PARAMETER :: mair   = 28.9644_jprb
  REAL(jprb), PARAMETER :: mh2o   = 18.01528_jprb
  REAL(jprb), PARAMETER :: mo3    = 47.9982_jprb
  REAL(jprb), PARAMETER :: mco2   = 44.0095_jprb
  REAL(jprb), PARAMETER :: mch4   = 16.04246_jprb
  REAL(jprb), PARAMETER :: mn2o   = 44.0128_jprb
  REAL(jprb), PARAMETER :: mco    = 28.0101_jprb
  REAL(jprb), PARAMETER :: mso2   = 64.064_jprb
  REAL(jprb), PARAMETER :: mo2    = 31.9988_jprb
  REAL(jprb), PARAMETER :: mno    = 30.0061_jprb
  REAL(jprb), PARAMETER :: mno2   = 46.0055_jprb
  REAL(jprb), PARAMETER :: mhno3  = 63.0128_jprb
  REAL(jprb), PARAMETER :: mocs   = 60.075_jprb
  REAL(jprb), PARAMETER :: mn2    = 28.0134_jprb
  REAL(jprb), PARAMETER :: mccl4  = 153.823_jprb
  REAL(jprb), PARAMETER :: mcfc11 = 137.368_jprb ! CFCl3
  REAL(jprb), PARAMETER :: mcfc12 = 120.914_jprb ! CF2Cl2
  REAL(jprb), PARAMETER :: mcfc14 = 88.0043_jprb ! CF4

  ! Ozone density at S.T.P in kg/m^2
  ! Reference: https://software.ecmwf.int/wiki/pages/viewpage.action?pageId=61121586
  ! [Ozone density at S.T.P.: 4.4615E-04 mol/m^2] * [O3 molecular weight: 47.9982E-03 kg/mol] = 2.1414E-05 kg/m^2
  REAL(jprb), PARAMETER :: rho_ozone_stp = 2.1414E-05_jprb

  ! Density of liquid water [kg/kg]
  REAL(jprb), PARAMETER :: rho_water = 1.E3_jprb

#ifdef _RTTOV12
  ! Fundamental constants taken from http://physics.nist.gov/cuu/index.html
  ! P.J. Mohr, B.N. Taylor, and D.B. Newell (2015), "The 2014 CODATA Recommended
  ! Values of the Fundamental Physical Constants" (Web Version 7.0)
  ! Values as of 11/12/2015

  ! Avogadro constant (mol^-1)
  REAL(jprb), PARAMETER :: na = 6.022140857E23_jprb

  ! Acceleration due to gravity (m/s^2, exact)
  REAL(jprb), PARAMETER :: gravity = 9.80665_jprb

  ! c1 = 2hc**2; c2 = hc/k
  REAL(jprb), PARAMETER :: planck_c1 = .00001191042953_jprb ! * mW/(m2 sr cm-4)
  REAL(jprb), PARAMETER :: planck_c2 = 1.4387774_jprb       ! cm K
  REAL(jprb), PARAMETER :: speedl = 29979245800.0_jprb      ! Speed of light cm/s (exact)

  ! Universal gas constant (J/mol/K)
  REAL(jprb), PARAMETER :: rgc = 8.3144598_jprb

  ! Boltzman constant (J K-1)  http://physics.nist.gov February 2017
  REAL(jprb), PARAMETER :: kb = 1.38064852E-23_jprb
#else
  ! Fundamental constants taken from http://physics.nist.gov/constants
  ! Eite Tiesinga, Peter J. Mohr, David B. Newell, and Barry N. Taylor (2019),
  ! "The 2018 CODATA Recommended Values of the Fundamental Physical Constants"
  ! (Web Version 8.0).

  ! Avogadro constant (mol^-1, exact)
  REAL(jprb), PARAMETER :: na = 6.02214076E23_jprb

  ! Universal gas constant (J/mol/K, exact)
  REAL(jprb), PARAMETER :: rgc = 8.314462618_jprb

  ! Acceleration due to gravity (m/s^2, exact)
  REAL(jprb), PARAMETER :: gravity = 9.80665_jprb

  ! Boltzman constant (J/K, exact)
  REAL(jprb), PARAMETER :: kb = 1.380649E-23_jprb

  ! Speed of light (cm/s, exact)
  REAL(jprb), PARAMETER :: speedl = 29979245800.0_jprb

  ! Planck constant: 6.62607015E-34 (J/Hz, exact, not used in RTTOV)

  ! c1 = 2hc**2; c2 = hc/k (h, c, k values are exact; c1, c2 truncated at 10sf)
  REAL(jprb), PARAMETER :: planck_c1 = 1.191042972E-05_jprb ! * mW/(m2 sr cm^-4)
  REAL(jprb), PARAMETER :: planck_c2 = 1.438776878_jprb     ! cm K
#endif

  ! Kaye & Laby latest library edition is 16e 1995, and gives
  ! * earth mean radius r= 6371.00 km (p191)
  !    [defined as [(r_equator)^2 (r_pole)]^1/3]
  REAL(jprb), PARAMETER :: pi            = 3.1415926536_jprb
  REAL(jprb), PARAMETER :: deg2rad       = pi / 180.0_jprb
  REAL(jprb), PARAMETER :: rad2deg       = 1._jprb / deg2rad
  REAL(jprb), PARAMETER :: earthradius   = 6371.00_jprb
  REAL(jprb), PARAMETER :: flatt         = 3.3528107E-3_jprb
  REAL(jprb), PARAMETER :: omega         = 7292115E-11_jprb
  REAL(jprb), PARAMETER :: eqrad         = 6378.137_jprb
  REAL(jprb), PARAMETER :: grave         = 9.7803267715_jprb
  REAL(jprb), PARAMETER :: pi_r          = 1._jprb / pi
  REAL(jprb), PARAMETER :: z4pi_r        = 0.25_jprb * pi_r
  REAL(jprb), PARAMETER :: sec_theta_eff = 1.743446796_jprb   ! theta_eff = 55 degrees

  ! The Cosmic Microwave Background Spectrum from the Full COBE FIRAS Data Set
  ! Fixsen D.J. et al Astrophysical Journal v.473, p.576 December 1996
  ! CMBR = 2.728 +- 0.004K
  REAL(jprb), PARAMETER :: tcosmic = 2.728_jprb

  ! zero temperature (K)
  REAL(jprb), PARAMETER :: t0 = 273.15_jprb

  ! standard pressure (hPa)
  REAL(jprb), PARAMETER :: p0 = 1013.25_jprb

  ! -------------------------------------------------------
  ! Satellite and instrument information
  ! -------------------------------------------------------
  ! Platform ID codes
  INTEGER(jpim), PARAMETER :: nplatforms = 71
  ! NB ID 22 was changed from "insat_3d" - use ID 40 for INSAT3
  !    IDs 41 and 45 are for ground-based and airborne sensors (experimental)
  INTEGER(jpim), PARAMETER :: &
         platform_id_noaa      = 1,  platform_id_dmsp      = 2,  platform_id_meteosat  = 3,  &
         platform_id_goes      = 4,  platform_id_gms       = 5,  platform_id_fy2       = 6,  &
         platform_id_trmm      = 7,  platform_id_ers       = 8,  platform_id_eos       = 9,  &
         platform_id_metop     = 10, platform_id_envisat   = 11, platform_id_msg       = 12, &
         platform_id_fy1       = 13, platform_id_adeos     = 14, platform_id_mtsat     = 15, &
         platform_id_coriolis  = 16, platform_id_jpss      = 17, platform_id_gifts     = 18, &
         platform_id_sentinel3 = 19, platform_id_meghatr   = 20, platform_id_kalpana   = 21, &
         platform_id_meteor    = 22, platform_id_fy3       = 23, platform_id_coms      = 24, &
         platform_id_meteorm   = 25, platform_id_gosat     = 26, platform_id_calipso   = 27, &
         platform_id_dummy     = 28, platform_id_gcomw     = 29, platform_id_nimbus    = 30, &
         platform_id_himawari  = 31, platform_id_mtg       = 32, platform_id_saral     = 33, &
         platform_id_metopsg   = 34, platform_id_landsat   = 35, platform_id_jason     = 36, &
         platform_id_gpm       = 37, platform_id_insat1    = 38, platform_id_insat2    = 39, &
         platform_id_insat3    = 40, platform_id_ground    = 41, platform_id_dscovr    = 42, &
         platform_id_clarreo   = 43, platform_id_ticfire   = 44, platform_id_aircraft  = 45, &
         platform_id_iss       = 46, platform_id_hj1       = 47, platform_id_gkompsat2 = 48, &
         platform_id_gcomc     = 49, platform_id_smos      = 50, platform_id_ors       = 51, &
         platform_id_fy4       = 52, platform_id_tropics   = 53, platform_id_gf5       = 54, &
         platform_id_hy2       = 55, platform_id_cloudsat  = 56, platform_id_cloudcore = 57, &
         platform_id_forum     = 58, platform_id_tempest   = 59, platform_id_jasoncs   = 60, &
         platform_id_electrol  = 61, platform_id_oms       = 62, platform_id_sentinel2 = 63, &
         platform_id_cloudctrl = 64, platform_id_trishna   = 65, platform_id_sbg       = 66, &
         platform_id_saphirhy  = 67, platform_id_aws       = 68, platform_id_cmimmw    = 69, &
         platform_id_earthcare = 70, platform_id_sentinl5p = 71

  ! Platform names
  INTEGER(jpim), PARAMETER :: len_platform = 9_jpim
  CHARACTER(LEN=len_platform), PARAMETER :: platform_name(nplatforms) = &
         (/ 'noaa     ', 'dmsp     ', 'meteosat ', 'goes     ', 'gms      ', &
            'fy2      ', 'trmm     ', 'ers      ', 'eos      ', 'metop    ', &
            'envisat  ', 'msg      ', 'fy1      ', 'adeos    ', 'mtsat    ', &
            'coriolis ', 'jpss     ', 'gifts    ', 'sentinel3', 'meghatr  ', &
            'kalpana  ', 'meteor   ', 'fy3      ', 'coms     ', 'meteor-m ', &
            'gosat    ', 'calipso  ', 'dummy    ', 'gcom-w   ', 'nimbus   ', &
            'himawari ', 'mtg      ', 'saral    ', 'metopsg  ', 'landsat  ', &
            'jason    ', 'gpm      ', 'insat1   ', 'insat2   ', 'insat3   ', &
            'ground   ', 'dscovr   ', 'clarreo  ', 'ticfire  ', 'aircraft ', &
            'iss      ', 'hj1      ', 'gkompsat2', 'gcom-c   ', 'smos     ', &
            'ors      ', 'fy4      ', 'tropics  ', 'gf5      ', 'hy2      ', &
            'cloudsat ', 'cloudcore', 'forum    ', 'tempest  ', 'jasoncs  ', &
            'electro-l', 'oms      ', 'sentinel2', 'cloudctrl', 'trishna  ', &
            'sbg      ', 'saphirhy ', 'aws      ', 'cmimmw   ', 'earthcare', &
            'sentinl5p' /)

  ! Instrument ID codes
  INTEGER(jpim), PARAMETER :: &
         inst_id_hirs       =   0, inst_id_msu      =   1, inst_id_ssu      =   2, inst_id_amsua    =   3, &
         inst_id_amsub      =   4, inst_id_avhrr    =   5, inst_id_ssmi     =   6, inst_id_vtpr1    =   7, &
         inst_id_vtpr2      =   8, inst_id_tmi      =   9, inst_id_ssmis    =  10, inst_id_airs     =  11, &
         inst_id_hsb        =  12, inst_id_modis    =  13, inst_id_atsr     =  14, inst_id_mhs      =  15, &
         inst_id_iasi       =  16, inst_id_amsre    =  17, inst_id_gmsim    =  18, inst_id_atms     =  19, &
         inst_id_mviri      =  20, inst_id_seviri   =  21, inst_id_goesim   =  22, inst_id_goessd   =  23, &
         inst_id_mtsatim    =  24, inst_id_vissr    =  25, inst_id_mvisr    =  26, inst_id_cris     =  27, &
         inst_id_crisfsr    =  28, inst_id_viirs    =  29, inst_id_windsat  =  30, inst_id_gifts    =  31, &
         inst_id_ssmt1      =  32, inst_id_ssmt2    =  33, inst_id_saphir   =  34, inst_id_madras   =  35, &
         inst_id_ssmisz     =  36, inst_id_vhrr     =  37, inst_id_insatim  =  38, inst_id_insatsd  =  39, &
         inst_id_mwts       =  40, inst_id_mwhs     =  41, inst_id_iras     =  42, inst_id_mwri     =  43, &
         inst_id_abi        =  44, inst_id_mi       =  45, inst_id_msumr    =  46, inst_id_tansofts =  47, &
         inst_id_iir        =  48, inst_id_mwr      =  49, inst_id_dummyir  =  50, inst_id_dummymw  =  51, &
         inst_id_dummyhi    =  52, inst_id_dummypo  =  53, inst_id_scams    =  54, inst_id_smmr     =  55, &
         inst_id_ahi        =  56, inst_id_irs      =  57, inst_id_altika   =  58, inst_id_iasing   =  59, &
         inst_id_tm         =  60, inst_id_fci      =  61, inst_id_amsr1    =  62, inst_id_amsr2    =  63, &
         inst_id_vissr2     =  64, inst_id_slstr    =  65, inst_id_tirs     =  66, inst_id_amr      =  67, &
         inst_id_oli        =  68, inst_id_iris     =  69, inst_id_ici      =  70, inst_id_gmi      =  71, &
         inst_id_mwts2      =  72, inst_id_mwhs2    =  73, inst_id_aster    =  74, inst_id_hatpro   =  75, &
         inst_id_mtvzagy    =  76, inst_id_metimage =  77, inst_id_mws      =  78, inst_id_mwi      =  79, &
         inst_id_epic       =  80, inst_id_mrir     =  81, inst_id_si       =  82, inst_id_mrfirs   =  83, &
         inst_id_mbfiri     =  84, inst_id_lhr      =  85, inst_id_ismar    =  86, inst_id_mersi1   =  87, &
         inst_id_mersi2     =  88, inst_id_ecostres =  89, inst_id_irmss    =  90, inst_id_olci     =  91, &
         inst_id_thir       =  92, inst_id_ami      =  93, inst_id_ikfs2    =  94, inst_id_li       =  95, &
         inst_id_sgli       =  96, inst_id_hiras    =  97, inst_id_giirs    =  98, inst_id_agri     =  99, &
         inst_id_pmr        = 100, inst_id_miras    = 101, inst_id_cowvr    = 102, inst_id_tropics  = 103, &
         inst_id_vims       = 104, inst_id_dpr      = 105, inst_id_hy2mwri  = 106, inst_id_cpr      = 107, &
         inst_id_saphirng   = 108, inst_id_forum    = 109, inst_id_tempest  = 110, inst_id_marss    = 111, &
         inst_id_virr       = 112, inst_id_sirs     = 113, inst_id_hirasfsr = 114, inst_id_hrir     = 115, &
         inst_id_amrc       = 116, inst_id_msugs    = 117, inst_id_hytes    = 118, inst_id_gems1    = 119, &
         inst_id_sntnl2_msi = 120, inst_id_hyms     = 121, inst_id_tir      = 122, inst_id_sbg      = 123, &
         inst_id_gome2      = 124, inst_id_saphirhy = 125, inst_id_aws      = 126, inst_id_cmimmw   = 127, &
         inst_id_ecare_msi  = 128, inst_id_tropomi  = 129, inst_id_sentinl5 = 130, inst_id_mersill  = 131, &
         inst_id_mwts3      = 132, inst_id_hiras2   = 133

  INTEGER(jpim), PARAMETER :: ninst = 134
  ! List of instruments - NB HIRS is number 0
  INTEGER(jpim), PARAMETER :: len_instrument = 8_jpim
  CHARACTER(LEN=len_instrument), PARAMETER :: inst_name(0:ninst-1) =     &
          (/ 'hirs    ', 'msu     ', 'ssu     ', 'amsua   ', 'amsub   ', &
             'avhrr   ', 'ssmi    ', 'vtpr1   ', 'vtpr2   ', 'tmi     ', &
             'ssmis   ', 'airs    ', 'hsb     ', 'modis   ', 'atsr    ', &
             'mhs     ', 'iasi    ', 'amsre   ', 'imager  ', 'atms    ', &
             'mviri   ', 'seviri  ', 'imager  ', 'sounder ', 'imager  ', &
             'vissr   ', 'mvisr   ', 'cris    ', 'cris-fsr', 'viirs   ', &
             'windsat ', 'gifts   ', 'ssmt1   ', 'ssmt2   ', 'saphir  ', &
             'madras  ', 'ssmisz  ', 'vhrr    ', 'imager  ', 'sounder ', &
             'mwts    ', 'mwhs    ', 'iras    ', 'mwri    ', 'abi     ', &
             'mi      ', 'msumr   ', 'tansofts', 'iir     ', 'mwr     ', &
             'dummyir ', 'dummymw ', 'dummyhi ', 'dummypo ', 'scams   ', &
             'smmr    ', 'ahi     ', 'irs     ', 'altika  ', 'iasing  ', &
             'tm      ', 'fci     ', 'amsr    ', 'amsr2   ', 'vissr   ', &
             'slstr   ', 'tirs    ', 'amr     ', 'oli     ', 'iris    ', &
             'ici     ', 'gmi     ', 'mwts2   ', 'mwhs2   ', 'aster   ', &
             'hatpro  ', 'mtvzagy ', 'metimage', 'mws     ', 'mwi     ', &
             'epic    ', 'mrir    ', 'si      ', 'mrfirs  ', 'mbfiri  ', &
             'lhr     ', 'ismar   ', 'mersi1  ', 'mersi2  ', 'ecostres', &
             'irmss   ', 'olci    ', 'thir    ', 'ami     ', 'ikfs2   ', &
             'li      ', 'sgli    ', 'hiras   ', 'giirs   ', 'agri    ', &
             'pmr     ', 'miras   ', 'cowvr   ', 'tropics ', 'vims    ', &
             'dpr     ', 'mwri    ', 'cpr     ', 'saphirng', 'forum   ', &
             'tempest ', 'marss   ', 'virr    ', 'sirs    ', 'hirasfsr', &
             'hrir    ', 'amrc    ', 'msugs   ', 'hytes   ', 'gems1   ', &
             'msi     ', 'hyms    ', 'tir     ', 'sbg     ', 'gome2   ', &
             'saphirhy', 'aws     ', 'cmimmw  ', 'msi     ', 'tropomi ', &
             'sentinl5', 'mersill ', 'mwts3   ', 'hiras2  ' /)

  ! Sensor type ID codes
  INTEGER(jpim), PARAMETER :: nsensors = 4
  INTEGER(jpim), PARAMETER :: &
         sensor_id_ir     = 1, &
         sensor_id_mw     = 2, &
         sensor_id_hi     = 3, &
         sensor_id_po     = 4

  ! Sensor type names
  CHARACTER(LEN=2), PARAMETER :: sensor_name(nsensors) = &
         (/ 'ir', 'mw', 'hi', 'po' /)

  ! Sensor types corresponding to entries in the inst_name array
  INTEGER(jpim), PARAMETER :: sensor_id(0:ninst-1) = (/ &
    sensor_id_ir, sensor_id_mw, sensor_id_ir, sensor_id_mw, sensor_id_mw, &
    sensor_id_ir, sensor_id_mw, sensor_id_ir, sensor_id_ir, sensor_id_mw, &
    sensor_id_mw, sensor_id_hi, sensor_id_mw, sensor_id_ir, sensor_id_ir, &
    sensor_id_mw, sensor_id_hi, sensor_id_mw, sensor_id_ir, sensor_id_mw, &
    sensor_id_ir, sensor_id_ir, sensor_id_ir, sensor_id_ir, sensor_id_ir, &
    sensor_id_ir, sensor_id_ir, sensor_id_hi, sensor_id_hi, sensor_id_ir, &
    sensor_id_po, sensor_id_hi, sensor_id_mw, sensor_id_mw, sensor_id_mw, &
    sensor_id_mw, sensor_id_mw, sensor_id_ir, sensor_id_ir, sensor_id_ir, &
    sensor_id_mw, sensor_id_mw, sensor_id_ir, sensor_id_mw, sensor_id_ir, &
    sensor_id_ir, sensor_id_ir, sensor_id_hi, sensor_id_ir, sensor_id_mw, &
    sensor_id_ir, sensor_id_mw, sensor_id_hi, sensor_id_po, sensor_id_mw, &
    sensor_id_mw, sensor_id_ir, sensor_id_hi, sensor_id_mw, sensor_id_hi, &
    sensor_id_ir, sensor_id_ir, sensor_id_mw, sensor_id_mw, sensor_id_ir, &
    sensor_id_ir, sensor_id_ir, sensor_id_mw, sensor_id_ir, sensor_id_hi, &
    sensor_id_mw, sensor_id_mw, sensor_id_mw, sensor_id_mw, sensor_id_ir, &
    sensor_id_mw, sensor_id_mw, sensor_id_ir, sensor_id_mw, sensor_id_mw, &
    sensor_id_ir, sensor_id_ir, sensor_id_hi, sensor_id_hi, sensor_id_ir, &
    sensor_id_hi, sensor_id_mw, sensor_id_ir, sensor_id_ir, sensor_id_ir, &
    sensor_id_ir, sensor_id_ir, sensor_id_ir, sensor_id_ir, sensor_id_hi, &
    sensor_id_ir, sensor_id_ir, sensor_id_hi, sensor_id_hi, sensor_id_ir, &
    sensor_id_ir, sensor_id_mw, sensor_id_po, sensor_id_mw, sensor_id_ir, &
    sensor_id_mw, sensor_id_mw, sensor_id_mw, sensor_id_mw, sensor_id_hi, &
    sensor_id_mw, sensor_id_mw, sensor_id_ir, sensor_id_ir, sensor_id_hi, &
    sensor_id_ir, sensor_id_mw, sensor_id_ir, sensor_id_hi, sensor_id_mw, &
    sensor_id_ir, sensor_id_mw, sensor_id_ir, sensor_id_ir, sensor_id_hi, &
    sensor_id_mw, sensor_id_mw, sensor_id_mw, sensor_id_ir, sensor_id_hi, &
    sensor_id_hi, sensor_id_ir, sensor_id_mw, sensor_id_hi /)

  ! -------------------------------------------------------
  ! ASCII coefficient file section names
  ! -------------------------------------------------------
  INTEGER(jpim), PARAMETER :: nsections = 48
  INTEGER(jpim), PARAMETER :: lensection = 34
  CHARACTER(LEN=lensection), PARAMETER :: section_types(nsections) = &
      (/ 'IDENTIFICATION                    ', 'LINE-BY-LINE                      ', &
         'FAST_MODEL_VARIABLES              ', 'FILTER_FUNCTIONS                  ', &
         'FUNDAMENTAL_CONSTANTS             ', 'SSIREM                            ', &
         'FASTEM                            ', 'REFERENCE_PROFILE                 ', &
         'PROFILE_LIMITS                    ', 'FAST_COEFFICIENTS                 ', &
         'COEF_SUB_FILES                    ', 'GAZ_UNITS                         ', &
         'DIMENSIONS                        ', 'FREQUENCIES                       ', &
         'HYDROMETEOR                       ', 'CONVERSIONS                       ', &
         'EXTINCTION                        ', 'ALBEDO                            ', &
         'ASYMMETRY                         ', 'GAS_SPECTRAL_INTERVAL             ', &
         'TRANSMITTANCE_TRESHOLD            ', 'SOLAR_SPECTRUM                    ', &
         'WATER_OPTICAL_CONSTANT            ', 'WAVE_SPECTRUM                     ', &
         'AEROSOLS                          ', 'WATER_CLOUD_OPAC                  ', &
         'WATER_CLOUD_DEFF                  ', 'ICE_CLOUD_BAUM                    ', &
         'PRINCOMP_PREDICTORS               ', 'PRINCOMP_EIGENVECTORS             ', &
         'PRINCOMP_COEFFICIENTS             ', 'EMISSIVITY_COEFFICIENTS           ', &
         'PC_REFERENCE_PROFILE              ', 'PC_PROFILE_LIMITS                 ', &
         'INSTRUMENT_NOISE                  ', 'PLANCK_WEIGHTED                   ', &
         'SOLAR_FAST_COEFFICIENTS           ', 'README_SPECTRAL_RESPONSE_FUNCTION ', &
         'NLTE_RADIANCE_COEFS               ', 'PRESSURE_MODULATED_CELL           ', &
         'IR_SEA_EMIS                       ', 'PROFILE_ENVELOPE                  ', &
         'MFASIS_GENERAL                    ', 'README_MFASIS                     ', &
         'MFASIS_DIMS                       ', 'MFASIS_LUTS                       ', &
         'REFLECTIVITY                      ', 'POLARISATION                      ' /)

  ! -------------------------------------------------------
  ! Gas information
  ! -------------------------------------------------------
  ! Gas ID codes
  INTEGER(jpim), PARAMETER :: ngases_max = 9
  INTEGER(jpim), PARAMETER :: &
      gas_id_mixed       =  1, gas_id_watervapour =  2, gas_id_ozone       =  3, &
      gas_id_wvcont      =  4, gas_id_co2         =  5, gas_id_n2o         =  6, &
      gas_id_co          =  7, gas_id_ch4         =  8, gas_id_so2         =  9!, &
      ! gas_id_o2          = 10, gas_id_no          = 11, gas_id_no2         = 12, &
      ! gas_id_hno3        = 13, gas_id_ocs         = 14, gas_id_n2          = 15, &
      ! gas_id_ccl4        = 16, gas_id_cfc11       = 17, gas_id_cfc12       = 18, &
      ! gas_id_cfc14       = 19

  ! Gas names
  CHARACTER(LEN=12), PARAMETER :: gas_name(ngases_max) = &
      (/ 'Mixed_gases ', 'Water_vapour', 'Ozone       ', 'WV_Continuum', &
         'CO2         ', 'N2O         ', 'CO          ', 'CH4         ', &
         'SO2         '/) !, 'O2          ','NO          ', 'NO2         ', &
         ! 'HNO3        ', 'OCS         ', 'N2          ', 'CCl4        ', &
         ! 'CFC-11      ', 'CFC-12      ', 'CFC-14      ' /)

  ! Convenient array of gas molecular masses
  REAL(jprb), PARAMETER :: gas_mass(ngases_max) = &
      (/ mair, mh2o, mo3, mh2o, mco2, mn2o, mco, mch4, mso2 /) !, mo2, mno, &
         ! mno2, mhno3, mocs, mn2, mccl4, mcfc11, mcfc12, mcfc14 /)

  ! -------------------------------------------------------
  ! Hard limits for control of input profile
  ! -------------------------------------------------------
  ! Temperature
  REAL(jprb), PARAMETER :: tskin_land_max = 1250.0_jprb ! K (for Tskin over land only - lavas, fires, etc)
  REAL(jprb), PARAMETER :: tmax           = 400.0_jprb  ! K (for all other input temperatures)
  REAL(jprb), PARAMETER :: tmin           = 90.0_jprb   ! K (for all input temperatures)
  ! Water Vapour
  REAL(jprb), PARAMETER :: qmax   = 0.60E+06_jprb    ! ppmv 0.373_jprb kg/kg
  REAL(jprb), PARAMETER :: qmin   = 0.1E-10_jprb     ! ppmv
  ! Ozone
  REAL(jprb), PARAMETER :: o3max  = 1000.0_jprb      ! ppmv  1.657E-3_jprb kg/kg
  REAL(jprb), PARAMETER :: o3min  = 0.1E-10_jprb     ! ppmv
  ! CO2
  REAL(jprb), PARAMETER :: co2max = 1000.0_jprb      ! ppmv
  REAL(jprb), PARAMETER :: co2min = 0.1E-10_jprb     ! ppmv
  ! CO
  REAL(jprb), PARAMETER :: comax  = 10.0_jprb        ! ppmv
  REAL(jprb), PARAMETER :: comin  = 0.1E-10_jprb     ! ppmv
  ! N2O
  REAL(jprb), PARAMETER :: n2omax = 10.0_jprb        ! ppmv
  REAL(jprb), PARAMETER :: n2omin = 0.1E-10_jprb     ! ppmv
  ! CH4
  REAL(jprb), PARAMETER :: ch4max = 50.0_jprb        ! ppmv
  REAL(jprb), PARAMETER :: ch4min = 0.1E-10_jprb     ! ppmv
  ! SO2
  REAL(jprb), PARAMETER :: so2max = 1000.0_jprb      ! ppmv
  REAL(jprb), PARAMETER :: so2min = 0.1E-10_jprb     ! ppmv
  ! Cloud Liquid Water
  REAL(jprb), PARAMETER :: clwmax = 1.0_jprb         ! kg/kg
  REAL(jprb), PARAMETER :: clwmin = 0.0_jprb         ! kg/kg
  ! Surface Pressure
  REAL(jprb), PARAMETER :: pmax   = 1100.0_jprb      ! surface pressure hPa
  REAL(jprb), PARAMETER :: pmin   = 400.0_jprb       ! hPa
  ! Surface Wind
  REAL(jprb), PARAMETER :: wmax   =  100.0_jprb      ! surface wind speed (m/s)
  ! Surface elevation
  REAL(jprb), PARAMETER :: elevmax = 20.0_jprb       ! km
  ! Zenith Angle
  REAL(jprb), PARAMETER :: zenmax = 75.0_jprb        ! zenith angle (Deg) = secant 3.86_jprb
  REAL(jprb), PARAMETER :: zenmaxv9 = 85.3_jprb      ! larger zenmax for v9 predictors = secant 12
  ! Cloud Top Pressure
  REAL(jprb), PARAMETER :: ctpmax = 1100.0_jprb      ! (hPa)
  REAL(jprb), PARAMETER :: ctpmin =   50.0_jprb      ! (hPa)
  ! Magnetic field strength
  REAL(jprb), PARAMETER :: bemax = 0.7_jprb          ! (Gauss)
  REAL(jprb), PARAMETER :: bemin = 0.2_jprb          ! (Guass)
  ! CLW diameter
  REAL(jprb), PARAMETER :: dgmin_clw = 2._jprb       ! (micron)
  REAL(jprb), PARAMETER :: dgmax_clw = 52._jprb      ! (micron)
  ! Ice crystal diameter
  REAL(jprb), PARAMETER :: dgmin_baum = 10._jprb     ! (micron)
  REAL(jprb), PARAMETER :: dgmax_baum = 120._jprb    ! (micron)
  ! Ice Water Content
  REAL(jprb), PARAMETER :: iwcmin_baum = 4.984E-5_jprb ! (g.m-3)
  REAL(jprb), PARAMETER :: iwcmax_baum = 0.1831_jprb   ! (g.m-3)

  ! -------------------------------------------------------
  ! Historical and default scaling of T and gas envelope to create limits
  ! -------------------------------------------------------
  REAL(jprb), PARAMETER :: tscale_v11 = 0.1_jprb
  REAL(jprb), PARAMETER :: gscale_v11 = 0.2_jprb
  REAL(jprb), PARAMETER :: tscale_def = 0.1_jprb
  REAL(jprb), PARAMETER :: gscale_def = 0.2_jprb

  ! -------------------------------------------------------
  ! Min/max optical depth and transmittance values
  ! -------------------------------------------------------
  ! maximum value of optical depth for transmittance calculation
  ! e(-30) -> 10**-14
  ! e(-50) -> 10**-22
#ifdef PARKIND1_SINGLE
  REAL(jprb), PARAMETER :: max_optical_depth = 30._jprb
#else
  REAL(jprb), PARAMETER :: max_optical_depth = 50._jprb
#endif
  REAL(jprb), PARAMETER :: min_tau = 1.0e-8_jprb
  REAL(jprb), PARAMETER :: min_od  = 1.0e-5_jprb         ! Also used by DOM

  ! -------------------------------------------------------
  ! Maximum solar zenith angle for which to apply solar calculation
  ! -------------------------------------------------------
  REAL(jprb), PARAMETER  :: max_sol_zen = 85.3_jprb ! = secant 12

  ! -------------------------------------------------------
  ! RTTOV auxiliary parameters
  ! -------------------------------------------------------
  ! Quality flag bits
  !   Optical depth prediction:
  INTEGER(jpim), PARAMETER  :: qflag_reg_limits             = 0  ! Gas op dep regression min/max limits exceeded
  INTEGER(jpim), PARAMETER  :: qflag_pc_aer_reg_limits      = 1  ! PC-RTTOV aerosol regression min/max limits exceeded
  !   MFASIS:
  INTEGER(jpim), PARAMETER  :: qflag_mfasis_zenangle        = 16 ! Max valid zenith angle exceeded
  INTEGER(jpim), PARAMETER  :: qflag_mfasis_sumzenangle     = 17 ! Max valid sum of zenith angles exceeded
  INTEGER(jpim), PARAMETER  :: qflag_mfasis_geometry_bounds = 20 ! Scattering angle out of LUT bounds
  INTEGER(jpim), PARAMETER  :: qflag_mfasis_opdpedia_bounds = 21 ! Some opdp/effdia dimension out of LUT bounds

  ! Gas units
  INTEGER(jpim), PARAMETER :: ngases_unit = 2
  INTEGER(jpim), PARAMETER :: &
          gas_unit_ppmvdry   = 0, &  ! volume mixing ratio [10^6 mol/mol] aka [ppmv] over dry air
          gas_unit_specconc  = 1, &  ! mass mixing ratio   [kg/kg]                   over wet air
          gas_unit_ppmv      = 2     ! volume mixing ratio [10^6 mol/mol] aka [ppmv] over wet air

  ! Interpolation modes
  !       MODE                     USER->COEF LEVEL          COEF->USER LEVEL
  !                            (profile interpolation)   (optical depth/weighting fn interpolation)
  ! interp_rochon:                    Rochon                 Rochon,     op dep
  ! interp_loglinear:                 Log-linear             Log-linear, op dep
  ! interp_rochon_loglinear:          Rochon                 Log-linear, op dep
  ! interp_rochon_wfn:                Rochon                 Rochon,     weighting fns
  ! interp_rochon_loglinear_wfn:      Rochon                 Log-linear, weighting fns
  INTEGER(jpim), PARAMETER :: ninterp_modes               = 5_jpim ! Number of valid interpolation options
  INTEGER(jpim), PARAMETER :: interp_rochon               = 1_jpim
  INTEGER(jpim), PARAMETER :: interp_loglinear            = 2_jpim
  INTEGER(jpim), PARAMETER :: interp_rochon_loglinear     = 3_jpim
  INTEGER(jpim), PARAMETER :: interp_rochon_wfn           = 4_jpim
  INTEGER(jpim), PARAMETER :: interp_rochon_loglinear_wfn = 5_jpim

  ! Surface types
  INTEGER(jpim), PARAMETER :: nsurftype = 2
  INTEGER(jpim), PARAMETER :: surftype_land = 0
  INTEGER(jpim), PARAMETER :: surftype_sea = 1
  INTEGER(jpim), PARAMETER :: surftype_seaice = 2

  ! Water types
  INTEGER(jpim), PARAMETER :: nwatertype = 1
  INTEGER(jpim), PARAMETER :: watertype_fresh_water = 0
  INTEGER(jpim), PARAMETER :: watertype_ocean_water = 1

  ! Emissivity model variables
  INTEGER(jpim), PARAMETER :: max_ir_sea_emis_model = 2 ! Highest IR sea emis model number available
  INTEGER(jpim), PARAMETER :: max_fastem_version = 6    ! Highest FASTEM version number available
  INTEGER(jpim), PARAMETER :: fastem_sp = 5             ! Number of fastem surface parameters

  ! Solar sea BRDF model variables
  INTEGER(jpim), PARAMETER :: max_solar_sea_brdf_model = 2 ! Highest solar sea BRDF model number available
  REAL(jprb),    PARAMETER :: min_windsp = 1.E-6_jprb      ! Minimum wind speed for sunglint model calculation

  ! Constants from Elfouhaily wave spectrum parameterisation
  ! (Numbers in parantheses refer to equations in the paper - see rttov_refsun.F90)
  INTEGER(jpim), PARAMETER :: nk_elf = 318                      ! Size of grid for Elfouhaily wave spectrum
  REAL(jprb),    PARAMETER :: cm_elf = 0.23_jprb                ! see just below (40)
  REAL(jprb),    PARAMETER :: km_elf = 2 * gravity / cm_elf**2  ! see just below (40) - this is 370 rad/m, near (24)
  REAL(jprb),    PARAMETER :: x0_elf = 2.2E4_jprb               ! (37) - for calculating kp
  REAL(jprb),    PARAMETER :: a0_elf = LOG(2._jprb) / 4         ! (59) - for calculating angular spread
  REAL(jprb),    PARAMETER :: ap_elf = 4                        ! (59) - for calculating angular spread

  ! MW CLW absorption options
  INTEGER(jpim), PARAMETER :: max_mw_clw_scheme = 3 ! Max valid value for MW clw_scheme
  INTEGER(jpim), PARAMETER :: mw_clw_scheme_liebe      = 1
  INTEGER(jpim), PARAMETER :: mw_clw_scheme_rosenkranz = 2
  INTEGER(jpim), PARAMETER :: mw_clw_scheme_tkc        = 3

  ! Liebe CLW absorption parameters
  REAL(jprb), PARAMETER :: dcoeff(8) = & ! Debye coefs
         (/ 17.1252_jprb, 134.2450_jprb, 310.2125_jprb,  5.667_jprb, &
           188.7979_jprb,  80.5419_jprb,   0.1157_jprb,  4.8417_jprb/)

  ! Coefficients for TKC permittivity parameterisation
  REAL(jprb), PARAMETER :: tkc_a_1 = 8.111E+1_jprb
  REAL(jprb), PARAMETER :: tkc_b_1 = 4.434E-3_jprb
  REAL(jprb), PARAMETER :: tkc_c_1 = 1.302E-13_jprb
  REAL(jprb), PARAMETER :: tkc_d_1 = 6.627E+2_jprb
  REAL(jprb), PARAMETER :: tkc_a_2 = 2.025E+0_jprb
  REAL(jprb), PARAMETER :: tkc_b_2 = 1.073E-2_jprb
  REAL(jprb), PARAMETER :: tkc_c_2 = 1.012E-14_jprb
  REAL(jprb), PARAMETER :: tkc_d_2 = 6.089E+2_jprb
  REAL(jprb), PARAMETER :: tkc_t_c = 1.342E+2_jprb
  REAL(jprb), PARAMETER :: tkc_s_0 = 8.7914E+1_jprb
  REAL(jprb), PARAMETER :: tkc_s_1 = -4.0440E-1_jprb
  REAL(jprb), PARAMETER :: tkc_s_2 = 9.5873E-4_jprb
  REAL(jprb), PARAMETER :: tkc_s_3 = -1.3280E-6_jprb

  ! Parameters to compute refractive index of air
  REAL(jprb), PARAMETER :: d1    = 8341.87_jprb
  REAL(jprb), PARAMETER :: d2    = 2405955.0_jprb
  REAL(jprb), PARAMETER :: d3    = 130.0_jprb
  REAL(jprb), PARAMETER :: d4    = 15996.0_jprb
  REAL(jprb), PARAMETER :: d5    = 38.9_jprb
  REAL(jprb), PARAMETER :: dco2  = 0.540_jprb
  REAL(jprb), PARAMETER :: ed1   = 96095.43_jprb
  REAL(jprb), PARAMETER :: ed2   = 0.601_jprb
  REAL(jprb), PARAMETER :: ed3   = 0.00972_jprb
  REAL(jprb), PARAMETER :: ed4   = 0.003661_jprb
  REAL(jprb), PARAMETER :: ew1   = 3.7345_jprb
  REAL(jprb), PARAMETER :: ew2   = 0.0401_jprb
  REAL(jprb), PARAMETER :: htop  = 100.0_jprb
  REAL(jprb), PARAMETER :: ctom  = 1.0E-4_jprb
  REAL(jprb), PARAMETER :: waver = 1700.0_jprb

  ! CO2 concentration assumed for atmospheric refractivity calculation
  REAL(jprb), PARAMETER :: co2_conc = 376._jprb

  ! Parameters for solar overcast radiance calculation
  REAL(jprb), PARAMETER :: overcast_albedo_wvn = 10000._jprb ! Wavenumber (cm-1) at which albedo changes
  REAL(jprb), PARAMETER :: overcast_albedo1    = 0.7_jprb    ! Overcast albedo for wvn > limit
  REAL(jprb), PARAMETER :: overcast_albedo2    = 0.6_jprb    ! Overcast albedo for wvn <= limit

  ! Parameters for Rayleigh scattering parameterization taken from Bucholtz 1995
  REAL(jprb), PARAMETER :: ray_ps = 1013.25_jprb, & ! Standard pressure (hPa) and temperature (K)
                           ray_ts = 288.15_jprb     !   at which parameterisation is computed

  REAL(jprb), PARAMETER :: ray_min_wlm = 0.2_jprb,         & ! Parameterisation valid for 0.2-4.0um
                           ray_max_wlm = 3.0_jprb,         & !   but we'll only allow Rayleigh up to 3um
                           ray_scs_wlm = 0.5_jprb,         & ! Wavelength limit: below 0.5um the
                           ray_scs_a1x = 3.01577E-28_jprb, & !   first set of parameters a1-d1 are used
                           ray_scs_a1b = 7.68246E-4_jprb,  & !   while above this a2-d2 are used.
                           ray_scs_b1  = -3.55212_jprb,    & !   a1x is for scattering cross-section
                           ray_scs_c1  = -1.35579_jprb,    & !   a1b is for scattering coefficient
                           ray_scs_d1  = -0.11563_jprb,    &
                           ray_scs_a2x = 4.01061E-28_jprb, &
                           ray_scs_a2b = 10.21675E-4_jprb, &
                           ray_scs_b2  = -3.99668_jprb,    &
                           ray_scs_c2  = -1.10298E-3_jprb, &
                           ray_scs_d2  = -2.71393E-2_jprb

  ! -------------------------------------------------------
  ! Polarisation definitions
  ! -------------------------------------------------------
  ! == pol_id +1
  !   1 average of vertical and horizontal
  !   2 nominal vertical at nadir, rotating with view angle
  !   3 nominal horizontal at nadir, rotating with view angle
  !   4 vertical
  !   5 horizontal
  !   6 + 45 minus -45 (3rd stokes vector)
  !   7 left circular - right circular (4th stokes vector)
  !   8 mixture of V and H according to angles pol_phi in rtcoef file:
  !       Tb = Tb_H * cos(pol_phi)^2 + Tb_V * sin(pol_phi)^2

  ! INTEGER(jpim), PARAMETER :: npolar_compute(7) = &
  !    (/ 2, 2, 2, 1, 1, 2, 4/)
  ! INTEGER(jpim), PARAMETER :: npolar_return(7) = &
  !    (/ 1, 1, 1, 1, 1, 2, 4/)

  ! pol_v and pol_h give proportion of v and h pol to use in emissivity calculation
  ! pol_s3 adds the 3rd/4th stokes vectors
  REAL(jprb), PARAMETER :: pol_v(3,8) = RESHAPE( &
      (/ 0.5_jprb, 0.0_jprb, 0.0_jprb, &
         0.0_jprb, 0.0_jprb, 1.0_jprb, &
         0.0_jprb, 1.0_jprb, 0.0_jprb, &
         1.0_jprb, 0.0_jprb, 0.0_jprb, &
         0.0_jprb, 0.0_jprb, 0.0_jprb, &
         0.0_jprb, 0.0_jprb, 0.0_jprb, &
         0.0_jprb, 0.0_jprb, 0.0_jprb, &
         1.0_jprb, 0.0_jprb, 0.0_jprb  /), (/3,8/) )
  REAL(jprb), PARAMETER :: pol_h(3,8) = RESHAPE( &
      (/ 0.5_jprb, 0.0_jprb, 0.0_jprb, &
         0.0_jprb, 1.0_jprb, 0.0_jprb, &
         0.0_jprb, 0.0_jprb, 1.0_jprb, &
         0.0_jprb, 0.0_jprb, 0.0_jprb, &
         1.0_jprb, 0.0_jprb, 0.0_jprb, &
         0.0_jprb, 0.0_jprb, 0.0_jprb, &
         0.0_jprb, 0.0_jprb, 0.0_jprb, &
         1.0_jprb, 0.0_jprb, 0.0_jprb  /), (/3,8/) )
  REAL(jprb), PARAMETER :: pol_s3(0:1,8) = RESHAPE( &
      (/ 0.0_jprb, 0.0_jprb, &
         0.0_jprb, 0.0_jprb, &
         0.0_jprb, 0.0_jprb, &
         0.0_jprb, 0.0_jprb, &
         0.0_jprb, 0.0_jprb, &
         1.0_jprb, 0.0_jprb, &
         0.0_jprb, 1.0_jprb, &
         0.0_jprb, 0.0_jprb  /), (/2,8/) )

  ! -------------------------------------------------------
  ! RTTOV-SCATT parameters
  ! -------------------------------------------------------
  ! Pressure of top level for line by line calculations (hPa)
  REAL(jprb), PARAMETER :: pressure_top = 0.004985_jprb
  ! Minimum single scattering albedo processed by rttov_scatt
  REAL(jprb), PARAMETER :: min_ssa = 1.0E-03_jprb
  ! Rain density (g.cm-3)
  REAL(jprb), PARAMETER :: rho_rain = 1.0_jprb
  ! Snow density (g.cm-3)
  REAL(jprb), PARAMETER :: rho_snow = 0.1_jprb
  ! Maximum optical depth
  REAL(jprb), PARAMETER :: max_scatt_optical_depth = 30._jprb
  ! (Deprecated) Flux to density conversions using RR = a * LWC^b, [RR]=mm/h, [LWC]=g/m^3
  ! See Geer et al. (2007) ECMWF Tech. Memo. 535, appendix B, for more info
  ! The variables store 1/a and 1/b for rain and snow and converted to kg/m2/s from mm/h:
  REAL(jprb), PARAMETER :: conv_rain(2) = (/ 3600.0_jprb/20.89_jprb , 1._jprb/1.15_jprb /)
  REAL(jprb), PARAMETER :: conv_sp  (2) = (/ 3600.0_jprb/29.51_jprb , 1._jprb/1.10_jprb /)
  ! Hydrometeor scattering optical properties filename stem
  CHARACTER(LEN=*),PARAMETER :: scattering_coef_filestem = 'hydrotable'
  ! NWP SAF hydrotable hydrometeor indices
  INTEGER(jpim), PARAMETER :: hydro_index_rain    = 1
  INTEGER(jpim), PARAMETER :: hydro_index_snow    = 2
  INTEGER(jpim), PARAMETER :: hydro_index_graupel = 3
  INTEGER(jpim), PARAMETER :: hydro_index_clw     = 4
  INTEGER(jpim), PARAMETER :: hydro_index_ciw     = 5
  ! Minimum allowable reflectivity output
  REAL(jprb), PARAMETER :: min_reflectivity = -999.0_jprb
  ! Output a tiny radiance as indicator that radar simulator radiances are not to be trusted
  REAL(jprb), PARAMETER :: min_radiance_radar = 1.0E-8_jprb
  ! Minimum allowable hydrometeor fraction in radar simulator
  REAL(jprb), PARAMETER :: min_cfrac_radar = 1.0E-4_jprb

  ! Flags to identify function in shared K/Adjoint routines
  INTEGER(jpim), PARAMETER :: adk_adjoint = 0
  INTEGER(jpim), PARAMETER :: adk_k       = 1

  ! -------------------------------------------------------
  ! Visible/IR scattering parameters
  ! -------------------------------------------------------
  ! Aerosols
  INTEGER(jpim), PARAMETER :: naer_opac = 13 ! Number of species in OPAC aerosol files
  INTEGER(jpim), PARAMETER :: naer_cams = 9  ! Number of species in CAMS aerosol files

  INTEGER(jpim), PARAMETER :: aer_id_opac = 1 ! ID for OPAC scaercoef files
  INTEGER(jpim), PARAMETER :: aer_id_cams = 2 ! ID for CAMS scaercoef files
  INTEGER(jpim), PARAMETER :: aer_id_user = 0 ! ID for user scaercoef files

  ! OPAC Aerosols
  CHARACTER(LEN=4), PARAMETER :: aer_opac_name(naer_opac) = &
    (/ 'inso', 'waso', 'soot', 'ssam', 'sscm', 'minm', 'miam', &
       'micm', 'mitr', 'suso', 'vola', 'vapo', 'asdu' /)

  ! OPAC MMR to number density conversion factors [g.m^-3]/[particle.cm^-3]
  ! NB the values stored in the scaercoef files are the reciprocal of these
  REAL(jprb), PARAMETER ::  aer_opac_confac(naer_opac) = &
    (/ 2.37E-05_jprb, 1.34E-09_jprb, 5.99E-11_jprb, 8.02E-07_jprb, 2.24E-04_jprb, 2.78E-08_jprb, &
       5.53E-06_jprb, 3.24E-04_jprb, 1.59E-05_jprb, 2.28E-08_jprb,   39.258_jprb,  13.431_jprb,  &
      1.473E-06_jprb /)

  ! CAMS Aerosols
  CHARACTER(LEN=4), PARAMETER :: aer_cams_name(naer_cams) = &
    (/ 'bcar', 'dus1', 'dus2', 'dus3', 'sulp', 'ssa1', 'ssa2', 'ssa3', 'omat' /)

  ! Constants for relative humidity calculation for aerosols
  REAL(jprb), PARAMETER :: e00 = 611.21_jprb
  REAL(jprb), PARAMETER :: t00 = 273.16_jprb
  REAL(jprb), PARAMETER :: ti  = t00 - 23.0_jprb

  ! Total number of cloud types (water and ice)
  INTEGER(jpim), PARAMETER :: ncldtyp = 6

  ! Liquid water clouds
  INTEGER(jpim), PARAMETER :: nclw_scheme = 2     ! Max valid value for clw_scheme
  INTEGER(jpim), PARAMETER :: clw_scheme_opac = 1
  INTEGER(jpim), PARAMETER :: clw_scheme_deff = 2
  INTEGER(jpim), PARAMETER :: nclwde_param = 1    ! Max valid value for clwde_param

  ! OPAC liquid water clouds
  INTEGER(jpim), PARAMETER :: nwcl_max = 5

  INTEGER(jpim), PARAMETER :: &
          wcl_id_stco       = 1, &
          wcl_id_stma       = 2, &
          wcl_id_cucc       = 3, &
          wcl_id_cucp       = 4, &
          wcl_id_cuma       = 5

  CHARACTER(LEN=4), PARAMETER :: wcl_opac_name(nwcl_max) = &
          (/ 'stco', 'stma', 'cucc', 'cucp', 'cuma' /)

  REAL(jprb), PARAMETER :: wcl_opac_deff(nwcl_max) = &  ! OPAC effective diameters (microns)
          (/ 14.66_jprb, 22.58_jprb, 11.54_jprb, 8.00_jprb, 25.36_jprb /)

  ! CLW Deff parameterisation (Martin et al, 1994)
  ! This minimum is applied by DWD, it doesn't come from the Martin et al paper
  ! but this is only applied to the parameterised Deff (cf dgmin_clw)
  REAL(jprb), PARAMETER :: martin_clwde_min = 10._jprb  ! microns
  ! k parameters for continental and maritime clouds [unitless]
  REAL(jprb), PARAMETER :: martin_k_land    = 0.67_jprb
  REAL(jprb), PARAMETER :: martin_k_sea     = 0.80_jprb
  ! Cloud particle number densities for continental and maritime clouds [m^-3]
  ! These values correspond approximately to the average ntot values in Fig 8c
  ! of Martin et al (300 and 100 cm^-3 respectively).
  REAL(jprb), PARAMETER :: martin_ntot_land = 3.0E8_jprb
  REAL(jprb), PARAMETER :: martin_ntot_sea  = 1.0E8_jprb

  ! Ice clouds
  INTEGER(jpim), PARAMETER :: nice_scheme = 3      ! Max valid value for ice_scheme
  INTEGER(jpim), PARAMETER :: ice_scheme_baum      = 1
  INTEGER(jpim), PARAMETER :: ice_scheme_baran2014 = 2
  INTEGER(jpim), PARAMETER :: ice_scheme_baran2018 = 3
  INTEGER(jpim), PARAMETER :: nicede_param = 4     ! Max valid value for icede_param
  INTEGER(jpim), PARAMETER :: baran_ngauss = 1000  ! Size of quadrature when computing Baran Leg. coefs

  ! Cloud overlap schemes
  INTEGER(jpim), PARAMETER :: ncloud_overlap = 2 ! Max valid value for cloud_overlap
  INTEGER(jpim), PARAMETER :: cloud_overlap_max_random = 1
  INTEGER(jpim), PARAMETER :: cloud_overlap_simple     = 2

  ! Minimum layer cloud fraction considered if grid_box_avg_cloud is true
  REAL(jprb),    PARAMETER :: cfrac_min = 1.E-6_jprb

  ! VIS/IR scattering models
  INTEGER(jpim), PARAMETER :: max_ir_scatt_model  = 2 ! Number of IR scattering models
  INTEGER(jpim), PARAMETER :: ir_scatt_dom        = 1
  INTEGER(jpim), PARAMETER :: ir_scatt_chou       = 2
  INTEGER(jpim), PARAMETER :: max_vis_scatt_model = 3 ! Number of solar scattering models
  INTEGER(jpim), PARAMETER :: vis_scatt_dom       = 1
  INTEGER(jpim), PARAMETER :: vis_scatt_single    = 2
  INTEGER(jpim), PARAMETER :: vis_scatt_mfasis    = 3

  ! Discrete Ordinates
  INTEGER(jpim), PARAMETER :: dom_min_nstr = 2     ! Min number of DOM streams allowed

  ! MFASIS
  INTEGER(jpim), PARAMETER :: mfasis_cld = 1      ! Particle type clouds (water and ice)
  INTEGER(jpim), PARAMETER :: mfasis_aer = 2      ! Particle type aersols

  INTEGER(jpim), PARAMETER :: mfasis_dim_kfourier = 0 ! Theta+ Fourier index (k) dimension type
  INTEGER(jpim), PARAMETER :: mfasis_dim_lfourier = 1 ! Theta- Fourier index (l) dimension type
  INTEGER(jpim), PARAMETER :: mfasis_dim_albedo   = 2 ! Albedo type dimension
  INTEGER(jpim), PARAMETER :: mfasis_dim_opdp     = 3 ! Optical depth type dimension
  INTEGER(jpim), PARAMETER :: mfasis_dim_effdia   = 4 ! Effective diameter type dimension
  INTEGER(jpim), PARAMETER :: mfasis_dim_scaangle = 5 ! Scattering angle type dimension
  ! Opposite angle conventions as DISORT/MFASIS paper
  REAL(jprb),    PARAMETER :: mfasis_maxzenangle    = 85._jprb  ! Max sat zenith angle
  REAL(jprb),    PARAMETER :: mfasis_maxsumzenangle = 150._jprb ! Max sum of zenith angles

  ! Phase function angle grids
  INTEGER(jpim), PARAMETER :: nphangle_lores = 208
  REAL(jprb), PARAMETER :: phangle_lores(nphangle_lores) = &
             (/   0.0_jprb,   0.1_jprb,   0.2_jprb,   0.3_jprb,   0.4_jprb,   0.5_jprb,   0.6_jprb, &
                  0.7_jprb,   0.8_jprb,   0.9_jprb,   1.0_jprb,   1.1_jprb,   1.2_jprb,   1.3_jprb, &
                  1.4_jprb,   1.5_jprb,   1.6_jprb,   1.7_jprb,   1.8_jprb,   1.9_jprb,   2.0_jprb, &
                  2.1_jprb,   2.2_jprb,   2.3_jprb,   2.4_jprb,   2.5_jprb,   2.6_jprb,   2.7_jprb, &
                  2.8_jprb,   2.9_jprb,   3.0_jprb,   4.0_jprb,   5.0_jprb,   6.0_jprb,   7.0_jprb, &
                  8.0_jprb,   9.0_jprb,  10.0_jprb,  11.0_jprb,  12.0_jprb,  13.0_jprb,  14.0_jprb, &
                 15.0_jprb,  16.0_jprb,  17.0_jprb,  18.0_jprb,  19.0_jprb,  20.0_jprb,  21.0_jprb, &
                 22.0_jprb,  23.0_jprb,  24.0_jprb,  25.0_jprb,  26.0_jprb,  27.0_jprb,  28.0_jprb, &
                 29.0_jprb,  30.0_jprb,  31.0_jprb,  32.0_jprb,  33.0_jprb,  34.0_jprb,  35.0_jprb, &
                 36.0_jprb,  37.0_jprb,  38.0_jprb,  39.0_jprb,  40.0_jprb,  41.0_jprb,  42.0_jprb, &
                 43.0_jprb,  44.0_jprb,  45.0_jprb,  46.0_jprb,  47.0_jprb,  48.0_jprb,  49.0_jprb, &
                 50.0_jprb,  51.0_jprb,  52.0_jprb,  53.0_jprb,  54.0_jprb,  55.0_jprb,  56.0_jprb, &
                 57.0_jprb,  58.0_jprb,  59.0_jprb,  60.0_jprb,  61.0_jprb,  62.0_jprb,  63.0_jprb, &
                 64.0_jprb,  65.0_jprb,  66.0_jprb,  67.0_jprb,  68.0_jprb,  69.0_jprb,  70.0_jprb, &
                 71.0_jprb,  72.0_jprb,  73.0_jprb,  74.0_jprb,  75.0_jprb,  76.0_jprb,  77.0_jprb, &
                 78.0_jprb,  79.0_jprb,  80.0_jprb,  81.0_jprb,  82.0_jprb,  83.0_jprb,  84.0_jprb, &
                 85.0_jprb,  86.0_jprb,  87.0_jprb,  88.0_jprb,  89.0_jprb,  90.0_jprb,  91.0_jprb, &
                 92.0_jprb,  93.0_jprb,  94.0_jprb,  95.0_jprb,  96.0_jprb,  97.0_jprb,  98.0_jprb, &
                 99.0_jprb, 100.0_jprb, 101.0_jprb, 102.0_jprb, 103.0_jprb, 104.0_jprb, 105.0_jprb, &
                106.0_jprb, 107.0_jprb, 108.0_jprb, 109.0_jprb, 110.0_jprb, 111.0_jprb, 112.0_jprb, &
                113.0_jprb, 114.0_jprb, 115.0_jprb, 116.0_jprb, 117.0_jprb, 118.0_jprb, 119.0_jprb, &
                120.0_jprb, 121.0_jprb, 122.0_jprb, 123.0_jprb, 124.0_jprb, 125.0_jprb, 126.0_jprb, &
                127.0_jprb, 128.0_jprb, 129.0_jprb, 130.0_jprb, 131.0_jprb, 132.0_jprb, 133.0_jprb, &
                134.0_jprb, 135.0_jprb, 136.0_jprb, 137.0_jprb, 138.0_jprb, 139.0_jprb, 140.0_jprb, &
                141.0_jprb, 142.0_jprb, 143.0_jprb, 144.0_jprb, 145.0_jprb, 146.0_jprb, 147.0_jprb, &
                148.0_jprb, 149.0_jprb, 150.0_jprb, 151.0_jprb, 152.0_jprb, 153.0_jprb, 154.0_jprb, &
                155.0_jprb, 156.0_jprb, 157.0_jprb, 158.0_jprb, 159.0_jprb, 160.0_jprb, 161.0_jprb, &
                162.0_jprb, 163.0_jprb, 164.0_jprb, 165.0_jprb, 166.0_jprb, 167.0_jprb, 168.0_jprb, &
                169.0_jprb, 170.0_jprb, 171.0_jprb, 172.0_jprb, 173.0_jprb, 174.0_jprb, 175.0_jprb, &
                176.0_jprb, 177.0_jprb, 178.0_jprb, 179.0_jprb, 180.0_jprb /)

  INTEGER(jpim), PARAMETER :: nphangle_hires = 498
  REAL(jprb), PARAMETER :: phangle_hires(nphangle_hires) = &
   (/ 0.00_jprb,   0.01_jprb,   0.02_jprb,   0.03_jprb,   0.04_jprb,   0.05_jprb,   0.06_jprb,   0.070_jprb, &
      0.08_jprb,   0.09_jprb,   0.10_jprb,   0.11_jprb,   0.12_jprb,   0.13_jprb,   0.14_jprb,   0.150_jprb, &
      0.16_jprb,   0.17_jprb,   0.18_jprb,   0.19_jprb,   0.20_jprb,   0.21_jprb,   0.22_jprb,   0.230_jprb, &
      0.24_jprb,   0.25_jprb,   0.26_jprb,   0.27_jprb,   0.28_jprb,   0.29_jprb,   0.30_jprb,   0.310_jprb, &
      0.32_jprb,   0.33_jprb,   0.34_jprb,   0.35_jprb,   0.36_jprb,   0.37_jprb,   0.38_jprb,   0.390_jprb, &
      0.40_jprb,   0.41_jprb,   0.42_jprb,   0.43_jprb,   0.44_jprb,   0.45_jprb,   0.46_jprb,   0.470_jprb, &
      0.48_jprb,   0.49_jprb,   0.50_jprb,   0.51_jprb,   0.52_jprb,   0.53_jprb,   0.54_jprb,   0.550_jprb, &
      0.56_jprb,   0.57_jprb,   0.58_jprb,   0.59_jprb,   0.60_jprb,   0.61_jprb,   0.62_jprb,   0.630_jprb, &
      0.64_jprb,   0.65_jprb,   0.66_jprb,   0.67_jprb,   0.68_jprb,   0.69_jprb,   0.70_jprb,   0.710_jprb, &
      0.72_jprb,   0.73_jprb,   0.74_jprb,   0.75_jprb,   0.76_jprb,   0.77_jprb,   0.78_jprb,   0.790_jprb, &
      0.80_jprb,   0.81_jprb,   0.82_jprb,   0.83_jprb,   0.84_jprb,   0.85_jprb,   0.86_jprb,   0.870_jprb, &
      0.88_jprb,   0.89_jprb,   0.90_jprb,   0.91_jprb,   0.92_jprb,   0.93_jprb,   0.94_jprb,   0.950_jprb, &
      0.96_jprb,   0.97_jprb,   0.98_jprb,   0.99_jprb,   1.00_jprb,   1.01_jprb,   1.02_jprb,   1.030_jprb, &
      1.04_jprb,   1.05_jprb,   1.06_jprb,   1.07_jprb,   1.08_jprb,   1.09_jprb,   1.10_jprb,   1.110_jprb, &
      1.12_jprb,   1.13_jprb,   1.14_jprb,   1.15_jprb,   1.16_jprb,   1.17_jprb,   1.18_jprb,   1.190_jprb, &
      1.20_jprb,   1.21_jprb,   1.22_jprb,   1.23_jprb,   1.24_jprb,   1.25_jprb,   1.26_jprb,   1.270_jprb, &
      1.28_jprb,   1.29_jprb,   1.30_jprb,   1.31_jprb,   1.32_jprb,   1.33_jprb,   1.34_jprb,   1.350_jprb, &
      1.36_jprb,   1.37_jprb,   1.38_jprb,   1.39_jprb,   1.40_jprb,   1.41_jprb,   1.42_jprb,   1.430_jprb, &
      1.44_jprb,   1.45_jprb,   1.46_jprb,   1.47_jprb,   1.48_jprb,   1.49_jprb,   1.50_jprb,   1.510_jprb, &
      1.52_jprb,   1.53_jprb,   1.54_jprb,   1.55_jprb,   1.56_jprb,   1.57_jprb,   1.58_jprb,   1.590_jprb, &
      1.60_jprb,   1.61_jprb,   1.62_jprb,   1.63_jprb,   1.64_jprb,   1.65_jprb,   1.66_jprb,   1.670_jprb, &
      1.68_jprb,   1.69_jprb,   1.70_jprb,   1.71_jprb,   1.72_jprb,   1.73_jprb,   1.74_jprb,   1.750_jprb, &
      1.76_jprb,   1.77_jprb,   1.78_jprb,   1.79_jprb,   1.80_jprb,   1.81_jprb,   1.82_jprb,   1.830_jprb, &
      1.84_jprb,   1.85_jprb,   1.86_jprb,   1.87_jprb,   1.88_jprb,   1.89_jprb,   1.90_jprb,   1.910_jprb, &
      1.92_jprb,   1.93_jprb,   1.94_jprb,   1.95_jprb,   1.96_jprb,   1.97_jprb,   1.98_jprb,   1.990_jprb, &
      2.00_jprb,   2.05_jprb,   2.10_jprb,   2.15_jprb,   2.20_jprb,   2.25_jprb,   2.30_jprb,   2.350_jprb, &
      2.40_jprb,   2.45_jprb,   2.50_jprb,   2.55_jprb,   2.60_jprb,   2.65_jprb,   2.70_jprb,   2.750_jprb, &
      2.80_jprb,   2.85_jprb,   2.90_jprb,   2.95_jprb,   3.00_jprb,   3.05_jprb,   3.10_jprb,   3.150_jprb, &
      3.20_jprb,   3.25_jprb,   3.30_jprb,   3.35_jprb,   3.40_jprb,   3.45_jprb,   3.50_jprb,   3.550_jprb, &
      3.60_jprb,   3.65_jprb,   3.70_jprb,   3.75_jprb,   3.80_jprb,   3.85_jprb,   3.90_jprb,   3.950_jprb, &
      4.00_jprb,   4.05_jprb,   4.10_jprb,   4.15_jprb,   4.20_jprb,   4.25_jprb,   4.30_jprb,   4.350_jprb, &
      4.40_jprb,   4.45_jprb,   4.50_jprb,   4.55_jprb,   4.60_jprb,   4.65_jprb,   4.70_jprb,   4.750_jprb, &
      4.80_jprb,   4.85_jprb,   4.90_jprb,   4.95_jprb,   5.00_jprb,   5.10_jprb,   5.20_jprb,   5.300_jprb, &
      5.40_jprb,   5.50_jprb,   5.60_jprb,   5.70_jprb,   5.80_jprb,   5.90_jprb,   6.00_jprb,   6.100_jprb, &
      6.20_jprb,   6.30_jprb,   6.40_jprb,   6.50_jprb,   6.60_jprb,   6.70_jprb,   6.80_jprb,   6.900_jprb, &
      7.00_jprb,   7.10_jprb,   7.20_jprb,   7.30_jprb,   7.40_jprb,   7.50_jprb,   7.60_jprb,   7.700_jprb, &
      7.80_jprb,   7.90_jprb,   8.00_jprb,   8.10_jprb,   8.20_jprb,   8.30_jprb,   8.40_jprb,   8.500_jprb, &
      8.60_jprb,   8.70_jprb,   8.80_jprb,   8.90_jprb,   9.00_jprb,   9.10_jprb,   9.20_jprb,   9.300_jprb, &
      9.40_jprb,   9.50_jprb,   9.60_jprb,   9.70_jprb,   9.80_jprb,   9.90_jprb,  10.00_jprb,  10.500_jprb, &
     11.00_jprb,  11.50_jprb,  12.00_jprb,  12.50_jprb,  13.00_jprb,  13.50_jprb,  14.00_jprb,  14.500_jprb, &
     15.00_jprb,  16.00_jprb,  17.00_jprb,  18.00_jprb,  19.00_jprb,  20.00_jprb,  21.00_jprb,  22.000_jprb, &
     23.00_jprb,  24.00_jprb,  25.00_jprb,  26.00_jprb,  27.00_jprb,  28.00_jprb,  29.00_jprb,  30.000_jprb, &
     31.00_jprb,  32.00_jprb,  33.00_jprb,  34.00_jprb,  35.00_jprb,  36.00_jprb,  37.00_jprb,  38.000_jprb, &
     39.00_jprb,  40.00_jprb,  41.00_jprb,  42.00_jprb,  43.00_jprb,  44.00_jprb,  45.00_jprb,  46.000_jprb, &
     47.00_jprb,  48.00_jprb,  49.00_jprb,  50.00_jprb,  51.00_jprb,  52.00_jprb,  53.00_jprb,  54.000_jprb, &
     55.00_jprb,  56.00_jprb,  57.00_jprb,  58.00_jprb,  59.00_jprb,  60.00_jprb,  61.00_jprb,  62.000_jprb, &
     63.00_jprb,  64.00_jprb,  65.00_jprb,  66.00_jprb,  67.00_jprb,  68.00_jprb,  69.00_jprb,  70.000_jprb, &
     71.00_jprb,  72.00_jprb,  73.00_jprb,  74.00_jprb,  75.00_jprb,  76.00_jprb,  77.00_jprb,  78.000_jprb, &
     79.00_jprb,  80.00_jprb,  81.00_jprb,  82.00_jprb,  83.00_jprb,  84.00_jprb,  85.00_jprb,  86.000_jprb, &
     87.00_jprb,  88.00_jprb,  89.00_jprb,  90.00_jprb,  91.00_jprb,  92.00_jprb,  93.00_jprb,  94.000_jprb, &
     95.00_jprb,  96.00_jprb,  97.00_jprb,  98.00_jprb,  99.00_jprb, 100.00_jprb, 101.00_jprb, 102.000_jprb, &
    103.00_jprb, 104.00_jprb, 105.00_jprb, 106.00_jprb, 107.00_jprb, 108.00_jprb, 109.00_jprb, 110.000_jprb, &
    111.00_jprb, 112.00_jprb, 113.00_jprb, 114.00_jprb, 115.00_jprb, 116.00_jprb, 117.00_jprb, 118.000_jprb, &
    119.00_jprb, 120.00_jprb, 121.00_jprb, 122.00_jprb, 123.00_jprb, 124.00_jprb, 125.00_jprb, 126.000_jprb, &
    127.00_jprb, 128.00_jprb, 129.00_jprb, 130.00_jprb, 131.00_jprb, 132.00_jprb, 133.00_jprb, 134.000_jprb, &
    135.00_jprb, 136.00_jprb, 137.00_jprb, 138.00_jprb, 139.00_jprb, 140.00_jprb, 141.00_jprb, 142.000_jprb, &
    143.00_jprb, 144.00_jprb, 145.00_jprb, 146.00_jprb, 147.00_jprb, 148.00_jprb, 149.00_jprb, 150.000_jprb, &
    151.00_jprb, 152.00_jprb, 153.00_jprb, 154.00_jprb, 155.00_jprb, 156.00_jprb, 157.00_jprb, 158.000_jprb, &
    159.00_jprb, 160.00_jprb, 161.00_jprb, 162.00_jprb, 163.00_jprb, 164.00_jprb, 165.00_jprb, 166.000_jprb, &
    167.00_jprb, 168.00_jprb, 169.00_jprb, 170.00_jprb, 171.00_jprb, 172.00_jprb, 173.00_jprb, 174.000_jprb, &
    175.00_jprb, 176.00_jprb, 176.25_jprb, 176.50_jprb, 176.75_jprb, 177.00_jprb, 177.25_jprb, 177.500_jprb, &
    177.75_jprb, 178.00_jprb, 178.25_jprb, 178.50_jprb, 178.75_jprb, 179.00_jprb, 179.25_jprb, 179.500_jprb, &
    179.75_jprb, 180.00_jprb /)

  ! -------------------------------------------------------
  ! Wavenumber intervals for RTTOV9 predictors
  ! -------------------------------------------------------
  REAL(jprb), PARAMETER :: rttov9_wv0690_50 =  690.50_jprb, &
                           rttov9_wv1050_00 = 1050.00_jprb, &
                           rttov9_wv1095_25 = 1095.25_jprb, &
                           rttov9_wv1100_25 = 1100.25_jprb, &
                           rttov9_wv1350_25 = 1350.25_jprb, &
                           rttov9_wv1750_25 = 1750.25_jprb, &
                           rttov9_wv1900_25 = 1900.25_jprb, &
                           rttov9_wv1995_00 = 1995.00_jprb, &
                           rttov9_wv2000_00 = 2000.00_jprb, &
                           rttov9_wv2250_00 = 2250.00_jprb, &
                           rttov9_wv2295_25 = 2295.25_jprb, &
                           rttov9_wv2360_00 = 2360.00_jprb, &
                           rttov9_wv2380_25 = 2380.25_jprb, &
                           rttov9_wv2660_25 = 2660.25_jprb, &
                           rttov9_wv2760_25 = 2760.25_jprb

END MODULE rttov_const
