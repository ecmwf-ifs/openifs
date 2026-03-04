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

MODULE YOMSATS

USE PARKIND1   , ONLY : JPIM, JPRB
USE RTTOV_CONST, ONLY : NINST

IMPLICIT NONE

SAVE

!*

!     /YOMSATS/  J. PAILLEUX 90/11/30

! Rose Munro (Dec 1997) Add pointers for multiple satellite series, 
!                       id's and subtypes 
! Rose Munro (Jun 1998) Add QMAXRT, O3MAXRT and O3MINRT for new RT model
! Dick Dee   2004-03-19 Support for variational sat bias correction
! L.v.Bremen (Apr 2003) Add MSG and MODIS in character list CHGEOSATN for AMV 
! M Szyndel  (March 2005) Add MTSAT, FY-2 and update MSG radiances
!     TABLES DEFINING ALL THE LOGISTICS FOR TOVS RADIANCE COMPUTATION.
! A.Benedetti 2005-Mar-22 Add MODIS in SATEM sensor list 
! N Bormann (Feb 2006)  Add METOP
! P Bauer/G Kelly/N Bormann (Feb 2006) Add AMSRE, SSMIS, TMI
! V Guidard (March 2007) Add IASI (after A Collard)
! A Collard (Oct 2006) Add IASI
! N Bormann (Nov 2007) Adjust CHGEOSATN
! N Bormann (Dec 2007) Update for RTTOV-9
! A Inness  (Aug 2009) Add varno to SATREO3GRP_TABLE
! J Munoz Sabater (May 2009) Introduce SMOS
! V Guidard (Feb 2009) Add LCO2_DIAG_AIRS
! W Bell (Feb 2009) Add Coriolis-Windsat
! Q Lu      (Mar 2009) Add FY3A
! N Bormann (Mar 2011) Add ntoplevels to satrgrp structure
! N Bormann (Mar 2011) Add ATMS
! N Bormann (Jul 2011) Remove unnecessary RT-subtype and related variables
! T McNally (Oct 2011) Add CRIS
! M Kazumori(May 2013) Add GCOM-W1 AMSR2
! R Eresmaa (Feb 2013) Add AVHRR and capability to a collocated instrument
! S Massart (Nov 2012) Add SCIAMACHY data
! LF Meunier, P Chambon (Nov 2013) Add SAPHIR
! S Migliorini (Nov 2013) Add logical mask for IASI all-sky channels
! P Lean    (Apr 2014) Add GPM/GMI
! LF Meunier (May 2014) O3 climatological profil defined as a parameter
! C Payan (Sep 2014) *GEOSATN variables renamed *AMV_VC, base 40_op2.03
!   definitions revision cleaner and more logical?
! H Lawrence (Feb 2015) Add index for AMSU-A instruments to SATGRP_TABLE, 
!                       for clearsky observation error definition
! A Geer (July 2015) Pre-OOPS cleaning
! C.Lupu (Feb 2017) : RTTOV 12.1 updates
! C Burrows (Jul 2019) Add GIIRS
! R Eresmaa (Apr 2020) Add HIRAS and IKFS2

!   COMMENTS:
!   --------

!     RADPRE  : PRESSURE LEVELS USED IN THE RAD. TRANSF. MODEL
!     VDZCOR  : SATEM THICKNESS CORRELATION MATRIX
!     OERRDZ  : SATEM THICKNESS OBSERVATION ERRORS
!     OERRPWCREL: SATEM PWC OBSERVATION ERROR, relative to SPWC
!     ORERPWC : SATEM PWC ROUND-OFF OBSERVATION ERROR
!     NJPPF   : Number of profiles that RTTOV is processing in one go
!     O3CLRT_N: OZONE climatology profile (numer of levels).
!     O3CLRT_P: OZONE climatology profile (pressure).
!     O3CLRT  : OZONE climatology profile.
!     NVATOVINDX: Index to TOVS data in control variable
!     NVATOVOFF:  Offset to TOVS data in control variable
!     NVATOVLEN:  Local length TOVS data in control variable
!     NVATOVLENP:  Length TOVS data in control variable all procs
!     TVTSBGE_LAND : TOVS Ts background error, land
!     TVTSBGE_SEA  : TOVS Ts background error, sea
!     TVTSBGE_ICE  : TOVS Ts background error, ice
!     AIRSCHAN: used AIRS channels
!     AIRSCHAN_BACK: AIRS channels
!     IASICHAN: used IASI channels
!     IASICHAN_BACK: IASI channels
!     IASICHAN: IASI channels
!     LIASICHAN: IASI channels all-sky mask
!     LALLSKY_IASI_BACK: IASI channels all-sky mask
!     CRISCHAN: CRIS channels
!     CRISCHAN_BACK: CRIS channels index
!     HIRASCHAN: HIRAS channels
!     HIRASCHAN_BACK: HIRAS channels index
!     IKFS2CHAN: IKFS2 channels
!     IKFS2CHAN_BACK: IKFS2 channels index

!  VARIABLES AND PARAMETERS USED AS DIMENSIONS IN THESE ARRAYS ARE
!  RESPECTIVELY DEFINED IN COMDECKS "YOMDIMO" AND "PARDIMO".

REAL(KIND=JPRB),ALLOCATABLE:: RADPRE(:)

INTEGER(KIND=JPIM),ALLOCATABLE:: NVATOVINDX(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NVATOVOFF(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NVATOVLENP(:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NVALCHAN(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NVALCHAN_MTS(:,:)
INTEGER(KIND=JPIM) :: NVATOVLEN
INTEGER(KIND=JPIM) :: NVATOVLENG

INTEGER(KIND=JPIM) :: NJPPF = 0

! Define SATOB methods numbering
INTEGER(KIND=JPIM), PARAMETER :: MXMETHOD         = 5
INTEGER(KIND=JPIM), PARAMETER :: MMETHOD_WVMW     = 5
CHARACTER(LEN=8), PARAMETER :: CHMETHOD(0:MXMETHOD) = &
 & (/'WVCL    ','IR      ','VIS     ','WVMIX   ','COM_CHAN','WVMW    '/)  

! Define AMV/SATOB Virtual-Constellation (:,0)=>polar (:,1)=>geo
! Based on COMMON CODE TABLES TO BINARY AND ALPHANUMERIC CODES
! http://www.wmo.int/pages/prog/www/WMOCodes/WMO306_vI2/LatestVERSION/WMO306_vI2_CommonTable_en.pdf
!CP probably still not enough flexible => MXSERIES approach as for radiances?
INTEGER(KIND=JPIM), PARAMETER :: MXAMV_VC         = 11
CHARACTER(LEN=9), PARAMETER :: CHAMV_VC(0:MXAMV_VC,0:1) = RESHAPE(&
  & (/'METOP    ','JAPAN-POL','NOAA     ','METEOR   ','INDIA-POL','FY-POLAR ',&
  &   'METOP    ','MODIS    ','POL-undef','Reserved ','dual-MTOP','LEO-GEO  ',&
  &   'METEOSAT ','HIMAWARI ','GOES     ','ELECTRO  ','INSAT    ','FY-GEO   ',&
  &   'METEOSAT ','GOES     ','COMS     ','Reserved ','Undefined','Undefined'/),&
  &   (/MXAMV_VC+1,2/))

! Diagnostic logical: T - RTSETUP has been called
LOGICAL  :: LRTSETUP


! Define sensor numbering for instruments that are not in the RTTOV framework
! MSG_HR is using inst_id=29, attributed in rttov_const to viirs
INTEGER(KIND=JPIM), PARAMETER :: MSENSOR_MSG_HR   =29
INTEGER(KIND=JPIM), PARAMETER :: NSENSOR_MERIS    =174
INTEGER(KIND=JPIM), PARAMETER :: NSENSOR_MOPITT   =090
INTEGER(KIND=JPIM), PARAMETER :: MSENSOR_SCIAMACHY=175
INTEGER(KIND=JPIM), PARAMETER :: MSENSOR_GOSAT    =516
INTEGER(KIND=JPIM), PARAMETER :: MSENSOR_OCO2     =432

! Define channel numbering offsets, by sensor.
INTEGER(KIND=JPIM), PARAMETER :: MSENSOR_CHANOFFL(0:NINST-1)= &
 & (/0,20,24,27,42,0,0,0,8,0,0,0,0,0,0,42,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
 &   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
 &   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
 &   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
 &   0,0,0,0/)  

! Define channel counts by sensor (AIRS and IASI are channel subsets): 
INTEGER(KIND=JPIM), PARAMETER :: M__NUMAIRSCHANS = 324
INTEGER(KIND=JPIM), PARAMETER :: M__MAXAIRSCHANS = 2378
INTEGER(KIND=JPIM)  :: M__NUMIASICHANS = -999
INTEGER(KIND=JPIM)  :: M__NUMCRISCHANS = -999
INTEGER(KIND=JPIM)  :: M__NUMHIRASCHANS = -999
INTEGER(KIND=JPIM)  :: M__NUMIKFS2CHANS = -999
INTEGER(KIND=JPIM)  :: M__NUMGIIRSCHANS = -999
INTEGER(KIND=JPIM), PARAMETER :: M__MAXIASICHANS = 8461
INTEGER(KIND=JPIM), PARAMETER :: M__MAXCRISCHANS = 2211
INTEGER(KIND=JPIM), PARAMETER :: M__MAXHIRASCHANS = 2275
INTEGER(KIND=JPIM), PARAMETER :: M__MAXIKFS2CHANS = 2701
INTEGER(KIND=JPIM), PARAMETER :: M__MAXGIIRSCHANS = 1650
INTEGER(KIND=JPIM)  :: MSENSOR_CHANCOUNT(0:NINST-1)= &
 & (/20,4,3,15,5,3,7,8,8,9,24,M__NUMAIRSCHANS,0,0,0,5,-999, &
 &   12,18,22,2,8,4,0,4,4,0,-999,0,8,16, &
 &   0,0,0,6,0,0,0,0,0,4,5,20,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,14, &
 &   0,0,0,0,0,0,0,13,13,15,0,0,29,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-999,0,0, &
 &   -999,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
 &   17,0/)

! Define usable AIRS channels
INTEGER(KIND=JPIM), PARAMETER :: AIRSCHAN(1:M__NUMAIRSCHANS) = &
 & (/1,6,7,10,11,15,16,17,20,21,22,24,27,28,30,36,39,40,42,51,52,54,55,56,&
 & 59,62,63,68,69,71,72,73,74,75,76,77,78,79,80,82,83,84,86,92,93,98,99,&
 & 101,104,105,108,110,111,113,116,117,123,124,128,129,138,139,144,145,&
 & 150,151,156,157,159,162,165,168,169,170,172,173,174,175,177,179,180,&
 & 182,185,186,190,192,193,198,201,204,207,210,213,215,216,218,221,224,&
 & 226,227,232,239,248,250,251,252,253,256,257,261,262,267,272,295,299,&
 & 300,305,308,309,310,318,321,325,333,338,355,362,375,453,475,484,497,&
 & 528,587,672,787,791,843,870,914,950,1003,1012,1019,1024,1030,1038,1048,&
 & 1069,1079,1082,1083,1088,1090,1092,1095,1104,1111,1115,1116,1119,1120,&
 & 1123,1130,1138,1142,1178,1199,1206,1221,1237,1252,1260,1263,1266,1278,&
 & 1285,1290,1301,1304,1329,1371,1382,1400,1401,1402,1403,1415,1424,1449,&
 & 1455,1466,1471,1477,1479,1488,1500,1519,1520,1538,1545,1565,1574,1583,&
 & 1593,1614,1627,1636,1644,1652,1669,1674,1681,1694,1708,1717,1723,1740,&
 & 1748,1751,1756,1763,1766,1771,1777,1780,1783,1794,1800,1803,1806,1812,&
 & 1826,1843,1852,1865,1866,1867,1868,1869,1872,1873,1875,1876,1877,1881,&
 & 1882,1883,1884,1897,1901,1911,1917,1918,1921,1923,1924,1928,1937,1938,&
 & 1939,1941,1946,1947,1948,1958,1971,1973,1988,1995,2084,2085,2097,2098,&
 & 2099,2100,2101,2103,2104,2106,2107,2108,2109,2110,2111,2112,2113,2114,&
 & 2115,2116,2117,2118,2119,2120,2121,2122,2123,2128,2134,2141,2145,2149,&
 & 2153,2164,2189,2197,2209,2226,2234,2280,2318,2321,2325,2328,2333,2339,&
 & 2348,2353,2355,2357,2363,2370,2371,2377/)  

INTEGER(KIND=JPIM) :: AIRSCHAN_BACK(M__MAXAIRSCHANS)

INTEGER(KIND=JPIM), ALLOCATABLE :: IASICHAN(:)

LOGICAL, ALLOCATABLE :: LIASICHAN(:)

INTEGER(KIND=JPIM) :: IASICHAN_BACK(M__MAXIASICHANS)

LOGICAL            :: LALLSKY_IASI_BACK(M__MAXIASICHANS)

INTEGER(KIND=JPIM), ALLOCATABLE :: CRISCHAN(:)

INTEGER(KIND=JPIM) :: CRISCHAN_BACK(M__MAXCRISCHANS)

INTEGER(KIND=JPIM), ALLOCATABLE :: HIRASCHAN(:)

INTEGER(KIND=JPIM) :: HIRASCHAN_BACK(M__MAXHIRASCHANS)

INTEGER(KIND=JPIM), ALLOCATABLE :: IKFS2CHAN(:)

INTEGER(KIND=JPIM) :: IKFS2CHAN_BACK(M__MAXIKFS2CHANS)

INTEGER(KIND=JPIM), ALLOCATABLE :: GIIRSCHAN(:)

INTEGER(KIND=JPIM) :: GIIRSCHAN_BACK(M__MAXGIIRSCHANS)

! Background error for TOVS Ts sink control variable
REAL(KIND=JPRB) :: TVTSBGE_LAND(0:NINST-1)
REAL(KIND=JPRB) :: TVTSBGE_SEA(0:NINST-1)
REAL(KIND=JPRB) :: TVTSBGE_ICE(0:NINST-1)

! Define table for administration of sat-ID, sensor and codetype...
TYPE SATGRP_T
INTEGER(KIND=JPIM)         :: BUFRID                   ! BUFR Satellite Identifier
INTEGER(KIND=JPIM)         :: CODETYPE                 ! Obs codetype (as in ODB)
INTEGER(KIND=JPIM)         :: SENSOR                   ! Sensor (as defined above)
INTEGER(KIND=JPIM)         :: SENSOR_SAFE              ! Sensor ID within range 0:NINST-1 for use with e.g. VarQC arrays
INTEGER(KIND=JPIM)         :: SATGROUP                 ! Satellite data group number
INTEGER(KIND=JPIM)         :: RETRTYPE                 ! Retrieval type
INTEGER(KIND=JPIM)         :: RTID                     ! SatID used by RTTOV
INTEGER(KIND=JPIM)         :: RTSERIES                 ! Series (as defined above)
INTEGER(KIND=JPIM)         :: RTCOEF_POS               ! Pointer to RT-coefficients
LOGICAL                    :: CLD_RTCALC_SCREEN        ! Perform cloudy RT computations in screening
LOGICAL                    :: CLD_RTCALC_ASSIM         ! Perform cloudy RT computations in the assimilation
INTEGER(KIND=JPIM)         :: NRTCHANNELS              ! No. of channels, this sensor
INTEGER(KIND=JPIM),POINTER :: RTCHAN_LIST(:)=>NULL()   ! Chans that RTTOV can compute
INTEGER(KIND=JPIM),POINTER :: VARBC_IX(:)=>NULL()      ! Bias parameter table index
INTEGER(KIND=JPIM)         :: NTOPLEVELS               ! Number of RTTOV levels above IFS model top
INTEGER(KIND=JPIM)         :: RTCOEF_POS_COLLOC        ! RT coefficient index of a collocated sensor
INTEGER(KIND=JPIM)         :: IOBSERR_IX               ! index for observation error definition 
CHARACTER(LEN=32)          :: SNAME                    ! name of sub-table
INTEGER(KIND=JPIM)         :: TREAT_2DGOM_ID           ! ID indicating 2d-GOM treatment (see gom2d_support)
END TYPE SATGRP_T
TYPE(SATGRP_T), ALLOCATABLE :: SATGRP_TABLE(:)         ! Sat data group admin table
INTEGER(KIND=JPIM) :: NSATGRP                          ! Size of SATGRP_TABLE

! Define table for administration of SATOB sat-ID, codetype...
TYPE SATOBGRP_T
INTEGER(KIND=JPIM)         :: BUFRID                   ! BUFR Satellite Identifier
INTEGER(KIND=JPIM)         :: CODETYPE                 ! Obs codetype (as in ODB)
INTEGER(KIND=JPIM)         :: METHOD                   ! Computational method
INTEGER(KIND=JPIM)         :: SATGROUP                 ! Satellite data group number
INTEGER(KIND=JPIM)         :: IPCOR                    ! Index of correction model
LOGICAL                    :: ASYM_FGCHECK             ! Asymmetric FG-check on/off
LOGICAL                    :: LOW_SPD_CHECK            ! Check to exclude low wind speeds on/off
CHARACTER(LEN=32)          :: OBS_OPER                 ! Obs operator method 
REAL(KIND=JPRB)            :: SIGMA                    ! Sigma used for obs operator [Pa]
REAL(KIND=JPRB)            :: ZREJMOD_FACTOR           ! Factor for fg check tuning
REAL(KIND=JPRB)            :: OERR_FACTOR              ! Factor for obs error tuning

INTEGER(KIND=JPIM)         :: RTID                     ! SatID used by RTTOV
INTEGER(KIND=JPIM)         :: RTSERIES                 ! Series (as defined above)
INTEGER(KIND=JPIM)         :: RTCOEF_POS               ! Pointer to RT-coefficients
END TYPE SATOBGRP_T
TYPE(SATOBGRP_T), ALLOCATABLE :: SATOBGRP_TABLE(:)     ! SATOB data group table
INTEGER(KIND=JPIM) :: NSATOBGRP                        ! Size of SATGRP_TABLE

! Define table for administration of sat-ID, sensor and codetype...
TYPE SATREO3GRP_T
INTEGER(KIND=JPIM)         :: BUFRID                   ! BUFR Satellite Identifier
INTEGER(KIND=JPIM)         :: SENSOR                   ! Sensor (as defined above)
INTEGER(KIND=JPIM)         :: VARNO                    ! Variable number
INTEGER(KIND=JPIM)         :: SATGROUP                 ! Satellite data group number
INTEGER(KIND=JPIM)         :: RETRTYPE                 ! Retrieval type
INTEGER(KIND=JPIM)         :: PRODTYPE                 ! Product type
END TYPE SATREO3GRP_T
TYPE(SATREO3GRP_T), ALLOCATABLE :: SATREO3GRP_TABLE(:)         ! Sat data group admin table
INTEGER(KIND=JPIM) :: NSATREO3GRP                          ! Size of SATREO3GRP_TABLE

! Define a O3 climatological profil that will be used if no Ozone field is
! availlable (MF only)

INTEGER(KIND=JPIM), PARAMETER:: O3CLRT_N = 44

! Note: Well below the ground level, a fake pressure level is added in 
!       order simplify the interpolation algo.
REAL(KIND=JPRB), PARAMETER:: O3CLRT_P(1:O3CLRT_N+1)=(/&
 & 0.0050_JPRB  , &
 & 0.1000_JPRB   ,0.2900_JPRB   ,0.6900_JPRB   ,1.4200_JPRB, &
 & 2.6110_JPRB   ,4.4070_JPRB   ,6.9500_JPRB   ,10.3700_JPRB, &
 & 14.8100_JPRB  ,20.4000_JPRB  ,27.2600_JPRB  ,35.5100_JPRB, &
 & 45.2900_JPRB  ,56.7300_JPRB  ,69.9700_JPRB  ,85.1800_JPRB, &
 & 102.0500_JPRB ,122.0400_JPRB ,143.8400_JPRB ,167.9500_JPRB, &
 & 194.3600_JPRB ,222.9400_JPRB ,253.7100_JPRB ,286.6000_JPRB, &
 & 321.5000_JPRB ,358.2800_JPRB ,396.8100_JPRB ,436.9500_JPRB, &
 & 478.5400_JPRB ,521.4600_JPRB ,565.5400_JPRB ,610.6000_JPRB, &
 & 656.4300_JPRB ,702.7300_JPRB ,749.1200_JPRB ,795.0900_JPRB, &
 & 839.9500_JPRB ,882.8000_JPRB ,922.4600_JPRB ,957.4400_JPRB, &
 & 985.8800_JPRB ,1005.4300_JPRB,1013.2500_JPRB, &
 & 99999.0000_JPRB /) * 100._JPRB ! (Pa)

REAL(KIND=JPRB), PARAMETER:: O3CLRT(1:O3CLRT_N)=(/&
 & 0.969339E-05_JPRB,                                                         &
 & 0.969339E-05_JPRB, 0.100043E-04_JPRB, 0.101194E-04_JPRB, 0.101751E-04_JPRB,&
 & 0.102181E-04_JPRB, 0.102463E-04_JPRB, 0.102127E-04_JPRB, 0.102454E-04_JPRB,&
 & 0.101497E-04_JPRB, 0.935746E-05_JPRB, 0.809846E-05_JPRB, 0.672055E-05_JPRB,&
 & 0.519183E-05_JPRB, 0.372400E-05_JPRB, 0.258175E-05_JPRB, 0.172204E-05_JPRB,&
 & 0.119041E-05_JPRB, 0.845131E-06_JPRB, 0.649756E-06_JPRB, 0.526723E-06_JPRB,&
 & 0.412770E-06_JPRB, 0.303262E-06_JPRB, 0.210869E-06_JPRB, 0.156281E-06_JPRB,&
 & 0.123483E-06_JPRB, 0.107569E-06_JPRB, 0.100696E-06_JPRB, 0.959905E-07_JPRB,&
 & 0.916880E-07_JPRB, 0.891324E-07_JPRB, 0.846609E-07_JPRB, 0.811542E-07_JPRB,&
 & 0.778113E-07_JPRB, 0.757306E-07_JPRB, 0.711752E-07_JPRB, 0.663465E-07_JPRB,&
 & 0.616040E-07_JPRB, 0.567638E-07_JPRB, 0.521222E-07_JPRB, 0.479273E-07_JPRB,&
 & 0.443969E-07_JPRB, 0.419678E-07_JPRB, 0.409984E-07_JPRB /)  

!     ------------------------------------------------------------------

END MODULE YOMSATS

