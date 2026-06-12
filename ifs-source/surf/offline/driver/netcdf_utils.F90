MODULE NETCDF_UTILS
! (C) Copyright 1995- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
USE PARKIND1  ,ONLY : JPIM,JPRB 
USE NETCDF
USE MPL_MODULE

IMPLICIT NONE
SAVE

TYPE TYPE_NCDF_VAR
  CHARACTER(LEN=20)  :: CNAME   ! VARIABLE NAME 
  CHARACTER(LEN=50)  :: CLNAME  ! ATTRIBUTE: LONG NAME
  CHARACTER(LEN=50)  :: CUNITS  ! ATTRIBUTE: UNITS
  CHARACTER(LEN=50)  :: CCOORD  ! ATTRIBUTE: coordinates
  CHARACTER(LEN=200)  :: COMMENT  ! ATTRIBUTE: coordinates
  INTEGER(KIND=JPIM) :: IACCUR  ! ACCURACY
  INTEGER(KIND=JPIM) :: IDLEVEL  ! Deflation LEVEL
  LOGICAL            :: LFILL     ! TRUE IF IT CAN BE MASKED
END TYPE TYPE_NCDF_VAR

TYPE TYPE_NCDF_DIMSID
  INTEGER(KIND=JPIM) :: NLATID  ! Latitude ID
  INTEGER(KIND=JPIM) :: NLONID  ! Longitude ID
  INTEGER(KIND=JPIM) :: NLEVSID ! soil levels ID
  INTEGER(KIND=JPIM) :: NTILID  ! tile ID
  INTEGER(KIND=JPIM) :: NVTID   ! vegetation type ID
  INTEGER(KIND=JPIM) :: NMONID  ! month id
  INTEGER(KIND=JPIM) :: NTIMID  ! Time id
  INTEGER(KIND=JPIM) :: NLEVSNID ! snow levels ID 
END TYPE TYPE_NCDF_DIMSID 
  
!       ----------------------------------------------------------------
!*    ** *NETCDF_UTILS* - I/O NETCDF UTILITIES 
!       ----------------------------------------------------------------

CONTAINS

SUBROUTINE INIT_NCDF_VAR(TVOUT,CNAME,IACCUR_IN,CCOORD_IN)
USE YOMLOG1S , ONLY : NACCUR   ,NDIMCDF, NDLEVEL


IMPLICIT NONE
TYPE(TYPE_NCDF_VAR),INTENT(OUT) :: TVOUT 
CHARACTER(LEN=*),INTENT(IN)     :: CNAME   ! VARIABLE NAME 
!! OPTIONAL 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)     :: IACCUR_IN  ! ACCURACY
CHARACTER(LEN=*),OPTIONAL,INTENT(IN)     :: CCOORD_IN  ! coordinates 

! local variables 
CHARACTER(LEN=50)  :: CLNAME  ! ATTRIBUTE: LONG NAME
CHARACTER(LEN=50)  :: CUNITS  ! ATTRIBUTE: UNITS
CHARACTER(LEN=50)  :: CCOORD  ! ATTRIBUTE: coordinates
CHARACTER(LEN=200)  :: COMMENT  ! ATTRIBUTE: comment
INTEGER(KIND=JPIM) :: IACCUR  ! ACCURACY
INTEGER(KIND=JPIM) :: IDLEVEL  ! Deflation LEVEL
LOGICAL            :: LFILL     ! TRUE IF IT CAN BE MASKED

!* -- set defaults here 
 CLNAME=CNAME
 CUNITS="-"
 CCOORD="NONE"
 COMMENT="NONE"
IF (NACCUR == 1) THEN
  IACCUR=NF90_FLOAT
ELSEIF (NACCUR == 2) THEN
  IACCUR=NF90_DOUBLE
ENDIF
IDLEVEL=NDLEVEL
LFILL = .TRUE.

!* -- change for particular cases.... 
SELECT CASE(CNAME)
  !* grid definitions
  CASE('lat')
    CLNAME='latitude'
    CUNITS='degrees_north'
    LFILL = .FALSE.
    IACCUR=NF90_DOUBLE
  CASE('lon')
    CLNAME='longitude'
    CUNITS='degrees_east'
    LFILL = .FALSE.
    IACCUR=NF90_DOUBLE
  CASE('nlevs')
    CLNAME='soil level '
    CUNITS='level'
    LFILL = .FALSE.
    IACCUR=NF90_INT
  CASE('nlevsn')
    CLNAME='snow level '
    CUNITS='level'
    LFILL = .FALSE.
    IACCUR=NF90_INT
  CASE('tile')
    CLNAME='tile number'
    CUNITS='level'
    IACCUR=NF90_INT
    LFILL = .FALSE.
    COMMENT='1=water;2=sea ice;3=interception;4=low veg;5=shaded snow;6=high veg;7=exposed snow;8=baresoil;9=lakes'
  CASE('vtype')
    CLNAME='vegetation type number'
    CUNITS='level'
    IACCUR=NF90_INT
    LFILL = .FALSE.
    COMMENT='1=low veg;2=high veg'
  ! temporal definitions 
  CASE('time')
    CLNAME='time'
    IACCUR=NF90_DOUBLE
    LFILL = .FALSE.
  CASE('timestp')
    CLNAME='model time step'
    IACCUR=NF90_INT
    LFILL = .FALSE.
  CASE('month')
    CLNAME='months'
    CUNITS='level'
    IACCUR=NF90_INT
    LFILL = .FALSE.
  CASE('Mask')
    CLNAME='Catchment mask'
    CUNITS='-'
    IACCUR=NF90_DOUBLE
    LFILL = .FALSE.
    COMMENT='1=run point;0=do not run point'
  CASE('z0m')
    CLNAME='momentum roughness'
    CUNITS='m'
  CASE('lz0h')
    CLNAME='log roughness length for heat'
    CUNITS='m'
  CASE('landsea')
    CLNAME='land sea mask'
    CUNITS='-'
  CASE('geopot')
    CLNAME='surface geopotential'
    CUNITS='m2 m-2'
  CASE('cvl')
    CLNAME='low vegetation coverage'
    CUNITS='-'
  CASE('cvh')
    CLNAME='high vegetation coverage'
    CUNITS='-'
  CASE('cu')
    CLNAME='urban coverage'
    CUNITS='-'
  CASE('tvl')
    CLNAME='low vegetation type'
    CUNITS='-'
  CASE('Ctype')
    CLNAME='C3/C4 type of photosynthetic pathway'
    CUNITS='-'
  CASE('tvh')
    CLNAME='high vegetation type'
    CUNITS='-'
  CASE('sotype')
    CLNAME='soil type'
    CUNITS='-'
  CASE('sdor')
    CLNAME='standard deviation of orography'
    CUNITS='m'
  CASE('sst')
    CLNAME='sea surface temperature'
    CUNITS='K'
  CASE('seaice')
    CLNAME='sea ice fraction'
    CUNITS='-'
  CASE('CanopInt')
    CLNAME='canopy interception depth'
    CUNITS='kg m-2'
  CASE('SWE')
    CLNAME='snow water equivalent'
    CUNITS='kg m-2'
  CASE('SWEML')
    CLNAME='snow water equivalent'
    CUNITS='kg m-2'
  CASE('SnowT')
    CLNAME='snow temperature'
    CUNITS='K'
  CASE('SnowTML')
    CLNAME='snow temperature'
    CUNITS='K'
  CASE('SAlbedo')
    CLNAME='Snow albedo'
    CUNITS='-'
  CASE('snowdens')
    CLNAME='Snow density'
    CUNITS='kg m-3'
  CASE('snowdensML')
    CLNAME='Snow density'
    CUNITS='kg m-3'
  CASE('slw')
    CLNAME='snow liquid water'
    CUNITS='kg m-2'
  CASE('slwML')
    CLNAME='snow liquid water'
    CUNITS='kg m-2'
  CASE('AvgSurfT','SkinT')
    CLNAME='skin temperature'
    CUNITS='K'
  CASE('SkinC')
    CLNAME='skin thermal conductivity'
    CUNITS='W m-2 k-1'
  CASE('TLICE')
    CLNAME='lake ice temperature'
    CUNITS='K'
  CASE('TLMNW')
    CLNAME='lake mean water temperature'
    CUNITS='K'
  CASE('TLWML')
    CLNAME='lake mixed layer temperature'
    CUNITS='K'
  CASE('TLBOT')
    CLNAME='lake bottom temperature'
    CUNITS='K'
  CASE('TLSF')
    CLNAME='lake temperature shape factor'
    CUNITS='-'
  CASE('HLICE')
    CLNAME='lake ice thickness'
    CUNITS='m'
  CASE('HLML')
    CLNAME='lake mixed layer thickness'
    CUNITS='m'
  CASE('LDEPTH')
    CLNAME='lake depth'
    CUNITS='m'
  CASE('CLAKE')
    CLNAME='lake cover'
    CUNITS='-'
  CASE('CLAKEF')
    CLNAME='lake and floodplain cover'
    CUNITS='-'
  CASE('SoilTemp')
    CLNAME='soil temperature'
    CUNITS='K'
  CASE('SoilMoist')
    CLNAME='soil moisture'
    CUNITS='kg m-2'
  CASE('icetemp')
    CLNAME='sea ice temperature'
    CUNITS='K'
  CASE('Malbedo')
    CLNAME='climatological surface albedo'
    CUNITS='-'
  CASE('Mlail')
    CLNAME='climatological low vegetation lai'
    CUNITS='m2 m-2'
  CASE('Mlaih')
    CLNAME='climatological high vegetation lai'
    CUNITS='m2 m-2'
  CASE('r0vt')
    CLNAME='reference respiration'
    CUNITS='kg m-2 s-1'
  CASE('fwet')
    CLNAME='climatological wetland fraction'
    CUNITS='-'
  CASE('Anday','An')
    CLNAME='net CO2 assimilation'
    CUNITS='kg m-2 s-1'
  CASE('Anfm')
    CLNAME='maximum leaf assimilation per vegtype'
    CUNITS='kg m-2 s-1'
  CASE('Rstr')
    CLNAME='respiration of above ground structural biomass'
    CUNITS='kg m-2 '
  CASE('Rstr2')
    CLNAME='respiration of below ground structural biomass'
    CUNITS='kg m-2 '
  CASE('Blast')
    CLNAME='active biomass of previous day'
    CUNITS='kg m-2 '
  CASE('laivt','lai')
    CLNAME='leaf area index'
    CUNITS='m2 m-2 '
  CASE('Biomstr')
    CLNAME='above ground structural biomass'
    CUNITS='kg m-2 '
  CASE('Biomstr2')
    CLNAME='below ground structural biomass'
    CUNITS='kg m-2 '
  CASE('Momu')
    CLNAME='momentum flux East-West'
    CUNITS='kg m-1 s-2'
  CASE('Momv')
    CLNAME='momentum flux North-South'
    CUNITS='kg m-1 s-2'
  CASE('SWnet')
    CLNAME='net shortwave radiation'
    CUNITS='W m-2'
  CASE('LWnet')
    CLNAME='net longwave radiation'
    CUNITS='W m-2'
  CASE('Qle','LEvt')
    CLNAME='latent heat flux'
    CUNITS='W m-2'
  CASE('Qh')
    CLNAME='sensible heat flux'
    CUNITS='W m-2'
  CASE('Qg')
    CLNAME='soil heat flux'
    CUNITS='W m-2'
  CASE('Qgsn')
    CLNAME='snow basal heat flux'
    CUNITS='W m-2'
  CASE('Qf')
    CLNAME='soil fusion heat flux'
    CUNITS='W m-2'
  CASE('Qfsn')
    CLNAME='snow fusion heat flux'
    CUNITS='W m-2'
  CASE('DelSoilHeat')
    CLNAME='soil heat content change'
    CUNITS='J m-2'
  CASE('DelColdCont')
    CLNAME='snow heat content change'
    CUNITS='J m-2'
  CASE('LWup')
    CLNAME='upward longwave flux'
    CUNITS='W m-2'
  CASE('SWdown')
    CLNAME='downward shortwave flux'
    CUNITS='W m-2'
  CASE('LWdown')
    CLNAME='downward longwave flux'
    CUNITS='W m-2'
  CASE('Snowf')
    CLNAME='snowfall rate'
    CUNITS='kg m-2 s-1'
  CASE('Rainf')
    CLNAME='rainfall rate'
    CUNITS='kg m-2 s-1'
  CASE('Intercept')
    CLNAME='Interception of rainfall rate'
    CUNITS='kg m-2 s-1'
  CASE('Evap')
    CLNAME='total evapotranspiration'
    CUNITS='kg m-2 s-1'
  CASE('Qs')
    CLNAME='surface runoff'
    CUNITS='kg m-2 s-1'
  CASE('Qsb')
    CLNAME='subsurface runoff'
    CUNITS='kg m-2 s-1'
  CASE('DelSoilMoist')
    CLNAME='soil moisture content change'
    CUNITS='kg m-2'
  CASE('DelSWE')
    CLNAME='snow water content change'
    CUNITS='kg m-2'
  CASE('DelIntercept')
    CLNAME='Insterception content change'
    CUNITS='kg m-2'
  CASE('Qsm')
    CLNAME='snowmelt'
    CUNITS='kg m-2 s*1'
  CASE('T2m')
    CLNAME='2 meters temperature'
    CUNITS='K'
  CASE('RH2m')
    CLNAME='2 meters relative humidity'
    CUNITS='-'
  CASE('D2m')
    CLNAME='2 meters dewpoint temperature'
    CUNITS='K'
  CASE('VegT')
    CLNAME='vegetation temperature'
    CUNITS='K'
  CASE('BaresoilT')
    CLNAME='temperature for bare soil'
    CUNITS='K'
  CASE('RadT')
    CLNAME='surface radiative temperature'
    CUNITS='K'
  CASE('Albedo')
    CLNAME='surface albedo'
    CUNITS='K'
  CASE('ECanop')
    CLNAME='interception evaporation'
    CUNITS='kg m-2 s-1'
  CASE('TVeg')
    CLNAME='vegetation transpiration'
    CUNITS='kg m-2 s-1'
  CASE('Conds')
    CLNAME='excess condensation'
    CUNITS='kg m-2 s-1'
  CASE('ESoil')
    CLNAME='bare soil evaporation'
    CUNITS='kg m-2 s-1'
  CASE('RootMoist')
    CLNAME='root zone soil moisture'
    CUNITS='kg m-2 s-1'
  CASE('SubSnow')
    CLNAME='snow sublimation'
    CUNITS='kg m-2 s-1'
  CASE('PotEvapI')
    CLNAME='Potential evaporation irrigation'
    CUNITS='kg m-2 s-1'
  CASE('PotEvapWA')
    CLNAME='Potential evaporation open water using air temp'
    CUNITS='kg m-2 s-1'
  CASE('PotEvapU')
    CLNAME='Potential evaporation unstressed'
    CUNITS='kg m-2 s-1'
  CASE('EWater')
    CLNAME='open water evaporation'
    CUNITS='kg m-2 s-1'
  CASE('SnowFrac')
    CLNAME='snow cover fraction'
    CUNITS='-'
  CASE('IceFrac')
    CLNAME='ice covered gridbox fraction'
    CUNITS='-'
  CASE('Fdepth')
    CLNAME='Frozen soil depth'
    CUNITS='m'
  CASE('Tdepth')
    CLNAME='Depth to soil thaw'
    CUNITS='m'
  CASE('SnowDepth')
    CLNAME='snow depth'
    CUNITS='m'
  CASE('SoilThick')
    CLNAME='soil layer thicknesses'
    CUNITS='m'
  CASE('SoilFC')
    CLNAME='soil field capacity'
    CUNITS='m3 m-3'
  CASE('SoilWilt')
    CLNAME='soil wilting point'
    CUNITS='m3 m-3'
  CASE('SoilSat')
    CLNAME='soil porosity'
    CUNITS='m3 m-3'
  CASE('RootDist')
    CLNAME='root fraction'
    CUNITS='-' 
  CASE('Ag')
    CLNAME='Gross CO2 assimilation'
    CUNITS='kg m-2 s-1' 
  CASE('Rd')
    CLNAME='Dark respiration'
    CUNITS='kg m-2 s-1' 
  CASE('Rsoil_str')
    CLNAME='Respiration from soil and structural biomass'
    CUNITS='kg m-2 s-1'  
  CASE('Reco')
    CLNAME='Ecosystem Respiration'
    CUNITS='kg m-2 s-1'  
  CASE('CO2flux')
    CLNAME='Net CO2 flux'
    CUNITS='kg m-2 s-1'
  CASE('RnoQ10')
    CLNAME='Respiration without Q10-dependency'
    CUNITS='kg m-2 s-1'
  CASE('CH4flux')
    CLNAME='Net CH4 flux'
    CUNITS='kg m-2 s-1'
  CASE('biomass')
    CLNAME='Active biomass'
    CUNITS='kg m-2'
  CASE('Bloss')
    CLNAME='Active biomass loss'
    CUNITS='kg m-2'
  CASE('Bgain')
    CLNAME='Active biomass gain'
    CUNITS='kg m-2'
  CASE('rc')
    CLNAME='canopy resistance'
    CUNITS='s m-1'
  CASE('ra')
    CLNAME='aerodynamic resistance'
    CUNITS='s m-1'
  CASE('f2vt')
    CLNAME='stress function F2'
    CUNITS='s m-1'
  CASE('dsvt')
    CLNAME='Specific humidity deficit at canopy level'
    CUNITS='kg kg-1'
  CASE('dmaxvt')
    CLNAME='Maximum specific humidity deficit'
    CUNITS='kg kg-1'
  CASE('tifr')
    CLNAME='tile fraction'
    CUNITS='-'
  CASE('vtfr')
    CLNAME='vegetation type fraction'
    CUNITS='-'
  CASE('intfr')
    CLNAME='Interception fraction'
    CUNITS='-'
  CASE('Tair')
    CLNAME='Air temperature (forcing)'
    CUNITS='K'
  CASE('Qair')
    CLNAME='Specific humidity (forcing)'
    CUNITS='kg/kg'
  CASE('CO2air')
    CLNAME='Atmospheric CO2 (forcing)'
    CUNITS='ppm'
  CASE('Uwind')
    CLNAME='U wind component (forcing)'
    CUNITS='m/s'
  CASE('Vwind')
    CLNAME='V wind component (forcing)'
    CUNITS='m/s'
  CASE('psurf')
    CLNAME='surface pressure (forcing)'
    CUNITS='Pa'
  CASE('x')
    CLNAME='gaussian reduced index'
    CUNITS='-'
  CASE('fldfrc')
    CLNAME='Flooded area fraction'
    CUNITS='-'
  CASE DEFAULT
    WRITE(*,*) CNAME, ' Not defined in INIT_NCDF_VAR'
    CALL ABOR1('INIT_NCDF_VAR')

END SELECT

IF ( PRESENT(IACCUR_IN) ) THEN
  IACCUR = IACCUR_IN
ENDIF
IF ( PRESENT(CCOORD_IN) ) THEN
  CCOORD = CCOORD_IN
ENDIF

!* -- INITIALIZE A NCDF_VAR TYPE 
TVOUT%CNAME   =  TRIM(CNAME)
TVOUT%CLNAME  =  TRIM(CLNAME)
TVOUT%CUNITS  =  TRIM(CUNITS)
TVOUT%IACCUR  =  IACCUR
TVOUT%IDLEVEL =  IDLEVEL
TVOUT%LFILL   =  LFILL
TVOUT%CCOORD  =  CCOORD
TVOUT%COMMENT  =  COMMENT
END SUBROUTINE INIT_NCDF_VAR

SUBROUTINE NC_DEF_VAR(NPOS,TVAR,DIMIDS,VARID)
USE YOMLUN1S , ONLY : RMISS,IMISS,ZMISS,NCTYPE,DMISS
IMPLICIT NONE 

!* -- define a netcdf variable 

INTEGER(KIND=JPIM),INTENT(IN)  :: NPOS
TYPE(TYPE_NCDF_VAR),INTENT(IN) :: TVAR
INTEGER(KIND=JPIM),INTENT(IN)  :: DIMIDS(:)
INTEGER(KIND=JPIM),INTENT(OUT) :: VARID

!* -- LOCAL 
LOGICAL :: SHUFFLE

SHUFFLE = .TRUE.
IF ( TVAR%IDLEVEL == 0 ) SHUFFLE = .FALSE.
! WRITE(*,*),'Creating:',NPOS,TRIM(TVAR%CNAME)
! print*,DIMIDS
IF ( NCTYPE == NF90_NETCDF4  .OR. NCTYPE == NF90_CLASSIC_MODEL+NF90_NETCDF4  ) THEN
  CALL NCERROR( NF90_DEF_VAR(NPOS,TRIM(TVAR%CNAME),TVAR%IACCUR,DIMIDS,VARID,&
                DEFLATE_LEVEL=TVAR%IDLEVEL,SHUFFLE=SHUFFLE,cache_preemption=0,cache_nelems=1,cache_size=10),&
                "Creating variable "//TRIM(TVAR%CNAME))
ELSEIF ( NCTYPE == NF90_64BIT_OFFSET .OR. NCTYPE == NF90_CLOBBER .OR. NCTYPE == NF90_SHARE  ) THEN
  CALL NCERROR( NF90_DEF_VAR(NPOS,TRIM(TVAR%CNAME),TVAR%IACCUR,DIMIDS,VARID),&
                "Creating variable "//TRIM(TVAR%CNAME))
ENDIF

CALL NCERROR( NF90_PUT_ATT(NPOS,VARID,'long_name',TRIM(TVAR%CLNAME)) )

CALL NCERROR( NF90_PUT_ATT(NPOS,VARID,'units',TRIM(TVAR%CUNITS)) )

IF ( TVAR%CCOORD .NE. 'NONE' ) THEN
  CALL NCERROR( NF90_PUT_ATT(NPOS,VARID,'coordinates',TRIM(TVAR%CCOORD)) )
ENDIF
IF ( TVAR%COMMENT .NE. 'NONE' ) THEN
  CALL NCERROR( NF90_PUT_ATT(NPOS,VARID,'comment',TRIM(TVAR%COMMENT)) )
ENDIF

IF ( TVAR%LFILL ) THEN
  IF (TVAR%IACCUR == NF90_DOUBLE ) THEN 
    CALL NCERROR( NF90_PUT_ATT(NPOS,VARID,'_FillValue',DMISS))
  ELSEIF (TVAR%IACCUR == NF90_INT ) THEN 
    CALL NCERROR( NF90_PUT_ATT(NPOS,VARID,'_FillValue',IMISS))
  ELSE
    CALL NCERROR( NF90_PUT_ATT(NPOS,VARID,'_FillValue',ZMISS))
  ENDIF
ENDIF

END SUBROUTINE NC_DEF_VAR

SUBROUTINE DEF_NC_FILE(NPOS,TDIMSID)

USE YOMDPHY  , ONLY : NCSS, NLAT,NLON,NTILES,NVHILO,NCLIMDAT,NCSNEC
USE YOMLOG1S , ONLY : CMODID   ,CVERID   ,NDIMCDF
USE YOMRIP   , ONLY : NINDAT   ,NSSSSS
USE YOMLUN1S , ONLY : NPOSRES  ,CTUNITS,NPOSCLM
USE YOEPHY   , ONLY : LESNML
IMPLICIT NONE 

!! General definition of netcdf files: dimensions and global attributes
INTEGER(KIND=JPIM),INTENT(IN) :: NPOS ! netcdf unit
TYPE(TYPE_NCDF_DIMSID),INTENT(OUT) :: TDIMSID

!* -- OTHERS
TYPE(TYPE_NCDF_VAR)          :: YD_VARINFO
INTEGER(KIND=JPIM)           :: VARID

ASSOCIATE( NLATID=>TDIMSID%NLATID, NLONID=>TDIMSID%NLONID, NLEVSID=>TDIMSID%NLEVSID, &
         & NTILID=>TDIMSID%NTILID, NVTID=>TDIMSID%NVTID  , NMONID=>TDIMSID%NMONID  , &
         & NTIMID=>TDIMSID%NTIMID, NLEVSNID=>TDIMSID%NLEVSNID )

!* -- 0 Initialize DIMIDS to -1
NLATID=-1
NLONID=-1
NLEVSID=-1
NTILID=-1
NVTID=-1
NMONID=-1
NTIMID=1

!* -- 1. create dimensions
  
IF (NDIMCDF == 2 ) THEN
  CALL NCERROR( NF90_DEF_DIM(NPOS,'lat',NLAT,NLATID))
  CALL NCERROR( NF90_DEF_DIM(NPOS,'lon',NLON,NLONID))
ELSE
  CALL NCERROR( NF90_DEF_DIM(NPOS,'x',NLON,NLONID))
  NLATID = NLONID
ENDIF
CALL NCERROR( NF90_DEF_DIM(NPOS,'nlevs',NCSS,NLEVSID))
CALL NCERROR( NF90_DEF_DIM(NPOS,'tile',NTILES,NTILID))
CALL NCERROR( NF90_DEF_DIM(NPOS,'vtype',NVHILO,NVTID))
IF ( NPOS == NPOSRES ) THEN
  CALL NCERROR( NF90_DEF_DIM(NPOS,'month',NCLIMDAT,NMONID))
  CALL NCERROR( NF90_DEF_DIM(NPOS,'time',1,NTIMID))
ELSEIF ( NPOS /= NPOSCLM ) THEN
  CALL NCERROR( NF90_DEF_DIM(NPOS,'time',NF90_UNLIMITED,NTIMID))
ENDIF
IF (LESNML) THEN
  CALL NCERROR( NF90_DEF_DIM(NPOS,'nlevsn',NCSNEC,NLEVSNID))
ENDIF

!* -- 3 global attributes

CALL NCERROR( NF90_PUT_ATT(NPOS, NF90_GLOBAL, 'modelID', CMODID) )
CALL NCERROR( NF90_PUT_ATT(NPOS, NF90_GLOBAL, 'versionID', CVERID) )
CALL NCERROR( NF90_PUT_ATT(NPOS, NF90_GLOBAL, 'start_day', NINDAT) )
CALL NCERROR( NF90_PUT_ATT(NPOS, NF90_GLOBAL, 'start_hour', NSSSSS) )
CALL NCERROR( NF90_PUT_ATT(NPOS, NF90_GLOBAL, 'SurfSgn_convention', 'Mathematical') )

!* -- 4. create variables 

!* =====================================
!*       DIMENSIONS 

!* -- latitude
CALL INIT_NCDF_VAR(YD_VARINFO,'lat')
CALL NC_DEF_VAR(NPOS,YD_VARINFO,(/NLATID/),VARID)

!* -- longitude
CALL INIT_NCDF_VAR(YD_VARINFO,'lon')
CALL NC_DEF_VAR(NPOS,YD_VARINFO,(/NLONID/),VARID)

!* -- soil levels
CALL INIT_NCDF_VAR(YD_VARINFO,'nlevs')
CALL NC_DEF_VAR(NPOS,YD_VARINFO,(/NLEVSID/),VARID)

!* -- Tiles
CALL INIT_NCDF_VAR(YD_VARINFO,'tile')
CALL NC_DEF_VAR(NPOS,YD_VARINFO,(/NTILID/),VARID)

!* -- Vegetation types
CALL INIT_NCDF_VAR(YD_VARINFO,'vtype')
CALL NC_DEF_VAR(NPOS,YD_VARINFO,(/NVTID/),VARID)

IF ( NPOS /= NPOSCLM ) THEN 
  !* -- Time
  CALL INIT_NCDF_VAR(YD_VARINFO,'time')
  YD_VARINFO%CUNITS=CTUNITS
  CALL NC_DEF_VAR(NPOS,YD_VARINFO,(/NTIMID/),VARID)
  
  !* -- timestep
  CALL INIT_NCDF_VAR(YD_VARINFO,'timestp')
  CALL NC_DEF_VAR(NPOS,YD_VARINFO,(/NTIMID/),VARID)
ENDIF 

IF ( NPOS == NPOSRES ) THEN
  !* -- months 
  CALL INIT_NCDF_VAR(YD_VARINFO,'month')
  CALL NC_DEF_VAR(NPOS,YD_VARINFO,(/NMONID/),VARID)
ENDIF

IF (LESNML) THEN
  CALL INIT_NCDF_VAR(YD_VARINFO,'nlevsn')
  CALL NC_DEF_VAR(NPOS,YD_VARINFO,(/NLEVSNID/),VARID)
ENDIF

! !* --  Grid masks
! IF(NDIMCDF == 1) THEN
!   CALL INIT_NCDF_VAR(YD_VARINFO,'Grid_Mask_x')
!   CALL NC_DEF_VAR(NPOS,YD_VARINFO,(/NLONID/),VARID)
! 
!   CALL INIT_NCDF_VAR(YD_VARINFO,'Grid_Mask_y')
!   CALL NC_DEF_VAR(NPOS,YD_VARINFO,(/NLONID/),VARID)
! ENDIF

END ASSOCIATE
END SUBROUTINE DEF_NC_FILE

SUBROUTINE NCERROR(STATUS,STRING,LABORT)
IMPLICIT NONE 

! HANDLE NETCDF ERRORS 
! USAGE: CALL NCERROR( NF90_(),STRING,LABORT)
! STATUS = RETURN OF NETCDF CALL FUNCTIONS
! STRING - EXTRA STRIN TO PRINT - OPTIONAL
! LABORT - IF TRUE (DEFAULT) ABORTS THE PROGRAM 

INTEGER(KIND=JPIM),INTENT(IN) :: STATUS
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: STRING
LOGICAL,INTENT(IN),OPTIONAL     :: LABORT

IF ( STATUS /= 0 ) THEN
  PRINT*, TRIM(NF90_STRERROR(STATUS))
  IF( PRESENT(STRING) ) WRITE(*,*) TRIM(STRING)
  IF ( PRESENT(LABORT) ) THEN
    IF (LABORT) CALL ABOR1('NCERROR:')
  ENDIF
  CALL ABOR1('NCERROR:')

ENDIF 

END SUBROUTINE NCERROR

END MODULE NETCDF_UTILS 
