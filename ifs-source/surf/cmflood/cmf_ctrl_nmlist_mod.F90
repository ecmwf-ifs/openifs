MODULE CMF_CTRL_NMLIST_MOD
!==========================================================
!* PURPOSE: read CaMa-Flood Model configulations from namelist ("input_flood.nam" as default)
!
!* CONTAINS:
! -- CMF_CONFIG_NAMELIST  : read namelist for CaMa-Flood 
! -- CMF_CONFIG_CHECK     : check config conflict
!
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  29Jul 2019
!             Adapted mostly from CMF v362 CONTROL0.F90
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
! shared variables in module
USE PARKIND1,                ONLY: JPIM, JPRB, JPRM
USE YOS_CMF_INPUT,           ONLY: LOGNAM
IMPLICIT NONE
CONTAINS
!####################################################################
! -- CMF_CONFIG_NAMELIST  : read namelist for CaMa-Flood 
! -- CMF_CONFIG_CHECK     : check config conflict
!
!####################################################################
SUBROUTINE CMF_CONFIG_NMLIST
USE YOS_CMF_INPUT,      ONLY: TMPNAM,   NSETFILE,   CSETFILE
! run version
USE YOS_CMF_INPUT,      ONLY: LADPSTP,  LFPLAIN,  LKINE,    LFLDOUT,  LPTHOUT,  LDAMOUT,  &
                            & LROSPLIT, LGDWDLY,  LSLPMIX,  LMEANSL,  LSEALEV,  LOUTPUT,  &
                            & LRESTART, LSTOONLY, LGRIDMAP, LLEAPYR,  LMAPEND,  LBITSAFE, &
                            & LSLOPEMOUTH,LWEVAP,LWEVAPFIX, LWEXTRACTRIV
! dimention & time
USE YOS_CMF_INPUT,      ONLY: CDIMINFO, DT,       NX,NY,    NLFP,     NXIN,NYIN,    INPN, &
                            & IFRQ_INP, DTIN,     WEST,EAST,NORTH,SOUTH
! parameters
USE YOS_CMF_INPUT,      ONLY: PMANRIV,  PMANFLD,  PDSTMTH,  PMINSLP,  PGRV, PCADP, &
                            & IMIS, RMIS, DMIS,   CSUFBIN,  CSUFVEC,  CSUFPTH,  CSUFCDF
USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID
IMPLICIT NONE
!
NAMELIST/NRUNVER/  LADPSTP,  LFPLAIN,  LKINE,    LFLDOUT,  LPTHOUT,  LDAMOUT, &
                   LROSPLIT, LGDWDLY,  LSLPMIX,  LMEANSL,  LSEALEV,  LOUTPUT, &
                   LRESTART, LSTOONLY, LGRIDMAP, LLEAPYR,  LMAPEND,  LBITSAFE,&
                   LSLOPEMOUTH,LWEVAP,LWEVAPFIX, LWEXTRACTRIV

NAMELIST/NDIMTIME/ CDIMINFO, DT, IFRQ_INP

NAMELIST/NPARAM/   PMANRIV, PMANFLD, PGRV,    PDSTMTH, PCADP,   PMINSLP, &
                   IMIS, RMIS, DMIS, CSUFBIN, CSUFVEC, CSUFPTH, CSUFCDF
!================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!--------------------"

! *** 0. SET INPUT UNIT AND OPEN FILE 
NSETFILE=INQUIRE_FID()               !!  for namelist
OPEN(NSETFILE,FILE=CSETFILE,STATUS="OLD")
WRITE(LOGNAM,*) "CMF::CONFIG_NMLIST: namelist opened: ", TRIM(CSETFILE), NSETFILE 

!============================
!*** 1. basic simulation run version

! * defaults
LADPSTP  = .TRUE.            !! true: use adaptive time step
LFPLAIN  = .TRUE.            !! true: consider floodplain (false: only river channel)
LKINE    = .FALSE.           !! true: use kinematic wave
LFLDOUT  = .TRUE.            !! true: floodplain flow (high-water channel flow) active
LPTHOUT  = .FALSE.           !! true: activate bifurcation scheme
LDAMOUT  = .FALSE.           !! true: activate dam operation (under development)

LROSPLIT = .FALSE.           !! true: input if surface (Qs) and sub-surface (Qsb) runoff
LWEVAP   = .FALSE.           !! true: input evaporation to extract from river 
LWEVAPFIX= .FALSE.           !! true: water balance closure extracting water from evap when available
LGDWDLY  = .FALSE.           !! true: Activate ground water reservoir and delay
LSLPMIX  = .FALSE.           !! true: activate mixed kinematic and local inertia based on slope
LWEXTRACTRIV=.FALSE.         !! true: also extract water from rivers 

LMEANSL  = .FALSE.           !! true : boundary condition for mean sea level
LSEALEV  = .FALSE.           !! true : boundary condition for variable sea level

LRESTART = .FALSE.           !! true: initial condition from restart file
LSTOONLY = .FALSE.           !! true: storage only restart (mainly for data assimilation)
LOUTPUT  = .TRUE.            !! true: use standard output (to file)

LGRIDMAP = .TRUE.            !! true: for standard XY gridded 2D map
LLEAPYR  = .TRUE.            !! true: neglect leap year (Feb29 skipped)
LMAPEND  = .FALSE.           !! true: for map data endian conversion
LBITSAFE = .FALSE.           !! true: for Bit Identical simulation (avoid OSM ATOMIC)

LSLOPEMOUTH = .FALSE.        !! true: prescribe water level slope == elevation slope on river month

!* change
REWIND(NSETFILE)
READ(NSETFILE,NML=NRUNVER)

WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "=== NAMELIST, NRUNVER ==="
WRITE(LOGNAM,*) "LADPSTP ",  LADPSTP
WRITE(LOGNAM,*) "LFPLAIN ",  LFPLAIN
WRITE(LOGNAM,*) "LKINE   ",  LKINE
WRITE(LOGNAM,*) "LFLDOUT ",  LFLDOUT
WRITE(LOGNAM,*) "LPTHOUT ",  LPTHOUT
WRITE(LOGNAM,*) "LDAMOUT ",  LDAMOUT
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "LROSPLIT ", LROSPLIT
WRITE(LOGNAM,*) "LWEVAP   ", LWEVAP
WRITE(LOGNAM,*) "LWEVAPFIX", LWEVAPFIX
WRITE(LOGNAM,*) "LWEXTRACTRIV", LWEXTRACTRIV
WRITE(LOGNAM,*) "LGDWDLY  ", LGDWDLY
WRITE(LOGNAM,*) "LSLPMIX  ", LSLPMIX
WRITE(LOGNAM,*) "LMEANSL: ", LSEALEV
WRITE(LOGNAM,*) "LSEALEV: ", LSEALEV
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "LRESTART ", LRESTART
WRITE(LOGNAM,*) "LSTOONLY ", LSTOONLY
WRITE(LOGNAM,*) "LOUTPUT  ", LOUTPUT
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "LGRIDMAP ", LGRIDMAP
WRITE(LOGNAM,*) "LLEAPYR  ", LLEAPYR
WRITE(LOGNAM,*) "LMAPEND  ", LMAPEND
WRITE(LOGNAM,*) "LBITSAFE ", LBITSAFE
WRITE(LOGNAM,*) "LSLOPEMOUTH ", LSLOPEMOUTH

!============================
!*** 2. set model dimention & time

!* defaults (from namelist)
CDIMINFO ="NONE"
DT       = 24*60*60          !! dt = 1day (automatically set by adaptive time step)
IFRQ_INP = 24                !! daily (24h) input

!* change
REWIND(NSETFILE)
READ(NSETFILE,NML=NDIMTIME)

DTIN  = IFRQ_INP*60*60       !! hour -> second

WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "=== NAMELIST, NCONF ==="
WRITE(LOGNAM,*) "CDIMINFO  ", TRIM(CDIMINFO)
WRITE(LOGNAM,*) "DT        ", DT
WRITE(LOGNAM,*) "DTIN      ", DTIN
WRITE(LOGNAM,*) "IFRQ_INP  ", IFRQ_INP

!==========
!* default (from diminfo)
NX    = 1440                 !! 15 minute resolution
NY    = 720
NLFP  = 10                   !! 10 floodplain layer
NXIN  = 360                  !! 1 degree input
NYIN  = 180
INPN  = 1                    !! maximum number of input grids corresponding to one CaMa-Flood grid
WEST  = -180.D0              !! west, east, north, south edges of the domain
EAST  =  180.D0
NORTH =  90.D0
SOUTH = -90.D0

!* value from CDIMINFO
IF( CDIMINFO/="NONE" )THEN
  WRITE(LOGNAM,*) "CMF::CONFIG_NMLIST: read DIMINFO ", TRIM(CDIMINFO)

  TMPNAM=INQUIRE_FID()
  OPEN(TMPNAM,FILE=CDIMINFO,FORM='FORMATTED')
  READ(TMPNAM,*) NX
  READ(TMPNAM,*) NY
  READ(TMPNAM,*) NLFP
  READ(TMPNAM,*) NXIN
  READ(TMPNAM,*) NYIN
  READ(TMPNAM,*) INPN
  READ(TMPNAM,*) 
  IF( LGRIDMAP )THEN
    READ(TMPNAM,*) WEST
    READ(TMPNAM,*) EAST
    READ(TMPNAM,*) NORTH
    READ(TMPNAM,*) SOUTH
  ENDIF
  CLOSE(TMPNAM)
ENDIF

!* check
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "=== DIMINFO ==="
WRITE(LOGNAM,*) "NX,NY,NLFP     ", NX,  NY,  NLFP
WRITE(LOGNAM,*) "NXIN,NYIN,INPN ", NXIN,NYIN,INPN
IF( LGRIDMAP ) THEN
  WRITE(LOGNAM,*) "WEST,EAST,NORTH,SOUTH ", WEST,EAST,NORTH,SOUTH
ENDIF

!============================
!*** 3. set PARAM: parameters
! * defaults
PMANRIV=0.03D0                              !! manning coefficient river
PMANFLD=0.10D0                              !! manning coefficient floodplain
PGRV   =9.8D0                               !! gravity accerelation
PDSTMTH=10000.D0                            !! downstream distance at river mouth [m]
PCADP  =0.7D0                                !! CFL coefficient
PMINSLP=1.D-5                               !! minimum slope (kinematic wave)

IMIS=-9999_JPIM
RMIS=1.E20_JPRM
DMIS=1.E20_JPRB

CSUFBIN='.bin'
CSUFVEC='.vec'
CSUFPTH='.pth'
CSUFCDF='.nc'

! * change
REWIND(NSETFILE)
READ(NSETFILE,NML=NPARAM)

WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "=== NAMELIST, NPARAM ==="
WRITE(LOGNAM,*) "PMANRIV  ", PMANRIV
WRITE(LOGNAM,*) "PMANRIV  ", PMANFLD
WRITE(LOGNAM,*) "PGRV     ", PGRV
WRITE(LOGNAM,*) "PDSTMTH  ", PDSTMTH
WRITE(LOGNAM,*) "PCADP    ", PCADP
WRITE(LOGNAM,*) "PMINSLP  ", PMINSLP
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "IMIS     ", IMIS
WRITE(LOGNAM,*) "RMIS     ", RMIS
WRITE(LOGNAM,*) "DMIS     ", DMIS
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "CSUFBIN  ", TRIM(CSUFBIN)
WRITE(LOGNAM,*) "CSUFVEC  ", TRIM(CSUFVEC)
WRITE(LOGNAM,*) "CSUFPTH  ", TRIM(CSUFPTH)
WRITE(LOGNAM,*) "CSUFCDF  ", TRIM(CSUFCDF)

!===============================
!*** CLOSE FILE 
CLOSE(NSETFILE)

WRITE(LOGNAM,*) "CMF::CONFIG_NMLIST: end "

WRITE(LOGNAM,*) "--------------------!"
WRITE(LOGNAM,*) ""

END SUBROUTINE CMF_CONFIG_NMLIST
!####################################################################




!####################################################################
SUBROUTINE CMF_CONFIG_CHECK
USE YOS_CMF_INPUT,      ONLY: LADPSTP,  LFPLAIN,  LKINE,    LPTHOUT,     &
                            & LROSPLIT, LGDWDLY,  LSEALEV, &
                            & LWEVAP, LWEVAPFIX,LWEXTRACTRIV
USE YOS_CMF_INPUT,      ONLY: DT, DTIN, DTSL
IMPLICIT NONE
!================================================

WRITE(LOGNAM,*) "CMF::CONFIG_CHECK: check setting conflicts"

!*** 1. check for time step
IF ( DT<60 .or. MOD( INT(DT),60 )/=0 ) THEN
  WRITE(LOGNAM,*) "DT= ", DT
  WRITE(LOGNAM,*) "DT should be multiple of 60. CaMa-Flood controls time by MINUTE"
  WRITE(LOGNAM,*) "stop"
  STOP 9
ENDIF

IF ( MOD( INT(DTIN), INT(DT) )/=0 ) THEN
  WRITE(LOGNAM,*) "DTIN, DT= ", DTIN, DT
  WRITE(LOGNAM,*) "DTIN should be multiple of DT"
  WRITE(LOGNAM,*) "stop"
  STOP 9
ENDIF

IF ( LSEALEV .and. MOD( INT(DTSL), INT(DT) )/=0 ) THEN
  WRITE(LOGNAM,*) "DTSL, DT= ", DTIN, DT
  WRITE(LOGNAM,*) "DTSL should be multiple of DT"
  WRITE(LOGNAM,*) "stop"
  STOP 9
ENDIF

!*** 2. check for physics options

IF ( .not.LFPLAIN .AND. .not.LKINE ) THEN
  WRITE(LOGNAM,*) "LFPLAIN=.false. & LKINE=.false."
!   WRITE(LOGNAM,*) "NO FLOODPLAIN OPTION only available with kinematic wave (LKINE=.true.)"
!   WRITE(LOGNAM,*) "stop"
!   STOP 9
  ! E. Dutra Allow this option for testing 
  WRITE(LOGNAM,*) "NO FLOODPLAIN OPTION should only be active with kinematic wave (LKINE=.true.)"
ENDIF

IF ( LKINE .AND. LADPSTP ) THEN
  WRITE(LOGNAM,*) "LKINE=.true. & LADPSTP=.true."
  WRITE(LOGNAM,*) "adaptive time step only available with local inertial equation (LKINE=.false.)"
  WRITE(LOGNAM,*) "STOP"
  STOP 9
ENDIF

IF ( LKINE .AND. LPTHOUT ) THEN
  WRITE(LOGNAM,*) "LKINE=.true. & LPATHOUT=.true."
  WRITE(LOGNAM,*) "bifurcation channel flow only available with local inertial equation (LKINE=.false.)"
  WRITE(LOGNAM,*) "STOP"
  STOP 9
ENDIF

IF ( LGDWDLY .AND. .NOT. LROSPLIT ) THEN
  WRITE(LOGNAM,*) "LGDWDLY=true and LROSPLIT=false"
  WRITE(LOGNAM,*) "Ground water reservoir can only be active when runoff splitting is on"
ENDIF

IF ( LWEVAPFIX .AND. .NOT. LWEVAP ) THEN
  WRITE(LOGNAM,*) "LWEVAPFIX=true and LWEVAP=false"
  WRITE(LOGNAM,*) "LWEVAPFIX can only be active if LWEVAP is active"
ENDIF 
IF ( LWEXTRACTRIV .AND. .NOT. LWEVAP ) THEN
  WRITE(LOGNAM,*) "LWEXTRACTRIV=true and LWEVAP=false"
  WRITE(LOGNAM,*) "LWEXTRACTRIV can only be active if LWEVAP is active"
ENDIF 

WRITE(LOGNAM,*) "CMF::CONFIG_CHECK: end"

END SUBROUTINE CMF_CONFIG_CHECK
!####################################################################

END MODULE CMF_CTRL_NMLIST_MOD
