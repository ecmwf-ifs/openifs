MODULE CMF_CTRL_BOUNDARY_MOD
!==========================================================
!* PURPOSE: Manage CaMa-Flood sea level boundary
!
!* CONTAINS:
! -- CMF_BOUNDARY_NMLIST : Read setting from namelist
! -- CMF_BOUNDARY_INIT   : Initialize boundary data file
! -- CMF_BOUNDARY_UPDATE : Update sea level boundary from file
! -- CMF_BOUNDARY_END    : Finalize   boundary data file
!
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
USE PARKIND1,                ONLY: JPIM, JPRB, JPRM
USE YOS_CMF_INPUT,           ONLY: LOGNAM, IFRQ_SL, DTSL
!============================
IMPLICIT NONE
SAVE
!*** NAMELIST/NBOUND/
! configulation
LOGICAL                         :: LSEALEVCDF          !! true : netCDF sea level boundary
! plain binary data
CHARACTER(LEN=256)              :: CSEALEVDIR          !! Sea level boundary DIRECTORY
CHARACTER(LEN=256)              :: CSEALEVPRE          !! Sea level boundary PREFIX
CHARACTER(LEN=256)              :: CSEALEVSUF          !! Sea level boundary SUFFIX
! netCDF data
CHARACTER(LEN=256)              :: CSEALEVCDF          !! Sea level netCDF file name
CHARACTER(LEN=256)              :: CVNSEALEV           !! Sea Level netCDF variable name
!
INTEGER(KIND=JPIM)              :: SYEARSL             !! Start YEAR  of netCDF sea level
INTEGER(KIND=JPIM)              :: SMONSL              !! Start MONTH of netCDF sea level
INTEGER(KIND=JPIM)              :: SDAYSL              !! Start DAY   of netCDF sea level
INTEGER(KIND=JPIM)              :: SHOURSL             !! Start DAY   of netCDF sea level
! for interporlation (netCDF only)
INTEGER(KIND=JPIM)              :: NSTATIONS           !! Number of ser level station
CHARACTER(LEN=256)              :: CSLMAP              !! Conversion table (Sta -> XY)

NAMELIST/NBOUND/   LSEALEVCDF, CSEALEVDIR, CSEALEVPRE, CSEALEVSUF,&
                     CSEALEVCDF, CVNSEALEV, SYEARSL, SMONSL, SDAYSL, SHOURSL, &
                     NSTATIONS,  CSLMAP
!*** local variables
#ifdef UseCDF_CMF
TYPE TYPESL
CHARACTER(LEN=256)              :: CNAME       !! Netcdf file name
CHARACTER(LEN=256)              :: CVAR        !! Netcdf variable name 
INTEGER(KIND=JPIM)              :: NCID        !! Netcdf file     ID
INTEGER(KIND=JPIM)              :: NVARID      !! Netcdf variable ID
INTEGER(KIND=JPIM)              :: NSTAID      !! Netcdf station  ID
INTEGER(KIND=JPIM)              :: NSTART      !! start date of netcdf (KMIN)
INTEGER(KIND=JPIM)              :: NSTEP       !! steps in netCDF
END TYPE TYPESL
TYPE(TYPESL)                    :: SLCDF  !!  Derived type for Sea Level boundary 

REAL(KIND=JPRM),ALLOCATABLE     :: R1SLIN(:)  ! 1D input boundary condition (m)
INTEGER(KIND=JPIM),ALLOCATABLE  :: I2SLMAP(:,:)
#endif

CONTAINS
!####################################################################
! -- CMF_BOUNDARY_NMLIST : Read setting from namelist
! -- CMF_BOUNDARY_INIT   : Initialize boundary data file
! -- CMF_BOUNDARY_UPDATE : Update sea level boundary from file
! -- CMF_BOUNDARY_END    : Finalize   boundary data file
!####################################################################
SUBROUTINE CMF_BOUNDARY_NMLIST
! reed setting from namelist
! -- Called from CMF_DRV_NMLIST
USE YOS_CMF_INPUT,      ONLY: CSETFILE,NSETFILE,LSEALEV
USE YOS_CMF_TIME,       ONLY: YYYY0
USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"

!*** 1. open namelist
NSETFILE=INQUIRE_FID()
OPEN(NSETFILE,FILE=CSETFILE,STATUS="OLD")
WRITE(LOGNAM,*) "CMF::BOUNDARY_NMLIST: namelist OPEN in unit: ", TRIM(CSETFILE), NSETFILE 

!*** 2. default value
LSEALEVCDF=.FALSE.            !! true for netcdf sea level data

CSEALEVDIR="./sealev/"
CSEALEVPRE="sealev"
CSEALEVSUF=".bin"

CSEALEVCDF="./sealev/"
CVNSEALEV="variable"
CSLMAP="./sealev/"

SYEARSL=0
SMONSL=0 
SDAYSL=0 
SHOURSL=0

IFRQ_SL  = 9999              !! default: dynamic sea level not used

!*** 3. read namelist
REWIND(NSETFILE)
READ(NSETFILE,NML=NBOUND)

IF( LSEALEV )THEN
  WRITE(LOGNAM,*)   "=== NAMELIST, NBOUNDARY ==="
  WRITE(LOGNAM,*)   "LSEALEVCDF: ", LSEALEVCDF
  IF( LSEALEVCDF )THEN
    WRITE(LOGNAM,*) "CSEALEVCDF: ", TRIM(CSEALEVCDF)
    WRITE(LOGNAM,*) "CVNSEALEV:  ", TRIM(CVNSEALEV)
    WRITE(LOGNAM,*) "SYEARSL:    ", SYEARSL
    WRITE(LOGNAM,*) "SMONSL:     ", SMONSL
    WRITE(LOGNAM,*) "SDAYSL:     ", SDAYSL
    WRITE(LOGNAM,*) "SHOURSL:    ", SHOURSL
    WRITE(LOGNAM,*) "NSTATIONS:  ", NSTATIONS
    WRITE(LOGNAM,*) "CSLMAP:     ", TRIM(CSLMAP)
  ELSE
    WRITE(LOGNAM,*) "CSEALEVDIR: ", TRIM(CSEALEVDIR)
    WRITE(LOGNAM,*) "CSEALEVPRE: ", TRIM(CSEALEVPRE)
    WRITE(LOGNAM,*) "CSEALEVSUF: ", TRIM(CSEALEVSUF)
  ENDIF
  WRITE(LOGNAM,*) "IFRQ_SL   ", IFRQ_SL
ENDIF

CLOSE(NSETFILE)

!*** 4. modify base date (shared for KMIN)
IF( LSEALEV .and. LSEALEVCDF )THEN
  YYYY0=MIN(YYYY0,SYEARSL)
  YYYY0=MAX(YYYY0,0)
ENDIF

DTSL  = IFRQ_SL *60          !! min  -> second

WRITE(LOGNAM,*) "CMF::BOUNDARY_NMLIST: end" 

END SUBROUTINE CMF_BOUNDARY_NMLIST
!####################################################################





!####################################################################
SUBROUTINE CMF_BOUNDARY_INIT
IMPLICIT NONE
!####################################################################
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"
WRITE(LOGNAM,*) "CMF::BOUNDARY_INIT: initialize boundary" 


IF( LSEALEVCDF )THEN
#ifdef UseCDF_CMF
  CALL CMF_BOUNDARY_INIT_CDF    !! initialize sea level boundary (netCDF only)
#endif
ENDIF

WRITE(LOGNAM,*) "CMF::BOUNDARY_INIT: end" 



#ifdef UseCDF_CMF
CONTAINS
!==========================================================
!+ CMF_BOUNDARY_INIT_CDF
!==========================================================
SUBROUTINE CMF_BOUNDARY_INIT_CDF
USE YOS_CMF_INPUT,           ONLY: TMPNAM,  DTSL, DT
USE YOS_CMF_TIME,            ONLY: KMINSTASL, KMINSTART, KMINEND
USE CMF_UTILS_MOD,           ONLY: NCERROR, INQUIRE_FID, DATE2MIN
USE NETCDF
IMPLICIT NONE
!* Local Variables 
INTEGER(KIND=JPIM)              :: NTIMEID,NCDFSTP
INTEGER(KIND=JPIM)              :: KMINENDSL
! ===============================================
!*** 1. calculate KMINSTASL (START KMIN for boundary)
KMINSTASL = DATE2MIN(SYEARSL*10000+SMONSL*100+SDAYSL,SHOURSL*100)

!*** 2. Initialize Type for sea level CDF
SLCDF%CNAME=TRIM(CSEALEVCDF)
SLCDF%CVAR=TRIM(CVNSEALEV)
SLCDF%NSTART=KMINSTASL
WRITE(LOGNAM,*) "CMF::BOUNRARY_INIT_CDF:", SLCDF%CNAME, SLCDF%NSTART

!*** Open netCDF sea level File 
CALL NCERROR( NF90_OPEN(SLCDF%CNAME,NF90_NOWRITE,SLCDF%NCID),'OPENING :'//SLCDF%CNAME )
CALL NCERROR( NF90_INQ_VARID(SLCDF%NCID,SLCDF%CVAR,SLCDF%NVARID) )
CALL NCERROR( NF90_INQ_DIMID(SLCDF%NCID,'time',NTIMEID),'GETTING TIME ID Sea Level Boundary')
CALL NCERROR( NF90_INQUIRE_DIMENSION(NCID=SLCDF%NCID,DIMID=NTIMEID,LEN=NCDFSTP),'GETTING TIME LENGTH')

CALL NCERROR( NF90_INQ_DIMID(SLCDF%NCID, 'stations', SLCDF%NSTAID ), 'GETTING STATION ID' ) 
CALL NCERROR( NF90_INQUIRE_DIMENSION(SLCDF%NCID, DIMID=SLCDF%NSTAID, LEN=NSTATIONS ), 'GETTING STATION NUMBER' )

WRITE(LOGNAM,*) "CMF::BOUNDARY_INIT_CDF: CNAME,NCID,VARID", TRIM(SLCDF%CNAME),SLCDF%NCID,SLCDF%NVARID
WRITE(LOGNAM,*) "STATION NUMBER:", NSTATIONS

!*** 4. check sealev forcing time 
IF ( KMINSTART .LT. KMINSTASL ) THEN 
  WRITE(LOGNAM,*) "Run start earlier than boundary data", TRIM(SLCDF%CNAME), KMINSTART, KMINSTASL
  STOP 9
ENDIF

KMINENDSL=KMINSTASL + NCDFSTP*INT(DTSL/60,JPIM)
IF ( KMINEND .GT. KMINENDSL  ) THEN 
  WRITE(LOGNAM,*) "Run end later than sealev data", TRIM(SLCDF%CNAME), KMINEND, KMINENDSL
  STOP 9
ENDIF

!*** 4. conversion table
ALLOCATE( R1SLIN(NSTATIONS)    ) ! 1D input boundary condition (m)
ALLOCATE( I2SLMAP(2,NSTATIONS) ) ! conversion matrix (Station -> XY )

TMPNAM=INQUIRE_FID()
OPEN(TMPNAM,FILE=CSLMAP,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*2*NSTATIONS)
READ(TMPNAM,REC=1) I2SLMAP
CLOSE(TMPNAM)


END SUBROUTINE CMF_BOUNDARY_INIT_CDF
#endif
!==========================================================

END SUBROUTINE CMF_BOUNDARY_INIT
!####################################################################





!####################################################################
SUBROUTINE CMF_BOUNDARY_UPDATE
! read runoff from file
USE YOS_CMF_INPUT,           ONLY: LMEANSL,  LSEALEV,  IFRQ_SL
USE YOS_CMF_TIME,            ONLY: IMIN, IYYYYMMDD, IHHMM
USE YOS_CMF_MAP,             ONLY: D2DWNELV, D2ELEVTN, D2SEALEV, D2MEANSL
IMPLICIT NONE
!* local variable
INTEGER(KIND=JPIM)              :: IUPDATE
!================================================
IUPDATE=0
IF( MOD( INT(IMIN),IFRQ_SL)==0 )THEN
  IUPDATE=1
ENDIF


IF( LSEALEV .and. IUPDATE==1 ) THEN
  WRITE(LOGNAM,*) "CMF::BOUNDARY_UPDATE: update at time: ", IYYYYMMDD, IHHMM


  IF( LSEALEVCDF )THEN
#ifdef UseCDF_CMF
    CALL CMF_BOUNDARY_GET_CDF
#endif
  ELSE
    CALL CMF_BOUNDARY_GET_BIN
  ENDIF
ENDIF

IF( LMEANSL ) THEN
  D2DWNELV(:,:)=D2ELEVTN(:,:) + D2MEANSL(:,:)
ELSE
  D2DWNELV(:,:)=D2ELEVTN(:,:)
ENDIF

IF( LSEALEV ) THEN
  D2DWNELV(:,:)=D2DWNELV(:,:) + D2SEALEV(:,:)
ENDIF

CONTAINS
!==========================================================
!+ CMF_BOUNDARY_GET_BIN
!+ CMF_BOUNDARY_GET_CDF
!==========================================================
SUBROUTINE CMF_BOUNDARY_GET_BIN
USE YOS_CMF_INPUT,           ONLY: TMPNAM, NX, NY
USE YOS_CMF_TIME,            ONLY: IYYYY, IMM, IDD, IHHMM
USE YOS_CMF_MAP,             ONLY: D2SEALEV
USE CMF_UTILS_MOD,           ONLY: MAP2VEC,INQUIRE_FID
IMPLICIT NONE
CHARACTER(LEN=256)              :: CIFNAME             !! INPUT FILE
CHARACTER(LEN=256)              :: CDATE               !!
REAL(KIND=JPRM)                 :: R2TMP(NX,NY)
!====================
!*** 1. set file name
WRITE(CDATE,'(I4.4,I2.2,I2.2,I4.4)') IYYYY,IMM,IDD,IHHMM
CIFNAME = TRIM(CSEALEVDIR)//'/'//TRIM(CSEALEVPRE)//TRIM(CDATE)//TRIM(CSEALEVSUF)
WRITE(LOGNAM,*) "CMF::BOUNDARY_GET_BIN: read sealev:",TRIM(CIFNAME)

!*** open & read sea level
TMPNAM=INQUIRE_FID()
OPEN(TMPNAM,FILE=CIFNAME,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TMP
CLOSE(TMPNAM)
CALL MAP2VEC(R2TMP,D2SEALEV)

END SUBROUTINE CMF_BOUNDARY_GET_BIN
!==========================================================
!+
!+
!+
!==========================================================
#ifdef UseCDF_CMF
SUBROUTINE CMF_BOUNDARY_GET_CDF
USE YOS_CMF_INPUT,           ONLY: DTSL, NX,NY
USE YOS_CMF_TIME,            ONLY: KMIN
USE YOS_CMF_MAP,             ONLY: D2SEALEV
USE CMF_UTILS_MOD,           ONLY: NCERROR, MAP2VEC
USE NETCDF
IMPLICIT NONE
!* Local variables
INTEGER(KIND=JPIM)              :: IRECSL, IX, IY, ISTAT
REAL(KIND=JPRM)                 :: R2TMP(NX,NY)
!===============
!*** 1. calculate irec
IRECSL = ( KMIN-SLCDF%NSTART )*60_JPIM / INT(DTSL,JPIM) + 1     !! (second from netcdf start time) / (input time step)
WRITE(LOGNAM,*) "CMF::BOUNDARY_GET_CDF:", TRIM(SLCDF%CNAME), IRECSL

!*** 2. read sea level
CALL NCERROR( NF90_GET_VAR(SLCDF%NCID,SLCDF%NVARID,R1SLIN,(/1,IRECSL/),(/NSTATIONS,1/)),'READING SEA LEVEL' )

!*** 3. convert 1D station data -> 2D map
R2TMP(:,:)=0.E0
DO ISTAT = 1, NSTATIONS
  IX = I2SLMAP(1,ISTAT)
  IY = I2SLMAP(2,ISTAT)
  IF( (IX .GT. 0) .AND. (IY .GT. 0) ) R2TMP(IX,IY) = R1SLIN(ISTAT)
END DO
CALL MAP2VEC(R2TMP,D2SEALEV)

END SUBROUTINE CMF_BOUNDARY_GET_CDF
#endif
!==========================================================

END SUBROUTINE CMF_BOUNDARY_UPDATE
!####################################################################





!####################################################################
SUBROUTINE CMF_BOUNDARY_END
#ifdef UseCDF_CMF
USE CMF_UTILS_MOD,         ONLY: NCERROR
USE NETCDF
#endif
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) "CMF::BOUNDARY_END: Finalize boundary module"

IF( LSEALEVCDF )THEN
#ifdef UseCDF_CMF
  CALL NCERROR( NF90_CLOSE(SLCDF%NCID))
  WRITE(LOGNAM,*) "Input netcdf sealev closed:",SLCDF%NCID
#endif
ENDIF 

WRITE(LOGNAM,*) "CMF::BOUNDARY_END: end"

END SUBROUTINE CMF_BOUNDARY_END
!####################################################################

END MODULE CMF_CTRL_BOUNDARY_MOD
