MODULE CMF_DRV_CONTROL_MOD
!==========================================================
!* PURPOSE: Initialize/Finalize CaMa-Flood Model 
!
!* CONTAINS:
! -- CMF_DRV_INPUT    : Set namelist & logfile
! -- CMF_DRV_INIT     : Initialize        CaMa-Flood
! -- CMF_DRV_END      : Finalize          CaMa-Flood
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
!** shared variables in module
USE PARKIND1,                ONLY: JPIM, JPRB, JPRM
USE YOS_CMF_INPUT,           ONLY: LOGNAM
IMPLICIT NONE
!** local variables
SAVE
REAL(KIND=JPRB)                 :: ZTT0, ZTT1, ZTT2   ! Time elapsed related 
!==========================================================
CONTAINS
!####################################################################
! -- CMF_DRV_INPUT    : Set namelist & logfile
! -- CMF_DRV_INIT     : Initialize        CaMa-Flood
! -- CMF_DRV_END      : Finalize          CaMa-Flood
!
!####################################################################
SUBROUTINE CMF_DRV_INPUT
! Read setting from namelist ("input_flood.nam" as default)
! -- Called from CMF_DRV_INIT
USE YOS_CMF_INPUT,           ONLY: LLOGOUT, LOGNAM, CLOGOUT, CSETFILE, LSEALEV, LOUTPUT
USE CMF_CTRL_NMLIST_MOD,     ONLY: CMF_CONFIG_NMLIST, CMF_CONFIG_CHECK
USE CMF_CTRL_TIME_MOD,       ONLY: CMF_TIME_NMLIST
USE CMF_CTRL_FORCING_MOD,    ONLY: CMF_FORCING_NMLIST
USE CMF_CTRL_BOUNDARY_MOD,   ONLY: CMF_BOUNDARY_NMLIST
USE CMF_CTRL_RESTART_MOD,    ONLY: CMF_RESTART_NMLIST
USE CMF_CTRL_OUTPUT_MOD,     ONLY: CMF_OUTPUT_NMLIST
USE CMF_CTRL_MAPS_MOD,       ONLY: CMF_MAPS_NMLIST
USE CMF_UTILS_MOD,           ONLY: INQUIRE_FID
IMPLICIT NONE
!================================================

!*** 0a. Set log file & namelist
! Preset in YOS_INPUT:  LLOGOUT=.TRUE.   CLOGOUT='./log_CaMa.txt'
! It can be modified in MAIN program before DRV_INPUT
IF( LLOGOUT )THEN
  LOGNAM=INQUIRE_FID()
  OPEN(LOGNAM,FILE=CLOGOUT,FORM='FORMATTED')  
ELSE
  LOGNAM=6  !! use standard output
  CLOGOUT="NONE"
ENDIF

WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!******************************"
WRITE(LOGNAM,*) "CMF::DRV_INPUT: log file:            ", TRIM(CLOGOUT), LOGNAM

!*** 0b. Input namelist filename
! Preset in YOS_INPUT:  CSETFILE="input_cmf.nam"
! It can be modified in MAIN program before DRV_INPUT
WRITE(LOGNAM,*) "CMF::DRV_INPUT: input namelist:      ", TRIM(CSETFILE)

!*** 1. CaMa-Flood configulation namelist
CALL CMF_CONFIG_NMLIST

CALL CMF_TIME_NMLIST

CALL CMF_MAPS_NMLIST

!*** 2. read namelist for each module
CALL CMF_FORCING_NMLIST

IF( LSEALEV )THEN
  CALL CMF_BOUNDARY_NMLIST
ENDIF

CALL CMF_RESTART_NMLIST

IF( LOUTPUT )THEN
  CALL CMF_OUTPUT_NMLIST
ENDIF

WRITE(LOGNAM,*) "CMF::DRV_INPUT: end reading namelist"

!*** 3. check configulation conflicts
CALL CMF_CONFIG_CHECK

WRITE(LOGNAM,*) "CMF::DRV_INPUT: finished"
WRITE(LOGNAM,*) "******************************!"
WRITE(LOGNAM,*) ""

END SUBROUTINE CMF_DRV_INPUT
!####################################################################





!####################################################################
SUBROUTINE CMF_DRV_INIT
! Initialize CaMa-Flood
! -- Called from CMF_DRV_INIT
USE YOS_CMF_INPUT,           ONLY: LRESTART, LSTOONLY, LOUTPUT, LSEALEV
USE CMF_CTRL_TIME_MOD,       ONLY: CMF_TIME_INIT
USE CMF_CTRL_MAPS_MOD,       ONLY: CMF_RIVMAP_INIT,  CMF_TOPO_INIT
USE CMF_CTRL_VARS_MOD,       ONLY: CMF_PROG_INIT,    CMF_DIAG_INIT
USE CMF_CTRL_FORCING_MOD,    ONLY: CMF_FORCING_INIT
USE CMF_CTRL_BOUNDARY_MOD,   ONLY: CMF_BOUNDARY_INIT
USE CMF_CTRL_OUTPUT_MOD,     ONLY: CMF_OUTPUT_INIT
USE CMF_CTRL_RESTART_MOD,    ONLY: CMF_RESTART_INIT
USE CMF_CALC_FLDSTG_MOD,     ONLY: CMF_CALC_FLDSTG
USE CMF_OPT_OUTFLW_MOD,      ONLY: CMF_CALC_OUTPRE
USE CMF_UTILS_MOD,           ONLY: INQUIRE_FID
!$ USE OMP_LIB    
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!******************************!"
WRITE(LOGNAM,*) "CMF::DRV_INIT: initialization start"

!*** 0b. get start time
CALL CPU_TIME(ZTT0)
!$ ZTT0=OMP_GET_WTIME()

!================================================
WRITE(LOGNAM,*) "CMF::DRV_INIT: (1) Set Time"

!*** 1a. Set time related 
CALL CMF_TIME_INIT

!================================================
WRITE(LOGNAM,*) "CMF::DRV_INIT: (2) Set River Map & Topography"

!*** 2a. Read input river map 
CALL CMF_RIVMAP_INIT

!*** 2b. Set topography 
CALL CMF_TOPO_INIT

!================================================
WRITE(LOGNAM,*) "CMF::DRV_INIT: (3) Set output & forcing modules"

!*** 3a. Create Output files 
IF( LOUTPUT )THEN
  CALL CMF_OUTPUT_INIT
ENDIF

!*** 3b. Initialize forcing data
CALL CMF_FORCING_INIT

!*** 3b. Initialize dynamic sea level boundary data
IF( LSEALEV )THEN
  CALL CMF_BOUNDARY_INIT
ENDIF
!================================================
WRITE(LOGNAM,*) "CMF::DRV_INIT: (4) Allocate prog&diag vars & initialize"

!*** 4a. Set initial prognostic variables 
CALL CMF_PROG_INIT

!*** 4b. Initialize (allocate) diagnostic arrays
CALL CMF_DIAG_INIT

IF( LRESTART )THEN
  CALL CMF_RESTART_INIT
ENDIF
!================================================
WRITE(LOGNAM,*) "CMF::DRV_INIT: (5) set flood stage at initial condition"

!*** 5a. Set flood stage
CALL CMF_CALC_FLDSTG

!*** 5b. reconstruct previous t-step flow 
IF( LRESTART .AND. LSTOONLY )THEN
  CALL CMF_CALC_OUTPRE
ENDIF
!================================================

!*** get initialization end time time
CALL CPU_TIME(ZTT1)
!$ ZTT1=OMP_GET_WTIME()

WRITE(LOGNAM,*) "CMF::DRV_INIT: initialization finished in:",ZTT1-ZTT0,' Seconds'
WRITE(LOGNAM,*) "CMF::DRV_INIT: end"

END SUBROUTINE CMF_DRV_INIT
!####################################################################




!####################################################################
SUBROUTINE CMF_DRV_END
! Finalize CaMa-Flood
USE YOS_CMF_INPUT,           ONLY: LOUTPUT, LSEALEV
USE CMF_CTRL_OUTPUT_MOD,     ONLY: CMF_OUTPUT_END
USE CMF_CTRL_FORCING_MOD,    ONLY: CMF_FORCING_END
USE CMF_CTRL_BOUNDARY_MOD,   ONLY: CMF_BOUNDARY_END
!$ USE OMP_LIB    
IMPLICIT NONE 
!==========================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!******************************!"
WRITE(LOGNAM,*) "CMF::DRV_END: finalize forcing & output modules"
CALL CMF_FORCING_END
IF( LOUTPUT )THEN
  CALL CMF_OUTPUT_END
ENDIF
IF( LSEALEV ) THEN
  CALL CMF_BOUNDARY_END
ENDIF

!*** get simulation end time
CALL CPU_TIME(ZTT2)
!$ ZTT2=OMP_GET_WTIME()
WRITE(LOGNAM,*) "CMF::DRV_END: simulation finished in:",ZTT2-ZTT0,' Seconds'

WRITE(LOGNAM,*) "CMF::DRV_END: close logfile"
WRITE(LOGNAM,*) "CMF::===== CALCULATION END ====="
CLOSE(LOGNAM)

END SUBROUTINE CMF_DRV_END
!####################################################################

END MODULE CMF_DRV_CONTROL_MOD
