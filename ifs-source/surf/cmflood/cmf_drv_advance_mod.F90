MODULE CMF_DRV_ADVANCE_MOD
!==========================================================
!* PURPOSE: Advance CaMa-Flood time integration  
!
!* CONTAINS:
! -- CMF_DRV_ADVANCE : Advance integration for KSPETS (given as argument)
!
!* INTERFACE:
! -- Called from "Main Program" or "Coupler"
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
USE PARKIND1,                ONLY: JPIM, JPRM, JPRB
USE YOS_CMF_INPUT,           ONLY: LOGNAM
IMPLICIT NONE
CONTAINS
!####################################################################
! -- CMF_DRV_ADVANCE : Advance integration for KSPETS
!
!
!####################################################################
SUBROUTINE CMF_DRV_ADVANCE(KSTEPS)
USE YOS_CMF_INPUT,           ONLY: LOUTPUT, LSEALEV, IFRQ_OUT
USE YOS_CMF_TIME,            ONLY: KSTEP, JYYYYMMDD, JHHMM, JHOUR, JMIN

!
USE CMF_CTRL_TIME_MOD,       ONLY: CMF_TIME_NEXT, CMF_TIME_UPDATE
USE CMF_CTRL_PHYSICS_MOD,    ONLY: CMF_ADVANCE_PHYSICS
USE CMF_CTRL_RESTART_MOD,    ONLY: CMF_RESTART_WRITE
USE CMF_CTRL_OUTPUT_MOD,     ONLY: CMF_OUTPUT_WRITE, LOUTCDF
USE CMF_CALC_DIAG_MOD,       ONLY: CMF_DIAG_AVERAGE, CMF_DIAG_RESET
USE CMF_CTRL_BOUNDARY_MOD,   ONLY: CMF_BOUNDARY_UPDATE
!$ USE OMP_LIB
IMPLICIT NONE 
SAVE
! Input argument 
INTEGER(KIND=JPIM)              :: KSTEPS             !! Number of timesteps to advance 
!* Local variables 
INTEGER(KIND=JPIM)              :: ISTEP              !! Time Step
REAL(KIND=JPRB)                 :: ZTT0, ZTT1         !! Time elapsed related 
!$ INTEGER(KIND=JPIM)           :: NTHREADS           !! OpenMP thread number
!==========================================================

!*** get OMP thread number
!$OMP PARALLEL
!$ NTHREADS=OMP_GET_MAX_THREADS()
!$OMP END PARALLEL 

!================================================
!*** START: time step loop
DO ISTEP=1,KSTEPS
  !============================
  !*** 0. get start CPU time
  CALL CPU_TIME(ZTT0)
  !$ ZTT0=OMP_GET_WTIME()
  
  !============================
  !*** 0.1 Write initial conditions in netcdf mode
  IF ( KSTEP == 0 .AND. LOUTPUT .AND. LOUTCDF ) THEN
    CALL CMF_OUTPUT_WRITE
  ENDIF
   
  !============================
  !*** 1. Set next time
  CALL CMF_TIME_NEXT               !! set KMINNEXT, JYYYYMMDD, JHHMM

  !*** (optional)
  IF( LSEALEV )THEN
    CALL CMF_BOUNDARY_UPDATE
  ENDIF

  !============================
  !*** 2. Advance model integration 
  CALL CMF_ADVANCE_PHYSICS

  !============================
  !*** 3. Write output file (when needed)
  IF( LOUTPUT .and. MOD(JHOUR,IFRQ_OUT)==0 .and. JMIN==0 )then
    !*** average variable
    CALL CMF_DIAG_AVERAGE

    !*** write output data
    CALL CMF_OUTPUT_WRITE

    !*** reset variable
    CALL CMF_DIAG_RESET
  ENDIF

  !============================ 
  !*** 4. Write restart file 
  CALL CMF_RESTART_WRITE

  !============================ 
  !*** 5. Update current time      !! Update KMIN, IYYYYMMDD, IHHMM (to KMINNEXT, JYYYYMMDD, JHHMM)
  CALL CMF_TIME_UPDATE

  !============================
  !*** 6. Check CPU time 
  CALL CPU_TIME(ZTT1)
  !$ ZTT1=OMP_GET_WTIME()
  WRITE(LOGNAM,*) "CMF::DRV_ADVANCE END: KSTEP, time (end of Tstep):", KSTEP, JYYYYMMDD, JHHMM
  WRITE(LOGNAM,*) "Elapsed cpu time", ZTT1-ZTT0,"Seconds"
  CALL FLUSH(LOGNAM)

ENDDO
!*** END:time step loop
!================================================

END SUBROUTINE CMF_DRV_ADVANCE
!####################################################################

END MODULE CMF_DRV_ADVANCE_MOD
