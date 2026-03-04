PROGRAM MAIN_cmf
!==========================================================
!* PURPOSE: CaMa-Flood default stand-alone driver
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
USE PARKIND1,                ONLY: JPRB, JPRM, JPIM
USE YOS_CMF_INPUT,           ONLY: NXIN, NYIN, DT,DTIN
USE YOS_CMF_TIME,            ONLY: NSTEPS
!
USE CMF_DRV_CONTROL_MOD,     ONLY: CMF_DRV_INPUT,   CMF_DRV_INIT,    CMF_DRV_END
USE CMF_DRV_ADVANCE_MOD,     ONLY: CMF_DRV_ADVANCE
USE CMF_CTRL_FORCING_MOD,    ONLY: CMF_FORCING_GET, CMF_FORCING_PUT
!$ USE OMP_LIB
IMPLICIT NONE
!** local variables
INTEGER(KIND=JPIM)              :: ISTEP              ! total time step
INTEGER(KIND=JPIM)              :: ISTEPADV           ! time step to be advanced within DRV_ADVANCE
REAL(KIND=JPRB),ALLOCATABLE     :: ZBUFF(:,:,:)       ! Buffer to store forcing runoff
!================================================

!*** 1a. Namelist handling
CALL CMF_DRV_INPUT
!*** 1b. INITIALIZATION
CALL CMF_DRV_INIT

!*** 1c. allocate data buffer for input forcing
ALLOCATE(ZBUFF(NXIN,NYIN,3))

!============================
!*** 2. MAIN TEMPORAL LOOP / TIME-STEP (NSTEPS calculated by DRV_INIT)

ISTEPADV=INT(DTIN/DT,JPIM)
DO ISTEP=1,NSTEPS,ISTEPADV

  !*  2a Read forcing from file, This is only relevant in Stand-alone mode 
  CALL CMF_FORCING_GET(ZBUFF(:,:,:))

  !*  2b Interporlate runoff & send to CaMa-Flood 
  CALL CMF_FORCING_PUT(ZBUFF(:,:,:))
 
  !*  2c  Advance CaMa-Flood model for ISTEPADV
  CALL CMF_DRV_ADVANCE(ISTEPADV)

ENDDO
!============================

!*** 3a. finalize CaMa-Flood 
DEALLOCATE(ZBUFF)
CALL CMF_DRV_END


!================================================

END PROGRAM MAIN_cmf
!####################################################################
