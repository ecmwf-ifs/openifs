MODULE CMF_CALC_DIAG_MOD
!==========================================================
!* PURPOSE: Manage prognostic/diagnostic variables in CaMa-Flood
!
!* CONTAINS:
! -- CMF_PROG_INIT      : Initialize Prognostic variables (include restart data handling)
! -- CMF_DIAG_INIT      : Initialize Diagnostic variables
! -- CMF_DIAG_AVERAGE   : Calculate time-average of Diagnostic Variables
! -- CMF_DIAG_RESET     : Reset Diagnostic Variables (Average & Maximum )
!!
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
! -- CMF_DIAG_AVE_MAX   : Add / Max of diagnostic variables at time step
! -- CMF_DIAG_AVERAGE   : Calculate time-average of Diagnostic Variables
! -- CMF_DIAG_RESET     : Reset Diagnostic Variables (Average & Maximum )
!
!####################################################################
SUBROUTINE CMF_DIAG_AVEMAX
USE YOS_CMF_INPUT,      ONLY: DT, LPTHOUT
USE YOS_CMF_MAP,        ONLY: NSEQALL,      NPTHOUT
USE YOS_CMF_PROG,       ONLY: D2RIVOUT,     D2FLDOUT,     D1PTHFLW,     D2GDWRTN,     D2RUNOFF,&
                            & D2ROFSUB
USE YOS_CMF_DIAG,       ONLY: D2OUTFLW,     D2RIVVEL,     D2PTHOUT,     D2PTHINF,     D2RIVDPH,&
                            & D2STORGE,     D2WEVAPEX, &
                            & D2RIVOUT_AVG, D2FLDOUT_AVG, D1PTHFLW_AVG, D2GDWRTN_AVG, D2RUNOFF_AVG,&
                            & D2ROFSUB_AVG, D2WEVAPEX_AVG, &
                            & D2OUTFLW_AVG, D2RIVVEL_AVG, D2PTHOUT_AVG, &
                            & NADD, &
                            & D2OUTFLW_MAX, D2RIVDPH_MAX, D2STORGE_MAX
IMPLICIT NONE
INTEGER(KIND=JPIM)            ::  ISEQ, IPTH
!$ SAVE                           ISEQ, IPTH
!====================
NADD=NADD+DT
!$OMP PARALLEL DO
DO ISEQ=1, NSEQALL
  D2RIVOUT_AVG(ISEQ,1)=D2RIVOUT_AVG(ISEQ,1)+D2RIVOUT(ISEQ,1)*DT
  D2FLDOUT_AVG(ISEQ,1)=D2FLDOUT_AVG(ISEQ,1)+D2FLDOUT(ISEQ,1)*DT
  D2RIVVEL_AVG(ISEQ,1)=D2RIVVEL_AVG(ISEQ,1)+D2RIVVEL(ISEQ,1)*DT
  D2OUTFLW_AVG(ISEQ,1)=D2OUTFLW_AVG(ISEQ,1)+D2OUTFLW(ISEQ,1)*DT

  D2PTHOUT_AVG(ISEQ,1)=D2PTHOUT_AVG(ISEQ,1)+D2PTHOUT(ISEQ,1)*DT-D2PTHINF(ISEQ,1)*DT

  D2GDWRTN_AVG(ISEQ,1)=D2GDWRTN_AVG(ISEQ,1)+D2GDWRTN(ISEQ,1)*DT
  D2RUNOFF_AVG(ISEQ,1)=D2RUNOFF_AVG(ISEQ,1)+D2RUNOFF(ISEQ,1)*DT
  D2ROFSUB_AVG(ISEQ,1)=D2ROFSUB_AVG(ISEQ,1)+D2ROFSUB(ISEQ,1)*DT
  D2WEVAPEX_AVG(ISEQ,1)= D2WEVAPEX_AVG(ISEQ,1) +D2WEVAPEX(ISEQ,1)*DT

  D2OUTFLW_MAX(ISEQ,1)=max( D2OUTFLW_MAX(ISEQ,1), abs(D2OUTFLW(ISEQ,1)) )
  D2RIVDPH_MAX(ISEQ,1)=max( D2RIVDPH_MAX(ISEQ,1),     D2RIVDPH(ISEQ,1)  )
  D2STORGE_MAX(ISEQ,1)=max( D2STORGE_MAX(ISEQ,1),     D2STORGE(ISEQ,1)  )
END DO
!$OMP END PARALLEL DO

IF( LPTHOUT )THEN
  !$OMP PARALLEL DO
  DO IPTH=1, NPTHOUT
    D1PTHFLW_AVG(IPTH,:)=D1PTHFLW_AVG(IPTH,:)+D1PTHFLW(IPTH,:)*DT
  END DO
  !$OMP END PARALLEL DO
ENDIF

END SUBROUTINE CMF_DIAG_AVEMAX
!####################################################################





!####################################################################
SUBROUTINE CMF_DIAG_AVERAGE
USE YOS_CMF_TIME,       ONLY: JYYYYMMDD, JHHMM
USE YOS_CMF_DIAG,       ONLY: D2DIAG_AVG, D1PTHFLW_AVG, NADD
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) "CMF::DIAG_AVERAGE: time-average", NADD, JYYYYMMDD, JHHMM
D2DIAG_AVG(:,:,:) = D2DIAG_AVG(:,:,:) /DBLE(NADD)
D1PTHFLW_AVG(:,:) = D1PTHFLW_AVG(:,:) /DBLE(NADD)
END SUBROUTINE CMF_DIAG_AVERAGE
!####################################################################





!####################################################################
SUBROUTINE CMF_DIAG_RESET
USE YOS_CMF_TIME,       ONLY: JYYYYMMDD, JHHMM
USE YOS_CMF_DIAG,       ONLY: D2DIAG_AVG, D1PTHFLW_AVG, D2DIAG_MAX, NADD
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) "CMF::DIAG_AVERAGE: reset", JYYYYMMDD, JHHMM
NADD=0
D2DIAG_AVG(:,:,:) = 0._JPRB
D1PTHFLW_AVG(:,:) = 0._JPRB 
D2DIAG_MAX(:,:,:) = 0._JPRB
END SUBROUTINE CMF_DIAG_RESET
!####################################################################

END MODULE CMF_CALC_DIAG_MOD
