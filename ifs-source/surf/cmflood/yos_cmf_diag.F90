MODULE YOS_CMF_DIAG
!==========================================================
!* PURPOSE: Shared diagnostic variables for CaMa-Flood
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
USE PARKIND1, ONLY: JPIM, JPRB, JPRM
IMPLICIT NONE
SAVE
!================================================
!*** Inst. diagnostics 
REAL(KIND=JPRB),ALLOCATABLE,TARGET :: D2DIAG(:,:,:)        !! Storage array for average diagnostics

REAL(KIND=JPRB),POINTER         :: D2RIVINF(:,:)           !! river      inflow   [m3/s] (from upstream)
REAL(KIND=JPRB),POINTER         :: D2RIVDPH(:,:)           !! river      depth    [m]
REAL(KIND=JPRB),POINTER         :: D2RIVVEL(:,:)           !! flow velocity       [m/s]

REAL(KIND=JPRB),POINTER         :: D2FLDINF(:,:)           !! floodplain inflow   [m3/s]
REAL(KIND=JPRB),POINTER         :: D2FLDDPH(:,:)           !! floodplain depth    [m]
REAL(KIND=JPRB),POINTER         :: D2FLDFRC(:,:)           !! flooded    fractipn [m2/m2]
REAL(KIND=JPRB),POINTER         :: D2FLDARE(:,:)           !! flooded    area     [m2]

REAL(KIND=JPRB),POINTER         :: D2PTHOUT(:,:)           !! flood path outflow   [m3/s]
REAL(KIND=JPRB),POINTER         :: D2PTHINF(:,:)           !! flood path inflow   [m3/s]

REAL(KIND=JPRB),POINTER         :: D2SFCELV(:,:)           !! water surface elev  [m]    (elevtn - rivhgt + rivdph)
REAL(KIND=JPRB),POINTER         :: D2OUTFLW(:,:)           !! total outflow       [m3/s] (rivout + fldout)
REAL(KIND=JPRB),POINTER         :: D2STORGE(:,:)           !! total storage       [m3]   (rivsto + fldsto)

REAL(KIND=JPRB),POINTER         :: D2WEVAPEX(:,:)          !! Evaporation water extracted

INTEGER(KIND=JPIM)              :: N2DIAG                  !! number of 2D diagnostics

!================================================
!*** Average diagnostics 
REAL(KIND=JPRB),ALLOCATABLE,TARGET :: D2DIAG_AVG(:,:,:)    !! Storage array for average diagnostics (nseqmax,1,variable)

REAL(KIND=JPRB),POINTER         :: D2RIVOUT_AVG(:,:)       !! average river       discharge
REAL(KIND=JPRB),POINTER         :: D2OUTFLW_AVG(:,:)       !! average total outflow       [m3/s] (rivout + fldout)  !! bugfix v362
REAL(KIND=JPRB),POINTER         :: D2FLDOUT_AVG(:,:)       !! average floodplain  discharge
REAL(KIND=JPRB),POINTER         :: D2RIVVEL_AVG(:,:)       !! average flow velocity
REAL(KIND=JPRB),POINTER         :: D2PTHOUT_AVG(:,:)       !! flood pathway net outflow (2D)

REAL(KIND=JPRB),POINTER         :: D2GDWRTN_AVG(:,:)       !! average ground water return flow
REAL(KIND=JPRB),POINTER         :: D2RUNOFF_AVG(:,:)       !! average input runoff
REAL(KIND=JPRB),POINTER         :: D2ROFSUB_AVG(:,:)       !! average input sub-surface runoff
REAL(KIND=JPRB),POINTER         :: D2WEVAPEX_AVG(:,:)      !! average extracted water evaporation

INTEGER(KIND=JPIM)              :: N2DIAG_AVG              !! Number of 2D diagnostics averages
REAL(KIND=JPRB)                 :: NADD                    !! sum DT to calculate average
!*** Average diagnostics (1D)
REAL(KIND=JPRB),ALLOCATABLE  :: D1PTHFLW_AVG(:,:)          !! bifurcation channel flow (1D, not 2D variable)

!================================================
!*** Daily max diagnostics 
REAL(KIND=JPRB),ALLOCATABLE,TARGET :: D2DIAG_MAX(:,:,:)    !! Storage array for maximum diagnostics (nseqmax,1,variable)

REAL(KIND=JPRB),POINTER         :: D2OUTFLW_MAX(:,:)       !! max total outflow       [m3/s] (rivout + fldout)
REAL(KIND=JPRB),POINTER         :: D2STORGE_MAX(:,:)       !! max total outflow       [m3/s] (rivout + fldout)
REAL(KIND=JPRB),POINTER         :: D2RIVDPH_MAX(:,:)       !! max total outflow       [m3/s] (rivout + fldout)

INTEGER(KIND=JPRB)              :: N2DIAG_MAX              !! Number of 2D diagnostics maximum

!================================================
!*** Global total
REAL(KIND=JPRB)                 :: DGLBRIVSTO              !! global river storage      [m3]
REAL(KIND=JPRB)                 :: DGLBFLDSTO              !! global floodplain storage [m3]
REAL(KIND=JPRB)                 :: DGLBSTOPRE              !! global water storage      [m3] (befre flow calculation)
REAL(KIND=JPRB)                 :: DGLBSTONXT              !! global water storage      [m3] (after flow calculation)
REAL(KIND=JPRB)                 :: DGLBSTONEW              !! global water storage      [m3] (after runoff input)
REAL(KIND=JPRB)                 :: DGLBRIVINF              !! global inflow             [m3] (rivinf + fldinf)
REAL(KIND=JPRB)                 :: DGLBRIVOUT              !! global outflow            [m3] (rivout + fldout)
REAL(KIND=JPRB)                 :: DGLBFLDARE              !! global flooded area       [m2]

END MODULE YOS_CMF_DIAG
