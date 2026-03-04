MODULE YOS_CMF_PROG
!==========================================================
!* PURPOSE: Shared prognostic variables for CaMa-Flood
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  2Aug 2019
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
!*** prognostics / state variables initial conditions
REAL(KIND=JPRB),ALLOCATABLE,TARGET :: D2PROG(:,:,:)      !! Array to store 2D variables 
INTEGER(KIND=JPIM)              :: ND2PROG

!================================================
!*** input runoff (interporlated)
REAL(KIND=JPRB),POINTER         :: D2RUNOFF(:,:)         !! input runoff             [m3/s]
REAL(KIND=JPRB),POINTER         :: D2ROFSUB(:,:)         !! input sub-surface runoff [m3/s]
REAL(KIND=JPRB),POINTER         :: D2WEVAP(:,:)          !! input Evaporation [m3/s]

!================================================
!*** river & floodpain
REAL(KIND=JPRB),POINTER         :: D2RIVSTO(:,:)         !! river      storage [m3]
REAL(KIND=JPRB),POINTER         :: D2FLDSTO(:,:)         !! floodplain storage [m3]
REAL(KIND=JPRB),POINTER         :: D2RIVOUT(:,:)         !! river      outflow [m3/s]
REAL(KIND=JPRB),POINTER         :: D2FLDOUT(:,:)         !! floodplain outflow [m3/s]

!================================================
!*** for implicit schemes of the local inertial equation
REAL(KIND=JPRB),POINTER         :: D2RIVOUT_PRE(:,:)     !! river      outflow [m3/s] (prev t-step)
REAL(KIND=JPRB),POINTER         :: D2RIVDPH_PRE(:,:)     !! river      depth   [m]    (prev t-step)
REAL(KIND=JPRB),POINTER         :: D2FLDOUT_PRE(:,:)     !! floodplain outflow [m3/s] (prev t-step)
REAL(KIND=JPRB),POINTER         :: D2FLDSTO_PRE(:,:)     !! floodplain storage [m3]   (prev t-step)

!================================================
!*** Groundwater Delay
REAL(KIND=JPRB),POINTER         :: D2GDWSTO(:,:)         !! ground water storage  [m3]
REAL(KIND=JPRB),POINTER         :: D2GDWRTN(:,:)         !! Ground water return flow [m3/s]

!================================================
!*** These have a different share, not part of the D2PROG array
REAL(KIND=JPRB),ALLOCATABLE     :: D1PTHFLW(:,:)         !! flood path outflow [m3/s]
REAL(KIND=JPRB),ALLOCATABLE     :: D1PTHFLW_PRE(:,:)     !! flood path outflow [m3/s] (prev t-step)

END MODULE YOS_CMF_PROG
