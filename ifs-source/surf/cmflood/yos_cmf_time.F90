MODULE YOS_CMF_TIME
!==========================================================
!* PURPOSE: Shared time-related variables
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
IMPLICIT NONE
!======================================
SAVE
! simulation time step
INTEGER(KIND=JPIM)              :: KSTEP              !! time step since start
INTEGER(KIND=JPIM)              :: NSTEPS             !! total time step (from start to end), given in CMF_TIME_INIT
! elapsed minute from base date (YYYY0,MM0,DD0)
INTEGER(KIND=JPIM)              :: KMIN               !! KMIN at the start of time step
INTEGER(KIND=JPIM)              :: KMINNEXT           !! KMIN at the end   of time step
!
INTEGER(KIND=JPIM)              :: KMINSTART          !! KMIN at the start of simulation
INTEGER(KIND=JPIM)              :: KMINEND            !! KMIN at the end   of simulation
!
INTEGER(KIND=JPIM)              :: KMINSTAIN          !! KMIN at the start of forcing runoff  data (netCDF)
INTEGER(KIND=JPIM)              :: KMINSTASL          !! KMIN at the start of boundary sealev data (netCDF)
! simulation start date:hour (KMINSTART)
INTEGER(KIND=JPIM)              :: ISYYYYMMDD         !! date     at simulation start time
INTEGER(KIND=JPIM)              :: ISHHMM             !! hour+min at simulation start time
INTEGER(KIND=JPIM)              :: ISYYYY
INTEGER(KIND=JPIM)              :: ISMM
INTEGER(KIND=JPIM)              :: ISDD
INTEGER(KIND=JPIM)              :: ISHOUR
INTEGER(KIND=JPIM)              :: ISMIN
! simulation end   date:hour (KMINEND)
INTEGER(KIND=JPIM)              :: IEYYYYMMDD         !! date     of simulation end time
INTEGER(KIND=JPIM)              :: IEHHMM             !! hour+min of simulation end time
INTEGER(KIND=JPIM)              :: IEYYYY
INTEGER(KIND=JPIM)              :: IEMM
INTEGER(KIND=JPIM)              :: IEDD
INTEGER(KIND=JPIM)              :: IEHOUR
INTEGER(KIND=JPIM)              :: IEMIN
!*** date:hour at START of time steop (KMIN)
INTEGER(KIND=JPIM)              :: IYYYYMMDD          !! date     at the start of time-step
INTEGER(KIND=JPIM)              :: IYYYY              !! year     at the start of time-step
INTEGER(KIND=JPIM)              :: IMM                !! month    at the start of time-step
INTEGER(KIND=JPIM)              :: IDD                !! day      at the start of time-step
INTEGER(KIND=JPIM)              :: IHHMM              !! hour+min at the start of time-step
INTEGER(KIND=JPIM)              :: IHOUR              !! hour     at the start of time-step
INTEGER(KIND=JPIM)              :: IMIN               !! min      at the start of time-step
!*** date:hour at END   of time steop (KMINNEXT)
INTEGER(KIND=JPIM)              :: JYYYYMMDD          !! date     at the end   of time-step
INTEGER(KIND=JPIM)              :: JYYYY              !! year     at the end   of time-step
INTEGER(KIND=JPIM)              :: JMM                !! month    at the end   of time-step
INTEGER(KIND=JPIM)              :: JDD                !! day      at the end   of time-step
INTEGER(KIND=JPIM)              :: JHHMM              !! hour+min at the end   of time-step
INTEGER(KIND=JPIM)              :: JHOUR              !! hour     at the end   of time-step
INTEGER(KIND=JPIM)              :: JMIN               !! min      at the end   of time-step

!*** base time to define kmin
INTEGER(KIND=JPIM)              :: YYYY0              !! base year
INTEGER(KIND=JPIM)              :: MM0                !! base month
INTEGER(KIND=JPIM)              :: DD0                !! base day
!==========================================================
END MODULE YOS_CMF_TIME
