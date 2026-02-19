MODULE CMF_CALC_DAMOUT_MOD
!==========================================================
!* PURPOSE: CaMa-Flood reservoir operation scheme (under development)
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
CONTAINS
!####################################################################
! -- CMF_CALC_DAMOUT
! --
! --
!####################################################################
SUBROUTINE CMF_CALC_DAMOUT
USE PARKIND1,           ONLY: JPIM, JPRB
USE YOS_CMF_MAP,        ONLY: I1NEXT,   I2VECTOR
USE YOS_CMF_PROG,       ONLY: D2RIVOUT, D2FLDOUT
USE YOS_CMF_DIAG,       ONLY: D2RIVINF, D2FLDINF
IMPLICIT NONE
!*** Local
INTEGER(KIND=JPIM)         :: IX, IY, ISEQ, JSEQ
REAL(KIND=JPRB)            :: DOUTFLW, DEXCESS     !! Total outflw, Excess flow, 
REAL(KIND=JPRB)            :: DRIVRED, DFLDRED     !! river flow reduction, floodplain flow reduction
!*** DAM Parameter
REAL(KIND=JPRB)            :: DMAXOUT
!================================================
!! CALC_DAMOUT is desinged to replace rivout/fldout by operation rule
!! Below is just an example. 
!================================================

!! for Stung Toreng (E105,94,N13.55)
IX=160   !! please check (ix,iy) of the target point
IY=215
ISEQ=I2VECTOR(IX,IY)  ! convert 2D-Map (ix,iy) to 1D-vector (iseq,1) 
JSEQ=I1NEXT(ISEQ)

!! dam operation (modify rivout & fldout)
DOUTFLW=D2RIVOUT(ISEQ,1)+D2FLDOUT(ISEQ,1)
DMAXOUT=40000.

IF( DOUTFLW > DMAXOUT )THEN
  DEXCESS=DOUTFLW-DMAXOUT
  IF( D2FLDOUT(ISEQ,1)>DEXCESS )THEN
    D2FLDOUT(ISEQ,1)=D2FLDOUT(1,ISEQ)-DEXCESS
    D2FLDINF(JSEQ,1)=D2FLDINF(JSEQ,1)-DEXCESS
  ELSE
    DRIVRED=D2RIVOUT(ISEQ,1)-DMAXOUT
    DFLDRED=D2FLDOUT(ISEQ,1)

    D2RIVOUT(ISEQ,1)=DMAXOUT
    D2FLDOUT(ISEQ,1)=0.0

    D2RIVINF(JSEQ,1)=D2RIVINF(JSEQ,1)-DRIVRED
    D2FLDINF(JSEQ,1)=D2FLDINF(JSEQ,1)-DFLDRED
  ENDIF
ENDIF

END SUBROUTINE CMF_CALC_DAMOUT
!####################################################################

END MODULE CMF_CALC_DAMOUT_MOD
