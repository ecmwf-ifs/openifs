MODULE CMF_CTRL_VARS_MOD
!==========================================================
!* PURPOSE: Manage prognostic/diagnostic variables in CaMa-Flood
!
!* CONTAINS:
! -- CMF_PROG_INIT      : Initialize Prognostic variables (include restart data handling)
! -- CMF_DIAG_INIT      : Initialize Diagnostic variables
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
! -- CMF_PROG_INIT      : Initialize Prognostic variables (include restart data handling)
! -- CMF_DIAG_INIT      : Initialize Diagnostic variables
!
!####################################################################
SUBROUTINE CMF_PROG_INIT
USE YOS_CMF_MAP,             ONLY: NSEQMAX, NPTHOUT, NPTHLEV
USE YOS_CMF_PROG,            ONLY: ND2PROG,      D2PROG, &
                                 & D2RUNOFF,     D2ROFSUB,     D2GDWSTO,     D2GDWRTN, &
                                 & D2RIVSTO,     D2FLDSTO,     D2RIVOUT,     D2FLDOUT, &
                                 & D2RIVOUT_PRE, D2FLDOUT_PRE, D2RIVDPH_PRE, D2FLDSTO_PRE, &
                                 & D1PTHFLW,     D1PTHFLW_PRE, D2WEVAP
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"

WRITE(LOGNAM,*) "CMF::PROG_INIT: prognostic variable initialization"

!*** 1. ALLOCATE 
ND2PROG=13
ALLOCATE( D2PROG(NSEQMAX,1,ND2PROG)     )
ALLOCATE( D1PTHFLW(NPTHOUT,NPTHLEV)     )
ALLOCATE( D1PTHFLW_PRE(NPTHOUT,NPTHLEV) )

D2RUNOFF     => D2PROG(:,:,1)
D2ROFSUB     => D2PROG(:,:,2)
D2RIVSTO     => D2PROG(:,:,3)
D2FLDSTO     => D2PROG(:,:,4)
D2RIVOUT     => D2PROG(:,:,5)
D2FLDOUT     => D2PROG(:,:,6)
D2RIVOUT_PRE => D2PROG(:,:,7)
D2FLDOUT_PRE => D2PROG(:,:,8)
D2RIVDPH_PRE => D2PROG(:,:,9)
D2FLDSTO_PRE => D2PROG(:,:,10)
D2GDWSTO     => D2PROG(:,:,11)
D2GDWRTN     => D2PROG(:,:,12)
D2WEVAP      => D2PROG(:,:,13)
D2PROG(:,:,:) = 0._JPRB

!============================
!*** 2a. set to zero 
D2RUNOFF(:,:)     = 0._JPRB
D2ROFSUB(:,:)     = 0._JPRB
D2RIVSTO(:,:)     = 0._JPRB
D2FLDSTO(:,:)     = 0._JPRB
D2RIVOUT(:,:)     = 0._JPRB
D2FLDOUT(:,:)     = 0._JPRB
D2RIVOUT_PRE(:,:) = 0._JPRB
D2FLDOUT_PRE(:,:) = 0._JPRB
D2RIVDPH_PRE(:,:) = 0._JPRB
D2FLDSTO_PRE(:,:) = 0._JPRB
D2GDWRTN(:,:)     = 0._JPRB
D2GDWSTO(:,:)     = 0._JPRB
D1PTHFLW(:,:)     = 0._JPRB
D1PTHFLW_PRE(:,:) = 0._JPRB

!***  2b. set initial water surface elevation to sea surface level
WRITE(LOGNAM,*) 'PROG_INIT: fill channels below downstream boundary'
CALL STORAGE_SEA_SURFACE

WRITE(LOGNAM,*) "CMF::PROG_INIT: end"

CONTAINS
!==========================================================
!+ STORAGE_SEA_SURFACE: set initial storage, assuming water surface not lower than downstream sea surface elevation
!+
!+
! ==================================================
SUBROUTINE STORAGE_SEA_SURFACE
! set initial storage, assuming water surface not lower than downstream sea surface elevation
USE YOS_CMF_MAP,  ONLY: NSEQALL,  I1NEXT
USE YOS_CMF_MAP,  ONLY: D2DWNELV, D2RIVELV,D2RIVHGT,D2RIVWTH,D2RIVLEN
IMPLICIT NONE
! local variables
INTEGER(KIND=JPIM) ::  ISEQ, JSEQ, KSEQ
REAL(KIND=JPRB)    ::  DSEAELV, DDPH
!!=================
DO ISEQ=1, NSEQALL
  JSEQ=ISEQ
  DO WHILE( I1NEXT(JSEQ)>0 )
    KSEQ=JSEQ
    JSEQ=I1NEXT(KSEQ)
  END DO

  DSEAELV=D2DWNELV(JSEQ,1) !! downstream boundary elevation

  !! set initial water level to sea level if river bed is lower than sea level
  DDPH=MAX( DSEAELV-D2RIVELV(ISEQ,1),0._JPRB )
  DDPH=MIN( DDPH,D2RIVHGT(ISEQ,1) )
  D2RIVSTO(ISEQ,1)=DDPH*D2RIVLEN(ISEQ,1)*D2RIVWTH(ISEQ,1)
END DO
    
END SUBROUTINE STORAGE_SEA_SURFACE
! ==================================================

END SUBROUTINE CMF_PROG_INIT
!####################################################################






!####################################################################
SUBROUTINE CMF_DIAG_INIT
USE YOS_CMF_MAP,        ONLY: NSEQMAX,NPTHOUT,NPTHLEV
USE YOS_CMF_DIAG,       ONLY: N2DIAG, D2DIAG, &
                            &   D2RIVINF, D2RIVDPH, D2RIVVEL, D2FLDINF, D2FLDDPH, D2FLDFRC, D2FLDARE, &
                            &   D2PTHOUT, D2PTHINF, D2SFCELV, D2OUTFLW, D2STORGE, D2WEVAPEX
USE YOS_CMF_DIAG,       ONLY: N2DIAG_AVG, D2DIAG_AVG, NADD, &
                            &   D2RIVOUT_AVG, D2FLDOUT_AVG, D2OUTFLW_AVG, D2RIVVEL_AVG, D2PTHOUT_AVG, &
                            &   D2GDWRTN_AVG, D2RUNOFF_AVG, D2ROFSUB_AVG, D1PTHFLW_AVG, D2WEVAPEX_AVG 
USE YOS_CMF_DIAG,       ONLY: N2DIAG_MAX, D2DIAG_MAX, &
                            &   D2STORGE_MAX, D2OUTFLW_MAX, D2RIVDPH_MAX
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"

WRITE(LOGNAM,*) "CMF::DIAG_INIT: initialize diagnostic variables"

!*** 1. snapshot 2D diagnostics
N2DIAG=13
ALLOCATE(D2DIAG(NSEQMAX,1,N2DIAG))
D2DIAG(:,:,:) = 0._JPRB
D2RIVINF => D2DIAG(:,:,1)
D2RIVDPH => D2DIAG(:,:,2)
D2RIVVEL => D2DIAG(:,:,3)
D2FLDINF => D2DIAG(:,:,4)
D2FLDDPH => D2DIAG(:,:,5)
D2FLDFRC => D2DIAG(:,:,6)
D2FLDARE => D2DIAG(:,:,7)
D2PTHOUT => D2DIAG(:,:,8)
D2PTHINF => D2DIAG(:,:,9)
D2SFCELV => D2DIAG(:,:,10)
D2OUTFLW => D2DIAG(:,:,11)
D2STORGE => D2DIAG(:,:,12)
D2WEVAPEX => D2DIAG(:,:,13)

!============================
!*** 2a. time-average 2D diagnostics
N2DIAG_AVG=9
ALLOCATE(D2DIAG_AVG(NSEQMAX,1,N2DIAG_AVG))
D2DIAG_AVG(:,:,:) = 0._JPRB 
D2RIVOUT_AVG => D2DIAG_AVG(:,:,1)
D2FLDOUT_AVG => D2DIAG_AVG(:,:,2)
D2OUTFLW_AVG => D2DIAG_AVG(:,:,3)
D2RIVVEL_AVG => D2DIAG_AVG(:,:,4)
D2PTHOUT_AVG => D2DIAG_AVG(:,:,5)
D2GDWRTN_AVG => D2DIAG_AVG(:,:,6)
D2RUNOFF_AVG => D2DIAG_AVG(:,:,7)
D2ROFSUB_AVG => D2DIAG_AVG(:,:,8)
D2WEVAPEX_AVG  => D2DIAG_AVG(:,:,9) 
NADD=0

!*** 2b time-average 1D Diagnostics (bifurcation channel)
ALLOCATE(D1PTHFLW_AVG(NPTHOUT,NPTHLEV))
D1PTHFLW_AVG(:,:) = 0._JPRB 

!============================
!*** 3. Maximum 2D Diagnostics 
N2DIAG_MAX=3
ALLOCATE(D2DIAG_MAX(NSEQMAX,1,N2DIAG_MAX))
D2DIAG_MAX(:,:,:) = 0._JPRB 
D2STORGE_MAX => D2DIAG_MAX(:,:,1)
D2OUTFLW_MAX => D2DIAG_MAX(:,:,2)
D2RIVDPH_MAX => D2DIAG_MAX(:,:,3)

WRITE(LOGNAM,*) "CMF::DIAG_INIT: end"

END SUBROUTINE CMF_DIAG_INIT
!####################################################################

END MODULE CMF_CTRL_VARS_MOD
