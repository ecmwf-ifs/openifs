MODULE CMF_CTRL_PHYSICS_MOD
!==========================================================
!* PURPOSE: call CaMa-Flood physics
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
! -- CMF_ADVANCE_PHYSICS
! --
! --
!####################################################################
SUBROUTINE CMF_ADVANCE_PHYSICS
USE PARKIND1,              ONLY: JPIM,   JPRB,    JPRM
USE YOS_CMF_INPUT,         ONLY: LOGNAM, DT,      LADPSTP
USE YOS_CMF_INPUT,         ONLY: LKINE,  LSLPMIX, LFLDOUT,   LPTHOUT,   LDAMOUT
USE YOS_CMF_PROG,          ONLY: D2FLDOUT
!
USE CMF_CALC_FLDSTG_MOD,   ONLY: CMF_CALC_FLDSTG
USE CMF_CALC_OUTFLW_MOD,   ONLY: CMF_CALC_OUTFLW
USE CMF_CALC_PTHOUT_MOD,   ONLY: CMF_CALC_PTHOUT
USE CMF_CALC_DAMOUT_MOD,   ONLY: CMF_CALC_DAMOUT
USE CMF_CALC_STONXT_MOD,   ONLY: CMF_CALC_STONXT
USE CMF_OPT_OUTFLW_MOD,    ONLY: CMF_CALC_OUTFLW_KINEMIX, CMF_CALC_OUTFLW_KINE
USE CMF_CALC_DIAG_MOD,     ONLY: CMF_DIAG_AVEMAX
IMPLICIT NONE
!! LOCAL
INTEGER(KIND=JPIM)            ::  IT, NT
REAL(KIND=JPRB)               ::  DT_DEF
!================================================
DT_DEF=DT

NT=1
IF( LADPSTP )THEN    ! adoptive time step
  CALL CALC_ADPSTP
ENDIF

!! ==========
DO IT=1, NT
!=== 1. Calculate river discharge 
  IF ( LKINE ) THEN
    CALL CMF_CALC_OUTFLW_KINE       !!  OPTION: kinematic
  ELSEIF( LSLPMIX ) THEN
    CALL CMF_CALC_OUTFLW_KINEMIX    !!  OPTION: mix local-inertial & kinematic based on slope
  ELSE
    CALL CMF_CALC_OUTFLW            !!  Default: Local inertial
  ENDIF

  IF( .not. LFLDOUT ) D2FLDOUT(:,:)=0.D0    !! OPTION: no high-water channel flow

! ---
  IF( LPTHOUT )THEN
    CALL CMF_CALC_PTHOUT            !! bifurcation channel flow
  ENDIF
! ---
  IF ( LDAMOUT ) THEN
    CALL CMF_CALC_DAMOUT            !! reservoir operation
  ENDIF

!=== 2.  Calculate the storage in the next time step in FTCS diff. eq.
  CALL CMF_CALC_STONXT

!=== 3. calculate river and floodplain staging
  CALL CMF_CALC_FLDSTG

!=== 4.  write water balance monitoring to IOFILE
  CALL CALC_WATBAL(IT)

!=== 5. calculate averages, maximum
  CALL CMF_DIAG_AVEMAX

END DO

DT=DT_DEF   !! reset DT

CONTAINS
!==========================================================
!+ CALC_ADPSTP
!+ CALC_WATBAL(IT)
!+
!==========================================================
SUBROUTINE CALC_ADPSTP
USE YOS_CMF_INPUT,      ONLY: PGRV, PDSTMTH, PCADP
USE YOS_CMF_MAP,        ONLY: D2NXTDST
USE YOS_CMF_MAP,        ONLY: NSEQALL,NSEQRIV
USE YOS_CMF_DIAG,       ONLY: D2RIVDPH
IMPLICIT NONE
!$ SAVE
INTEGER(KIND=JPIM)         :: ISEQ
REAL(KIND=JPRB)            :: DT_MIN
REAL(KIND=JPRB)            :: DDPH, DDST
!$OMP THREADPRIVATE(DDPH,DDST)
!================================================

DT_MIN=DT_DEF
!$OMP PARALLEL DO REDUCTION(MIN:DT_MIN)
DO ISEQ=1, NSEQRIV
  DDPH=MAX(D2RIVDPH(ISEQ,1),0.01D0 )
  DDST=D2NXTDST(ISEQ,1)
  DT_MIN=min( DT_MIN, PCADP*DDST * (PGRV*DDPH)**(-0.5) )
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO REDUCTION(MIN:DT_MIN)
DO ISEQ=NSEQRIV+1, NSEQALL
  DDPH=MAX(D2RIVDPH(ISEQ,1),0.01D0 )
  DDST=PDSTMTH
  DT_MIN=min( DT_MIN, PCADP*DDST * (PGRV*DDPH)**(-0.5) )
END DO
!$OMP END PARALLEL DO
NT=INT( DT_DEF * DT_MIN**(-1.) -0.01 )+1
DT=DT_DEF * REAL(NT)**(-1.)

IF( NT>=2 ) WRITE(LOGNAM,'(A15,I4,3F10.2)') "ADPSTP: NT=",NT, DT_DEF, DT_MIN, DT

END SUBROUTINE CALC_ADPSTP
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE CALC_WATBAL(IT)
USE YOS_CMF_TIME,            ONLY: KMIN
USE YOS_CMF_DIAG,            ONLY: DGLBRIVSTO,DGLBFLDSTO,DGLBSTOPRE,DGLBSTONXT,DGLBSTONEW,DGLBRIVINF,DGLBRIVOUT,DGLBFLDARE
USE CMF_UTILS_MOD,           ONLY: MIN2DATE,SPLITDATE,SPLITHOUR
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)   :: IT        !! step in adaptive time loop
!*** LOCAL
REAL(KIND=JPRB)                 :: DMISSING !! water ballance error [kg]
!*** local physics time
INTEGER(KIND=JPIM)              :: PKMIN
INTEGER(KIND=JPIM)              :: PYEAR, PMON, PDAY, PHOUR, PMIN
INTEGER(KIND=JPIM)              :: PYYYYMMDD, PHHMM
!*** PARAMETER
REAL(KIND=JPRB)                 ::  DORD
PARAMETER                          (DORD=1.D-9)
! ================================================
PKMIN=INT ( KMIN + IT*DT/60_JPRB )
CALL MIN2DATE(PKMIN,PYYYYMMDD,PHHMM)
CALL SPLITDATE(PYYYYMMDD,PYEAR,PMON,PDAY)
CALL SPLITHOUR(PHHMM,PHOUR,PMIN)

DMISSING   = DGLBSTOPRE - DGLBSTONXT + DGLBRIVINF - DGLBRIVOUT
WRITE(LOGNAM,'(I4.4,4(A1,I2.2),I6,3F12.3,2x,2F12.3,2x,2F12.3,G12.3,F12.3)') &
  PYEAR, '/', PMON, '/', PDAY, '_', PHOUR, ':', PMIN, IT,   &
  DGLBSTOPRE*DORD, DGLBSTONXT*DORD, DGLBSTONEW*DORD ,DGLBRIVSTO*DORD, DGLBFLDSTO*DORD, &
  DGLBRIVINF*DORD, DGLBRIVOUT*DORD, DMISSING*DORD,   DGLBFLDARE*DORD

END SUBROUTINE CALC_WATBAL
!==========================================================

END SUBROUTINE CMF_ADVANCE_PHYSICS
!###############################################################


END MODULE CMF_CTRL_PHYSICS_MOD
