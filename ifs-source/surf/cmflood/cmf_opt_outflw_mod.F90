MODULE CMF_OPT_OUTFLW_MOD
!==========================================================
!* PURPOSE: CaMa-Flood optional physics for discharge 
!
!* CONTAINS:
! -- CMF_CALC_OUTFLW_KINE:    calculate outflow by kinematic wave
! -- CMF_CALC_OUTFLW_KINEMIX: calculate outflow by kinematic wave + local inercial using slope threshold
! -- CMF_CALC_OUTPRE:         reconstruct previous time-step outflow by diffusion waev, for LSTOONLY restart
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
! CMF_CALC_OUTFLW_KINE:    calculate outflow by kinematic wave
! CMF_CALC_OUTFLW_KINEMIX: calculate outflow by kinematic wave + local inercial using slope threshold
! CMF_CALC_OUTPRE:         reconstruct previous time-step outflow by diffusion waev, for LSTOONLY restart
!
!####################################################################
SUBROUTINE CMF_CALC_OUTFLW_KINE
! Calculate discharge, mix kinematic & local inertial, depending on slope
USE PARKIND1,           ONLY: JPIM, JPRB
USE YOS_CMF_INPUT,      ONLY: DT,       LBITSAFE, PMANFLD,  PMINSLP,&
                             &LSLOPEMOUTH
USE YOS_CMF_MAP,        ONLY: I1NEXT,   NSEQALL,  NSEQRIV
USE YOS_CMF_MAP,        ONLY: D2RIVELV, D2ELEVTN, D2NXTDST, D2RIVWTH
USE YOS_CMF_MAP,        ONLY: D2RIVLEN, D2RIVMAN, D2ELEVSLOPE
USE YOS_CMF_PROG,       ONLY: D2RIVSTO, D2RIVOUT, D2FLDSTO, D2FLDOUT
USE YOS_CMF_DIAG,       ONLY: D2RIVDPH, D2RIVVEL, D2RIVINF, D2FLDDPH, D2FLDINF, D2SFCELV
IMPLICIT NONE
!$ SAVE
      INTEGER(KIND=JPIM)   ::  ISEQ, JSEQ
      REAL(KIND=JPRB)      ::  DSLOPE,   DAREA , DVEL,   DSLOPE_F, DARE_F, DVEL_F
!$OMP THREADPRIVATE     (JSEQ, DSLOPE,   DAREA , DVEL,   DSLOPE_F, DARE_F, DVEL_F)
!================================================

!*** 0. calculate surface water elevation, reset inflow

!$OMP PARALLEL DO
DO ISEQ=1, NSEQALL
  D2SFCELV(ISEQ,1)     = D2RIVELV(ISEQ,1) + D2RIVDPH(ISEQ,1)
  D2RIVINF(ISEQ,1) = 0.D0
  D2FLDINF(ISEQ,1) = 0.D0
END DO
!$OMP END PARALLEL DO

!============================
!*** 1a. discharge for usual river grid

!$OMP PARALLEL DO
DO ISEQ=1, NSEQRIV
  JSEQ   = I1NEXT(ISEQ)
! === river flow
  DSLOPE = (D2ELEVTN(ISEQ,1)-D2ELEVTN(JSEQ,1)) * D2NXTDST(ISEQ,1)**(-1.D0)
  DSLOPE = max(DSLOPE,PMINSLP)
  DVEL   = D2RIVMAN(ISEQ,1)**(-1.D0) * DSLOPE**0.5D0 * D2RIVDPH(ISEQ,1)**(2.D0/3.D0)
  DAREA  = D2RIVWTH(ISEQ,1) * D2RIVDPH(ISEQ,1)

  D2RIVVEL(ISEQ,1) = DVEL
  D2RIVOUT(ISEQ,1) = DAREA * DVEL
  D2RIVOUT(ISEQ,1) = MIN(  D2RIVOUT(ISEQ,1), D2RIVSTO(ISEQ,1)/DT )
!=== floodplain flow
  DSLOPE_F = min( 0.005D0,DSLOPE )    !! set min [instead of using weir equation for efficiency]
  DVEL_F   = PMANFLD**(-1.D0) * DSLOPE_F**0.5D0 * D2FLDDPH(ISEQ,1)**(2.D0/3.D0)
  DARE_F   = D2FLDSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.D0)
  DARE_F   = MAX( DARE_F - D2FLDDPH(ISEQ,1)*D2RIVWTH(ISEQ,1), 0.D0 )   !! remove above river channel area

  D2FLDOUT(ISEQ,1) = DARE_F * DVEL_F
  D2FLDOUT(ISEQ,1) = MIN(  D2FLDOUT(ISEQ,1), D2FLDSTO(ISEQ,1)/DT )
END DO
!$OMP END PARALLEL DO

IF( LBITSAFE )THEN  !! Aboid OMP ATOMIC for bit identical simulation
  DO ISEQ=1, NSEQRIV ! for normal pixels
    JSEQ=I1NEXT(ISEQ)
    D2RIVINF(JSEQ,1) = D2RIVINF(JSEQ,1) + D2RIVOUT(ISEQ,1)             !! total inflow to a grid (from upstream)
    D2FLDINF(JSEQ,1) = D2FLDINF(JSEQ,1) + D2FLDOUT(ISEQ,1)
  END DO
ELSE
!$OMP PARALLEL DO
  DO ISEQ=1, NSEQRIV ! for normal pixels
    JSEQ=I1NEXT(ISEQ)
!$OMP ATOMIC
    D2RIVINF(JSEQ,1) = D2RIVINF(JSEQ,1) + D2RIVOUT(ISEQ,1)             !! total inflow to a grid (from upstream)
!$OMP ATOMIC
    D2FLDINF(JSEQ,1) = D2FLDINF(JSEQ,1) + D2FLDOUT(ISEQ,1)
  END DO
!$OMP END PARALLEL DO
ENDIF

!============================
!*** 1b. discharge for river mouth grids

!$OMP PARALLEL DO
DO ISEQ=NSEQRIV+1, NSEQALL
!=== Kinematic approach, river mouth flow
  IF ( LSLOPEMOUTH ) THEN
    ! prescribed slope 
    DSLOPE = D2ELEVSLOPE(ISEQ,1)
  ELSE
    DSLOPE = PMINSLP
  ENDIF
  DVEL   = D2RIVMAN(ISEQ,1)**(-1.) * DSLOPE**0.5 * D2RIVDPH(ISEQ,1)**(2./3.)
  DAREA  = D2RIVWTH(ISEQ,1) * D2RIVDPH(ISEQ,1)

  D2RIVVEL(ISEQ,1) = DVEL
  D2RIVOUT(ISEQ,1) = DAREA * DVEL
  D2RIVOUT(ISEQ,1) = MIN(  D2RIVOUT(ISEQ,1), D2RIVSTO(ISEQ,1)/DT )

!=== kinematic, floodplain mouth flow
  DSLOPE_F = min( 0.005D0,DSLOPE )    !! set min [instead of using weir equation for efficiency]
  DVEL_F  = PMANFLD**(-1.D0) * DSLOPE_F**0.5D0 * D2FLDDPH(ISEQ,1)**(2.D0/3.D0)
  DARE_F   = D2FLDSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.D0)
  DARE_F   = MAX( DARE_F - D2FLDDPH(ISEQ,1)*D2RIVWTH(ISEQ,1), 0.D0 )   !! remove above river channel area

  D2FLDOUT(ISEQ,1) = DARE_F * DVEL_F
  D2FLDOUT(ISEQ,1) = MIN(  D2FLDOUT(ISEQ,1), D2FLDSTO(ISEQ,1)/DT )
END DO
!$OMP END PARALLEL DO

END SUBROUTINE CMF_CALC_OUTFLW_KINE
!####################################################################





!####################################################################
SUBROUTINE CMF_CALC_OUTFLW_KINEMIX
! Calculate discharge, mix kinematic & local inertial, depending on slope
USE PARKIND1,           ONLY: JPIM, JPRB
USE YOS_CMF_INPUT,      ONLY: DT,       LBITSAFE
USE YOS_CMF_INPUT,      ONLY: PDSTMTH,  PMANFLD,  PGRV , PMINSLP
USE YOS_CMF_MAP,        ONLY: I1NEXT,   NSEQALL,  NSEQRIV,  NSEQMAX
USE YOS_CMF_MAP,        ONLY: D2RIVELV, D2ELEVTN, D2NXTDST, D2RIVWTH, D2RIVHGT
USE YOS_CMF_MAP,        ONLY: D2RIVLEN, D2RIVMAN, I2MASK, D2DWNELV
USE YOS_CMF_PROG,       ONLY: D2RIVSTO, D2RIVOUT, D2FLDSTO, D2FLDOUT
USE YOS_CMF_PROG,       ONLY: D2RIVOUT_PRE, D2RIVDPH_PRE, D2FLDOUT_PRE, D2FLDSTO_PRE
USE YOS_CMF_DIAG,       ONLY: D2RIVDPH, D2RIVVEL, D2RIVINF, D2FLDDPH, D2FLDINF, D2SFCELV
IMPLICIT NONE
!*** Local
REAL(KIND=JPRB)            :: D2SFCELV_PRE(NSEQMAX,1)                  !! water surface elevation (t-1) [m]
REAL(KIND=JPRB)            :: D2FLDDPH_PRE(NSEQMAX,1)                  !! floodplain depth (t-1)        [m]
REAL(KIND=JPRB)            :: D2STOOUT(NSEQMAX,1)                      !! total outflow from a grid     [m3]
REAL(KIND=JPRB)            :: D2RATE(NSEQMAX,1)                        !! outflow correction

INTEGER(KIND=JPIM),SAVE    :: ISEQ, JSEQ
REAL(KIND=JPRB),SAVE       :: DSLOPE,   DOUT_PRE,   DFLW,   DFLW_PRE,   DFLW_IMP,   DAREA , DVEL
REAL(KIND=JPRB),SAVE       :: DSLOPE_F, DOUT_PRE_F, DFLW_F, DFLW_PRE_F, DFLW_IMP_F, DARE_F, DVEL_F, DARE_PRE_F, DARE_IMP_F
REAL(KIND=JPRB),SAVE       :: OUT_R1, OUT_R2, OUT_F1, OUT_F2, DIUP, DIDW, DSFCMAX, DSFCMAX_PRE
!$OMP THREADPRIVATE    (JSEQ, DSLOPE,   DOUT_PRE,   DFLW,   DFLW_PRE,   DFLW_IMP,   DAREA,  DVEL)
!$OMP THREADPRIVATE    (      DSLOPE_F, DOUT_PRE_F, DFLW_F, DFLW_PRE_F, DFLW_IMP_F, DARE_F, DVEL_F, DARE_PRE_F, DARE_IMP_F)
!$OMP THREADPRIVATE    (      OUT_R1, OUT_R2, OUT_F1, OUT_F2, DIUP, DIDW, DSFCMAX, DSFCMAX_PRE)
!================================================

!$OMP PARALLEL DO
DO ISEQ=1, NSEQALL
  D2SFCELV(ISEQ,1)     = D2RIVELV(ISEQ,1) + D2RIVDPH(ISEQ,1)
  D2SFCELV_PRE(ISEQ,1) = D2RIVELV(ISEQ,1) + D2RIVDPH_PRE(ISEQ,1)
  D2FLDDPH_PRE(ISEQ,1) = MAX( D2RIVDPH_PRE(ISEQ,1)-D2RIVHGT(ISEQ,1), 0.D0 )

  D2RIVINF(ISEQ,1) = 0.D0
  D2FLDINF(ISEQ,1) = 0.D0
  D2STOOUT(ISEQ,1) = 0.D0
  D2RATE(ISEQ,1) = 1._JPRB
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
DO ISEQ=1, NSEQRIV                                                    !! for normal cells
  JSEQ=I1NEXT(ISEQ) ! next cell's pixel
  
  IF (I2MASK(ISEQ,1) == 0 ) THEN

    DSFCMAX    =MAX( D2SFCELV(ISEQ,1),    D2SFCELV(JSEQ,1) )
    DSFCMAX_PRE=MAX( D2SFCELV_PRE(ISEQ,1),D2SFCELV_PRE(JSEQ,1) )
    DSLOPE = ( D2SFCELV(ISEQ,1)-D2SFCELV(JSEQ,1) ) * D2NXTDST(ISEQ,1)**(-1.D0)
    DSLOPE_F = MAX( -0.005D0, min( 0.005D0,DSLOPE ))    !! set max&min [instead of using weir equation for efficiency]

  !=== River Flow ===
    DFLW   = DSFCMAX - D2RIVELV(ISEQ,1)        !!  flow cross-section depth
    DAREA  = D2RIVWTH(ISEQ,1) * DFLW                                            !!  flow cross-section area

    DFLW_PRE=DSFCMAX_PRE - D2RIVELV(ISEQ,1)
    DFLW_IMP=MAX( (DFLW*DFLW_PRE)**0.5D0 ,1.D-6 )                                            !! semi implicit flow depth

    IF( DFLW_IMP>1.D-5 .and. DAREA>1.D-5 )THEN 
      DOUT_PRE= D2RIVOUT_PRE(ISEQ,1) * D2RIVWTH(ISEQ,1)**(-1.D0)                         !! outflow (t-1) [m2/s] (unit width)
      D2RIVOUT(ISEQ,1) = D2RIVWTH(ISEQ,1) * ( DOUT_PRE + PGRV*DT*DFLW_IMP*DSLOPE ) &
                               * ( 1.D0 + PGRV*DT*D2RIVMAN(ISEQ,1)**2.D0*abs(DOUT_PRE)*DFLW_IMP**(-7.D0/3.D0) )**(-1.D0)
      D2RIVVEL(ISEQ,1) = D2RIVOUT(ISEQ,1) * DAREA**(-1.D0)
    ELSE
      D2RIVOUT(ISEQ,1) = 0.D0
      D2RIVVEL(ISEQ,1) = 0.D0
    ENDIF

  !=== Floodplain Flow ===
    DFLW_F   = MAX( DSFCMAX-D2ELEVTN(ISEQ,1), 0.D0 )
    DARE_F   = D2FLDSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.D0)
    DARE_F   = MAX( DARE_F - D2FLDDPH(ISEQ,1)*D2RIVWTH(ISEQ,1), 0.D0 )   !! remove above river channel area

    DFLW_PRE_F = DSFCMAX_PRE - D2ELEVTN(ISEQ,1)
    DFLW_IMP_F = MAX( (MAX(DFLW_F*DFLW_PRE_F,0._JPRB))**0.5D0, 1.D-6 )

    DARE_PRE_F = D2FLDSTO_PRE(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.D0)
    DARE_PRE_F = MAX( DARE_PRE_F - D2FLDDPH_PRE(ISEQ,1)*D2RIVWTH(ISEQ,1), 1.D-6 )   !! remove above river channel area
    DARE_IMP_F = max( (DARE_F*DARE_PRE_F)**0.5D0, 1.D-6 )

    IF( DFLW_IMP_F>1.D-5 .and. DARE_IMP_F>1.D-5 )THEN 
      DOUT_PRE_F = D2FLDOUT_PRE(ISEQ,1)
      D2FLDOUT(ISEQ,1) = ( DOUT_PRE_F + PGRV*DT*DARE_IMP_F*DSLOPE_F ) &
                        * (1.D0 + PGRV*DT*PMANFLD**2.D0*abs(DOUT_PRE_F)*DFLW_IMP_F**(-4.D0/3.D0)*DARE_IMP_F**(-1.D0) )**(-1.D0)
    ELSE
      D2FLDOUT(ISEQ,1) = 0.D0
    ENDIF

    IF( D2FLDOUT(ISEQ,1)*D2RIVOUT(ISEQ,1)<0.D0 ) D2FLDOUT(ISEQ,1)=0.D0  !! stabilization
    
  ELSE
    ! Kinematic wave, river flow
    DSLOPE = (D2ELEVTN(ISEQ,1)-D2ELEVTN(JSEQ,1)) * D2NXTDST(ISEQ,1)**(-1.)
    DSLOPE = max(DSLOPE,PMINSLP)
    DVEL   = D2RIVMAN(ISEQ,1)**(-1.) * DSLOPE**0.5 * D2RIVDPH(ISEQ,1)**(2./3.)
    DAREA  = D2RIVWTH(ISEQ,1) * D2RIVDPH(ISEQ,1)

    D2RIVVEL(ISEQ,1) = DVEL
    D2RIVOUT(ISEQ,1) = DAREA * DVEL
    D2RIVOUT(ISEQ,1) = MIN(  D2RIVOUT(ISEQ,1), D2RIVSTO(ISEQ,1)/DT )
    !! kinematic wave, floodplain flow
    DSLOPE_F = min( 0.005D0,DSLOPE )    !! set max&min [instead of using weir equation for efficiency]
    DVEL_F  = PMANFLD**(-1.D0) * DSLOPE_F**0.5D0 * D2FLDDPH(ISEQ,1)**(2.D0/3.D0)
    DARE_F   = D2FLDSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.D0)
    DARE_F   = MAX( DARE_F - D2FLDDPH(ISEQ,1)*D2RIVWTH(ISEQ,1), 0.D0 )   !! remove above river channel area
  
    D2FLDOUT(ISEQ,1) = DARE_F * DVEL_F
    D2FLDOUT(ISEQ,1) = MIN(  D2FLDOUT(ISEQ,1), D2FLDSTO(ISEQ,1)/DT )
  ENDIF 
    
END DO
!$OMP END PARALLEL DO

IF( LBITSAFE )THEN  !! Avoid OMP ATOMIC for bit identical simulation
  DO ISEQ=1, NSEQRIV                                                    !! for normal cells
    JSEQ=I1NEXT(ISEQ) ! next cell's pixel
    OUT_R1 = max(  D2RIVOUT(ISEQ,1),0.D0 )
    OUT_R2 = max( -D2RIVOUT(ISEQ,1),0.D0 )
    OUT_F1 = max(  D2FLDOUT(ISEQ,1),0.D0 )
    OUT_F2 = max( -D2FLDOUT(ISEQ,1),0.D0 )
    DIUP=(OUT_R1+OUT_F1)*DT
    DIDW=(OUT_R2+OUT_F2)*DT
    D2STOOUT(ISEQ,1) = D2STOOUT(ISEQ,1) + DIUP
    D2STOOUT(JSEQ,1) = D2STOOUT(JSEQ,1) + DIDW
  END DO
ELSE
!$OMP PARALLEL DO
  DO ISEQ=1, NSEQRIV                                                    !! for normal cells
    JSEQ=I1NEXT(ISEQ) ! next cell's pixel
    OUT_R1 = max(  D2RIVOUT(ISEQ,1),0.D0 )
    OUT_R2 = max( -D2RIVOUT(ISEQ,1),0.D0 )
    OUT_F1 = max(  D2FLDOUT(ISEQ,1),0.D0 )
    OUT_F2 = max( -D2FLDOUT(ISEQ,1),0.D0 )
    DIUP=(OUT_R1+OUT_F1)*DT
    DIDW=(OUT_R2+OUT_F2)*DT
!$OMP ATOMIC
    D2STOOUT(ISEQ,1) = D2STOOUT(ISEQ,1) + DIUP 
!$OMP ATOMIC
    D2STOOUT(JSEQ,1) = D2STOOUT(JSEQ,1) + DIDW 
  END DO
!$OMP END PARALLEL DO
ENDIF


!$OMP PARALLEL DO                                                     !! for river mouth grids
DO ISEQ=NSEQRIV+1, NSEQALL
  IF (I2MASK(ISEQ,1) == 0 ) THEN

    DSLOPE = ( D2SFCELV(ISEQ,1) - D2DWNELV(ISEQ,1) ) * PDSTMTH ** (-1.D0)
    DSLOPE_F = MAX( -0.005D0, min( 0.005D0,DSLOPE ))    !! set max&min [instead of using weir equation for efficiency]
  !=== river mouth flow ===

    DFLW   = D2RIVDPH(ISEQ,1)
    DAREA  = D2RIVWTH(ISEQ,1) * DFLW

    DFLW_PRE=D2RIVDPH_PRE(ISEQ,1)
    DFLW_IMP=MAX( (DFLW*DFLW_PRE)**0.5D0, 1.D-6 )                                    !! semi implicit flow depth

    IF( DFLW_IMP>1.D-5 .and. DAREA>1.D-5 )THEN 
      DOUT_PRE = D2RIVOUT_PRE(ISEQ,1) * D2RIVWTH(ISEQ,1)**(-1.D0)
      D2RIVOUT(ISEQ,1) = D2RIVWTH(ISEQ,1) * ( DOUT_PRE + PGRV*DT*DFLW_IMP*DSLOPE ) &
                               * ( 1.D0 + PGRV*DT*D2RIVMAN(ISEQ,1)**2.D0*abs(DOUT_PRE)*DFLW_IMP**(-7.D0/3.D0) )**(-1.D0)
      D2RIVVEL(ISEQ,1) = D2RIVOUT(ISEQ,1) * DAREA**(-1.D0)
    ELSE
      D2RIVOUT(ISEQ,1) = 0.D0
      D2RIVVEL(ISEQ,1) = 0.D0
    ENDIF

  !=== floodplain mouth flow ===

    DFLW_F   = D2SFCELV(ISEQ,1)-D2ELEVTN(ISEQ,1)

    DARE_F   = D2FLDSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.D0)
    DARE_F   = MAX( DARE_F - D2FLDDPH(ISEQ,1)*D2RIVWTH(ISEQ,1), 0.D0 )   !! remove above river channel area

    DFLW_PRE_F = D2SFCELV_PRE(ISEQ,1)-D2ELEVTN(ISEQ,1)
    DFLW_IMP_F = MAX( (MAX(DFLW_F*DFLW_PRE_F,0._JPRB))**0.5D0, 1.D-6 )

    DARE_PRE_F = D2FLDSTO_PRE(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.D0)
    DARE_PRE_F = MAX( DARE_PRE_F - D2FLDDPH_PRE(ISEQ,1)*D2RIVWTH(ISEQ,1), 1.D-6 )   !! remove above river channel area
    DARE_IMP_F = max( (DARE_F*DARE_PRE_F)**0.5D0, 1.D-6 )

    IF( DFLW_IMP_F>1.D-5 .and. DARE_IMP_F>1.D-5 )THEN 
      DOUT_PRE_F = D2FLDOUT_PRE(ISEQ,1)
      D2FLDOUT(ISEQ,1) = ( DOUT_PRE_F + PGRV*DT*DARE_IMP_F*DSLOPE_F ) &
                        * (1.D0 + PGRV*DT*PMANFLD**2.D0*abs(DOUT_PRE_F)*DFLW_IMP_F**(-4.D0/3.D0)*DARE_IMP_F**(-1.D0) )**(-1.D0)
    ELSE
      D2FLDOUT(ISEQ,1) = 0.D0
    ENDIF

    IF( D2FLDOUT(ISEQ,1)*D2RIVOUT(ISEQ,1)<0.D0 ) D2FLDOUT(ISEQ,1)=0.D0  !! stabilization
  
  ELSE
    ! Kinematic approach, river channel flow
    DSLOPE = PMINSLP
    DVEL   = D2RIVMAN(ISEQ,1)**(-1.) * DSLOPE**0.5 * D2RIVDPH(ISEQ,1)**(2./3.)
    DAREA  = D2RIVWTH(ISEQ,1) * D2RIVDPH(ISEQ,1)

    D2RIVVEL(ISEQ,1) = DVEL
    D2RIVOUT(ISEQ,1) = DAREA * DVEL
    D2RIVOUT(ISEQ,1) = MIN(  D2RIVOUT(ISEQ,1), D2RIVSTO(ISEQ,1)/DT )

    !! kinematic wave, floodplain flow
    DSLOPE_F = min( 0.005D0,DSLOPE )    !! set max&min [instead of using weir equation for efficiency]
    DVEL_F  = PMANFLD**(-1.D0) * DSLOPE_F**0.5D0 * D2FLDDPH(ISEQ,1)**(2.D0/3.D0)
    DARE_F   = D2FLDSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.D0)
    DARE_F   = MAX( DARE_F - D2FLDDPH(ISEQ,1)*D2RIVWTH(ISEQ,1), 0.D0 )   !! remove above river channel area
  
    D2FLDOUT(ISEQ,1) = DARE_F * DVEL_F
    D2FLDOUT(ISEQ,1) = MIN(  D2FLDOUT(ISEQ,1), D2FLDSTO(ISEQ,1)/DT )
  ENDIF 
 
!=== check outflow ===
  OUT_R1 = max(  D2RIVOUT(ISEQ,1),0.D0 )
  OUT_F1 = max(  D2FLDOUT(ISEQ,1),0.D0 )
  D2STOOUT(ISEQ,1) = D2STOOUT(ISEQ,1) + OUT_R1*DT + OUT_F1*DT
END DO
!$OMP END PARALLEL DO


!! === ourflow modification to avoid negative storage ======

!$OMP PARALLEL DO                                                     !! outflow correcttion if total outflow > storage
DO ISEQ=1, NSEQALL
  IF ( D2STOOUT(ISEQ,1) > 1.D-8 .AND. I2MASK(ISEQ,1) == 0 ) THEN
    D2RATE(ISEQ,1)   = min( (D2RIVSTO(ISEQ,1)+D2FLDSTO(ISEQ,1)) * D2STOOUT(ISEQ,1)**(-1.D0), 1.D0 )
  ENDIF
END DO
!$OMP END PARALLEL DO

IF( LBITSAFE )THEN  !! Aboid OMP ATOMIC for bit identical simulation
  DO ISEQ=1, NSEQRIV ! for normal pixels
    JSEQ=I1NEXT(ISEQ)
    IF( D2RIVOUT(ISEQ,1) >= 0.D0 )THEN
      D2RIVOUT(ISEQ,1) = D2RIVOUT(ISEQ,1)*D2RATE(ISEQ,1)
      D2FLDOUT(ISEQ,1) = D2FLDOUT(ISEQ,1)*D2RATE(ISEQ,1)
    ELSE
      D2RIVOUT(ISEQ,1) = D2RIVOUT(ISEQ,1)*D2RATE(JSEQ,1)
      D2FLDOUT(ISEQ,1) = D2FLDOUT(ISEQ,1)*D2RATE(JSEQ,1)
    ENDIF
    D2RIVINF(JSEQ,1) = D2RIVINF(JSEQ,1) + D2RIVOUT(ISEQ,1)             !! total inflow to a grid (from upstream)
    D2FLDINF(JSEQ,1) = D2FLDINF(JSEQ,1) + D2FLDOUT(ISEQ,1)

    D2RIVOUT_PRE(ISEQ,1)=D2RIVOUT(ISEQ,1)                              !! save outflow (t)
    D2RIVDPH_PRE(ISEQ,1)=D2RIVDPH(ISEQ,1)                              !! save depth   (t)
    D2FLDOUT_PRE(ISEQ,1)=D2FLDOUT(ISEQ,1)                              !! save outflow (t)
    D2FLDSTO_PRE(ISEQ,1)=D2FLDSTO(ISEQ,1)
  END DO
ELSE
!$OMP PARALLEL DO
  DO ISEQ=1, NSEQRIV ! for normal pixels
    JSEQ=I1NEXT(ISEQ)
    IF( D2RIVOUT(ISEQ,1) >= 0.D0 )THEN
      D2RIVOUT(ISEQ,1) = D2RIVOUT(ISEQ,1)*D2RATE(ISEQ,1)
      D2FLDOUT(ISEQ,1) = D2FLDOUT(ISEQ,1)*D2RATE(ISEQ,1)
    ELSE
      D2RIVOUT(ISEQ,1) = D2RIVOUT(ISEQ,1)*D2RATE(JSEQ,1)
      D2FLDOUT(ISEQ,1) = D2FLDOUT(ISEQ,1)*D2RATE(JSEQ,1)
    ENDIF
!$OMP ATOMIC
    D2RIVINF(JSEQ,1) = D2RIVINF(JSEQ,1) + D2RIVOUT(ISEQ,1)             !! total inflow to a grid (from upstream)
!$OMP ATOMIC
    D2FLDINF(JSEQ,1) = D2FLDINF(JSEQ,1) + D2FLDOUT(ISEQ,1)

    D2RIVOUT_PRE(ISEQ,1)=D2RIVOUT(ISEQ,1)                              !! save outflow (t)
    D2RIVDPH_PRE(ISEQ,1)=D2RIVDPH(ISEQ,1)                              !! save depth   (t)
    D2FLDOUT_PRE(ISEQ,1)=D2FLDOUT(ISEQ,1)                              !! save outflow (t)
    D2FLDSTO_PRE(ISEQ,1)=D2FLDSTO(ISEQ,1)
  END DO
!$OMP END PARALLEL DO
ENDIF

!$OMP PARALLEL DO
DO ISEQ=NSEQRIV+1, NSEQALL ! for river mouth
  D2RIVOUT(ISEQ,1) = D2RIVOUT(ISEQ,1)*D2RATE(ISEQ,1)
  D2FLDOUT(ISEQ,1) = D2FLDOUT(ISEQ,1)*D2RATE(ISEQ,1)
END DO
!$OMP END PARALLEL DO

D2RIVOUT_PRE(:,1)=D2RIVOUT(:,1)                              !! save outflow (t)
D2RIVDPH_PRE(:,1)=D2RIVDPH(:,1)                              !! save depth   (t)
D2FLDOUT_PRE(:,1)=D2FLDOUT(:,1)                              !! save outflow (t)
D2FLDSTO_PRE(:,1)=D2FLDSTO(:,1)

END SUBROUTINE CMF_CALC_OUTFLW_KINEMIX
!####################################################################





!####################################################################
SUBROUTINE CMF_CALC_OUTPRE
! to Calculate discharge, diffusive wave, initialization for storage only restart
USE PARKIND1,           ONLY: JPIM, JPRB
USE YOS_CMF_INPUT,      ONLY: LFLDOUT,  LPTHOUT
USE YOS_CMF_INPUT,      ONLY: PMANFLD,  PDSTMTH
USE YOS_CMF_MAP,        ONLY: I1NEXT,   NSEQALL,  NSEQRIV,  NPTHOUT
USE YOS_CMF_MAP,        ONLY: D2RIVELV, D2ELEVTN, D2NXTDST, D2RIVWTH, D2RIVLEN, D2RIVMAN, D2DWNELV
USE YOS_CMF_MAP,        ONLY: NPTHOUT,  NPTHLEV, PTH_UPST, PTH_DOWN, PTH_DST, PTH_ELV, PTH_WTH, PTH_MAN
USE YOS_CMF_PROG,       ONLY: D2RIVOUT_PRE, D2FLDOUT_PRE, D1PTHFLW_PRE ,D2RIVDPH_PRE         !! output
USE YOS_CMF_PROG,       ONLY: D2FLDSTO                                                       !! input
USE YOS_CMF_DIAG,       ONLY: D2RIVDPH, D2SFCELV, D2FLDDPH                                   !! input
IMPLICIT NONE
!$ SAVE
      INTEGER(KIND=JPIM) ::  ISEQ, JSEQ, IPTH, ILEV, ISEQP, JSEQP
      REAL(KIND=JPRB)    ::  DSFCMAX, DSLOPE, DAREA, DFLW, DSLOPE_F, DARE_F, DFLW_F
!$OMP THREADPRIVATE         (DSFCMAX, DSLOPE, DAREA, DFLW, DSLOPE_F, DARE_F, DFLW_F, JSEQ, ILEV, ISEQP, JSEQP)
!================================================

!$OMP PARALLEL DO
DO ISEQ=1, NSEQALL
  D2SFCELV(ISEQ,1)     = D2RIVELV(ISEQ,1) + D2RIVDPH(ISEQ,1)
  D2RIVDPH_PRE(ISEQ,1) = D2RIVDPH(ISEQ,1)                           !! bugfix v362
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
DO ISEQ=1, NSEQRIV                                                    !! for normal cells
  JSEQ=I1NEXT(ISEQ)

  DSFCMAX    =MAX( D2SFCELV(ISEQ,1),    D2SFCELV(JSEQ,1) )
  DSLOPE = ( D2SFCELV(ISEQ,1)-D2SFCELV(JSEQ,1) ) * D2NXTDST(ISEQ,1)**(-1.D0)
  DSLOPE_F = MAX( -0.005D0, min( 0.005D0,DSLOPE ))    !! set max&min [instead of using weir equation for efficiency]

!=== River Flow ===
  DFLW  = DSFCMAX - D2RIVELV(ISEQ,1)
  DAREA = D2RIVWTH(ISEQ,1) * DFLW

  IF( DAREA>1.D-5 )THEN
    D2RIVOUT_PRE(ISEQ,1) = DAREA * ( D2RIVMAN(ISEQ,1)**(-1.D0) * DFLW**(2.D0/3.D0) * abs(DSLOPE)**(0.5D0) )
    IF( DSLOPE<0.D0 ) D2RIVOUT_PRE(ISEQ,1)=-D2RIVOUT_PRE(ISEQ,1)
  ELSE
    D2RIVOUT_PRE(ISEQ,1) = 0.D0
  ENDIF

!=== Floodplain Flow ===
  DFLW_F   = MAX( DSFCMAX-D2ELEVTN(ISEQ,1), 0.D0 )
  DARE_F   = D2FLDSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.D0)
  DARE_F   = MAX( DARE_F - D2FLDDPH(ISEQ,1)*D2RIVWTH(ISEQ,1), 0.D0 )   !! remove above river channel area

  IF( DARE_F>1.D-5 )THEN 
    D2FLDOUT_PRE(ISEQ,1) = DARE_F * ( PMANFLD**(-1.D0) * DFLW_F**(2.D0/3.D0) * abs(DSLOPE_F)**(0.5D0) )
    IF( DSLOPE_F<0.D0 ) D2FLDOUT_PRE(ISEQ,1)=-D2FLDOUT_PRE(ISEQ,1)
  ELSE
    D2FLDOUT_PRE(ISEQ,1) = 0.D0
  ENDIF
END DO
!$OMP END PARALLEL DO


!$OMP PARALLEL DO                                                     !! for river mouth grids
DO ISEQ=NSEQRIV+1, NSEQALL
  DSLOPE = ( D2SFCELV(ISEQ,1)-D2DWNELV(ISEQ,1) ) * PDSTMTH**(-1.D0)
  DSLOPE_F = MAX( -0.005D0, min( 0.005D0,DSLOPE ))    !! set max&min [instead of using weir equation for efficiency]

!=== river mouth flow ===
  DFLW   = D2RIVDPH(ISEQ,1)
  DAREA  = D2RIVWTH(ISEQ,1) * DFLW
  IF( DAREA>1.D-5 )THEN
    D2RIVOUT_PRE(ISEQ,1) = DAREA * ( D2RIVMAN(ISEQ,1)**(-1.D0) * DFLW**(2.D0/3.D0) * abs(DSLOPE)**(0.5D0) )
    IF( DSLOPE<0.D0 ) D2RIVOUT_PRE(ISEQ,1)=-D2RIVOUT_PRE(ISEQ,1)
  ELSE
    D2RIVOUT_PRE(ISEQ,1) = 0.D0
  ENDIF
!=== floodplain mouth flow ===
  DFLW_F   = D2SFCELV(ISEQ,1)-D2ELEVTN(ISEQ,1)
  DARE_F   = D2FLDSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.D0)
  DARE_F   = MAX( DARE_F - D2FLDDPH(ISEQ,1)*D2RIVWTH(ISEQ,1), 0.D0 )   !! remove above river channel area
  IF( DARE_F>1.D-5 )THEN
    D2FLDOUT_PRE(ISEQ,1) = DARE_F * ( PMANFLD**(-1.D0) * DFLW_F**(2.D0/3.D0) * abs(DSLOPE_F)**(0.5D0) )
    IF( DSLOPE_F<0.D0 ) D2FLDOUT_PRE(ISEQ,1)=-D2FLDOUT_PRE(ISEQ,1)
  ELSE
    D2FLDOUT_PRE(ISEQ,1) = 0.D0
  ENDIF
END DO
!$OMP END PARALLEL DO


IF( LPTHOUT )THEN
 !$OMP PARALLEL DO
 DO IPTH=1, NPTHOUT  
  ISEQP=PTH_UPST(IPTH)
  JSEQP=PTH_DOWN(IPTH)

  DSLOPE  = (D2SFCELV(ISEQP,1)-D2SFCELV(JSEQP,1)) * PTH_DST(IPTH)**(-1.D0)
  DO ILEV=1, NPTHLEV
    DFLW = MAX(D2SFCELV(ISEQP,1),D2SFCELV(JSEQP,1)) - PTH_ELV(IPTH,ILEV) 
    DFLW = MAX(DFLW,0.D0)

    IF( DFLW>1.D-5 )THEN
      D1PTHFLW_PRE(IPTH,ILEV) = PTH_WTH(IPTH,ILEV) * DFLW * ( PTH_MAN(ILEV)**(-1.D0) * DFLW**(2.D0/3.D0) * abs(DSLOPE)**(0.5D0) )
      IF( DSLOPE<0.D0 ) D1PTHFLW_PRE(IPTH,ILEV)=-D1PTHFLW_PRE(IPTH,ILEV)
    ELSE
      D1PTHFLW_PRE(IPTH,ILEV) = 0.D0
    ENDIF
  END DO
 END DO
 !$OMP END PARALLEL DO
ENDIF

!! when high-water flow is not defined
IF( .not. LFLDOUT )THEN
  D2FLDOUT_PRE(:,:)=0.D0
ENDIF

END SUBROUTINE CMF_CALC_OUTPRE
!####################################################################

END MODULE CMF_OPT_OUTFLW_MOD
