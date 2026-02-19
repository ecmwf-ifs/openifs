MODULE CMF_CALC_PTHOUT_MOD
!==========================================================
!* PURPOSE: subroutine for bifurcation channel flow
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
! -- CMF_CALC_PTHOUT
! --
! --
!####################################################################
SUBROUTINE CMF_CALC_PTHOUT
USE PARKIND1,           ONLY: JPIM, JPRB
USE YOS_CMF_INPUT,      ONLY: DT, PGRV, LBITSAFE
USE YOS_CMF_MAP,        ONLY: NSEQMAX, NPTHOUT, NPTHLEV, PTH_UPST, PTH_DOWN, PTH_DST, &
                            & PTH_ELV, PTH_WTH, PTH_MAN, I2MASK
USE YOS_CMF_MAP,        ONLY: D2RIVELV
USE YOS_CMF_PROG,       ONLY: D2RIVSTO, D2FLDSTO, D1PTHFLW, D2RIVOUT, D2FLDOUT
USE YOS_CMF_PROG,       ONLY: D1PTHFLW_PRE, D2RIVDPH_PRE
USE YOS_CMF_DIAG,       ONLY: D2PTHOUT, D2PTHINF, D2RIVINF, D2FLDINF, D2SFCELV
IMPLICIT NONE
!*** Local
REAL(KIND=JPRB)    ::  D2SFCELV_PRE(NSEQMAX,1)                  !! water surface elev (t-1)[m] (for floodplain)
REAL(KIND=JPRB)    ::  D2RATE(NSEQMAX,1)                        !! outflow correction

INTEGER(KIND=JPIM),SAVE ::  IPTH, ILEV, ISEQ, ISEQP, JSEQP
REAL(KIND=JPRB),SAVE    ::  DSLOPE, DFLW, DOUT_PRE, DFLW_PRE, DFLW_IMP, DSTO_TMP
!$OMP THREADPRIVATE         (DSLOPE, DFLW, DOUT_PRE, DFLW_PRE, DFLW_IMP, DSTO_TMP, ILEV, ISEQP, JSEQP)

!================================================
!$OMP PARALLEL DO
DO ISEQ=1, NSEQMAX
  D2SFCELV_PRE(ISEQ,1) = D2RIVELV(ISEQ,1)+D2RIVDPH_PRE(ISEQ,1)
  D2PTHOUT(ISEQ,1) = 0.D0
  D2PTHINF(ISEQ,1) = 0.D0
  D2RATE(ISEQ,1)=-999.
END DO
!$OMP END PARALLEL DO

D1PTHFLW(:,:) = 0._JPRB 
!$OMP PARALLEL DO
DO IPTH=1, NPTHOUT  
  ISEQP=PTH_UPST(IPTH)
  JSEQP=PTH_DOWN(IPTH)
  !! Avoid calculation outside of domain
  IF (ISEQP == 0 .OR. JSEQP== 0 ) CYCLE
  IF (I2MASK(ISEQP,1) == 1 .OR. I2MASK(JSEQP,1) == 1 ) CYCLE
  
  DSLOPE  = (D2SFCELV(ISEQP,1)-D2SFCELV(JSEQP,1)) * PTH_DST(IPTH)**(-1.D0)
  DSLOPE = max(-0.005D0,min(0.005D0,DSLOPE))                                    !! v390 stabilization

  DO ILEV=1, NPTHLEV

    DFLW = MAX(D2SFCELV(ISEQP,1),D2SFCELV(JSEQP,1)) - PTH_ELV(IPTH,ILEV) 
    DFLW = MAX(DFLW,0.D0)

    DFLW_PRE = MAX(D2SFCELV_PRE(ISEQP,1),D2SFCELV_PRE(JSEQP,1)) - PTH_ELV(IPTH,ILEV)
    DFLW_PRE = MAX(DFLW_PRE,0.D0)

    DFLW_IMP = (DFLW*DFLW_PRE)**0.5D0                                       !! semi implicit flow depth
    IF( DFLW_IMP<=0.D0 ) DFLW_IMP=DFLW

    IF( DFLW_IMP>1.D-5 )THEN                         !! local inertial equation, see [Bates et al., 2010, J.Hydrol.]
      DOUT_PRE = D1PTHFLW_PRE(IPTH,ILEV) * PTH_WTH(IPTH,ILEV)**(-1.D0)                         !! outflow (t-1) [m2/s] (unit width)
      D1PTHFLW(IPTH,ILEV) = PTH_WTH(IPTH,ILEV) * ( DOUT_PRE + PGRV*DT*DFLW_IMP*DSLOPE ) &
                         * ( 1.D0 + PGRV*DT*PTH_MAN(ILEV)**2.D0*abs(DOUT_PRE)*DFLW_IMP**(-7.D0/3.D0) )**(-1.D0)
    ELSE
      D1PTHFLW(IPTH,ILEV) = 0.D0
    ENDIF
  END DO
END DO
!$OMP END PARALLEL DO

IF( LBITSAFE )THEN  !! Aboid OMP ATOMIC for bit identical simulation
  DO IPTH=1, NPTHOUT  
    ISEQP=PTH_UPST(IPTH)
    JSEQP=PTH_DOWN(IPTH)
    DO ILEV=1, NPTHLEV
      IF( D1PTHFLW(IPTH,ILEV) >= 0.D0 )THEN                                  !! total outflow from each grid
        D2PTHOUT(ISEQP,1) = D2PTHOUT(ISEQP,1) + D1PTHFLW(IPTH,ILEV)
      ELSE
        D2PTHOUT(JSEQP,1) = D2PTHOUT(JSEQP,1) - D1PTHFLW(IPTH,ILEV)
      ENDIF
    END DO
  END DO
ELSE
!$OMP PARALLEL DO
  DO IPTH=1, NPTHOUT  
    ISEQP=PTH_UPST(IPTH)
    JSEQP=PTH_DOWN(IPTH)
    DO ILEV=1, NPTHLEV
      IF( D1PTHFLW(IPTH,ILEV) >= 0.D0 )THEN                                  !! total outflow from each grid
!$OMP ATOMIC
        D2PTHOUT(ISEQP,1) = D2PTHOUT(ISEQP,1) + D1PTHFLW(IPTH,ILEV)
      ELSE
!$OMP ATOMIC
        D2PTHOUT(JSEQP,1) = D2PTHOUT(JSEQP,1) - D1PTHFLW(IPTH,ILEV)
      ENDIF
    END DO
  END DO
!$OMP END PARALLEL DO
ENDIF

!$OMP PARALLEL DO                                              !! calculate total outflow from a grid
DO ISEQ=1, NSEQMAX
  IF( D2PTHOUT(ISEQ,1) > 1.D-10 )THEN
    DSTO_TMP = ( D2RIVSTO(ISEQ,1)+D2FLDSTO(ISEQ,1) ) &
                  - D2RIVOUT(ISEQ,1)*DT + D2RIVINF(ISEQ,1)*DT - D2FLDOUT(ISEQ,1)*DT + D2FLDINF(ISEQ,1)*DT
    D2RATE(ISEQ,1) = MIN( DSTO_TMP * (D2PTHOUT(ISEQ,1)*DT)**(-1.), 1.D0 )
  ELSE
    D2RATE(ISEQ,1) = 1.D0
  ENDIF
  D2PTHOUT(ISEQ,1) = D2PTHOUT(ISEQ,1) * D2RATE(ISEQ,1)
END DO
!$OMP END PARALLEL DO


IF( LBITSAFE )THEN  !! Aboid OMP ATOMIC for bit identical simulation
  DO IPTH=1, NPTHOUT
    ISEQP=PTH_UPST(IPTH)
    JSEQP=PTH_DOWN(IPTH)
    DO ILEV=1, NPTHLEV
      IF( D1PTHFLW(IPTH,ILEV) >= 0.D0 )THEN
        D1PTHFLW(IPTH,ILEV) = D1PTHFLW(IPTH,ILEV)*D2RATE(ISEQP,1)
        D2PTHINF(JSEQP,1) = D2PTHINF(JSEQP,1) + D1PTHFLW(IPTH,ILEV)             !! total inflow [m3/s] (from upstream)
      ELSE
        D1PTHFLW(IPTH,ILEV) = D1PTHFLW(IPTH,ILEV)*D2RATE(JSEQP,1)
        D2PTHINF(ISEQP,1) = D2PTHINF(ISEQP,1) - D1PTHFLW(IPTH,ILEV)             !! total inflow [m3/s] (from upstream)
      ENDIF
      D1PTHFLW_PRE(IPTH,ILEV)=D1PTHFLW(IPTH,ILEV)
    END DO
  END DO
ELSE
!$OMP PARALLEL DO   
  DO IPTH=1, NPTHOUT
    ISEQP=PTH_UPST(IPTH)
    JSEQP=PTH_DOWN(IPTH)
    DO ILEV=1, NPTHLEV
      IF( D1PTHFLW(IPTH,ILEV) >= 0.D0 )THEN
        D1PTHFLW(IPTH,ILEV) = D1PTHFLW(IPTH,ILEV)*D2RATE(ISEQP,1)
!$OMP ATOMIC
        D2PTHINF(JSEQP,1) = D2PTHINF(JSEQP,1) + D1PTHFLW(IPTH,ILEV)             !! total inflow [m3/s] (from upstream)
      ELSE
        D1PTHFLW(IPTH,ILEV) = D1PTHFLW(IPTH,ILEV)*D2RATE(JSEQP,1)
!$OMP ATOMIC
        D2PTHINF(ISEQP,1) = D2PTHINF(ISEQP,1) - D1PTHFLW(IPTH,ILEV)             !! total inflow [m3/s] (from upstream)
      ENDIF
      D1PTHFLW_PRE(IPTH,ILEV)=D1PTHFLW(IPTH,ILEV)
    END DO
  END DO
!$OMP END PARALLEL DO
ENDIF

END SUBROUTINE CMF_CALC_PTHOUT
!####################################################################

END MODULE CMF_CALC_PTHOUT_MOD
