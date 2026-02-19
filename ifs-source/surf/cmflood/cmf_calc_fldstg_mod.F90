MODULE CMF_CALC_FLDSTG_MOD
!==========================================================
!* PURPOSE: calculate river and floodplain staging
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
! -- CMF_CALC_FLDSTG
! --
! --
!####################################################################
SUBROUTINE CMF_CALC_FLDSTG
USE PARKIND1,           ONLY: JPIM, JPRB
USE YOS_CMF_INPUT,      ONLY: NLFP
USE YOS_CMF_MAP,        ONLY: NSEQALL
USE YOS_CMF_MAP,        ONLY: D2GRAREA, D2RIVLEN, D2RIVWTH, D2RIVELV, D2RIVSTOMAX, D2FLDSTOMAX, D2FLDGRD, DFRCINC
USE YOS_CMF_PROG,       ONLY: D2RIVSTO, D2FLDSTO
USE YOS_CMF_DIAG,       ONLY: D2RIVDPH, D2FLDDPH, D2FLDFRC, D2FLDARE, D2SFCELV, DGLBRIVSTO, DGLBFLDSTO, DGLBFLDARE
IMPLICIT NONE

!*** LOCAL
!$ SAVE
INTEGER(KIND=JPIM)         :: ISEQ, I
REAL(KIND=JPRB)            :: DSTOALL, DSTONOW, DSTOPRE, DWTHNOW, DWTHPRE, DDPHPRE, DWTHINC
!$OMP THREADPRIVATE        (I,DSTOALL, DSTONOW, DSTOPRE, DWTHNOW, DWTHPRE, DDPHPRE, DWTHINC)
!================================================
DGLBRIVSTO=0.D0
DGLBFLDSTO=0.D0
DGLBFLDARE=0.D0

!$OMP PARALLEL DO REDUCTION(+:DGLBRIVSTO,DGLBFLDSTO,DGLBFLDARE)
DO ISEQ=1, NSEQALL
!
  DSTOALL = D2RIVSTO(ISEQ,1) + D2FLDSTO(ISEQ,1)
  IF( DSTOALL > D2RIVSTOMAX(ISEQ,1) )THEN
    I=1
    DSTOPRE = D2RIVSTOMAX(ISEQ,1)
    DWTHPRE = D2RIVWTH(ISEQ,1)
    DDPHPRE = 0.D0
    DWTHINC = D2GRAREA(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * DFRCINC
    DO WHILE( DSTOALL > D2FLDSTOMAX(ISEQ,1,I) .AND. I<=NLFP)
      DSTOPRE = D2FLDSTOMAX(ISEQ,1,I)
      DWTHPRE = DWTHPRE + DWTHINC
      DDPHPRE = DDPHPRE + D2FLDGRD(ISEQ,1,I) * DWTHINC
      I=I+1
      IF( I>NLFP ) EXIT
    END DO
    IF( I>NLFP )THEN
      DSTONOW = DSTOALL - DSTOPRE
      DWTHNOW = 0.D0
      D2FLDDPH(ISEQ,1) = DDPHPRE + DSTONOW * DWTHPRE**(-1.) * D2RIVLEN(ISEQ,1)**(-1.)
    ELSE
      DSTONOW =  DSTOALL - DSTOPRE
      DWTHNOW = -DWTHPRE + ( DWTHPRE**2. + 2.D0 * DSTONOW * D2RIVLEN(ISEQ,1)**(-1.) * D2FLDGRD(ISEQ,1,I)**(-1.) )**0.5
      D2FLDDPH(ISEQ,1) = DDPHPRE + D2FLDGRD(ISEQ,1,I) * DWTHNOW
    ENDIF
    D2RIVSTO(ISEQ,1) = D2RIVSTOMAX(ISEQ,1) + D2RIVLEN(ISEQ,1) * D2RIVWTH(ISEQ,1) * D2FLDDPH(ISEQ,1)
    D2RIVDPH(ISEQ,1) = D2RIVSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * D2RIVWTH(ISEQ,1)**(-1.)
!
    D2FLDSTO(ISEQ,1) = DSTOALL - D2RIVSTO(ISEQ,1)
    D2FLDSTO(ISEQ,1) = MAX( D2FLDSTO(ISEQ,1), 0.D0 )
    D2FLDFRC(ISEQ,1) = (-D2RIVWTH(ISEQ,1) + DWTHPRE + DWTHNOW ) * (DWTHINC*10.D0)**(-1.)
    D2FLDFRC(ISEQ,1) = MAX( D2FLDFRC(ISEQ,1),0.D0)
    D2FLDFRC(ISEQ,1) = MIN( D2FLDFRC(ISEQ,1),1.D0)
    D2FLDARE(ISEQ,1) = D2GRAREA(ISEQ,1)*D2FLDFRC(ISEQ,1)
  ELSE
    D2RIVSTO(ISEQ,1) = DSTOALL
    D2RIVDPH(ISEQ,1) = DSTOALL * D2RIVLEN(ISEQ,1)**(-1.) * D2RIVWTH(ISEQ,1)**(-1.)
    D2RIVDPH(ISEQ,1) = MAX( D2RIVDPH(ISEQ,1), 0.D0 )
    D2FLDSTO(ISEQ,1) = 0.D0
    D2FLDDPH(ISEQ,1) = 0.D0
    D2FLDFRC(ISEQ,1) = 0.D0
    D2FLDARE(ISEQ,1) = 0.D0
  ENDIF
  D2SFCELV(ISEQ,1)     = D2RIVELV(ISEQ,1) + D2RIVDPH(ISEQ,1)

  DGLBRIVSTO      = DGLBRIVSTO + D2RIVSTO(ISEQ,1)
  DGLBFLDSTO      = DGLBFLDSTO + D2FLDSTO(ISEQ,1)
  DGLBFLDARE      = DGLBFLDARE + D2FLDARE(ISEQ,1)

END DO
!$OMP END PARALLEL DO

END SUBROUTINE CMF_CALC_FLDSTG
!####################################################################

END MODULE CMF_CALC_FLDSTG_MOD
