! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

MODULE STOPH_MIX

USE PARKIND1 , ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------

!*    control parameters for stochastic physics

!      Variables for forced S.V.s
!      LFORCENL   : add perturbations to model tendencies during forecast
!                   period NFORCESTART and NFORCEEND (in hours)

!      Variables for stochastic physics: option PDT
!      NFRSTOPH_SPBS   : time frequency (number of time steps) between
!                        each call to spectral backscatter
!      NFRSTOPH_SPBS_PAT   : time frequency (number of time steps) between
!                        each call to spectral pattern update
!      NFRSTOPH_VC     : time frequency (number of time steps) between
!                        each call to vorticity confinement
!      NSTOCHOPT  : =1 dissipation rate based on massflux-formulation;
!      NSTOCHOPT  : =2 dissipation rate based on updraft-formulation
!      NSTOCHOPT  : =3 dissipation rate based on revised updraft-formulation
!      SQRTCORR   : sqrt of vertical correlation matrix for SPBS perturbations
!      Variables for stochastic physics: common to CASBS and SPBS
!      GPSTREAM        : Gridpoint field for streamfunction perturbations
!      GPTOTDISS       : Gridpoint field for total dissipation rates
!      GPVORTGRAD      : gridpoint field for grad(vorticity) (used in numerical
!                        dissipation rate calculation and vorticity confinement
!      RSMOOTH         : Spectral filter for dissipation rate field
!      ALPHA_DEEP_CONV : Entrainment cloud fraction for deep convection
!      ALPHA_SHAL_CONV : Entrainment cloud fraction for shallow and mid-level convection
!      LEXTRAFIELDS    : .T. to write extrafields
!      LSPBSNORM       ; .T. to calculate and print out energy input
!      LSPBSDISS       ; .T. to calculate and print out global mean dissipationrate
!      SPDP(:)         : average level thickness in Pa (used for energy input calculation)

!      Variables for stochastic physics: option SPBS
!      SPSTREAM        : Spectral field for streamfunction perturbations
!      SPSTREAM_FORC   : Spectral field for streamfunction forcing field
!      SPVELPOT        : Spectral field for velocity potential pattern field
!      SPVELPOT_FORC   : Spectral field for velocity potential forcing field
!      SPG_AMP         : Noise amplitude of streamfunction forcing 
!      ALPHA_STO       : Autoregressive parameter 
!      ONEMINALPHA_NFRSPBS : (1-Autoregressive parameter) if pattern is only updated every NFRSTOPH_SPBS timestep
!      RATIO_BACKSCAT  : Backscatter ratio
!      RSPBS_TAU       : decorrelation time of pattern
!      LSTOPH_SPBS     : .T. to add spectral backscatter streamfunction perturbations
!      LSTOPH_SPBS_FAST: .T. pattern is only updated every NFRSTOPH_SPBS timesteps 
!                        (bigger steps -> should be equivalent statistically)
!      LSTOPH_SPBS_VORT: .T. use vorticity-based forcing  ansatz
!      LSTOPH_TAPER    : .T. to suppress stream function forcing in the PBL
!      LSTOPH_JBCOR    : .T. to make SPBS perturbations vertically correlated similar to Jb statistics for vorticity
!      LSTOPH_UNCORR   : .T. to make SPBS perturbations vertically uncorrelated 
!      LSTOPH_UNIFORM  : .T. to make SPBS perturbations uniformly distributed, instead 
!                        of Gaussian distributed
!      LSTOPH_RVP      : .T. to use random vertical profiles in SPBS
!      LSTOPH_RVPOLD   : .T. to use old RVP (p-dependency only appropriate for 91 and 62 levels)
!      LSTOPH_VARALPHA : .T. to use scale-dependent decorrelation time in SPBS
!      LSTOPH_GAUSS    : .T. to use gaussian for dissipation smoothing
!      LSPBS_DISSGW    : .T. orographic GWD contributes to dissipation rate estimate
!      LSPBS_DISSNUM   : .T. horizontal diffusion contributes to dissipation rate estimate
!      LSPBS_DISSNUM_CT: .T. use horizontal diffusion in SPBS consistent with hor. diff in model
!      LSPBS_DISSCU    : .T. deep convection contributes to dissipation rate estimate
!      RSMOOTHSCALE    : "decorrelation" scale for smoothing
!      RSIGMA2_EPS     : variance of complex random numbers (cancels out for gaussian distributed random numbers in the AR1 process
!                        needs to be 2/12 when uniform distributed random numbers are used)
!      NIMRAN          : address in a spectral array of (m, n=m) for m in increasing order
!                        ([m=0,n=0], [0,1], [0,2],...,[1,1], [1,2],...,[2,2], [2,3],..,[NSMAX,NSMAX])
!      NSMAXSPBS       : maximum wavenumber forced by SPBS

!      Variables for stochastic physics: option SPBS with random vertical profiles
!      RVP_MULMIN      : parameter for wavemode dependent vertical correlation
!      RVP_MULMAX      : parameter for wavemode dependent vertical correlation
!      RVP_MULEXP      : parameter for wavemode dependent vertical correlation
!      RVP_MULNSMAX    : parameter for wavemode dependent vertical correlation
!      RVP_MUL_A       : parameter for pressure dependent vertical correlation
!      RVP_MUL_B       : parameter for pressure dependent vertical correlation
!      RVP_MUL_C       : parameter for pressure dependent vertical correlation
!      RVP_MUL_D       : parameter for pressure dependent vertical correlation
!      RVP_MUL_1       : parameter for new pressure dependent vertical correlation
!      RVP_MUL_2       : parameter for new pressure dependent vertical correlation
!      RVP_MUL(:)      : Multiplicator for pressure dependent vertical correlation
     
!      Variables for stochastic physics: option CASBS
!      LSTOPH_CASBS    : .T. to add cellular automaton streamfunction perturbations
!      MCELL           : Lat-lon grid holding the 'number of lives' for each CA cell
!      RWGHT           : Weighted and smoothed CA pattern
!      AMAGSTOPH_CASBS : Magnitude of forcing 
!      ADLATSTOPH_CA   : Gridsize of cellular automaton in zonal direction
!      ADLONSTOPH_CA   : Gridsize of cellular automaton in meridional direction
!      RFLUX_DET_CLIP  : upper limit for convective mass flux detrainment rate

!      Variables for vorticity confinement 
!      LVORTCON        : .T. to switch on vorticity confinement
!      VC_CON          :  value of the 'epsilon' parameter in Steinhoff's type 1 vorticity confinement

TYPE :: TSTOPH
REAL(KIND=JPRB),ALLOCATABLE :: RSTOPHCA(:)

REAL(KIND=JPRB),ALLOCATABLE :: SQRTCORR(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SPSTREAM(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SPVELPOT(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SPSTREAM_FORC(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SPVELPOT_FORC(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SPG_AMP(:)
REAL(KIND=JPRB),ALLOCATABLE :: ALPHA_STO(:)
REAL(KIND=JPRB),ALLOCATABLE :: ONEMINALPHA_NFRSPBS(:)
REAL(KIND=JPRB),ALLOCATABLE :: RSMOOTH(:)

REAL(KIND=JPRB),ALLOCATABLE :: GPSTREAM(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: GPVELPOT(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: GPTOTDISS(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: GPTOTDISS_SMOOTH(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: GPVORTGRAD(:,:,:)

INTEGER(KIND=JPIM),ALLOCATABLE :: MCELL(:,:)
REAL(KIND=JPRB),   ALLOCATABLE :: RWGHT(:,:)

INTEGER(KIND=JPIM) :: NFRSTOPH_SPBS
INTEGER(KIND=JPIM) :: NFRSTOPH_SPBS_PAT
INTEGER(KIND=JPIM) :: NFRSTOPH_VC
INTEGER(KIND=JPIM) :: NFORCESTART
INTEGER(KIND=JPIM) :: NFORCEEND
INTEGER(KIND=JPIM) :: NSTOCHOPT

REAL(KIND=JPRB) :: ALPHA_DEEP_CONV
REAL(KIND=JPRB) :: ALPHA_SHAL_CONV
REAL(KIND=JPRB) :: SLDISSFAC
REAL(KIND=JPRB) :: RATIO_BACKSCAT
REAL(KIND=JPRB) :: RSPBS_TAU
REAL(KIND=JPRB) :: RATIO_BACKSCAT_CON2NUM
REAL(KIND=JPRB) :: REXPONENT
REAL(KIND=JPRB) :: VC_CON
REAL(KIND=JPRB) :: RFLUX_DET_CLIP
REAL(KIND=JPRB) :: BIHARM
REAL(KIND=JPRB) :: RSMOOTHSCALE
REAL(KIND=JPRB) :: RSIGMA2_EPS

INTEGER(KIND=JPIM) :: INCUT

REAL(KIND=JPRB) :: AMAGSTOPH_CASBS
REAL(KIND=JPRB) :: ADLATSTOPH_CA
REAL(KIND=JPRB) :: ADLONSTOPH_CA

LOGICAL :: LSTOPH_SPBS,LSTOPH_SPBS_FAST,LEXTRAFIELDS,LSTOPH_JBCOR,LSTOPH_UNCORR,LSTOPH_UNIFORM
LOGICAL :: LSTOPH_SPBS_VORT
LOGICAL :: LSTOPH_RVP,LSTOPH_RVPOLD
LOGICAL :: LSTOPH_TAPER
LOGICAL :: LSTOPH_INI
LOGICAL :: LSTOPH_CASBS
LOGICAL :: LVORTCON
LOGICAL :: LFORCENL
LOGICAL :: LSTOPH_VARALPHA
LOGICAL :: LSTOPH_GAUSS
LOGICAL :: LSPBS_DISSGW, LSPBS_DISSNUM, LSPBS_DISSCU 
LOGICAL :: LSPBS_DISSNUM_CT
LOGICAL :: LSPBSNORM
LOGICAL :: LSPBSDISS
REAL(KIND=JPRB),ALLOCATABLE::   SPDP(:)

INTEGER(KIND=JPIM),ALLOCATABLE:: NIMRAN(:)

INTEGER(KIND=JPIM) :: NSMAXSPBS

REAL(KIND=JPRB) ::    RVP_MULMIN
REAL(KIND=JPRB) ::    RVP_MULMAX
REAL(KIND=JPRB) ::    RVP_MULEXP
REAL(KIND=JPRB) ::    RVP_MULNSMAX
REAL(KIND=JPRB),ALLOCATABLE ::  RVP_MULFACT(:)
REAL(KIND=JPRB) ::    RVP_MUL_A, RVP_MUL_B, RVP_MUL_C, RVP_MUL_D
REAL(KIND=JPRB) ::    RVP_MUL_1, RVP_MUL_2
REAL(KIND=JPRB),ALLOCATABLE ::   RVP_MUL(:)
REAL(KIND=JPRB) ::    TAPER_SIGMATOP, TAPER_SIGMABOT, TAPER0, TAPER1,&
                    & TAPER2, TAPER3
REAL(KIND=JPRB),ALLOCATABLE::   TAPER_FACT(:)

!variables for T-Backscatter

LOGICAL :: LSTOPH_SPBS_T
REAL(KIND=JPRB) :: REXPONENT_T
REAL(KIND=JPRB) :: RATIO_APE2KE
REAL(KIND=JPRB),ALLOCATABLE :: SPTEMP(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SPTEMP_FORC(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: SPG_AMP_T(:)
REAL(KIND=JPRB),ALLOCATABLE :: ALPHA_STO_T(:)
REAL(KIND=JPRB),ALLOCATABLE :: ONEMINALPHA_NFRSPBS_T(:)
REAL(KIND=JPRB),ALLOCATABLE :: GPTEMP(:,:,:)

REAL(KIND=JPRB) ::    RVP_MULMIN_T
REAL(KIND=JPRB) ::    RVP_MULMAX_T
REAL(KIND=JPRB) ::    RVP_MULEXP_T
REAL(KIND=JPRB) ::    RVP_MULNSMAX_T
REAL(KIND=JPRB),ALLOCATABLE ::    RVP_MULFACT_T(:)
REAL(KIND=JPRB) ::    RVP_MUL_A_T, RVP_MUL_B_T, RVP_MUL_C_T, RVP_MUL_D_T
REAL(KIND=JPRB) ::    RVP_MUL_1_T, RVP_MUL_2_T
REAL(KIND=JPRB),ALLOCATABLE ::   RVP_MUL_T(:)
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TSTOPH
!============================================================================

!     ------------------------------------------------------

CONTAINS

SUBROUTINE INITIALIZE_CELLS(YD_RANDOM_STREAMS,KCELL,PWGHT,KIIP,KJJP)
USE YOMRANDOM_STREAMS , ONLY : TRANDOM_STREAMS
USE RANDOM_NUMBERS_MIX, ONLY : UNIFORM_DISTRIBUTION
USE YOMLUN            , ONLY : NULOUT
USE YOMHOOK           , ONLY : LHOOK,   DR_HOOK, JPHOOK 

IMPLICIT NONE
TYPE(TRANDOM_STREAMS) , INTENT(INOUT) :: YD_RANDOM_STREAMS
INTEGER(KIND=JPIM)    , INTENT(IN)    :: KIIP
INTEGER(KIND=JPIM)    , INTENT(IN)    :: KJJP
INTEGER(KIND=JPIM)    , INTENT(INOUT) :: KCELL(KIIP*4,KJJP*4)
REAL(KIND=JPRB)       , INTENT(INOUT) :: PWGHT(KIIP,KJJP)

INTEGER(KIND=JPIM)            :: IC,JC
INTEGER(KIND=JPIM), PARAMETER :: INCLUMPS=200
INTEGER(KIND=JPIM), PARAMETER :: ILIVES=10
INTEGER(KIND=JPIM)            :: I,ICLUMP,JCLUMP
REAL(KIND=JPRB)               :: ZRAND_NOS(2*INCLUMPS+10)   ! define the random numbers
REAL(KIND=JPHOOK)               :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('STOPH_MIX:INITIALIZE_CELLS',0,ZHOOK_HANDLE)

IC=KIIP*4
JC=KJJP*4

CALL UNIFORM_DISTRIBUTION(ZRAND_NOS,YD_RANDOM_STREAMS%STOCHPHYS_CABS)

KCELL=0    ! array initialization
WRITE(NULOUT,*) 'initializing KCELLs'

DO I=1,2*INCLUMPS,2
  ICLUMP= 2 + NINT((IC-3)*ZRAND_NOS(I))
  JCLUMP= 2 + NINT((JC-3)*ZRAND_NOS(I+1))
  KCELL(ICLUMP,JCLUMP)= ILIVES
  KCELL(ICLUMP+1,JCLUMP)= ILIVES
  KCELL(ICLUMP-1,JCLUMP)= ILIVES
  KCELL(ICLUMP,JCLUMP+1)= ILIVES
  KCELL(ICLUMP,JCLUMP-1)= ILIVES
  KCELL(ICLUMP+1,JCLUMP-1)= ILIVES
  KCELL(ICLUMP-1,JCLUMP-1)= ILIVES
ENDDO

CALL FLUSH(NULOUT)
PWGHT=0.0  ! array initialization

IF (LHOOK) CALL DR_HOOK('STOPH_MIX:INITIALIZE_CELLS',1,ZHOOK_HANDLE)
END SUBROUTINE INITIALIZE_CELLS

!     ------------------------------------------------------
SUBROUTINE WEIGHTING_FIELD(KCELL,PWGHT,KIIP,KJJP)
!     ------------------------------------------------------

!     square the KCELL values (on range 0->31); average over
!     4x4 KCELL blocks, and then normalize

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK 
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KIIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KJJP
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCELL(KIIP*4,KJJP*4)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWGHT(KIIP,KJJP)
INTEGER(KIND=JPIM), PARAMETER    :: ILIVES=10
INTEGER(KIND=JPIM), PARAMETER    :: ISMOOTH=1
INTEGER(KIND=JPIM)               :: I,J,II,JJ
REAL(KIND=JPRB)                  :: ZSUM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('STOPH_MIX:WEIGHTING_FIELD',0,ZHOOK_HANDLE)

ZSUM= 0.0
DO II=0,KIIP-1
   DO JJ=0,KJJP-1
      PWGHT(II+1,JJ+1)= 0.0
      DO I=1,4
        DO J=1,4
          PWGHT(II+1,JJ+1)= PWGHT(II+1,JJ+1) + KCELL(II*4+I,JJ*4+J)
        ENDDO
      ENDDO
       ZSUM= ZSUM + PWGHT(II+1,JJ+1)
   ENDDO
ENDDO

DO I=1,ISMOOTH
  CALL SMOOTH121(PWGHT,KIIP,KJJP)
ENDDO

DO II=1,KIIP
  DO JJ=1,KJJP
    PWGHT(II,JJ)=(KIIP*KJJP)*PWGHT(II,JJ)/ZSUM - 1._JPRB
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('STOPH_MIX:WEIGHTING_FIELD',1,ZHOOK_HANDLE)
END SUBROUTINE WEIGHTING_FIELD

SUBROUTINE SMOOTH121(PFLD,KIIP,KJJP)

!    apply 1-2-1 smoothing in both directions, first in i followed by j then
!    j followed by i... average the resulting values for symmetry

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK 
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KIIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KJJP
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFLD(KIIP,KJJP)
INTEGER(KIND=JPIM)               :: I,J,IM1,IP1,JM1,JP1
REAL(KIND=JPRB), DIMENSION(KIIP,KJJP) :: ZFLD_SMOOTH
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('STOPH_MIX:SMOOTH121',0,ZHOOK_HANDLE)

DO I=1,KIIP
  IM1= I - 1
  IP1= I + 1
  IF (I==1)   IM1= KIIP
  IF (I==KIIP) IP1= 1
  DO J=1,KJJP
    JM1= J - 1
    JP1= J + 1
    IF (J==1)   JM1= KJJP
    IF (J==KJJP) JP1= 1
    ZFLD_SMOOTH(I,J)= 0.0625*(&
    &      PFLD(IM1,JP1) + PFLD(IP1,JP1) + PFLD(IM1,JM1) + PFLD(IP1,JM1) +&
    &      2*( PFLD(I,JP1)   + PFLD(IM1,J)   + PFLD(IP1,J)   + PFLD(I,JM1) ) +&
    &  4*  PFLD(I,J) )
  ENDDO
ENDDO

PFLD= ZFLD_SMOOTH
IF (LHOOK) CALL DR_HOOK('STOPH_MIX:SMOOTH121',1,ZHOOK_HANDLE)
END SUBROUTINE SMOOTH121

SUBROUTINE UPDATE_CELLS(KCELL,KCELLS,KIIP,KJJP)
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK 
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KIIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KJJP
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCELL(KIIP*4,KJJP*4)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCELLS(0:KIIP*4+1,0:KJJP*4+1)
INTEGER(KIND=JPIM)               :: IC,JC
INTEGER(KIND=JPIM), PARAMETER    :: ILIVES=10
INTEGER(KIND=JPIM)               :: I,J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('STOPH_MIX:UPDATE_CELLS',0,ZHOOK_HANDLE)

IC=KIIP*4
JC=KJJP*4
DO I=1,IC
  DO J=1,JC
    KCELLS(I,J)= KCELL(I,J)
  ENDDO
ENDDO

CALL WRAP_CELLS(KCELLS,KIIP,KJJP)

IF (LHOOK) CALL DR_HOOK('STOPH_MIX:UPDATE_CELLS',1,ZHOOK_HANDLE)
END SUBROUTINE UPDATE_CELLS

SUBROUTINE WRAP_CELLS(KCELLS,KIIP,KJJP)
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK 
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KIIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KJJP
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCELLS(0:KIIP*4+1,0:KJJP*4+1)
INTEGER(KIND=JPIM)               :: IC,JC
INTEGER(KIND=JPIM), PARAMETER    :: ILIVES=10
INTEGER(KIND=JPIM)               :: I,J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('STOPH_MIX:WRAP_CELLS',0,ZHOOK_HANDLE)

IC=KIIP*4
JC=KJJP*4
DO I= 0,IC+1
  KCELLS(I,0)= KCELLS(I,JC)
  KCELLS(I,JC+1)= KCELLS(I,1)
ENDDO

DO J= 0,JC+1
  KCELLS(0,J)= KCELLS(IC,J)
  KCELLS(IC+1,J)= KCELLS(1,J)
ENDDO

IF (LHOOK) CALL DR_HOOK('STOPH_MIX:WRAP_CELLS',1,ZHOOK_HANDLE)
END SUBROUTINE WRAP_CELLS

FUNCTION NEIGHBOUR_COUNT(KCOL,KROW,KCELLS,KIIP,KJJP)
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK 
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KIIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KJJP
INTEGER(KIND=JPIM),INTENT(IN)    :: KCOL
INTEGER(KIND=JPIM),INTENT(IN)    :: KROW
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCELLS(0:KIIP*4+1,0:KJJP*4+1)
INTEGER(KIND=JPIM), PARAMETER    :: ILIVES=10
INTEGER(KIND=JPIM)               :: ICOUNT,II,JJ
INTEGER(KIND=JPIM)               :: NEIGHBOUR_COUNT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('STOPH_MIX:NEIGHBOUR_COUNT',0,ZHOOK_HANDLE)

ICOUNT= 0
DO II= -1,1
  DO JJ= -1,1
    IF (.NOT.(II==0.AND.JJ==0)) THEN
      IF (KCELLS(KCOL-II,KROW-JJ) /= 0) ICOUNT=  ICOUNT+1
    ENDIF
  ENDDO
ENDDO

NEIGHBOUR_COUNT= ICOUNT

IF (LHOOK) CALL DR_HOOK('STOPH_MIX:NEIGHBOUR_COUNT',1,ZHOOK_HANDLE)
END FUNCTION NEIGHBOUR_COUNT

FUNCTION FERTILE_NEIGHBOUR_COUNT(KCOL,KROW,KCELLS,KIIP,KJJP)

!  count those neighbours that are youngest i.e. those whose KCELL value is ILIVES

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK 
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KIIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KJJP
INTEGER(KIND=JPIM),INTENT(IN)    :: KCOL
INTEGER(KIND=JPIM),INTENT(IN)    :: KROW
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCELLS(0:KIIP*4+1,0:KJJP*4+1)
INTEGER(KIND=JPIM), PARAMETER    :: ILIVES=10
INTEGER(KIND=JPIM)               :: ICOUNT,II,JJ
INTEGER(KIND=JPIM)               :: FERTILE_NEIGHBOUR_COUNT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('STOPH_MIX:FERTILE_NEIGHBOUR_COUNT',0,ZHOOK_HANDLE)

ICOUNT= 0
DO II= -1,1
  DO JJ= -1,1
     IF (.NOT.(II==0.AND.JJ==0)) THEN
       IF (KCELLS(KCOL+II,KROW+JJ)==ILIVES) ICOUNT= ICOUNT+1
     ENDIF
  ENDDO
ENDDO

FERTILE_NEIGHBOUR_COUNT= ICOUNT

IF (LHOOK) CALL DR_HOOK('STOPH_MIX:FERTILE_NEIGHBOUR_COUNT',1,ZHOOK_HANDLE)
END FUNCTION FERTILE_NEIGHBOUR_COUNT


FUNCTION SIGTOFACT(PSIG)

!**** *SIGTOFACT* - Compute parameter in Laplace distribution from desired standard deviation

!     Purpose.
!     --------
!        Compute the parameter for a Laplace distribution from the 
!        desired standard deviation of a phaseshift
!        Since phaseshifts are periodic the standard deviation of
!        the phaseshifts is different from the standard deviation
!        of the distribution

!**   Interface.
!     ----------
!        =SIGTOFACT(PSIG)

!        Explicit arguments :
!        --------------------

!        PSIG     : Desired standard deviation of phaseshifts

!        Implicit arguments :  none.
!        --------------------

!     Method.
!     -------
!        The standard deviation for a Laplace distribution, considering periodicity
!        on multiples of 2*PI, is given by:
!          SIG=SQRT(2*b*(b-PI/SINH(PI/b)))
!        If the probabilty density function is given by:
!          p=EXP(-ABS(x-mean)/b)/(2*b)
!        There is no analytical expression for b as function of SIG.
!        Hence, a least squares fit of the form b**2=VAR/ALOG(1+VARMAX-VAR)+a1*VAR
!        is used. VARMAX is the maximum variance the phaseshifts can have, i.e. the
!        variance of a uniform distribution between -PI and +PI = PI**2/3
!     Externals.
!     ----------
!      none

!     Reference.
!     ----------
!        none

!     Author.
!     -------
!        Martin Steinheimer  *ECMWF*

!     Modifications.
!     --------------
!        Original : 31-03-2010
!     ------------------------------------------------------------------
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK 
USE YOMCST   , ONLY : RPI
IMPLICIT NONE
REAL(KIND=JPRB),INTENT(IN) :: PSIG
REAL(KIND=JPRB) :: SIGTOFACT
REAL(KIND=JPRB) :: ZA1, ZVAR, ZOFF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('STOPH_MIX:SIGTOFACT',0,ZHOOK_HANDLE)

ZA1=-0.25_JPRB
ZOFF=1+RPI**2/3._JPRB
ZVAR=PSIG**2
IF (ZVAR >= ZOFF-1.0001_JPRB) ZVAR=ZOFF-1.0001_JPRB  !avoid infinite LOG(ZOFF-ZVAR)

SIGTOFACT=SQRT(ZA1*ZVAR+ZVAR/LOG(ZOFF-ZVAR))

IF (LHOOK) CALL DR_HOOK('STOPH_MIX:SIGTOFACT',1,ZHOOK_HANDLE)
END FUNCTION SIGTOFACT




SUBROUTINE SPNORMEI(YDLAP,YDDIM,PX,PXFORC,PSM,KLEV)
!**** *SPNORMEI* - Compute energy input norms in spectral space

!     Purpose.
!     --------
!        Compute the norm in spectral of the energy input due to
!        vorticity and/or divergence forcing at a given m
!        (based on SPNORMBM)

!**   Interface.
!     ----------
!        *CALL* *SPNORMEI(PX,PXFORC,PSM,KFLEV)

!        Explicit arguments :
!        --------------------
!        PX      : Input array
!        PXFORC  : Input array of forcing
!        PSM     : spectrum at constant m
!        KLEV    : number of levels of computation

!        Implicit arguments :  none.
!        --------------------

!     Method.
!     -------

!     Externals. None
!     ----------
!      Called by SPNORMBS.

!     Reference.
!     ----------
!        none

!     Author.
!     -------
!        Martin Steinheimer  *ECMWF*

!     Modifications.
!     --------------
!        Original : 17-02-2010
!     ------------------------------------------------------------------


USE YOMLAP   , ONLY : TLAP
USE YOMDIM   , ONLY : TDIM
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE
TYPE(TLAP)        ,INTENT(IN)  :: YDLAP
TYPE(TDIM)        ,INTENT(IN)  :: YDDIM
REAL(KIND=JPRB)   ,INTENT(IN)  :: PX(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PXFORC(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PSM(:,:)
INTEGER(KIND=JPIM),INTENT(IN)  :: KLEV

INTEGER(KIND=JPIM) :: INM, ISP, JN, JLEV, JNML
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('STOPH_MIX:SPNORMEI',0,ZHOOK_HANDLE)
ASSOCIATE(NUMP=>YDDIM%NUMP, NSMAX=>YDDIM%NSMAX, &
 & MYMS=>YDLAP%MYMS, NASM0=>YDLAP%NASM0, RLAPIN=>YDLAP%RLAPIN)
!$OMP PARALLEL DO SCHEDULE(STATIC,1)  PRIVATE(JNML,JLEV,INM,ISP,JN)
DO JNML=1,NUMP
  INM=MYMS(JNML)
  DO JLEV=1,KLEV

!*       1.    COMPUTE for a given n (=KNM). !not implemented yet 
!                                             -> more or less useless anyway 
!                                             because only working for single processor run
!              -----------------------------

!     ------------------------------------------------------------------

!*       2.    COMPUTE for a given m (=KNM).
!              -----------------------------

    PSM(JLEV,JNML)=0.0_JPRB
    IF(INM == 0)THEN
      DO JN=0,NSMAX
        ISP=NASM0(0)+JN*2
        PSM(JLEV,JNML)=PSM(JLEV,JNML)-RLAPIN(JN)*(&
                    &  2.0_JPRB*PX(JLEV,ISP)*PXFORC(JLEV,ISP)+PXFORC(JLEV,ISP)**2)
      ENDDO
    ELSE
      DO JN=INM,NSMAX
        ISP=NASM0(INM)+(JN-INM)*2
        PSM(JLEV,JNML)=PSM(JLEV,JNML)-2.0_JPRB*RLAPIN(JN)*&
         & ( 2.0_JPRB*PX(JLEV,ISP)*PXFORC(JLEV,ISP)+2.0_JPRB*PX(JLEV,ISP+1)*PXFORC(JLEV,ISP+1)&
         &  +PXFORC(JLEV,ISP)**2+PXFORC(JLEV,ISP+1)**2)  
      ENDDO
    ENDIF

  ENDDO
ENDDO
!$OMP END PARALLEL DO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('STOPH_MIX:SPNORMEI',1,ZHOOK_HANDLE)
END SUBROUTINE SPNORMEI

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TSTOPH), INTENT(IN) :: SELF
INTEGER      , INTENT(IN) :: KDEPTH
INTEGER      , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC

IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_stoch%yrstoph : '
IF (ALLOCATED(SELF%RSTOPHCA)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSTOPHCA ALLOCATED OF SHAPE ',& 
 &        SHAPE(SELF%RSTOPHCA), ' SUM = ', SUM(SELF%RSTOPHCA)
IF (ALLOCATED(SELF%SQRTCORR)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SQRTCORR ALLOCATED OF SHAPE ',& 
 &        SHAPE(SELF%SQRTCORR), ' SUM = ', SUM(SELF%SQRTCORR)
IF (ALLOCATED(SELF%SPSTREAM)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SPSTREAM ALLOCATED OF SHAPE ',& 
 &        SHAPE(SELF%SPSTREAM), ' SUM = ', SUM(SELF%SPSTREAM)
IF (ALLOCATED(SELF%SPVELPOT)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SPVELPOT ALLOCATED OF SHAPE ',& 
 &        SHAPE(SELF%SPVELPOT), ' SUM = ', SUM(SELF%SPVELPOT)
IF (ALLOCATED(SELF%SPSTREAM_FORC)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SPSTREAM_FORC ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%SPSTREAM_FORC), ' SUM = ', SUM(SELF%SPSTREAM_FORC)
IF (ALLOCATED(SELF%SPVELPOT_FORC)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SPVELPOT_FORC ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%SPVELPOT_FORC), ' SUM = ', SUM(SELF%SPVELPOT_FORC)
IF (ALLOCATED(SELF%SPG_AMP)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SPG_AMP ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%SPG_AMP), ' SUM = ', SUM(SELF%SPG_AMP)
IF (ALLOCATED(SELF%ALPHA_STO)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'ALPHA_STO ALLOCATED OF SHAPE ', SHAPE(SELF%ALPHA_STO), ' SUM = ', SUM(SELF%ALPHA_STO)
IF (ALLOCATED(SELF%ONEMINALPHA_NFRSPBS)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'ONEMINALPHA_NFRSPBS ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%ONEMINALPHA_NFRSPBS), ' SUM = ', SUM(SELF%ONEMINALPHA_NFRSPBS)
IF (ALLOCATED(SELF%RSMOOTH)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSMOOTH ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%RSMOOTH), ' SUM = ', SUM(SELF%RSMOOTH)
IF (ALLOCATED(SELF%GPSTREAM)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GPSTREAM ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%GPSTREAM), ' SUM = ', SUM(SELF%GPSTREAM)
IF (ALLOCATED(SELF%GPVELPOT)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GPVELPOT ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%GPVELPOT), ' SUM = ', SUM(SELF%GPVELPOT)
IF (ALLOCATED(SELF%GPTOTDISS)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GPTOTDISS ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%GPTOTDISS), ' SUM = ', SUM(SELF%GPTOTDISS)
IF (ALLOCATED(SELF%GPTOTDISS_SMOOTH)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GPTOTDISS_SMOOTH ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%GPTOTDISS_SMOOTH), ' SUM = ', SUM(SELF%GPTOTDISS_SMOOTH)
IF (ALLOCATED(SELF%GPVORTGRAD)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GPVORTGRAD ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%GPVORTGRAD), ' SUM = ', SUM(SELF%GPVORTGRAD)
IF (ALLOCATED(SELF%MCELL)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'MCELL ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%MCELL), ' SUM = ', SUM(SELF%MCELL)
IF (ALLOCATED(SELF%RWGHT)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RWGHT ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%RWGHT), ' SUM = ', SUM(SELF%RWGHT)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFRSTOPH_SPBS = ', SELF%NFRSTOPH_SPBS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFRSTOPH_SPBS_PAT = ', SELF%NFRSTOPH_SPBS_PAT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFRSTOPH_VC = ', SELF%NFRSTOPH_VC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFORCESTART = ', SELF%NFORCESTART
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NFORCEEND = ', SELF%NFORCEEND
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSTOCHOPT = ', SELF%NSTOCHOPT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'ALPHA_DEEP_CONV = ', SELF%ALPHA_DEEP_CONV
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'ALPHA_SHAL_CONV = ', SELF%ALPHA_SHAL_CONV
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SLDISSFAC = ', SELF%SLDISSFAC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RATIO_BACKSCAT = ', SELF%RATIO_BACKSCAT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSPBS_TAU = ', SELF%RSPBS_TAU
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RATIO_BACKSCAT_CON2NUM = ', SELF%RATIO_BACKSCAT_CON2NUM
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REXPONENT = ', SELF%REXPONENT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'VC_CON = ', SELF%VC_CON
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RFLUX_DET_CLIP = ', SELF%RFLUX_DET_CLIP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'BIHARM = ', SELF%BIHARM
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSMOOTHSCALE = ', SELF%RSMOOTHSCALE
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSIGMA2_EPS = ', SELF%RSIGMA2_EPS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'INCUT = ', SELF%INCUT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'AMAGSTOPH_CASBS = ', SELF%AMAGSTOPH_CASBS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'ADLATSTOPH_CA = ', SELF%ADLATSTOPH_CA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'ADLONSTOPH_CA = ', SELF%ADLONSTOPH_CA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSTOPH_SPBS = ', SELF%LSTOPH_SPBS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSTOPH_SPBS_FAST = ', SELF%LSTOPH_SPBS_FAST
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LEXTRAFIELDS = ', SELF%LEXTRAFIELDS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSTOPH_JBCOR = ', SELF%LSTOPH_JBCOR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSTOPH_UNCORR = ', SELF%LSTOPH_UNCORR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSTOPH_UNIFORM = ', SELF%LSTOPH_UNIFORM
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSTOPH_SPBS_VORT = ', SELF%LSTOPH_SPBS_VORT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSTOPH_RVP = ', SELF%LSTOPH_RVP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSTOPH_RVPOLD = ', SELF%LSTOPH_RVPOLD
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSTOPH_TAPER = ', SELF%LSTOPH_TAPER
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSTOPH_INI = ', SELF%LSTOPH_INI
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSTOPH_CASBS = ', SELF%LSTOPH_CASBS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LVORTCON = ', SELF%LVORTCON
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LFORCENL = ', SELF%LFORCENL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSTOPH_VARALPHA = ', SELF%LSTOPH_VARALPHA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSTOPH_GAUSS = ', SELF%LSTOPH_GAUSS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSPBS_DISSGW = ', SELF%LSPBS_DISSGW
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSPBS_DISSNUM = ', SELF%LSPBS_DISSNUM
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSPBS_DISSCU = ', SELF%LSPBS_DISSCU
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSPBS_DISSNUM_CT = ', SELF%LSPBS_DISSNUM_CT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSPBSNORM = ', SELF%LSPBSNORM
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSPBSDISS = ', SELF%LSPBSDISS
IF (ALLOCATED(SELF%SPDP)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SPDP ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%SPDP), ' SUM = ', SUM(SELF%SPDP)
IF (ALLOCATED(SELF%NIMRAN)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NIMRAN ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%NIMRAN), ' SUM = ', SUM(SELF%NIMRAN)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSMAXSPBS = ', SELF%NSMAXSPBS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MULMIN = ', SELF%RVP_MULMIN
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MULMAX = ', SELF%RVP_MULMAX
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MULEXP = ', SELF%RVP_MULEXP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MULNSMAX = ', SELF%RVP_MULNSMAX
IF (ALLOCATED(SELF%RVP_MULFACT)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MULFACT ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%RVP_MULFACT), ' SUM = ', SUM(SELF%RVP_MULFACT)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MUL_A = ', SELF%RVP_MUL_A
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MUL_B = ', SELF%RVP_MUL_B
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MUL_C = ', SELF%RVP_MUL_C
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MUL_D = ', SELF%RVP_MUL_D
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MUL_1 = ', SELF%RVP_MUL_1
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MUL_2 = ', SELF%RVP_MUL_2
IF (ALLOCATED(SELF%RVP_MUL)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MUL ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%RVP_MUL), ' SUM = ', SUM(SELF%RVP_MUL)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TAPER_SIGMATOP = ', SELF%TAPER_SIGMATOP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TAPER_SIGMABOT = ', SELF%TAPER_SIGMABOT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TAPER0 = ', SELF%TAPER0
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TAPER1 = ', SELF%TAPER1
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TAPER2 = ', SELF%TAPER2
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TAPER3 = ', SELF%TAPER3
IF (ALLOCATED(SELF%TAPER_FACT)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TAPER_FACT ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%TAPER_FACT), ' SUM = ', SUM(SELF%TAPER_FACT)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSTOPH_SPBS_T = ', SELF%LSTOPH_SPBS_T
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'REXPONENT_T = ', SELF%REXPONENT_T
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RATIO_APE2KE = ', SELF%RATIO_APE2KE
IF (ALLOCATED(SELF%SPTEMP)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SPTEMP ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%SPTEMP), ' SUM = ', SUM(SELF%SPTEMP)
IF (ALLOCATED(SELF%SPTEMP_FORC)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SPTEMP_FORC ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%SPTEMP_FORC), ' SUM = ', SUM(SELF%SPTEMP_FORC)
IF (ALLOCATED(SELF%SPG_AMP_T)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SPG_AMP_T ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%SPG_AMP_T), ' SUM = ', SUM(SELF%SPG_AMP_T)
IF (ALLOCATED(SELF%ALPHA_STO_T)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'ALPHA_STO_T ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%ALPHA_STO_T), ' SUM = ', SUM(SELF%ALPHA_STO_T)
IF (ALLOCATED(SELF%ONEMINALPHA_NFRSPBS_T)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'ONEMINALPHA_NFRSPBS_T ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%ONEMINALPHA_NFRSPBS_T), ' SUM = ', SUM(SELF%ONEMINALPHA_NFRSPBS_T)
IF (ALLOCATED(SELF%GPTEMP)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'GPTEMP ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%GPTEMP), ' SUM = ', SUM(SELF%GPTEMP)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MULMIN_T = ', SELF%RVP_MULMIN_T
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MULMAX_T = ', SELF%RVP_MULMAX_T
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MULEXP_T = ', SELF%RVP_MULEXP_T
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MULNSMAX_T = ', SELF%RVP_MULNSMAX_T
IF (ALLOCATED(SELF%RVP_MULFACT_T)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MULFACT_T ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%RVP_MULFACT_T), ' SUM = ', SUM(SELF%RVP_MULFACT_T)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MUL_A_T = ', SELF%RVP_MUL_A_T
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MUL_B_T = ', SELF%RVP_MUL_B_T
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MUL_C_T = ', SELF%RVP_MUL_C_T
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MUL_D_T = ', SELF%RVP_MUL_D_T
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MUL_1_T = ', SELF%RVP_MUL_1_T
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MUL_2_T = ', SELF%RVP_MUL_2_T
IF (ALLOCATED(SELF%RVP_MUL_T)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVP_MUL_T ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%RVP_MUL_T), ' SUM = ', SUM(SELF%RVP_MUL_T)

END SUBROUTINE PRINT_CONFIGURATION

END MODULE STOPH_MIX
