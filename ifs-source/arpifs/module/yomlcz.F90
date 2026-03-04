! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE YOMLCZ

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE CONTROL_VECTORS_MOD, ONLY: CONTROL_VECTOR

IMPLICIT NONE
SAVE

!     -----------------------------------------------------------------
!*    ** *YOMLCZ* - CONTROL PARAMETERS FOR LANCZOS ALGORITHM
!     -----------------------------------------------------------------

#ifdef RS6K
! xlf90 bug
TYPE(CONTROL_VECTOR), SAVE :: YVAZX0, YVAZG0, YSPFORCE
#else
TYPE(CONTROL_VECTOR) :: YVAZX0, YVAZG0, YSPFORCE
#endif
REAL(KIND=JPRB), ALLOCATABLE :: RLANBUF(:,:), RLANBUF_OBSCOR(:,:)
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: GPFORCEU(:,:,:,:),GPFORCEV(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: GPFORCET(:,:,:,:),GPFORCEQ(:,:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: GPFORCESP(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: RLRAIN(:) 
REAL(KIND=JPRB), ALLOCATABLE :: RITZVALS(:) 
INTEGER(KIND=JPIM) :: MEMBFGS
INTEGER(KIND=JPIM) :: NITERL
INTEGER(KIND=JPIM) :: NITERL_OBSCOR
INTEGER(KIND=JPIM) :: NWEIGL
INTEGER(KIND=JPIM) :: NLANCOUNT
INTEGER(KIND=JPIM) :: NLANTYPE
INTEGER(KIND=JPIM) :: NLANNORM
INTEGER(KIND=JPIM) :: NEIGEVO
INTEGER(KIND=JPIM) :: NLEVMIN
INTEGER(KIND=JPIM) :: NLEVMAX
INTEGER(KIND=JPIM) :: NWTRMIN0
INTEGER(KIND=JPIM) :: NWTRMAX0
INTEGER(KIND=JPIM) :: NWTRMIN1
INTEGER(KIND=JPIM) :: NWTRMAX1
INTEGER(KIND=JPIM) :: NEWNORMT0
INTEGER(KIND=JPIM) :: NINNER
INTEGER(KIND=JPIM) :: NJDSTOP
INTEGER(KIND=JPIM) :: NOPMSTOP
INTEGER(KIND=JPIM) :: NRAINSTART
INTEGER(KIND=JPIM) :: N_DIM_SUBSPACE
INTEGER(KIND=JPIM) :: NRITZREAD
INTEGER(KIND=JPIM) :: NRITZNUMB
INTEGER(KIND=JPIM) :: NRITZGBH1
INTEGER(KIND=JPIM) :: NRITZGBH2
REAL(KIND=JPRB) :: GREDBFGS
REAL(KIND=JPRB) :: XKAPA
REAL(KIND=JPRB) :: XMIN_RITZ
REAL(KIND=JPRB) :: ALAT1
REAL(KIND=JPRB) :: ALON1
REAL(KIND=JPRB) :: ALAT3
REAL(KIND=JPRB) :: ALON3
REAL(KIND=JPRB) :: COEQTERM
REAL(KIND=JPRB) :: COENEWQTERM
REAL(KIND=JPRB) :: TSTEP_STATE_4D
LOGICAL :: LANCZOS
LOGICAL :: LEVOLC
LOGICAL :: LSCALC
LOGICAL :: LOCNORM
LOGICAL :: LSPTRLC0
LOGICAL :: LSPTRLC1
LOGICAL :: LSELU
LOGICAL :: LSELV
LOGICAL :: LSELT
LOGICAL :: LSELQ
LOGICAL :: LSELSP
LOGICAL :: LNEWNORMT0
LOGICAL :: LCHSYMEIG
LOGICAL :: LSYMCHECK
LOGICAL :: LFORCE
LOGICAL :: LFORCEWR
LOGICAL :: LRAIN
LOGICAL :: LRENORMALIZE
LOGICAL :: L_USE_CONGRAD
LOGICAL :: L_SOS
LOGICAL :: L_EOFS
LOGICAL :: L_BALANCED_REDUCTION
LOGICAL :: L_SUBSPACE_SVS
LOGICAL :: LSMSSIG_SVINI

TYPE (CONTROL_VECTOR), ALLOCATABLE :: YSTATE_VECTOR_4D(:)
TYPE (CONTROL_VECTOR), DIMENSION(:), ALLOCATABLE :: YV_SUBSPACE
INTEGER(KIND=JPIM) :: NSTEPS_PER_STATE

!     ------------------------------------------------------------------

!*     Variables for eigensystem evaluation with Lanczos code

!      LANCZOS : .T. to compute the unstable perturbations, false to
!                 read them from an input file and compute the TL evolution
!      LFORCE  : .T. to compute Forcing SVs: unstable perturbations for
!                     model tendencies
!      L_USE_CONGRAD : .T. to use CONGRAD Lanczos code. .F. to use LANDR.

!      L_SOS   : .T. to compute (finite time integral) stochastic optimals
!                    i.e. eigenvectors of Q=\int_0^T M^* M dt (where T is NSTOP)
!      L_EOFS  : .T. to compute (finite time integral) EOFs
!                    i.e. eigenvectors of P=\int_0^T M M^* dt (where T is NSTOP)
!      L_BALANCED_REDUCTION : .T. to compute balancing vectors X and Y for model reduction
!      L_SUBSPACE_SVS : .T. to compute SVs in a subspace orthogonal to a 
!                       given subspace, identified by N_DIM_SUBSPACE vectors 
!                       in file 'sv_subspace'
!      TSTEP_STATE_4D : time interval between states in the 4d state vector.

!      LEVOLC   : .T. to write the evolution of the basic state and
!                 of some eigenvectors (TL evol)
!      LSCALC   : .T. to write in the unit NULUSR3 the scalar product
!                 weight factors SCALP
!      LOCNORM   : .T. to localize the norm computation in gg space
!      LSPTRLC   : .T. to truncate the state vector in the spectral space
!                           0 = at initial time
!                           1 = at final time
!      LNEWNORMT0: .T. to re-define the norm at initial time t0
!      LRENORMALIZE: .T. to normalize SVs with TE metric at initial time
!      LCHSYMEIG:  .T. to check eigenvectors and the symmetry of the Hessian
!      NITERL    : maximum number of Lanczos/Jacobi-Davidson(outer) iterations
!      NINNER    : number of Jacobi-Davidson inner iterations
!      NJDSTOP   : value of NSTOP used for TL/AD integrations
!      NOPMSTOP  : value of NSTOP used in Hessian calculation
!      NRAINSTART: value of NSTEP when accumulation of rain starts
!      NWEIGL    : maximum number of wanted eigenvalues
!      NLANCOUNT : counter for the Lanczos iteration
!      N_DIM_SUBSPACE : number of used vectors in 'sv_subspace'
!      XKAPA     : acceptable precision for the eigenvalues
!      XMIN_RITZ : smallest accepted eigenvalue for CONGRAD
!      COEQTERM  : multiplication coefficient in front of the q-term
!                  in the total energy norm scalar weights
!      NLANTYPE  : type of integration:
!                           1 = for TL and AD integration
!                           2 = for TL and AD with covariance matrix constraint
!                           3 = multiply by covariance matrix only
!                 ...      99 = no integration, OBS ERR correlations          
!      NLANNORM  : norm control: 1 = for Total Energy
!                                2 = for Kinetic Energy
!      ALAT1,ALON1,ALAT3,ALON3 : coordinate of the local area used for
!                                norm computation
!      NLEVMIN, NLEVMAX : min and max level for the local area norm comp
!      NEWNORMT0: control the re-definition of the norm at t0:
!             1 to re-define the norm at initial time to be total-energy
!             2                                             kin-energy
!             3                                             vor-2
!             4                                             psi-2
!             5                                             rot-kin-energy
!      COENEWQTERM:multiplication coefficient in front of the q-term
!                  in the total energy norm scalar weights
!      NWTRMIN, NWTRMAX : min and max total wavenumber of the spectral
!                 components not to be set to zero by the spectral truncation
!                           0 = at initial time
!                           1 = at final time
!      LSELU  : .T. to select the U-wind in the local norm
!      LSELV  : .T. to select the V-wind in the local norm
!      LSELT  : .T. to select the temperature in the local norm
!      LSELQ  : .T. to select the humidity in the local norm
!      LSELSP : .T. to select the surface pressure in the local norm
!      LRAIN  : .T. to select rain in the final norm
!      RLANBUF: memory buffer used by Lanczos alg.
!      MEMBFGS: maximum number of (y,s) pairs used to form the BFGS
!               preconditioner
!      GREDBFGS:required reduction in gradient norm for PCGBFGS
!      YVAZX0  : Control variable at reference point (used in OPM).
!      YVAZG0  : Gradient at reference point (used in OPM).

!      YSTATE_VECTOR_4D : Array of state vectors (used for balanced reduction).
!      LSMSSIG_SVINI : set SMS event i when initial SVs have been computed

END MODULE YOMLCZ
