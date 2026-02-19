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

MODULE YOMVAR

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Switches for variational assimilation
! LTEST    : .T. = TEST OF THE GRADIENT BEFORE THE MINIMIZATION
! LREPRO4DVAR: .T. = Runs 4D-Var in bit reproducible mode (slower)
! LSLADREP   : CONTROL BIT REPRODUCIBILITY
!            : T - SLAD IS BIT REPRODUCIBLE WHEN NUMBER OF PROCESSORS
!            :     OR PARTITITIONING IS CHANGED
!            : F - SLAD IS NOT BIT REPRODUCIBLE WHEN NUMBER OF PROCESSORS
!            :     OR PARTITITIONING IS CHANGED
! L_CHECK_CONVERGENCE : .T. = ABORT IF THE MINIMIZATION FAILS
! R_NORM_REDUCTION_ABORT_LEVEL: = ABORT IF THE RATIO OF FINAL/INITIAL
!                            GRADIENT NORM SQUARED EXCEEDS THIS FACTOR
! L_CHECK_GRADIENT: .T. => CONGRAD checks its gradient each iteration.
!                          (This performs an additional call to SIM4D each iteration.
!                           If a problem is found, the gradient is written out and
!                           the minimization aborts.)
! RTOL_CHECK_GRADIENT: Gradient error norm (relative to gradient norm) required
!                      to trigger an abort if L_CHECK_GRADIENT==.T..
! LFCOBS  : .F. = Switch to compute the forecast sensitivity with respect to
! LFCOBSTEST  : .F. = Switch to compute the spectral and grid point norms during the temporal loop in cnt4tl/ad
! LCLDSINK : .T. = SWITCH FOR CLOUD SINK VARIABLE
! CTOPBGE : 5.= BG ERROR   CLOUD SINK VARIABLE
! CAMTBGE : 0.01 = BG ERROR   CLOUD SINK VARIABLE
! LTOVSCV : .T. = SWITCH FOR TOVS CONTROL VARIABLE
! LTOVSREP: .T. = Reproducible TOVS scalar product
! LJC     : .T. = ADD A JC TERM IN THE COST FUNCTION
! LJCDFI  : .T. = compute a digital formulation of a Jc term
! LUSEJCDFI: .T. = a Jc-DFI term is added to the cost function when
! LINITCV : .T. = include initial condition in 4D-Var control variable.
! LMODERR : .T. = include model error term in 4D-Var.
! LVARBC  : .T. = include bias parameters in control variable.
! LPERTMRPAS : .T. = blacklist passive observations in perturbed EDA members
! LTRREF  : Controls what is read/write in array SPA5 or SPA7
!          .T.: reference trajectory  ---> SPA7
!          .F.: current trajectory    ---> SPA5
! LREFINC : Informs on the nature of the reference SPA7 content (with LFCOBS).
!          .T. SPA7 contains an increment
!              (differences between 2 historic fieldsets)
!          .F. SPA7 contains an historic fieldset.
! LZOWA   : .T. = DIAGNOSTIC OF COST FUNCTION BY ZONAL WAVE NUMBER
!           .F. = DIAGNOSTIC OF COST FUNCTION BY TOTAL WAVE NUMBER
! LTWANA  : .T. = WRITE ON FA FILES THE CURRENT ANALYSIS
! LTWGRA  : .T. = WRITE ON FA FILES THE CURRENT GRADIENT
! LTWLCZ  : .T. = WRITE ON FA FILES THE CURRENT SV
! LTWINC  : .T. = WRITE ON FA FILES THE CURRENT INCREMENT
! LTWBGV  : .T. = WRITE ON FA FILES A RANDOM BACKGROUND VECTOR
! LTWCGL  : .T. = WRITE ON FA FILES AN EIGENVECTOR OF THE HESSIAN
! LGRASCAL: .T. = WRITE GRADIENT WITH RESPECT TO SCALP INNER PRODUCT.
! LSKIPMIN: .T. = FORCE SKIPPING OF MINIMIZATION (For diagnostics)
! L_ABS_CONVERGENCE: .T. => convergence criterion is the absolute
!                           value of the final gradient.
! L_INFO_CONVERGENCE: .T. => convergence criterion is relative increase
!                            in information (DFI) per iteration.
! NITER   : MAX NUMBER OF ITERATIONS
! NSIMU   : MAX NUMBER OF SIMULATIONS
! NITER_MIN : MINIMUM NUMBER OF ITERATIONS (NB: CONGRAD ONLY)
! NINFRA  : EXPECTED FRACTION OF THE DIMINUTION OF THE COST FUNCTION
! RCVGE   : Convergence of the minimization is judged to have been achieved
!           when:
!             IF (L_ABS_CONVERGENCE) THEN
!               The norm of the gradient is less than RCVGE.
!              (NB: The value specified by RCVGE is not the true norm of the
!               proper cost function. It is the square of twice the norm.
!               I.e. it is the value reported as FINAL GRADIENT.)
!             ELSEIF (L_INFO_CONVERGENCE) THEN
!               The increase in Jb compared with the preceding iteration
!               is less than RCVGE. (Note that Jb represents "degrees
!               of freedom for signal", DFI, so the change in Jb per
!               iteration is a measure of the amount of information added
!               per iteration.
!             ELSE
!               The norm of the gradient is reduced by a factor of RCVGE.
!             ENDIF
! NMIMP   : CONTROLS THE IMPRESSIONS OF THE OPTIMIZER
! NSIM4D  : COUNTER (number of calls to SIM4D)
! NITER4D : COUNTER (number of iterations of the minimisation)
! NSIM4DL : VALUE OF NSIM4D AT LAST SIMULATION
! NDIAG   : Diagnostic management (switched on only for M1GCX)
! RDX     : starting point for GRTEST (generally 1.E-15)
! RXMIN   : minimum distance in the sup-norm distinguishable by the optimizer
! ALPHAG  : weighting coefficient for the weak constraint term Jc.
! ALPHAV  : weighting coefficient to scale the Jr cost function & gradient
! NOANEF  : Number of analysis error fields
! NBGVECS : Number of random background vectors to be generated.
! LBGTRUNC: .T. ===> filter the randomization estimate of sigma_b.
! NBGTRUNC: wavenumber above which sigma_b coefficients are zero.
! LBGOBS  : .T. ===> Compute Bg errors for observed quantities.
! LBGM    : .T. ===> Propagate bg-errors using TL model.
! LANOBS  : .T. ===> Compute and Store in ODB the analysis sensitivity to observations
! LWREINI : .T. ===> write perturbed initial file for LELAM and not LTEST
! LWRIBVEC: .T. ===> write random background vectors to files
! LWRIBVEC_FULL: .T. ===> write random background vectors to gribfull files
! LREABVEC: .T. ===> read  random background vectors from files
! LWRIEVEC: .T. ===> write eigenvectors of the Hessian to files
! NWRIEVEC: Maximum number of Hessian eigenvectors to write
! LEVECCNTL : .T. ===> write eigenvectors in control space (before CHAVARIN)
! LEVECGRIB : .T. ===> write eigenvectors in GRIB format
! N_DIAGS_CONVERGENCE: Convergence diagnostics level
!                      0 - No diagnostics
!                      1 - Print diagnostics if bad convergence
!                      2 - Print diagnostics
! N_DIAGS_EIGENVECS  : Save eigenvectors in GRIB (could replace LEVECGRIB)
!                      0 - Don't save eigenvectors
!                      1 - Save eigenvectors if bad convergence
!                      2 - Save eigenvectors
! LWRISIGB: .T. ===> write standard deviation of background error
! LWRISIGA: .T. ===> write standard deviation of analysis   error
! LWRISIGF: .T. ===> write standard deviation of forecast   error
! CFNSIGB : file to which standard deviation of background error is written
! CFNSIGA : file to which standard deviation of analysis   error is written
! CFNSIGF : file to which standard deviation of forecast   error is written
! MBGVEC  : COUNTER (number of random vector being constructed).
! NFGFCLEN: Length of first guess forecast (hours)
! LN1CG1  : .T. ===> Minimize using N1CG1
! LCONGRAD: .T. ===> Minimize using CONGRAD (cost function must be quadratic)
! L3DFGAT : .T. ===> Use 3d FGAT (4dVar with TL model=identity in minimization)
! NPRECO  : 0 ===> N1CG1 : unpreconditioned CG
!         : 2 ===> N1CG1 : L_BFGS preconditioning
! NBFGSB  : 0 ===> N1CG1 : don't build any preconditioner
!         : 2 ===> N1CG1 : build a L_BFGS preconditioner
! N1IMP   : Printing level
! NSELECT : Selection of the pairs to build the L_BFGS preconditioner
!         : 0 ===> FIFO strategy (last pairs saved)
!         : 1 ===> Uniform selection (pairs are distributed uniformly through CG run)
!         : 2 ===> Selection by the Rayleigh quotient
! ZEPSNEG : N1CG1 : control the positivity of the Hessian during the minimization
! LAMV_REASSIGN_PASSIVE : Passive calculation of predictors for AMV bias correction
! LAMV_HEIGHT_ADJUST: Reassign heights of AMV diagnosed above the model cloud
! LREO3_BCOR: Switch for ozone bias correction
! LCH4_BCOR: Switch for methane bias correction

! L_GUESS_RUNTIME     : .T. => estimate when the minimization will end.
! NITER_GUESS_RUNTIME : Estimate when the job will end at iteration NITER_GUESS_RUNTIME...
! NFREQ_GUESS_RUNTIME : ...and every NFREQ_GUESS_RUNTIME iterations thereafter.
! NMEM_GUESS_RUNTIME  : Base the estimate on the NMEM_GUESS_RUNTIME most recent iterations.

! FILTERFACTOR : Filtering factor applied to trajectory fields
! FILTEREXPO: Filtering exponent applied to trajectory fields
! FILTERRESOL: Filtering of trajectory fields only applied for this resolution and above

! LJBIMPACT: Obs impact on Jb+Jq
! LJCIMPACT: Obs impact on Jc
! LMONITOR_FCDEPAR: configuration for forecast departure monitoring
! ======== FORECAST DEPARTURE MONITORING============================
! NUPTRA_RANGE : index of forecat timerange
!========= SIMULATED OBSERVATION EVENTS MANAGEMENT==================

! NFRREF  : frequency of reference observation events
! NREFTS  : array containing observation events steps
!     EXPLANATION :
!     1) IF NREFTS(0)=0 ACTION IF MOD(JSTEP,NFRREF)=0
!     2) IF NREFTS(0)>0 NREFTS(0) SIGNIFICANT NUMBERS IN
! NREFTS ARE THEN CONSIDERED AND :
!       ACTION FOR JSTEP=NREFTS(.)*NFRREF

!========= WRITES CURRENT ANALYSIS ON FILES =================

! NFRANA  : frequency of writes, relative to the number of simulations
! NANATS  : array containing write events steps
!     EXPLANATION :
!     1) IF NANATS(0)=0 ACTION IF MOD(JSTEP,NFRANA)=0
!     2) IF NANATS(0)>0 NANATS(0) SIGNIFICANT NUMBERS IN
! NANATS ARE THEN CONSIDERED AND :
!       ACTION FOR NSIM4D=NANATS(.)*NFRANA
! MSIME   : number of the current emsemble member

!========= WRITES CURRENT GRADIENT ON FILES =================

! NFRGRA  : frequency of writes, relative to the number of simulations
! NGRATS  : array containing write events steps
!     EXPLANATION : As for NFRANA and NANATS above

!========= ASSIMILATION WITH THE T.L. MODEL (nconf 131) ==========
! NUPTRA  : number of updates of the trajectory during the minimisation.
! MUPTRA  : Maximum number of updates of the trajectory during the minimisation

!========= Combined conjugate-gradient and Lanczos algorithm =====
! LAVCGL  : .T. ====> use combined conjugate-gradient / Lanczos scheme
! LMPCGL  : .T. ====> precondition conjugate-gradient minimization
! R_MAX_CNUM_PC : Maximum allowed condition number for the preconditioner
! CFNPCV  : filename prefix for files containing preconditioner vectors
! NPCVECS : number of vectors which make up the preconditioner
! EVBCGL  : max relative error in eigenvalues of written-out vectors
! MCGLVEC : COUNTER identifies the vector being written
! LGCV    : .T. => Calculate the Generalized Cross Validation function
! NITERGCV: Number of iterations for trace calculation
! LGCVJO  : (Internal) calculate Jo gradient for random departues.
! GCVJO   : (Internal) Jo at analysis point
! LFDBERR :  .T. => write bg and an error estimates to FDB

! L_INFO_CONTENT      : .T. => Calculate information content (degrees of freedom for signal, etc.)
! N_INFO_CONTENT_METHOD : 1 => Use Bai et al's algorthm
!                         2 => Use analysis-difference method
! N_INFO_CONTENT_SEED : Random number seed for randomized trace estimate.
! LPROPTL: .T. => 3DFGAT through 4DVAR looping with backward TL propagation
! LAEOLUSAMD  : .F. = include Aeolus production of auxiliary met data
! LBACKGE    : .T. = compute background error variance from an ensemble for full variables
! LBACKGECV  : .T. = compute background error variance from an ensemble for control variables
! LBACKGERENORM : .F. = compute renormalisation coefficients induced by C matrix in wavelet space
! LCONSTANTZFCE : .F. = specify that the variance is constant (and equal to 1) when doing a randomisation of B matrix
! LUSEWAVRENORM : .F. = use renormalisation coefficients to renormalise the correlation matrix
! LWRISB_VPROF:.T. = write vertical profile of sigmab from writesd/fltbgerr (in bgvecs and bgevecs)
! LFAMEMBERS : .T. = read ensemble members in FA format (no need for femars/readvec anymore)
! LTOY42     : Small resolution 4DVAR
! LSUSPQLIM  : Disable SUSPQLIM calculation
! LSPINT     : .T. = the low resolution fg and analyses will be directly read into
!                    the high resolution run and padded with zeros (ECMWF only)
! NDATE_TIME_WINDOW_END  integer date and time of 4dvar window end
! LBGPERT  : .T. = RANDBG perturbations
! DELTA : scaling factor of the size of random perturbations
! LCHRESINCR : .F. : Read control vector, change resolution in control vector
!                  space, transform to analysis incremnt, add to background
!                  and write out (HIRLAM)
! LINC_TOVSCV : update skin temperature estimate starting from estimate in previous outer loop
! LUSE_EDA_SKT: use skin temperature background error calculated from the EDA
! LECV : extended control variable
! L_FGNOTBG : .T. = First minimization starts with first guess not equal to background
! L_BGREC : .T. = Re-centre background (in EDA, re-centre perturbed background on control background)
! L_FGREC : .T. = Re-centre first-guess (in EDA, add total increment from control minimisation NFGREC_MIN)
!                 Re-centre by adding increments from control, cannot use control first-guess directly
!                 because background of vontrol and EDa members differ.
! NFGREC_MIN : Number of control minimisation to take first-guess from. We take the sum
!              of the control analysis increments 0 to NFGREC_MIN.
! L_COMPBGDEP : .T. = Background departures calculated from a separate background trajectory
! LWRICV    : .T. = Write out control vector to file
INTEGER(KIND=JPIM), PARAMETER :: JPNRFT=40
INTEGER(KIND=JPIM),PROTECTED :: NREFTS(0:JPNRFT)
INTEGER(KIND=JPIM),PROTECTED :: NANATS(0:JPNRFT)
INTEGER(KIND=JPIM),PROTECTED :: NGRATS(0:JPNRFT)
CHARACTER (LEN = 80),PROTECTED ::  CFNSIGB
CHARACTER (LEN = 80),PROTECTED ::  CFNSIGA
CHARACTER (LEN = 80),PROTECTED ::  CFNSIGF
CHARACTER (LEN = 80),PROTECTED ::  CFNPCV
INTEGER(KIND=JPIM),PROTECTED :: NITER
INTEGER(KIND=JPIM),PROTECTED :: NITER_MIN
INTEGER(KIND=JPIM),PROTECTED :: NSIMU
INTEGER(KIND=JPIM) :: NINFRA
INTEGER(KIND=JPIM),PROTECTED :: NMIMP
!INTEGER(KIND=JPIM),PROTECTED :: NSIM4D
INTEGER(KIND=JPIM) :: NDIAG
INTEGER(KIND=JPIM),PROTECTED :: NFRREF
INTEGER(KIND=JPIM),PROTECTED :: NFRANA
INTEGER(KIND=JPIM),PROTECTED :: NFRGRA
INTEGER(KIND=JPIM),PROTECTED :: NSIM4DL
INTEGER(KIND=JPIM) :: MSIME
!INTEGER(KIND=JPIM),PROTECTED :: NUPTRA
INTEGER(KIND=JPIM) :: NOBSCURE_NUMBER ! The name explains all
INTEGER(KIND=JPIM),PROTECTED :: MUPTRA
INTEGER(KIND=JPIM) :: NPCVECS ! Changed in PREPPCM
INTEGER(KIND=JPIM) :: NBGVECS ! Changed in BGEVECS
INTEGER(KIND=JPIM) :: NHEVECS ! Changed in XFORMEV
INTEGER(KIND=JPIM) :: NSSBGV_COUNT ! Counter
INTEGER(KIND=JPIM) :: NSSHEV_COUNT ! Counter
INTEGER(KIND=JPIM) :: MCGLVEC !Counter
INTEGER(KIND=JPIM) :: MBGVEC ! Counter
INTEGER(KIND=JPIM) :: NOANEF ! Set in SUANEBUF
INTEGER(KIND=JPIM),PROTECTED :: NFGFCLEN
!INTEGER(KIND=JPIM),PROTECTED :: NITER4D
INTEGER(KIND=JPIM),PROTECTED :: NBGTRUNC
INTEGER(KIND=JPIM),PROTECTED :: NITERGCV
INTEGER(KIND=JPIM),PROTECTED :: NWRIEVEC
INTEGER(KIND=JPIM),PROTECTED :: N_DIAGS_CONVERGENCE
INTEGER(KIND=JPIM),PROTECTED :: N_DIAGS_EIGENVECS
INTEGER(KIND=JPIM),PROTECTED :: NDATE_TIME_WINDOW_END
REAL(KIND=JPRB),PROTECTED :: RDX
REAL(KIND=JPRB),PROTECTED :: ALPHAG
REAL(KIND=JPRB),PROTECTED :: ALPHAV
REAL(KIND=JPRB),PROTECTED :: RXMIN
REAL(KIND=JPRB),PROTECTED :: EVBCGL
REAL(KIND=JPRB),PROTECTED :: RCVGE
REAL(KIND=JPRB),PROTECTED :: R_NORM_REDUCTION_ABORT_LEVEL
REAL(KIND=JPRB) :: GCVJO
REAL(KIND=JPRB),PROTECTED :: ZEPSNEG
REAL(KIND=JPRB),PROTECTED :: R_MAX_CNUM_PC
REAL(KIND=JPRB),PROTECTED :: RTOL_CHECK_GRADIENT
REAL(KIND=JPRB),PROTECTED :: FILTERFACTOR
REAL(KIND=JPRB),PROTECTED :: FILTEREXPO
REAL(KIND=JPRB),PROTECTED :: FILTERRESOL
REAL(KIND=JPRB),PROTECTED :: CTOPBGE
REAL(KIND=JPRB),PROTECTED :: CAMTBGE
REAL(KIND=JPRB),PROTECTED :: RCOEFCO2         ! Cooef. std. CO2 MACC
REAL(KIND=JPRB),PROTECTED :: RCOEFCH4         ! Cooef. std. CH4 MACC
REAL(KIND=JPRB),PROTECTED :: RCOEFCO          ! Cooef. std. CO CHEM
REAL(KIND=JPRB),PROTECTED :: RCOEFNO2         ! Cooef. std. NO2 CHEM
REAL(KIND=JPRB),PROTECTED :: RCOEFGO3         ! Cooef. std. O3 CHEM
REAL(KIND=JPRB),PROTECTED :: DELTA
LOGICAL :: LFCOBS
LOGICAL,PROTECTED :: LFCOBSTEST
LOGICAL,PROTECTED :: LTEST
LOGICAL,PROTECTED :: LREPRO4DVAR
LOGICAL,PROTECTED :: LSLADREP
LOGICAL :: LTRREF
LOGICAL,PROTECTED :: LREFINC
LOGICAL,PROTECTED :: LZOWA
LOGICAL :: LJC
LOGICAL,PROTECTED :: LINITCV
LOGICAL :: LENSCV
LOGICAL,PROTECTED :: LMODERR
LOGICAL,PROTECTED :: LVARBC
LOGICAL,PROTECTED :: LPERTMRPAS
LOGICAL :: LTWANA
LOGICAL :: LTWGRA
LOGICAL :: LTWLCZ
LOGICAL,PROTECTED :: LAVCGL
LOGICAL :: LMPCGL
LOGICAL :: LTWINC
LOGICAL,PROTECTED :: LGRASCAL
LOGICAL,PROTECTED :: LSKIPMIN
LOGICAL,PROTECTED :: L_ABS_CONVERGENCE
LOGICAL,PROTECTED :: L_INFO_CONVERGENCE
LOGICAL :: LTWBGV ! Set in BGVECS
LOGICAL :: LTWCGL ! Set in XFORMEV
LOGICAL,PROTECTED :: LWRIBVEC
LOGICAL,PROTECTED :: LWRIBVEC_FULL
LOGICAL,PROTECTED :: LREABVEC
LOGICAL,PROTECTED :: LWRIEVEC
LOGICAL,PROTECTED :: LEVECCNTL
LOGICAL,PROTECTED :: LEVECGRIB
LOGICAL,PROTECTED :: LWRISIGB
LOGICAL,PROTECTED :: LWRISIGA
LOGICAL,PROTECTED :: LWRISIGF
LOGICAL,PROTECTED :: LBGTRUNC
LOGICAL,PROTECTED :: LBGOBS
LOGICAL,PROTECTED :: LBGM
LOGICAL,PROTECTED :: LANOBS
LOGICAL :: L_CHECK_CONVERGENCE
LOGICAL,PROTECTED :: L_CHECK_GRADIENT
LOGICAL :: LJCDFI ! Turned off and on like a jojo
LOGICAL,PROTECTED :: LUSEJCDFI
LOGICAL,PROTECTED :: LTOVSCV
LOGICAL,PROTECTED :: LTOVSREP
LOGICAL,PROTECTED :: LGCV
LOGICAL,PROTECTED :: L_FGNOTBG
LOGICAL,PROTECTED :: L_BGREC
LOGICAL,PROTECTED :: L_FGREC
INTEGER(KIND=JPIM),PROTECTED :: NFGREC_MIN
LOGICAL,PROTECTED :: L_COMPBGDEP
LOGICAL,PROTECTED :: LWRICV
LOGICAL :: LGCVJO ! FORECAST_ERROR
LOGICAL,PROTECTED :: LWREINI
LOGICAL,PROTECTED :: LN1CG1
LOGICAL,PROTECTED :: LCONGRAD
LOGICAL :: L3DFGAT ! In TLPROP
LOGICAL,PROTECTED :: LFDBERR
LOGICAL,PROTECTED :: LAMV_REASSIGN_PASSIVE
LOGICAL,PROTECTED :: LAMV_HEIGHT_ADJUST
LOGICAL,PROTECTED :: LREO3_BCOR
LOGICAL,PROTECTED :: LCH4_BCOR
LOGICAL,PROTECTED :: L_INFO_CONTENT
LOGICAL,PROTECTED :: LCLDSINK
LOGICAL,PROTECTED :: L_GUESS_RUNTIME
LOGICAL,PROTECTED :: LPROPTL
LOGICAL,PROTECTED :: LAEOLUSAMD
LOGICAL :: LBACKGE
LOGICAL,PROTECTED :: LBACKGECV
LOGICAL,PROTECTED :: LBACKGERENORM
LOGICAL,PROTECTED :: LUSEWAVRENORM
LOGICAL,PROTECTED :: LCONSTANTZFCE
LOGICAL,PROTECTED :: LWRISB_VPROF
LOGICAL :: LDIAG_LCT ! SET in BGVECS and BGEVECS
LOGICAL,PROTECTED :: LFAMEMBERS
LOGICAL,PROTECTED :: LTOY42
LOGICAL,PROTECTED :: LSUSPQLIM
LOGICAL,PROTECTED :: LSPINT
LOGICAL,PROTECTED :: LJBIMPACT = .FALSE.
LOGICAL,PROTECTED :: LJCIMPACT = .FALSE.
LOGICAL :: LMONITOR_FCDEPAR = .FALSE. ! Changed in shuffle used by Odbtools in ODB - uncomprehensible
LOGICAL,PROTECTED :: LENDA = .FALSE.
LOGICAL,PROTECTED :: LBGPERT = .false.
LOGICAL,PROTECTED :: LCHRESINCR
LOGICAL,PROTECTED :: LINC_TOVSCV
LOGICAL,PROTECTED :: LUSE_EDA_SKT= .false.
LOGICAL,PROTECTED :: LECV
INTEGER(KIND=JPIM),PROTECTED :: NUPTRA_RANGE = 0
INTEGER(KIND=JPIM),PROTECTED :: N_INFO_CONTENT_METHOD
INTEGER(KIND=JPIM),PROTECTED :: N_INFO_CONTENT_SEED
INTEGER(KIND=JPIM),PROTECTED :: N1IMP
INTEGER(KIND=JPIM),PROTECTED :: NPRECO
INTEGER(KIND=JPIM),PROTECTED :: NBFGSB
INTEGER(KIND=JPIM),PROTECTED :: NSELECT
INTEGER(KIND=JPIM),PROTECTED :: NITER_GUESS_RUNTIME
INTEGER(KIND=JPIM),PROTECTED :: NFREQ_GUESS_RUNTIME
INTEGER(KIND=JPIM),PROTECTED :: NMEM_GUESS_RUNTIME
!     ------------------------------------------------------------------
CONTAINS
SUBROUTINE SETUP_VAR(KULOUT)

!**** *SETUP_VAR*   - Routine to initialize variational flags common

!     Purpose.
!     --------
!        Initialize variational control modules: YOMSENS, YOMVRTL, YOMCOSJB, YOMVAR,
!         and also some variables of YOMJQ, YOMMODERR, YOMVCGL, YOMDIAGVAR
!        Reads namelist NAMSENS, NAMVRTL and NAMVAR.

!**   Interface.
!     ----------
!        *CALL* *SETUP_VAR(...)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!         output modules: YOMSENS, YOMVRTL, YOMCOSJB, YOMVAR
!                         a subset of YOMJQ, YOMMODERR,  YOMVCGL, YOMDIAGVAR

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Authors.
!     -------
!      Philippe Courtier  and Jean-Noel Thepaut *DMN/ECMWF*  90-12-01

!     Modifications.
!     --------------
!      Modified by N. Bormann   : 01-09-21 Initialise lamv_reassign_passive
!      Modified by A. Dethof    : 03-04240 Initialise lreo3_bcor
!      Modified by D. Dee       : 03-10-17 Initialise lvarbc
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Y.Tremolet    01-Apr-2004 Setup model error stats
!      G. Radnoti    03-12-2004 : initialize lproptl
!      D. Tan        14-Mar-2005: Initialise laeolusamd, laeolusl2bp
!      A. Benedetti  14-Mar-2005: Initialise laerod
!      G. Desroziers and K. Yessad (sept 2005):
!       - split option LTRAJHR into LTRAJHR_ALTI and LTRAJHR_SURF.
!       - adapt option LTRAJHR to METEO-FRANCE configurations.
!      Modified by R. Engelen    : Mar 2008 Initialise lch4_bcor
!      G. Desroziers 28-Feb-2008: Possible computation of ensemble bg error variance
!      Y.Tremolet    27-Nov-2008 Setup for long windows
!      G. Desroziers 22-Dec-2008: Enable transf. of ARPEGE file in GRIB format (to be used in femars)
!      H.Varella     15-Nov-2011: Option LWRIBVEC_FULL included
!      M. Rennie     13-Apr-2012: Remove LAEOLUSL2BP
!      K. Yessad (oct 2013): cleaning, re-write in a readable way.
!      LF. Meunier   29-Oct-2013: Remove L_OPENMP_CV
!      K. Yessad (July 2014): use L4DVAR, in order to avoid use of NSTOP.
!      Y. Michel     10-Mar-2015 Local Correlation Tensor
!      V.Chabot      27-Janv-2016 Renormalisation coefficient for wavelet correlation matrix
!      N. Bormann    01-Aug-2016 CVarBC
!      M. Hamrud     May 2018 module routine setup_var, with contents previously in suvar
!      K. Lean       10-Jan-2020 Option LAMV_HEIGHT_ADJUST added
!     -------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULNAM
USE YOMCT0   , ONLY : NCONF, LOBSC1, LBACKG, LFDBOP, LECMWF, LELAM, L4DVAR, L_OOPS
USE YOMMP0   , ONLY : LOPT_SCALAR
USE YOMDIAGVAR,ONLY : DIAG_4DVAR
USE YOMVCGL  , ONLY : NSAVEEV, NSAVEPC
USE YOMMODERR, ONLY : N_COUPLED_WINDOWS
USE YOMJQ    , ONLY : LSTATMERR
USE YOMSENS  , ONLY : LGRVOL   ,NJROPT   ,LBSENS
USE YOMVRTL  , ONLY : L131TL   ,LTLINT   ,LOBSTL   ,LDRYTL   ,LIDMODEL
USE YOMCOSJB , ONLY : LJBZERO, LJPZERO, LJHZERO, LJLZERO, LJTZERO
USE ALGORITHM_STATE_MOD, ONLY : SETUP_ALGORITHM_STATE,SET_NUPTRA,SET_MUPTRA,GET_NSIM4D,SET_NSIM4D,&
 &                              GET_NUPTRA,GET_MUPTRA,GET_ALGOR_TYPE
USE YOMVRTLX, ONLY : LMINI, L801TL
!     -------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KULOUT

!     -------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: J
INTEGER(KIND=JPIM) :: NUPTRA
REAL(KIND=JPRB) :: Z
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -------------------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"




!     -------------------------------------------------------------------------

#include "namvar.nam.h"
#include "namvrtl.nam.h"
#include "namsens.nam.h"

!     -------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('YOMVAR:SETUP_VAR',0,ZHOOK_HANDLE)

!      ----------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!              -------------------

!        1.1 YOMJQ variables.

LSTATMERR = .FALSE.
!       random perturbations
LBGPERT = .FALSE.
DELTA   = 1.
MSIME   = 0



!        1.3 YOMMODERR variables.

N_COUPLED_WINDOWS=0

!        1.4 YOMVCGL variables.

NSAVEEV=5
NSAVEPC=3

!        1.5 YOMSENS variables.

LGRVOL=.FALSE.
NJROPT=1
LBSENS=.FALSE.

!        1.6 YOMVRTL variables.

L131TL=.FALSE.
LOBSTL=.FALSE.
LTLINT=.FALSE.
LDRYTL=.FALSE.
LIDMODEL=.FALSE.

!        1.7 YOMCOSJB variables.

LJBZERO=.FALSE.
LJPZERO=.FALSE.
LJHZERO=.FALSE.
LJLZERO=.FALSE.
LJTZERO=.FALSE.

!        1.8 YOMVAR variables linked to minimisation.

! used for M1QN3:
NMIMP=4
RXMIN=SQRT(EPSILON(Z))

! used for N1CG1:
LN1CG1=.FALSE.
ZEPSNEG=1.E-8_JPRB ! not DOCTOR compliant name
N1IMP=5
NSELECT=1
NPRECO=0
NBFGSB=2

! used for CONGRAD:
LCONGRAD=.TRUE.
NPCVECS=10
LWRIEVEC=.FALSE.
NWRIEVEC=10
N_DIAGS_CONVERGENCE=1
N_DIAGS_EIGENVECS=1
LEVECCNTL=.TRUE.
LEVECGRIB=.FALSE.
L_CHECK_GRADIENT=.FALSE.
RTOL_CHECK_GRADIENT=1.E-2_JPRB
LMPCGL=.FALSE.
EVBCGL=1.E-1_JPRB
LAVCGL=.FALSE.

! used for all minimizers:
IF (LECMWF) THEN
  NITER_MIN=10
  RCVGE=0.03_JPRB
  L_INFO_CONVERGENCE=.TRUE.
ELSE
  NITER_MIN=0
  RCVGE=1.E-2_JPRB
  L_INFO_CONVERGENCE=.FALSE.
ENDIF
NITER=50
NSIMU=60
NINFRA=10 !  NINFRA WILL BE CHANGED IN SUIOMI ACCORDING TO LWARM
R_NORM_REDUCTION_ABORT_LEVEL=0.01_JPRB
L_CHECK_CONVERGENCE=.TRUE.
L_ABS_CONVERGENCE=.FALSE.
L_INFO_CONTENT=.FALSE.
N_INFO_CONTENT_METHOD=1
N_INFO_CONTENT_SEED=0

!        1.9 YOMVAR variables: other ones.

LTEST =.FALSE.
LSLADREP=.FALSE.
LREPRO4DVAR =.FALSE.

IF (LECMWF) THEN
  IF (NCONF/100 == 0 .OR. NCONF == 302 .OR. (NCONF==131 .AND. L4DVAR)) THEN
    LJC=.FALSE.
  ELSE
    LJC=.TRUE.
  ENDIF
ELSE
  IF (NCONF/100 == 0.OR.NCONF == 302) THEN
    LJC=.FALSE.
  ELSE
    LJC=.TRUE.
  ENDIF
ENDIF
IF (LECMWF) THEN
  IF (NCONF==131 .AND. L4DVAR) THEN
    LJCDFI=.TRUE.
    LUSEJCDFI=.TRUE.
  ELSE
    LJCDFI=.FALSE.
    LUSEJCDFI=.FALSE.
  ENDIF
ELSE
  LJCDFI=.FALSE.
  LUSEJCDFI=.FALSE.
ENDIF
LJCIMPACT=.FALSE.
LJBIMPACT=.FALSE.

LFCOBS=.FALSE.
LFCOBSTEST=.FALSE.
LREFINC=.FALSE.

LINITCV=.TRUE.
LMODERR=.FALSE.
LVARBC=.FALSE.
LPERTMRPAS=.FALSE.

L_FGNOTBG=.FALSE.
L_BGREC=.FALSE.
L_FGREC=.FALSE.
NFGREC_MIN=0
L_COMPBGDEP=.FALSE.
LWRICV=.FALSE.

LAEOLUSAMD=.FALSE.
FILTERFACTOR=1.0_JPRB
FILTEREXPO=2.0_JPRB
FILTERRESOL=255.0_JPRB
LTRREF=.FALSE.
LWRIBVEC=.FALSE.
LWRIBVEC_FULL=.FALSE.
LREABVEC=.FALSE.
LWRISIGB=.FALSE.
LWRISIGA=.FALSE.
LWRISIGF=.FALSE.
LWRISB_VPROF=.FALSE.
LZOWA=.TRUE.
LWREINI=.FALSE.
LTWANA=.FALSE.
LTWGRA=.FALSE.
LTWINC=.FALSE.
LMINI=.FALSE.
LSKIPMIN=.FALSE.
R_MAX_CNUM_PC=10.0_JPRB
NFGFCLEN=6
NSIM4DL=999
NDIAG=1
RDX=1.0E-12_JPRB
CALL SET_NSIM4D(0)
IF (LECMWF) THEN
  ALPHAG=100._JPRB
ELSE
  ALPHAG=1.E+6_JPRB
ENDIF
ALPHAV=1.0_JPRB

NFRREF=1
NREFTS(0:JPNRFT)=0
NFRANA=10
NANATS(0:JPNRFT)=0
NFRGRA=1000
IF (LECMWF) THEN
  NGRATS(0)=1
  NGRATS(1)=1
  NGRATS(2:JPNRFT)=0
ELSE
  NGRATS(0:JPNRFT)=0
ENDIF

LGRASCAL=.FALSE.
NUPTRA=99999
MUPTRA=0 ! 3DVAR value, must be updated for 4DVAR
NUPTRA_RANGE=0
NOBSCURE_NUMBER = 0

IF(NCONF == 131) THEN
  NBGVECS=50
  LBGTRUNC=.TRUE.
  NBGTRUNC=1000
  LBGOBS=.TRUE.
  LBGM=.FALSE.
  LANOBS=.FALSE.
ELSE
  NBGVECS=0
  LBGTRUNC=.FALSE.
  NBGTRUNC=1000
  LBGOBS=.FALSE.
  LBGM=.FALSE.
  LANOBS=.FALSE.
ENDIF
NHEVECS=0

NSSBGV_COUNT=0
NSSHEV_COUNT=0
LTOVSCV=(LOBSC1.OR.NCONF/100 == 1)
LTOVSREP=.FALSE.
LGCV=.FALSE.
LGCVJO=.FALSE.
NITERGCV=40

CFNSIGB='sigma_b'
CFNSIGA='sigma_a'
CFNSIGF='sigma_f'
CFNPCV ='precon'

! Variables of YOMVRTL/NAMVRTL
L131TL=.FALSE.
LOBSTL=.FALSE.
LTLINT=.FALSE.
LDRYTL=.FALSE.
L801TL=.FALSE.
LGRVOL=.FALSE.
NJROPT=1
LBSENS=.FALSE.

LIDMODEL=.FALSE.
L3DFGAT=.FALSE.

L_GUESS_RUNTIME=.TRUE.
NITER_GUESS_RUNTIME=10
NFREQ_GUESS_RUNTIME=5
NMEM_GUESS_RUNTIME =10

IF (LECMWF) THEN
  LSPINT=.TRUE.
ELSE
  LSPINT=.FALSE.
ENDIF

LAMV_REASSIGN_PASSIVE=.FALSE.
LAMV_HEIGHT_ADJUST=.FALSE.
LREO3_BCOR=.FALSE.
LCH4_BCOR=.FALSE.
LPROPTL=.FALSE.
LCLDSINK=.TRUE.
CTOPBGE=3.0_JPRB ! not DOCTOR compliant name
CAMTBGE=0.01_JPRB ! not DOCTOR compliant name
RCOEFCO2=1.0_JPRB
RCOEFCH4=1.0_JPRB
RCOEFCO=1.0_JPRB
RCOEFNO2=1.0_JPRB
RCOEFGO3=1.0_JPRB
LBACKGECV=.FALSE.
LBACKGE=.FALSE.
LBACKGERENORM=.FALSE.
LUSEWAVRENORM=.FALSE.
LCONSTANTZFCE=.FALSE.
LDIAG_LCT=.FALSE.
LFAMEMBERS=.TRUE.
LTWLCZ=.FALSE.
LTOY42=.FALSE.
LSUSPQLIM=(.NOT.LELAM)
NDATE_TIME_WINDOW_END=0
LMONITOR_FCDEPAR=.FALSE.
LENDA=.FALSE.
LCHRESINCR = .FALSE.
LINC_TOVSCV = .TRUE.
LECV = .FALSE.

!      ----------------------------------------------------------------

!*       2.    Modifies default values: read namelists.
!              ----------------------------------------

CALL POSNAM(NULNAM,'NAMVAR')
READ(NULNAM,NAMVAR)

CALL POSNAM(NULNAM,'NAMVRTL')
READ(NULNAM,NAMVRTL)

CALL POSNAM(NULNAM,'NAMSENS')
READ(NULNAM,NAMSENS)

!      ----------------------------------------------------------------

!*       3.    Checkings and resettings.
!              -------------------------


!IF (.NOT.LSLAG) THEN
!  LVECADIN=.TRUE.
!ENDIF

!        3.2 YOMVRTL variables.

! resetting
IF (NCONF/100 == 8) THEN
  L131TL=.FALSE.
ENDIF

! check: OBSTL requires L131TL
IF (LOBSTL .AND. .NOT.L131TL) CALL ABOR1 (' SETUP_VAR: LOBSTL requires L131TL')

!        3.3 YOMVAR variables linked to minimisation.

! check: LN1CG1 and LCONGRAD mutually exclusive
IF (LN1CG1 .AND. LCONGRAD) CALL ABOR1(' SETUP_VAR: LN1CG1 and LCONGRAD are mutually exclusive')

! check: L_.._CONVERGENCE
IF (L_INFO_CONVERGENCE .AND..NOT. LCONGRAD) CALL ABOR1 (' SETUP_VAR: L_INFO_CONVERGENCE requires LCONGRAD')
IF (L_ABS_CONVERGENCE .AND. L_INFO_CONVERGENCE) &
 & CALL ABOR1 (' SETUP_VAR: L_ABS_CONVERGENCE and L_INFO_CONVERGENCE mutually exclusive')

! check: LAVCGL requires config 131
IF (LAVCGL .AND. NCONF /= 131) CALL ABOR1 (' SETUP_VAR: LAVCGL requires config 131')

! check: LAVCGL requires CONGRAD
IF (LAVCGL .AND..NOT. LCONGRAD) CALL ABOR1 (' SETUP_VAR: LAVCGL requires LCONGRAD')

!        3.4 YOMVAR variables: other variables.

! check: testings and resettings on LREFINC
!  We can notice that if (LECMWF .AND. LFCOBS .AND. LTEST) there is no possible value for LREFINC
IF (LTEST) THEN
  WRITE(KULOUT,'('' SETUP_VAR: LTEST = T => LREFINC = F'')')
  LREFINC=.FALSE.
ENDIF

IF (NCONF == 131 .AND. LBGPERT) LBACKG=.TRUE.
! ky: one must choose between these two redundant testings.
! resetting
IF (LECMWF .AND. LFCOBS) THEN
  WRITE(KULOUT,'('' SETUP_VAR: LECMWF and LFCOBS => LREFINC = T'')')
  LREFINC=.TRUE.
ENDIF
! check: standard use of forecast sensitivity to observations at ECMWF
IF (LECMWF .AND. LFCOBS .AND. (.NOT.LREFINC)) THEN
  CALL ABOR1(' SETUP_VAR: LECMWF=T and LFCOBS requires LREFINC=T.')
ENDIF

IF (LBGPERT .AND. NCONF /= 131) CALL ABOR1('LBGPERT requires c131')

IF (LBGPERT .AND. NBGVECS == 0)&
 & CALL ABOR1('LBGPERT requires randomozation, but NBGVECS=0')


! resetting: make sure AN diagnostics switches don't interfere with other configurations.
IF (NCONF/100 /= 1) THEN
  NFRANA=1
  NANATS(0:JPNRFT)=0
  NFRGRA=1
  NGRATS(0:JPNRFT)=0
ENDIF

! check: LBACKG (in YOMCT0/NAMCT0) must be T in some cases:
IF (NCONF/100 == 8 .AND. LBSENS .AND. (.NOT.LBACKG)) THEN
  CALL ABOR1(' SETUP_VAR: Conf 8xx with LBSENS=T requires LBACKG=T.')
ENDIF

! check: NUPTRA must be specified in the namelist for security
IF ( (LOBSC1.OR.(NCONF==131)) .AND. (NUPTRA==99999) )&
 & CALL ABOR1(' SETUP_VAR: NUPTRA must be specified in NAMVAR (e.g.-1)')
IF ( (LOBSC1.OR.(NCONF==131)) .AND. (NUPTRA>MUPTRA) )&
 & CALL ABOR1(' SETUP_VAR: NUPTRA must be smaller or equal to MUPTRA')

! check: LBGOBS and LBGM require randomisation
IF(LBGOBS .AND. NBGVECS == 0)&
 & CALL ABOR1(' SETUP_VAR: LBGOBS requires randomisation, but NBGVECS=0')
IF(LBGM .AND. NBGVECS == 0)&
 & CALL ABOR1(' SETUP_VAR: LBGM requires randomisation, but NBGVECS=0')
!*    check: so far LBGOBS or LBGM not implemented under LDIAG_LCT
IF (LDIAG_LCT .AND. LBGOBS) &
     & CALL ABOR1 ('LDIAG_LCT and LBGOBS are mutually exclusive')
IF (LDIAG_LCT .AND. LBGM) &
     & CALL ABOR1 ('LDIAG_LCT and LBGM are mutually exclusive')

! check: all components needed to make 4D-Var reproducible
LSLADREP = LSLADREP .OR. LREPRO4DVAR
LTOVSREP = LTOVSREP .OR. LREPRO4DVAR

! additional set-up
LFDBERR = LAVCGL.AND.LECMWF.AND.LFDBOP

!        3.5 YOMDIAGVAR variables.

DIAG_4DVAR%NOBS=0
DIAG_4DVAR%NCOST=0
DIAG_4DVAR%JO=0.0_JPRB
DIAG_4DVAR%JB=0.0_JPRB
DIAG_4DVAR%JC=0.0_JPRB
DIAG_4DVAR%JQ=0.0_JPRB
DIAG_4DVAR%JP=0.0_JPRB
DIAG_4DVAR%JH=0.0_JPRB
DIAG_4DVAR%JCVARBC=0.0_JPRB

! New way of tracking state of Algorithm (shared with OOPS)
IF (L_OOPS) THEN
  CALL SETUP_ALGORITHM_STATE('OOPS')
ELSE
  CALL SETUP_ALGORITHM_STATE('IFS')
ENDIF
CALL SET_NUPTRA(NUPTRA)
CALL SET_MUPTRA(MUPTRA)

WRITE(KULOUT,*) ' SETUP_VAR ALGORITHM: ', GET_ALGOR_TYPE()
WRITE(KULOUT,*) ' SETUP_VAR    NUPTRA: ', GET_NUPTRA()
WRITE(KULOUT,*) ' SETUP_VAR    MUPTRA: ', GET_MUPTRA()

!      -----------------------------------------------------------

!*       4.    Print final values.
!              -------------------

WRITE(KULOUT,'(''  '')')
WRITE(KULOUT,'('' --- PRINTINGS IN SETUP_VAR: '')')

!        4.1 YOMJQ variables.

WRITE(KULOUT,'('' MODULE YOMJQ'')')
WRITE(KULOUT,'(''  LSTATMERR= '',L2)') LSTATMERR

!        4.2 YOMMODERR variables.

WRITE(KULOUT,'('' MODULE YOMMODERR'')')
WRITE(KULOUT,'(''  N_COUPLED_WINDOWS = '',I6)') N_COUPLED_WINDOWS


!        4.4 YOMVCGL variables.
WRITE(KULOUT,'('' MODULE YOMVCGL'')')
WRITE(KULOUT,'(''  NSAVEEV = '',I6)') NSAVEEV
WRITE(KULOUT,'(''  NSAVEPC = '',I6)') NSAVEPC

!        4.5 YOMSENS variables.

WRITE(KULOUT,'('' MODULE YOMSENS'')')
WRITE(KULOUT,'(''  LGRVOL = '',L2,'' NJROPT = '',I6)') LGRVOL,NJROPT
WRITE(KULOUT,'(''  LBSENS = '',L2)')LBSENS

!        4.6 YOMVRTL variables.

WRITE(KULOUT,'('' MODULE YOMVRTL '')')
WRITE(KULOUT,*) ' L131TL=',L131TL,' LOBSTL=',LOBSTL
WRITE(KULOUT,*) ' LDRYTL=',LDRYTL
WRITE(KULOUT,*) ' LTLINT=',LTLINT
WRITE(KULOUT,*) ' LIDMODEL=',LIDMODEL

!        4.7 YOMCOSJB variables.

WRITE(KULOUT,'('' MODULE YOMCOSJB '')')
WRITE(KULOUT,'(''  LJBZERO = '',L2)') LJBZERO
WRITE(KULOUT,'(''  LJPZERO = '',L2)') LJPZERO
WRITE(KULOUT,'(''  LJHZERO = '',L2)') LJHZERO
WRITE(KULOUT,'(''  LJTZERO = '',L2)') LJTZERO
WRITE(KULOUT,'(''  LJLZERO = '',L2)') LJLZERO

!        4.8 YOMVAR variables linked to minimisation.

WRITE(KULOUT,'('' MODULE YOMVAR: variables linked to minimisation '')')

! used for M1QN3:
WRITE(KULOUT,'(''  NMIMP = '',I6)') NMIMP
WRITE(KULOUT,'(''  RXMIN = '',E13.7)') RXMIN

! used for N1CG1:
WRITE(KULOUT,'(''  LN1CG1 = '',L2)') LN1CG1
WRITE(KULOUT,'(''  ZEPSNEG = '',E13.7)') ZEPSNEG
WRITE(KULOUT,'(''  N1IMP = '',I6)') N1IMP
WRITE(KULOUT,'(''  NSELECT = '',I6)') NSELECT
WRITE(KULOUT,'(''  NPRECO = '',I6)') NPRECO
WRITE(KULOUT,'(''  NBFGSB = '',I6)') NBFGSB

! used for CONGRAD:
WRITE(KULOUT,'(''  LCONGRAD = '',L2)') LCONGRAD
WRITE(KULOUT,'(''  NPCVECS = '',I6)') NPCVECS
WRITE(KULOUT,'(''  LWRIEVEC = '',L2)') LWRIEVEC
WRITE(KULOUT,'(''  NWRIEVEC = '',I9)') NWRIEVEC
WRITE(KULOUT,'(''  N_DIAGS_CONVERGENCE = '',I6)') N_DIAGS_CONVERGENCE
WRITE(KULOUT,'(''  N_DIAGS_EIGENVECS = '',I6)') N_DIAGS_EIGENVECS
WRITE(KULOUT,'(''  LEVECCNTL = '',L2)') LEVECCNTL
WRITE(KULOUT,'(''  LEVECGRIB = '',L2)') LEVECGRIB
WRITE(KULOUT,'(''  L_CHECK_GRADIENT = '',L2)') L_CHECK_GRADIENT
WRITE(KULOUT,'(''  RTOL_CHECK_GRADIENT = '',E13.7)') RTOL_CHECK_GRADIENT
WRITE(KULOUT,'(''  LMPCGL = '',L2)') LMPCGL
WRITE(KULOUT,'(''  EVBCGL = '',E13.7)') EVBCGL
WRITE(KULOUT,'(''  LAVCGL = '',L2)') LAVCGL

! used for all minimizers:
WRITE(KULOUT,'(''  NITER = '',I6)') NITER
WRITE(KULOUT,'(''  NITER_MIN = '',I6)') NITER_MIN
WRITE(KULOUT,'(''  NSIMU = '',I6)') NSIMU
WRITE(KULOUT,'(''  NINFRA = '',I6)') NINFRA
WRITE(KULOUT,'(''  RCVGE = '',E13.7)') RCVGE
WRITE(KULOUT,'(''  R_NORM_REDUCTION_ABORT_LEVEL = '',E13.7)') R_NORM_REDUCTION_ABORT_LEVEL
WRITE(KULOUT,'(''  L_ABS_CONVERGENCE = '',L2)') L_ABS_CONVERGENCE
WRITE(KULOUT,'(''  L_INFO_CONVERGENCE = '',L2)') L_INFO_CONVERGENCE
WRITE(KULOUT,'(''  L_CHECK_CONVERGENCE = '',L2)') L_CHECK_CONVERGENCE
WRITE(KULOUT,'(''  L_INFO_CONTENT = '',L2)') L_INFO_CONTENT
WRITE(KULOUT,'(''  N_INFO_CONTENT_METHOD = '',I6)') N_INFO_CONTENT_METHOD
WRITE(KULOUT,'(''  N_INFO_CONTENT_SEED = '',I6)') N_INFO_CONTENT_SEED

!        4.9 YOMVAR variables: other ones.

WRITE(KULOUT,'('' MODULE YOMVAR: other variables '')')
WRITE(KULOUT,'(''  LTWANA = '',L2)') LTWANA
WRITE(KULOUT,'(''  LTWGRA = '',L2)') LTWGRA
WRITE(KULOUT,'(''  LTWINC = '',L2)') LTWINC
WRITE(KULOUT,'(''  R_MAX_CNUM_PC = '',E13.7)') R_MAX_CNUM_PC
WRITE(KULOUT,'(''  NFGFCLEN = '',I6)') NFGFCLEN
WRITE(KULOUT,'(''  NHEVECS = '',I6)') NHEVECS
WRITE(KULOUT,'(''  NSSBGV_COUNT = '',I6)') NSSBGV_COUNT
WRITE(KULOUT,'(''  NSSHEV_COUNT = '',I6)') NSSHEV_COUNT
WRITE(KULOUT,'(''  LGCV = '',L2)') LGCV
WRITE(KULOUT,'(''  LGCVJO = '',L2)') LGCVJO
WRITE(KULOUT,'(''  NITERGCV = '',I6)') NITERGCV
WRITE(KULOUT,'(''  CFNSIGB = '',A)') CFNSIGB
WRITE(KULOUT,'(''  CFNSIGA = '',A)') CFNSIGA
WRITE(KULOUT,'(''  CFNSIGF = '',A)') CFNSIGF
WRITE(KULOUT,'(''  CFNPCV = '',A)') CFNPCV
WRITE(KULOUT,'(''  L3DFGAT = '',L2)') L3DFGAT
WRITE(KULOUT,'(''  L_GUESS_RUNTIME = '',L2)') L_GUESS_RUNTIME
WRITE(KULOUT,'(''  NITER_GUESS_RUNTIME = '',I6)') NITER_GUESS_RUNTIME
WRITE(KULOUT,'(''  NFREQ_GUESS_RUNTIME = '',I6)') NFREQ_GUESS_RUNTIME
WRITE(KULOUT,'(''  NMEM_GUESS_RUNTIME = '',I6)') NMEM_GUESS_RUNTIME
WRITE(KULOUT,'(''  LAMV_REASSIGN_PASSIVE = '',L2)') LAMV_REASSIGN_PASSIVE
WRITE(KULOUT,'(''  LAMV_HEIGHT_ADJUST = '',L2)') LAMV_HEIGHT_ADJUST
WRITE(KULOUT,'(''  LREO3_BCOR = '',L2)') LREO3_BCOR
WRITE(KULOUT,'(''  LCH4_BCOR = '',L2)') LCH4_BCOR
WRITE(KULOUT,'(''  LCLDSINK = '',L2)') LCLDSINK
WRITE(KULOUT,'(''  CTOPBGE = '',E13.7)') CTOPBGE
WRITE(KULOUT,'(''  CAMTBGE = '',E13.7)') CAMTBGE
WRITE(KULOUT,'(''  LTWLCZ = '',L2)') LTWLCZ
WRITE(KULOUT,'(''  LTOY42 = '',L2)') LTOY42
WRITE(KULOUT,'(''  LSUSPQLIM = '',L2)') LSUSPQLIM
WRITE(KULOUT,'(''  LFDBERR = '',L2)') LFDBERR
WRITE(KULOUT,'(''  LENDA = '',L2)') LENDA
WRITE(KULOUT,'(''  LTEST = '',L2,'' LTRREF = '',L2,'' LZOWA = '',L2,'' LWREINI = '',L2)')&
 & LTEST,LTRREF,LZOWA,LWREINI
WRITE(KULOUT,'(''  LJC = '',L2, &
 & '' LJCDFI = '',L2,'' LUSEJCDFI = '',L2,'' LJCIMPACT = '',L2,'' LJBIMPACT = '',L2, &
 & '' ALPHAG = '',E13.7,'' ALPHAV = '',E13.7)') &
 & LJC,LJCDFI,LUSEJCDFI,LJCIMPACT,LJBIMPACT,ALPHAG,ALPHAV
WRITE(KULOUT,'(''  LGRASCAL = '',L2)') LGRASCAL
WRITE(KULOUT,'(''  LWRISIGB= '',L2,'' LWRISIGA= '',L2,'' LWRISIGF= '',L2,'' LWRISB_VPROF= '',L2, &
 & '' LWRIBVEC= '',L2,'' LWRIBVEC_FULL= '',L2,'' LREABVEC= '',L2)')&
 & LWRISIGB,LWRISIGA,LWRISIGF,LWRISB_VPROF,LWRIBVEC,LWRIBVEC_FULL,LREABVEC
WRITE(KULOUT,'(''  NSIM4D = '',I6,'' NDIAG = '',I6,'' NSIM4DL = '',I6)') GET_NSIM4D(),NDIAG,NSIM4DL
WRITE(KULOUT,'(''  RDX  = '',E13.7)') RDX
WRITE(KULOUT,'(''  NFRREF = '',I6,'' NFRANA = '',I6,'' NFRGRA = '',I6)') NFRREF,NFRANA,NFRGRA
WRITE(KULOUT,*) ' NREFTS = ',NREFTS(0),(NREFTS(J),J=1,ABS(NREFTS(0)))
WRITE(KULOUT,*) ' NANATS = ',NANATS(0),(NANATS(J),J=1,ABS(NANATS(0)))
WRITE(KULOUT,*) ' NGRATS = ',NGRATS(0),(NGRATS(J),J=1,ABS(NGRATS(0)))
WRITE(KULOUT,'(''  NUPTRA = '',I6,'' MUPTRA = '',I6,'' NBGVECS= '',I6, &
 & '' LBGTRUNC= '',L2,'' NBGTRUNC= '',I6, &
 & '' LBGOBS= '',L2,'' LBGM= '',L2,'' LANOBS= '',L2)')&
 & NUPTRA,MUPTRA,NBGVECS,LBGTRUNC,NBGTRUNC,LBGOBS,LBGM,LANOBS
WRITE(KULOUT,'(''  NUPTRA_RANGE = '',I6)') NUPTRA_RANGE
WRITE(KULOUT,'(''  LTOVSCV = '',L2)') LTOVSCV
WRITE(KULOUT,'(''  FILTERFACTOR = '',E13.7)') FILTERFACTOR
WRITE(KULOUT,'(''  FILTEREXPO = '',E13.7)') FILTEREXPO
WRITE(KULOUT,'(''  FILTERRESOL = '',E13.7)') FILTERRESOL
WRITE(KULOUT,'(''  LINITCV= '',L2,'' LMODERR= '',L2)') LINITCV,LMODERR
WRITE(KULOUT,'(''  LVARBC= '',L2)') LVARBC
WRITE(KULOUT,'(''  LPERTMRPAS= '',L2)') LPERTMRPAS
WRITE(KULOUT,'(''  LAEOLUSAMD= '',L2)') LAEOLUSAMD
WRITE(KULOUT,'(''  LREPRO4DVAR= '',L2,'' LSLADREP= '',L2,'' LTOVSREP= '',L2)') LREPRO4DVAR,LSLADREP,LTOVSREP
WRITE(KULOUT,'(''  LSKIPMIN= '',L2)') LSKIPMIN
WRITE(KULOUT,'(''  LFCOBS= '',L2,'' LFCOBSTEST= '',L2)') LFCOBS, LFCOBSTEST
WRITE(KULOUT,'(''  LREFINC= '',L2)') LREFINC
WRITE(KULOUT,'(''  LSPINT = '',L2)') LSPINT
WRITE(KULOUT,'(''  LPROPTL = '',L2)') LPROPTL
WRITE(KULOUT,'(''  LBACKGE = '',L2)') LBACKGE
WRITE(KULOUT,'(''  LBACKGECV = '',L2)') LBACKGECV
WRITE(KULOUT,'(''  LDIAG_LCT = '',L2)') LDIAG_LCT
WRITE(KULOUT,'(''  L_FGNOTBG= '',L2,'' LWRICV= '',L2,''  L_BGREC= '',L2,''  L_FGREC= '',L2,''  L_COMPBGDEP= '',L2)') L_FGNOTBG, LWRICV, L_BGREC, L_FGREC, L_COMPBGDEP
WRITE(KULOUT,'(''  NFGREC_MIN = '',I6)') NFGREC_MIN
WRITE(KULOUT,*) ' DATE AND TIME OF 4DVAR WINDOW END=',NDATE_TIME_WINDOW_END
WRITE(KULOUT,*) ' LMONITOR_FCDEPAR=',LMONITOR_FCDEPAR
IF (LMONITOR_FCDEPAR) WRITE(KULOUT,*) ' Forecast departure configuration LMONITOR_FCDEPAR ON'

!      -----------------------------------------------------------

!*       5.    Additional set-up.
!              ------------------

! The call to SUIOMI could be done there, because it sets-up YOMIOMI variables
! containing additional variables about minimisations.





!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('YOMVAR:SETUP_VAR',1,ZHOOK_HANDLE)
END SUBROUTINE SETUP_VAR

END MODULE YOMVAR
