! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWPHYS

!*    ** *YOWPHYS* - PARAMETERS FOR WAVE PHYSICS PARAMETERISATION 

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!     *XKAPPA*    VON KARMAN CONSTANT.
      REAL(KIND=JWRB), PARAMETER :: XKAPPA=0.40_JWRB
!     *XNLEV*     REAL      WIND SPEED REFERENCE HEIGHT
      REAL(KIND=JWRB), PARAMETER :: XNLEV=10.0_JWRB
!     *RNU*      KINEMATIC AIR VISCOSITY (when coupled get from the atmospheric model)
      REAL(KIND=JWRB) :: RNU
!     *RNUM*     REDUCED KINEMATIC AIR VISCOSITY FOR MOMENTUM TRANSFER (as RNU)
      REAL(KIND=JWRB) :: RNUM
!     *PRCHAR*   DEFAULT VALUE FOR CHARNOCK
      REAL(KIND=JWRB) :: PRCHAR



!     WIND INPUT ::
!     ==========
!     *BETAMAX*  PARAMETER FOR WIND INPUT.
      REAL(KIND=JWRB) :: BETAMAX

!     *BETAMAXOXKAPPA2*  BETAMAX/XKAPPA**2 
      REAL(KIND=JWRB) :: BETAMAXOXKAPPA2

!     *BMAXOKAP*         DELTA_THETA_RN * BETAMAXOXKAPPA2 /XKAPPA
      REAL(KIND=JWRB) :: BMAXOKAP

!     *BMAXOKAPDTH*      BMAXOKAP * DELTH
      REAL(KIND=JWRB) :: BMAXOKAPDTH

!     *GAMNCONST*           DELTA_THETA*0.5_JWRB*ZPI**4*GM1*3*BETAMAXOXKAPPA2/XKAPPA
      REAL(KIND=JWRB) :: GAMNCONST

!     *ZALP*      SHIFTS GROWTH CURVE.
      REAL(KIND=JWRB) :: ZALP

!     *ALPHA*     MINIMUM CHARNOCK CONSTANT WITH NO WAVES.
      REAL(KIND=JWRB) :: ALPHA
      REAL(KIND=JWRB) :: ALPHAMIN
!     MAXIMUM CHARNOCK
      REAL(KIND=JWRB), PARAMETER :: ALPHAMAX=0.11_JWRB

!     *CHNKMIN_U* WIND THRESHOLD USED TO REDUCED MIN CHARNOCK (see *CHNKMIN*)
      REAL(KIND=JWRB) :: CHNKMIN_U

!     *TAUWSHELTER* SHELTERING COEFFICIENT in Ardhuin et al. PHYSICS
      REAL(KIND=JWRB) :: TAUWSHELTER

!     MAXIMUM PHILLIPS PARAMETER USED TO CONTROL MAXIMUM STEEPNESS
      REAL(KIND=JWRB) :: ALPHAPMAX

!     MINIMUM PHILLIPS PARAMETER ALLOWED :  ALPHAPMINFAC/WAVE_AGE
      REAL(KIND=JWRB), PARAMETER :: ALPHAPMINFAC=0.1_JWRB
      REAL(KIND=JWRB) :: FLMINFAC 

!     DIRECTIONALITY CORRECTION FACTORS IN THE GOWTH RATE RENORMALISATION (SEE JANSSEN ECMWF TECH MEMO 845)
      REAL(KIND=JWRB) :: DELTA_THETA_RN
      REAL(KIND=JWRB) :: RN1_RN
!     if LLCAPCHNK, DELTA_THETA_RN is enhanced by factor 1+DTHRN_A*(1+TANH(U10-DTHRN_U))
!     This is intended to model the impact of unrepresented effects on the drag for very high winds
!     (i.e. spray, foam, etc...)
      REAL(KIND=JWRB) :: DTHRN_A
      REAL(KIND=JWRB) :: DTHRN_U

!     COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!     FREQUENCIES LE MAX(TAILFACTOR*MAX(FMNWS,FM),TAILFACTOR_PM*FPM),
!     WHERE FPM IS THE PIERSON-MOSKOWITZ FREQUENCY BASED ON FRICTION
!     VELOCITY. (FPM=G/(FRIC*ZPI*USTAR))
      REAL(KIND=JWRB) :: TAILFACTOR
      REAL(KIND=JWRB) :: TAILFACTOR_PM

!     FOR THE GRAVITY-CAPILLARY MODEL, ONE NEEDS TO SPECIFY ANGULAR ADJUSTMENT ANG_GC (SEE *SETWAVPHYS*)
!     ANG_GC = ANG_GC_A + ANG_GC_B * TANH(ANG_GC_C * USTAR**2) 
      REAL(KIND=JWRB) :: ANG_GC_A, ANG_GC_B, ANG_GC_C

!     Negative wind input, ARDHUIN et al. 2010:
      REAL(KIND=JWRB), PARAMETER :: SWELLF = 0.66_JWRB ! controls the turbulent swell dissipation
      REAL(KIND=JWRB), PARAMETER :: SWELLF2 = -0.018_JWRB
      REAL(KIND=JWRB), PARAMETER :: SWELLF3 = 0.022_JWRB
      REAL(KIND=JWRB) :: SWELLF4  !! See *SETWAVPHYS*
      REAL(KIND=JWRB), PARAMETER :: SWELLF5 = 1.2_JWRB ! controls the viscous swell dissipation
      REAL(KIND=JWRB), PARAMETER :: SWELLF6 = 1.0_JWRB
      REAL(KIND=JWRB) :: SWELLF7  !! See *SETWAVPHYS*
      REAL(KIND=JWRB) :: SWELLF7M1 !!!! set it to 1 if you decide to have SWELLF7=0
      REAL(KIND=JWRB) :: Z0RAT   !! See *SETWAVPHYS*
      REAL(KIND=JWRB) :: Z0TUBMAX   !! See *SETWAVPHYS*

      REAL(KIND=JWRB), PARAMETER :: ABMIN = 0.3_JWRB
      REAL(KIND=JWRB), PARAMETER :: ABMAX = 8.0_JWRB 


!     WHITECAP DISSIPATION ::
!     ====================


!     Whitecap dissipation WAM cycle 4:

      REAL(KIND=JWRB) :: CDIS
      REAL(KIND=JWRB) :: DELTA_SDIS
      REAL(KIND=JWRB) :: CDISVIS


!     Whitecap dissipation, ARDHUIN et al. 2010:
!     TEST 473:
!     Br:
      REAL(KIND=JWRB), PARAMETER :: SDSBR = 9.0E-4_JWRB

!     Saturation dissipation coefficient
      INTEGER(KIND=JWIM), PARAMETER :: ISDSDTH = 80_JWIM 
      INTEGER(KIND=JWIM), PARAMETER :: ISB=2_JWIM
      INTEGER(KIND=JWIM), PARAMETER :: IPSAT=2_JWIM

      REAL(KIND=JWRB), PARAMETER :: SSDSC2 = -2.2E-5_JWRB
      REAL(KIND=JWRB), PARAMETER :: SSDSC4 = 1.0_JWRB
      REAL(KIND=JWRB), PARAMETER :: SSDSC6 = 0.3_JWRB 
      REAL(KIND=JWRB), PARAMETER :: MICHE = 1.0_JWRB 


!     Cumulative dissipation coefficient
!!!      REAL(KIND=JWRB), PARAMETER :: SSDSC3 = -0.40344_JWRB 
!!!   This is quite an expensive computation. Setting it to 0 will disable its calculation.
!!!   It was found that if the high frequency tail is prescribed (see frcutindex),
!!!   then the results are very similar.
!!!   Converserly, it will be required, in particular for the wave modified fluxes to NEMO
!!!   when the high frequency tail is not prescribed and used in all calculation
      REAL(KIND=JWRB), PARAMETER :: SSDSC3 = 0.0_JWRB 
      REAL(KIND=JWRB), PARAMETER :: SSDSBRF1   = 0.5_JWRB
!     28.16 = 22.0 * 1.6² * 1/2 with  
!     22.0 (Banner & al. 2000, figure 6) 
!     1.6  the coefficient that transforms  SQRT(B) to Banner et al. (2000)'s epsilon
!     1/2  factor to correct overestimation of Banner et al. (2000)'s breaking probability due to zero-crossing analysis
      REAL(KIND=JWRB), PARAMETER :: BRKPBCOEF=28.16_JWRB  

!     Wave-turbulence interaction coefficient 
      REAL(KIND=JWRB) :: SSDSC5  !! See *SETWAVPHYS*

!     NSDSNTH is the number of directions on both used to compute the spectral saturation  
      INTEGER(KIND=JWIM) :: NSDSNTH
!     NDIKCUMUL is the  integer difference in frequency bands
      INTEGER(KIND=JWIM) :: NDIKCUMUL

      INTEGER(KIND=JWIM), ALLOCATABLE :: INDICESSAT(:,:)
      REAL(KIND=JWRB), ALLOCATABLE :: SATWEIGHTS(:,:)
! ----------------------------------------------------------------------
      END MODULE YOWPHYS
