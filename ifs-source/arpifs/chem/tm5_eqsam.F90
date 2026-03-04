SUBROUTINE TM5_EQSAM(PT,PRH,PNH,PNO3,PSO4,PYEQ)


!**   DESCRIPTION 
!     ----------
!
!   Part of TM5 routines for IFS chemistry: 
!--------------------------------------------------------------------------
!   Eulerian backward Iteration
!   Chemistry solver for the CBM4 scheme 
!--------------------------------------------------------------------------
!
!
!
!**   INTERFACE.
!     ----------
!          *TM5_EQSAM* IS CALLED FROM *CHEM_tm5*.

! INPUTS:
! -------
! PT                     : Temperature[K]
! PRH                    : relative humidity [0-1]
! PNH                    : initial concentrations of NH3+NH4          (mol/m3)
! PNO3                   : initial concentrations of HNO3+NO3_A       (mol/m3)
! PSO4                   : initial concentrations of SO4              (mol/m3)
! 
!
!
! OUTPUTS:
! -------
! PYEQ (4)               : final   concentrations of selected tracers  ( mol/m3)
!
! LOCAL:
! -------
!
!
!     AUTHOR.
!     -------
!        Swen Metzger *MPCH-Mainz*
!        Implemented in C-IFS by Vincent Huijnen   *KNMI* 
!        TM5-community    
!
!     MODIFICATIONS.
!     --------------
!        First implementation : 2011-06-16


!______________________________________________________________________________________________________
!      Originally written by Swen Metzger 3/11/99. Modified 2002, 2003.
!
!      Department of Atmospheric Chemistry, Max-Planck-Institute for Chemistry.
!      email: metzger@mpch-mainz.mpg.de
!      http://www.mpch-mainz.mpg.de/~metzger
!
!      COPYRIGHT 1999-2003
!
!      purpose
!      -------
!      EQSAM (EQuilibrium Simplified Aerosol Model) is a new and computationally efficient thermodynamic
!      aerosol composition model that allows to calculate the gas/aerosol equilibrium partitioning,
!      including aerosol water, sufficiently fast and accurate for global (or even regional) modeling.
!      EQSAM is based on a number of parameterizations, including single solute molalities and activity 
!      coefficients (AC). The thermodynamic framework (domains and subdomains, internally mixed aerosols) 
!      is the same as of more sophisticated thermodynamic equilibrium models (EQMs), e.g. of ISORROPIA 
!      (Nenes et al., 1998). Details are given in the references below (and the references therein).
!
!      The main assumption on which EQSAM/EQMs are based is thermodynamical and chemical equilibrium. 
!      From this assumption it directly follows that the aerosol water activity (aw) equals the ambient 
!      relative humidity (RH), if the water vapor pressure is sufficiently larger than the partial vapor
!      pressure of the aerosol compounds. This is approximately true for tropospheric aerosols. Given the 
!      large amount of water vapor present, water vapor and aerosol water equilibrate relatively faster 
!      compared to all other aerosol compounds. This is subsequently also true for single aerosol compounds.
!      The water activity of single solutes must also equal RH under this assumption. Therefore, the so 
!      called ZSR-relation is (and can be) used to calculate the aerosol associated water mass (simply
!      from the sum of all water mass fractions that are derived from measured single solute molalities). 
!
!      In contrast to other EQMs, EQSAM utilizes the fact that the RH fixes the water activity 
!      (under the above assumptions) and the consequence that any changes in RH also causes changes in 
!      the aerosol water mass and, hence, aerosol activity (including activity coefficients). Thus, an decrease
!      (increase) in RH decrease (increases) the aerosol water mass (and water activity). This can change the
!      aerosol composition, e.g. due to condensation (evaporation/crystallization), because the vapor pressure 
!      above the aerosol reduces (increases). In turn, a vapor pressure reduction (increase) due to changes
!      in the aerosol composition is compensated by an associated condensation (evaporation) of water vapor 
!      to maintain the aerosol molality to remain constant (because aw=RH). Furthermore, the aerosol water 
!      mainly depends on the aerosol mass and the type of solute, so that parameterizations of single solute 
!      molalities and activity coefficients can be defined, only depending on the type of solute and RH.
!      The advantage of using such parameterizations is that the entire aerosol equilibrium composition 
!      can be solved analytically, i.e. non-iteratively, which considerably reduces the amount of CPU time 
!      that is usually need for aerosol thermodynamic calculations (especially if an EQM is incorporated in
!      an aerosol dynamical model that is in turn embedded in a high resolution regional or global model).
!
!      However, EQSAM should still be regarded as a starting point for further developments. There is still 
!      room for improvements. For instance, this code is not yet numerically optimized (vectorized) and a 
!      number of improvements with respect to an explicit treatment of additional equilibrium reactions,
!      missing (or only implicit) dissociation, and a basic parameterization of the water uptake. 
!      
!      Note that EQSAM was originally developed to calculate the gas/aerosol equilibrium partitioning of the 
!      ammonium-sulfate-nitrate-water system for climate models, excluding solid compounds. 
!      This version (eqsam_v03d.f90) is extended with respect to sea salt. Solids/hysteresis are treated in a 
!      simplified manner. Results of a box model comparison with ISORROPIA will be available from the web page.
!      Please also note that the water uptake is based on additional (unpublished) parameterizations for single 
!      solute molalities, which are derived from tabulated measurements used in ISORROPIA. Note further that 
!      this extended version (eqsam_v03d.f90) is not yet published. A publication is in progress.
!
!
! Version History:
!
!  eqsam_CIFS (Vincent Huijnen, KNMI, June 2011):
!   - Na, Po and Cl dependency is removed
!   - redundant calculations have been removed
!   - input / output arrays have been made as small as possible 
!
!  eqsam_v03d.f90 (MPI-CH, June 2003): 
!   - gama parameterizations now according to Metzger 2002 (JGR Appendix)
!   - improved pH calculations (still restricted to strong acids)
!   - removed bug that lead to too high nitrate formation at dry and cold regions (UT/LS) 
!   - removed bug in solid/hysteresis calculations 
!     (both bugs introduced in eqsam_v03b.f90 by cleaning up eqsam_v02a.f90)
!   
!  eqsam_v03c.f90 (MPI-CH, April 2003):
!   - more accurate paramterizations of single solute molalities (Na, Cl species)
!   - cleanded up RHD subdomain structure
!   - improved water uptake (Na, Cl species)
!
!      
!   
!      method
!      ------
!      equilibrium / internal mixture assumption / aw=rh
!      System: NH3,NH4+/H2SO4+,HSO4-,SO4--/HNO3,NO3-, HCl,Cl-/Na+, H2O 
!        (K+,Ca++,Mg++)
!      
!      references
!      ---------
!      Swen Metzger Ph.D Thesis, University Utrecht, 2000.
!   http://www.library.uu.nl/digiarchief/dip/diss/1930853/inhoud.htm
!
!      Metzger, S. M., F. J. Dentener, J. Lelieveld, and S. N. Pandis, 
!   GAS/AEROSOL PARTITIONING I: A COMPUTATIONALLY EFFICIENT MODEL, 
!   J Geophys. Res., 107, D16, 10.1029/2001JD001102, 2002
!   http://www.agu.org/journals/jd/jd0216/2001JD001102/index.html
!      Metzger, S. M., F. J. Dentener, A. Jeuken, and M. Krol, J. Lelieveld, 
!   GAS/AEROSOL PARTITIONING II: GLOBAL MODELING RESULTS, 
!   J Geophys. Res., 107, D16, 10.1029/2001JD001103, 2002.
!   http://www.agu.org/journals/jd/jd0216/2001JD001103/index.html
!_________________________________________________________________________________________





USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
! USE YOMLUN   , ONLY : NULOUT, NULERR 

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------


REAL(KIND=JPRB),INTENT(IN)    :: PT,PRH,PNH,PNO3,PSO4
REAL(KIND=JPRB),INTENT(OUT)   :: PYEQ(4)    ! final concentrations

! * LOCAL 

REAL(KIND=JPRB),PARAMETER    :: ZT0=298.15,ZR=82.0567E-6 ! in cu.m*atm/deg/mole
REAL(KIND=JPRB),PARAMETER    :: ZMWH20=55.51*18.01, ZERO=0.0_JPRB
REAL(KIND=JPRB),PARAMETER    :: ZGF1=0.25,ZK=2.    ! exponents of AC-RH functions
!
REAL(KIND=JPRB)           :: ZT0T,ZKAN,ZGAMA,ZGG,ZGF,ZGFN
REAL(KIND=JPRB)           :: ZX0,ZX1,ZX2,ZX3,ZX4,ZX5,ZX6,ZXK10,ZXK6
REAL(KIND=JPRB)           :: ZFLAG,ZZKAN,ZCOEF
REAL(KIND=JPRB)           :: ZPNH4,ZPNO3,ZGNO3,ZGNH3
REAL(KIND=JPRB),DIMENSION(2) :: ZRHDZ ! RHD / MRHD arrays for different aerosol types
! salt solutes:
!   1 = NH4NO3, 2 = 2H-SO4
REAL(KIND=JPRB),DIMENSION(2),PARAMETER :: ZRHDA = (/0.67500, 0.4000  /)
REAL(KIND=JPRB),DIMENSION(2),PARAMETER :: ZRHDE = (/ 262.000, 384.00 /)



REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TM5_EQSAM',0,ZHOOK_HANDLE )

ZFLAG=1._JPRB

! SULFATE RICH

IF(PNH <= (2._JPRB*PSO4))  ZFLAG=3._JPRB

! SULFATE VERY RICH CASE if NH4/SO4 < 1

IF(PNH <= PSO4) ZFLAG=4._JPRB

! SULFATE NEUTRAL CASE

!VH IF(PNH > PSO4) ZFLAG=2._JPRB
!VH Fix on 20 Jan 2017. 
IF(PNH > 2_JPRB*PSO4) ZFLAG=2._JPRB



! CALCULATE TEMPERATURE DEPENDENCY FOR SOME RHDs

ZRHDZ(:)=ZRHDA(:)*EXP(ZRHDE(:)*(1._JPRB/PT-1._JPRB/ZT0))

! ACCOUNT FOR VARIOUS AMMOMIUM/SODIUM SULFATE SALTS ACCORDING TO MEAN VALUE AS OF ISORROPIA

ZGG=2.0     ! (Na)2SO4 / (NH4)2SO4 IS THE PREFFERED SPECIES FOR SULFATE DEFICIENT CASES
IF(ZFLAG == 3._JPRB) THEN
   IF(PRH <= ZRHDZ(2_JPIM)) THEN  ! ACCOUNT FOR MIXTURE OF (NH4)2SO4(s) & NH4HSO4(s) & (NH4)3H(SO4)2(s) 
      ZGG=1.677_JPRB    !    (Na)2SO4 &  NaHSO4
      !    ZGG=1.5
   ELSEIF(PRH > ZRHDZ(2_JPIM).AND.PRH <= ZRHDZ(1)) THEN ! MAINLY (Na)2SO4 / (NH4)2SO4(s) & (NH4)3H(SO4)2(s)
      ZGG=1.75_JPRB
      !    GG=1.5
   ELSEIF(PRH > ZRHDZ(1_JPIM)) THEN   ! (NH4)2SO4(S) & NH4HSO4(S) & SO4-- & HSO4-
      ZGG=1.5_JPRB   !  (Na)2SO4 &  NaHSO4
   ENDIF
ENDIF
IF(ZFLAG == 4._JPRB) ZGG=1.0_JPRB   ! IF SO4 NEUTRALIZED, THEN ONLY AS NaHSO4 / NH4HSO4(S) OR  HSO4- / H2SO4

!VH Note: ZRHD=PRH

! GET WATER ACTIVITIES ACCORDING TO METZGER, 2000.
! FUNCTION DERIVED FROM ZSR RELATIONSHIP DATA (AS USED IN ISORROPIA)

! CALCULATE TEMPERATURE DEPENDENT EQUILIBRIUM CONSTANTS

ZT0T=ZT0/PT
ZCOEF=1.0+LOG(ZT0T)-ZT0T

! EQUILIBRIUM CONSTANT NH4NO3(s) <==> NH3(g) + HNO3(g) [atm^2] (ISORROPIA)

ZXK10= 5.746E-17_JPRB * EXP(-74.38*(ZT0T-1.0_JPRB) + 6.120_JPRB*ZCOEF)
ZKAN = ZXK10/(ZR*PT)/(ZR*PT)

! EQUILIBRIUM CONSTANT  NH4CL(s) <==> NH3(g) + HCL(g) [atm^2] (ISORROPIA)

!ZXK6  = 1.086e-16
ZXK6 = 1.086E-16_JPRB * EXP(-71.00_JPRB*(ZT0T-1.0_JPRB) + 2.400_JPRB*ZCOEF)

! GET MEAN MOLAL IONIC ACTIVITY COEFF ACCORDING TO METZGER, 2002.

ZGAMA=(PRH**ZFLAG/(1000._JPRB/ZFLAG*(1._JPRB-PRH)+ZFLAG))
ZGAMA = ZGAMA**ZGF1     ! ONLY ZGAMA TYPE OF NH4NO3, NaCl, etc. NEEDED SO FAR

ZGFN=ZK*ZK        ! ZK=2, i.e. condensation of 2 water molecules per 1 mole ion pair
ZGF=ZGFN*ZGF1       ! = GFN[=Nw=4] * ZGF1[=(1*1^1+1*1^1)/2/Nw=1/4] = 1
! ONLY ZGAMA TYPE OF NH4NO3, NH4Cl, etc. needed so far

ZGAMA=PRH**ZGF/((ZGFN*ZMWH20*(1._JPRB/PRH-1._JPRB)))**ZGF1
ZGAMA = MIN(MAX(0._JPRB,ZGAMA),1.0_JPRB)     ! FOCUS ON 0-1 SCALE
ZGAMA = (1._JPRB-ZGAMA)**ZK      ! transplate into aqueous phase equillibrium and account for 
! enhanced uptake of aerosol precursor gases with increasing RH
! (to match the results of ISORROPIA)

! CALCULATE ZRHD DEPENDENT EQ: IF RH <  ZRHD => NH4NO3(s) <==> NH3 (g) + HNO3(g) (ISORROPIA)
!         IF RH >> ZRHD => HNO3  (g)   -> NO3 (aq)

ZX0   = MAX(ZERO,MIN(PNH,ZGG*PSO4))      ! MAX AMMOMIUM SULFATE
ZX1   = MAX(ZERO,MIN(PNH-ZX0,PNO3))      ! MAX AMMOMIUM NITRATE
!

ZX2   = MAX(PNH-ZX1-ZX0,ZERO)    ! INTERIM RESIDUAL NH3
ZX3   = MAX(PNO3-ZX1,ZERO)      ! INTERIM RESIDUAL HNO3
!
ZZKAN=2._JPRB*ZGAMA

ZX4   = ZX2 + ZX3
ZX5   = SQRT(ZX4*ZX4+ZKAN*ZZKAN*ZZKAN)
ZX6   = 0.5_JPRB*(-ZX4+ZX5)
ZX6   = MIN(ZX1,ZX6)

ZGNH3 = ZX2 + ZX6     ! INTERIM RESIDUAl NH3
ZGNO3 = ZX3 + ZX6     ! RESIDUAl HNO3

ZGNH3=MAX(0._JPRB,MIN(ZGNH3,PNH))
ZGNO3=MAX(0._JPRB,MIN(ZGNO3,PNO3))

!VH : quantity not needed to be evaluated:  ANH4 = PNH -  ZGNH3 ! aqueous phase ammonium [mol/m^3]
!VH : quantity not needed to be evaluated:  ANO3 = PNO3 - ZGNO3 ! aqueous phase nitrate  [mol/m^3]

! total particulate matter (neutralized)
ZPNH4=PNH-ZGNH3                    ! particulate ammonium [mol/m^3]
ZPNO3=MAX(0._JPRB,PNO3-ZGNO3)           ! particulate nitrate    [mol/m^3]

! solid matter - not needed in output, so not calculated. Note that of SNH4/SNO3 is positive
!this means a loss term in the concentrations!!!
!VH SNH4=PNH4-ANH4        ! solid phase ammonium   [mol/m^3]
!VH SNO3=ZPNO3-ANO3       ! solid phase nitrate    [mol/m^3]

PYEQ(1) = ZGNO3     ! residual HNO3  (g)      [mol/m^3]
PYEQ(2) = ZGNH3     ! residual NH3   (g)      [mol/m^3]
PYEQ(3) = ZPNH4    ! particulate ammonium (p=a+s)     [mol/m^3]
PYEQ(4) = ZPNO3    ! particulate nitrate  (p=a+s)     [mol/m^3]

IF (LHOOK) CALL DR_HOOK('TM5_EQSAM',1,ZHOOK_HANDLE )
END SUBROUTINE TM5_EQSAM





