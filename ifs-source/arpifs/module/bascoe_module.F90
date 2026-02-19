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

MODULE BASCOE_MODULE

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

  INTEGER(KIND=JPIM), PARAMETER  :: ntrace = 64  ! All tracers for chemistry
  !
  ! components numbers
  !
  INTEGER(KIND=JPIM), PARAMETER :: iO3      =1
  INTEGER(KIND=JPIM), PARAMETER :: iH2O2    =2 
  INTEGER(KIND=JPIM), PARAMETER :: iCH4     =3
  INTEGER(KIND=JPIM), PARAMETER :: iCO      =4
  INTEGER(KIND=JPIM), PARAMETER :: iHNO3    =5
  INTEGER(KIND=JPIM), PARAMETER :: iCH3OOH  =6 
  INTEGER(KIND=JPIM), PARAMETER :: iCH2O    =7
  INTEGER(KIND=JPIM), PARAMETER :: iNO      =8
  INTEGER(KIND=JPIM), PARAMETER :: iHO2     =9
  INTEGER(KIND=JPIM), PARAMETER :: iCH3    =10
  INTEGER(KIND=JPIM), PARAMETER :: iCH3O   =11
  INTEGER(KIND=JPIM), PARAMETER :: iHCO    =12
  INTEGER(KIND=JPIM), PARAMETER :: iCH3O2  =13 
  INTEGER(KIND=JPIM), PARAMETER :: iOH     =14
  INTEGER(KIND=JPIM), PARAMETER :: iNO2    =15
  INTEGER(KIND=JPIM), PARAMETER :: iN2O5   =16
  INTEGER(KIND=JPIM), PARAMETER :: iHO2NO2 =17
  INTEGER(KIND=JPIM), PARAMETER :: iNO3    =18
  INTEGER(KIND=JPIM), PARAMETER :: iN2O    =19
  INTEGER(KIND=JPIM), PARAMETER :: iH2O    =20
  INTEGER(KIND=JPIM), PARAMETER :: iOCLO   =21
  INTEGER(KIND=JPIM), PARAMETER :: iHCL    =22
  INTEGER(KIND=JPIM), PARAMETER :: iCLONO2 =23
  INTEGER(KIND=JPIM), PARAMETER :: iHOCL   =24
  INTEGER(KIND=JPIM), PARAMETER :: iCL2    =25
  INTEGER(KIND=JPIM), PARAMETER :: iHBR    =26
  INTEGER(KIND=JPIM), PARAMETER :: iBRONO2 =27 
  INTEGER(KIND=JPIM), PARAMETER :: iCL2O2  =28
  INTEGER(KIND=JPIM), PARAMETER :: iHOBR   =29
  INTEGER(KIND=JPIM), PARAMETER :: iBRCL   =30
  INTEGER(KIND=JPIM), PARAMETER :: iCFC11  =31
  INTEGER(KIND=JPIM), PARAMETER :: iCFC12  =32
  INTEGER(KIND=JPIM), PARAMETER :: iCFC113 =33
  INTEGER(KIND=JPIM), PARAMETER :: iCFC114 =34
  INTEGER(KIND=JPIM), PARAMETER :: iCFC115 =35
  INTEGER(KIND=JPIM), PARAMETER :: iCCL4   =36
  INTEGER(KIND=JPIM), PARAMETER :: iCLNO2  =37
  INTEGER(KIND=JPIM), PARAMETER :: iCH3CCL3=38
  INTEGER(KIND=JPIM), PARAMETER :: iCH3CL  =39
  INTEGER(KIND=JPIM), PARAMETER :: iHCFC22 =40
  INTEGER(KIND=JPIM), PARAMETER :: iCH3BR  =41
  INTEGER(KIND=JPIM), PARAMETER :: iHF     =42
  INTEGER(KIND=JPIM), PARAMETER :: iHA1301 =43
  INTEGER(KIND=JPIM), PARAMETER :: iHA1211 =44
  INTEGER(KIND=JPIM), PARAMETER :: iCHBR3  =45
  INTEGER(KIND=JPIM), PARAMETER :: iCLOO   =46
  INTEGER(KIND=JPIM), PARAMETER :: iO      =47
  INTEGER(KIND=JPIM), PARAMETER :: iO1D    =48
  INTEGER(KIND=JPIM), PARAMETER :: iN      =49
  INTEGER(KIND=JPIM), PARAMETER :: iCLO    =50
  INTEGER(KIND=JPIM), PARAMETER :: iCL     =51
  INTEGER(KIND=JPIM), PARAMETER :: iBR     =52
  INTEGER(KIND=JPIM), PARAMETER :: iBRO    =53
  INTEGER(KIND=JPIM), PARAMETER :: iH      =54
  INTEGER(KIND=JPIM), PARAMETER :: iH2     =55
  INTEGER(KIND=JPIM), PARAMETER :: iCO2    =56
  INTEGER(KIND=JPIM), PARAMETER :: iBR2    =57
  INTEGER(KIND=JPIM), PARAMETER :: ICH2BR2 =58
  INTEGER(KIND=JPIM), PARAMETER :: iStratAer=59
  INTEGER(KIND=JPIM), PARAMETER :: IOCS    =60
  INTEGER(KIND=JPIM), PARAMETER :: ISO2    =61
  INTEGER(KIND=JPIM), PARAMETER :: ISO3    =62
  INTEGER(KIND=JPIM), PARAMETER :: ISO4    =63
  INTEGER(KIND=JPIM), PARAMETER :: IH2SO4  =64

  ! heterogeneous reactions
  !
  INTEGER(KIND=JPIM), PARAMETER  :: NHET = 9

!-----------------------------------------------------------------------
!  Variables for Aerosols (STS) & PSC's (STS,SAT,NAT,ICE).
! NBINS  = Number of particle size bins - MUST be >=36 (tested at v3s10)
! NTYPES = Number of types of particles modelled in PSCBOX
! NWORK  = Dimension of work array
! NCLASS = Number of size classes in culumated size distributions
! rmin,rmax = min and max radius of particles (in microns)
! Notice HUGE array psss where size distributions are stored.
!-----------------------------------------------------------------------
  INTEGER(KIND=JPIM), PARAMETER :: nbins  = 36, &
     &                             ntypes = 4,  &
     &                             nwork  = 13, &
     &                             nclass = 13

  REAL(KIND=JPRB), PARAMETER :: rmin=0.002, rmax=36.0  ! microns

  REAL(KIND=JPRB), DIMENSION(NBINS,8) :: PTSIZE
! related to aerosol size distribution
  REAL(KIND=JPRB)                     :: D1

! -- Variable used in aero setup using values from 2D model of S. Bekki
!VH  INTEGER(KIND=JPIM), PARAMETER   :: NAER = 7
!VH Limit dimension to 2: Only SAD and NTOT are currently in use...
  INTEGER(KIND=JPIM), PARAMETER   :: NAER = 2

!-----------------------------------------------------------------------
!     indexes declaration for aerosols/psc fields
!-----------------------------------------------------------------------
  INTEGER(KIND=JPIM), PARAMETER :: iaer_sad   =1
  INTEGER(KIND=JPIM), PARAMETER :: iaer_ntot  =2
!VH  INTEGER(KIND=JPIM), PARAMETER :: iaer_h2so4 =3
!VH  INTEGER(KIND=JPIM), PARAMETER :: iaer_vsed  =4
!VH  INTEGER(KIND=JPIM), PARAMETER :: iaer_psc   =5
!VH  INTEGER(KIND=JPIM), PARAMETER :: iaer_nat   =6
!VH  INTEGER(KIND=JPIM), PARAMETER :: iaer_ice   =7

!-----------------------------------------------------------------------
!     aerosols/psc fields short names
!-----------------------------------------------------------------------
  CHARACTER(len=5), PARAMETER, DIMENSION(NAER) :: aer_name= &
     & (/'  sad' &
     & , ' ntot' /)
!VH     & , 'h2so4' &
!VH     & , ' vsed' &
!VH     & , '  psc' &
!VH     & , '  nat' &
!VH     & , '  ice' /)

!-----------------------------------------------------------------------
!     aerosols/psc fields units
!-----------------------------------------------------------------------
  CHARACTER(len=11), PARAMETER, DIMENSION(NAER) :: aer_units= &
     & (/'micron2/cm3' &
     &  ,' partic/cm3' /)
!VH     &  ,'        vmr' &
!VH     &  ,'        m/s' &
!VH     &  ,'       none' &
!VH     &  ,'      hours' &
!VH     &  ,'      hours' /)


!     Data requierd to calculate aerosol number densities consistent
!     with SAGE II averaged extinction measurements, assuming constant
!     lognormal size distribution median radius and geometric standard
!     deviation.
!     Data source: Hitchman et al.,  JGR 99, 20689, 1994.
!

INTEGER(KIND=JPIM), PARAMETER :: ILAT_SAGE=19,IALT_SAGE=14

REAL(KIND=JPRB), DIMENSION (ILAT_SAGE), PARAMETER :: LAT_SAGE = (/ &
     &-90.,-80.,-70.,-60.,-50.,-40.,-30.,-20.,-10., &
     & 0.,10.,20.,30.,40.,50.,60.,70.,80.,90. /)
REAL(KIND=JPRB), DIMENSION(IALT_SAGE), PARAMETER :: ALT_SAGE = (/ &
     & 8.,10.,12.,14.,16.,18.,20.,22.,24.,26.,28.,30.,32.,34./)
REAL(KIND=JPRB), DIMENSION(ILAT_SAGE), PARAMETER :: ZTROP_SAGE = (/ &
     &8.,8.,8.,8.,8.,8.,10.,16.,16., &
     &16.,16.,14.,10.,8.,8.,8.,8.,8.,8. /)
REAL(KIND=JPRB), DIMENSION(ILAT_SAGE,IALT_SAGE), PARAMETER :: E_SAGE= &
     & RESHAPE ((/ &
     &3.9,3.9,5.0,5.1,5.0,6.7,3.4,2.1,2.1, &
     &2.7,2.5,1.7,4.0,6.0,5.8,6.6,5.9,6.1,6.1, &
     &5.3,5.3,5.6,5.7,5.3,5.1,3.4,2.1,2.1, &
     &2.7,2.5,1.7,4.0,7.1,7.7,8.6,7.9,7.8,7.8, &
     &5.6,5.6,5.3,4.6,4.5,3.8,2.6,2.1,2.1, &
     &2.7,2.5,1.7,3.0,5.0,6.0,5.7,7.1,6.6,6.6, &
     &4.8,4.8,4.5,3.9,4.2,3.3,2.5,2.1,2.1, &
     &2.7,2.5,1.7,2.8,4.2,5.3,4.6,5.8,5.4,5.4, &
     &2.1,2.1,2.8,3.1,4.2,3.8,2.8,2.1,2.1, &
     &2.7,2.5,2.3,3.0,4.4,4.8,4.1,4.4,4.1,4.1, &
     &1.0,1.0,1.7,2.0,3.3,3.7,3.9,3.3,3.2, &
     &2.9,2.8,3.2,4.1,4.1,3.6,3.1,2.6,2.4,2.4, &
     &.62,.62,.91,1.1,2.0,2.5,3.1,4.0,5.5, &
     &8.1,6.7,4.4,3.2,2.6,2.1,1.7,1.6,1.2,1.2, &
     &.31,.31,.47,.56,1.0,1.3,1.7,2.6,4.3, &
     &5.5,4.7,3.0,1.8,1.3,1.0,.76,.70,.54,.54, &
     &.15,.15,.25,.28,.48,.66,.87,1.5,2.6, &
     &3.1,2.5,1.5,.88,.60,.47,.33,.34,.28,.28, &
     &.08,.08,.13,.14,.20,.29,.42,.74,1.3, &
     &1.5,1.3,.80,.42,.26,.21,.13,.16,.16,.16, &
     &.05,.05,.08,.07,.08,.12,.18,.34,.56, &
     &.66,.57,.37,.18,.11,.09,.06,.08,.10,.10, &
     &.05,.05,.05,.04,.04,.05,.08,.14,.22, &
     &.24,.20,.14,.08,.05,.05,.03,.05,.06,.06, &
     &.01,.01,.02,.02,.03,.03,.04,.06,.09, &
     &.09,.07,.06,.04,.03,.03,.02,.02,.02,.02, &
     &.01,.01,.01,.02,.02,.02,.02,.03,.04, &
     &.04,.03,.03,.02,.02,.02,.02,.01,.01,.01 /),(/ILAT_SAGE,IALT_SAGE/))

REAL(KIND=JPRB), DIMENSION(ILAT_SAGE,IALT_SAGE) :: AREA_LAT
REAL(KIND=JPRB), DIMENSION(ILAT_SAGE,IALT_SAGE) :: Y2AREA_LAT

LOGICAL    :: LREAD_AEROCLIM = .TRUE. 


!----------------------------------------------------------------------------------------------
! J. Debosscher: declarations related to Aerosol SAD climatology
INTEGER(KIND=JPIM)           :: NMONTH_CLIM, NLAT_CLIM, NLEV_CLIM, MONTHNUM
REAL(KIND=JPRB), ALLOCATABLE :: SAD_CLIM(:,:,:),P_CLIM(:,:,:),LAT_CLIM(:)
!----------------------------------------------------------------------------------------------


! Tropopause climatology from SOCRATES
INTEGER(KIND=JPIM),PARAMETER                    :: NLAT_PTROPO=180
REAL(KIND=JPRB), DIMENSION(NLAT_PTROPO)         :: PTROP_SOC_B, LATS_PTROPO

  ! Boundary conditions at surface for stratospheric species
  ! Which currently don't feature emissions (!) 
  ! Just a simple parameter, taken from BASCOE field for 1 April 2008
  ! Units in [mole/mole]
  ! Note that we also introduce BC for O3, to prevent a blow-up in the troposphere
  ! also include CH4.
  INTEGER(KIND=JPIM), PARAMETER  :: NBC = 21
  INTEGER(KIND=JPIM),DIMENSION(NBC), PARAMETER ::BASCOE_BC=(/ &
  & IN2O  , ICFC11  , ICFC12 , ICFC113, ICFC114, &
  & ICFC115, ICH2BR2, ICCL4 , ICH3CCL3, IHCFC22, &
  & IHA1301, IHA1211, ICH3Br, ICHBR3  , ICH3CL , &
  & ICO2   , IO3    , ICH4  , IBRO    , IHBR   , &
  & IOCS /)
  CHARACTER (LEN = 12), DIMENSION(NBC), PARAMETER :: BASCOE_BCNAME = &
     &      (/"n2o         " , &
     &        "cfc11       " , &
     &        "cfc12       " , &
     &        "cfc113      " , &
     &        "cfc114      " , &
     &        "cfc115      " , &
     &        "ch2br2      " , &
     &        "ccl4        " , &
     &        "ch3ccl3     " , &
     &        "hcfc22      " , &
     &        "ha1301      " , &
     &        "ha1211      " , &
     &        "ch3br       " , &
     &        "chbr3       " , &
     &        "ch3cl       " , &
     &        "co2         " , &
     &        "o3          " , &
     &        "ch4         " , &
     &        "bro         " , &
     &        "hbr         " , &
     &        "ocs         " /)

  REAL(KIND=JPRB) , DIMENSION(NBC), PARAMETER ::BASCOE_BCVAL=(/ &
  &  3.27E-7 , 2.33E-10, 5.20E-10, 7.27E-11, 1.63E-11, &
  &  8.43E-12, 1.20E-12, 83.2E-12, 3.70E-12, 2.29E-10, &
  &  3.30E-12, 3.75E-12, 6.65E-12, 1.17E-12, 5.34E-10, &
  &  400.E-6 , 3.00E-8 , 1.74E-6 , 0.10E-12, 0.10E-12, &
  &  266.E-12 /)
  ! *WAERNING* o3 here is arbitrary !
  REAL(KIND=JPRB) , DIMENSION(NBC), PARAMETER ::BASCOE_BCVAL_1750=(/ &
  &  2.74E-7 , 0.00E-10, 0.00E-10, 0.00E-11, 0.00E-12, &
  &  0.00E-12, 1.20E-12, 2.94E-14, 0.00E-11, 0.00E-10, &
  &  0.00E-12, 0.00E-12, 5.30E-12, 1.20E-12, 4.57E-10, &
  &  277.E-6 , 3.00E-8 , 0.731E-6, 0.00E-12, 0.00E-12, &
  &  266.E-12 /)

  REAL(KIND=JPRB) , DIMENSION(NBC), PARAMETER ::BASCOE_BCVAL_1850=(/ &
  &  2.73E-7 , 0.00E-10, 0.00E-10, 0.00E-11, 0.00E-12, &
  &  0.00E-10, 1.20E-12, 2.94E-14, 0.00E-11, 0.00E-10, &
  &  0.00E-12, 0.00E-12, 5.30E-12, 1.20E-12, 4.57E-10, &
  &  284.E-6 , 3.00E-8 , 0.808E-6, 0.00E-12, 0.00E-12, &
  &  266.E-12 /)


END MODULE BASCOE_MODULE
