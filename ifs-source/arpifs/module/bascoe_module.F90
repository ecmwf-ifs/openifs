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


! parameters that define the NAT/ICE PSC parameterization in BASCOE-based stratospheric chemistry:
     !NAT sedimentation lifetime (days): 
   REAL(KIND=JPRB), PARAMETER :: RCHEM_SEDIM_NAT=20.
     !Surface Area Density of NAT PSC as soon as they appear
   REAL(KIND=JPRB), PARAMETER :: RCHEM_SAD_NAT_PSC=2.E-7
     !Supersaturation ratio for NAT PSC
   REAL(KIND=JPRB), PARAMETER :: RCHEM_SSR_NAT_PSC=10.

     !ICE sedimentation lifetime (days): 
   REAL(KIND=JPRB), PARAMETER :: RCHEM_SEDIM_ICE=9.
     !Surface Area Density of ICE PSC as soon as they appear
   REAL(KIND=JPRB), PARAMETER :: RCHEM_SAD_ICE_PSC=2.E-6
     !Supersaturation ratio for ICE PSC
   REAL(KIND=JPRB), PARAMETER :: RCHEM_SSR_ICE_PSC=1.

!----------------------------------------------------------------------------------------------
! J. Debosscher: declarations related to Aerosol SAD climatology
INTEGER(KIND=JPIM)           :: NSYEAR_CLIM, NEYEAR_CLIM, NYEAR_CLIM, NMONTH_CLIM, NLAT_CLIM, NLEV_CLIM
REAL(KIND=JPRB), ALLOCATABLE :: SAD_CLIM(:,:,:,:),P_CLIM(:,:,:,:),LAT_CLIM(:)
!----------------------------------------------------------------------------------------------


! Tropopause climatology from SOCRATES
INTEGER(KIND=JPIM),PARAMETER                    :: NLAT_PTROPO=180
REAL(KIND=JPRB), DIMENSION(NLAT_PTROPO)         :: PTROP_SOC_B, LATS_PTROPO



END MODULE BASCOE_MODULE
