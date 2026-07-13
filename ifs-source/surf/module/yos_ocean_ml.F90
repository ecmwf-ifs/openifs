MODULE YOS_OCEAN_ML

! YOS_OCEAN_ML : Data Module for the Ocean Mixed Layer Model.

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!---------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     , JPRB

IMPLICIT NONE

SAVE

INTEGER(KIND=JPIM),PARAMETER :: NSCLRO = 2 ! Number of scalars, 2 for T & S 
INTEGER(KIND=JPIM),PARAMETER :: NVELO = 2  ! Number of velocity, 2 for U & V

TYPE :: TOCEAN_ML
INTEGER(KIND=JPIM) :: NDIFFO
INTEGER(KIND=JPIM) :: NZTBLO          ! Dimension of M-O table
INTEGER(KIND=JPIM) :: NUTBLO          ! ...
INTEGER(KIND=JPIM) :: NOCNSTEP_KPP    ! No. of ocean time step
INTEGER(KIND=JPIM) :: NJERLOV_KPP     ! Jerlov water type
INTEGER(KIND=JPIM) :: NITERMAX_KPP    ! Max. num. of iteration 
                                      ! of Euler-backward scheme  

LOGICAL :: LEOCML               ! Ocean mixed layer switch on/off
LOGICAL :: LOCDEPTH_KPP         ! .TRUE. : ocean topo.,  .FALSE. : fixed depth
LOGICAL :: LINIT_WSCALE         ! initialize universal function table
LOGICAL :: LCDIAG_KPP           ! .TRUE. : compute diagnostic variables
LOGICAL :: LVL_NEW_KPP          ! new vertical level definition
LOGICAL :: LGEO_KPP             ! geostrophic flow decomposition 

REAL(KIND=JPRB) :: RDEPTH_KPP   ! used for fixed ocean depth (>0) [m]
REAL(KIND=JPRB) :: RDSCALE_KPP  ! scale parameter, = -lambda in LMD94 Appendix D.
                                ! =  0.0 : constant
                                ! > 0.0 more layer concentrated near the surface

! new vertical level definition
! PDO = ZA*ZX1 + ZB*(ZX2**RVC_KPP) + RVD_KPP*(ZX2**RVE_KPP)  
REAL(KIND=JPRB) :: RVC_KPP      
REAL(KIND=JPRB) :: RVD_KPP    
REAL(KIND=JPRB) :: RVE_KPP
REAL(KIND=JPRB) :: RVDEL_KPP

! bldepth
REAL(KIND=JPRB) :: REPS_KPP     ! small number to avoid divided by 0
REAL(KIND=JPRB) :: RICR_KPP     ! Critical Richardson Number for BL depth
REAL(KIND=JPRB) :: REPSILON_KPP ! epsilon: lower limit of sigma
REAL(KIND=JPRB) :: RCEKMAN_KPP  ! Constans for estimating Ekman depth
REAL(KIND=JPRB) :: RCMONOB_KPP  ! Tuning parameter(?) for Ekman depth
REAL(KIND=JPRB) :: RCV_KPP      ! constants for vertical turbulent shear
REAL(KIND=JPRB) :: RVONK        ! von Kalman constant
REAL(KIND=JPRB) :: RVTC         

! blmix
REAL(KIND=JPRB) :: RCSTAR_KPP   ! constant for convective velocity scale

! ddmix
REAL(KIND=JPRB) :: RRRHO0       ! Rp=(alpha*delT)/(beta*delS)
REAL(KIND=JPRB) :: RRRHO0_D06   ! Rp  Danabasoglu et al. (2006)
REAL(KIND=JPRB) :: RDSFMAX      ! salt fingering diffusivity in eq.(31) LMD94
REAL(KIND=JPRB) :: RDMOL        ! molecular viscosity

! ri_iwmix
INTEGER(KIND=JPIM) :: NMRI_KPP  ! number of vertical smoothing passes
REAL(KIND=JPRB) :: RRIINFTY     ! critical Ri number for shear instability  
REAL(KIND=JPRB) :: RDIFM_MAX    ! max visc due to shear instability  (m^2/s)
REAL(KIND=JPRB) :: RDIFM_IW     ! background/internal waves visc(m^2/s)
REAL(KIND=JPRB) :: RDIFS_MAX    ! max diff due to shear instability  (m^2/s)
REAL(KIND=JPRB) :: RDIFS_IW     ! background/internal waves diff(m^2/s)

! swfrac_opt
REAL(KIND=JPRB) :: RFAC_JERLOV(5)
REAL(KIND=JPRB) :: RA1_JERLOV(5)
REAL(KIND=JPRB) :: RA2_JERLOV(5)
REAL(KIND=JPRB) :: RMIN_JERLOV

! tridrhs
INTEGER(KIND=JPIM) :: NPD_KPP

! vmix 
LOGICAL :: LDD_KPP
LOGICAL :: LSF_NEW_KPP
LOGICAL :: LKPP_KPP
LOGICAL :: LNBFLX_KPP
LOGICAL :: LRI_KPP

REAL(KIND=JPRB) :: RDIFM_BOT    ! Diffusivity of the bottom boundary
REAL(KIND=JPRB) :: RDIFS_BOT    ! ...
REAL(KIND=JPRB) :: RSICE        ! Salinity of Sea Ice
REAL(KIND=JPRB) :: RTICE        ! Temperature below Sea Ice
REAL(KIND=JPRB) :: RRAMSICE     ! Parameter for sea-ice flux

! wscale
REAL(KIND=JPRB),ALLOCATABLE :: RWMT(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: RWST(:,:)
REAL(KIND=JPRB) :: RZETMIN_KPP
REAL(KIND=JPRB) :: RZETMAX_KPP
REAL(KIND=JPRB) :: RUMIN_KPP
REAL(KIND=JPRB) :: RUMAX_KPP

REAL(KIND=JPRB) :: RC1_KPP  ! constant for universal function
REAL(KIND=JPRB) :: RC2_KPP  ! ...
REAL(KIND=JPRB) :: RC3_KPP  ! ...
REAL(KIND=JPRB) :: RZETAM   ! ...
REAL(KIND=JPRB) :: RZETAS   ! ...
REAL(KIND=JPRB) :: RAM_KPP  ! ...
REAL(KIND=JPRB) :: RAS_KPP  ! ...
REAL(KIND=JPRB) :: RCM_KPP  ! ...
REAL(KIND=JPRB) :: RCS_KPP  ! constants for vertical turbulent shear

! semi-implicit time integration
REAL(KIND=JPRB) :: RLAMBDA_KPP
REAL(KIND=JPRB) :: RTOLFAC_KPP
END TYPE TOCEAN_ML

!---------------------------------------------------------------------

END MODULE YOS_OCEAN_ML
