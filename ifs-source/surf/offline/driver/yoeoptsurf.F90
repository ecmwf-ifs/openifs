MODULE YOEOPTSURF
! (C) Copyright 2021- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

USE PARKIND1  ,ONLY : JPRB, JPIM

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOEOPTSURF* - OPTIMISED PARAMETERS RELATED TO LAND SURFACE
!     -----------------------------------------------------------------

!        * SURFACE PARAMETERS USED IN OSM OPTIMIZATION  *


! Include vegetationn surface parameters that are optimized with observations (passed via namelist to surface model)

REAL(KIND=JPRB), ALLOCATABLE :: RVR0VT(:,:) ! reference ecosystem respiration [Kg CO2 m-2 s-1]
REAL(KIND=JPRB), ALLOCATABLE :: RVCMAX25(:,:) ! Maximum rate of Rubisco activity-limited carboxylation at 25°C (\mu mol.m^{-2}.s^{-1})
REAL(KIND=JPRB), ALLOCATABLE :: RHUMREL(:,:) !scaling factor for soil moisture stress (optimized per PFT with FLUXNET data)
REAL(KIND=JPRB), ALLOCATABLE :: RA1(:,:)   ! Empirical factor involved in the calculation of fvpd (-)
                                       ! See Table 2 of Yin et al. (2009)
REAL(KIND=JPRB), ALLOCATABLE :: RB1(:,:)   ! Empirical factor involved in the calculation of fvpd (-)
                                       ! See Table 2 of Yin et al. (2009)
REAL(KIND=JPRB), ALLOCATABLE :: RG0(:,:)  ! Residual stomatal conductance when irradiance approaches zero (mol CO2 m−2 s−1 bar−1)
REAL(KIND=JPRB), ALLOCATABLE :: RGM25(:,:) ! Mesophyll diffusion conductance at 25Â°C (mol m-2 s-1 bar-1) 
REAL(KIND=JPRB), ALLOCATABLE :: RE_VCMAX(:,:) ! Energy of activation for Vcmax (J mol-1)
                                          ! See Table 2 of Yin et al. (2009) for C4 plants
                                          ! and Kattge & Knorr (2007) for C3 plants (table 3)
REAL(KIND=JPRB), ALLOCATABLE :: RE_JMAX(:,:) ! Energy of activation for Jmax (J mol-1)
                                         ! See Table 2 of Yin et al. (2009) for C4 plants
                                         ! and Kattge & Knorr (2007) for C3 plants (table 3)

! Optimized parameters read from namelist (1D version)
REAL(KIND=JPRB), ALLOCATABLE :: OVR0VT(:) ! reference ecosystem respiration [Kg CO2 m-2 s-1]
REAL(KIND=JPRB), ALLOCATABLE :: OVCMAX25(:) ! Maximum rate of Rubisco activity-limited carboxylation at 25°C (\mu mol.m^{-2}.s^{-1})
REAL(KIND=JPRB), ALLOCATABLE :: OHUMREL(:) !scaling factor for soil moisture stress (optimized per PFT with FLUXNET data)
REAL(KIND=JPRB), ALLOCATABLE :: OA1(:)   ! Empirical factor involved in the calculation of fvpd (-)
                                       ! See Table 2 of Yin et al. (2009)
REAL(KIND=JPRB), ALLOCATABLE :: OB1(:)   ! Empirical factor involved in the calculation of fvpd (-)
                                       ! See Table 2 of Yin et al. (2009)
REAL(KIND=JPRB), ALLOCATABLE :: OG0(:)  ! Residual stomatal conductance when irradiance approaches zero (mol CO2 m−2 s−1 bar−1)
REAL(KIND=JPRB), ALLOCATABLE :: OGM25(:) ! Mesophyll diffusion conductance at 25Â°C (mol m-2 s-1 bar-1) 
REAL(KIND=JPRB), ALLOCATABLE :: OE_VCMAX(:) ! Energy of activation for Vcmax (J mol-1)
                                          ! See Table 2 of Yin et al. (2009) for C4 plants
                                          ! and Kattge & Knorr (2007) for C3 plants (table 3)
REAL(KIND=JPRB), ALLOCATABLE :: OE_JMAX(:) ! Energy of activation for Jmax (J mol-1)
                                         ! See Table 2 of Yin et al. (2009) for C4 plants
                                         ! and Kattge & Knorr (2007) for C3 plants (table 3)

!


!
!     REFERENCE.
!     ----------

!     Anna Agusti-Panareda       E.C.M.W.F.      02/06/2021

!     MODIFICATIONS
!     -------------

!     
!     ------------------------------------------------------------------

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! Vegetation model parameters optimized with observations
! RVR0VT: REAL PFT ARRAY : reference ecosystem respiration [Kg CO2 m-2 s-1]! 
! RVCMAX25 : REAL PFT ARRAY: MAX RATE of Rubisco activity-limited carboxylation at 25°C (\mu mol.m^{-2}.s^{-1})
! RHUMREL : REAL PFT ARRAY: SCALING FACTOR for soil moisture stress (optimized per PFT with FLUXNET data)
! RA1: REAL ARRAY (PHOTOSYNTHESIS PATHWAY): Empirical factor in the calculation of fvpd (-)
! RB1: REAL ARRAY (PHOTOSYNTHESIS PATHWAY): Empirical factor in the calculation of fvpd (-)
! RG0: REAL ARRAY (PHOTOSYNTHESIS PATHWAY): Residual stomatal conductance when irradiance approaches zero (mol CO2 m−2 s−1 bar−1)
! RGM25 REAL ARRAY (PHOTOSYNTHESIS PATHWAY): Mesophyll diffusion conductance at 25oC (mol m-2 s-1 bar-1) 
! RE_VCMAX: REAL ARRAY (PHOTOSYNTHESIS PATHWAY): Energy of activation for Vcmax (J mol-1)
! RE_JMAX: REAL ARRAY (PHOTOSYNTHESIS PATHWAY): Energy of activation for Jmax (J mol-1)
!     -----------------------------------------------------------------
END MODULE YOEOPTSURF
