MODULE YOS_URB

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!------------------------------------------------------------------------------

!     PURPOSE
!     -------
! Declare parameters used for urban tile

!     INTERFACE.
!     ----------
!     CALLLED FROM *SUSURF*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     Urban scheme based on Macdonald 1998, Harman 2004, Oleson 2008 and Porson 2010

!     MODIFICATIONS
!     -------------
!==============================================================================

! Modules used : 
USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

TYPE :: TURB

REAL (KIND = JPRB):: RBUIZ0M        ! Roughness length for momentum for building material 
REAL (KIND = JPRB):: RWALTHK        ! Thickness of wall (m)  
REAL (KIND = JPRB):: RROOTHK        ! Thickness of roof (m)
!REAL (KIND = JPRB):: RFLOTHK        ! Thickness of floor (m)
REAL (KIND = JPRB):: RROATHK        ! Thickness of road (m)
REAL (KIND = JPRB):: RWALALB        ! Wall albedo
REAL (KIND = JPRB):: RROOALB        ! Roof albedo
REAL (KIND = JPRB):: RROAALB        ! Road albedo
REAL (KIND = JPRB):: RWALEMIS       ! Wall emissivity
REAL (KIND = JPRB):: RROOEMIS       ! Roof emissivity
REAL (KIND = JPRB):: RROAEMIS       ! Road emissivity
REAL (KIND = JPRB):: RSOIEMIS       ! Soil emissivity !NOTE THIS CAN BE TAKEN FROM MODEL
!REAL (KIND = JPRB):: RFLOEMIS       ! Floor emissivity
REAL (KIND = JPRB):: RWALVHC        ! Volumetric heat capactiy of wall (K m-3K-1)
REAL (KIND = JPRB):: RROOVHC        ! Volumetric heat capactiy of roof (K m-3K-1)
REAL (KIND = JPRB):: RROAVHC        ! Volumetric heat capactiy of road (K m-3K-1)
REAL (KIND = JPRB):: RCANVHC        ! Volumetric heat capactiy of canyon (K m-3K-1)
!REAL (KIND = JPRB):: RFLOVHC        ! Volumetric heat capactiy of floor (K m-3K-1)
REAL (KIND = JPRB):: RAIRVHC        ! Volumetric heat capactiy of air (K m-3K-1)
REAL (KIND = JPRB):: RSOITMP        ! Sub-urban soil temperature (K)
REAL (KIND = JPRB):: RSOITC         ! Soil thermal conductivity (W m-1K-1)   !NOTE THIS CAN BE TAKEN FROM MODEL
REAL (KIND = JPRB):: RWALTC         ! Wall thermal conductivity (W m-2K-1)
REAL (KIND = JPRB):: RROOTC         ! Roof thermal conductivity (W m-2K-1)
REAL (KIND = JPRB):: RCANTC         ! Canyon thermal conductivity (W m-2K-1)
REAL (KIND = JPRB):: RROATC         ! Road thermal conductivity (W m-2K-1)
REAL (KIND = JPRB):: RSATSH         ! Saturated specific humidity (kg kg-1)

REAL (KIND = JPRB):: RCDG           ! Drag Coefficient (Macdonald 1998)
REAL (KIND = JPRB):: RCDA           ! Coefficient For Canyon (A) (Macdonald 1998)
REAL (KIND = JPRB):: RCHIS          ! Fraction of SW radiation scatter by sky
REAL (KIND = JPRB):: RSBCONS        ! Stefan-Boltzmann Constant (W m-2K-4)
REAL (KIND = JPRB):: RGACC          ! Acceleration due to gravity (m s-2)
REAL (KIND = JPRB):: RAIRRHO        ! Density of air (Kg m-3)    !NOTE THIS CAN BE TAKEN FROM MODEL
REAL (KIND = JPRB):: RMODZ          ! Lowest model height (m)    !NOTE THIS CAN BE TAKEN FROM MODEL
REAL (KIND = JPRB):: RVKSQ          ! Von Karman constant squared
REAL (KIND = JPRB):: REXPDR         ! Exponential decay of recirculation

REAL (KIND = JPRB):: RHWR           ! Height-to-roadwidth (aspect) ratio
REAL (KIND = JPRB):: RHGT           ! Average building height (m)
REAL (KIND = JPRB):: RWRR           ! Road-to-building+road ratio (width)
REAL (KIND = JPRB):: RCANALB        ! Canyon Albedo
REAL (KIND = JPRB):: RCANEMIS       ! Canyon Emissivity
REAL (KIND = JPRB):: RURBEMIS       ! Urban Emissivity
REAL (KIND = JPRB):: RCANZTM        ! Canyon Roughness Length for Momentum
REAL (KIND = JPRB):: RCANZTH        ! Canyon Roughness Length for Heat
REAL (KIND = JPRB):: RCANHC         ! Canyon Heat Coefficient
REAL (KIND = JPRB):: RROOZTM        ! Roof Roughness Length for Momentum
REAL (KIND = JPRB):: RROOZTH        ! Roof Roughness Length for Heat
REAL (KIND = JPRB):: RROOHC         ! Roof Heat Coefficient
REAL (KIND = JPRB):: RCANRES        ! Canyon Resistance
REAL (KIND = JPRB):: RROORES        ! Roof Resistance

! AVERAGED URBAN VALUES

REAL (KIND = JPRB):: RURBZTM        ! Urban Roughness Length for Momentum
REAL (KIND = JPRB):: RURBZTH        ! Urban Roughness Length for Heat
REAL (KIND = JPRB):: RURBHC         ! Urban Heat Coefficient
REAL (KIND = JPRB):: RURBRES        ! Urban Resistance
REAL (KIND = JPRB):: RURBTC         ! Urban Thermal Conductivity
REAL (KIND = JPRB):: RURBVHC        ! Urban Volumetric Heat Capacity
REAL (KIND = JPRB):: RURBTC1        ! Urban Thermal Conductivity (level-1)

! MOISTURE VARIABLES

REAL (KIND = JPRB):: RURBALP        ! URBAN ALPHA-FACTOR IN VAN GENUCHTEN
REAL (KIND = JPRB):: RURBCON        ! URBAN WATER CONDUCTIVITY AT SATURATION
REAL (KIND = JPRB):: RURBLAM        ! URBAN LAMBDA-FACTOR IN VAN GENUCHTEN
REAL (KIND = JPRB):: RURBSAT        ! URBAN SOIL WATER CONTENT AT SATURATION
REAL (KIND = JPRB):: RURBSRES       ! URBAN SOIL WATER CONTENT RESIDUAL



LOGICAL           :: LURBAN         ! True when using urban scheme for energy balance (not used)
LOGICAL           :: LURBUI         ! True when using building heating scheme for energy balance (not used)

END TYPE TURB

!==============================================================================

END MODULE YOS_URB
