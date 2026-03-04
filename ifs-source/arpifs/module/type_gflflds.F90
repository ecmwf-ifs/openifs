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

MODULE TYPE_GFLFLDS

USE PARKIND1  ,ONLY : JPIM
USE YOM_YGFL  ,ONLY : JPAERO,JPGHG,JPCHEM,JPLIMA,JPNOGW,JPUVP,JPAEROUT,JPERA40, &
 & JPGFL,JPNAMED_GFL,JPFORC,JPEZDIAG,JPSLDIA,JPCH4S,JPAEROCLIM,JPPHYCTY,JPFSD,&
 & JPEDRP
USE YOMCOMPO    , ONLY : NPEMIS3D
USE YOE_AERODIAG, ONLY : NPAERAOT, NPAERLISI

IMPLICIT NONE
SAVE

!-------------------------------------------------------------------------
! Derived types for describing the list of variables in the GFL structure:
! The descriptors themselves
! (YGFL and YCOMP) can be found in module yom_ygfl.F90.
!-------------------------------------------------------------------------
! Modifications:
!   F. Vana & M. Kharoutdinov 06-Feb-2015  Super-parametrization scheme.
!   2017-11-11 M Ahlgrimm - add cloud heterogeneity FSD

! ky: to be moved later in YOM_YGFL

INTEGER(KIND=JPIM),PARAMETER :: JPEXT=JPGFL-JPNAMED_GFL-JPGHG-JPCHEM-JPLIMA-JPFORC- &
 & JPEZDIAG-JPAERO-JPERA40-JPNOGW-JPSLDIA-JPCH4S-NPAERAOT-NPAERLISI-&
 & JPAEROUT-JPUVP-JPAEROCLIM-JPPHYCTY-JPFSD-JPEDRP

TYPE TYPE_IGFLFLD
! List of GFL variables without derivatives for integer variables.
! * historical (advectable) variables:
INTEGER(KIND=JPIM) :: IQ        ! Specific humidity
INTEGER(KIND=JPIM) :: IL        ! Liquid water
INTEGER(KIND=JPIM) :: ILCONV    ! Liquid water (CONV. PART)
INTEGER(KIND=JPIM) :: IICONV    ! Ice    water (CONV. PART)
INTEGER(KIND=JPIM) :: IRCONV    ! Rain         (CONV. PART)
INTEGER(KIND=JPIM) :: ISCONV    ! Snow         (CONV. PART)
INTEGER(KIND=JPIM) :: ILRAD     ! Total liquid water for radiation
INTEGER(KIND=JPIM) :: II        ! Ice
INTEGER(KIND=JPIM) :: IIRAD     ! Total ice for radiation
INTEGER(KIND=JPIM) :: IA        ! Cloud fraction
INTEGER(KIND=JPIM) :: IO3       ! Ozone
INTEGER(KIND=JPIM) :: IS        ! Snow
INTEGER(KIND=JPIM) :: IRR       ! Rain
INTEGER(KIND=JPIM) :: IG        ! Graupels
INTEGER(KIND=JPIM) :: IH        ! Hail
INTEGER(KIND=JPIM) :: ITKE      ! TKE (total kinetic entrainment)
INTEGER(KIND=JPIM) :: IAERO     ! Aerosols
INTEGER(KIND=JPIM) :: IEXT      ! Extra-GFL

INTEGER(KIND=JPIM) :: ILIMA     ! LIMA group
! * other variables:
INTEGER(KIND=JPIM) :: IDAL      ! downdraught mesh fraction (for ALARO physics)
INTEGER(KIND=JPIM) :: IDOM      ! downdraught vertical velocity (for ALARO physics)
INTEGER(KIND=JPIM) :: IUAL      ! updraught mesh fraction (for ALARO physics)
INTEGER(KIND=JPIM) :: IUOM      ! updraught vertical velocity (for ALARO physics)
INTEGER(KIND=JPIM) :: IUEN      ! pseudo-historic entrainment (for ALARO physics)
INTEGER(KIND=JPIM) :: IUNEBH    ! pseudo-historic convective cloudiness (for ALARO physics)
INTEGER(KIND=JPIM) :: ITTE      ! total turbulent energy
INTEGER(KIND=JPIM) :: IMXL      ! prognostic mixing lentgh
INTEGER(KIND=JPIM) :: ISHTUR    ! shear source term for turbulence            
INTEGER(KIND=JPIM) :: IFQTUR    ! flux form source term for turbulence -moisture
INTEGER(KIND=JPIM) :: IFSTUR    ! flux form source term for turbulence -enthalpy
INTEGER(KIND=JPIM) :: IRKTH     ! Rasch-Kristjansson enthalpy tendency.
INTEGER(KIND=JPIM) :: IRKTQV    ! Rasch-Kristjansson water vapour tendency.
INTEGER(KIND=JPIM) :: IRKTQC    ! Rasch-Kristjansson condensates tendency.
INTEGER(KIND=JPIM) :: ISRC      ! Second-order flux (for AROME)
INTEGER(KIND=JPIM) :: ICPF      ! Convective precipitation flux
INTEGER(KIND=JPIM) :: ISPF      ! Stratiform precipitation flux
INTEGER(KIND=JPIM) :: ICVGQ     ! Moisture convergence (for MF physics)
INTEGER(KIND=JPIM) :: ISDSAT    ! Standard deviation of the saturation depression
INTEGER(KIND=JPIM) :: ICVV      ! Convective vertical velocity
INTEGER(KIND=JPIM) :: IQVA      ! Total humidity variation (for HIRLAM physics)
INTEGER(KIND=JPIM) :: ILRCH4    ! CH4 (methane) loss rate
INTEGER(KIND=JPIM) :: IEMIS3D   ! 3D emissions for atmospheric composition
INTEGER(KIND=JPIM) :: IFORC     ! Forcings (1D model)
INTEGER(KIND=JPIM) :: IEZDIAG   ! Easy diagnostics (for AROME physics)
INTEGER(KIND=JPIM) :: IGHG      ! Greenhouse gases (ECMWF)
INTEGER(KIND=JPIM) :: ICHEM     ! Chemical species (ECMWF)
INTEGER(KIND=JPIM) :: IERA40    ! ERA40 reanalysis fields
INTEGER(KIND=JPIM) :: INOGW     ! Diagnostic fields for NORO GWD scheme (ECMWF)
INTEGER(KIND=JPIM) :: IAERAOT   ! Aerosol layer optical thicknesses
INTEGER(KIND=JPIM) :: IAERLISI  ! Aerosol lidar simulator (extinction and backscatter)
INTEGER(KIND=JPIM) :: IAEROUT   ! Output aerosols for diagnostics
INTEGER(KIND=JPIM) :: IAERCLIM ! 3D aerosol input
INTEGER(KIND=JPIM) :: IUVP      ! Output fields for UV processor
INTEGER(KIND=JPIM) :: ISLDIA    ! Semi-Lagrangian dynamics diagnostic fields
INTEGER(KIND=JPIM) :: IPHYCTY   ! Specific mass change from physics (s-1) for CTY equation
INTEGER(KIND=JPIM) :: IFSD      ! cloud heterogenity FSD
INTEGER(KIND=JPIM) :: IEDRP     ! turbulence EDR Parameter diagnostics
END TYPE TYPE_IGFLFLD

TYPE TYPE_IGFLFLDD
! List of GFL variables with derivatives for integer variables.
! * historical (advectable) variables:
INTEGER(KIND=JPIM) :: IQ        ! Specific humidity
INTEGER(KIND=JPIM) :: IQL       !
INTEGER(KIND=JPIM) :: IQM       !
INTEGER(KIND=JPIM) :: IL        ! Liquid water
INTEGER(KIND=JPIM) :: ILL       !
INTEGER(KIND=JPIM) :: ILM       !
INTEGER(KIND=JPIM) :: ILCONV    ! Liquid water (CONV. PART)
INTEGER(KIND=JPIM) :: ILCONVL   !
INTEGER(KIND=JPIM) :: ILCONVM   !
INTEGER(KIND=JPIM) :: IICONV    ! Ice    water (CONV. PART)
INTEGER(KIND=JPIM) :: IICONVL   !
INTEGER(KIND=JPIM) :: IICONVM   !
INTEGER(KIND=JPIM) :: IRCONV    ! Rain         (CONV. PART)
INTEGER(KIND=JPIM) :: IRCONVL   !
INTEGER(KIND=JPIM) :: IRCONVM   !
INTEGER(KIND=JPIM) :: ISCONV    ! Snow         (CONV. PART)
INTEGER(KIND=JPIM) :: ISCONVL   !
INTEGER(KIND=JPIM) :: ISCONVM   !
INTEGER(KIND=JPIM) :: ILRAD     ! Total liquid water for radiation
INTEGER(KIND=JPIM) :: ILRADL    !
INTEGER(KIND=JPIM) :: ILRADM    !
INTEGER(KIND=JPIM) :: II        ! Ice
INTEGER(KIND=JPIM) :: IIL       !
INTEGER(KIND=JPIM) :: IIM       !
INTEGER(KIND=JPIM) :: IIRAD     ! Total ice for radiation
INTEGER(KIND=JPIM) :: IIRADL    !
INTEGER(KIND=JPIM) :: IIRADM    !
INTEGER(KIND=JPIM) :: IA        ! Cloud fraction
INTEGER(KIND=JPIM) :: IAL       !
INTEGER(KIND=JPIM) :: IAM       !
INTEGER(KIND=JPIM) :: IO3       ! Ozone
INTEGER(KIND=JPIM) :: IO3L      !
INTEGER(KIND=JPIM) :: IO3M      !
INTEGER(KIND=JPIM) :: IS        ! Snow
INTEGER(KIND=JPIM) :: ISL       !
INTEGER(KIND=JPIM) :: ISM       !
INTEGER(KIND=JPIM) :: IRR       ! Rain
INTEGER(KIND=JPIM) :: IRRL      !
INTEGER(KIND=JPIM) :: IRRM      !
INTEGER(KIND=JPIM) :: IG        ! Graupels
INTEGER(KIND=JPIM) :: IGL       !
INTEGER(KIND=JPIM) :: IGM       !
INTEGER(KIND=JPIM) :: IH        ! Hail
INTEGER(KIND=JPIM) :: IHL       !
INTEGER(KIND=JPIM) :: IHM       !
INTEGER(KIND=JPIM) :: ITKE      ! TKE (total kinetic entrainment)
INTEGER(KIND=JPIM) :: ITKEL     ! 
INTEGER(KIND=JPIM) :: ITKEM     !
INTEGER(KIND=JPIM) :: IAERO     ! Aerosols
INTEGER(KIND=JPIM) :: IAEROL    ! 
INTEGER(KIND=JPIM) :: IAEROM    !
INTEGER(KIND=JPIM) :: IEXT      ! Extra-GFL
INTEGER(KIND=JPIM) :: IEXTL     ! 
INTEGER(KIND=JPIM) :: IEXTM     !

INTEGER(KIND=JPIM) :: ILIMA  ! LIMA group
INTEGER(KIND=JPIM) :: ILIMAL ! 
INTEGER(KIND=JPIM) :: ILIMAM !
! * other variables (irrelevant horizontal derivatives):
INTEGER(KIND=JPIM) :: IDAL      ! downdraught mesh fraction (for ALARO physics)
INTEGER(KIND=JPIM) :: IDOM      ! downdraught vertical velocity (for ALARO physics)
INTEGER(KIND=JPIM) :: IUAL      ! updraught mesh fraction (for ALARO physics)
INTEGER(KIND=JPIM) :: IUOM      ! updraught vertical velocity (for ALARO physics)
INTEGER(KIND=JPIM) :: IUEN      ! pseudo-historic entrainment (for ALARO physics)
INTEGER(KIND=JPIM) :: IUNEBH    ! pseudo-historic convective cloudiness (for ALARO physics)
INTEGER(KIND=JPIM) :: ITTE      ! total turbulent energy
INTEGER(KIND=JPIM) :: IMXL      ! prognostic mixing lentgh
INTEGER(KIND=JPIM) :: ISHTUR    ! shear source term for turbulence           
INTEGER(KIND=JPIM) :: IFQTUR    ! flux form source term for turbulence -moisture
INTEGER(KIND=JPIM) :: IFSTUR    ! flux form source term for turbulence -enthalpy
INTEGER(KIND=JPIM) :: IRKTH     ! Rasch-Kristjansson enthalpy tendency.
INTEGER(KIND=JPIM) :: IRKTQV    ! Rasch-Kristjansson water vapour tendency.
INTEGER(KIND=JPIM) :: IRKTQC    ! Rasch-Kristjansson condensates tendency.
INTEGER(KIND=JPIM) :: ISRC      ! Second-order flux (for AROME)
INTEGER(KIND=JPIM) :: ICPF      ! Convective precipitation flux
INTEGER(KIND=JPIM) :: ISPF      ! Stratiform precipitation flux
INTEGER(KIND=JPIM) :: ICVGQ     ! Moisture convergence (for MF physics)
INTEGER(KIND=JPIM) :: ISDSAT    ! Standard deviation of the saturation depression
INTEGER(KIND=JPIM) :: ICVV      ! Convective vertical velocity
INTEGER(KIND=JPIM) :: IQVA      ! Total humidity variation (for HIRLAM physics)
INTEGER(KIND=JPIM) :: ILRCH4    ! CH4 (methane) loss rate
INTEGER(KIND=JPIM) :: IEMIS3D   ! 3D emissions for atmospheric composition
INTEGER(KIND=JPIM) :: IFORC     ! Forcings (1D model)
INTEGER(KIND=JPIM) :: IEZDIAG   ! Easy diagnostics (for AROME physics)
INTEGER(KIND=JPIM) :: IGHG      ! Greenhouse gases (ECMWF)
INTEGER(KIND=JPIM) :: ICHEM     ! Chemistry fields (ECMWF)
INTEGER(KIND=JPIM) :: IERA40    ! ERA40 reanalysis fields
INTEGER(KIND=JPIM) :: INOGW     ! Diagnostic fields for NORO GWD scheme (ECMWF)
INTEGER(KIND=JPIM) :: IAERAOT   ! Aerosol layer optical thicknesses
INTEGER(KIND=JPIM) :: IAERLISI  ! Aerosol lidar simulator (extinction and backscatter)
INTEGER(KIND=JPIM) :: IAEROUT   ! Output aerosols for diagnostics
INTEGER(KIND=JPIM) :: IAERCLIM ! 3D aerosol input
INTEGER(KIND=JPIM) :: IUVP      ! Output fields for UV processor
INTEGER(KIND=JPIM) :: ISLDIA    ! Semi-Lagrangian dynamics diagnostic fields
INTEGER(KIND=JPIM) :: IPHYCTY   ! Specific mass change from physics (s-1) for CTY equation
INTEGER(KIND=JPIM) :: IFSD      ! cloud heterogeneity FSD
INTEGER(KIND=JPIM) :: IEDRP     ! turbulence EDR Parameter diagnostics
END TYPE TYPE_IGFLFLDD

TYPE TYPE_LGFLFLDD
! List of GFL variables with derivatives for logical variables.
! * historical (advectable) variables:
LOGICAL :: LLQ        ! Specific humidity
LOGICAL :: LLQL       !
LOGICAL :: LLQM       !
LOGICAL :: LLL        ! Liquid water
LOGICAL :: LLLL       !
LOGICAL :: LLLM       !
LOGICAL :: LLLCONV    ! Liquid water (CONV. PART)
LOGICAL :: LLLCONVL   !
LOGICAL :: LLLCONVM   !
LOGICAL :: LLICONV    ! Ice    water (CONV. PART)
LOGICAL :: LLICONVL   !
LOGICAL :: LLICONVM   !
LOGICAL :: LLRCONV    ! Rain         (CONV. PART)
LOGICAL :: LLRCONVL   !
LOGICAL :: LLRCONVM   !
LOGICAL :: LLSCONV    ! Snow         (CONV. PART)
LOGICAL :: LLSCONVL   !
LOGICAL :: LLSCONVM   !
LOGICAL :: LLLRAD     ! Total liquid water for radiation
LOGICAL :: LLLRADL    !
LOGICAL :: LLLRADM    !
LOGICAL :: LLI        ! Ice
LOGICAL :: LLIL       !
LOGICAL :: LLIM       !
LOGICAL :: LLIRAD     ! Total ice for radiation
LOGICAL :: LLIRADL    !
LOGICAL :: LLIRADM    !
LOGICAL :: LLA        ! Cloud fraction
LOGICAL :: LLAL       !
LOGICAL :: LLAM       !
LOGICAL :: LLO3       ! Ozone
LOGICAL :: LLO3L      !
LOGICAL :: LLO3M      !
LOGICAL :: LLS        ! Snow
LOGICAL :: LLSL       !
LOGICAL :: LLSM       !
LOGICAL :: LLRR       ! Rain
LOGICAL :: LLRRL      !
LOGICAL :: LLRRM      !
LOGICAL :: LLG        ! Graupels
LOGICAL :: LLGL       !
LOGICAL :: LLGM       !
LOGICAL :: LLH        ! Hail
LOGICAL :: LLHL       !
LOGICAL :: LLHM       !
LOGICAL :: LLTKE      ! TKE (total kinetic entrainment)
LOGICAL :: LLTKEL     ! 
LOGICAL :: LLTKEM     !
LOGICAL :: LLAERO(JPAERO)    ! Aerosols
LOGICAL :: LLAEROL(JPAERO)   ! 
LOGICAL :: LLAEROM(JPAERO)   !
LOGICAL :: LLEXT(JPEXT)      ! Extra-GFL
LOGICAL :: LLEXTL(JPEXT)     ! 
LOGICAL :: LLEXTM(JPEXT)     !

LOGICAL :: LLLIMA(JPLIMA)     ! LIMA group
LOGICAL :: LLLIMAL(JPLIMA)    ! 
LOGICAL :: LLLIMAM(JPLIMA)    !
! * other variables (irrelevant horizontal derivatives):
LOGICAL :: LLDAL      ! downdraught mesh fraction (for ALARO physics)
LOGICAL :: LLDOM      ! downdraught vertical velocity (for ALARO physics)
LOGICAL :: LLUAL      ! updraught mesh fraction (for ALARO physics)
LOGICAL :: LLUOM      ! updraught vertical velocity (for ALARO physics)
LOGICAL :: LLUEN      ! pseudo-historic entrainment (for ALARO physics)
LOGICAL :: LLUNEBH    ! pseudo-historic convective cloudiness (for ALARO physics)
LOGICAL :: LLTTE      ! total turbulent energy
LOGICAL :: LLMXL      ! prognostic mixing length
LOGICAL :: LLSHTUR    ! shear source term for turbulence           
LOGICAL :: LLFQTUR    ! flux form source term for turbulence -moisture
LOGICAL :: LLFSTUR    ! flux form source term for turbulence -enthalpy
LOGICAL :: LLRKTH     ! Rasch-Kristjansson enthalpy tendency.
LOGICAL :: LLRKTQV    ! Rasch-Kristjansson water vapour tendency.
LOGICAL :: LLRKTQC    ! Rasch-Kristjansson condensates tendency.
LOGICAL :: LLSRC      ! Second-order flux (for AROME)
LOGICAL :: LLCPF      ! Convective precipitation flux
LOGICAL :: LLSPF      ! Stratiform precipitation flux
LOGICAL :: LLCVGQ     ! Moisture convergence (for MF physics)
LOGICAL :: LLSDSAT    ! Standard deviation of the saturation depression
LOGICAL :: LLCVV      ! Convective vertical velocity
LOGICAL :: LLQVA      ! Total humidity variation (for HIRLAM physics)
LOGICAL :: LLLRCH4    ! CH4 (methane) loss rate
LOGICAL :: LLEMIS3D(NPEMIS3D) ! 3D emissions for atmospheric composition
LOGICAL :: LLFORC(JPFORC)     ! Forcings (1D model)
LOGICAL :: LLEZDIAG(JPEZDIAG) ! Easy diagnostics (for AROME physics)
LOGICAL :: LLGHG(JPGHG)       ! Greenhouse gases (ECMWF)
LOGICAL :: LLCHEM(JPCHEM)     ! Chemistry fields (ECMWF) 
LOGICAL :: LLERA40(JPERA40)   ! ERA40 reanalysis fields
LOGICAL :: LLNOGW(JPNOGW)     ! Diagnostic fields for NORO GWD scheme (ECMWF)
LOGICAL :: LLAERAOT(NPAERAOT) ! Aerosol layer optical thicknesses
LOGICAL :: LLAERLISI(NPAERLISI)!Aerosol lidar simulator (extinction and backscatter)
LOGICAL :: LLAEROUT(JPAEROUT) ! Output aerosols for diagnostics
LOGICAL :: LLAEROCLIM(JPAEROCLIM) ! 3D aerosol inputs
LOGICAL :: LLUVP(JPUVP)       ! Output fields for UV processor
LOGICAL :: LLSLDIA(JPSLDIA)   ! Semi-Lagrangian dynamics diagnostic fields
LOGICAL :: LLPHYCTY   ! Specific mass change from physics (s-1) for CTY equation
LOGICAL :: LLEDRP(JPEDRP)     ! Diagnostic EDR Parameter turbulence fields (ECMWF)
END TYPE TYPE_LGFLFLDD

!-------------------------------------------------------------------------
END MODULE TYPE_GFLFLDS
