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

MODULE YOMCTESSELDIM

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!*

!   CTESSEL: Carbon model (for PASSIVE introduction)

!   some usefull dimensions

!  ---------------------------------------------------------------------
INTEGER(KIND=JPIM), PARAMETER  :: IKVTYPES=2
INTEGER(KIND=JPIM), PARAMETER  :: IKDHVBIOS=0
INTEGER(KIND=JPIM), PARAMETER  :: IKDHFBIOS=5  ! leaf biomass, leaf biomass loss, leaf biomass gain, above ground
   ! structural biomass, below ground structural biomass
INTEGER(KIND=JPIM), PARAMETER  :: IKDHVCO2S=0
INTEGER(KIND=JPIM), PARAMETER  :: IKDHFCO2S=6  ! gross assimilation, dark respiration, net assimilation,
   ! soil+structural biomass respiration, ecosystem respiration, net CO2 flux
INTEGER(KIND=JPIM), PARAMETER  :: IKDHVVEGS=7  ! fraction, LAI (m2/m2),
   ! Canopy conductance, Aerodynamic conductance, stress function F2,
   ! Ds (spec. humidity deficit at canopy level),Dmax
INTEGER(KIND=JPIM), PARAMETER  :: IKDHFVEGS=1  ! Latent heat flux (W/m2)


!  ---------------------------------------------------------------------


END MODULE YOMCTESSELDIM
