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

MODULE PTRGFU

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------

!*    Pointers for cumulated fluxes diagnostics (CFU).

TYPE TCFUPTR

INTEGER(KIND=JPIM) :: MSTRDU  ! pointer for "U-wind gravity wave stress" flux
INTEGER(KIND=JPIM) :: MSTRDV  ! pointer for "V-wind gravity wave stress" flux
INTEGER(KIND=JPIM) :: MSTRCU  ! pointer for "contribution of convection to U" flux
INTEGER(KIND=JPIM) :: MSTRCV  ! pointer for "contribution of convection to V" flux
INTEGER(KIND=JPIM) :: MFDICQ  ! pointer for "contribution of convection to q" flux
INTEGER(KIND=JPIM) :: MFDICS  ! pointer for "contribution of convection to (cp T)" flux
INTEGER(KIND=JPIM) :: MSTRTU  ! pointer for "contribution of turbulence to U" flux
INTEGER(KIND=JPIM) :: MSTRTV  ! pointer for "contribution of turbulence to V" flux
INTEGER(KIND=JPIM) :: MFDITQ  ! pointer for "contribution of turbulence to T" flux
INTEGER(KIND=JPIM) :: MFDITS  ! pointer for "contribution of turbulence to (cp T)" flux
INTEGER(KIND=JPIM) :: MFPLCL  ! pointer for "convective precipitation"
INTEGER(KIND=JPIM) :: MFPLSL  ! pointer for "large scale precipitation (stratiform)"
INTEGER(KIND=JPIM) :: MFPLCN  ! pointer for "convective snow fall"
INTEGER(KIND=JPIM) :: MFPLSN  ! pointer for "large scale snow fall (stratiform)"
INTEGER(KIND=JPIM) :: MFPLCG  ! pointer for "convective graupel fall"
INTEGER(KIND=JPIM) :: MFPLSG  ! pointer for "large scale graupel fall (stratiform)"
INTEGER(KIND=JPIM) :: MFPLCH  ! pointer for "convective hail fall"
INTEGER(KIND=JPIM) :: MFPLSH  ! pointer for "large scale hail fall (stratiform)"
INTEGER(KIND=JPIM) :: MFCCQL  ! pointer for "liquid condensation due to convection" flux
INTEGER(KIND=JPIM) :: MFCCQN  ! pointer for "solid condensation due to convection" flux
INTEGER(KIND=JPIM) :: MFCSQL  ! pointer for "large scale liquid condensation" flux
INTEGER(KIND=JPIM) :: MFCSQN  ! pointer for "large scale solid condensation" flux
INTEGER(KIND=JPIM) :: MFRSO   ! pointer for "solar radiation" flux
INTEGER(KIND=JPIM) :: MFRTH   ! pointer for "surface radiation" flux
INTEGER(KIND=JPIM) :: MFNEB   ! pointer for "cloud cover"
INTEGER(KIND=JPIM) :: MWS     ! pointer for "soil moisture"
INTEGER(KIND=JPIM) :: MSNS    ! pointer for "snow mass"
INTEGER(KIND=JPIM) :: MQTOT   ! pointer for "total precipitable water"
INTEGER(KIND=JPIM) :: MFNEBT  ! pointer for "total cloudiness"
INTEGER(KIND=JPIM) :: MOZONT  ! pointer for "total ozone"
INTEGER(KIND=JPIM) :: MSPRES  ! pointer for "surface pressure"
INTEGER(KIND=JPIM) :: MFCL    ! pointer for "latent heat" flux
INTEGER(KIND=JPIM) :: MFCS    ! pointer for "sensible heat" flux
INTEGER(KIND=JPIM) :: MFDISH  ! pointer for "surface enthalpy" flux (due to the dissipation of kinetic energy)
INTEGER(KIND=JPIM) :: MFTOPH  ! pointer for "top mesospheric enthalpy" flux (+ dissipation).
INTEGER(KIND=JPIM) :: MQICE   ! pointer for "solid specific moisture"
INTEGER(KIND=JPIM) :: MQLI    ! pointer for "liquid specific moisture"
INTEGER(KIND=JPIM) :: MFRSOC0 ! pointer for "top clear sky shortwave radiative" flux
INTEGER(KIND=JPIM) :: MFRTHC0 ! pointer for "top clear sky longwave radiative" flux
INTEGER(KIND=JPIM) :: MFRSOC1 ! pointer for "surface clear sky shortwave radiative" flux
INTEGER(KIND=JPIM) :: MFRTHC1 ! pointer for "surface clear sky longwave radiative" flux
INTEGER(KIND=JPIM) :: MFRSDNI ! pointer for "surface direct normal irradiance" flux
INTEGER(KIND=JPIM) :: MFRSGNI ! pointer for "surface global normal irradiance" flux
INTEGER(KIND=JPIM) :: MFRSOPS ! pointer for "surface parallel solar" flux
INTEGER(KIND=JPIM) :: MFRSOPT ! pointer for "top parallel solar" flux
INTEGER(KIND=JPIM) :: MFRSOLU ! pointer for "surface downward moon radiation" flux
INTEGER(KIND=JPIM) :: MFRSODS ! pointer for "surface down solar" flux
INTEGER(KIND=JPIM) :: MFRTHDS ! pointer for "surface down thermic" flux
INTEGER(KIND=JPIM) :: MFCLL   ! pointer for "liquid latent heat" flux
INTEGER(KIND=JPIM) :: MFCLN   ! pointer for "solid latent heat" flux
INTEGER(KIND=JPIM) :: MFEVL   ! pointer for "liquid evaporation" flux
INTEGER(KIND=JPIM) :: MFEVN   ! pointer for "snow evaporation" flux
INTEGER(KIND=JPIM) :: MFONTE  ! pointer for "melt snow" flux
INTEGER(KIND=JPIM) :: MFCHSP  ! pointer for "heat flux in soil" flux
INTEGER(KIND=JPIM) :: MFEVV   ! pointer for "evapotranspiration" flux
INTEGER(KIND=JPIM) :: MFLWSP  ! pointer for "water flux in soil" flux
INTEGER(KIND=JPIM) :: MFTR    ! pointer for "transpiration" flux
INTEGER(KIND=JPIM) :: MRUISL  ! pointer for "interception soil layer runoff" flux
INTEGER(KIND=JPIM) :: MRUISP  ! pointer for "deep soil runoff" flux
INTEGER(KIND=JPIM) :: MRUISS  ! pointer for "surface soil runoff" flux
INTEGER(KIND=JPIM) :: MNEBCON ! pointer for "convective cloud cover"
INTEGER(KIND=JPIM) :: MNEBHAU ! pointer for "high cloud cover"
INTEGER(KIND=JPIM) :: MNEBMOY ! pointer for "medium cloud cover"
INTEGER(KIND=JPIM) :: MNEBBAS ! pointer for "low cloud cover"
INTEGER(KIND=JPIM) :: MFGEL   ! pointer for "deep frost" flux
INTEGER(KIND=JPIM) :: MFGELS  ! pointer for "surface frost" flux
INTEGER(KIND=JPIM) :: MFDUTP  ! pointer for "filtered duration of total precipitations" flux
INTEGER(KIND=JPIM) :: MFLASH  ! pointer for "cumulative lighntning diagnotics"

END TYPE TCFUPTR

!     -----------------------------------------------------------------
END MODULE PTRGFU
