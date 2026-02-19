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

MODULE PTRXFU

USE PARKIND1  ,ONLY : JPIM    

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------

!*    Pointers for instantaneous fluxes diagnostics (XFU).

TYPE TXFUPTR

INTEGER(KIND=JPIM) :: MXTRDU  ! pointer for "U-wind gravity wave stress" flux 
INTEGER(KIND=JPIM) :: MXTRDV  ! pointer for "V-wind gravity wave stress" flux
INTEGER(KIND=JPIM) :: MXTRCU  ! pointer for "contribution of convection to U" flux
INTEGER(KIND=JPIM) :: MXTRCV  ! pointer for "contribution of convection to V" flux
INTEGER(KIND=JPIM) :: MXDICQ  ! pointer for "contribution of convection to q" flux
INTEGER(KIND=JPIM) :: MXDICS  ! pointer for "contribution of convection to (cp T)" flux
INTEGER(KIND=JPIM) :: MXTRTU  ! pointer for "contribution of turbulence to U" flux
INTEGER(KIND=JPIM) :: MXTRTV  ! pointer for "contribution of turbulence to V" flux
INTEGER(KIND=JPIM) :: MXDITQ  ! pointer for "contribution of turbulence to T" flux
INTEGER(KIND=JPIM) :: MXDITS  ! pointer for "contribution of turbulence to (cp T)" flux
INTEGER(KIND=JPIM) :: MXPLCL  ! pointer for "convective precipitation"
INTEGER(KIND=JPIM) :: MXPLSL  ! pointer for "large scale precipitation (stratiform)"
INTEGER(KIND=JPIM) :: MXPLCN  ! pointer for "convective snow fall"
INTEGER(KIND=JPIM) :: MXPLSN  ! pointer for "large scale snow fall (stratiform)"
INTEGER(KIND=JPIM) :: MXPLCG  ! pointer for "convective graupel fall"
INTEGER(KIND=JPIM) :: MXPLSG  ! pointer for "large scale graupel fall (stratiform)"
INTEGER(KIND=JPIM) :: MXPLCH  ! pointer for "convective hail fall"
INTEGER(KIND=JPIM) :: MXPLSH  ! pointer for "large scale hail fall (stratiform)"
INTEGER(KIND=JPIM) :: MXRSO   ! pointer for "solar radiation" flux
INTEGER(KIND=JPIM) :: MXRTH   ! pointer for "surface radiation" flux
INTEGER(KIND=JPIM) :: MXQICE  ! pointer for "ice water"
INTEGER(KIND=JPIM) :: MXQLI   ! pointer for "liquid water"
INTEGER(KIND=JPIM) :: MXNEB   ! pointer for "cloud cover"
INTEGER(KIND=JPIM) :: MXNEBT  ! pointer for "total cloudiness"
INTEGER(KIND=JPIM) :: MXUCLS  ! pointer for "U-component of wind at 10 meters (pbl)"
INTEGER(KIND=JPIM) :: MXVCLS  ! pointer for "V-component of wind at 10 meters (pbl)"
INTEGER(KIND=JPIM) :: MXNUCLS ! pointer for "U-component of neutral wind at 10 meters (pbl)"
INTEGER(KIND=JPIM) :: MXNVCLS ! pointer for "V-component of neutral wind at 10 meters (pbl)"
INTEGER(KIND=JPIM) :: MXTCLS  ! pointer for "temperature at 2 meters (pbl)
INTEGER(KIND=JPIM) :: MXQCLS  ! pointer for "specific humidity at 2 meters (pbl)"
INTEGER(KIND=JPIM) :: MXRHCLS ! pointer for "relative humidity at 2 meters (pbl)"
INTEGER(KIND=JPIM) :: MXTPWCLS! pointer for "wet bulb temperature at 2 meters (pbl)"
INTEGER(KIND=JPIM) :: MXTN    ! pointer for "minimum temperature at 2 meters"
INTEGER(KIND=JPIM) :: MXTX    ! pointer for "maximum temperature at 2 meters"
INTEGER(KIND=JPIM) :: MXMRT   ! pointer for "mean radiant temperature at 2 meters"
INTEGER(KIND=JPIM) :: MXHUX   ! pointer for "maximum relative humidity at 2 meters (pbl)"
INTEGER(KIND=JPIM) :: MXHUN   ! pointer for "minimum relative humidity at 2 meters (pbl)"
INTEGER(KIND=JPIM) :: MXNBCON ! pointer for "convective cloud cover"
INTEGER(KIND=JPIM) :: MXNBHAU ! pointer for "high cloud cover"
INTEGER(KIND=JPIM) :: MXNBMOY ! pointer for "medium cloud cover"
INTEGER(KIND=JPIM) :: MXNBBAS ! pointer for "low cloud cover"
INTEGER(KIND=JPIM) :: MXFONTE ! pointer for "melt snow" flux
INTEGER(KIND=JPIM) :: MXFCHSP ! pointer for "heat flux in soil" flux
INTEGER(KIND=JPIM) :: MXFLWSP ! pointer for "water flux in soil" flux
INTEGER(KIND=JPIM) :: MXFEVL  ! pointer for "liquid evaporation" flux
INTEGER(KIND=JPIM) :: MXRUISP ! pointer for "deep soil runoff" flux
INTEGER(KIND=JPIM) :: MXRUISS ! pointer for "surface soil runoff" flux
INTEGER(KIND=JPIM) :: MXFEVV  ! pointer for "evapotranspiration" flux
INTEGER(KIND=JPIM) :: MXFTR   ! pointer for "transpiration" flux
INTEGER(KIND=JPIM) :: MXRUISL ! pointer for "interception soil layer runoff" flux
INTEGER(KIND=JPIM) :: MXCAPE  ! pointer for "CAPE (convective available potential energy)" flux
INTEGER(KIND=JPIM) :: MXCTOP  ! pointer for "pressure of top of deep convection" flux
INTEGER(KIND=JPIM) :: MXMOCON ! pointer for "moisture convergence" flux
INTEGER(KIND=JPIM) :: MXCLPH  ! pointer for "height of the PBL out of the model" flux
INTEGER(KIND=JPIM) :: MXVEIN  ! pointer for "ventilation index in PBL" flux
INTEGER(KIND=JPIM) :: MXUGST  ! pointer for "U-momentum of gusts out of the model" flux
INTEGER(KIND=JPIM) :: MXVGST  ! pointer for "V-momentum of gusts out of the model" flux
INTEGER(KIND=JPIM) :: MXUGST2 ! pointer for "U-momentum of gusts2
INTEGER(KIND=JPIM) :: MXVGST2 ! pointer for "V-momentum of gusts2 
INTEGER(KIND=JPIM) :: MXTHW   ! pointer for "theta'_w" surface flux
INTEGER(KIND=JPIM) :: MXXDIAGH! pointer for hail diagnostic
INTEGER(KIND=JPIM) :: MXVISICLD! pointer for visibility due to water and/or ice cloud
INTEGER(KIND=JPIM) :: MXVISIHYD! pointer for visibility due to precipitations
INTEGER(KIND=JPIM) :: MXCLWC  ! pointer for maximum of CLWC
INTEGER(KIND=JPIM) :: MXPTYPE ! pointer for frequent precipitations type
INTEGER(KIND=JPIM) :: MXPTYPESEV! pointer for severe precipitations type
! ****2 = idem for another period
INTEGER(KIND=JPIM) :: MXVISICLD2! pointer for visibility due to water and/or ice cloud
INTEGER(KIND=JPIM) :: MXVISIHYD2! pointer for visibility due to precipitations
INTEGER(KIND=JPIM) :: MXCLWC2  ! pointer for maximum of CLWC
INTEGER(KIND=JPIM) :: MXPTYPE2 ! pointer for frequent precipitations type
INTEGER(KIND=JPIM) :: MXPTYPESEV2! pointer for severe precipitations type

END TYPE TXFUPTR
!     -----------------------------------------------------------------
END MODULE PTRXFU
