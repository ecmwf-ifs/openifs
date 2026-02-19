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

!****--------------------------------------------------------------------
!**** NAMCUMF : Namelist for the massflux scheme
!**** 
!****  Author   : M. Leutbecher  ECMWF 01/2006 
!****
!****  Modifications:
!****    P. Bechtold 21-09-2011  Adding more parameters for parameter
!****                            optimisation (project using EPPES)
!**** 
!****
!****--------------------------------------------------------------------
!***
!***  RMFSOLUV: implicitness factor for mass flux solver for momentum
!***  RMFSOLTQ: implicitness factor for mass flux solver for T and q
!***  RMFSOLCT: implicitness factor for chemical tracers
!***  ENTRORG : entrainment for positively buoyant deep convection 1/(m)
!***  ENTSHALP: shallow entrainment defined as ENTSHALP*ENTRORG
!***  DETRPEN : detrainment rate for penetrative convection (1/m)
!***  RMFDEPS : fractional massflux for downdrafts at LFS
!***  RTAUA   : proportionality constant for adjsutment time
!***  RPRCON  : coefficient for determining conversion from cloud water to rain
!***  ENTRDD  : entrainment rate for cumulus downdrafts
!***  RDEPTHS : maximum cloud depth (Pa) for shallow convection
!***  LMFPEN  : switch  for deep convection
!***  LMFCUCA : switch to modulate base massflux by CA or 2d advect field
!***
!***---------------------------------------------------------------------
NAMELIST / NAMCUMF / RMFSOLUV, RMFSOLTQ, RMFSOLCT, RMFCFL, RMFSOLRHS, ENTRORG, ENTSHALP, &
 & DETRPEN, RTAUA, RMFDEPS, RPRCON, ENTRDD, RDEPTHS, LMFPEN, LMFSCV, LMFDSNOW, LMFCUCA, &
 & RCAPDCYCL, RMINCAPE, RHEBC, RBASE0, RMINCIN, ENTSTPC1 ,ENTSTPC2, ENTSTPC3, NJKT7, RMFLIA,&
 & LMFWETB, LMFGLAC, LSCVLIQ, LMFDUDV, LMFENTHCONS, RUVPER, RMFADVW, RMFADVWDD, RCPECONS, RCAPQADV
