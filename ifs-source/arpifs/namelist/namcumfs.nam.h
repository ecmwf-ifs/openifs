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
!**** NAMCUMFS : Namelist for the simplified convection scheme
!**** 
!****  Author   : P. Lopez  ECMWF 01/2002 
!****
!****  Modifications:
!**** 
!****  P. Lopez   06-Oct-2005   Removed microphysics parameters
!****  P. Lopez   22-Oct-2007   Added implcitness factors
!****  P. Lopez   10-May-2016   Added switches for rain glaciation and wet bulb 
!****                           in snow melting as well as RHS implicitness factor
!****  P. Lopez   9-Apr-2019    Added switch for reducing mass flux CFL criterion 
!****                           for short time steps 
!****
!****--------------------------------------------------------------------
!***
!***  LECUMFS  : SWITCH ON SIMPLIFIED CONVECTIVE MASS-FLUX SCHEME IN TRAJECTORY
!***  LREGCV   : SWITCH ON REGULARIZATIONS OF SIMPLIFIED CONVECTION SCHEME 
!***             IN TANGENT-LINEAR AND ADJOINT CALCULATIONS
!***  LMFDUDV2 : SWITCH ON CONVECTIVE MOMENTUM TRANSPORT CALCULATIONS
!***  LMFDD2   : SWITCH ON CONVECTIVE DOWNDRAFT COMPUTATIONS
!***  LMFWETB2 : USE WET BULB TEMPERATURE IN SNOW MELTING COMPUTATIONS
!***  LMFGLAC2 : ADD RAIN GLACIATION CONTRIBUTION
!***  LMFCFL2_SHSTEP : APPLY REDUCTION OF MASS FLUX CFL FOR SHORT TIME STEPS
!***  RMFSOLUV2: implicitness factor for mass flux solver for momentum
!***  RMFSOLTQ2: implicitness factor for mass flux solver for T and q
!***  RMFSOLCT2: implicitness factor for chemical tracers
!***  RMFSOLRHS2: implicitness factor for model tendency in RHS of implicit solver
!***  RMFCFL2  : Mass-flux limiter for T and q (times CFL criterion)
!***
!***---------------------------------------------------------------------
NAMELIST / NAMCUMFS / LECUMFS, LREGCV, LMFDUDV2, LMFDD2, LMFWETB2, LMFGLAC2, LMFCFL2_SHSTEP, &
                    & RMFSOLUV2, RMFSOLTQ2, RMFSOLCT2, RMFSOLRHS2, RMFCFL2, &
                    & RCAPDCYCL2, RTAU02, RMFLIA2



