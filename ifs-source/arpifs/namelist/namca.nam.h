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
!**** NAMCA : Namelist for the Cellular Automaton
!**** 
!****  Author   : M. Steinheimer  ECMWF 06/2009 
!****
!****  Modifications:
!**** 
!****
!****--------------------------------------------------------------------

!***
!***  LCUCONV_CA: switch for activateing the CA
!***  LCA_ADVECT: switch for "kind of" semi-Lagrangian advection
!***  LCA_GLOBAL: switch for global CA
!***  LCA_SMOOTH: swith for smoothing CA on large scale grid
!***  LCA_TEST: switch for initialize CA at single point
!***  LCA_RANTROP  : switch for random coupling of CA to deep convection (concerns initialization)
!***  NTESTPROC: set on which processor the single point should lie
!***  NTESTGP: set which gridpoint the single point is on
!***  CA_FORC: switch to choose convective forcing
!***  CA_WIND: switch to choose CA-wind (real/idealized)
!***  CA_PROB: switch to choose probabilities
!***  NFRCASEED: switch to choose frequency of seeding new cells in CA
!***  NLIVES: switch to choose number of lives (scaled by CAPE in callpar)
!***  NFERTYRS: switch to choose max number of steps a cell can be fertile
!***  RCA_SEEDPROB: Probability of random seeding for global CA
!***  RCADELX : Length of individual CA cell (meters)
!***  NDXUNREAL: set X component of unreal wind
!***  NDYUNREAL: set Y component of unreal wind
!***  NSPINUP  : set number of spin-up cycles for global pattern
!***  LCA_ADVTEST  : switch for advection test (CA only evolved first 5 steps)
!***  LCA_EXTRACT : switch for extraction of CA fields (to plot output)
!***  LCA_NBDEBUG  : switch for neighbour debug output
!***---------------------------------------------------------------------

NAMELIST / NAMCA / LCUCONV_CA, LCA_ADVECT, LCA_GLOBAL, LCA_SMOOTH, LCA_TEST, &
 & CA_FORC, CA_WIND, CA_PROB, NFRCASEED, NLIVES, NFERTYRS, RCA_SEEDPROB, &
 & NTESTPROC,NTESTGP,NDXUNREAL,NDYUNREAL,NSPINUP,LCA_RANTROP,&
 & LCA_ADVTEST,LCA_EXTRACT,LCA_NBDEBUG, RCADELX
