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
!**** NAMVDF : Namelist for the vert diffusion scheme
!****
!****  Author   : P. Bechtold ECMWF 12/2014
!****
!***  RLAM   : ASYMPTOTIC MIXING LENGTH FOR MOMENTUM
!***  RVDIFTS: FACTOR FOR TIME STEP WEIGHTING IN *VDF....*
!***  LWDS   : .T. for Wood/Diamantakis/Staniforth scheme      
!***  NSUBST : Number of substeps in VDF           

!****
!****--------------------------------------------------------------------
NAMELIST / NAMVDF / RLAM, NSUBST, LWDS, RVDIFTS, RTOFDALPHA, REISTHSC


