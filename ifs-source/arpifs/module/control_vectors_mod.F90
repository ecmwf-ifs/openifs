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

MODULE CONTROL_VECTORS_MOD

!   Purpose.
!   --------
!     Define CONTROL_VECTOR type and related operations.
!
!   Author.
!   -------
!     Y. Tremolet
!
!   Modifications.
!   -------------- 
!     Original   02-08-01  
!     Y.Tremolet 02-08-20  Split in several files.
!     C.Fischer  03-03-05  Aladin spectral treatment added.
!      N.B.1: the mean wind profiles, 2*nflev values, resp. mean U and mean V,
!             are treated as an extra appended stripe to the control vector
!             (addressed between ibgn_la and iend_la, as proposed by yannick).
!             however, they do actually belong to the Aladin state vector as well,
!             so that I decided to include them inside the spectral_fields_mod
!             f90 structure (spec%sp1d). this somewhat mixed situation looks
!             uggly, but appears to me as the best compromise between 2 conflicting
!             properties of these bloody dangerous coefficients: they are part
!             of the spectral state, but they are not at all distributed the
!             same way as the usual (horizontal) spectral fields.
!             Actually, the mean wind profiles are assumed to be duplicated on
!             all processors (the Aladin/IFS code should make sure that SPA1
!             is indeed present and equal on all procs, as output from the direct
!             spectral transforms).
!      N.B.2: the 1D spectral profiles in the control vector are controlled
!             by two distinct variables:
!             * key yd_cva%lam1d  tells you that this appendix indeed should be
!               present in the control
!             * integer yd_cva%ns1d tells you how many fields it contains. This is
!               usually 2 (or 0 if absent) for zonal and meridional mean winds
!             Note that the CONTROL_VECTORS_MOD module now handles two distinct
!             logicals for aladin: one simply tells that we must handle the
!             Aladin horizontal spectra (key LELAM in CONTROL_VECTOR_DATA) while
!             the other tells that the mean wind is or not in the cv (key yd_cva%lam1d
!             which used to be yd_cva%lelam in Yannick's original version).
!      N.B.3: for Ald, very strictly coding, one should also modify the
!             nserr truncation into nserr+nmserr to level up with nsmax+nmsmax.
!             the present use of only one nserr dictates that the model error
!             control variable only works for square domains.
!      N.B.4: the mean wind component would also be missing in the model error
!             term since it belongs to the state vector.
!     M.Fisher 2003-11-27  Wavelet Jb and coding-norm compliance
!        Y.Tremolet    18-Mar-2004 Add CTLVEC_NORM
!        Y.Tremolet    23-Mar-2004 Made initial condition optional
!        C.Soci 16-Oct-2004  initialize nmsmax in setup_ctlvec
!        Y.Tremolet    25-Jul-2004 Split and some clean-up
!        M.Fisher      15-Feb-2013 Add CTLVEC_STRUCT
!      Y. Michel, MF, June 2018 Extention of the control variable for sqrt EnVar
! ------------------------------------------------------------------

USE CONTROL_VECTORS_BASE_MIX, ONLY: CONTROL_VECTOR, &
                              & ALLOCATE_CTLVEC, DEALLOCATE_CTLVEC

USE CONTROL_VECTORS_OPER_MOD, ONLY: ASSIGNMENT(=), INVERSE

USE CONTROL_VECTORS_PARA_MOD, ONLY: DOT_PRODUCT, MAXVAL, MINVAL, SUM, &
                              & CTLVEC_SQNORM, CTLVEC_NORM

USE CONTROL_VECTORS_DATA_MIX, ONLY: SETUP_CTLVEC, CTLVEC_STRUCT,CTLVEC_STRUCT_ENS, CONTROL_VECTOR_DATA_STRUCT

IMPLICIT NONE
PRIVATE
PUBLIC CONTROL_VECTOR, SETUP_CTLVEC, ALLOCATE_CTLVEC, DEALLOCATE_CTLVEC, &
     & CONTROL_VECTOR_DATA_STRUCT, CTLVEC_STRUCT, CTLVEC_STRUCT_ENS,&
     & ASSIGNMENT(=), INVERSE, &
     & DOT_PRODUCT, MAXVAL, MINVAL, SUM, CTLVEC_SQNORM, CTLVEC_NORM

! ------------------------------------------------------------------

END MODULE CONTROL_VECTORS_MOD
