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


!     Initialize variables that control :
!     -> Trajectory management at t-dt
!     -> Physics in tangent-linear and adjoint

!     ------------------------------------------------------------------
NAMELIST/NAMTRAJP/ LETRAJP,LETRAJPT,&
      &NLWFR,NLOOPLW,&
      &LERADI2,LERADS2,LERADSW2,LERADLW2,LERADN2,LERADFL2,&
      &LOPTLWPR,LWLCLHR,LEPCLD2,&
      &LEDCLD2,LENCLD2,LEVAPLS2,LEVDIF2,LEGWDG2,LECUMF2,&
      &LECOND2,LEGWWMS2,LREGWWMS,LEQNGT2,LESURF2,LREGSF,&
      &LEKPERT,LEKPERTS,LENOPERT,LNCLIN,LREGCL,&
      &LTRACLNPH
!     ------------------------------------------------------------------

