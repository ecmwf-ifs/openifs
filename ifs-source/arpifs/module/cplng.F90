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

MODULE CPLNG

   USE CPLNG_TYPES_MOD
   USE CPLNG_DATA_MOD
   USE CPLNG_INIT_MOD
   USE CPLNG_EXCHANGE_MOD
   USE CPLNG_FINALIZE_MOD

   IMPLICIT NONE

   PRIVATE

   PUBLIC CPLNG_IS_ACTIVE

   PUBLIC CPLNG_FLD_TYPE

   PUBLIC CPLNG_ADD_FLD
   PUBLIC CPLNG_ADD_FLD_COMPLETED
   PUBLIC CPLNG_IDX
    
   PUBLIC CPLNG_FLD_TYPE_GRIDPOINT
   PUBLIC CPLNG_FLD_TYPE_SPECTRAL

   PUBLIC CPLNG_INIT
   PUBLIC CPLNG_EXCHANGE
   PUBLIC CPLNG_FINALIZE

   PUBLIC CPL_IN
   PUBLIC CPL_OUT
   PUBLIC CPL_OUTINST

END MODULE CPLNG
