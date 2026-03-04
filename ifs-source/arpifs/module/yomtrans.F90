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

MODULE YOMTRANS

USE PARKIND1  ,ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    SPECTRAL TRANSFORMS CONTROL VARIABLES

!     NPROMATR : max. number of fields to be transformed at each call in
!                arpege/ifs inverse transforms
!                = 0 : all fields
!                > 0 : such number of fields
!     NEPROMATR : max. number of fields to be transformed at each call in
!                aladin inverse transforms
!                = 0 : all fields
!                > 0 : such number of fields
!     NMAX_RESOL: number of computed resolutions
!     RDISTR_E  : between 0. and 1. (debault 0.) part of E zone that will
!                 be distributed
!     LUSEFLT,LUSERPNM,LKEEPRPNM - see TRANS package
!     LALLOPERM : allocate certain arrays permanently
!     LFFTW     : Use FFTW if true (see TRANS package)


INTEGER(KIND=JPIM) :: NPROMATR
INTEGER(KIND=JPIM) :: NEPROMATR
INTEGER(KIND=JPIM) :: NMAX_RESOL
REAL(KIND=JPRB) :: RDISTR_E 
LOGICAL         :: LMONO_TRANS
LOGICAL         :: LUSEFLT
LOGICAL         :: LUSERPNM
LOGICAL         :: LKEEPRPNM
LOGICAL         :: LALLOPERM
LOGICAL         :: LFFTW
!     ------------------------------------------------------------------
END MODULE YOMTRANS
