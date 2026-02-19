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

MODULE YOMANEB

USE PARKIND1  ,ONLY : JPIM, JPRB

USE YOM_YGFL, ONLY : JPGHG_ASSIM, JPCHEM_ASSIM

USE YOMJG, ONLY : JPAEROCV

IMPLICIT NONE

SAVE

!       H.Varella : 2010-04-20 Divergence included
!       Y. Michel : 2014-11-18 Buffers for Local Correlation Tensor
!     ------------------------------------------------------------------

!*    Buffer for gridpoint workfile for analysis and new forecast errors

! ANEBUF                : buffer
! NANEBU                : locates u errors in buffer
! NANEBV                : locates v errors in buffer
! NANEBZ                : locates geopotential errors in buffer
! NANEBQ                : locates Q errors in buffer
! NANEBR                : locates RH errors in buffer
! NANEBT                : locates T errors in buffer
! NANEBDIV              : locates divergence errors in buffer
! NANEBVO               : locates vorticity errors in buffer
! NANEBRAD              : locates TOVS radiance errors in buffer
! NANEBTCW              : locates Total Column Water errors in buffer
! NANEBRH2              : locates 2M Relative Humidity errors in buffer
! NANEBU10              : locates 10M u-component errors in buffer
! NANEBV10              : locates 10M v-component errors in buffer
! NANEBT2               : locates 2M temperature errors in buffer
! NANEBO3               : locates ozone errors in buffer
! NANEBAERO             : locates AERO in buffer
! NANEBTO3              : locates total column ozone errors in buffer
! NANEBSP               : locates surface pressure errors in buffer
! NANEVARS              : number of analysis error variables
! NANERADS              : number of analysis error radiances
! NANERADL              : list of analysis error radiances channels

! NGRBCHEM(JPCHEM_ASSIM)    : Gribcodes of assimilated CHEM fields
! NGRBGHGASSIM(JPGHG_ASSIM) : Gribcodes of assimilated GHG fields

!*    Buffer for gridpoint workfile for Local Correlation Tensor

! DERZANEBUF            : buffer
! DERMANEBUF            : buffer
! DERZMANEBUF           : buffer

INTEGER(KIND=JPIM) :: NANEBU
INTEGER(KIND=JPIM) :: NANEBV
INTEGER(KIND=JPIM) :: NANEBZ
INTEGER(KIND=JPIM) :: NANEBQ
INTEGER(KIND=JPIM) :: NANEBR
INTEGER(KIND=JPIM) :: NANEBT
INTEGER(KIND=JPIM) :: NANEBDIV
INTEGER(KIND=JPIM) :: NANEBVO
INTEGER(KIND=JPIM) :: NANEBRAD
INTEGER(KIND=JPIM) :: NANEBTCW
INTEGER(KIND=JPIM) :: NANEBRH2
INTEGER(KIND=JPIM) :: NANEBU10
INTEGER(KIND=JPIM) :: NANEBV10
INTEGER(KIND=JPIM) :: NANEBT2
INTEGER(KIND=JPIM) :: NANEBO3
INTEGER(KIND=JPIM), DIMENSION(JPGHG_ASSIM)  :: NANEBGHG
INTEGER(KIND=JPIM), DIMENSION(JPCHEM_ASSIM) :: NANEBCHEM
INTEGER(KIND=JPIM), DIMENSION(JPAEROCV)     :: NANEBAERO
INTEGER(KIND=JPIM) :: NANEBTO3
INTEGER(KIND=JPIM) :: NANEBSP
INTEGER(KIND=JPIM) :: NANEVARS
INTEGER(KIND=JPIM) :: NANERADS

REAL(KIND=JPRB),    ALLOCATABLE :: ANEBUF (:,:,:)
REAL(KIND=JPRB),    ALLOCATABLE :: DERZANEBUF (:,:,:), DERMANEBUF (:,:,:), DERZMANEBUF (:,:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: NANERADL(:)
INTEGER(KIND=JPIM), DIMENSION(JPCHEM_ASSIM) :: NGRBCHEM
INTEGER(KIND=JPIM), DIMENSION(JPGHG_ASSIM)  :: NGRBGHGASSIM

!     ------------------------------------------------------------------
END MODULE YOMANEB

