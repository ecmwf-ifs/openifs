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

MODULE YOMCT1

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Control variables for the model - changed at level 1 during ex.

!     N1POS  : OVER-RIDING SWITCH FOR POST-PROCESSING (0=FALSE)
!     N1HIS  : OVER-RIDING SWITCH FOR HISTORY WRITE-UP (0=FALSE)
!     N1GDI  : OVER-RIDING SWITCH FOR GRID-POINT DIAGNOSTICS (0=FALSE)
!     N1SDI  : OVER-RIDING SWITCH FOR SPECTRAL DIAGNOSTICS (0=FALSE)
!     N1DHP  :     "         "     "  DIAGNOSTIC PRINTS DDH (0=F)
!     N1DHFG :     "         "     "  DIAGNOSTIC FILE OUTPUTS
!                                                GLOBAL DDH (0=F)
!     N1DHFZ :     "         "     "  DIAGNOSTIC FILE OUTPUTS
!                                                ZONAL MEANS DDH (0=F)
!     N1DHFD :     "         "     "  DIAGNOSTIC FILE OUTPUTS
!                                                LIMITED AREAS DDH (0=F)
!     N1CFU  :     "         "     "  ACCUMULATED FLUX WRITE-UPS
!     N1XFU  :     "         "     "  INSTANTANEOUS FLUX WRITE-UPS
!     N1RES  :     "         "     "  RESTART WRITE-UP
!     N1SFXHIS :   "         "     " SURFEX HISTORY WRITE-UP (0=FALSE)
!     N1MASSCON :   "         "     " MASS CONSERVATION (0=FALSE)
!     LRFILAF: .T. = CALL LFILAF : CATALOG FILE LFI
!     LWRSPEC: .T. = write out spectral data on historical file

INTEGER(KIND=JPIM) :: N1POS
INTEGER(KIND=JPIM) :: N1HIS
INTEGER(KIND=JPIM) :: N1GDI
INTEGER(KIND=JPIM) :: N1SDI
INTEGER(KIND=JPIM) :: N1DHP
INTEGER(KIND=JPIM) :: N1RES
INTEGER(KIND=JPIM) :: N1XFU
INTEGER(KIND=JPIM) :: N1DHFG
INTEGER(KIND=JPIM) :: N1DHFZ
INTEGER(KIND=JPIM) :: N1DHFD
INTEGER(KIND=JPIM) :: N1CFU
INTEGER(KIND=JPIM) :: N1SFXHIS
INTEGER(KIND=JPIM) :: N1MASSCON
LOGICAL :: LRFILAF
LOGICAL :: LWRSPEC

!     ------------------------------------------------------------------
END MODULE YOMCT1
