MODULE YOMLUN1S
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD,JPRM
IMPLICIT NONE
SAVE

!     ------------------------------------------------------------------

!*    Logical units used by code

!     RMISS  :   missing data value
!     NULOUT :   output unit
!     NULNAM :   unit number for namelist
!     NULFOR :   unit number for forcing data
!     NPOSRES:   unit number for restart file
!     NPOSGG :   unit number for post-processing of prognostic variables
!     NPOSEFL:   unit number for post-processing of energy fluxes
!     NPOSWAT:   unit number for post-processing of water budget
!     NPOSSUS:   unit number for post-processing of surface state
!     NPOSSUB:   unit number for post-processing of sub-surface state
!     NPOSEVA:   unit number for post-processing of evaporation comp
!     NPOSCLD:   unit number for cold-season processes
!     NPOSCLM:   unit number for fixed climate fields
!     NPOSDFO:   unit number for post-processing of diagnostic variables
!                             (forcing)
!     NPOSDBD:   unit number for post-processing of diagnostic variables
!                             (forcing)
!     NPOSDTIk:   unit number for post-processing of diagnostic variables
!                             (tile k budget)
!     NPOSDST:   unit number for post-processing of diagnostic variables
!                             (snow T budget)
!     NPOSDT0:   unit number for post-processing of diagnostic variables
!                             (skin temperature budget)
!     NPOSDT1:   unit number for post-processing of diagnostic variables
!                             (temperature layer 1 budget)
!     NPOSDT2:   unit number for post-processing of diagnostic variables
!                             (temperature layer 2 budget)
!     NPOSDT3:   unit number for post-processing of diagnostic variables
!                             (temperature layer 3 budget)
!     NPOSDT4:   unit number for post-processing of diagnostic variables
!                             (temperature layer 4 budget)
!     NPOSDSW:   unit number for post-processing of diagnostic variables
!                             (snow budget)
!     NPOSDW0:   unit number for post-processing of diagnostic variables
!                             (interception water budget)
!     NPOSDW1:   unit number for post-processing of diagnostic variables
!                             (soil water layer 1 budget)
!     NPOSDW2:   unit number for post-processing of diagnostic variables
!                             (soil water layer 2 budget)
!     NPOSDW3:   unit number for post-processing of diagnostic variables
!                             (soil water layer 3 budget)
!     NPOSDW4:   unit number for post-processing of diagnostic variables
!                             (soil water layer 4 budget)
!     NPOSDMO:   unit number for post-processing of diagnostic variables
!                             (momentum fluxes)
!     NULGP0:    unit number for initial prognostic fields (yomgp1s.h)
!     NULGPD0:   unit number for initial geographic fields (yomgpd1s.h)
!     NPOSRC:   unit number for post-processing of diagnostic variables
!                             (canopy resistance)
!     NPOSOCP:   unit number for ocean mixed layer model prog. variables.
!     NPOSOCD:   unit number for ocean mixed layer model diag. variables.

!     NPOSCO2:   unit number for post-processing of CO2 fluxes
!     NPOSBIO:   unit number for post-processing of biomass
!     NPOSVEG:   unit number for post-processing of vegetation variables
!     NPOSFRA:   unit number for post-processing of fraction
!     NPOSCLM:   unit number for post-processing of fixed climate fields
!     NPOSTIL:   unit number for post-processing of tiled output
!     NPOSVTY:   unit number for post-processing of output per vegetation type
!     NPOSFOR:   unit number for forcing fields
!     NPOSD2M:   unit number for post-processing of diagnostics as 2m and 10m 
!     NPOSGGD:   unit number for post-processing of prognostics mean


INTEGER(KIND=JPIM) :: NPOSRES
INTEGER(KIND=JPIM) :: NULOUT
INTEGER(KIND=JPIM) :: NULNAM
INTEGER(KIND=JPIM) :: NULFOR
INTEGER(KIND=JPIM) :: NPOSGG
INTEGER(KIND=JPIM) :: NPOSGGD
INTEGER(KIND=JPIM) :: NPOSEFL
INTEGER(KIND=JPIM) :: NPOSWAT
INTEGER(KIND=JPIM) :: NPOSSUS
INTEGER(KIND=JPIM) :: NPOSSUB
INTEGER(KIND=JPIM) :: NPOSEVA
INTEGER(KIND=JPIM) :: NPOSCLD
INTEGER(KIND=JPIM) :: NPOSCLM
INTEGER(KIND=JPIM) :: NPOSDFO
INTEGER(KIND=JPIM) :: NPOSDBD
INTEGER(KIND=JPIM) :: NPOSDTI1
INTEGER(KIND=JPIM) :: NPOSDTI2
INTEGER(KIND=JPIM) :: NPOSDTI3
INTEGER(KIND=JPIM) :: NPOSDTI4
INTEGER(KIND=JPIM) :: NPOSDTI5
INTEGER(KIND=JPIM) :: NPOSDTI6
INTEGER(KIND=JPIM) :: NPOSDTI7
INTEGER(KIND=JPIM) :: NPOSDTI8
INTEGER(KIND=JPIM) :: NPOSDT0
INTEGER(KIND=JPIM) :: NPOSDST
INTEGER(KIND=JPIM) :: NPOSDT1
INTEGER(KIND=JPIM) :: NPOSDT2
INTEGER(KIND=JPIM) :: NPOSDT3
INTEGER(KIND=JPIM) :: NPOSDT4
INTEGER(KIND=JPIM) :: NPOSDSW
INTEGER(KIND=JPIM) :: NPOSDW0
INTEGER(KIND=JPIM) :: NPOSDW1
INTEGER(KIND=JPIM) :: NPOSDW2
INTEGER(KIND=JPIM) :: NPOSDW3
INTEGER(KIND=JPIM) :: NPOSDW4
INTEGER(KIND=JPIM) :: NPOSDMO
INTEGER(KIND=JPIM) :: NULGP0
INTEGER(KIND=JPIM) :: NULGPD0
INTEGER(KIND=JPIM) :: NPOSRC

INTEGER(KIND=JPIM) :: NPOSLKE  ! E. DUTRA  LAKE DIGANGOSTICS FILE (NETCDF)
INTEGER(KIND=JPIM) :: NPOSGGL  ! E. DUTRA LAKE PROGNOSTIC VARIABLES 
INTEGER(KIND=JPIM) :: NPOSDTI9 ! E. DUTRA  TILE 9 SKIN DIAGNOSTICS 

INTEGER(KIND=JPIM) :: NPOSOCP  !KPP
INTEGER(KIND=JPIM) :: NPOSOCD  !KPP
INTEGER(KIND=JPIM) :: NPOSRESO !KPP

INTEGER(KIND=JPIM) :: NPOSCO2
INTEGER(KIND=JPIM) :: NPOSBIO
INTEGER(KIND=JPIM) :: NPOSVEG
INTEGER(KIND=JPIM) :: NPOSEXT
INTEGER(KIND=JPIM) :: NPOSTIL
INTEGER(KIND=JPIM) :: NPOSVTY
INTEGER(KIND=JPIM) :: NPOSFOR
INTEGER(KIND=JPIM) :: NPOSD2M


REAL(KIND=JPRB) :: RMISS  ! missing value - single/double
REAL(KIND=JPRD) :: DMISS  ! real double missing values 
REAL(KIND=JPRM) :: ZMISS  ! missing value - real 
INTEGER*4 :: IMISS ! missing value - integer 

CHARACTER(LEN=100)     :: CTUNITS  !! time information for netcdf file 
INTEGER(KIND=JPIM)     :: NCTYPE   !! TYPE of NETCDF FILE 
!     ------------------------------------------------------------------
END MODULE YOMLUN1S
