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

MODULE YOMMTS

IMPLICIT NONE

SAVE

! Handling of simulated satellite images (former yommts)
! (Cles d'activation de la production de temperatures de brillance)

! Overall key :
LOGICAL :: LMTS = .FALSE.
! Keys per satellites :
LOGICAL :: LCHAN_MSAT(7:11,1:8)   = .FALSE.
LOGICAL :: LCHAN_GOES(11:17,1:10) = .FALSE.
LOGICAL :: LCHAN_MTSAT(1:1,1:4)   = .FALSE.
LOGICAL :: LCHAN_HIMA(8:8,1:10)   = .FALSE.

END MODULE YOMMTS
