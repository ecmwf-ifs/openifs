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

MODULE YOMRANDOM_STREAMS

!**** *YOMRANDOM_STREAMS*  - Globally-available random number streams 

!  Original :      2004-05-06 - Mike Fisher
!  Modifications: 07-09-2006 - Judith Berner
!                              added random number streams for stochastic physics
!                  Jan-2016   - SJ Lock : removed LSTOPH option   
                  

USE RANDOM_NUMBERS_MIX, ONLY: RANDOMNUMBERSTREAM

IMPLICIT NONE

SAVE

TYPE :: TRANDOM_STREAMS
TYPE (RANDOMNUMBERSTREAM) :: SCAN2MTL
TYPE (RANDOMNUMBERSTREAM) :: STOCHPHYS_CABS
TYPE (RANDOMNUMBERSTREAM) :: STOCHPHYS_SPBS
TYPE (RANDOMNUMBERSTREAM) :: STOCHPHYS_RVP
TYPE (RANDOMNUMBERSTREAM) :: STOPH_SPBS_T
TYPE (RANDOMNUMBERSTREAM) :: STOPH_RVP_T
END TYPE TRANDOM_STREAMS

!!TYPE(TRANDOM_STREAMS), POINTER :: YR_RANDOM_STREAMS => NULL()

END MODULE YOMRANDOM_STREAMS
