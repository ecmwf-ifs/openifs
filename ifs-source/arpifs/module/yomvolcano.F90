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

MODULE YOMVOLCANO

! volcano emissions, dates and source strength

USE PARKIND1  ,ONLY : JPRB, JPIM

IMPLICIT NONE
SAVE
! volcano emissions
INTEGER(KIND=JPIM),PARAMETER :: JVOCDAT=200
INTEGER(KIND=JPIM)           ::  NVOCDATES
INTEGER(KIND=JPIM), DIMENSION(JVOCDAT) :: IVOCSTART ! start end date for specific source strength&height
LOGICAL            :: LVOCENS   ! tracer ensemble

REAL(KIND=JPRB),DIMENSION(JVOCDAT) :: SLVOC1, SLVOC2 ! injection height range in km  ! namelist input 
REAL(KIND=JPRB),DIMENSION(10) :: SLVOCES1, SLVOCES2 ! injection height range in km for 10 tracers ensemble
REAL(KIND=JPRB) , DIMENSION(JVOCDAT) :: SEMIVOC   ! kg/s
REAL(KIND=JPRB) , DIMENSION(JVOCDAT) :: SEMIVOCFLX   ! kg/m**2
INTEGER(KIND=JPIM),DIMENSION(JVOCDAT) :: ILVOCJB ! JB injection height model levels 
REAL(KIND=JPRB)  :: SEMIVOCENS   ! kg/s
REAL(KIND=JPRB)   :: SEMIVOCFLXENS   ! kg/m**2s
! volcano location
LOGICAL            :: LVOCMP   ! volcano on this mpi task
INTEGER(KIND=JPIM) :: IVOCGP   ! index of volcano location betwen 1 NPROMA
INTEGER(KIND=JPIM) :: IVOCBLK  ! block position of volcanoe in NGPTOT  
REAL(KIND=JPRB)    :: SVOCLAT, SVOCLON ! lat lon of volcano

!     -------------------------------------------------------------------------
END MODULE YOMVOLCANO
