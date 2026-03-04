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

MODULE YOMCVMNH

!Module to define the parameters of the KFB
!convection parametrization of MNH.
! Modified by E. Bazile 09.05.2007
!    parameters for KFB in namelist

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

TYPE TCVMNH
!--------------------------------------------

!LDEEP : Control Key to deep convection scheme
!LSHALLOW : Control key to shallow convection scheme
!LDIAGCONV : Control key to compute diagnostics on convection
!LSETTADJ : logical to set convective adjustment time by user
!LREFRESH_ALL ! refresh or not tendencies at every call
!LDOWN ! take or not convective downdrafts into account 

!OTADJD : user specified time step to call deep convection
!OTADJS : user specified time step to call shallow convection

!NSETENS : use of several calls to deep_conv (0 only one call)
!NIICE : USe of Ice in convection 1 yes 0 no

LOGICAL :: LDEEP, LSHALLOW, LDIAGCONV
LOGICAL :: LSETTADJ, LREFRESH_ALL, LDOWN, LSMOOTH
REAL(KIND=JPRB) :: OTADJD
REAL(KIND=JPRB) :: OTADJS
INTEGER(KIND=JPIM) :: NSETENS, NIICE
REAL(KIND=JPRB) :: XA25
REAL(KIND=JPRB) :: XCRAD
REAL(KIND=JPRB) :: XCDEPTH
REAL(KIND=JPRB) :: XCDEPTH_D
REAL(KIND=JPRB) :: XDTPERT
REAL(KIND=JPRB) :: XENTR
REAL(KIND=JPRB) :: XZLCL
REAL(KIND=JPRB) :: XZPBL
REAL(KIND=JPRB) :: XWTRIG
REAL(KIND=JPRB) :: XNHGAM
REAL(KIND=JPRB) :: XTFRZ1
REAL(KIND=JPRB) :: XTFRZ2 
REAL(KIND=JPRB) :: XSTABT
REAL(KIND=JPRB) :: XSTABC
REAL(KIND=JPRB) :: XAW
REAL(KIND=JPRB) :: XBW
REAL(KIND=JPRB) :: XATPERT
REAL(KIND=JPRB) :: XBTPERT

END TYPE TCVMNH

!!TYPE(TCVMNH), POINTER :: YRCVMNH => NULL()

END MODULE YOMCVMNH
