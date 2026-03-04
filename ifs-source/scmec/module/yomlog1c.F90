! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

MODULE YOMLOG1C

!    SWITCHES RELATED TO ONE-COLUMN MODEL

! NAME      TYPE       DEFAULT   PURPOSE
! ------- : -------- : ------- : ----------------------------------------------
! LDYNFOR : LOGICAL  : TRUE    : TRUE - if dynamical calculation done producing
!                                       dynamical tendencies.
! LUGVG   : LOGICAL  : TRUE    : TRUE - if ug and vg are used.
! LVERVEL : LOGICAL  : TRUE    : TRUE - if omega ("vervel") is used (has to be 
!                                on even if LETADOT, as w is used in callpar).
! LETADOT : LOGICAL  : FALSE   : TRUE - if etadot * dp/deta is used for 
!                                vertical advection.
! LUPWIND : LOGICAL  : TRUE    : TRUE - if upwind and forward in time 
!                                differencing is used for vertical advection.
!                                FALSE- if the Lax forward in time, centered
!                                in space scheme is used for vert. advection.
! LWADV   : LOGICAL  : TRUE    : TRUE - if vertical advection is used.
! LWADVCLD: LOGICAL  : FALSE   : TRUE - if vertical advection of clouds (a,l,i) is used.
! LTADV   : LOGICAL  : FALSE   : TRUE - if hor.adv. of t is used.
! LQADV   : LOGICAL  : FALSE   : TRUE - if hor.adv. of q is used.
! LUVADV  : LOGICAL  : FALSE   : TRUE - if hor.adv. of u,v is used.
! LVARFOR : LOGICAL  : FALSE   : TRUE - if varying large-scale forcing.
! LVARSST : LOGICAL  : FALSE   : TRUE - if varying lat, lon, SST, surface 
!                                pressure.
! LRELAX  : LOGICAL  : FALSE   : TRUE - if relaxation is used.
! LUVREL  : LOGICAL  : TRUE    : TRUE - if relaxation for u and v is used.
! LTQREL  : LOGICAL  : TRUE    : TRUE - if TRUE - if relaxation for t and q.
! RUVREL_PLEV : REAL : 2000hPa : only pressure levels above this are relaxed UV
! RTQREL_PLEV : REAL : 2000hPa : only pressure levels above this are relaxed TQ
! NFRFOR  : INTEGER  : 1       : frequency of large-scale forcing input 
!                                [steps].
! NFRSST  : INTEGER  : 1       : frequency of changing SST [steps].
! NFROBS  : INTEGER  : 1       : frequency of changing relaxation values 
!                                [steps].
! NSTRTINI: INTEGER  : 1       : # of first step in ini. cond./forc. data used.
! NPOSPRG : INTEGER  :         : output unit number (prognostic variables)
! NPOSDIA : INTEGER  :         : output unit number (diagnostic variables)
! NPOSDIA2: INTEGER  :         : output unit number (diagnostic variables - 2)
! NPOSASC : INTEGER  :         : output unit number (ascii output)
! RDTRELAX: REAL     : 10800   : dt for relaxation [s].
! OUTFORM : CHARACTER: netcdf  : identifies the output data format ("ascii" or 
!                                "netcdf")
! CMODID  : CHARACTER:         : model code identification
! CSIMID  : CHARACTER:         : simulation identification

!     -----------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

LOGICAL LDYNFOR
LOGICAL LUGVG
LOGICAL LVERVEL
LOGICAL LETADOT
LOGICAL LUPWIND
LOGICAL LWADV
LOGICAL LWADVCLD
LOGICAL LTADV
LOGICAL LQADV
LOGICAL LUVADV
LOGICAL LVARFOR
LOGICAL LVARSST
LOGICAL LRELAX
LOGICAL LUVREL
LOGICAL LTQREL
INTEGER(KIND=JPIM) :: NFRFOR
INTEGER(KIND=JPIM) :: NFRSST
INTEGER(KIND=JPIM) :: NFROBS
INTEGER(KIND=JPIM) :: NSTRTINI
INTEGER(KIND=JPIM) :: NPOSPRG
INTEGER(KIND=JPIM) :: NPOSDIA
INTEGER(KIND=JPIM) :: NPOSDIA2
INTEGER(KIND=JPIM) :: NPOSASC
REAL(KIND=JPRB) :: RDTRELAX
REAL(KIND=JPRB) :: RUVREL_PLEV
REAL(KIND=JPRB) :: RTQREL_PLEV
CHARACTER*100 OUTFORM
CHARACTER*100 CMODID
CHARACTER*100 CSIMID

!     -----------------------------------------------------------------

END MODULE YOMLOG1C
