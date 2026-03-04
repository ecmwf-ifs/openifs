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

MODULE YOMNUD

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE PARDIM, ONLY : JPMXLE

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------

! LNUDG  : USE OF NUDGING TERMS

! LNUDTE : nudging switch for temperature (altitude)
! LNUDSH : nudging switch for specific humidity (altitude)
! LNUDDI : nudging switch for divergence (altitude)
! LNUDVO : nudging switch for vorticity (altitude)
! LNUDSV : nudging switch for scalar variables (altitude)
! LNUDLP : nudging switch for log(surf pressure)
! LNUDST : nudging switch for surface temperature
! LNUDSM : nudging switch for surface moisture
! LNUDRM : nudging switch for root moisture
! LNUDSD : nudging switch for snow depth

! LNDRIV : nudging switch for using horizontal mask 

! LNUDTG : nudging switch for grid point temperature
! LNUDQG : nudging switch for grid point specific humidity
! LNUDUG : nudging switch for grid point U-wind
! LNUDVG : nudging switch for grid point V-wind
! LNUDPG : nudging switch for grid point surface pressure

! LWNUDG : nudging switch for time variable coefficient


! NFNUDG : number of analysis time steps in memory
! NFRNUDG: frequency of nudging (time-steps)
! NTNUDG : first vertical level where nudging is active
! NTOTFNUDG3 : number of 3D fields when grid point nudging
! NTOTFNUDG2 : number of 2D fields when grid point nudging

! NSPNU1 : wave number below which full nudging is applied
! NSPNU2 : wave number beyond which no nudging is applied

! XNUDTE : nudging coefficient for temperature (altitude)
! XNUDSH : nudging coefficient for specific humidity (altitude)
! XNUDDI : nudging coefficient for divergence (altitude)
! XNUDVO : nudging coefficient for vorticity (altitude)
! XNUDSV : nudging coefficient for scalar variables (altitude)
! XNUDLP : nudging coefficient for log(surf pressure)
! XNUDST : nudging coefficient for surface temperature
! XNUDSM : nudging coefficient for surface moisture
! XNUDRM : nudging coefficient for root moisture
! XNUDSD : nudging coefficient for snow depth

! XNUDTG : nudging coefficient for grid point temperature
! XNUDQG : nudging coefficient for grid point specific humidity
! XNUDUG : nudging coefficient for grid point U-wind
! XNUDVG : nudging coefficient for grid point V-wind
! XNUDPG : nudging coefficient for grid point surface pressure

! XNUVERT(NFLEVG) : vertical profile of nudging coefficient

INTEGER(KIND=JPIM) :: NFNUDG
INTEGER(KIND=JPIM) :: NFRNUDG
INTEGER(KIND=JPIM) :: NTNUDG
INTEGER(KIND=JPIM) :: NTOTFNUDG3
INTEGER(KIND=JPIM) :: NTOTFNUDG2
INTEGER(KIND=JPIM) :: NSPNU1 
INTEGER(KIND=JPIM) :: NSPNU2
LOGICAL :: LNUDG
LOGICAL :: LNUDDI
LOGICAL :: LNUDLP
LOGICAL :: LNUDRM
LOGICAL :: LNUDSD
LOGICAL :: LNUDSH
LOGICAL :: LNUDSM
LOGICAL :: LNUDST
LOGICAL :: LNUDSV
LOGICAL :: LNUDTE
LOGICAL :: LNUDVO
LOGICAL :: LNDRIV 

LOGICAL :: LNUDTG 
LOGICAL :: LNUDQG 
LOGICAL :: LNUDUG 
LOGICAL :: LNUDVG 
LOGICAL :: LNUDPG 

LOGICAL :: LWNUDG

REAL(KIND=JPRB) :: XNUVERT(1:JPMXLE)
REAL(KIND=JPRB) :: XNUDDI
REAL(KIND=JPRB) :: XNUDLP
REAL(KIND=JPRB) :: XNUDRM
REAL(KIND=JPRB) :: XNUDSD
REAL(KIND=JPRB) :: XNUDSH
REAL(KIND=JPRB) :: XNUDSM
REAL(KIND=JPRB) :: XNUDST
REAL(KIND=JPRB) :: XNUDSV
REAL(KIND=JPRB) :: XNUDTE
REAL(KIND=JPRB) :: XNUDVO
REAL(KIND=JPRB) :: XNUDTG 
REAL(KIND=JPRB) :: XNUDQG 
REAL(KIND=JPRB) :: XNUDUG 
REAL(KIND=JPRB) :: XNUDVG 
REAL(KIND=JPRB) :: XNUDPG 

REAL(KIND=JPRB),ALLOCATABLE :: XVUST(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: XVUSM(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: XVURM(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: XVUSD(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: XVUTG(:,:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: XVUQG(:,:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: XVUUG(:,:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: XVUVG(:,:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: XVUPG(:,:,:)


END MODULE YOMNUD
