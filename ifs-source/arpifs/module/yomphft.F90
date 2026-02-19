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

MODULE YOMPHFT
!     Purpose.
!     --------
!     Declaration of: 
!     - type TYPE_APFT, used for descriptions in DDH
!     - variables and arrays needed to transport physical fluxes
!       and tendencies from physical routines to DDH 

!     Author.
!     -------
!      Tomislav Kovacic

!     Modifications.
!     --------------
!      Original : 2006-03-23
!      Modifications:
!      2007-03-16 T.Kovacic APFT four dimensinal
!-------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM,JPRB

IMPLICIT NONE

SAVE

TYPE TYPE_APFT
  CHARACTER(LEN=1)   :: CFT
  CHARACTER(LEN=2)   :: CVAR
  CHARACTER(LEN=10)  :: CNAME
END TYPE TYPE_APFT


INTEGER(KIND=JPIM)              :: NPROGVAR ! number of prognostic variables
                                            ! changed by physical processes
INTEGER(KIND=JPIM), ALLOCATABLE :: MJJ1(:) !first index in APFT for variable 'i'
INTEGER(KIND=JPIM), ALLOCATABLE :: MJJ2(:) !last index in APFT for variable 'i'
INTEGER(KIND=JPIM)              :: NAPHFT  !number of fluxes/tendencies in APFT
INTEGER(KIND=JPIM)              :: NDDHFT  !number of fluxes/tendencies in YAPFT

TYPE(TYPE_APFT),    ALLOCATABLE :: YAPFT(:)!descriptions of fluxes & tendencies

REAL(KIND=JPRB), ALLOCATABLE :: APFT(:,:,:,:)  ! array with fluxes & tendencies 
                                               !1. dimension: horisontal extend
                                               !2. dimansion: vertical extend
                                               !3. dimension: NAPHFT
                                               !4. dimension: thread

!  Positions of momentum turbulent fluxes
INTEGER(KIND=JPIM)  :: NUUTURFT ! third index in APFT for U turb. flux
INTEGER(KIND=JPIM)  :: NVVTURFT ! third index in APFT for V turb. flux

!  Positions of convective fluxes
INTEGER(KIND=JPIM)  :: NQLCONV ! third index in APFT for convective rain 
INTEGER(KIND=JPIM)  :: NQNCONV ! third index in APFT for convective snow

INTEGER(KIND=JPIM)  :: NCDPIPHFT ! number of fluxes/tendencies in CDPI

!  Positions of fluxes for CDPI
                               ! third index in APFT for
INTEGER(KIND=JPIM)  :: NQVQL1 ! Pl'
INTEGER(KIND=JPIM)  :: NQVQI1 ! Pi'
INTEGER(KIND=JPIM)  :: NQLQR2 ! Pl''
INTEGER(KIND=JPIM)  :: NQIQS2 ! Pi''
INTEGER(KIND=JPIM)  :: NQRQV3 ! Pl'''
INTEGER(KIND=JPIM)  :: NQSQV3 ! Pi'''
INTEGER(KIND=JPIM)  :: NQGQV3 ! Pg'''
INTEGER(KIND=JPIM)  :: NQR0   ! Pr
INTEGER(KIND=JPIM)  :: NQI0   ! Pi
INTEGER(KIND=JPIM)  :: NQS0   ! Ps
INTEGER(KIND=JPIM)  :: NQG0   ! Pg


END MODULE YOMPHFT
