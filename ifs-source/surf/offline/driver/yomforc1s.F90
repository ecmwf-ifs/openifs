MODULE YOMFORC1S
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD

IMPLICIT NONE
SAVE

!     -----------------------------------------------

!*    ATMOSPHERIC FORCING DATA FOR 1D SURFACE SCHEME, PILPS VERSION

!     -----------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPTYFC=12

REAL(KIND=JPRB),ALLOCATABLE,TARGET:: GFOR(:,:,:)
REAL(KIND=JPRB),POINTER :: UFI(:,:)
REAL(KIND=JPRB),POINTER :: VFI(:,:)
REAL(KIND=JPRB),POINTER :: TFI(:,:)
REAL(KIND=JPRB),POINTER :: QFI(:,:)
REAL(KIND=JPRB),POINTER :: CO2FI(:,:)
REAL(KIND=JPRB),POINTER :: PSFI(:,:)
REAL(KIND=JPRB),POINTER :: SRFFI(:,:)
REAL(KIND=JPRB),POINTER :: TRFFI(:,:)
REAL(KIND=JPRB),POINTER :: R30FI(:,:)
REAL(KIND=JPRB),POINTER :: S30FI(:,:)
REAL(KIND=JPRB),POINTER :: R30FI_C(:,:)
REAL(KIND=JPRB),POINTER :: S30FI_C(:,:)
REAL(KIND=JPRD) :: DTIMFC                ! Forcing frequency 
REAL(KIND=JPRD) :: RTSTFC                ! first step of loaded forcing ref. time 
INTEGER(KIND=JPIM) :: NSTPFC
INTEGER(KIND=JPIM) :: JPSTPFC           ! Dimension of the forcing set in the namelist
INTEGER(KIND=JPIM) :: DIMFORC            ! forcing dimension 1/2
! REAL(KIND=JPRB) :: PHISTA

!------------------------------------------------------------     
END MODULE YOMFORC1S
