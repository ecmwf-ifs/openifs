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

MODULE YOMIO_SERV_MAP_PLAN

!**** *YOMIO_SERV_MAP_PLAN*  - IO server plans (send and recv)

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE* 
!      Original : 01-10-2014


USE PARKIND1, ONLY : JPIM, JPRB
USE IOMULTIBUF_MOD, ONLY : IOMULTIBUF

IMPLICIT NONE

TYPE IO_SERV_SEND_PLAN
  TYPE (IOMULTIBUF),   POINTER :: YLBUFA (:) => NULL () ! Full
  TYPE (IOMULTIBUF),   POINTER :: YLBUFD (:) => NULL () ! Field data
  INTEGER (KIND=JPIM), POINTER :: IPRIO  (:) => NULL () ! IO server ranks
  INTEGER (KIND=JPIM)          :: INIO  = 0_JPIM        ! Number of IO ranks involved
  INTEGER (KIND=JPIM)          :: IFNUM = 0_JPIM        ! Number of fields
END TYPE IO_SERV_SEND_PLAN

TYPE IO_SERV_RECV_PLAN
  TYPE (IOMULTIBUF),   POINTER :: YLBUFA (:)   => NULL () ! Headers
  REAL (KIND=JPRB),    POINTER :: ZLBUFD (:,:) => NULL () ! Field data
  INTEGER (KIND=JPIM)          :: INIO  = 0_JPIM          ! Number of IO ranks involved
  INTEGER (KIND=JPIM)          :: IFNUM = 0_JPIM          ! Number of fields
  INTEGER (KIND=JPIM), POINTER :: IREQA (:) => NULL ()    ! Requests IDs for YLBUFA
  INTEGER (KIND=JPIM), POINTER :: IREQD (:) => NULL ()    ! Requests IDs for ZLBUFD
END TYPE IO_SERV_RECV_PLAN

SAVE

END MODULE YOMIO_SERV_MAP_PLAN

